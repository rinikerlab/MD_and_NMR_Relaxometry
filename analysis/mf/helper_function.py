#!/usr/bin/env python3

from pyDR.Selection import select_tools as selt
import fitting_functions as fit

import MDAnalysis as mda

import copy
import numpy as np
import pandas as pd

import emcee


def find_slice_indices(topo, traj, block_size):
    """
    This function will return a list of indices from where to start/stop each block
    to perform an iRED analysis averaged over many blocks. The block size should be chosen 
    to be ~ 5 times the tumbling time of the protein in solution. 
    
    Note: The function may also be used for block averaged Model Free Analysis 

    https://pubs.acs.org/doi/pdf/10.1021/ct500181v
    
    The functions figures this out by creating an MDAnalysis trajectory and gathering the information from it. 
    
    Arguments
    ---------
        topo:
            path to topology (.gro)
        traj:
            path to trajectory (.traj)
        block_size: 
            size of the block in [ns] 
    
    Returns
    --------
        indices_blocks: List [Tuple(int, int)]
            list of indices of the frames for which we want to start and end each block
    """
    u = mda.Universe(topo, traj)
    
    n_frames_per_block = int(block_size * 1000 / u.trajectory.dt)
    n_frames_tot = u.trajectory.n_frames    
    
    return [(n_frames_per_block*i , n_frames_per_block*(i+1)) for i in range(int(n_frames_tot/n_frames_per_block))]


def find_cutoff_index(time, t_cutoff):
    """
    Find the index where time matches t_cutoff (ns in Andy's frames object)
    which allows to select the cutoff for the correlation function
    """
    return np.argmin(time < t_cutoff)


def perform_different_fits(time, ct):
    """
    This function performs the two types of fits for the correlation function given. 
    """
    
    # Initialize them all to None
    soln_bh = soln_mcmc = None
        
    # 1 : Perform basin hopping fit
    
    print ('Start Basin Hopping')
    try:
        
        soln_bh = fit.fit_like_ferrage_2018(time, ct)
        ct_bh = fit.extended_model_free(time, *soln_bh.x)
        
        soln_bh.labels = ['Ss2', 'tau_s', 'Sf2', 'tau_f']
        
        # Get the "error" of the fit
        error = np.sqrt(np.power(ct_bh - ct, 2))
    except:
        print ('Basin Hopping fit returned an error')
    
    # 2: MCMC
   
    sampler = None
    ct_mcmc = None

    #try:
    nWalkers = 50
    nParam = len(soln_bh.x) + 1
    nMcmc = 2500

    bounds = fit.gen_bounds_mcmc(2)
    pos = fit.intialize_mcmc_walkers(bounds, nWalkers)
    
    print ('Start MCMC') 
    sampler = emcee.EnsembleSampler(nWalkers, nParam, fit.log_probability, args=(time, ct, error))
    sampler.run_mcmc(pos, nMcmc, progress=True)

    # Extract results
    discard_first_n = int(nMcmc / 2)
    samples = sampler.chain[:, discard_first_n:, :].reshape((-1, nParam))
    mean = np.percentile(samples, 50, axis=0)
    sixteenth = np.percentile(samples, 16, axis=0)
    heightyfourth = np.percentile(samples, 84, axis=0)

    soln_mcmc = pd.DataFrame(data = np.array([mean, sixteenth, heightyfourth]).T, 
                            columns = ['mean', 'sixteenth', 'heightyforth'], 
                            index = ['Ss2', 'tau_s', 'Sf2', 'tau_f', 'log_f'] )
    
    ct_mcmc = fit.extended_model_free(time, *mean[:-1])

    #except:
    #    print ('Error during the MCMC')
    
    return soln_bh, ct_bh, soln_mcmc, ct_mcmc


def perform_monoexponential_fits(time, ct):
    """
    This function performs the two types of fits for the correlation function given. 
    """

    # Initialize them all to None
    soln_bh = soln_mcmc = None

    # 1 : Perform basin hopping fit

    print ('Start Basin Hopping')
    try:

        soln_bh = fit.fit_like_ferrage_monoexp(time, ct)
        ct_bh = fit.simple_model_free(time, *soln_bh.x)

        soln_bh.labels = ['S2', 'tau']

        # Get the "error" of the fit
        error = np.sqrt(np.power(ct_bh - ct, 2))
    except:
        print ('Basin Hopping fit returned an error')

    # 2: MCMC
    """
    sampler = None
    ct_mcmc = None

    #try:
    nWalkers = 50
    nParam = len(soln_bh.x) + 1
    nMcmc = 2500

    bounds = fit.gen_bounds_mcmc(1)
    pos = fit.intialize_mcmc_walkers(bounds, nWalkers)

    print ('Start MCMC')
    sampler = emcee.EnsembleSampler(nWalkers, nParam, fit.log_probability, args=(time, ct, error))
    sampler.run_mcmc(pos, nMcmc, progress=True)

    # Extract results
    discard_first_n = int(nMcmc / 2)
    samples = sampler.chain[:, discard_first_n:, :].reshape((-1, nParam))
    mean = np.percentile(samples, 50, axis=0)
    sixteenth = np.percentile(samples, 16, axis=0)
    heightyfourth = np.percentile(samples, 84, axis=0)

    soln_mcmc = pd.DataFrame(data = np.array([mean, sixteenth, heightyfourth]).T,
                            columns = ['mean', 'sixteenth', 'heightyforth'],
                            index = ['S2', 'tau', 'log_f'] )

    ct_mcmc = fit.simple_model_free(time, *mean[:-1])
    """
    #except:
    #    print ('Error during the MCMC')

    return soln_bh #, ct_bh, soln_mcmc, ct_mcmc



