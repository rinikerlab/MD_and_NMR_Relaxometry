#!/usr/bin/env python3

import pyDR # this notebook worked well with git commit: 15b82d57ec10193603d9486500925ff5249747f0
from pyDR.Selection import select_tools as selt

import emcee
from scipy import optimize

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import MDAnalysis as mda

import pickle, copy, sys, os

def save_pickle(out_path, data):
    with open(out_path, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_pickle(pickle_path):
    with open(pickle_path, 'rb') as handle:
        return pickle.load(handle)

from my_selection import select_caxis_iles, selection_all_methyls

# These were useful when I was compared the effect on the fit.
#BLOCK_SIZE = int(sys.argv[2])  we found chunks of 1 mus were good
#CUTOFF_TIME = int(sys.argv[3]) we found 50 ns was good

def plotCorrelationFunctions(t, cts_raw, cts_fit, cutoff, labels, title, out_path):
    """
    Plots the correlation functions, to make sure the fit is ok
    
    t: time
    
    cts_raw: list of correlation functions directly calculated from MD
    cts_fit: list of correlation functions fit with the MCMC
    
    labels: labels of the isoleucines (or other residue)
    title: title to give the plot
    
    """
    from matplotlib.pyplot import cm

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = [8, 6])
    color = iter(cm.rainbow(np.linspace(0, 1, len(cts_raw))))

    for i, (ct, ct_fit) in enumerate(zip(cts_raw, cts_fit)):

        c = next(color)

        ax.plot(t[:cutoff], ct[:cutoff], label = f'res {labels[i]}', c=c, alpha = 0.5)
        ax.plot(t[:cutoff], ct_fit, ls = '--', c=c, lw = 1)
        # break

    ax.set_ylim([0, 1])

    ax.set_ylabel('C(t) [a.u.]')
    ax.set_xlabel('time [ns]')

    ax.set_title(title)
    ax.legend()

    fig.savefig(out_path)
    plt.close()
    return 
    

def blockaveraged_ModelFree(top_path, traj, block_size, step=10, cutoff_time = 50, out_dir = None, label_ff = ''):
    """
    This function will perform a model free analysis by dividing the trajectory into different chunks, 
    
    1) extracting the correlation functions with a specific cutoff, 
    2) fitting them with the Basin-Hopping / MCMC approach used by Fabien
    3) reporting results with average / std edviations among those blocks. 
    
        Arguments
    ---------
        top_path:
            path to topology (.gro)
        traj:
            path to trajectory (.traj)
        block_size: 
            size of the block in [ns]
        step: 
            stride over the data (makes it faster)
        cutoff_time:
            time [ns] at which to cutoff the correlation functions before the fit
        out_dir:
            direcotry to print the results in
        ff_name:
            name of the ff 
        
        
    
    Returns
    --------
        Write it here
    
    """
    
    
    # 1: Find how to separate our data:
    from helper_function import find_slice_indices, perform_different_fits, find_cutoff_index
    
    idx_blocks = find_slice_indices(top_path, traj, block_size=block_size)

    # 2: Perform an iRED analysis on each of those blocks

    data_all_blocks = []

    for i, (start, end) in enumerate(idx_blocks):
        print (str(start) + '\t --> \t' + str(end))
        
        molsys = pyDR.Selection.MolSys.MolSys(top_path, traj, t0 = start, tf = end, step = step) 
        molsel = pyDR.Selection.MolSys.MolSelect(molsys)
        
        # 1) Extract correlation functions
        
        selection_all_methyls(molsys, molsel)

        frames = pyDR.Frames.FrameObj(molsel, )

        try:
            ct_data = frames.md2data()
        except:
            pass
        
        
        # 2) Fit correlation functions
        
        cut_idx = find_cutoff_index(frames.t, 50)

        full_data = {}

        cts_fit = []
        
        for k, (label, ct_full) in enumerate(zip(ct_data.label, frames.Ct['ct'])):
            time = frames.t[:cut_idx] * 1e-9 # convert to seconds
            ct = ct_full[:cut_idx]

            soln_bh, ct_bh, soln_mcmc, ct_mcmc = perform_different_fits(time, ct)
            
            cts_fit.append(ct_mcmc)
            full_data[label] = {'basin_hop':soln_bh,
                                     'mcmc':soln_mcmc, 
                                     }
            
            save_pickle(out_dir + str(i+1) + '.pkl', full_data)

        # 4) Append the data (containing fit parameters)
        data_all_blocks.append(full_data)
        
        # 5) Save the data so we can keep track of it. 
        save_pickle(out_dir + str(i+1) + '.pkl', full_data)
         
    return data_all_blocks


def main():
    
    
    NAME = sys.argv[1]
    
    print (f'We will work with: {NAME}')

    OUT_DIR_PREFIX = '/cluster/home/cchampion/work/NMR/ubiquitin/MF_all_methyls'

    # 1: Define paths for topology and trajectory
    root_dir = '/cluster/work/igc/mlehner/nmr_project/1ubq/'

    top_path = root_dir + '/ubq_desolv.gro'
    
    traj = f'{root_dir}/{NAME}/run_001/traj_4us_PBC_fit.xtc'
    
    # or for charmm
    #traj = f'{root_dir}/{NAME}/production/traj_4us_PBC_fit.xtc'
    print (traj) 
    out_dir = f'{OUT_DIR_PREFIX}/{NAME}_block'
    
    print ('quit before execution')
    exit() # comment this out to actually refit everything

    data_all_blocks = blockaveraged_ModelFree(top_path, 
                                              traj, 
                                              block_size = 1000, 
                                              out_dir = out_dir, 
                                              label_ff = NAME,
                                              cutoff_time = 50
                                             )

if __name__ == '__main__':
    main()

