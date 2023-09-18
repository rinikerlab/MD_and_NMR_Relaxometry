# This file will contain all functions 
# performing the fits, to lighten the jupyter notebooks.
# so we can just import functions here

import numpy as np
import scipy.optimize as opt

from random import uniform

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.ticker import MaxNLocator

# My bounds
UPPER_BOUND_TAU_S = 1e-7
LOWER_BOUND_TAU_S = 1e-10

UPPER_BOUND_TAU_F = 0.5 * 1e-9
LOWER_BOUND_TAU_F = 1e-14

# Ferrage bounds:
#UPPER_BOUND_TAU_S = 1e-8
#LOWER_BOUND_TAU_S = 1e-10

#UPPER_BOUND_TAU_F = 1e-10
#LOWER_BOUND_TAU_F = 1e-14

# Fitting to a set of 6 exponentials as done by Hoffmann et al.
# Make the total S^2 not be a parameter (to ensure a sum of 1)

def decaying_exp1(t, A, tau_a):
    return (1-A) + A*np.exp(-t/tau_a)

def decaying_exp2(t, A, tau_a, B, tau_b):
    return (1-A-B) + A*np.exp(-t/tau_a) + B*np.exp(-t/tau_b)

def decaying_exp3(t, A, tau_a, B, tau_b, C, tau_c):
    return (1-A-B-C) + A*np.exp(-t/tau_a) + B*np.exp(-t/tau_b) + C*np.exp(-t/tau_c)

def decaying_exp4(t, A, tau_a, B, tau_b, C, tau_c, D, tau_d):
    cor_t = (1-A-B-C-D) + A*np.exp(-t/tau_a) + B*np.exp(-t/tau_b)
    cor_t += C*np.exp(-t/tau_c) + D*np.exp(-t/tau_d)
    return cor_t 

def decaying_exp5(t, A, tau_a, B, tau_b, C, tau_c, D, tau_d, E, tau_e):
    cor_t = (1-A-B-C-D-E) + A*np.exp(-t/tau_a) + B*np.exp(-t/tau_b)
    cor_t+= C*np.exp(-t/tau_c) + D*np.exp(-t/tau_d) + E*np.exp(-t/tau_e)
    return cor_t

def decaying_exp6(t, A, tau_a, B, tau_b, C, tau_c, D, tau_d, E, tau_e, F, tau_f):
    cor_t = (1-A-B-C-D-E-F) + A*np.exp(-t/tau_a) + B*np.exp(-t/tau_b) + C*np.exp(-t/tau_c) 
    cor_t += D*np.exp(-t/tau_d) + E*np.exp(-t/tau_e) + F*np.exp(-t/tau_f)
    return cor_t

#
# Functions to automatically generate guess parameters 
# for the linear-least squares fit
#

def gen_guess_params(num_exponentials, previous_guess=None):
    """
        This function will generate guess parameters 
        to use with the N exponential decay functions, 
        or with the extended model free function
    
        Guess performed based on results of a prior fit if such a fit
        was attempted. 

        This function needs to generate initial guesses so the timescales go from
        slower -> faster.

        First guess will be tau_slow

    """
    
    if previous_guess is not None:
        guess = np.append(previous_guess, previous_guess[-2:])
        guess[-1] = guess[-1] / 10 # not sure why this didn't work anymore...
    else:
        guess = np.tile([1, 1e-10], num_exponentials) # first guess will be 1 ns, had 100 before chosing to work in s.
        if num_exponentials == 1: return guess

        for i, g in enumerate(guess):
            if i % 2 == 1: 
                guess[i] = guess[i] / 10**(int(i/2))
    return guess

def gen_bounds(num_exponentials):
    """
        This function will generate  
        bounds, so all order parameters are in [0,1]
        and correlation times are positive [0, np.inf]
    """
    
    return (0, np.tile([1, np.inf], num_exponentials))
#
# Functions to check quality of the fit
#

def out_of_bounds(params):
    """
        At the moment this assumes paramaters given in the 
        same orders as in the 6 functions above.
    """
    if len(params) == 1:
        return False
    elif len(params) % 2 == 0:
        s = sum(params[0::2])
        print ('sum of values = ' + str(s))
        return (s>1)
    else:
        s = params[0]+sum(params[1::2]) # Add S2_long + A, B, C etc.
        print ('sum of values = ' + str(s))
        return (s>1)

def calc_chi2(y, y_fit, dy = []):
    """
    Calculate the goodness of fit based on the chi parameter
    
    y: data fit to
    
    y_fit: fit to that data
    
    d_y: standard deviations of y
    
    """
    if dy != []:
        return np.sum( (y-y_fit)**2.0/dy )/len(y)
    else:
        return np.sum( (y-y_fit)**2.0 )/len(y)

    
#
# Reproducing work from Ferrage 2018 et al.
# Rewriting the code myself based on: https://emcee.readthedocs.io/en/stable/tutorials/line/     
# 

def extended_model_free(t, Ss2, tau_s, Sf2, tau_f):
    """
    Extended model free type function, with a slow (Ss2, tau_s)
    and fast Sf2 and tau_f motion.
    """
    return (Ss2 * Sf2) + (1 - Sf2)*np.exp(-t/tau_f) + Sf2*(1-Ss2)*np.exp(-t/tau_s)

def pseudoChi2(params, time, ct_raw):
    """
    This function will calculate the chi2 (cost function for basin hopping) 
    based on the current parameters in the optimization procedure
    """
    Ss2 = params[0]
    tau_s = params[1]
    Sf2 = params[2]
    tau_f = params[3]
    return np.sum((ct_raw - extended_model_free(time, Ss2, tau_s, Sf2, tau_f))**2)


def fit_like_ferrage_2018(time, ct_raw):
    """
    Re-implementation of code from Fabien et al.
    
    Here we really call it with the same arguments to reproduce
    """
   
    niter = 1000
    initial_guess = gen_guess_params(2)
    
    # From Ferrage 2018
    #bounds = [[0.0, 1.0], [1e-10, 1e-8], # this second one is tau_s (in s)
    #          [0.0, 1.0], [1e-14, 1e-10]] # this last one is tau_f (in s)
    
    # Mine:
    bounds = [[0.0, 1.0], [LOWER_BOUND_TAU_S, UPPER_BOUND_TAU_S],
              [0.0, 1.0], [LOWER_BOUND_TAU_F, UPPER_BOUND_TAU_F]]

    minimizer_kwargs = dict(args = (time, ct_raw), method = "L-BFGS-B", bounds=bounds, jac=False)

    popt = opt.basinhopping(pseudoChi2, x0 = initial_guess, niter=niter,
                            minimizer_kwargs=minimizer_kwargs,
                            accept_test=None, disp=False)
    return popt

#
# Functions to perform and plot results of the MCMC fits. 
#

def intialize_mcmc_walkers(bounds, nWalkers):
    """
    This function makes the initial parameters for each walker of the MCMC run by 
    selecting randomly the intial position of every parameter (from uniform distrib in given bounds)
    
    """
    
    init_thetas = np.zeros([nWalkers, len(bounds)])
    for w in range(nWalkers):
        for p, bound in enumerate(bounds):
            init_thetas[w][p] = uniform(bound[0], bound[1])
            
    return init_thetas
        
def gen_bounds_mcmc(num_exponentials, popt = None):
    """
        This function will generate  
        bounds for the MCMC, so all order parameters are in [0,1]
        and correlation times are in the correct time ranges (same as Ferrage 2018)
        
        Note we have n+1 parameters as we need to add the log_f term (bounded by [0.000001, 1.0])
        
    """
    
    # Get bounds much closer to the actual results.
    if num_exponentials == 2:
        return [ [0, 1],[LOWER_BOUND_TAU_S, UPPER_BOUND_TAU_S], [0, 1], [LOWER_BOUND_TAU_F, UPPER_BOUND_TAU_F], [0.00001, 1.0]]
    return None

def log_likelihood(theta, x, y, yerr):
    Ss2, tau_s, Sf2, tau_f, log_f = theta
    model = extended_model_free(x, Ss2, tau_s, Sf2, tau_f)
    sigma2 = yerr ** 2 + model ** 2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))


def log_prior(theta): # here we will give in seconds for simplicity 
    Ss2, tau_s, Sf2, tau_f, log_f = theta    
    if 0.0 < Sf2 < 1.0 and 0.0 < Ss2 < 1.0 and LOWER_BOUND_TAU_F < tau_f < UPPER_BOUND_TAU_F and  LOWER_BOUND_TAU_S < tau_s < UPPER_BOUND_TAU_S  and -15 < log_f < 1.0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

#
# Plotting:
#

def plot_mcmc_parameter_search(samples, out_path, title):
    """ 
        Plots the convergence of the parameter search by the MCMC walkers.
    """    
    
    labels = [r"$S_s^2$", r"$\tau_s$", "$S_f^2$", r"$\tau_f$", r"log(f)"]
    num_params = len(samples[0][0])
    fig, axes = plt.subplots(num_params, figsize=(10, 7), sharex=True)
    
    for i in range(num_params):
        ax = axes[i]
        
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number");    
    
    plt.suptitle(title)

    plt.savefig(out_path)
    plt.close()

    return None    

def plot_mcmc_parameter_distribution(flat_samples, out_path, title):
    """ 
        Plots the parameter probability distrubtion functions based 
        on the MCMC search.

        flat_samples are obtained by calling the emcee function sampler.get_chain() with 
        appropriate arguments.        

    """

    labels = [r"$S_s^2$", r"$\tau_s$", "$S_f^2$", r"$\tau_f$", r"log(f)"]
    xlabel = [r"$S_s^2$", r"$\tau_s$ [s]", "$S_f^2$", r"$\tau_f$ [s]", r"log(f)"]
    xlims = [[0, 1], None, [0, 1], None]

    num_params = len(flat_samples[0]) 

    fig, axes = plt.subplots(num_params-1, figsize=(4, 10))

    name = "tab20"
    cmap = get_cmap(name)  # type: matplotlib.colors.ListedColormap
    colors = cmap.colors  # type: list

    for i, param in enumerate(flat_samples.T):
        if i == num_params-1: break # we don't show log(f)
        ax = axes[i] 
        if i == 0 or i == 2:
            text = 'mean: ' + f'{np.mean(param):.2f}'
        else:
            text = 'mean: ' + f'{np.mean(param):.2e}'
        
        bins = np.histogram_bin_edges(param, bins = 50)
        ax.hist(param, bins=bins, label = text, color = colors[i])
        ax.set_title('Probability distribution for: ' + labels[i])
        ax.set_xlabel(xlabel[i])
        ax.legend()
    
    plt.suptitle(title)

    fig.tight_layout()
    plt.savefig(out_path)
    plt.close()

    return None

