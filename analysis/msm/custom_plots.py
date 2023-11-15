import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import pyemma

def plot_energized_data(dihedrals, energies, title = None):
    
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize =[7,6])
    
    chi1, chi2 = dihedrals.T[0], dihedrals.T[1]
    
    #colors = [cmap[int(i-1)] for i in clustered_dat]
    im = ax.scatter(dihedrals.T[0], dihedrals.T[1], c=energies, cmap='viridis', s =1)
    
    ax.set_ylim([0, 360])
    ax.set_xlim([0, 360])
    
    ax.set_ylabel('$\chi_2$')
    ax.set_xlabel('$\chi_1$')
    
    ax.legend(loc="lower left", ncol = 3)
    
    fig.colorbar(im, orientation='vertical')

    return fig, ax

def plot_clustered_data(dihedrals, clustered_dat, title = None):
    
    from matplotlib import cm
    active_qualitative_map_mligs = lambda n: plt.cm.viridis(np.linspace(0,1, n+1))

    n_clusters = max(clustered_dat)
    cmap = active_qualitative_map_mligs(n_clusters)
    
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize =[6, 6])
    
    chi1, chi2 = dihedrals.T[0], dihedrals.T[1]
    
    #colors = [cmap[int(i-1)] for i in clustered_dat]
    #ax.scatter(dihedrals.T[0], dihedrals.T[1], color = colors, s =1)
    
    for i in range(n_clusters+1):
        idx =np.where(clustered_dat == i)
    
        sub_chi1 = [chi1[j] for j in idx]
        sub_chi2 = [chi2[j] for j in idx]
        
        ax.scatter(sub_chi1, sub_chi2, color = cmap[i], s =1, label = f'cluster {i+1}')
        
    ax.set_ylim([0, 360])
    ax.set_xlim([0, 360])
    
    ax.set_ylabel('$\chi_2$')
    ax.set_xlabel('$\chi_1$')
    
    ax.legend(loc="lower left", ncol = 3)
    
    return fig, ax

def create_custom_cmap(base_color, reverse=True):
    """
    Will create a default colormap from white to base_color
    
    """
    tmp_name = 'test'
    if reverse:
        return LinearSegmentedColormap.from_list(tmp_name, [base_color, '#ffffff'], N=50)
    else:
        return LinearSegmentedColormap.from_list(tmp_name, ['#ffffff', base_color], N=50)

def nviridis(n):
    """
    return list of n colors extracted from viridis
    """
    return plt.cm.viridis(np.linspace(0,1,n))

def get_histogram(
        xall, yall, nbins=100,
        weights=None, avoid_zero_count=False, histrange=None):
    """Compute a two-dimensional histogram.

    Parameters
    ----------
    xall : ndarray(T)
        Sample x-coordinates.
    yall : ndarray(T)
        Sample y-coordinates.
    nbins : int, optional, default=100
        Number of histogram bins used in each dimension.
    weights : ndarray(T), optional, default=None
        Sample weights; by default all samples have the same weight.
    avoid_zero_count : bool, optional, default=True
        Avoid zero counts by lifting all histogram elements to the
        minimum value before computing the free energy. If False,
        zero histogram counts would yield infinity in the free energy.

    Returns
    -------
    x : ndarray(nbins, nbins)
        The bins' x-coordinates in meshgrid format.
    y : ndarray(nbins, nbins)
        The bins' y-coordinates in meshgrid format.
    z : ndarray(nbins, nbins)
        Histogram counts in meshgrid format.

    """
    z, xedge, yedge = np.histogram2d(xall, yall, bins=nbins, weights=weights, range=histrange)
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    if avoid_zero_count:
        z = np.maximum(z, np.min(z[z.nonzero()]))
    return x, y, z.T # transpose to match x/y-directions

def _to_density(z, specific_dens=None):
    """Normalize histogram counts.

    Parameters
    ----------
    z : ndarray(T)
        Histogram counts.

    """
    if specific_dens is not None:
        return z / float(specific_dens.sum())    
    else:
        return z / float(z.sum())


def _to_free_energy(z, minener_zero=False, specific_dens=None):
    """Compute free energies from histogram counts.

    Parameters
    ----------
    z : ndarray(T)
        Histogram counts.
    minener_zero : boolean, optional, default=False
        Shifts the energy minimum to zero.

    Returns
    -------
    free_energy : ndarray(T)
        The free energy values in units of kT.

    """
    pi = _to_density(z, specific_dens=specific_dens)
    free_energy = np.inf * np.ones(shape=z.shape)
    nonzero = pi.nonzero()
    free_energy[nonzero] = -np.log(pi[nonzero])

    if minener_zero:
        minim = np.min(free_energy[nonzero])
        free_energy[nonzero] -= minim
        return free_energy, minim
    return free_energy

def get_offsets(dihedrals, mtraj, n_states, specific_dens):
    """
    Here we align w.r.t. the entire density, to ensure the minima of the free energy are respected
    """ 
    
    mins = np.zeros(n_states)
    
    for i in range(n_states):
        idx = np.where(mtraj == i)
        _x, _y, _z = get_histogram(dihedrals[0][idx], dihedrals[1][idx], avoid_zero_count=False, histrange=[[0, 360], [0, 360]])
        _f, mi = _to_free_energy(_z, minener_zero=True, specific_dens=specific_dens) * 1
        
        mins[i] = mi
    
    mins -= np.min(mins)    
    return mins

def plot_msm_final(dihedrals, mtraj, alpha = 1, title = None, colors=None, with_leg=True):
   
    import matplotlib.ticker as ticker

    fontsize = 22

    fig, ax = plt.subplots(figsize=(15, 9))

    n_states = np.max(mtraj) +1
   
    if colors is None:
            colors = nviridis(n_states)
   
    cmaps = [create_custom_cmap(c) for c in colors]

    ax.set_ylim([0, 360])
    ax.set_xlim([0, 360])

    levels = np.arange(11)

    # To get the scale includng all 3 of them with proper location of the minima:
    x, y, z = get_histogram(dihedrals[0]+1000, dihedrals[1]+1000, avoid_zero_count=False, histrange=[[1000, 1360], [1000, 1360]])
    f, minim = _to_free_energy(z, minener_zero=True) * 1
    cs = ax.contourf(x, y, f, cmap='Greys_r', alpha=1, levels = levels)
   
    # To get coloring right we need to do this step twice:
    offsets = get_offsets(dihedrals, mtraj, n_states, z)
    
    contours = []

    # Now plot the actual data
    for i in range(n_states):
        idx = np.where(mtraj == i)
        _x, _y, _z = get_histogram(dihedrals[0][idx], dihedrals[1][idx], avoid_zero_count=False, histrange=[[0, 360], [0, 360]])
        _f, _ = _to_free_energy(_z, minener_zero=True) * 1

        _f += offsets[i] # gets the right colors for everything (minima is 0 only for lowest macrostate)

        percent_pop = np.rint(np.count_nonzero(mtraj == i) / len(mtraj) * 100)
        
        ax.contourf(_x, _y, _f, cmap=cmaps[i], alpha=1, levels = levels)
        contour = ax.contour(_x, _y, _f, colors='black', alpha=1, levels = levels, linewidths=0.5) 
         
        # To try circling it 
        # contour = ax.contour(_x, _y, _f, colors=[colors[i]], alpha=1, levels = [8, 10], linewidths=4)
        #vec = contour.allsegs[7][-1]
        #ax.text(np.max(vec.T[0]), np.max(vec.T[1]), s=f'{i+1}', fontsize = 25)

        if percent_pop == 0:
            ax.scatter(-100, -100, marker = 's', label = f'$S_{i+1}$: < 1%', color = colors[i])
        else:
            ax.scatter(-100, -100, marker = 's', label = f'$S_{i+1}$: {percent_pop:.0f}%', color = colors[i])
   
    ax.set_box_aspect(1)
   
    cbar = fig.colorbar(cs, shrink = 0.85)
    cbar.ax.set_ylabel('Free Energy / kT', fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    ax.set_xlabel(r'$\chi{}_1$ [degrees]', fontsize = fontsize)
    ax.set_ylabel(r'$\chi{}_2$ [degrees]', fontsize = fontsize)

    if title is not None:
        ax.set_title(title, fontsize = fontsize)

    ax.set_xticks(np.arange(60, 361, 120))
    ax.set_yticks(np.arange(60, 361, 120))
   
    ax.tick_params(axis="both",direction="in")
    ax.set_xticklabels(np.arange(60, 361, 120) , fontsize=fontsize)
    ax.set_yticklabels(np.arange(60, 361, 120) , fontsize=fontsize)
   
    if with_leg:
        #ax.legend(loc='lower center', fontsize = 0.6* fontsize, ncols = n_states, edgecolor='black')
        ax.legend(loc='center left', fontsize = fontsize, ncols = 1 , edgecolor='black', bbox_to_anchor=(-0.55, 0.5))

    return fig, ax, contours

def plot_msm_final2(dihedrals, mtraj, alpha = 1, title = None, colors=None, with_leg=True):
    
    from matplotlib.colors import LinearSegmentedColormap
    
    def add_labels(ax, fontsize):
        ax.set_box_aspect(1)
        ax.set_xlabel(r'$\chi{}_1$ [degrees]', fontsize = fontsize)
        ax.set_ylabel(r'$\chi{}_2$ [degrees]', fontsize = fontsize)

        ax.set_xticks(np.arange(60, 361, 120))
        ax.set_yticks(np.arange(60, 361, 120))

        ax.tick_params(axis="both",direction="in", width = 2, length = 10)
        ax.set_xticklabels(np.arange(60, 361, 120) , fontsize=fontsize)
        ax.set_yticklabels(np.arange(60, 361, 120) , fontsize=fontsize)
        
        # also set limits
        ax.set_ylim([0, 360])
        ax.set_xlim([0, 360])
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(2)

        
        
    import matplotlib.ticker as ticker

    
    def find_color(x, y):
        """
        finds an appropriate color from the tab20b
        """
        colors = [ 
                '#cedb9c', # light green
                '#e7cb94', # light orange   # FDD5B0   # old one: F99B45
                '#9c9ede', # light blue
                '#637939', # dark green
                '#bd9e39', # dark orange
                '#393b79', # dark blue
                '#7b4173', # dark purple 
                '#de9ed6', # light purple
                '#d6616b', # pink
                ]
        
        
        idx_x = (np.argmin(np.abs([60, 180, 360] - np.average(x))))
        idx_y = (np.argmin(np.abs([60, 180, 360] - np.average(y))))
        
        code = 10*idx_y + idx_x
        mapping = [0, 1, 2, 10, 11, 12, 20, 21, 22]
        return colors[mapping.index(code)]
        
    fontsize = 28
    
    fig = plt.figure(layout="constrained", figsize=(24, 9))
    mosaic = """AB"""
    ax_dict = fig.subplot_mosaic(mosaic)

    ax = ax_dict['A']
    n_states = np.max(mtraj) +1

    levels = np.arange(11)

    # To get the scale includng all 3 of them with proper location of the minima:
    x, y, z = get_histogram(dihedrals[0], dihedrals[1], avoid_zero_count=False, histrange=[[0, 360], [0, 360]])
    f, minim = _to_free_energy(z, minener_zero=True) * 1
    #cs = ax.contourf(x, y, f, cmap='Greys_r', alpha=1, levels = levels)
    cs = ax.contourf(x, y, f, cmap=create_custom_cmap('#141414'), alpha=1, levels = levels)
    
    contour = ax.contour(x, y, f, colors='black', alpha=1, levels = levels, linewidths=0.5)
    
    
    # To get coloring right we need to do this step twice:
    offsets = get_offsets(dihedrals, mtraj, n_states, z)
    
    add_labels(ax, fontsize)
    
    cbar = fig.colorbar(cs, shrink = 0.85)
    cbar.ax.set_ylabel('Free Energy / kT', fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    
    contours = []
    
    ax = ax_dict['B']

    # Now plot the actual data
    for i in range(n_states):
        idx = np.where(mtraj == i)
        percent_pop = np.rint(np.count_nonzero(mtraj == i) / len(mtraj) * 100)
        c = find_color(dihedrals[0][idx], dihedrals[1][idx])
        
        # BEFORE WE HAD
        _x, _y, _z = get_histogram(dihedrals[0][idx], dihedrals[1][idx], avoid_zero_count=False, histrange=[[0, 360], [0, 360]])
        _f, _ = _to_free_energy(_z, minener_zero=True) * 1

        _f += offsets[i] # gets the right colors for everything (minima is 0 only for lowest macrostate)
        ax.contourf(_x, _y, _f, cmap= LinearSegmentedColormap.from_list('test', [c, c], N=50), 
                    alpha=1, levels = levels)
        contour = ax.contour(_x, _y, _f, colors=c, alpha=1, levels = levels, linewidths=2)
        
        # simple scatter plot
        #ax.scatter(dihedrals[0][idx], dihedrals[1][idx], color=c, s = 10)
        
        #break
        
        if percent_pop == 0:
            ax.scatter(-100, -100, marker = 's', label = f'$S_{i+1}$: < 1%', color = c, s = 144)
        else:
            ax.scatter(-100, -100, marker = 's', label = f'$S_{i+1}$: {percent_pop:.0f}%', color = c, s = 144)
    
    add_labels(ax, fontsize)
    

    
    if title is not None:
        fig.suptitle(title, fontsize = fontsize*1.25)
    


    if with_leg:
        #ax.legend(loc='lower center', fontsize = 0.6* fontsize, ncols = n_states, edgecolor='black')
        legend = ax.legend(loc='center left', fontsize = fontsize, ncols = 1 , edgecolor='black', bbox_to_anchor=(1, 0.5))
        frame = legend.get_frame()
        frame.set_linewidth(2)

    return fig, ax, contours


def draw_arrow(ax, pos_1, pos_2, label="", width=1.0, arrow_curvature=1.0, color="black",
                patchA=None, patchB=None, shrinkA=2, shrinkB=1, arrow_label_size=None, arrow_label_location=.55):
    r""" Draws a slightly curved arrow from (x1,y1) to (x2,y2). 
         Will allow the given patches at start and end.
         copied from deeptime ! 
     """
    from matplotlib import patches
    
    # set arrow properties
    dist = np.linalg.norm(pos_2 - pos_1)
    arrow_curvature *= 0.075  # standard scale
    rad = arrow_curvature / dist
    tail_width = width
    head_width = max(2., 2 * width)

    if width < 1:
        color = 'lightgrey'

    arrow_style = patches.ArrowStyle.Simple(head_length=head_width, head_width=head_width, tail_width=tail_width)
    connection_style = patches.ConnectionStyle.Arc3(rad=-rad)
    arr = patches.FancyArrowPatch(posA=pos_1, posB=pos_2, arrowstyle=arrow_style,
                                  connectionstyle=connection_style, color=color,
                                  shrinkA=shrinkA, shrinkB=shrinkB, patchA=patchA, patchB=patchB,
                                  transform=ax.transData, 
                                    )
    ax.add_patch(arr)
    
    # Bezier control point
    control_vertex = np.array(arr.get_connectionstyle().connect(pos_1, pos_2).vertices[1])
    # quadratic Bezier at slightly shifted midpoint t = arrow_label_location
    t = arrow_label_location  # shorthand
    ptext = (1 - t) ** 2 * pos_1 + 2 * (1 - t) * t * control_vertex + t ** 2 * pos_2

    ax.text(*ptext, label, size=arrow_label_size, horizontalalignment='center', verticalalignment='center',
            zorder=3, transform=ax.transData)


#
#
# Helper methods to plot to MSM data
#
#

def plot_eigenmodes(msm, dihedrals, num=4):

    eigvec = msm.eigenvectors_right()
    print('The first eigenvector is one: {} (min={}, max={})'.format(
        np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))
    
    dtraj = msm.dtrajs_full[0] # works for us because we hacve a single traj
    
    
    ncol = 2 if num >= 2 else 1 
    nrow = int(np.rint(num/ncol))
    
    figsize = [4.5*ncol, 4.5*nrow]
    
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
    for i, ax in enumerate(axes.flat):
        if i == num:
            break
        pyemma.plots.plot_contour(
            *dihedrals.T,
            eigvec[dtraj, i + 1],
            ax=ax,
            cmap='PiYG',
            cbar_label='{}. right eigenvector'.format(i + 2),
            mask=True)
        ax.set_xlabel('chi 1')
        ax.set_ylabel('chi 2')
        ax.set_ylim([0, 360])
        ax.set_xlim([0, 360])
    fig.tight_layout()

    return fig, axes

def plot_its(msm, nits=5):
    def its_separation_err(ts, ts_err):
        """
        Error propagation from ITS standard deviation to timescale separation.
        """
        return ts[:-1] / ts[1:] * np.sqrt(
            (ts_err[:-1] / ts[:-1])**2 + (ts_err[1:] / ts[1:])**2)

    timescales_mean = msm.sample_mean('timescales', k=nits)
    timescales_std = msm.sample_std('timescales', k=nits)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].errorbar(
        range(1, nits + 1),
        timescales_mean,
        yerr=timescales_std,
        fmt='.', markersize=10)

    axes[1].errorbar(
        range(1, nits),
        timescales_mean[:-1] / timescales_mean[1:],
        yerr=its_separation_err(
            timescales_mean,
            timescales_std),
        fmt='.',
        markersize=10,
        color='C0')

    for i, ax in enumerate(axes):
        ax.set_xticks(range(1, nits + 1))
        ax.grid(True, axis='x', linestyle=':')

    axes[0].axhline(msm.lag * 0.01, lw=1.5, color='k')
    axes[0].axhspan(0, msm.lag * 0.01, alpha=0.3, color='k')
    axes[0].set_xlabel('implied timescale index')
    axes[0].set_ylabel('implied timescales / ns')
    axes[1].set_xticks(range(1, nits))
    #axes[1].set_xticklabels(
    #    ["{:d}/{:d}".format(k, k + 1) for k in range(1, nits + 2)],
    #    rotation=45)
    axes[1].set_xlabel('implied timescale indices')
    axes[1].set_ylabel('timescale separation')
    fig.tight_layout()

    return fig, axes

def plot_metastable(msm, dihedrals):
    
    dtraj = msm.dtrajs_full[0]

    fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flat):
        pyemma.plots.plot_contour(
            *dihedrals.T,
            msm.metastable_distributions[i][dtraj],
            ax=ax,
            cmap='viridis',
            mask=True,
            cbar_label='metastable distribution {}'.format(i + 1))
        ax.set_xlabel(' chi 1')
        ax.set_ylabel(' chi 2')
    fig.tight_layout()
    
    return fig, axes

def plot_pcca(msm, dihedrals, nstates):
    
    dtraj = msm.dtrajs_full[0]
    metastable_traj = msm.metastable_assignments[dtraj]
    
    fig, ax = plt.subplots(figsize=(5, 4))
    _, _, misc = pyemma.plots.plot_state_map(
        *dihedrals.T, metastable_traj, ax=ax)

    ax.set_xlabel('chi 1')
    ax.set_ylabel('chi 2')

    misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
                                 for i in range(nstates)])
    fig.tight_layout()

    return fig, ax
