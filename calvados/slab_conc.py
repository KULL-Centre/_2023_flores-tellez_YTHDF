import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from numba import jit
import string
from scipy.ndimage import gaussian_filter1d
from mdtraj import element
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from scipy.optimize import least_squares
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
from matplotlib.colors import LogNorm
import warnings
import itertools
warnings.filterwarnings('ignore')
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
from statsmodels.tsa.stattools import acf

import sys

def calcProfile(seq,m,T,L,value,error,tmin=1200,tmax=None,fbase='.',
    plot=False,pairs=False,X=None,blockingpath='/Users/sobuelow/software/BLOCKING',
    pden=3.,pdil=6.):
    sys.path.append(blockingpath)
    from main import BlockAnalysis
    if not pairs:
        X = m
    h = np.load(f'{fbase}/{m}/{T:d}/{X}_{T:d}.npy')
    # print(h.shape)
    N = len(seq)
    conv = 100/6.022/N/L/L*1e3
    h = h[tmin:tmax]*conv # corresponds to (binned) concentration in mM
    lz = h.shape[1]+1
    edges = np.arange(-lz/2.,lz/2.,1)/10
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    profile = lambda x,a,b,c,d : .5*(a+b)+.5*(b-a)*np.tanh((np.abs(x)-c)/d) # hyperbolic function, parameters correspond to csat etc.
    residuals = lambda params,*args : ( args[1] - profile(args[0], *params) )
    hm = np.mean(h,axis=0)
    z1 = z[z>0]
    h1 = hm[z>0]
    z2 = z[z<0]
    h2 = hm[z<0]
    p0=[1,1,1,1]
    res1 = least_squares(residuals, x0=p0, args=[z1, h1], bounds=([0]*4,[100]*4)) # fit to hyperbolic function
    res2 = least_squares(residuals, x0=p0, args=[z2, h2], bounds=([0]*4,[100]*4))

    cutoffs1 = [res1.x[2]-pden*res1.x[3],-res2.x[2]+pden*res2.x[3]] # position of interface - half width
    cutoffs2 = [res1.x[2]+pdil*res1.x[3],-res2.x[2]-pdil*res2.x[3]] # get far enough from interface for dilute phase calculation

    bool1 = np.logical_and(z<cutoffs1[0],z>cutoffs1[1])
    bool2 = np.logical_or(z>cutoffs2[0],z<cutoffs2[1])

    dilarray = np.apply_along_axis(lambda a: a[bool2].mean(), 1, h) # concentration in range [bool2]
    denarray = np.apply_along_axis(lambda a: a[bool1].mean(), 1, h)

    if (np.abs(cutoffs2[1]/cutoffs2[0]) > 2) or (np.abs(cutoffs2[1]/cutoffs2[0]) < 0.5): # ratio between right and left should be close to 1
        print('NOT CONVERGED',m,cutoffs1,cutoffs2)
        print(res1.x,res2.x)

    if plot:
        N = 10
        fig, ax = plt.subplots(1,2,figsize=(8,4))
        ax[0].plot(z,hm)
        for c1,c2 in zip(cutoffs1,cutoffs2):
            ax[0].axvline(c1,color='gray')
            ax[0].axvline(c2,color='black')
        ax[0].set(xlabel='z [nm]', ylabel='Concentration [mM]')
        ax[0].set(yscale='log')
        ax[0].set(title=m)
        xs = np.arange(len(h))
        dilrunavg = np.convolve(dilarray, np.ones(N)/N, mode='same')
        ax[1].plot(xs,dilarray)
        # ax[1].plot(xs,dilrunavg)
        ax[1].set(xlabel='Timestep',ylabel='cdil [mM]')
        fig.tight_layout()

    dil = hm[bool2].mean() # average concentration
    den = hm[bool1].mean()

    block_dil = BlockAnalysis(dilarray)
    block_den = BlockAnalysis(denarray)
    block_dil.SEM()
    block_den.SEM()

    # if pairs:
    #     mX = f'{m}_{X}'
    # else:
    # mX = m

    if pairs:
        value.loc[m,f'{X}_dil'] = block_dil.av 
        value.loc[m,f'{X}_den'] = block_den.av 

        error.loc[m,f'{X}_dil'] = block_dil.sem 
        error.loc[m,f'{X}_den'] = block_den.sem
    else:
        value.loc[m,'{:d}_dil'.format(T)] = block_dil.av 
        value.loc[m,'{:d}_den'.format(T)] = block_den.av 

        error.loc[m,'{:d}_dil'.format(T)] = block_dil.sem 
        error.loc[m,'{:d}_den'.format(T)] = block_den.sem

    # if pairs:
    #     value.loc[m, f'{X}_dil'] = block_dil.av
    #     error.loc[m, 'X'] = X
    # else:
    #     value.loc[m, 'X'] = m
    #     error.loc[m, 'X'] = m
    return(value, error)