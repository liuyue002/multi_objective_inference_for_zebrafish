'''
    This file is adapted from the code by Cleveland et al.
    This script is meant to be called from within MATLAB, 
    specifically by params_to_tda.m,
    and won't run by itself in Python.
    The following variables need to be passed in from MATLAB:
    D: a matrix of pairwise distances between cells
    dim: dimension to compute persistent homology (usually 0 or 1)
    p_max: maximum persistence of a barcode to consider
    pland_step: number of grid points to compute persistence landscape
'''
import sys

# Note: if you get errors due to Python not finding installed packages,
# try manually add the packages's install location to the search path like so:
#sys.path.extend(['', '/apps/external/conda/2025.02/lib/python312.zip', '/apps/external/conda/2025.02/lib/python3.12', '/apps/external/conda/2025.02/lib/python3.12/lib-dynload', '/home/liu4194/.local/lib/python3.12/site-packages', '/apps/external/conda/2025.02/lib/python3.12/site-packages'])

import numpy as np
from ripser import ripser
from persim import PersLandscapeApprox

dim=int(dim)
pland_step=int(pland_step)
if D.size == 0:
    barcodes = np.zeros((0,2))
    pland=np.zeros((0,pland_step))
else:
    bars = ripser(D,  distance_matrix=True, maxdim=dim,thresh=650)
    barcodes = bars['dgms'][dim]
    pla = PersLandscapeApprox(dgms=bars['dgms'],hom_deg=dim,start=0, stop=p_max, num_steps=pland_step)
    pland=pla.values
    if (pland.size == 1) and (pland[0] == 'empty'):
        pland=np.zeros((0,pland_step))
