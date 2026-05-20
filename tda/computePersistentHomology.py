'''
Compute persistent homology and construct persistence landscape
given pairwise distance from a point-cloud dataset.

Copyright (C) Yue Liu, Alexandria Volkening, 2026
Adapted from Cleveland, Zhu, Sandstede, & Volkening (2023)

This program is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public 
License along with this program. If not, see 
https://www.gnu.org/licenses/.
'''


'''
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
