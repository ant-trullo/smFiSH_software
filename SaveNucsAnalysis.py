"""This function saves the nucli analysis.

It writes a '.npy' for each of the matrices defined.
"""

import numpy as np


class SaveNucsAnalysis:
    """Main class, does all the job"""
    def __init__(self, folder2write, nucs_2d_det, nucs_3d_det, nucs_ellips):
        
        np.save(folder2write + '/nucs_segm.npy', nucs_2d_det)
        np.save(folder2write + '/nucs_piled.npy', nucs_3d_det)
        np.save(folder2write + '/nucs_ellips.npy', nucs_ellips)
    
