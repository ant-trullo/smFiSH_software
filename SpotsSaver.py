"""This function saves the results of spots detection for smiFISh analysis."""

import numpy as np
from skimage.measure import regionprops


class SpotsSaver:
    def __init__(self, spts_segm):

        spts_ctrs  =  np.zeros((3, spts_segm.max() + 200))
        rgp        =  regionprops(spts_segm)
        
        for k in range(len(rgp)):
            spts_ctrs[:, k]  =  rgp[k]['centroid'][0], rgp[k]['centroid'][1], rgp[k]['centroid'][2]

        spts_ctrs  =  np.delete(spts_ctrs, np.where(spts_ctrs.sum(0) == 0)[0], axis=1)
