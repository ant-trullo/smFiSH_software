"""This function splits overlapping ellipsoids.

Given 2 ellipsoids that overlapps, this function splits the
overlapping region between the 2. The splitting is gained 
dilating alternatively z by z the slides of both the ellpsoids.
Of course the dilation must happen inside the overlapping region
only.
"""

import numpy as np
from scipy.ndimage.morphology import binary_dilation



class SplitEllipsoids:
    """Main class, does all the job"""
    def __init__(self, joint_ellipsoids):

        frames2work  =  np.where(np.sum((joint_ellipsoids == 3), axis=(1,2)) != 0)[0]             # z frames with nuclei to work on
        ellipsoids   =  np.zeros(joint_ellipsoids.shape)                                          # initialization of the matrix with splitted ellipsoids
        
        for z in frames2work:                                                                     # all the operation are performed z by z
            mask_tot   =  np.sign(joint_ellipsoids[z]).astype(np.int)                             # total mask of the overlapping slides
#            mask_ovrl  =  (joint_ellipsoids[z] == 3) * 1                                         # mask of the overlapping
            mask1      =  (joint_ellipsoids[z] == 1) * 1                                          # mask of the first ellipsoid
            mask2      =  (joint_ellipsoids[z] == 2) * 1                                          # mask of the second ellipsoid

            while (mask_tot - mask1 - mask2).sum() != 0:                                          # check if the job is done (overlapping surfaces reduces to 0)
                mask1  =  binary_dilation(mask1) * (1 - mask2) * mask_tot                         # mask1 dilation is masked with the other mask to avoid overlapping and with the total mask to avoid that the ellipsoid grows out of the starting overlapping
                mask2  =  binary_dilation(mask2) * (1 - mask1) * mask_tot                         # same thing
            
            ellipsoids[z]  =  mask1 * 1 + mask2 * 2                                               # dilated ellipsoids slices are updated in the output matrix
            
        self.ellipsoids  =  ellipsoids + joint_ellipsoids * (1 - np.sign(ellipsoids))             # add to the final result the ellipdoids in the frames in which there where no overlapping
