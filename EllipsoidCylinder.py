"""This function expands ellipsoid.

The biggest ellips of the ellipsoid is copied in all the z
overlapping the center of the original ellips giving an
ellipsoid. Than the function takes an average between the
two solids using skeletonize algorithm.
"""

import numpy as np
from skimage.measure import regionprops
from skimage.morphology import skeletonize
from scipy.ndimage.morphology import binary_fill_holes


class EllipsoidCylinder:
    """Main function, does all the job"""
    def __init__(self, ellipsoid):

        xlen, ylen    =  ellipsoid.shape[1:]
        cylinder      =  np.zeros(ellipsoid.shape, dtype=np.int32)                # initialize cylinder matrix
        j             =  np.argmax(ellipsoid.sum(2).sum(1))                       # find the biggest ellipse (2D) of the ellipsoid
        msk_cyl       =  ellipsoid[j].astype(np.int32)                            # define the mask with the ellipse
        [i, j]        =  np.where(msk_cyl != 0)  
        msk_cyl       =  msk_cyl[i.min():i.max() + 1, j.min():j.max() + 1]        # crop the mask
        xsize, ysize  =  msk_cyl.shape
        frames2work   =  np.where(ellipsoid.sum(2).sum(1) != 0)[0]                # z-frames with the ellipsoid

        for k in frames2work:
            ctrs  =  np.round(regionprops(ellipsoid[k].astype(np.int))[0]['centroid']).astype(np.int)                                        # for each z measure the centrois position of the ellipsoid slide
            x0    =  ctrs[0] - xsize // 2
            x1    =  ctrs[0] - xsize // 2 + xsize
            y0    =  ctrs[1] - ysize // 2
            y1    =  ctrs[1] - ysize // 2 + ysize

            cylinder[k, np.max([0, x0]):np.min([x1, xlen]), np.max([0, y0]):np.min([y1, ylen])]  =  msk_cyl[np.max([0, -x0]):np.min([xlen, xlen - x0]), np.max([0, -y0]):np.min([ylen, ylen - y0])]     # add the mask in the proper location

        semi_cyl  =  np.zeros(cylinder.shape)                                   # initialize the buffer matrix
        for r in frames2work:
            crown     =  1 * ((cylinder[r] - ellipsoid[r]) == 1)                # define the crown, cylinder - ellipse
            semi_bff  =  binary_fill_holes(skeletonize(crown))
            if np.sum(semi_bff) > 500:
                semi_cyl[r]  =  semi_bff                                        # if the crown is ok, we take and skeleton and and it to have the average between ellipse of the ellipsoid and the one of the cylinder
            else:
                # print(r)
                semi_cyl[r]  =  ellipsoid[r]                                    # if the crown is not ok, we take the ellipse of the ellipsoide

        self.semi_cyl  =  semi_cyl
