"""This function clusters spots togheter to get factories.

Starting from the segmented spots, we cluster them giving
the same tag depending on their distance.
"""


import numpy as np
from skimage.measure import label, regionprops_table
from scipy.ndimage import binary_dilation
from cv2 import dilate
# import time


class ClusterInFactoriesEven:
    """Cluster spots, very fast, but works with even distances (0, 2, 4, 6, ...)"""
    def __init__(self, spts_lbls, dist_thr):

        spts_nosegm  =  label(np.sign(spts_lbls), connectivity=1).astype(np.int32)              # remove the segmentation: all the connected components are considered single objects
        spts_dil  =  np.zeros(spts_lbls.shape, dtype=np.int32)                                  # dilation to merge several connected components

        if dist_thr == 0:
            self.spts_fact  =  spts_nosegm
        else:

            for z in range(spts_lbls.shape[0]):
                spts_dil[z]  =  binary_dilation(np.sign(spts_nosegm[z]), iterations=dist_thr // 2)  # dilation happen frame by frame to avoid the growth in z (super fast)

            spts_dil_lbls  =  label(spts_dil, connectivity=1).astype(np.int32)                      # 3D label on the dilated spots
            spts_fact      =  spts_dil_lbls * np.sign(spts_nosegm)                                  # multiplication to make spots inherit labels from dilated matrix

            self.spts_fact  =  spts_fact



class ClusterInFactoriesOdd:
    """Cluster spots one by one, very slow but works for distance odd sitances too (1, 2, 3, ...)"""
    def __init__(self, spts_lbls, dist_thr):

        spts_nosegm       =  label(np.sign(spts_lbls))                                         # remove the segmentation: all the connected components are considered single objects

        spts_dil  =  np.zeros(spts_nosegm.shape, dtype=np.int32)                               # dilation to merge several connected components
        for z in range(spts_nosegm.shape[0]):
            spts_dil[z]  =  binary_dilation(np.sign(spts_nosegm[z]), iterations=dist_thr)      # dilation happen frame by frame to avoid the growth in z (super fast)

        rgp_dil  =  regionprops_table(label(spts_dil, connectivity=1), spts_nosegm, properties=["label", "intensity_image", "coords"])   # regionprops of dilated spots, using spts_nosegm as intensity image

        idxs2work  =  []
        for k in range(rgp_dil["label"].size):
            if np.unique(rgp_dil["intensity_image"][k]).size > 2:                           # if a label of dilated spots has more than one value in "fake" intensity, means that 2 spots are involved, we need to work on
                idxs2work.append(k)

        mtx2work  =  np.zeros(spts_dil.shape, dtype=bool)
        for kk in idxs2work:
            for kk_in in range(rgp_dil["coords"][kk].shape[0]):
                mtx2work[rgp_dil["coords"][kk][kk_in][0], rgp_dil["coords"][kk][kk_in][1], rgp_dil["coords"][kk][kk_in][2]]  =  1    # isolate the spots to work on building the 3D matrix by coordinates

        tags2work  =  np.unique(mtx2work * spts_nosegm)[1:]                                                                          # tags of spts_nosegm to work on
        rgp2work   =  regionprops_table(spts_nosegm, properties=["label", "coords"])                                                 # regionprops of these spts_nosegm
        idxs       =  np.where(np.in1d(rgp2work["label"], tags2work) == 1)[0]                                                        # indexes (no tags) to work on

        kernel        =  np.ones((3, 3), dtype=np.uint8)
        kernel[0, 0]  =  0
        kernel[0, 2]  =  0
        kernel[0, 2]  =  0
        kernel[2, 0]  =  0

        while idxs.size > 0:                                                                                                       # while the list of tags to check is still populated

            idx  =  idxs[0]
            print(idxs.size)
            spt  =  np.zeros(spts_nosegm.shape, dtype=np.int32)                                                                     # initilize single spot matrix

            spt[rgp2work["coords"][idx][:, 0], rgp2work["coords"][idx][:, 1], rgp2work["coords"][idx][:, 2]]  =  1                  # build the single spot by coordinates

            rgp_loc    =  regionprops_table(spt, properties=["coords"])
            zz_coords  =  np.unique(rgp_loc["coords"][0][:, 0])
            for zz in zz_coords:                                                                                                    # z-coordinates of the spot
                # spt[zz]  =  dilate(spt[zz].astype(np.uint8), kernel, iterations=dist_thr)                                           # dilation of the spot z by z
                spt[zz]  =  binary_dilation(spt[zz].astype(np.uint8), iterations=dist_thr)                                           # dilation of the spot z by z

                
            idxs_zeros  =  np.where(spt == 1)
            tags_loc    =  np.unique(spts_nosegm[idxs_zeros])[1:]

            print(tags_loc)

            if tags_loc.size > 1:                                                                                                   # if there is overlapping, you have 2 tags
                tags_loc  =  np.delete(tags_loc, np.argwhere(tags_loc == rgp2work["label"][idx]))                                   # remove the tag of the overlapping spot
                for tag_loc in tags_loc:
                    spts_nosegm[spts_nosegm == tag_loc]  =  rgp2work["label"][idx]                                                  # clustering: the overlapped spot takes the same tag of the overlapping one
                    idxs                                 =  np.delete(idxs, np.argwhere(rgp2work["label"] == tag_loc))              # remove the tag of the overlapped spot from the list of the tags to check
                idxs = np.delete(idxs, 0)                                                                                           # remove the tag ov the checked spot from the list of the tags to check

            else:
                idxs  = np.delete(idxs, 0)                                                                                          # if the dilated spot does not overlap just remove the tag from the list

        self.spts_fact  =  spts_nosegm




class ClusterInFactoriesOdd_Stunt:
    """Cluster spots one by one, very slow but works for distance odd sitances too (1, 2, 3, ...)"""
    def __init__(self, spts_lbls, dist_thr):

        kernel        =  np.ones((3, 3), dtype=np.uint8)
        kernel[0, 0]  =  0
        kernel[0, 2]  =  0
        kernel[2, 2]  =  0
        kernel[2, 0]  =  0

        spts_nosegm       =  label(np.sign(spts_lbls)).astype(np.int32)                                         # remove the segmentation: all the connected components are considered single objects

        spts_dil  =  np.zeros(spts_nosegm.shape, dtype=np.int32)                               # dilation to merge several connected components
        for z in range(spts_nosegm.shape[0]):
            spts_dil[z]  =  dilate(np.sign(spts_nosegm[z]).astype(np.uint8), kernel, iterations=dist_thr)      # dilation happen frame by frame to avoid the growth in z (super fast)

        rgp_dil  =  regionprops_table(label(spts_dil, connectivity=1).astype(np.int32), spts_nosegm, properties=["label", "intensity_image", "coords"])   # regionprops of dilated spots, using spts_nosegm as intensity image

        idxs2work  =  []
        for k in range(rgp_dil["label"].size):
            if np.unique(rgp_dil["intensity_image"][k]).size > 2:                           # if a label of dilated spots has more than one value in "fake" intensity, means that 2 spots are involved, we need to work on
                idxs2work.append(k)

        mtx2work  =  np.zeros(spts_dil.shape, dtype=bool)
        for kk in idxs2work:
            mtx2work[rgp_dil["coords"][kk][:, 0], rgp_dil["coords"][kk][:, 1], rgp_dil["coords"][kk][:, 2]]  =  1    # isolate the spots to work on building the 3D matrix by coordinates

        
        tags2work  =  np.unique(mtx2work * spts_nosegm)[1:]                                                                          # tags of spts_nosegm to work on
        print(tags2work)
        rgp2work   =  regionprops_table(spts_nosegm, properties=["label", "coords"])                                                 # regionprops of these spts_nosegm
        idxs       =  np.where(np.in1d(rgp2work["label"], tags2work) == 1)[0]                                                     # indexes (no tags) to work on

        while idxs.size > 0:                                                                                                       # while the list of tags to check is still populated

            idx  =  idxs[0]
            print(idxs.size)
            spt  =  np.zeros(spts_nosegm.shape, dtype=np.int32)                                                                     # initilize single spot matrix

            spt[rgp2work["coords"][idx][:, 0], rgp2work["coords"][idx][:, 1], rgp2work["coords"][idx][:, 2]]  =  1                  # build the single spot by coordinates

            rgp_loc  =  regionprops_table(spt, properties=["coords"])

            zz_coords  =  np.unique(rgp_loc["coords"][0][:, 0])
            for zz in zz_coords:                                                                        # z-coordinates of the spot
                spt[zz]  =  dilate(spt[zz].astype(np.uint8), kernel, iterations=dist_thr)                                                           # dilation of the spot z by z

                
            idxs_zeros  =  np.where(spt == 1)
            tags_loc    =  np.unique(spts_nosegm[idxs_zeros])[1:]

            print(tags_loc)

            if tags_loc.size > 1:                                                                                                   # if there is overlapping, you have 2 tags
                tags_loc  =  np.delete(tags_loc, np.argwhere(tags_loc == rgp2work["label"][idx]))                                   # remove the tag of the overlapping spot
                for tag_loc in tags_loc:
                    spts_nosegm[spts_nosegm == tag_loc]  =  rgp2work["label"][idx]                                                  # clustering: the overlapped spot takes the same tag of the overlapping one
                    idxs                                 =  np.delete(idxs, np.argwhere(rgp2work["label"] == tag_loc))              # remove the tag of the overlapped spot from the list of the tags to check
                idxs = np.delete(idxs, 0)                                                                                           # remove the tag ov the checked spot from the list of the tags to check

            else:
                idxs  = np.delete(idxs, 0)                                                                                          # if the dilated spot does not overlap just remove the tag from the list

        self.spts_fact  =  spts_nosegm
