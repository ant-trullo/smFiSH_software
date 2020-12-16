"""This function selects spots whos centroid as a certain z-value.

Given as input the spots matrix and the z boundaries, the function
measures centroid positions and selects them accordingly with
boundaties.
"""


import numpy as np
from skimage.measure import regionprops_table, regionprops


class RemoveCentralSpots:
    """Remove central detected spots"""
    def __init__(self, spts_detected, first_last_frame):

        rgp      =  regionprops_table(spts_detected, properties=["label", "coords", "centroid"])     # properties of spots organized in a dictionary
        rgp_cc   =  regionprops(spts_detected)                                                       # regionprops to measure the centroids coordinates (the prvious gives approximation of it)
        rgp_arr  =  []
        for k in range(len(rgp["label"])):
            rgp_arr.append([rgp_cc[k]["label"], int(np.round(rgp_cc[k]["centroid"][0])), int(np.round(rgp_cc[k]["centroid"][1])), int(np.round(rgp_cc[k]["centroid"][2]))])    # array with not approximate centroids-coordinates and label

        rgp_arr  =  np.asarray(rgp_arr)                                                                                                                                         # convert to array
        
        for k in range(rgp["label"].size):                                                                                                                                      # substitute the approximated centroids-coordinates with the non approximated ones
            idx2subst                                                         =  np.where(rgp_arr[:, 0] == rgp["label"][k])[0][0]                                               # check the indexes of the label
            rgp["centroid-0"][k], rgp["centroid-1"][k], rgp["centroid-2"][k]  =  rgp_arr[idx2subst][1], rgp_arr[idx2subst][2], rgp_arr[idx2subst][3]                            # substitute centroids-coordinates with the proper not approximated ones

        idxs_dwn  =  np.where(rgp["centroid-0"] <= first_last_frame[0])                                                                                                         # select spots  whos centroid is below a certain z-threshold
        idxs_up   =  np.where(rgp["centroid-0"] >= first_last_frame[1])                                                                                                         # select spots whos centroid is above a certain z-threshold

        idxs2keep      =  np.append(idxs_dwn, idxs_up)                                                                                                                          # join the 2 array of indexes
        spts_selected  =  np.zeros(spts_detected.shape, dtype=np.int32)                                                                                                         # initialize the output matrix
        
        for kk in idxs2keep:
            for kk_in in range(rgp["coords"][kk].shape[0]):
                spts_selected[rgp["coords"][kk][kk_in][0], rgp["coords"][kk][kk_in][1], rgp["coords"][kk][kk_in][2]]  =  rgp["label"][kk]      # insert the spots in the output matrix by coordinates witht the proper tag

        self.spts_selected  =  spts_selected


class RemoveCentralSegmentedSpots:
    """Remove central segmented spots"""
    def __init__(self, spts_segmented, first_last_frame):

        rgp      =  regionprops_table(spts_segmented, properties=["label", "coords", "centroid"])     # properties of spots organized in a dictionary
        rgp_cc   =  regionprops(spts_segmented)                                                       # regionprops to measure the centroids coordinates (the prvious gives approximation of it)
        rgp_arr  =  []
        for k in range(len(rgp["label"])):
            rgp_arr.append([rgp_cc[k]["label"], int(np.round(rgp_cc[k]["centroid"][0])), int(np.round(rgp_cc[k]["centroid"][1])), int(np.round(rgp_cc[k]["centroid"][2]))])    # array with not approximate centroids-coordinates and label

        rgp_arr  =  np.asarray(rgp_arr)                                                                                                                                        # convert to array
        
        for k in range(rgp["label"].size):                                                                                                                                     # substitute the approximated centroids-coordinates with the non approximated ones
            idx2subst                                                         =  np.where(rgp_arr[:, 0] == rgp["label"][k])[0][0]                                              # check the indexes of the label
            rgp["centroid-0"][k], rgp["centroid-1"][k], rgp["centroid-2"][k]  =  rgp_arr[idx2subst][1], rgp_arr[idx2subst][2], rgp_arr[idx2subst][3]                           # substitute centroids-coordinates with the proper not approximated ones

        idxs_dwn  =  np.where(rgp["centroid-0"] <= first_last_frame[0])                                                                                                        # select spots  whos centroid is below a certain z-threshold
        idxs_up   =  np.where(rgp["centroid-0"] >= first_last_frame[1])                                                                                                        # select spots whos centroid is above a certain z-threshold

        idxs2rem  =  np.intersect1d(idxs_dwn, idxs_up)                                                                                                                         # join the 2 array of labels

        msk  =  np.zeros(spts_segmented.shape, dtype=np.int32)                                                                                                                 # initialize matrix to remove from the input
        for kk in idxs2rem:
            for kk_in in range(rgp["coords"][kk].shape[0]):
                msk[rgp["coords"][kk][kk_in][0], rgp["coords"][kk][kk_in][1], rgp["coords"][kk][kk_in][2]]  =  1      # insert the spots by coordinates: it is B&W

        self.spts_selected  =  spts_segmented * (1 - msk)


