"""This function desciminates TS fromm single molecules.

"""


import numpy as np
import kmeans1d
from skimage.measure import regionprops_table


class SearchTS:
    """Main class, does all the job"""
    def __init__(self, spts, raw_data_spts):

        rgp  =  regionprops_table(spts, raw_data_spts, properties=["label", "intensity_image", "area"])

        lbl_vol_ints_avints  =  []
        for k in range(len(rgp["label"])):
            lbl_vol_ints_avints.append([rgp["label"][k], rgp["area"][k], rgp["intensity_image"][k].sum(), rgp["intensity_image"][k].sum() / rgp["area"][k]])

        lbl_vol_ints_avints  =  np.asarray(lbl_vol_ints_avints)

        clusters_vol, ctrds     =  kmeans1d.cluster(lbl_vol_ints_avints[:, 1], 3)
        clusters_ints, ctrds    =  kmeans1d.cluster(lbl_vol_ints_avints[:, 2], 3)
        clusters_avints, ctrds  =  kmeans1d.cluster(lbl_vol_ints_avints[:, 3], 3)

        clusters_vol     =  np.asarray(clusters_vol)
        clusters_ints    =  np.asarray(clusters_ints)
        clusters_avints  =  np.asarray(clusters_avints)

        idxs  =  np.where(clusters_vol + clusters_ints + clusters_avints == 6)[0]

        ts_mtx  =  np.zeros(spts.shape, dtype=np.int32)
        for k in idxs:
            ts_mtx  +=  int(lbl_vol_ints_avints[k, 0]) * (spts == int(lbl_vol_ints_avints[k, 0]))

        self.ts_mtx  =  ts_mtx
