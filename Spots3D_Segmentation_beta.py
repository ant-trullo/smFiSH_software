"""This function filters the segmented spots and eventually corrects them.

It takes as input the 3D matrix of detected spots and removes from the
small ones. All the remaining are studied in terms of sphericity:
if sphericity is smaller than a certain threshold, the spot is considered as
several spots together which are than segmented with a clustering algorithm.
"""


from importlib import reload
import multiprocessing
import numpy as np
from skimage.measure import regionprops
from skimage.morphology import label
from sklearn.cluster import KMeans
from PyQt5 import QtWidgets

import SearchBadSpotsUtilitySupport
import ClusterInFactories
import SearchTS



def inertia_sphericity(rr):                                                               # sphericity is calculated from the eigenvalues of inertia tensor: it is given by  
    """From regionprops results (intertia momentum eigenvalues) calculate a real number parameter"""
    return np.sqrt((rr[0] - rr[1]) ** 2 + (rr[0] - rr[2]) ** 2 + (rr[1] - rr[2]) ** 2)    # the square root of the sum of the square of the differencies of the 3 components



class Spots3D_Segmentation:
    """Main class, coordinates the action in the multiprocessing"""
    def __init__(self, spts, raw_data_spts, int_thr, vol_thr, min_numb_z, sph_thr, segm_clstr_flag):

        int_thr     =  np.int32(int_thr)                 # all these variable are converted in np.int32 to metch with the cython function
        vol_thr     =  np.int32(vol_thr)
        min_numb_z  =  np.int32(min_numb_z)
        sph_thr     =  np.int32(sph_thr)
                                
        pbar    =  ProgressBar()                          # progressbar widget
        pbar.show()
        pbar.update_progressbar(1)

        rgp        =  regionprops(spts)                                 # regionporops of the detected spots
        cpu_owe    =  multiprocessing.cpu_count()                       # number of cpu
        idxs_chop  =  np.array_split(np.arange(len(rgp)), cpu_owe)      # chop the indexes for multiprocessing
        job_args   =  []
        for n in range(cpu_owe):
            job_args.append([spts, rgp, idxs_chop[n], int_thr, vol_thr, min_numb_z])       # organize input arguments in lists for multiprocessing

        pool     =  multiprocessing.Pool()                              # pool
        # results  =  pool.map(Spots3D_Segmentation_beta.SearchBadSpotsUtility, job_args)           # find spots that don't fit with selected parameters (bad spots)
        results  =  pool.map(SearchBadSpotsUtility, job_args)           # find spots that don't fit with selected parameters (bad spots)
        pool.close()

        pbar.update_progressbar(2)

        msk  =  results[0].msk
        for jj in range(1, len(results)):
            msk  +=  results[jj].msk                                    # merge pool results
        spts_new  =  spts * (1 - msk * 1).astype('uint16')              # remove bad spots
        del results

        pbar.update_progressbar(3)

        if segm_clstr_flag == "s":
            ts_spots  =  SearchTS.SearchTS(spts_new, raw_data_spts).ts_mtx

            spts_lbls  =  label(spts_new, connectivity=1).astype(np.int32)                   # relabel good spots to have only consecutive tags for the spots - this is useful after
            rgp_new    =  regionprops(spts_new)                                              # regionprops of the spots to identify which one must be segmented

            spts_bad   =  np.zeros(spts_lbls.shape, dtype=np.int32)

            for k in range(len(rgp_new)):                                                   # for each spot calculate sphericity
                sph_bff  =  inertia_sphericity(rgp_new[k]['inertia_tensor_eigvals'])
                if sph_bff > sph_thr:                                                       # if it is higher than a certain threshold, the spots is added to the bad spots matrix
                    spts_bad  +=  (spts_new == rgp_new[k]['label'])

            spts_bad  *=  1 - np.sign(ts_spots)

            spts_lbls  =  label(spts_lbls * (1 - spts_bad), connectivity=1).astype(np.int32)    # bad spots are removed
            rgp_bad    =  regionprops(label(spts_bad, connectivity=1))                          # regionprops of the bad spots

            if len(rgp_bad) < 5 * cpu_owe:                                                                                                      # arbitrary parameter to chose between multiprocess or not
                spts_fixed  =  SegmCluster([rgp_bad, np.zeros(spts_lbls.shape, dtype=np.int32), range(len(rgp_bad)), sph_thr]).spts_lbls        # spots segmentation working on single processor
                # spts_fixed  =  Spots3D_Segmentation_beta.SegmCluster([rgp_bad, np.zeros(spts_lbls.shape, dtype=np.int), range(len(rgp_bad)), sph_thr]).spts_lbls        # spots segmentation working on single processor

                spts_lbls       +=  (spts_lbls.max() + 1) * np.sign(spts_fixed) + spts_fixed                                                   # adding the "ex-bad" spots to the total matrix with proper tags

                self.spts_lbls   =  spts_lbls

            else:
                job_args2  =  []                                                                                                                  # data input preparation for multiprocessing purposes
                jjjs       =  np.array_split(range(len(rgp_bad)), cpu_owe)
                for k in range(cpu_owe):
                    job_args2.append([rgp_bad, np.zeros(spts_lbls.shape, dtype=np.int32), jjjs[k], sph_thr])                                        # input data organized in a list
                    
                pool2     =  multiprocessing.Pool()                                                                                              # pool
                results2  =  pool2.map(SegmCluster, job_args2)                                                                                   # spots segmentation in multiprocessing
                # results2  =  pool2.map(Spots3D_Segmentation_beta.SegmCluster, job_args2)                                                                                   # spots segmentation in multiprocessing
                pool2.close()

                for k in range(cpu_owe):
                    spts_lbls  +=  (spts_lbls.max() + 1) * np.sign(results2[k].spts_lbls) + results2[k].spts_lbls                                # adding the "ex-bad" spots to the total matrix with proper tags

                self.spts_lbls  =  spts_lbls

        else:
            reload(ClusterInFactories)
            if sph_thr - 2 * (sph_thr // 2) == 0:
                self.spts_lbls  =  ClusterInFactories.ClusterInFactoriesEven(spts_new, sph_thr).spts_fact
            else:
                print("dlfkgjl")
                self.spts_lbls  =  ClusterInFactories.ClusterInFactoriesOdd_Stunt(np.copy(spts_new), sph_thr).spts_fact




class SearchBadSpotsUtility:
    """Class to work in multiprocess; this selects the spots that will be removed"""
    def __init__(self, arg_vars):
 
        list_tags    =  []                                          # initialization of lists of regionprops to feed a cython function
        list_vols    =  []
        list_coords  =  []
        for j in arg_vars[2]:
            list_tags.append(np.int32(arg_vars[1][j]['label']))                      # labels (tags) list
            list_vols.append(np.int32(arg_vars[1][j]['area']))                       # volume list
            list_coords.append(arg_vars[1][j]['coords'][:, 0].astype(np.int32))      # coordinate list of the pixels of each spot

        self.msk  =  SearchBadSpotsUtilitySupport.SearchBadSpotsUtilitySupport(arg_vars[0], list_tags, list_vols, list_coords, arg_vars[4], arg_vars[5])     # cython function to search bad spots



class SegmCluster:
    """Class to work in multiprocess; this segments spots"""
    def __init__(self, args_in):

        rgp_bad    =  args_in[0] 
        spts_lbls  =  args_in[1]
        lbls2work  =  args_in[2]
        sph_thr    =  args_in[3]

        for kk in lbls2work:
            tag_coords  =  rgp_bad[kk]['coords']                                            # coordinates of all the pixels ub the spot
            kmeans      =  KMeans(n_clusters=2, random_state=0).fit(tag_coords)             # clustering coordinates of pixels by distance: at this stage we cluster in two groups
            segm_2      =  np.zeros(spts_lbls.shape, dtype=np.int32)                        # initialization of the matrix with the spot segmented in two

            for ll in range(tag_coords.shape[0]):
                segm_2[tag_coords[ll][0], tag_coords[ll][1], tag_coords[ll][2]]  =  kmeans.labels_[ll] + 1      # assign the sub tag to the pixels

            rgp_s2  =  regionprops(segm_2)                                                                      # regionprops of the two segmented spots

            if inertia_sphericity(rgp_s2[0]['inertia_tensor_eigvals']) > sph_thr or inertia_sphericity(rgp_s2[1]['inertia_tensor_eigvals']) > sph_thr:     # check sphericity of the new sub spots
                kmeans  =  KMeans(n_clusters=3, random_state=0).fit(tag_coords)                                 # if not good, clustering in three groups
                segm_3  =  np.zeros(spts_lbls.shape, dtype=np.int32)

                for lll in range(tag_coords.shape[0]):
                    segm_3[tag_coords[lll][0], tag_coords[lll][1], tag_coords[lll][2]]  =  kmeans.labels_[lll] + 1   # assign the sub tag to the pixels

                spts_lbls  +=  spts_lbls.max() * np.sign(segm_3) + segm_3                                            # add the new spots (clustering in 3) to the final matrix
                
            else:
                spts_lbls  +=  spts_lbls.max() * np.sign(segm_2) + segm_2                                            # add the new spots (clustering in 3) to the final matrix if they gave good results (line 130)

        self.spts_lbls  =  spts_lbls



class ProgressBar(QtWidgets.QWidget):
    """Simple progressbar widget"""
    def __init__(self, parent=None, total1=4):
        super(ProgressBar, self).__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        """Class to update progressbar1"""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()
