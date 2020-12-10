"""This function segments 3D nuclei z by z.

Raw data are 3D (z-stack) nuclei and we perform segmentation slice by slice.
"""

import numpy as np
from skimage.filters import median, gaussian, threshold_otsu
from skimage.morphology import remove_small_objects, label, disk   # , watershed
# from skimage.measure import regionprops
from scipy.ndimage.morphology import binary_fill_holes, binary_dilation, binary_erosion
# from scipy.ndimage import distance_transform_edt
from PyQt5 import QtWidgets


class NucsSegmenter:
    """Main class, does all the job"""
    def __init__(self, nucs):

        zlen    =  nucs.shape[0]
        nucs_f  =  np.zeros(nucs.shape)
        for z in range(zlen):
            nucs_f[z]  =  gaussian(nucs[z], 1.5)                                      # gaussian 2D filtering for each frame

        nucs_log   =  np.log(nucs_f)                                                  # logarithm is needed because of the non-homoheneity of the DAPI in the nucleus
        j          =  np.argmax(nucs_log.sum(2).sum(1))                               # otsu threshold based on the most intense frame: we need no detection in the empty frames
        otsu_thr   =  threshold_otsu(nucs_log[j])                                     # otsu threshold
        nucs_thr   =  nucs_log > otsu_thr                                             # thresholding

        pbar  =  ProgressBar(total1=zlen)
        pbar.show()
        pbar.update_progressbar1(0)

        nucs_lbls  =  np.zeros(nucs.shape, dtype=np.int)
        for z in range(zlen):
            pbar.update_progressbar1(z + 1)
            nucs_lbls[z]  =  binary_fill_holes(nucs_thr[z])                                # fill the holes in the nuclei shadow 
            nucs_lbls[z]  =  label(nucs_lbls[z])                                           # labels nuclei shadow
            nucs_lbls[z]  =  median(np.sign(nucs_lbls[z]).astype('uint8'), disk(5))        # median filter to compact the nuclei shadows after thresholding
            nucs_lbls[z]  =  binary_dilation(nucs_lbls[z], iterations=3)                   # fill the holes in the nuclei shadow 
            nucs_lbls[z]  =  binary_erosion(nucs_lbls[z], iterations=3)                    # fill the holes in the nuclei shadow 
#            nucs_lbls[z]  =  binary_fill_holes(nucs_lbls[z])                               # fill the holes in the nuclei shadow 
            nucs_lbls[z]  =  label(nucs_lbls[z])                                           # fill the holes in the nuclei shadow 
            nucs_lbls[z]  =  remove_small_objects(nucs_lbls[z], 250)                       # removes small connected components
            nucs_lbls[z]  =  remove_small_objects(nucs_lbls[z], 250)                       # removes small holes in connected components
            

        # rgps  =  []
        # for z in range(zlen):
        #     rgps.append(regionprops(nucs_lbls[z]))                      # regionprops to extract information about all the area of nuclei in each frame

        # aarr  =  []                                                     # list of all the areas
        # for z in range(zlen):
        #     for k in range(len(rgps[z])):        
        #         aarr.append(rgps[z][k]['area'])

        # aarr      =  np.asarray(aarr)
        # area_thr  =  aarr.mean() * 1.5                                  # area threshold to detect nuclei that must be segmented with watershed

        # nucs2correct  =  np.zeros(nucs.shape)
        # for z in range(zlen):
        #     for k in range(len(rgps[z])):
        #         if rgps[z][k]['area'] > area_thr:                                                               # selection with area threshold
        #             if rgps[z][k]['area']/rgps[z][k]['bbox_area'] < 0.6 or rgps[z][k]['eccentricity'] > 0.7:        # selection with eccentricity or ratio between area and bounding box area
        #                 nucs2correct[z]  +=  rgps[z][k]['label'] * (nucs_lbls[z] == rgps[z][k]['label']) 

        # frames2work  =  np.where(np.sum(nucs2correct, axis=(1,2)) != 0)[0]

        # pbar.progressbar2.setMaximum(frames2work.size)
        # pbar_idx2  =  0        

        # """To define the peaks for the watershed algorithm we work on the distance image:
        #     we perform recursively a thresholding on the distance image starting from a 
        #     very high value and decreasing it at each step. The number of connected components 
        #     varies at each step, in particular it increases at beginning, reaches a maximum and 
        #     than decreases when connected components start to touch each other. The selected value 
        #     for the thresholding is at the maximum number of components."""


        # for jj in frames2work:
        #     pbar.update_progressbar2(pbar_idx2 + 1)
        #     pbar_idx2     +=  1
        #     distance       =  distance_transform_edt(nucs2correct[jj])       # distance matrix
        #     dist_thr       =  distance.max()                                 # starting thresholding value
        #     num_obj        =  []                                                
        #     l              =  0
        #     distance_thrs  =  np.zeros((int(dist_thr) + 1,) + (distance.shape))   # array for the thresholded images
        
        #     while dist_thr > 0:
        #         mtx_bff  =  distance > dist_thr                     # thresholding
        #         mtx_bff  =  binary_erosion(mtx_bff)                 # erosion to avoid wrong results due to small touch
        #         mtx_bff  =  remove_small_objects(mtx_bff, 15)       # remove small objects to avoid wrong results
        #         num_obj.append(label(mtx_bff).max())                # number of objects
        #         distance_thrs[l]   =  mtx_bff                       # store the thresholded image
        #         dist_thr          -=  1
        #         l                 +=  1

        
        #     num_obj  =  np.asarray(num_obj)
        #     j_slctd  =  np.where(num_obj == num_obj.max())[0][-1]       # find maximum number of components

        #     rgp_bff   =  regionprops(label(distance_thrs[j_slctd]))     # detect centers
        #     peak_mtx  =  np.zeros(mtx_bff.shape)
        #     for ll in range(len(rgp_bff)):
        #         peak_mtx[rgp_bff[ll]['centroid'][0].astype(np.int), rgp_bff[ll]['centroid'][1].astype(np.int)]  =  1   # peak matrix

        #     markers         =  label(peak_mtx)
        #     new_nucs        =  watershed(-distance, markers, mask=np.sign(nucs2correct[jj]))     # watershed action
        #     nucs_lbls[jj]  *=  (1 - np.sign(new_nucs))
        #     nucs_lbls[jj]  +=  new_nucs + np.sign(new_nucs) * nucs_lbls[jj].max()                # substituing in the output matrix
            
        pbar.close()

        self.nucs_lbls      =  nucs_lbls
        # self.num_obj        =  num_obj
        # self.distance_thrs  =  distance_thrs



class ProgressBar(QtWidgets.QWidget):
    """Double progress bar widget"""
    def __init__(self, parent=None, total1=20, total2=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        self.progressbar2  =  QtWidgets.QProgressBar()
        self.progressbar2.setMinimum(1)
        self.progressbar2.setMaximum(total2)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)
        main_layout.addWidget(self.progressbar2, 1, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar1(self, val1):
        """First progress bar updater"""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()

    def update_progressbar2(self, val2):
        """Second progress bar updater"""
        self.progressbar2.setValue(val2)
        QtWidgets.qApp.processEvents()
