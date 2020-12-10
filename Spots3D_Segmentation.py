"""This function filters the segmented spots a,d eventually corrects it.

It takes as input the 3D matrix of segmented spots and removes from the
small and the 2D spots. All the remaining are studied in terms of solidity:
if solidity is smaller than a certain threshold, the spot is considered as
several spots together which are than segmented with a clustering algorithm.
"""

import multiprocessing
import numpy as np
from skimage.measure import regionprops, label
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from PyQt5 import QtWidgets

import SearchBadSpotsUtilitySupport

# import SingleSpot3D_Segmentation



class Spots3D_Segmentation:
    """Main class, coordinates the action in the multiprocessing"""
    def __init__(self, spts, raw_spts, int_thr, vol_thr, min_numb_z):

        pbar    =  ProgressBar()
        pbar.show()
        pbar.update_progressbar(1)

        raw_spts   =  raw_spts.astype(np.float)
#        rgp        =  regionprops(spts, raw_spts)
        rgp        =  regionprops(spts)

        cpu_owe    =  multiprocessing.cpu_count()
        idxs_chop  =  np.array_split(np.arange(len(rgp)), cpu_owe)
        job_args   =  []
        for n in range(cpu_owe):
            job_args.append([spts.astype(np.int), rgp, idxs_chop[n], int_thr, vol_thr, min_numb_z])       # organize input arguments in lists for multiprocessing

        pool     =  multiprocessing.Pool()
        results  =  pool.map(SearchBadSpotsUtility, job_args)
        pool.close()

        pbar.update_progressbar(2)
        
        msk  =  results[0].msk
        for jj in range(1, len(results)):
            msk  +=  results[jj].msk

        spts_new  =  spts * (1 - msk * 1).astype('uint16')
        results   =  0
        del results

        pbar.update_progressbar(3)

        distance        =  ndi.distance_transform_edt(np.sign(spts_new))
#        local_maxi      =  peak_local_max(distance, indices=False, footprint=np.ones((25, 25, 25)), exclude_border=0)
        local_maxi      =  peak_local_max(distance, indices=False, footprint=np.ones((25, 25, 25)), exclude_border=0)
        local_maxi_lbl  =  label(local_maxi, connectivity=1)
        rgp_local       =  regionprops(local_maxi_lbl)
        new_local_maxi  =  np.zeros(local_maxi.shape)

        for k in range(len(rgp_local)):
            new_local_maxi[int(rgp_local[k]['centroid'][0]), int(rgp_local[k]['centroid'][1]), int(rgp_local[k]['centroid'][2])] = 1

        markers    =  ndi.label(new_local_maxi)[0]
        spts_lbls  =  watershed(- distance, markers, mask=np.sign(spts_new))

        pbar.update_progressbar(4)

        spts_left   =  label(np.sign(spts_new) - np.sign(spts_lbls), connectivity=1)
        spts_lbls   +=  (spts_lbls.max() + 1) * np.sign(spts_left) + spts_left

        # rgp_final  =  regionprops(spts_lbls)
        # spts_ctrs  =  np.zeros((3, len(rgp_final)), dtype=int)
        # for k in range(len(rgp_final)):
        #     spts_ctrs[:, k]  =  np.round(rgp_final[k]['centroid'][0]), np.round(rgp_final[k]['centroid'][1]), np.round(rgp_final[k]['centroid'][2])


        # self.spts_ctrs  =  spts_ctrs
        self.spts_lbls  =  spts_lbls



class SearchBadSpotsUtility:
    """Class to work in multiprocess; this selects the spots that will be analyzed"""
    def __init__(self, arg_vars):
 
        # msk          =  np.zeros(arg_vars[0].shape)
        list_tags    =  []
        list_vols    =  []
        list_coords  =  []
        for j in arg_vars[2]:
            list_tags.append(arg_vars[1][j]['label'])
            list_vols.append(arg_vars[1][j]['area'])
            list_coords.append(arg_vars[1][j]['coords'][:, 0])

        self.msk  =  SearchBadSpotsUtilitySupport.SearchBadSpotsUtilitySupport(arg_vars[0], list_tags, list_vols, list_coords, arg_vars[4], arg_vars[5])


#        for j in iddxx:
#            list_tags.append(rgp_sel[j]['label'])
#            list_vols.append(rgp_sel[j]['area'])
#            list_coords.append(rgp_sel[j]['coords'][:, 0])



class ProgressBar(QtWidgets.QWidget):
    """Simple progressbar widget"""
    def __init__(self, parent=None, total1=4):
        super(ProgressBar, self).__init__(parent)
        self.name_line1 = QtWidgets.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()



# class Spots3D_Segmentation:
#     def __init__(self, spts):

#         rgp        =  regionprops(spts)
#         cpu_ow     =  multiprocessing.cpu_count()
#         idxs_chop  =  np.array_split(np.arange(1, rgp[len(rgp) - 1]['label']), cpu_ow)
#         job_args   =  []
#         for n in range(cpu_ow):
#             job_args.append([spts, rgp, idxs_chop[n]])

#         pool     =  multiprocessing.Pool()
#         results  =  pool.map(SearchBadSpotsUtility, job_args)
#         pool.close()
#         print('Flag1')
#         msk  =  results[0].msk
#         for jj in range(1, len(results)):
#             msk  +=  results[jj].msk

#         spts_new  =  spts * (1 - msk * 1).astype('uint16')
#         results   =  0
#         del results
#         print('Flag2')
#         distance        =  ndi.distance_transform_edt(np.sign(spts_new))
#         local_maxi      =  peak_local_max(distance, indices=False, footprint=np.ones((10, 10, 10)), min_distance=6, exclude_border=1)
#         local_maxi_lbl  =  label(local_maxi)
#         rgp_local       =  regionprops(local_maxi_lbl)

#         new_local_maxi  =  np.zeros(local_maxi.shape)

#         for k in range(len(rgp_local)):
#             new_local_maxi[int(rgp_local[k]['centroid'][0]), int(rgp_local[k]['centroid'][1]), int(rgp_local[k]['centroid'][2])] = 1

#         markers    =  ndi.label(new_local_maxi)[0]
#         spts_lbls  =  watershed(-distance, markers, mask=np.sign(spts_new))
#         print('Flag3')
#         # rgp_final  =  regionprops(spts_lbls)
#         # spts_ctrs  =  np.zeros((3, len(rgp_final)), dtype=int)
#         # for k in range(len(rgp_final)):
#         #     spts_ctrs[:, k]  =  np.round(rgp_final[k]['centroid'][0]), np.round(rgp_final[k]['centroid'][1]), np.round(rgp_final[k]['centroid'][2])


#         # self.spts_ctrs  =  spts_ctrs
#         self.spts_lbls  =  spts_lbls



# class SearchBadSpotsUtility:
#     def __init__(self, arg_vars):

#             msk  =  np.zeros(arg_vars[0].shape)

#             for j in arg_vars[2]:                                           # intertia_tensor is a 3x3 matrix with a row=[0, 0, 0] if the spots is present in only one z or x or y plane
#                 if arg_vars[1][j]['area'] < 6 or len(np.where(arg_vars[1][j]['inertia_tensor'].sum(1) == 0)[0]) > 0:        
#                     msk  +=  (arg_vars[0] == arg_vars[1][j]['label'])

#             self.msk  =  msk





# spt_sing    =  (spts == 7305)
# rgp_sing    =  regionprops(spt_sing*1)
# coord_sing  =  rgp_sing[0]['coords']
# spt_zoom    =  spt_sing[coord_sing[:, 0].min()-5:coord_sing[:, 0].max()+5, coord_sing[:, 1].min()-5:coord_sing[:, 1].max()+5, coord_sing[:, 2].min()-5:coord_sing[:, 2].max()+5]

# rgp_zoom    =  regionprops(spt_zoom*1)
# coord_zoom  =  rgp_zoom[0]['coords']


# HINTS: try to work with a 3D medianfileter to identify the centroids ana than associate things in a random way (not complitely, introduce some fucking criteria) and go on up
# to the moment you are satisfied by the solidity value of each components you find

# d2               =  np.empty(ant_tets[:, :, :, 0].shape + (4,), dtype=np.ubyte)
# d2[:, :, :, :3]  =  ant_tets
# d2[..., 3]       =  d2[..., 0]*0.3 + d2[..., 1]*0.3
# d2[..., 3]       =  (d2[..., 3].astype(float) / 255.) **2 * 255


# from pyqtgraph.Qt import QtCore, QtGui
# import pyqtgraph.opengl as gl

# app  =  QtGui.QApplication([])
# w    =  gl.GLViewWidget()
# w.opts['distance'] = 20
# w.show()
# # w.setWindowTitle('pyqtgraph example: GLVolumeItem')

# # g = gl.GLGridItem()
# # g.scale(10, 1, 10)
# # w.addItem(g)

# v = gl.GLVolumeItem(d2)
# v.translate(50,50,100)
# w.addItem(v)
