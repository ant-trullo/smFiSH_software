"""This function reads the data save by the analysis.

Each class takes care of a part of the analysis to load
it and respect the methods name to be sure you can rerun
algorithms in the main gui.
"""


import numpy as np
from skimage.color import label2rgb
from skimage.measure import label



class AnalysisLoaderNucsArea:
    """Load the nuclei dilation analysis"""
    def __init__(self, analysis_folder):

        nucs_dil  =  np.fromfile(analysis_folder + "/nucs_dil.bin", "uint16")
        nucs_dil  =  nucs_dil[2:].reshape((nucs_dil[1], nucs_dil[0]))

        self.nucs_dil  =  nucs_dil


class AnalysisLoaderNucsLabels:
    """Load nuclei segmentation analysis"""
    def __init__(self, analysis_folder, spts_mature):

        nucs_lbls  =  np.load(analysis_folder + "/nucs_segmented.npy")

        mycmap             =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3)) / 255.0
        self.nucs_lbls_3c  =  label2rgb(nucs_lbls, bg_label=0, bg_color=[0, 0, 0], colors=mycmap)

        spts  =  np.zeros(nucs_lbls.shape)
        for k in range(spts_mature.shape[1]):
            spts[spts_mature[3, k], spts_mature[4, k], spts_mature[5, k]]  =  1

        self.nucs_lbls_3c[:, :, :, 0]  *=  (1 - spts)
        self.nucs_lbls_3c[:, :, :, 1]  *=  (1 - spts)
        self.nucs_lbls_3c[:, :, :, 2]  *=  (1 - spts)
        self.nucs_lbls_3c[:, :, :, 1]  +=  spts

        self.nucs_lbls  =  nucs_lbls
        self.spts       =  spts



class AnalysisLoaderNucs:
    """Loads piled nuclei and fitting ellipsoids"""
    def __init__(self, analysis_folder):

        self.nucs_2d_det  =  np.load(analysis_folder + '/nucs_segm.npy')
        self.nucs_3d_det  =  np.load(analysis_folder + '/nucs_piled.npy')
        self.nucs_ellips  =  np.load(analysis_folder + '/nucs_ellips.npy')



class AnalysisLoaderSptsStudy:
    """Load the spots study"""
    def __init__(self, analysis_folder, a_b_flag):

        spts_g       =  np.load(analysis_folder + '/spts_g' + a_b_flag + '.npy')
        spts_others  =  np.load(analysis_folder + '/spts_other_vals' + a_b_flag + '.npy')

        self.spts_g    =  spts_g
        self.mu        =  spts_others[0]
        self.sigma     =  spts_others[1]
        self.spts_num  =  spts_others[2:]



class AnalysisLoaderSptsDetected:
    """"Load the detected (not segmented) spots"""
    def __init__(self, analysis_folder, a_b_flag):

        self.spts_detected  =  label(np.load(analysis_folder + '/spts_lbls' + a_b_flag + '.npy')).astype(np.int32)



class AnalysisLoaderSptsSegm:
    """Load the segmentation"""
    def __init__(self, analysis_folder, a_b_flag):

        self.spts_lbls  =  np.load(analysis_folder + '/spts_lbls' + a_b_flag + '.npy')
        if type(self.spts_lbls[1, 1, 1]) != np.int32:
            self.spts_lbls  =  self.spts_lbls.astype(np.int32)
            np.save(analysis_folder + '/spts_lbls' + a_b_flag + '.npy', self.spts_lbls)


class AnalysisLoaderParams:
    """Reads analysis parameters"""
    def __init__(self, analysis_folder, a_b_flag):

        file                   =  open(analysis_folder + '/params_text' + a_b_flag + '.txt', "r")
        a                      =  file.readlines()
        self.min_thr_value     =  int(a[0][20:-1])
        self.max_thr_value     =  int(a[1][20:-1])
        self.steps_thr_value   =  int(a[2][18:-1])
        self.sld_thr_current   =  int(a[3][25:])
        self.av_int_thr        =  int(a[4][30:])
        self.vol_thr           =  int(a[5][19:])
        self.min_numb_z_value  =  int(a[6][29:])
        self.sph_thr           =  int(a[7][23:])
        self.segm_clstr_flag   =  a[8][24:]
        self.thr_vals          =  np.linspace(self.min_thr_value, self.max_thr_value, self.steps_thr_value)
