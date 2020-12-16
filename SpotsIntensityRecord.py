"""This function gives a vector of all the intensity of the detected spots.

Given the 3D matrices of raw data and detected spots, the output is a 7xN
matrix in which the components are the spot label, volume, intensity, max 
value and the 3 coordinates (z, x, y).
"""


import numpy as np
from skimage.measure import regionprops_table, regionprops
from PyQt5 import QtWidgets



class SpotsIntensityRecord:
    """Main class, does all the job"""
    def __init__(self, spts_segm, raw_spts):                                         # input_data are: spts_segm, raw_spts

        if spts_segm.max() > 0:
            rgp      =  regionprops_table(spts_segm.astype(np.int32), raw_spts, properties=["label", "area", "intensity_image", "max_intensity", "centroid"])
            rgp_cc   =  regionprops(spts_segm.astype(np.int32))                                                       # regionprops to measure the centroids coordinates (the prvious gives approximation of it)
            rgp_arr  =  []
            for k in range(len(rgp["label"])):
                rgp_arr.append([rgp_cc[k]["label"], int(np.round(rgp_cc[k]["centroid"][0])), int(np.round(rgp_cc[k]["centroid"][1])), int(np.round(rgp_cc[k]["centroid"][2]))])    # array with not approximate centroids-coordinates and label

            rgp_arr  =  np.asarray(rgp_arr)                                                                                                                                         # convert to array
        
            for k in range(rgp["label"].size):                                                                                                                                      # substitute the approximated centroids-coordinates with the non approximated ones
                idx2subst                                                         =  np.where(rgp_arr[:, 0] == rgp["label"][k])[0][0]                                               # check the indexes of the label
                rgp["centroid-0"][k], rgp["centroid-1"][k], rgp["centroid-2"][k]  =  rgp_arr[idx2subst][1], rgp_arr[idx2subst][2], rgp_arr[idx2subst][3]                            # substitute centroids-coordinates with the proper not approximated ones

            spts_ints_ctrs  =  np.zeros((7, rgp["area"].size), dtype=np.int32)

            pbar     =  ProgressBar(total1=len(rgp))
            pbar.show()
            pbar.update_progressbar(0)

            for k in range(rgp["area"].size):
                pbar.update_progressbar(k)
                spts_ints_ctrs[:, k]  =  np.array([rgp['label'][k], rgp['area'][k], rgp['intensity_image'][k].sum(), rgp['max_intensity'][k], rgp['centroid-0'][k], rgp['centroid-1'][k], rgp['centroid-2'][k]])

            pbar.close()

        else:
            spts_ints_ctrs  =  np.array([[], []])
            
        self.spts_ints_ctrs  =  spts_ints_ctrs



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
        """Method to update progressbar"""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()
