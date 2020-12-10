"""This function detects spots (single molecules and TS) in 3D.

Inputs are raw data plus a threshold value.
"""



import numpy as np
from scipy.ndimage import filters
from scipy.stats import norm
from skimage.morphology import label
# from skimage.measure import regionprops
from PyQt5 import QtGui, QtWidgets
# import multiprocessing



class SpotsDetection3D_RecursiveDoG:
    """Detection is implemented as a difference of Gaussian filter and works recursively on a number of threshold:
        the output is the number of detected spots per thrshold value."""
    def __init__(self, green4d, thr_vals):

        steps  =  green4d.shape[0]
        pbar   =  ProgressBar(total1=thr_vals.size + steps)
        pbar.show()

        spts_g        =  np.zeros(green4d.shape)
        steps         =  green4d.shape[0]
        for t in range(steps):
            # print(t)
            pbar.update_progressbar(t)
            spts_g[t, :, :]        =  filters.gaussian_filter(green4d[t, :, :].astype(np.float), 1.2) - filters.gaussian_filter(green4d[t, :, :].astype(np.float), 2.2)

        (mu, sigma)  =  norm.fit(spts_g)                                          # histogram is fitted with a Gaussian function
        spts_num     =  np.zeros(thr_vals.size)


        for k in range(thr_vals.size):
            pbar.update_progressbar(k + steps)
            # print(k)
            spts_num[k]  =  label(spts_g > mu + thr_vals[k] * sigma, connectivity=1).max()

        pbar.close()

        self.spts_g    =  spts_g
        self.spts_num  =  spts_num
        self.mu        =  mu
        self.sigma     =  sigma



class SpotsDetection3D_SingleThr:
    """Same as previous class, but for a single threshold value"""
    def __init__(self, spts_gl, mu, sigma, thr):

        spts_gl_thr   =  spts_gl > mu + thr * sigma                           # thresholding on the histogram
        spts_gl_lbls  =  label(spts_gl_thr, connectivity=1)                   # labelling

        self.spts_detect  =  spts_gl_lbls.astype(np.int32)



class ProgressBar(QtGui.QWidget):
    """Simple progressbar widget"""
    def __init__(self, parent=None, total1=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1 = QtGui.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtGui.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()
