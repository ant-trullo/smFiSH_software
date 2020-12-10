"""This function fitts piled nuclei with a 3D ellipsoid.

Taking as input argument the z-stack of nuclei segmented
frame by frame and piled up, this function gives as output
the matrix of all the fitting ellipsoids (one per nucleus)
and 
"""


import numpy as np
from skimage.measure import regionprops
from scipy.spatial import ConvexHull
from scipy.ndimage.morphology import binary_erosion
from matplotlib.path import Path
from PyQt5 import QtWidgets

import EllipsoidCylinder
import EllipsoidTool
import SplitEllipsoids


class NucsEllipsoidFitting:
    """Main function, does all the job"""
    def __init__(self, nucs_piled, nucs_fitting_flag="Ellips"):

        mask4rmv              =  np.ones(nucs_piled[0].shape, dtype=np.int32)
        mask4rmv[4:-4, 4:-4]  =  0                                                  # border mask, to detect nuclei touching the border

        idxs2rmv  =  np.unique(mask4rmv * nucs_piled)[1:]                           # tags of nuclei touching the border
        for vv in idxs2rmv:
            nucs_piled  *=  (1 - (nucs_piled == vv))                                # removing nuclei touching the border

        zlen, xlen, ylen  =  nucs_piled.shape
        nucs_elips        =  np.zeros(nucs_piled.shape)                             # matrix with all the fitting ellipsoids

        idxs  =  np.unique(nucs_piled)[1:]
        pbar  =  ProgressBar(total1=idxs.size)
        pbar.show()
        pbar_idx  =  1

        a  =  EllipsoidTool.EllipsoidTool()                                  # initialize the fitting tool 
        for k in idxs:
            # pbar.update_progressbar(pbar_idx)
            pbar_idx  +=  1
            print(k)
            single     =  (nucs_piled == k)                                         # isolate the matrix of a single nucleus (z stack)
            single_er  =  binary_erosion(single, iterations=1)                      # erosion to take only the external surface
            mtx        =  single * 1  -  single_er * 1                              # external 3D surface (to speed up the fitting)
            rgp        =  regionprops(mtx)                                          # regionprops of the external surface to get the coordintes of all the points of the surface

            params   =  a.getMinVolEllipse(rgp[0]['coords'])                           # calculate fitting parameters
            brd_pts  =  np.asarray(a.plotEllipsoid(params[0], params[1], params[2]))   # array of coordinates points of the fitting ellipsoid   

            brd_pts                            =  np.abs(np.round(brd_pts).astype(np.int))    # approximate the coordinates to integers (they will be used as pixel coordinates)
            brd_pts[brd_pts[:, 0] >= zlen, 0]  =  zlen - 1
            brd_pts[brd_pts[:, 1] >= xlen, 1]  =  xlen - 1
            brd_pts[brd_pts[:, 2] >= ylen, 2]  =  ylen - 1


            """The result of the fitting tool is a set of coordinate of the fitting ellipsoid.
            These coordinates do not give a solid surface, but a fuzzy surface (don't know why),
            so we need to compact the points to get the final solid ellipsoid. This job is done
            using convex Hull surfaces.
            """

            brd_mtx  =  np.zeros(single.shape, dtype=np.int32)
            for ll in range(10000):                               
                brd_mtx[brd_pts[ll, 0], brd_pts[ll, 1], brd_pts[ll, 2]]  =  1           # 3D matrix with the 

            frames2work  =  np.where(np.sum(brd_mtx, axis=(1,2)) != 0)[0]               # select frames with overlapping
            brd_elpsd    =  np.zeros(brd_mtx.shape, dtype=np.int32)
            for z in frames2work:                                                       # frame by in z-stack
                xy_coords  =  np.transpose(np.asarray(np.where(brd_mtx[z] == 1)))       # coordinate of the ellipsoid border in the considered frame
                if xy_coords.shape[0] > 2:
                    hull          =  ConvexHull(xy_coords)                                  # convex hull surface
                    x, y          =  np.meshgrid(np.arange(xlen), np.arange(ylen))          # make a canvas with coordinates
                    x, y          =  x.flatten(), y.flatten()                               # adapting coordinates
                    points        =  np.vstack((y,x)).T                                     
                    p             =  Path(hull.points[hull.vertices])
                    grid          =  p.contains_points(points)                              # points inside the contour
                    brd_elpsd[z]  =  grid.reshape(xlen, ylen)                               # reshape (last commands come from a tutorial, I don't know details, but works)
                
            if nucs_fitting_flag == "Semi-Ellip":
                brd_elpsd  =  EllipsoidCylinder.EllipsoidCylinder(brd_elpsd).semi_cyl
                
            idxs2splt  =  np.unique(brd_elpsd * nucs_elips)[1:].astype(np.int)          # checking if the ellipsoid just found has overlapping with the other

            if idxs2splt.size > 0:                                                      # if there is overlapping, all the overlapping ellipsoids are organized in couples and on each couple the SplitEllipsoid class will work
                for ii in idxs2splt:
                    elips2splt   =  (nucs_elips == ii) * 1  +  brd_elpsd * 2                    # first couple of ellipsoids to split
                    new_elips    =  SplitEllipsoids.SplitEllipsoids(elips2splt).ellipsoids      # splitting
                    nucs_elips  *=  (1 - (nucs_elips == ii))                                    # the ellipsoid already done is removed from the output matrix
                    nucs_elips  +=  ii * (new_elips == 1)                                       # the ellipsoid already done is updated in the output matrix with the cut to remove overlapping
                    brd_elpsd    =  (new_elips == 2)                                            # the starting ellipsoid is updated with the cut. This goes on for all the ellipsoids involved
                
                nucs_elips  +=  k * (new_elips == 2)                                    # if no more ellipsoids needs to be cut, the starting ellipsoid is added to the output matrix

            else:
                nucs_elips  +=  k * brd_elpsd                                           # if there is no overlapping, the starting ellipsoid is directly added to the output with the proper tag
            pbar.close()


        self.nucs_elips  =  nucs_elips
                
                
                
class ProgressBar(QtWidgets.QWidget):
    """Simple progress bar widget"""
    def __init__(self, parent=None, total1=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar  =  QtWidgets.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total1)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        """update the progressbar"""
        self.progressbar.setValue(val1)
        QtWidgets.qApp.processEvents()
