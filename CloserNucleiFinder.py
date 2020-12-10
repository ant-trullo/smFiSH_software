"""This function finds the closer connected component to thestarting point.

Input is 'nuclei_lbl', which is the matrix-image with
connected labeled objects, and points, the coordinates of the point. It
returns the label value of the closest object.
To do that, it starts from a square centred in point and with a side of 1.
If the maximum in the square is zero, the side is increased and the maximum
is searched into a broader square. The side is kept smaller than 5.
"""

import numpy as np


class CloserNucleiFinder:
    def __init__(self, nuclei_lbl, point):
        self.nuclei_lbl  =  nuclei_lbl
        self.point       =  point

        mx_pt  =  0
        side   =  1
        edge_x, edge_y   =  self.nuclei_lbl.shape

        while mx_pt == 0 and side < 50:
            mx_pt  =   (self.nuclei_lbl[np.max((self.point[0] - side, 0)):np.min((self.point[0] + side, edge_x - 1)) + 1, np.max((self.point[1] - side, 0)):np.min((self.point[1] + side, edge_y - 1)) + 1]).max()
            side   +=  1

        self.mx_pt  =  mx_pt
