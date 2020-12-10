"""This function determines the area of the cylindroids of each nucleus.

Given a Z-X-Y matrix in which nuclei are segmented frame by frame, this
function sums up the segmentation in Z and than, working with dilation,
expands in 2D the calc of nuclei up to the filling of all the space. In
this way every 2D point can be associated to a nucleus.
"""

import numpy as np
import multiprocessing
from scipy.ndimage.morphology import binary_erosion
from skimage.morphology import label, remove_small_objects

import DilationFunction
import RemoveBorderNuclei


class NucleiArea:
    """Main class, does all the job"""
    # def __init__(self, nucs_lbls, spts_ctrs, spts):
    def __init__(self, nucs_lbls):

        zlen, xlen, ylen  =  nucs_lbls.shape

        nucs_er  =  np.zeros(nucs_lbls.shape, dtype=int)
        for z in range(zlen):
            nucs_er[z, :, :]  =  label(binary_erosion(nucs_lbls[z, :, :], iterations=6), connectivity=1)               # the matrix of segmented nuclei is eroded first
            nucs_er[z, :, :]  =  remove_small_objects(nucs_er[z, :, :], 120)                                           # small connected components are deleted

        j_max     =  np.argmax(np.sign(nucs_er).sum(2).sum(1))                                                          # find z with the highest surface of detected nuclei
        nucs_prj  =  np.copy(nucs_er[j_max, :, :])
        new_idx   =  1000
        for z in range(1, zlen):
            idxs  =  np.unique(nucs_er[z, :, :])[1:]                                                                   # 
            for k in idxs:
                smpl  =  ((nucs_er[z, :, :] == k) * nucs_prj).reshape((xlen * ylen))

                if smpl.sum() == 0:
                    nucs_prj  +=  new_idx * (nucs_er[z, :, :] == k)
                    new_idx   +=  1

                else:
                    smpl      =   np.delete(smpl, np.where(smpl == 0)[0])
                    smpl      =   np.median(smpl).astype(int)
                    nucs_prj  +=  smpl * (nucs_er[z, :, :] == k) * (1 - np.sign(nucs_prj))

        nucs_prj  =  np.copy(nucs_er[j_max, :, :])
        nucs_prj  =  label(nucs_prj)

        cpu_owe   =  multiprocessing.cpu_count()
        nucs_dil  =  np.copy(nucs_prj)

        while (nucs_dil == 0).sum() > 0:
            print((nucs_dil == 0).sum())
            idxs      =  np.unique(nucs_dil)[1:]
            tags      =  np.split(idxs, (int(idxs.size / cpu_owe + 1)) * np.arange(1, cpu_owe, dtype=np.int))
            args      =  []
            for jj in range(len(tags)):
                args.append([tags[jj], nucs_dil])
            pool     =  multiprocessing.Pool()
            results  =  pool.map(DilationFunction.DilationFunction, args)
            pool.close()
            for j in range(len(results)):
                nucs_dil  +=  results[j].msk * (1 - np.sign(nucs_dil))
            nucs_dil  =  label(nucs_dil)

        nucs_dil    =  RemoveBorderNuclei.RemoveBorderNuclei(nucs_dil).nucs_dil

        # ref_spts  =  np.zeros(spts_ctrs.shape[1])
        # for k in range(spts_ctrs.shape[1]):
        #     ref_spts[k]  =  nucs_dil[spts_ctrs[1, k], spts_ctrs[0, k]]

        # fls_clrd  =  np.zeros(np.append(nucs_prj.shape, 3))
        # for j in range(ref_spts.size):
        #     fls_clrd[:, :, 0]  +=  (nucs_prj == ref_spts[j]) * 1
        #
        # fls_clrd[:, :, 2]  =  np.sign(nucs_prj) * (1 - np.sign(fls_clrd[:, :, 0])) * fls_clrd[:, :, 0].max()
        # fls_clrd[:, :, 1]  =  spts.sum(0) * fls_clrd[:, :, 0].max()


        self.nucs_lbls  =  nucs_lbls
        self.nucs_dil   =  nucs_dil
        # self.fls_clrd   =  fls_clrd
