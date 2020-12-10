import numpy as np
from skimage.morphology import label
from skimage.morphology import binary_dilation
from skimage.measure import regionprops

import CloserNucleiFinder


class ProjectConnect:
    def __init__(self, nucs_det_lbls2, spts):

        steps, xlen, ylen  =  nucs_det_lbls2.shape

        nucs_det_lbls  =  np.zeros(nucs_det_lbls2.shape)
        for t in range(steps):                                                                                          # correction for drift during the acquisition
            nucs_det_lbls[t, :xlen - t*2, :]  =  nucs_det_lbls2[t, t*2:, :]

        idxs      =  np.unique(nucs_det_lbls)[1:]
        nucs_mip  =  np.zeros(nucs_det_lbls[0, :, :].shape)
        for k in idxs:                                                                                                  # mip of the nucs_det_lbls: it projects all nuclei in 2D leaving them the same tag as before.
            nucs_mip  +=  k * np.sign((nucs_det_lbls == k).sum(0)) * (1 - np.sign(nucs_mip))                            # The last multiplication is to avoid overlapping

        spts_lbls  =  label(spts)
        rgp        =  regionprops(spts_lbls)
        ref_spts   =  np.zeros(len(rgp))
        for k in range(spts_lbls.max()):
            ref_spts[k]  =  CloserNucleiFinder.CloserNucleiFinder(nucs_mip, np.asarray(rgp[k]['centroid'][1:], dtype=int)).mx_pt

        fls_clrd  =  np.zeros(np.append(nucs_mip.shape, 3))
        for j in range(ref_spts.size):
            fls_clrd[:, :, 0]  +=  (nucs_mip == ref_spts[j]) * 1

        fls_clrd[:, :, 2]  =  np.sign(nucs_mip) * (1 - np.sign(fls_clrd[:, :, 0])) * fls_clrd[:, :, 0].max()
        fls_clrd[:, :, 1]  =  spts.sum(0) * fls_clrd[:, :, 0].max()

        self.nucs_det_lbls  =  nucs_det_lbls
        self.nucs_mip  =  nucs_mip
        self.fls_clrd  =  fls_clrd



class ProjectConnectByDilation:
    def __init__(self, nucs_det_lbls, spts, spts_ctrs):


        idxs      =  np.unique(nucs_det_lbls)[1:]
        nucs_mip  =  np.zeros(nucs_det_lbls[0, :, :].shape)
        for k in idxs:                                                                                                  # mip of the nucs_det_lbls: it projects all nuclei in 2D leaving them the same tag as before.
            nucs_mip  +=  k * np.sign((nucs_det_lbls == k).sum(0)) * (1 - np.sign(nucs_mip))                            # The last multiplication is to avoid overlapping

        idxs  =  np.unique(nucs_mip)[1:]

        nucs_dil  =  np.copy(nucs_mip)
        while (nucs_dil == 0).sum() > 0:
            for k in idxs:
                nuc_bff   =   binary_dilation(nucs_dil == k)
                nucs_dil  *=  (1 - (nucs_dil == k))
                nucs_dil  +=  k * nuc_bff * (1 - np.sign(nucs_dil))


        # spts_lbls  =  label(spts)
        # rgp        =  regionprops(spts_lbls)
        ref_spts   =  np.zeros(spts_ctrs.shape[1])
        for k in range(spts_ctrs.shape[1]):
            ref_spts[k]  =  nucs_dil[spts_ctrs[1, k], spts_ctrs[0, k]]

        fls_clrd  =  np.zeros(np.append(nucs_mip.shape, 3))
        for j in range(ref_spts.size):
            fls_clrd[:, :, 0]  +=  (nucs_mip == ref_spts[j]) * 1

        fls_clrd[:, :, 2]  =  np.sign(nucs_mip) * (1 - np.sign(fls_clrd[:, :, 0])) * fls_clrd[:, :, 0].max()
        fls_clrd[:, :, 1]  =  spts.sum(0) * fls_clrd[:, :, 0].max()

        self.nucs_det_lbls  =  nucs_det_lbls
        self.nucs_dil       =  nucs_dil
        self.fls_clrd       =  fls_clrd






class ProjectConnectByDilationBackup:
    def __init__(self, nucs_det_lbls, spts, spts_ctrs):

        # steps, xlen, ylen  =  nucs_det_lbls.shape

        # nucs_det_lbls  =  np.zeros(nucs_det_lbls2.shape)
        # for t in range(steps):                                                                                          # correction for drift during the acquisition
        #     nucs_det_lbls[t, :xlen - t*2, :]  =  nucs_det_lbls2[t, t*2:, :]


        idxs      =  np.unique(nucs_det_lbls)[1:]
        nucs_mip  =  np.zeros(nucs_det_lbls[0, :, :].shape)
        for k in idxs:                                                                                                  # mip of the nucs_det_lbls: it projects all nuclei in 2D leaving them the same tag as before.
            nucs_mip  +=  k * np.sign((nucs_det_lbls == k).sum(0)) * (1 - np.sign(nucs_mip))                            # The last multiplication is to avoid overlapping

        idxs  =  np.unique(nucs_mip)[1:]

        nucs_dil  =  np.copy(nucs_mip)
        while (nucs_dil == 0).sum() > 0:
            for k in idxs:
                nuc_bff   =   binary_dilation(nucs_dil == k)
                nucs_dil  *=  (1 - (nucs_dil == k))
                nucs_dil  +=  k * nuc_bff * (1 - np.sign(nucs_dil))

        # borders  =  np.zeros(nucs_dil.shape)                                                                          # detects the line that separate nuclei. It is here just in case
        # for k in idxs:
        #     borders  +=  binary_dilation(nucs_dil == k)
        # borders  =  (borders > 1)


        spts_lbls  =  label(spts)
        rgp        =  regionprops(spts_lbls)
        ref_spts   =  np.zeros(spts_ctrs.shape[1])
        for k in range(spts_ctrs.shape[1]):
            ref_spts[k]  =  nucs_dil[spts_ctrs[1, k], spts_ctrs[0, k]]

        fls_clrd  =  np.zeros(np.append(nucs_mip.shape, 3))
        for j in range(ref_spts.size):
            fls_clrd[:, :, 0]  +=  (nucs_mip == ref_spts[j]) * 1

        fls_clrd[:, :, 2]  =  np.sign(nucs_mip) * (1 - np.sign(fls_clrd[:, :, 0])) * fls_clrd[:, :, 0].max()
        fls_clrd[:, :, 1]  =  spts.sum(0) * fls_clrd[:, :, 0].max()

        self.nucs_det_lbls  =  nucs_det_lbls
        self.nucs_dil       =  nucs_dil
        self.fls_clrd       =  fls_clrd
