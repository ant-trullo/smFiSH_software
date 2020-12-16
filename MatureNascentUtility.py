import numpy as np



class MatureNascent:
    def __init__(self, spts_lbls, spts_ints_ctrs):

        mtx2thr  =  np.copy(spts_lbls)

        for k in range(spts_ints_ctrs.shae[1]):
            mtx2thr[mtx2thr == spts_ints_ctrs[0, k]]  =  spts_ints_ctrs[1, k]

        self.mtx2thr  =  mtx2thr
