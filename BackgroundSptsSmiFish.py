"""This function defines cells around spots to evaluate the background for each spot."""


import multiprocessing
import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion
# from skimage.measure import regionprops_table




class BackgroundSptsSmiFish:
    """Main class, organizes input data for multiprocessing purposes"""
    def __init__(self, spts_lbls):

        idxs       =  np.unique(spts_lbls)[1:]
        cpu_ow     =  multiprocessing.cpu_count()
        jobs_idxs  =  np.array_split(idxs, cpu_ow)
        jobs_args  =  []
        for j in range(cpu_ow):
            jobs_args.append([spts_lbls, jobs_idxs[j]])

        pool     =  multiprocessing.Pool()
        results  =  pool.map(BackgroundSptsSmiFishUtility, jobs_args)
        pool.close()

        cages    =  np.zeros(spts_lbls.shape, dtype=np.int32)
        tot_msk  =  np.zeros(spts_lbls.shape, dtype=np.int8)
        for k in range(len(results)):
            cages    +=  results[k].cages
            tot_msk  +=  results[k].tot_msk

        cages  *=  (tot_msk == 1)

        self.cages  =  cages



class BackgroundSptsSmiFishUtility:
    """class Utility, to be implemented in a multicore pool"""
    def __init__(self, input_arg):

        cages    =  np.zeros(input_arg[0].shape, dtype=np.int32)
        tot_msk  =  np.zeros(input_arg[0].shape, dtype=np.int8)

        for k in input_arg[1]:
            sing_spt  =   (input_arg[0] == k)
            bffr      =   binary_dilation(sing_spt, iterations=5)
            cages     +=  k * (bffr ^ binary_erosion(bffr, iterations=2))
            tot_msk   +=  bffr * 1

        self.cages    =  cages
        self.tot_msk  =  tot_msk
