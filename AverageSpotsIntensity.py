"""This function calculates the average intensity of the detected spots.EnginioClient

Starting from spots (segmented and labeled) and the spots raw data , it records 
all the average intensities and gives as output the matrix of labeled spots with
their average intensity instead of the tag. It uses regionprops to do the job, 
plus multiprocessing for speed.

MAYBE WE CAN AVOID TO DEFINE A MATRIX WITH AVERAGE INTENSITY INSTEAD OF TAGS
AND GO FOR A LIST OF INTENSITY VALUES WHICH CAN BE FASTER FOR THE ALGORITHM.

"""


import numpy as np
import multiprocessing
from skimage.measure import regionprops



class AverageSpotsIntensity:
    def __init__(self, spts_lbls, raw_spts):

        rgp        =  regionprops(spts_lbls, raw_spts           # regionprops to collect info on spots
        idxs_tot   =  np.arange(1, spts_lbls.max() + 1)         # list of all the spot tags

        cpu_ow  =  multiprocessing.cpu_count()
        idx_gr  =  np.array_split(idxs_tot, cpu_ow)

        jobs_args  =  []
        for k in range(cpu_ow):
            jobs_args.append([idx_gr[k], rgp, spts_lbls])       # chop of input and organization of data  (setq evil-emacs-state-cursor '("green" bar))in list  for multiprocessing purposes

        av_int_map  =  np.zeros(spts_lbls.shape)

        pool     =  multiprocessing.Pool()
        results  =  pool.map(Func4Mp, jobs_args)
        pool.close()

        for jj in range(1, len(results)):                       # mergin results of multiprocessing
            av_int_map  +=  results[jj].int_map

        self.av_int_map  =  av_int_map



class Func4Mp:
    def __init__(self, args):

        int_map  =  np.zeros(args[2].shape)                             # matrix initialization

        for idx in args[0]:                                                                                             # The output matrix changes the label of the input with the average intensity of the spots.
            int_map  +=  (args[2] == args[1][idx - 1]['label']) * args[1][idx - 1]['intensity_image'].sum() / args[1][idx - 1]['area']  # Average intensity is callculated with 'intensity_image' divided by the 'area', all built-in of regionprops. 
 
        self.int_map  =  int_map
