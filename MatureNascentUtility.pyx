import numpy as np
cimport numpy as np

#import pyximport
#pyximport.install(reload_support=True)

DTYPE = np.int32
ctypedef np.int32_t DTYPE_t

cpdef MatureNascentUtility(np.ndarray[DTYPE_t, ndim=3] spts_lbls, np.ndarray[DTYPE_t, ndim=2] spts_ints_ctrs):
    cdef int zlen = spts_lbls.shape[0]
    cdef int xlen = spts_lbls.shape[1]
    cdef int ylen = spts_lbls.shape[2]
    cdef int nlbs = spts_ints_ctrs.shape[1]
    cdef np.ndarray mtx2thr = np.zeros([zlen, xlen, ylen], dtype=DTYPE)


    for z in range(zlen):
        print(z)
        for x in range(xlen):
            for y in range(ylen):
                # print(spts_lbls[z, x, y])
                if spts_lbls[z, x, y] != 0:
                    k                 =  np.where(spts_ints_ctrs[0, :] == spts_lbls[z, x, y])
                    # print(k, z, x, y, spts_lbls[z, x, y])
                    mtx2thr[z, x, y]  =  spts_ints_ctrs[1, k[0][0]]
                    
    return mtx2thr
