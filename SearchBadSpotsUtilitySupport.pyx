import numpy as np
cimport numpy as np
#import pyximport
#pyximport.install(reload_support=True)

DTYPE = np.int32
ctypedef np.int32_t DTYPE_t

 
cpdef unique_cython_int_1dim(np.ndarray[DTYPE_t, ndim=1] a):
    cdef int x
    cdef int xlen = a.size
    cdef set s = set()
    cdef list s_list = []
    
    for x in range(xlen):
        s.add(a[x])
    for s_el in s:
        s_list.append(s_el)

    return s_list



#cpdef SearchBadSpotsUtilitySupport(np.ndarray[long long, ndim=3] a, list l, list v, list c, int vol, int num_frms):
cpdef SearchBadSpotsUtilitySupport(np.ndarray[DTYPE_t, ndim=3] a, list l, list v, list c, int vol, int num_frms):
    cdef int len_b  =  len(l)
    cdef int zlen   =  a.shape[0]
    cdef int xlen   =  a.shape[1]
    cdef int ylen   =  a.shape[2]
    cdef int tag_ref
    cdef np.ndarray msk = np.zeros([zlen, xlen, ylen], dtype=DTYPE)

    for j in range(len_b):                                           
        if v[j] < vol or len(unique_cython_int_1dim(c[j])) < num_frms:      # the constrain on the coordinates asks to a spots to be 3D, not 2D
#            print(v[j])
            tag_ref  =  l[j]
            for z in range(zlen):
                for x in range(xlen):
                    for y in range(ylen):
                        if a[z, x, y] == tag_ref:
                            msk[z, x, y]  +=  1
    
    return msk
