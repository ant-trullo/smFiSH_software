"""This function loads raw data for smiFish analysis.

Raw data are .czi. The input is the file name and the channel numbers,
the output are the matrices for the multi z nuclei and the multi z
for spots channels, plus the mip of the nuclei channel. If the nuclei channel is
not asked, the nucs_mip is going to be a zeros array.
chs_spts_nucs is a list that has 3 channel numbers for spots_a, spots_b and nuclei 
respectively. If one of this number is -1, the relative channel is absent.
"""



import numpy as np
import czifile


class RawDataLoader:
    """Main class, does all the job"""
    def __init__(self, raw_data_fname, chs_spts_nucs):

        a      =  czifile.CziFile(raw_data_fname)                                                                                   # read info about pixel size
        b      =  a.metadata()

        start           =  b.find("ScalingZ")
        end             =  b[start + 9:].find("ScalingZ")
        self.pix_sizeZ  =  np.round(float(b[start + 9:start + 7 + end]) * 1000000, decimals=4) 

        start           =  b.find("ScalingX")
        end             =  b[start + 9:].find("ScalingX")
        self.pix_sizeX  =  np.round(float(b[start + 9:start + 7 + end]) * 1000000, decimals=4) 

        filedata  =  np.squeeze(czifile.imread(raw_data_fname))

        if chs_spts_nucs[0] != -1:
            self.spts_a  =  filedata[chs_spts_nucs[0], :, :, :]
            for z in range(self.spts_a.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                self.spts_a[z, :, :]  =  np.rot90(self.spts_a[z, :, ::-1])


        if chs_spts_nucs[1] != -1:
            self.spts_b    =  filedata[chs_spts_nucs[1], :, :, :]
            for z in range(self.spts_b.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                self.spts_b[z, :, :]  =  np.rot90(self.spts_b[z, :, ::-1])

        if chs_spts_nucs[2] != -1:
            self.nucs    =  filedata[chs_spts_nucs[2], :, :, :]
            for z in range(self.nucs.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
                self.nucs[z, :, :]  =  np.rot90(self.nucs[z, :, ::-1])

            self.nucs_mip  =  np.zeros(self.nucs.shape[1:])                 # nuclei mip image
            for x in range(self.nucs.shape[1]):
                self.nucs_mip[x, :]    =  filedata[chs_spts_nucs[2], :, x, :].max(0)
                
        else:                                                               # if nuclei channel is not asked, it will be a zeros matrix with the proper shape
            if chs_spts_nucs[0] != -1:
                self.nucs_mip  =  np.zeros(self.spts_a.shape[1:])
            elif chs_spts_nucs[1] != -1:
                self.nucs_mip  =  np.zeros(self.spts_b.shape[1:])





#class RawDataLoader:
#    def __init__(self, raw_data_fname, chs_spts_nucs):
#
#        a      =  czifile.CziFile(raw_data_fname)                                                                                   # read info about pixel size
#        b      =  a.metadata()
#
#        start           =  b.find("ScalingZ")
#        end             =  b[start + 9:].find("ScalingZ")
#        self.pix_sizeZ  =  np.round(float(b[start + 9:start + 7 + end]) * 1000000, decimals=4) 
#
#        start           =  b.find("ScalingX")
#        end             =  b[start + 9:].find("ScalingX")
#        self.pix_sizeX  =  np.round(float(b[start + 9:start + 7 + end]) * 1000000, decimals=4) 
#
#        filedata  =  np.squeeze(czifile.imread(raw_data_fname))
#
#        if chs_spts_nucs.size == 2:                             # check the number of channels of the raw data
#            spts_a    =  filedata[chs_spts_nucs[0], :, :, :]
#            nucs      =  filedata[chs_spts_nucs[1], :, :, :]
#
#            for z in range(nucs.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
#                spts_a[z, :, :]  =  np.rot90(spts_a[z, :, ::-1])
#                nucs[z, :, :]    =  np.rot90(nucs[z, :, ::-1])
#
#            self.spts_a  =  spts_a
#            self.nucs    =  nucs
#
#        if chs_spts_nucs.size == 3:
#            spts_a    =  filedata[chs_spts_nucs[0], :, :, :]
#            spts_b    =  filedata[chs_spts_nucs[1], :, :, :]
#            nucs      =  filedata[chs_spts_nucs[2], :, :, :]
#
#            for z in range(nucs.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
#                spts_a[z, :, :]  =  np.rot90(spts_a[z, :, ::-1])
#                spts_b[z, :, :]  =  np.rot90(spts_b[z, :, ::-1])
#                nucs[z, :, :]    =  np.rot90(nucs[z, :, ::-1])
#
#            self.spts_a     =  spts_a
#            self.spts_b     =  spts_b
#            self.nucs       =  nucs
##            self.pix_sizeX  =  pix_sizeX
##            self.pix_sizeZ  =  pix_sizeZ
#

# class RawDataLoader:
#     def __init__(self, raw_data_fname):

#         if raw_data_fname[-3:] == 'lsm' or raw_data_fname[-3:] == 'tif':
#             filedata   =  tifffile.imread(raw_data_fname)
#             spts       =  filedata[:, 1, :, :]
#             nucs       =  filedata[:, 0, :, :]

#         if raw_data_fname[-3:] == 'czi':
#             filedata  =  np.squeeze(czifile.imread(raw_data_fname))
#             spts      =  filedata[0, :, :, :]
#             nucs      =  filedata[0, :, :, :]

#         for z in range(spts.shape[0]):                          # both are rotated and mirrored to match ImageJ standards
#             spts[z, :, :]  =  np.rot90(spts[z, :, ::-1])
#             nucs[z, :, :]  =  np.rot90(nucs[z, :, ::-1])


#         self.spts  =  spts
#         self.nucs  =  nucs
