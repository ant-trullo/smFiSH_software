"""This function organizes information of mRNA splitting in above and below, single mol and clusters.
Inputs are the matrices of the segmented spots and the relative raw data matrix, plus bondaries of above and below
and intensity thresholds for both cases.
"""


from os import listdir
import numpy as np
from skimage.measure import regionprops_table
import xlwt
from xlrd import open_workbook

import RawDataLoader

# analysis_folder  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select the Directory with the Analysis"))
# raw_data_fname   =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data file of the analyze", filter='*.czi')[0])


class Quantify_mRNA:
    """Only class, does all the job"""
    def __init__(self, analysis_folder, raw_data_fname, int_thr_above, int_thr_below, vol_thr, clstr_dist):

        first_last_frame  =  np.load(analysis_folder + '/first_last_frame.npy')
        spts_red          =  np.load(analysis_folder + '/spts_lbls_b.npy')
        chs_spts_nucs     =  np.fromfile(analysis_folder + '/chs_spts_nucs.bin', 'uint16')
        raw_data          =  RawDataLoader.RawDataLoader(raw_data_fname, chs_spts_nucs)

        files  =  listdir(analysis_folder)                   # search in the analysis folder the xls files to read the background values
        for file in files:
            if file[-6:] == "_b.xls":                        # # select the .xls file with spot analysis results, red
                xls_b_name  =  file
                break

        book2read_b  =  open_workbook(analysis_folder + "/" + xls_b_name)     # load xls
        sheet_b      =  book2read_b.sheet_by_index(0)                         # copy tags from xls into a python-list
        tags_list_b  =  sheet_b.col_values(0)                                 # remove first 3 entries - are just titles
        tags_list_b  =  tags_list_b[3:]                                       # list of strings ("Spt_ ...")
        tags_list_b  =  [int(tag[4:])for tag in tags_list_b]                  # list of tags as integer numbers
        bkg_list_b   =  sheet_b.col_values(4)[3:]                             # copy background values in a python-list

        rgp_red  =  regionprops_table(spts_red, raw_data.spts_b, properties=["label", "intensity_image", "area", "centroid"])    # table of the needed properties of all the spots

        clstr_above  =  []                      # list for the properties of all the clustered spots above nucs 
        clstr_below  =  []                      # list for the properties of all the clustered spots below nucs
        sm_above     =  []                      # list for the properties of all the single mol spots above nucs
        sm_below     =  []                      # list for the properties of all the single mol spots below nucs


        for k in range(rgp_red["area"].size):
            if rgp_red["centroid-0"][k] <= first_last_frame[0]:
                if rgp_red["intensity_image"][k].sum() <= int_thr_above:
                    bkg      =  bkg_list_b[tags_list_b.index(rgp_red["label"][k])]
                    spt_int  =  rgp_red["intensity_image"][k].sum()
                    sm_above.append([int(rgp_red["label"][k]), int(spt_int), bkg, spt_int - bkg * rgp_red["area"][k], spt_int / bkg, int(rgp_red["area"][k]), int(rgp_red["centroid-0"][k]), int(rgp_red["centroid-1"][k]), int(rgp_red["centroid-2"][k])])
                else:
                    bkg      =  bkg_list_b[tags_list_b.index(rgp_red["label"][k])]
                    spt_int  =  rgp_red["intensity_image"][k].sum()
                    clstr_above.append([int(rgp_red["label"][k]), int(spt_int), bkg, spt_int - bkg * rgp_red["area"][k], spt_int / bkg, int(rgp_red["area"][k]), int(rgp_red["centroid-0"][k]), int(rgp_red["centroid-1"][k]), int(rgp_red["centroid-2"][k])])
            else:
                if rgp_red["intensity_image"][k].sum() <= int_thr_below:
                    bkg      =  bkg_list_b[tags_list_b.index(rgp_red["label"][k])]
                    spt_int  =  rgp_red["intensity_image"][k].sum()
                    sm_below.append([int(rgp_red["label"][k]), int(spt_int), bkg, spt_int - bkg * rgp_red["area"][k], spt_int / bkg, int(rgp_red["area"][k]), int(rgp_red["centroid-0"][k]), int(rgp_red["centroid-1"][k]), int(rgp_red["centroid-2"][k])])
                else:
                    bkg      =  bkg_list_b[tags_list_b.index(rgp_red["label"][k])]
                    spt_int  =  rgp_red["intensity_image"][k].sum()
                    clstr_below.append([int(rgp_red["label"][k]), int(spt_int), bkg, spt_int - bkg * rgp_red["area"][k], spt_int / bkg, int(rgp_red["area"][k]), int(rgp_red["centroid-0"][k]), int(rgp_red["centroid-1"][k]), int(rgp_red["centroid-2"][k])])

        book    =  xlwt.Workbook(encoding='utf-8')
        sheet1  =  book.add_sheet("Above")
        sheet2  =  book.add_sheet("Below")

        sheet1.write(0, 0, "Single Mol")
        sheet1.write(1, 0, "Spt_id")
        sheet1.write(1, 1, "Ints")
        sheet1.write(1, 2, "bkg")
        sheet1.write(1, 3, "Ints - bkg")
        sheet1.write(1, 4, "Ints / bkg")
        sheet1.write(1, 5, "Vols")
        sheet1.write(1, 6, "z-coord")
        sheet1.write(1, 7, "x-coord")
        sheet1.write(1, 8, "y-coord")

        sheet1.write(0, 10, "Clusters")
        sheet1.write(1, 10, "Spt_id")
        sheet1.write(1, 11, "Ints")
        sheet1.write(1, 12, "bkg")
        sheet1.write(1, 13, "Ints - bkg")
        sheet1.write(1, 14, "Ints / bkg")
        sheet1.write(1, 15, "Vols")
        sheet1.write(1, 16, "z-coord")
        sheet1.write(1, 17, "x-coord")
        sheet1.write(1, 18, "y-coord")

        sheet2.write(0, 0, "Single Mol")
        sheet2.write(1, 0, "Spt_id")
        sheet2.write(1, 1, "Ints")
        sheet2.write(1, 2, "bkg")
        sheet2.write(1, 3, "Ints - bkg")
        sheet2.write(1, 4, "Ints / bkg")
        sheet2.write(1, 5, "Vols")
        sheet2.write(1, 6, "z-coord")
        sheet2.write(1, 7, "x-coord")
        sheet2.write(1, 8, "y-coord")

        sheet2.write(0, 10, "Clusters")
        sheet2.write(1, 10, "Spt_id")
        sheet2.write(1, 11, "Ints")
        sheet2.write(1, 12, "bkg")
        sheet2.write(1, 13, "Ints - bkg")
        sheet2.write(1, 14, "Ints / bkg")
        sheet2.write(1, 15, "Vols")
        sheet2.write(1, 16, "z-coord")
        sheet2.write(1, 17, "x-coord")
        sheet2.write(1, 18, "y-coord")

        idx_s_a  =  2
        for k_s_a in range(len(sm_above)):
            sheet1.write(idx_s_a, 0, sm_above[k_s_a][0])
            sheet1.write(idx_s_a, 1, sm_above[k_s_a][1])
            sheet1.write(idx_s_a, 2, sm_above[k_s_a][2])
            sheet1.write(idx_s_a, 3, sm_above[k_s_a][3])
            sheet1.write(idx_s_a, 4, sm_above[k_s_a][4])
            sheet1.write(idx_s_a, 5, sm_above[k_s_a][5])
            sheet1.write(idx_s_a, 6, sm_above[k_s_a][6])
            sheet1.write(idx_s_a, 7, sm_above[k_s_a][7])
            sheet1.write(idx_s_a, 8, sm_above[k_s_a][8])

            idx_s_a  +=  1

        idx_c_a  =  2
        for k_c_a in range(len(clstr_above)):
            sheet1.write(idx_c_a, 10, clstr_above[k_c_a][0])
            sheet1.write(idx_c_a, 11, clstr_above[k_c_a][1])
            sheet1.write(idx_c_a, 12, clstr_above[k_c_a][2])
            sheet1.write(idx_c_a, 13, clstr_above[k_c_a][3])
            sheet1.write(idx_c_a, 14, clstr_above[k_c_a][4])
            sheet1.write(idx_c_a, 15, clstr_above[k_c_a][5])
            sheet1.write(idx_c_a, 16, clstr_above[k_c_a][6])
            sheet1.write(idx_c_a, 17, clstr_above[k_c_a][7])
            sheet1.write(idx_c_a, 18, clstr_above[k_c_a][8])

            idx_c_a  +=  1

        idx_s_b  =  2
        for k_s_b in range(len(sm_below)):
            sheet2.write(idx_s_b, 0, sm_below[k_s_b][0])
            sheet2.write(idx_s_b, 1, sm_below[k_s_b][1])
            sheet2.write(idx_s_b, 2, sm_below[k_s_b][2])
            sheet2.write(idx_s_b, 3, sm_below[k_s_b][3])
            sheet2.write(idx_s_b, 4, sm_below[k_s_b][4])
            sheet2.write(idx_s_b, 5, sm_below[k_s_b][5])
            sheet2.write(idx_s_b, 6, sm_below[k_s_b][6])
            sheet2.write(idx_s_b, 7, sm_below[k_s_b][7])
            sheet2.write(idx_s_b, 8, sm_below[k_s_b][8])

            idx_s_b  +=  1

        idx_c_b  =  2
        for k_c_b in range(len(clstr_below)):
            sheet2.write(idx_c_b, 10, clstr_below[k_c_b][0])
            sheet2.write(idx_c_b, 11, clstr_below[k_c_b][1])
            sheet2.write(idx_c_b, 12, clstr_below[k_c_b][2])
            sheet2.write(idx_c_b, 13, clstr_below[k_c_b][3])
            sheet2.write(idx_c_b, 14, clstr_below[k_c_b][4])
            sheet2.write(idx_c_b, 15, clstr_below[k_c_b][5])
            sheet2.write(idx_c_b, 16, clstr_below[k_c_b][6])
            sheet2.write(idx_c_b, 17, clstr_below[k_c_b][7])
            sheet2.write(idx_c_b, 18, clstr_below[k_c_b][8])

            idx_c_b  +=  1


        book.save(analysis_folder + "/" + xls_b_name[:-22] + "_mRNA_SM_Clusters" + str(vol_thr) + "_clstrDist_" + str(clstr_dist) + ".xls")
        # book.save('/home/atrullo/Desktop' + "/" + xls_b_name[:-22] + "_mRNA_SM_Clusters" + str(vol_thr) + "_clstrDist_" + str(clstr_dist) + ".xls")



