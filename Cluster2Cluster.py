"""This function tries to correlate the mRNA cluster intensity with the overlapping protein intensity.

The starting point is the analysis folder with generated xls. 
"""


from os import listdir
import numpy as np
import xlwt
from xlrd import open_workbook
from skimage.measure import regionprops_table


class Cluster2Cluster:
    """Only class, does all the job"""
    def __init__(self, analysis_folder):

        global xls_b_name
        spts_red    =  np.load(analysis_folder + '/spts_lbls_b.npy')                        # load segmented red spots
        spts_green  =  np.load(analysis_folder + '/spts_lbls_a.npy')                        # load segmented green spots
        rgp_red     =  regionprops_table(spts_red, properties=["label", "coords"])          # regionprops of red spots to get coordinates and labels


        first_last_frame  =  np.load(analysis_folder + '/first_last_frame.npy')             # first and last considered z-frames

        files  =  [f for f in listdir(analysis_folder) if not f.startswith('.')]            # list of files in analysis folder (exclude hidden files in case the xls file is opened gives troubles at least under linux)
        for file in files:
            if "Spots_intensity_b" in file:                                                 # select the .xls file with spot analysis results, red
                xls_b_name  =  file
                break

        book_b      =  open_workbook(analysis_folder + "/" + xls_b_name)                    # load xls file
        sheet_b     =  book_b.sheet_by_index(0)
        red_tags    =  sheet_b.col_values(0)[3:]                                            # import tags as a list (tags are in column 3, strating from the raw 3)
        red_tags    =  np.asarray([np.int32(tag[4:]) for tag in red_tags])                  # tags are strings like "Spt_45", extract 45 as a int32 number
        red_ints    =  sheet_b.col_values(2)[3:]                                            # import intensity of the spots from xls as a list
        red_ints    =  np.asanyarray([np.int32(inntt) for inntt in red_ints])               # convert intensities in integers
        red_zcoors  =  sheet_b.col_values(8)[3:]                                            # import z-coords as a list
        red_zcoors  =  np.asarray([np.int32(zcoord) for zcoord in red_zcoors])              # convert z-coords in integers

        idxs_above  =  np.where(red_zcoors <= first_last_frame[0])[0]                       # indexes of above spots
        idxs_below  =  np.where(red_zcoors >= first_last_frame[1])[0]                       # indexes of below spots

        red_above        =  np.zeros((len(idxs_above), 3), dtype=np.int32)                  # initialize a matrix to organize label, tags and z-coords for each spots above
        red_above[:, 0]  =  red_tags[list(idxs_above)]                                      # fill the matrix
        red_above[:, 1]  =  red_ints[list(idxs_above)]
        red_above[:, 2]  =  red_zcoors[list(idxs_above)]
        
        red_below        =  np.zeros((len(idxs_below), 3), dtype=np.int32)                  # initialize a matrix to organize label, tags and z-coords for each spots below
        red_below[:, 0]  =  red_tags[list(idxs_below)]                                      # fill the matrix
        red_below[:, 1]  =  red_ints[list(idxs_below)]
        red_below[:, 2]  =  red_zcoors[list(idxs_below)]

        for file in files:    
            if "Green_Overlapping" in file:                                                 # select the .xls file with green overlapping data
                xls_gr_overl_name  =  file
                break

        book_gr_overl         =  open_workbook(analysis_folder + '/' + xls_gr_overl_name)                               # load xls file
        sheet_gr_overl_above  =  book_gr_overl.sheet_by_index(0)
        sheet_gr_overl_below  =  book_gr_overl.sheet_by_index(1)

        int_thr_sm_above  =  sheet_gr_overl_above.col(24)[6].value + 3 * sheet_gr_overl_above.col(24)[7].value          # read average intensity and standard deviation of red sm-spots from the xls (above)
        int_thr_sm_below  =  sheet_gr_overl_below.col(24)[6].value + 3 * sheet_gr_overl_below.col(24)[7].value          # read average intensity and standard deviation of red sm-spots from the xls (below)

        num_bins_above   =  red_above[:, 1].max() / int_thr_sm_above + 1                                                # number of bins to categorize all the above spots (1 sm range, 2 sm range and so on)
        bins_above       =  np.arange(num_bins_above) * int_thr_sm_above                                                # declaring bins
        digit_above      =  np.digitize(red_above[:, 1], bins_above)                                                    # populating bins
        dgt_score_above  =  []

        num_bins_below   =  red_below[:, 1].max() / int_thr_sm_below + 1                                                # number of bins to categorize all the below spots (1 sm range, 2 sm range and so on)
        bins_below       =  np.arange(num_bins_below) * int_thr_sm_below                                                # declaring bins
        digit_below      =  np.digitize(red_below[:, 1], bins_below)                                                    # populating bins
        dgt_score_below  =  []

        for k in range(1, digit_above.max() + 1):
            dgt_score_above_partial  =  []                                                                              # initializing the list to check the overlapping
            idxs                     =  np.where(digit_above == k)[0]                                                   # indexes of tags of the spots in the bin
            for jj in idxs:
                idx_lbl  =  np.where(rgp_red["label"] == red_above[jj, 0])[0][0]                                        # coordinates of the selected spot in the matrix
                bff      =  [spts_green[tuple(rgp_red["coords"][idx_lbl][k_coord])] for k_coord in range(rgp_red["coords"][idx_lbl].shape[0])]     # values taken by the spts_green matrix to check the overlapping
                if np.sum(bff) == 0:
                    dgt_score_above_partial.append(0)                                                                   # if overlapping is not present, append 0 to the partial results list
                else:
                    dgt_score_above_partial.append(1)                                                                   # if overlapping is present, append 1 to the partial results list (its easy to check the percentage of overlapping)
            dgt_score_above.append(dgt_score_above_partial)                                                             # append partial result to global results list


        dgt_rslt_above  =  np.zeros((bins_above.size - 1, 2))                                                           # matrix with the overlapping results
        for aa in range(dgt_rslt_above.shape[0]):
            dgt_rslt_above[aa, :]  =  len(dgt_score_above[aa]), 100 * np.sum(dgt_score_above[aa]) / len(dgt_score_above[aa])    # measure the number of spots and the % of spots overlapped by green


        for k in range(1, digit_below.max() + 1):                                                                       # same as before, but for below
            dgt_score_below_partial  =  []
            idxs                     =  np.where(digit_below == k)[0]
            for jj in idxs:
                idx_lbl  =  np.where(rgp_red["label"] == red_below[jj, 0])[0][0]
                bff      =  [spts_green[tuple(rgp_red["coords"][idx_lbl][k_coord])] for k_coord in range(rgp_red["coords"][idx_lbl].shape[0])]
                if np.sum(bff) == 0:
                    dgt_score_below_partial.append(0)
                else:
                    dgt_score_below_partial.append(1)
            dgt_score_below.append(dgt_score_below_partial)

        dgt_rslt_below  =  np.zeros((bins_below.size - 1, 2))
        for aa in range(dgt_rslt_below.shape[0]):
            dgt_rslt_below[aa, :]  =  len(dgt_score_below[aa]), 100 * np.sum(dgt_score_below[aa]) / len(dgt_score_below[aa])

        book_perc_rslts  =  xlwt.Workbook(encoding='utf-8')
        sheet            =  book_perc_rslts.add_sheet("Sheet1")

        sheet.write(0, 0, "Above")
        sheet.write(1, 0, "min")
        sheet.write(1, 1, "max")
        sheet.write(1, 2, "numb of spots")
        sheet.write(1, 3, "% overlapp")

        sheet.write(0, 5, "Below")
        sheet.write(1, 5, "min")
        sheet.write(1, 6, "max")
        sheet.write(1, 7, "numb of spots")
        sheet.write(1, 8, "% overlapp")

        for bb in range(dgt_rslt_above.shape[0]):
            sheet.write(2 + bb, 0, bins_above[bb])
            sheet.write(2 + bb, 1, bins_above[bb + 1])
            sheet.write(2 + bb, 2, dgt_rslt_above[bb, 0])
            sheet.write(2 + bb, 3, dgt_rslt_above[bb, 1])

        for cc in range(dgt_rslt_below.shape[0]):
            sheet.write(2 + cc, 5, bins_below[cc])
            sheet.write(2 + cc, 6, bins_below[cc + 1])
            sheet.write(2 + cc, 7, dgt_rslt_below[cc, 0])
            sheet.write(2 + cc, 8, dgt_rslt_below[cc, 1])

        book_perc_rslts.save(analysis_folder + xls_gr_overl_name.replace("Green_Overlapping", "Percentage_mRNA_overlapping_by_intensity"))
        print(analysis_folder + xls_gr_overl_name.replace("Green_Overlapping", "Percentage_mRNA_overlapping_by_intensity"))
        self.dgt_rslt_above  =  dgt_rslt_above
        self.dgt_rslt_below  =  dgt_rslt_below
