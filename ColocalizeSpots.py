
"""This function organizes spots in terms of colocalization.

Given the results of an anlysis, this function loads the results and splits
green spots into 2 categories: the ones overlapping red spots and the ones not
overlapping. The same is done with the red spots. Finally it writes results in
an excel file, showing some plots too.
"""

from os import listdir
import datetime
import numpy as np
from skimage.measure import regionprops_table, regionprops
import pyqtgraph as pg
import pyqtgraph.exporters
import xlwt
from xlrd import open_workbook
from scipy.optimize import curve_fit

import RawDataLoader



class ColocalizeSpots:
    """Main class, does all the job"""
    def __init__(self, analysis_folder, raw_data_fname, soft_version, vol_thr, clstr_dist):

        spts_green        =  np.load(analysis_folder + '/spts_lbls_a.npy')
        spts_red          =  np.load(analysis_folder + '/spts_lbls_b.npy')
        first_last_frame  =  np.load(analysis_folder + '/first_last_frame.npy')


        files  =  listdir(analysis_folder)                   # search in the analysis folder the xls files to read the background values
        for file in files:
            if file[-6:] == "_a.xls":                        # select the .xls file with spot analysis results, green
                xls_a_name  =  file
                break
        for file in files:
            if file[-6:] == "_b.xls":                        # # select the .xls file with spot analysis results, red
                xls_b_name  =  file
                break

        book2read_a  =  open_workbook(analysis_folder + "/" + xls_a_name)              # load xls
        sheet_a      =  book2read_a.sheet_by_index(0)
        tags_list_a  =  sheet_a.col_values(0)                                          # copy tags from xls into a python-list
        tags_list_a  =  tags_list_a[3:]                                                # remove first 3 entries - are just titles
        tags_list_a  =  [int(tag[4:]) for tag in tags_list_a]
        bkg_list_a   =  sheet_a.col_values(4)
        bkg_list_a   =  bkg_list_a[3:]                                                 # copy background values in a python-list

        book2read_b  =  open_workbook(analysis_folder + "/" + xls_b_name)
        sheet_b      =  book2read_b.sheet_by_index(0)
        tags_list_b  =  sheet_b.col_values(0)
        tags_list_b  =  tags_list_b[3:]
        tags_list_b  =  [int(tag[4:])for tag in tags_list_b]
        bkg_list_b   =  sheet_b.col_values(4)
        bkg_list_b   =  bkg_list_b[3:]

        chs_spts_nucs  =  np.fromfile(analysis_folder + '/chs_spts_nucs.bin', 'uint16')       # read channels info
        raw_data       =  RawDataLoader.RawDataLoader(raw_data_fname, chs_spts_nucs)          # load raw data

        coloc_mtx              =  np.zeros(spts_green.shape + (3,), dtype=np.int32)
        coloc_mtx[:, :, :, 0]  =  np.sign(spts_green)
        coloc_mtx[:, :, :, 1]  =  np.sign(spts_red)

        # perc of red on green #
        tags_red_ongr  =  np.unique(spts_red * np.sign(spts_green))[1:]                                                                        # identify the label of transcription (b) overlapping green
        rgp_red        =  regionprops_table(spts_red, raw_data.spts_b, properties=["label", "coords", "intensity_image", "area", "centroid"])  # labels and coordinates of all the transcription spots

        rgp_cc   =  regionprops(spts_red)                                                       # regionprops to measure the centroids coordinates (the prvious gives approximation of it)
        rgp_arr  =  []
        for k in range(len(rgp_red["label"])):
            rgp_arr.append([rgp_cc[k]["label"], int(np.round(rgp_cc[k]["centroid"][0])), int(np.round(rgp_cc[k]["centroid"][1])), int(np.round(rgp_cc[k]["centroid"][2]))])    # array with not approximate centroids-coordinates and label

        rgp_arr  =  np.asarray(rgp_arr)                                                                                                                                        # convert to array

        for k in range(rgp_red["label"].size):                                                                                                                                 # substitute the approximated centroids-coordinates with the non approximated ones
            idx2subst                                                                     =  np.where(rgp_arr[:, 0] == rgp_red["label"][k])[0][0]                              # check the indexes of the label
            rgp_red["centroid-0"][k], rgp_red["centroid-1"][k], rgp_red["centroid-2"][k]  =  rgp_arr[idx2subst][1], rgp_arr[idx2subst][2], rgp_arr[idx2subst][3]               # substitute centroids-coordinates with the proper not approximated ones


        mrna_in_transl  =  np.round(tags_red_ongr.size / rgp_red["label"].size, decimals=3)                                                     # ratio of mrna ofverlapped by translational spots (red on green) 

        idxs_red_on_gr    =  np.in1d(rgp_red["label"], tags_red_ongr).nonzero()[0]                                                              # indexes (not labels) of red spots on green spots 
        ints_red_ongreen  =  []
        vols_red_ongreen  =  []
        img_red_ongreen   =  np.zeros(spts_green.shape, dtype=np.int32)

        for rr in idxs_red_on_gr:
            ints_red_ongreen.append(rgp_red["intensity_image"][rr].sum())                                                                  # only for on the overlapping labels, write the intensity values
            vols_red_ongreen.append(rgp_red["area"][rr])                                                                                   # only for on the overlapping labels, write the intensity values
            for rr_in in range(rgp_red["coords"][rr].shape[0]):
                img_red_ongreen[rgp_red["coords"][rr][rr_in][0], rgp_red["coords"][rr][rr_in][1], rgp_red["coords"][rr][rr_in][2]]  =  1   # spts matrix of the b spots overlapping the green 

        hist_ints_red_ongreen  =  np.histogram(ints_red_ongreen, bins=100)
        p1                     =  pg.plot(hist_ints_red_ongreen[1], hist_ints_red_ongreen[0], stepMode=True, pen='r', fillLevel=0, brush=(0, 0, 255, 150), title="Red overlapping on green (Intensity)")
        exporter               =  pg.exporters.ImageExporter(p1.plotItem)
        exporter.export(analysis_folder + '/Hist_Intensity_Red_On_Green.png')
        hist_vols_red_ongreen  =  np.histogram(vols_red_ongreen, bins=100)

        def gauss(x, mu, sigma, A):
            return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

        def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
            return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

        p2           =  pg.plot(hist_vols_red_ongreen[1], hist_vols_red_ongreen[0], stepMode=True, pen='r', fillLevel=0, brush=(0, 255, 255, 150), title="Red overlapping on green (Volume)")

        try:

            params, cov  =  curve_fit(bimodal, hist_vols_red_ongreen[1][:-1], hist_vols_red_ongreen[0], p0=[60, 10, 250, 210, 30, 50], maxfev=3000)
            yy           =  bimodal(hist_vols_red_ongreen[1][:-1], params[0], params[1], params[2], params[3], params[4], params[5])
            myPen        =  pg.mkPen('g', width=2)
            p2.plot(hist_vols_red_ongreen[1][:-1], yy, pen=myPen)
            text1        =  pg.TextItem("mean = " + str(np.round(params[0], decimals=3)))
            text2        =  pg.TextItem("std  = " + str(np.round(params[1], decimals=3)))
            p2.addItem(text1)
            p2.addItem(text2)
            text1.setPos(params[3] + 2 * params[4], params[-1] * 3)
            text2.setPos(params[3] + 2 * params[4], params[-1] * 3 - 10)

        except RuntimeError:
            pass

        exporter               =  pg.exporters.ImageExporter(p2.plotItem)
        exporter.export(analysis_folder + '/Hist_Volume_Red_On_Green.png')
        pg.image(img_red_ongreen + np.sign(spts_red), title="% mRNA in translation")


        # green on red
        tags_green_onred  =  np.unique(spts_green * np.sign(spts_red))[1:]                                                                          # identify the label of translation (a) overlapping red
        rgp_green         =  regionprops_table(spts_green, raw_data.spts_a, properties=["label", "coords", "intensity_image", "area", "centroid"])  # labels and coordinates of all the tranlation spots

        rgp_ff   =  regionprops(spts_green)                                                       # regionprops to measure the centroids coordinates (the prvious gives approximation of it)
        rgp_rra  =  []
        for k in range(len(rgp_green["label"])):
            rgp_rra.append([rgp_ff[k]["label"], int(np.round(rgp_ff[k]["centroid"][0])), int(np.round(rgp_ff[k]["centroid"][1])), int(np.round(rgp_ff[k]["centroid"][2]))])    # array with not approximate centroids-coordinates and label

        rgp_rra  =  np.asarray(rgp_rra)                                                                                                                                        # convert to array

        for k in range(rgp_green["label"].size):                                                                                                                                 # substitute the approximated centroids-coordinates with the non approximated ones
            idx2subst                                                                           =  np.where(rgp_rra[:, 0] == rgp_green["label"][k])[0][0]                        # check the indexes of the label
            rgp_green["centroid-0"][k], rgp_green["centroid-1"][k], rgp_green["centroid-2"][k]  =  rgp_rra[idx2subst][1], rgp_rra[idx2subst][2], rgp_rra[idx2subst][3]           # substitute centroids-coordinates with the proper not approximated ones

        green_onmrna      =  np.round(tags_green_onred.size / rgp_green["label"].size, decimals=3)                                     # ratio of mrna ofverlapped by translational spots (red on green)

        idxs_green_onred  =  np.in1d(rgp_green["label"], tags_green_onred).nonzero()[0]                                                           # indexes (not labels) of red spots on green spots 
        ints_green_onred  =  []
        vols_green_onred  =  []
        img_green_onred   =  np.zeros(spts_red.shape, dtype=np.int32)

        for gg in idxs_green_onred:
            ints_green_onred.append(rgp_green["intensity_image"][gg].sum())                                                       # only for on the overlapping labels, write the intensity values
            vols_green_onred.append(rgp_green["area"][gg])                                                                        # only for on the overlapping labels, write the intensity values
            for gg_in in range(rgp_green["coords"][gg].shape[0]):
                img_green_onred[rgp_green["coords"][gg][gg_in][0], rgp_green["coords"][gg][gg_in][1], rgp_green["coords"][gg][gg_in][2]]  =  1   # spts matrix of the b spots overlapping the green 


        hist_ints_green_onred  =  np.histogram(ints_green_onred, bins=100)
        p3                     =  pg.plot(hist_ints_green_onred[1], hist_ints_green_onred[0], stepMode=True, pen='g', fillLevel=0, brush=(0, 0, 255, 150), title="Green overlapping on red (Intensity)")
        exporter               =  pg.exporters.ImageExporter(p3.plotItem)
        exporter.export(analysis_folder + '/Hist_Intensity_Green_On_Red.png')
        hist_vols_green_onred  =  np.histogram(vols_green_onred, bins=100)
        p4                     =  pg.plot(hist_vols_green_onred[1], hist_vols_green_onred[0], stepMode=True, pen='g', fillLevel=0, brush=(0, 255, 255, 150), title="Green overlapping on red (Volume)")
        exporter               =  pg.exporters.ImageExporter(p4.plotItem)
        exporter.export(analysis_folder + '/Hist_Volume_Green_On_Red.png')
        pg.image(img_green_onred + np.sign(spts_green), title="% Protein in transcription")

        # red not overlapping
        idxs_red_not_overlp  =  np.logical_not(np.in1d(rgp_red["label"], tags_red_ongr)).nonzero()[0]                                           # indexes of red spots not overlapping green
        ints_red_not_overlp  =  []
        vols_red_not_overlp  =  []
        img_red_not_overlp   =  np.zeros(spts_red.shape, dtype=np.int32)

        for nr in idxs_red_not_overlp:
            ints_red_not_overlp.append(rgp_red["intensity_image"][nr].sum())                                            # only for on the not-overlapping labels, write the intensity values
            vols_red_not_overlp.append(rgp_red["area"][nr])                                                             # only for on the not-overlapping labels, write the intensity values
            for nr_in in range(rgp_red["coords"][nr].shape[0]):
                img_red_not_overlp[rgp_red["coords"][nr][nr_in][0], rgp_red["coords"][nr][nr_in][1], rgp_red["coords"][nr][nr_in][2]]  =  1        # spts matrix of the b spots not overlapping the green 

        red_not_overlp_clrs              =  np.zeros(spts_red.shape + (3,), dtype=np.int32)
        red_not_overlp_clrs[:, :, :, 0]  =  np.sign(img_red_not_overlp)
        red_not_overlp_clrs[:, :, :, 1]  =  np.sign(spts_green)
        pg.image(red_not_overlp_clrs, title="red NON overlapping")

        hist_ints_red_not_ovrl   =  np.histogram(ints_red_not_overlp, bins=100)
        p5                       =  pg.plot(hist_ints_red_not_ovrl[1], hist_ints_red_not_ovrl[0], stepMode=True, pen='r', fillLevel=0, brush=(0, 0, 255, 150), title="Red NOT overlapping (Intensity)")
        exporter                 =  pg.exporters.ImageExporter(p5.plotItem)
        exporter.export(analysis_folder + '/Hist_Intensity_Red_Not_Overl.png')
        hist_vols_red_not_ovrl   =  np.histogram(vols_red_not_overlp, bins=100)
        p6                       =  pg.plot(hist_vols_red_not_ovrl[1], hist_vols_red_not_ovrl[0], stepMode=True, pen='r', fillLevel=0, brush=(0, 255, 255, 150), title="Red NOT overlapping (Volume)")
        exporter                 =  pg.exporters.ImageExporter(p6.plotItem)
        exporter.export(analysis_folder + '/Hist_Volume_Red_Not_Overl.png')

        # green not overlapping
        idxs_green_not_overlp  =  np.logical_not(np.in1d(rgp_green["label"], tags_green_onred)).nonzero()[0]                                       # indexes of green spots not overlapping green
        ints_green_not_overlp  =  []
        vols_green_not_overlp  =  []
        img_green_not_overlp   =  np.zeros(spts_green.shape, dtype=np.int32)

        for ng in idxs_green_not_overlp:
            ints_green_not_overlp.append(rgp_green["intensity_image"][ng].sum())                                # only for on the not-overlapping labels, write the intensity values
            vols_green_not_overlp.append(rgp_green["area"][ng])                                                 # only for on the not-overlapping labels, write the intensity values
            for ng_in in range(rgp_green["coords"][ng].shape[0]):
                img_green_not_overlp[rgp_green["coords"][ng][ng_in][0], rgp_green["coords"][ng][ng_in][1], rgp_green["coords"][ng][ng_in][2]]  =  1   # matrix of the a spots not overlapping the green

        green_not_overlp_clrs              =  np.zeros(spts_green.shape + (3,), dtype=np.int32)
        green_not_overlp_clrs[:, :, :, 0]  =  np.sign(img_green_not_overlp)
        green_not_overlp_clrs[:, :, :, 1]  =  np.sign(spts_red)
        pg.image(green_not_overlp_clrs, title="green NON overlapping")

        hist_ints_green_not_ovrl  =  np.histogram(ints_green_not_overlp, bins=100)
        p7                        =  pg.plot(hist_ints_green_not_ovrl[1], hist_ints_green_not_ovrl[0], stepMode=True, pen='g', fillLevel=0, brush=(0, 0, 255, 150), title="Green NOT overlapping (Intensity)")
        exporter                  =  pg.exporters.ImageExporter(p7.plotItem)
        exporter.export(analysis_folder + '/Hist_Intensity_Green_Not_Overl.png')
        hist_green_not_ovrl_vols  =  np.histogram(vols_green_not_overlp, bins=100)
        p8                        =  pg.plot(hist_green_not_ovrl_vols[1], hist_green_not_ovrl_vols[0], stepMode=True, pen='g', fillLevel=0, brush=(0, 255, 255, 150), title="Green NOT overlapping (Volume)")
        exporter                  =  pg.exporters.ImageExporter(p8.plotItem)
        exporter.export(analysis_folder + '/Hist_Volume_Green_Not_Overl.png')

        book_green_not_overlapping    =  xlwt.Workbook(encoding='utf-8')
        sheet1_green_not_overlapping  =  book_green_not_overlapping.add_sheet("Above")
        sheet2_green_not_overlapping  =  book_green_not_overlapping.add_sheet("Below")
        sheet1_green_not_overlapping.write(0, 0, "Spt_id")
        sheet1_green_not_overlapping.write(0, 1, "Ints")
        sheet1_green_not_overlapping.write(0, 2, "bkg")
        sheet1_green_not_overlapping.write(0, 3, "Ints - bkg")
        sheet1_green_not_overlapping.write(0, 4, "Ints / bkg")
        sheet1_green_not_overlapping.write(0, 5, "Vols")
        sheet1_green_not_overlapping.write(0, 6, "z-coord")
        sheet1_green_not_overlapping.write(0, 7, "x-coord")
        sheet1_green_not_overlapping.write(0, 8, "y-coord")
        sheet2_green_not_overlapping.write(0, 0, "Spt_id")
        sheet2_green_not_overlapping.write(0, 1, "Ints")
        sheet2_green_not_overlapping.write(0, 2, "bkg")
        sheet2_green_not_overlapping.write(0, 3, "Ints - bkg")
        sheet2_green_not_overlapping.write(0, 4, "Ints / bkg")
        sheet2_green_not_overlapping.write(0, 5, "Vols")
        sheet2_green_not_overlapping.write(0, 6, "z-coord")
        sheet2_green_not_overlapping.write(0, 7, "x-coord")
        sheet2_green_not_overlapping.write(0, 8, "y-coord")

        book    =  xlwt.Workbook(encoding='utf-8')
        sheet1  =  book.add_sheet("Red Overlapping")
        sheet2  =  book.add_sheet("Green overlapping")
        sheet3  =  book.add_sheet("Red NOT Overlapping")
        sheet4  =  book.add_sheet("Green NOT overlapping")

        sheet1.write(0, 0, "Spt_id")
        sheet1.write(0, 1, "Ints")
        sheet1.write(0, 2, "bkg")
        sheet1.write(0, 3, "Ints - bkg")
        sheet1.write(0, 4, "Ints / bkg")
        sheet1.write(0, 5, "Vols")
        sheet1.write(0, 6, "z-coord")
        sheet1.write(0, 7, "x-coord")
        sheet1.write(0, 8, "y-coord")
        sheet1.write(0, 10, "% of mrna in transl")
        sheet1.write(1, 10, mrna_in_transl)
        sheet1.write(3, 10, "% of transl on mrna")
        sheet1.write(4, 10, green_onmrna)
        sheet1.write(0, 14, "date")
        sheet1.write(1, 14, datetime.datetime.now().strftime("%d-%b-%Y"))
        sheet1.write(0, 15, "soft_version")
        sheet1.write(1, 15, soft_version)

        sheet2.write(0, 0, "Spt_id")
        sheet2.write(0, 1, "Ints")
        sheet2.write(0, 2, "bkg")
        sheet2.write(0, 3, "Ints - bkg")
        sheet2.write(0, 4, "Ints / bkg")
        sheet2.write(0, 5, "Vols")
        sheet2.write(0, 6, "z-coord")
        sheet2.write(0, 7, "x-coord")
        sheet2.write(0, 8, "y-coord")

        sheet3.write(0, 0, "Spt_id")
        sheet3.write(0, 1, "Ints")
        sheet3.write(0, 2, "bkg")
        sheet3.write(0, 3, "Ints - bkg")
        sheet3.write(0, 4, "Ints / bkg")
        sheet3.write(0, 5, "Vols")
        sheet3.write(0, 6, "z-coord")
        sheet3.write(0, 7, "x-coord")
        sheet3.write(0, 8, "y-coord")

        sheet4.write(0, 0, "Spt_id")
        sheet4.write(0, 1, "Ints")
        sheet4.write(0, 2, "bkg")
        sheet4.write(0, 3, "Ints - bkg")
        sheet4.write(0, 4, "Ints / bkg")
        sheet4.write(0, 5, "Vols")
        sheet4.write(0, 6, "z-coord")
        sheet4.write(0, 7, "x-coord")
        sheet4.write(0, 8, "y-coord")

        sht_idxs  =  1
        for k_red in idxs_red_on_gr:
            bkg      =  bkg_list_b[tags_list_b.index(rgp_red["label"][k_red])]               # select background values ridden in the .xls file from the list
            spt_int  =  np.sum(rgp_red["intensity_image"][k_red])                            # define the intensity variable instead of calculating each time

            sheet1.write(sht_idxs, 0, int(rgp_red["label"][k_red]))
            sheet1.write(sht_idxs, 1, int(np.sum(rgp_red["intensity_image"][k_red])))
            sheet1.write(sht_idxs, 2, bkg)
            sheet1.write(sht_idxs, 3, spt_int - bkg * rgp_red["area"][k_red])
            sheet1.write(sht_idxs, 4, spt_int / bkg)
            sheet1.write(sht_idxs, 5, int(rgp_red["area"][k_red]))
            sheet1.write(sht_idxs, 6, int(rgp_red["centroid-0"][k_red]))
            sheet1.write(sht_idxs, 7, int(rgp_red["centroid-1"][k_red]))
            sheet1.write(sht_idxs, 8, int(rgp_red["centroid-2"][k_red]))
            sht_idxs  +=  1

        sht_idxs  =  1
        for k_gr in idxs_green_onred:
            bkg      =  bkg_list_a[tags_list_a.index(rgp_green["label"][k_gr])]
            spt_int  =  np.sum(rgp_green["intensity_image"][k_gr])

            sheet2.write(sht_idxs, 0, int(rgp_green["label"][k_gr]))
            sheet2.write(sht_idxs, 1, int(spt_int))
            sheet2.write(sht_idxs, 2, bkg)
            sheet2.write(sht_idxs, 3, spt_int - bkg * rgp_green["area"][k_gr])
            sheet2.write(sht_idxs, 4, spt_int / bkg)
            sheet2.write(sht_idxs, 5, int(rgp_green["area"][k_gr]))
            sheet2.write(sht_idxs, 6, int(rgp_green["centroid-0"][k_gr]))
            sheet2.write(sht_idxs, 7, int(rgp_green["centroid-1"][k_gr]))
            sheet2.write(sht_idxs, 8, int(rgp_green["centroid-2"][k_gr]))
            sht_idxs  +=  1

        sht_idxs  =  1
        for k_red2 in idxs_red_not_overlp:
            bkg      =  bkg_list_b[tags_list_b.index(rgp_red["label"][k_red2])]
            spt_int  =  np.sum(rgp_red["intensity_image"][k_red2])

            sheet3.write(sht_idxs, 0, int(rgp_red["label"][k_red2]))
            sheet3.write(sht_idxs, 1, int(np.sum(rgp_red["intensity_image"][k_red2])))
            sheet3.write(sht_idxs, 2, bkg)
            sheet3.write(sht_idxs, 3, spt_int - bkg * rgp_red["area"][k_red2])
            sheet3.write(sht_idxs, 4, spt_int / bkg)
            sheet3.write(sht_idxs, 5, int(rgp_red["area"][k_red2]))
            sheet3.write(sht_idxs, 6, int(rgp_red["centroid-0"][k_red2]))
            sheet3.write(sht_idxs, 7, int(rgp_red["centroid-1"][k_red2]))
            sheet3.write(sht_idxs, 8, int(rgp_red["centroid-2"][k_red2]))
            sht_idxs  +=  1

        sht_idxs   =  1
        sht_idxs1  =  1
        sht_idxs2  =  1

        for k_gr2 in idxs_green_not_overlp:
            bkg      =  bkg_list_a[tags_list_a.index(rgp_green["label"][k_gr2])]
            spt_int  =  np.sum(rgp_green["intensity_image"][k_gr2])

            sheet4.write(sht_idxs, 0, int(rgp_green["label"][k_gr2]))
            sheet4.write(sht_idxs, 1, int(np.sum(rgp_green["intensity_image"][k_gr2])))
            sheet4.write(sht_idxs, 2, bkg)
            sheet4.write(sht_idxs, 3, spt_int - bkg * rgp_green["area"][k_gr2])
            sheet4.write(sht_idxs, 4, spt_int / bkg)
            sheet4.write(sht_idxs, 5, int(rgp_green["area"][k_gr2]))
            sheet4.write(sht_idxs, 6, int(rgp_green["centroid-0"][k_gr2]))
            sheet4.write(sht_idxs, 7, int(rgp_green["centroid-1"][k_gr2]))
            sheet4.write(sht_idxs, 8, int(rgp_green["centroid-2"][k_gr2]))

            if int(rgp_green["centroid-0"][k_gr2]) <= first_last_frame[0]:
                sheet1_green_not_overlapping.write(sht_idxs1, 0, int(rgp_green["label"][k_gr2]))
                sheet1_green_not_overlapping.write(sht_idxs1, 1, int(np.sum(rgp_green["intensity_image"][k_gr2])))
                sheet1_green_not_overlapping.write(sht_idxs1, 2, bkg)
                sheet1_green_not_overlapping.write(sht_idxs1, 3, spt_int - bkg * rgp_green["area"][k_gr2])
                sheet1_green_not_overlapping.write(sht_idxs1, 4, spt_int / bkg)
                sheet1_green_not_overlapping.write(sht_idxs1, 5, int(rgp_green["area"][k_gr2]))
                sheet1_green_not_overlapping.write(sht_idxs1, 6, int(rgp_green["centroid-0"][k_gr2]))
                sheet1_green_not_overlapping.write(sht_idxs1, 7, int(rgp_green["centroid-1"][k_gr2]))
                sheet1_green_not_overlapping.write(sht_idxs1, 8, int(rgp_green["centroid-2"][k_gr2]))
                sht_idxs1  +=  1
            else:
                sheet2_green_not_overlapping.write(sht_idxs2, 0, int(rgp_green["label"][k_gr2]))
                sheet2_green_not_overlapping.write(sht_idxs2, 1, int(np.sum(rgp_green["intensity_image"][k_gr2])))
                sheet2_green_not_overlapping.write(sht_idxs2, 2, bkg)
                sheet2_green_not_overlapping.write(sht_idxs2, 3, spt_int - bkg * rgp_green["area"][k_gr2])
                sheet2_green_not_overlapping.write(sht_idxs2, 4, spt_int / bkg)
                sheet2_green_not_overlapping.write(sht_idxs2, 5, int(rgp_green["area"][k_gr2]))
                sheet2_green_not_overlapping.write(sht_idxs2, 6, int(rgp_green["centroid-0"][k_gr2]))
                sheet2_green_not_overlapping.write(sht_idxs2, 7, int(rgp_green["centroid-1"][k_gr2]))
                sheet2_green_not_overlapping.write(sht_idxs2, 8, int(rgp_green["centroid-2"][k_gr2]))
                sht_idxs2  +=  1

            sht_idxs  +=  1

        book.save(analysis_folder + "/" + xls_a_name[:-22] + "_Overlapping_results_volThr" + str(vol_thr) + "_clstrDist_" + str(clstr_dist) + ".xls")
        book_green_not_overlapping.save(analysis_folder + "/" + xls_a_name[:-22] + "_Green_Not_Overlapping_volThr" + str(vol_thr) + "_clstrDist_" + str(clstr_dist) + ".xls")

        ints_green_not_overl_above      =  []
        vols_green_not_overl_above      =  []
        ints_green_not_overl_below      =  []
        vols_green_not_overl_below      =  []
        ints_green_not_overl_abovembkg  =  []
        ints_green_not_overl_belowmbkg  =  []
        for k_gr2 in idxs_green_not_overlp:
            bkg      =  bkg_list_a[tags_list_a.index(rgp_green["label"][k_gr2])]
            spt_int  =  np.sum(rgp_green["intensity_image"][k_gr2])

            if int(rgp_green["centroid-0"][k_gr2]) <= first_last_frame[0]:
                ints_green_not_overl_above.append(int(np.sum(rgp_green["intensity_image"][k_gr2])))
                ints_green_not_overl_abovembkg.append(spt_int - bkg * rgp_green["area"][k_gr2])
                vols_green_not_overl_above.append(int(rgp_green["area"][k_gr2]))
            else:
                ints_green_not_overl_below.append(int(np.sum(rgp_green["intensity_image"][k_gr2])))
                ints_green_not_overl_belowmbkg.append(spt_int - bkg * rgp_green["area"][k_gr2])
                vols_green_not_overl_below.append(int(rgp_green["area"][k_gr2]))


        mtx_pair_above               =  []
        mtx_pair_below               =  []
        multi_over_green_tags_above  =  []
        multi_over_green_tags_below  =  []

        for kk_gr in idxs_green_onred:
            # print(kk_red)
            if int(rgp_green["centroid-0"][kk_gr]) <= first_last_frame[0]:

                line_list  =  []
                red_tag    =  []
                for jj in range(rgp_green["coords"][kk_gr].shape[0]):
                    red_tag.append(spts_red[rgp_green["coords"][kk_gr][jj][0], rgp_green["coords"][kk_gr][jj][1], rgp_green["coords"][kk_gr][jj][2]])
                red_tag  =  np.trim_zeros(np.unique(red_tag))
                if red_tag.size == 1:
                    red_idx  =  np.where(rgp_red["label"] == red_tag)[0][0]

                    bkg_g      =  bkg_list_a[tags_list_a.index(rgp_green["label"][kk_gr])]
                    spt_int_g  =  np.sum(rgp_green["intensity_image"][kk_gr])
                    bkg_r      =  bkg_list_b[tags_list_b.index(red_tag)]
                    spt_int_r  =  np.sum(rgp_red["intensity_image"][red_idx])

                    line_list.append([int(rgp_green["label"][kk_gr]), int(spt_int_g), bkg_g, spt_int_g - bkg_g * rgp_green["area"][kk_gr], spt_int_g / bkg_g, int(rgp_green["area"][kk_gr]), int(rgp_green["centroid-0"][kk_gr]), int(rgp_green["centroid-1"][kk_gr]), int(rgp_green["centroid-2"][kk_gr]), int(red_tag), int(spt_int_r), bkg_r, spt_int_r - bkg_r * rgp_red["area"][red_idx], spt_int_r / bkg_r, int(rgp_red["area"][red_idx]), int(rgp_red["centroid-0"][red_idx]), int(rgp_red["centroid-1"][red_idx]), int(rgp_red["centroid-2"][red_idx])])

                    mtx_pair_above.append(line_list)
                else:
                    multi_over_green_tags_above.append(kk_gr)

            if int(rgp_green["centroid-0"][kk_gr]) >= first_last_frame[1]:

                line_list  =  []
                red_tag    =  []
                for jj in range(rgp_green["coords"][kk_gr].shape[0]):
                    red_tag.append(spts_red[rgp_green["coords"][kk_gr][jj][0], rgp_green["coords"][kk_gr][jj][1], rgp_green["coords"][kk_gr][jj][2]])
                red_tag  =  np.trim_zeros(np.unique(red_tag))
                if red_tag.size == 1:
                    red_idx  =  np.where(rgp_red["label"] == red_tag)[0][0]

                    bkg_g      =  bkg_list_a[tags_list_a.index(rgp_green["label"][kk_gr])]
                    spt_int_g  =  np.sum(rgp_green["intensity_image"][kk_gr])
                    bkg_r      =  bkg_list_b[tags_list_b.index(red_tag)]
                    spt_int_r  =  np.sum(rgp_red["intensity_image"][red_idx])

                    line_list.append([int(rgp_green["label"][kk_gr]), int(spt_int_g), bkg_g, spt_int_g - bkg_g * rgp_green["area"][kk_gr], spt_int_g / bkg_g, int(rgp_green["area"][kk_gr]), int(rgp_green["centroid-0"][kk_gr]), int(rgp_green["centroid-1"][kk_gr]), int(rgp_green["centroid-2"][kk_gr]), int(red_tag), int(spt_int_r), bkg_r, spt_int_r - bkg_r * rgp_red["area"][red_idx], spt_int_r / bkg_r, int(rgp_red["area"][red_idx]), int(rgp_red["centroid-0"][red_idx]), int(rgp_red["centroid-1"][red_idx]), int(rgp_red["centroid-2"][red_idx])])

                    mtx_pair_below.append(line_list)
                else:
                    multi_over_green_tags_below.append(kk_gr)

        mtx_pair_above  =  np.squeeze(np.asarray(mtx_pair_above))
        mtx_pair_below  =  np.squeeze(np.asarray(mtx_pair_below))


        # red_ints = []
        # for k in range(rgp_red["area"].size):
        #     red_ints.append(rgp_red["intensity_image"][k].sum())

        # hist_red_ints         =  np.histogram(red_ints, bins=200)
        # params_red_ints, cov  =  curve_fit(bimodal, hist_red_ints[1][:-1], hist_red_ints[0], p0=[hist_red_ints[1][np.argmax(hist_red_ints[0])], 1000, hist_red_ints[0].max(), 40000, 30, 100], maxfev=3000)
        # vv_red_ints           =  bimodal(hist_red_ints[1][1:], params_red_ints[0], params_red_ints[1], params_red_ints[2], params_red_ints[3], params_red_ints[4], params_red_ints[5])
        # w                     =  pg.plot(hist_red_ints[1], hist_red_ints[0], stepMode=True, title="red_ints")
        # w.plot(hist_red_ints[1][1:], vv_red_ints, pen='r')

        # hist_red_vols  =  np.histogram(rgp_red["area"], bins=200)
        # params_red_vols, cov  =  curve_fit(bimodal, hist_red_vols[1][:-1], hist_red_vols[0], p0=[hist_red_vols[1][np.argmax(hist_red_vols[0])], 10, hist_red_vols[0].max(), 120, 30, 50], maxfev=3000)
        # vv_red_vols           =  bimodal(hist_red_vols[1][1:], params_red_vols[0], params_red_vols[1], params_red_vols[2], params_red_vols[3], params_red_vols[4], params_red_vols[5])
        # w                     =  pg.plot(hist_red_vols[1], hist_red_vols[0], stepMode=True, title="red_vols")
        # w.plot(hist_red_vols[1][1:], vv_red_vols, pen='r')

        red_ints_above  =  []
        red_ints_below  =  []
        red_vols_above  =  []
        red_vols_below  =  []
        for kk in range(rgp_red["area"].size):
            if rgp_red["centroid-0"][kk] <= first_last_frame[0]:
                red_ints_above.append(rgp_red["intensity_image"][kk].sum())
                red_vols_above.append(rgp_red["area"][kk])
            else:
                red_ints_below.append(rgp_red["intensity_image"][kk].sum())
                red_vols_below.append(rgp_red["area"][kk])

        # green_ints_above  =  []
        # green_ints_below  =  []
        # green_vols_above  =  []
        # green_vols_below  =  []
        # for kk in range(rgp_green["area"].size):
        #     if rgp_green["centroid-0"][kk] <= first_last_frame[0]:
        #         green_ints_above.append(rgp_green["intensity_image"][kk].sum())
        #         green_vols_above.append(rgp_green["area"][kk])
        #     else:
        #         green_ints_below.append(rgp_green["intensity_image"][kk].sum())
        #         green_vols_below.append(rgp_green["area"][kk])

        hist_red_ints_above         =  np.histogram(red_ints_above, bins=100)
        params_red_ints_above, cov  =  curve_fit(gauss, hist_red_ints_above[1][:-1], hist_red_ints_above[0], p0=[hist_red_ints_above[1][np.argmax(hist_red_ints_above[0])], hist_red_ints_above[1][np.argmax(hist_red_ints_above[0])] / 10, hist_red_ints_above[0].max()], maxfev=3000)
        vv_red_ints_above           =  gauss(hist_red_ints_above[1][1:], params_red_ints_above[0], params_red_ints_above[1], params_red_ints_above[2])
        w                           =  pg.plot(hist_red_ints_above[1], hist_red_ints_above[0], stepMode=True, title="red_ints_above")
        w.plot(hist_red_ints_above[1][1:], vv_red_ints_above, pen='r')

        # params_red_ints_above_bimod, cov  =  curve_fit(bimodal, hist_red_ints_above[1][:-1], hist_red_ints_above[0], p0=[hist_red_ints_above[1][np.argmax(hist_red_ints_above[0])], hist_red_ints_above[1][np.argmax(hist_red_ints_above[0])] / 2, hist_red_ints_above[0].max(), 40000, 10000, 25], maxfev=3000)
        # vv_red_ints_above_dimob           =  bimodal(hist_red_ints_above[1][1:], params_red_ints_above_bimod[0], params_red_ints_above_bimod[1], params_red_ints_above_bimod[2], params_red_ints_above_bimod[3], params_red_ints_above_bimod[4], params_red_ints_above_bimod[5])
        # w                           =  pg.plot(hist_red_ints_above[1], hist_red_ints_above[0], stepMode=True, title="red_ints_above, double gaussian")
        # w.plot(hist_red_ints_above[1][1:], vv_red_ints_above_dimob, pen='r')

        hist_red_ints_below         =  np.histogram(red_ints_below, bins=100)
        params_red_ints_below, cov  =  curve_fit(gauss, hist_red_ints_below[1][:-1], hist_red_ints_below[0], p0=[hist_red_ints_below[1][np.argmax(hist_red_ints_below[0])], hist_red_ints_below[1][np.argmax(hist_red_ints_below[0])] / 10, hist_red_ints_below[0].max()], maxfev=3000)
        vv_red_ints_below           =  gauss(hist_red_ints_below[1][1:], params_red_ints_below[0], params_red_ints_below[1], params_red_ints_below[2])
        w                           =  pg.plot(hist_red_ints_below[1], hist_red_ints_below[0], stepMode=True, title="red_ints_below")
        w.plot(hist_red_ints_below[1][1:], vv_red_ints_below, pen='r')

        # params_red_ints_below_bimo, cov  =  curve_fit(bimodal, hist_red_ints_below[1][:-1], hist_red_ints_below[0], p0=[hist_red_ints_below[1][np.argmax(hist_red_ints_below[0])], hist_red_ints_below[1][np.argmax(hist_red_ints_below[0])] / 2, hist_red_ints_below[0].max(), 40000, 1000, 80], maxfev=3000)
        # vv_red_ints_below_bimo           =  bimodal(hist_red_ints_below[1][1:], params_red_ints_below_bimo[0], params_red_ints_below_bimo[1], params_red_ints_below_bimo[2], params_red_ints_below_bimo[3], params_red_ints_below_bimo[4], params_red_ints_below_bimo[5])
        # w                           =  pg.plot(hist_red_ints_below[1], hist_red_ints_below[0], stepMode=True, title="red_ints_below")
        # w.plot(hist_red_ints_below[1][1:], vv_red_ints_below_bimo, pen='r')


        hist_red_vols_above         =  np.histogram(red_vols_above, bins=100)
        params_red_vols_above, cov  =  curve_fit(gauss, hist_red_vols_above[1][:-1], hist_red_vols_above[0], p0=[hist_red_vols_above[1][np.argmax(hist_red_vols_above[0])], hist_red_vols_above[1][np.argmax(hist_red_vols_above[0])] / 10, hist_red_vols_above[0].max()], maxfev=3000000)
        vv_red_vols_above           =  gauss(hist_red_vols_above[1][1:], params_red_vols_above[0], params_red_vols_above[1], params_red_vols_above[2])
        w                           =  pg.plot(hist_red_vols_above[1], hist_red_vols_above[0], stepMode=True, title="red_vols_above")
        w.plot(hist_red_vols_above[1][1:], vv_red_vols_above, pen='r')

        hist_red_vols_below         =  np.histogram(red_vols_below, bins=100)
        params_red_vols_below, cov  =  curve_fit(gauss, hist_red_vols_below[1][:-1], hist_red_vols_below[0], p0=[hist_red_vols_below[1][np.argmax(hist_red_vols_below[0])], hist_red_vols_below[1][np.argmax(hist_red_vols_below[0])] / 10, hist_red_vols_below[0].max()], maxfev=3000)
        vv_red_vols_below           =  gauss(hist_red_vols_below[1][1:], params_red_vols_below[0], params_red_vols_below[1], params_red_vols_below[2])
        w                           =  pg.plot(hist_red_vols_below[1], hist_red_vols_below[0], stepMode=True, title="red_vols_below")
        w.plot(hist_red_vols_below[1][1:], vv_red_vols_below, pen='r')

        red_ints_tot  =  []
        red_vols_tot  =  []
        for kk in range(rgp_red["area"].size):
            red_ints_tot.append(rgp_red["intensity_image"][kk].sum())
            red_vols_tot.append(rgp_red["area"][kk])

        hist_red_ints_tot         =  np.histogram(red_ints_tot, bins=100)
        params_red_ints_tot, cov  =  curve_fit(gauss, hist_red_ints_tot[1][:-1], hist_red_ints_tot[0], p0=[hist_red_ints_tot[1][np.argmax(hist_red_ints_below[0])], hist_red_ints_tot[1][np.argmax(hist_red_ints_tot[0])] / 10, hist_red_ints_tot[0].max()], maxfev=3000)
        vv_red_ints_tot           =  gauss(hist_red_ints_tot[1][1:], params_red_ints_tot[0], params_red_ints_tot[1], params_red_ints_tot[2])
        w                           =  pg.plot(hist_red_ints_tot[1], hist_red_ints_tot[0], stepMode=True, title="red_ints_tot")
        w.plot(hist_red_ints_tot[1][1:], vv_red_ints_tot, pen='r')

        hist_red_vols_tot         =  np.histogram(red_vols_tot, bins=100)
        params_red_vols_tot, cov  =  curve_fit(gauss, hist_red_vols_tot[1][:-1], hist_red_vols_tot[0], p0=[hist_red_vols_tot[1][np.argmax(hist_red_vols_tot[0])], hist_red_vols_tot[1][np.argmax(hist_red_vols_tot[0])] / 10, hist_red_vols_tot[0].max()], maxfev=3000)
        vv_red_vols_tot           =  gauss(hist_red_vols_tot[1][1:], params_red_vols_tot[0], params_red_vols_tot[1], params_red_vols_tot[2])
        w                         =  pg.plot(hist_red_vols_tot[1], hist_red_vols_tot[0], stepMode=True, title="red_vols_tot")
        w.plot(hist_red_vols_tot[1][1:], vv_red_vols_tot, pen='r')



        book_overl_pair  =  xlwt.Workbook(encoding='utf-8')
        sheet1_pair      =  book_overl_pair.add_sheet("Above")
        sheet2_pair      =  book_overl_pair.add_sheet("Below")
        sheet3_pair      =  book_overl_pair.add_sheet("Multiple Above")
        sheet4_pair      =  book_overl_pair.add_sheet("Multiple Below")

        sheet1_pair.write(0, 0, "Green")
        sheet1_pair.write(1, 0, "Spt_id")
        sheet1_pair.write(1, 1, "Ints")
        sheet1_pair.write(1, 2, "bkg")
        sheet1_pair.write(1, 3, "Ints - bkg")
        sheet1_pair.write(1, 4, "Ints / bkg")
        sheet1_pair.write(1, 5, "Vols")
        sheet1_pair.write(1, 6, "z-coord")
        sheet1_pair.write(1, 7, "x-coord")
        sheet1_pair.write(1, 8, "y-coord")
        sheet1_pair.write(0, 10, "RED")
        sheet1_pair.write(1, 10, "Spt_id")
        sheet1_pair.write(1, 11, "Ints")
        sheet1_pair.write(1, 12, "bkg")
        sheet1_pair.write(1, 13, "Ints - bkg")
        sheet1_pair.write(1, 14, "Ints / bkg")
        sheet1_pair.write(1, 15, "Vols")
        sheet1_pair.write(1, 16, "z-coord")
        sheet1_pair.write(1, 17, "x-coord")
        sheet1_pair.write(1, 18, "y-coord")

        sheet2_pair.write(0, 0, "GREEN")
        sheet2_pair.write(1, 0, "Spt_id")
        sheet2_pair.write(1, 1, "Ints")
        sheet2_pair.write(1, 2, "bkg")
        sheet2_pair.write(1, 3, "Ints - bkg")
        sheet2_pair.write(1, 4, "Ints / bkg")
        sheet2_pair.write(1, 5, "Vols")
        sheet2_pair.write(1, 6, "z-coord")
        sheet2_pair.write(1, 7, "x-coord")
        sheet2_pair.write(1, 8, "y-coord")
        sheet2_pair.write(0, 10, "RED")
        sheet2_pair.write(1, 10, "Spt_id")
        sheet2_pair.write(1, 11, "Ints")
        sheet2_pair.write(1, 12, "bkg")
        sheet2_pair.write(1, 13, "Ints - bkg")
        sheet2_pair.write(1, 14, "Ints / bkg")
        sheet2_pair.write(1, 15, "Vols")
        sheet2_pair.write(1, 16, "z-coord")
        sheet2_pair.write(1, 17, "x-coord")
        sheet2_pair.write(1, 18, "y-coord")

        print(mtx_pair_above.shape)
        for jjj_a in range(mtx_pair_above.shape[0]):
            sheet1_pair.write(2 + jjj_a, 0, mtx_pair_above[jjj_a, 0])
            sheet1_pair.write(2 + jjj_a, 1, mtx_pair_above[jjj_a, 1])
            sheet1_pair.write(2 + jjj_a, 2, mtx_pair_above[jjj_a, 2])
            sheet1_pair.write(2 + jjj_a, 3, mtx_pair_above[jjj_a, 3])
            sheet1_pair.write(2 + jjj_a, 4, mtx_pair_above[jjj_a, 4])
            sheet1_pair.write(2 + jjj_a, 5, mtx_pair_above[jjj_a, 5])
            sheet1_pair.write(2 + jjj_a, 6, mtx_pair_above[jjj_a, 6])
            sheet1_pair.write(2 + jjj_a, 7, mtx_pair_above[jjj_a, 7])
            sheet1_pair.write(2 + jjj_a, 8, mtx_pair_above[jjj_a, 8])
            sheet1_pair.write(2 + jjj_a, 10, mtx_pair_above[jjj_a, 9])
            sheet1_pair.write(2 + jjj_a, 11, mtx_pair_above[jjj_a, 10])
            sheet1_pair.write(2 + jjj_a, 12, mtx_pair_above[jjj_a, 11])
            sheet1_pair.write(2 + jjj_a, 13, mtx_pair_above[jjj_a, 12])
            sheet1_pair.write(2 + jjj_a, 14, mtx_pair_above[jjj_a, 13])
            sheet1_pair.write(2 + jjj_a, 15, mtx_pair_above[jjj_a, 14])
            sheet1_pair.write(2 + jjj_a, 16, mtx_pair_above[jjj_a, 15])
            sheet1_pair.write(2 + jjj_a, 17, mtx_pair_above[jjj_a, 16])
            sheet1_pair.write(2 + jjj_a, 18, mtx_pair_above[jjj_a, 17])

        for jjj_b in range(mtx_pair_below.shape[0]):
            sheet2_pair.write(2 + jjj_b, 0, mtx_pair_below[jjj_b, 0])
            sheet2_pair.write(2 + jjj_b, 1, mtx_pair_below[jjj_b, 1])
            sheet2_pair.write(2 + jjj_b, 2, mtx_pair_below[jjj_b, 2])
            sheet2_pair.write(2 + jjj_b, 3, mtx_pair_below[jjj_b, 3])
            sheet2_pair.write(2 + jjj_b, 4, mtx_pair_below[jjj_b, 4])
            sheet2_pair.write(2 + jjj_b, 5, mtx_pair_below[jjj_b, 5])
            sheet2_pair.write(2 + jjj_b, 6, mtx_pair_below[jjj_b, 6])
            sheet2_pair.write(2 + jjj_b, 7, mtx_pair_below[jjj_b, 7])
            sheet2_pair.write(2 + jjj_b, 8, mtx_pair_below[jjj_b, 8])
            sheet2_pair.write(2 + jjj_b, 10, mtx_pair_below[jjj_b, 9])
            sheet2_pair.write(2 + jjj_b, 11, mtx_pair_below[jjj_b, 10])
            sheet2_pair.write(2 + jjj_b, 12, mtx_pair_below[jjj_b, 11])
            sheet2_pair.write(2 + jjj_b, 13, mtx_pair_below[jjj_b, 12])
            sheet2_pair.write(2 + jjj_b, 14, mtx_pair_below[jjj_b, 13])
            sheet2_pair.write(2 + jjj_b, 15, mtx_pair_below[jjj_b, 14])
            sheet2_pair.write(2 + jjj_b, 16, mtx_pair_below[jjj_b, 15])
            sheet2_pair.write(2 + jjj_b, 17, mtx_pair_below[jjj_b, 16])
            sheet2_pair.write(2 + jjj_b, 18, mtx_pair_below[jjj_b, 17])

        sheet1_pair.write(3, 23, "Average Vols")
        sheet1_pair.write(3, 24, params_red_vols_above[0])
        sheet1_pair.write(4, 23, "Std Vols")
        sheet1_pair.write(4, 24, params_red_vols_above[1])
        sheet1_pair.write(6, 23, "Average Ints")
        sheet1_pair.write(6, 24, params_red_ints_above[0])
        sheet1_pair.write(7, 23, "Std Ints")
        sheet1_pair.write(7, 24, params_red_ints_above[1])

        sheet2_pair.write(3, 23, "Average Vols")
        sheet2_pair.write(3, 24, params_red_vols_below[0])
        sheet2_pair.write(4, 23, "Std Vols")
        sheet2_pair.write(4, 24, params_red_vols_below[1])
        sheet2_pair.write(6, 23, "Average Ints")
        sheet2_pair.write(6, 24, params_red_ints_below[0])
        sheet2_pair.write(7, 23, "Std Ints")
        sheet2_pair.write(7, 24, params_red_ints_below[1])


        # WITHOUT ABOVE-BELOW SPLITTING
        sheet1_pair.write(20, 23, "Average Vols TOT")
        sheet1_pair.write(20, 24, params_red_vols_tot[0])
        sheet1_pair.write(21, 23, "Std Vols TOT")
        sheet1_pair.write(21, 24, params_red_vols_tot[1])
        sheet1_pair.write(23, 23, "Average Ints TOT")
        sheet1_pair.write(23, 24, params_red_ints_tot[0])
        sheet1_pair.write(24, 23, "Std Ints TOT")
        sheet1_pair.write(24, 24, params_red_ints_tot[1])

        sheet2_pair.write(20, 23, "Average Vols TOT")
        sheet2_pair.write(20, 24, params_red_vols_tot[0])
        sheet2_pair.write(21, 23, "Std Vols TOT")
        sheet2_pair.write(21, 24, params_red_vols_tot[1])
        sheet2_pair.write(23, 23, "Average Ints TOT")
        sheet2_pair.write(23, 24, params_red_ints_tot[0])
        sheet2_pair.write(24, 23, "Std Ints TOT")
        sheet2_pair.write(24, 24, params_red_ints_tot[1])
        # ####################################

        sheet3_pair.write(0, 0, "Green")
        sheet3_pair.write(1, 0, "Spt_id")
        sheet3_pair.write(1, 1, "Ints")
        sheet3_pair.write(1, 2, "bkg")
        sheet3_pair.write(1, 3, "Ints - bkg")
        sheet3_pair.write(1, 4, "Ints / bkg")
        sheet3_pair.write(1, 5, "Vols")
        sheet3_pair.write(1, 6, "z-coord")
        sheet3_pair.write(1, 7, "x-coord")
        sheet3_pair.write(1, 8, "y-coord")

        sheet4_pair.write(0, 0, "Green")
        sheet4_pair.write(1, 0, "Spt_id")
        sheet4_pair.write(1, 1, "Ints")
        sheet4_pair.write(1, 2, "bkg")
        sheet4_pair.write(1, 3, "Ints - bkg")
        sheet4_pair.write(1, 4, "Ints / bkg")
        sheet4_pair.write(1, 5, "Vols")
        sheet4_pair.write(1, 6, "z-coord")
        sheet4_pair.write(1, 7, "x-coord")
        sheet4_pair.write(1, 8, "y-coord")

        idx_w3     =  0
        idx_title  =  0
        for kk_mm_a in multi_over_green_tags_above:
            line_list  =  []
            red_tag    =  []
            for jj in range(rgp_green["coords"][kk_mm_a].shape[0]):
                red_tag.append(spts_red[rgp_green["coords"][kk_mm_a][jj][0], rgp_green["coords"][kk_mm_a][jj][1], rgp_green["coords"][kk_mm_a][jj][2]])

            bkg      =  bkg_list_a[tags_list_a.index(rgp_green["label"][kk_mm_a])]
            spt_int  =  np.sum(rgp_green["intensity_image"][kk_mm_a])

            sheet3_pair.write(2 + idx_w3, 0, int(rgp_green["label"][kk_mm_a]))
            sheet3_pair.write(2 + idx_w3, 1, int(spt_int))
            sheet3_pair.write(2 + idx_w3, 2, bkg)
            sheet3_pair.write(2 + idx_w3, 3, spt_int - bkg * rgp_green["area"][kk_mm_a])
            sheet3_pair.write(2 + idx_w3, 4, spt_int / bkg)
            sheet3_pair.write(2 + idx_w3, 5, int(rgp_green["area"][kk_mm_a]))
            sheet3_pair.write(2 + idx_w3, 6, int(rgp_green["centroid-0"][kk_mm_a]))
            sheet3_pair.write(2 + idx_w3, 7, int(rgp_green["centroid-0"][kk_mm_a]))
            sheet3_pair.write(2 + idx_w3, 8, int(rgp_green["centroid-0"][kk_mm_a]))

            red_tag    =  np.trim_zeros(np.unique(red_tag))
            idx_col    =  1
            for k_m_a in range(red_tag.size):
                idx_title  =  max(idx_title, red_tag.size)
                bkg_r      =  bkg_list_b[tags_list_b.index(red_tag[k_m_a])]
                red_idx    =  np.where(rgp_red["label"] == red_tag[k_m_a])[0][0]
                spt_int_r  =  np.sum(rgp_red["intensity_image"][red_idx])

                sheet3_pair.write(2 + idx_w3, 0 + 10 * idx_col, int(red_tag[k_m_a]))
                sheet3_pair.write(2 + idx_w3, 1 + 10 * idx_col, int(spt_int_r))
                sheet3_pair.write(2 + idx_w3, 2 + 10 * idx_col, bkg_r)
                sheet3_pair.write(2 + idx_w3, 3 + 10 * idx_col, spt_int_r - bkg_r * rgp_red["area"][red_idx])
                sheet3_pair.write(2 + idx_w3, 4 + 10 * idx_col, spt_int_r / bkg_r)
                sheet3_pair.write(2 + idx_w3, 5 + 10 * idx_col, int(rgp_red["area"][red_idx]))
                sheet3_pair.write(2 + idx_w3, 6 + 10 * idx_col, int(rgp_red["centroid-0"][red_idx]))
                sheet3_pair.write(2 + idx_w3, 7 + 10 * idx_col, int(rgp_red["centroid-0"][red_idx]))
                sheet3_pair.write(2 + idx_w3, 8 + 10 * idx_col, int(rgp_red["centroid-0"][red_idx]))

                idx_col  +=  1

            idx_w3  +=  1

        for i_title in range(1, idx_title + 1):
            sheet3_pair.write(0, 10 * i_title, "RED")
            sheet3_pair.write(1, 10 * i_title, "Spt_id")
            sheet3_pair.write(1, 1 + 10 * i_title, "Ints")
            sheet3_pair.write(1, 2 + 10 * i_title, "bkg")
            sheet3_pair.write(1, 3 + 10 * i_title, "Ints - bkg")
            sheet3_pair.write(1, 4 + 10 * i_title, "Ints / bkg")
            sheet3_pair.write(1, 5 + 10 * i_title, "Vols")
            sheet3_pair.write(1, 6 + 10 * i_title, "z-coords")
            sheet3_pair.write(1, 7 + 10 * i_title, "x-coords")
            sheet3_pair.write(1, 8 + 10 * i_title, "y-coords")


        idx_title  =  0
        idx_w4     =  0
        for kk_mm_b in multi_over_green_tags_below:
            line_list  =  []
            red_tag    =  []
            for jj in range(rgp_green["coords"][kk_mm_b].shape[0]):
                red_tag.append(spts_red[rgp_green["coords"][kk_mm_b][jj][0], rgp_green["coords"][kk_mm_b][jj][1], rgp_green["coords"][kk_mm_b][jj][2]])

            bkg      =  bkg_list_a[tags_list_a.index(rgp_green["label"][kk_mm_b])]
            spt_int  =  np.sum(rgp_green["intensity_image"][kk_mm_b])

            sheet4_pair.write(2 + idx_w4, 0, int(rgp_green["label"][kk_mm_b]))
            sheet4_pair.write(2 + idx_w4, 1, int(spt_int))
            sheet4_pair.write(2 + idx_w4, 2, bkg)
            sheet4_pair.write(2 + idx_w4, 3, spt_int - bkg * rgp_green["area"][kk_mm_b])
            sheet4_pair.write(2 + idx_w4, 4, spt_int / bkg)
            sheet4_pair.write(2 + idx_w4, 5, int(rgp_green["area"][kk_mm_b]))
            sheet4_pair.write(2 + idx_w4, 6, int(rgp_green["centroid-0"][kk_mm_b]))
            sheet4_pair.write(2 + idx_w4, 7, int(rgp_green["centroid-0"][kk_mm_b]))
            sheet4_pair.write(2 + idx_w4, 8, int(rgp_green["centroid-0"][kk_mm_b]))

            red_tag    =  np.trim_zeros(np.unique(red_tag))
            idx_col    =  1
            for k_m_b in range(red_tag.size):
                idx_title  =  max(idx_title, red_tag.size)
                bkg_r      =  bkg_list_b[tags_list_b.index(red_tag[k_m_b])]
                red_idx    =  np.where(rgp_red["label"] == red_tag[k_m_b])[0][0]
                spt_int_r  =  np.sum(rgp_red["intensity_image"][red_idx])

                sheet4_pair.write(2 + idx_w4, 0 + 10 * idx_col, int(red_tag[k_m_b]))
                sheet4_pair.write(2 + idx_w4, 1 + 10 * idx_col, int(spt_int_r))
                sheet4_pair.write(2 + idx_w4, 2 + 10 * idx_col, bkg_r)
                sheet4_pair.write(2 + idx_w4, 3 + 10 * idx_col, spt_int_r - bkg_r * rgp_red["area"][red_idx])
                sheet4_pair.write(2 + idx_w4, 4 + 10 * idx_col, spt_int_r / bkg_r)
                sheet4_pair.write(2 + idx_w4, 5 + 10 * idx_col, int(rgp_red["area"][red_idx]))
                sheet4_pair.write(2 + idx_w4, 6 + 10 * idx_col, int(rgp_red["centroid-0"][red_idx]))
                sheet4_pair.write(2 + idx_w4, 7 + 10 * idx_col, int(rgp_red["centroid-0"][red_idx]))
                sheet4_pair.write(2 + idx_w4, 8 + 10 * idx_col, int(rgp_red["centroid-0"][red_idx]))

                idx_col  +=  1

            idx_w4  +=  1

        for i_title in range(1, idx_title + 1):
            sheet4_pair.write(0, 10 * i_title, "RED")
            sheet4_pair.write(1, 10 * i_title, "Spt_id")
            sheet4_pair.write(1, 1 + 10 * i_title, "Ints")
            sheet4_pair.write(1, 2 + 10 * i_title, "bkg")
            sheet4_pair.write(1, 3 + 10 * i_title, "Ints - bkg")
            sheet4_pair.write(1, 4 + 10 * i_title, "Ints / bkg")
            sheet4_pair.write(1, 5 + 10 * i_title, "Vols")
            sheet4_pair.write(1, 6 + 10 * i_title, "z-coords")
            sheet4_pair.write(1, 7 + 10 * i_title, "x-coords")
            sheet4_pair.write(1, 8 + 10 * i_title, "y-coords")


        book_overl_pair.save(analysis_folder + "/" + xls_a_name[:-22] + "_Green_Overlapping_volThr" + str(vol_thr) + "_clstrDist_" + str(clstr_dist) + ".xls")



        self.spts_green      =  spts_green
        self.spts_red        =  spts_red
        self.coloc_mtx       =  coloc_mtx
        self.mrna_in_transl  =  mrna_in_transl
        self.green_onmrna    =  green_onmrna







        #         sht1_idxs  +=  1
        #     if rgp_green["area"][kk_gr] <= vols_green_not_overl_below_mean_fttg + vols_green_not_overl_below_std_surf and int(rgp_green["centroid-0"][kk_gr]) >= first_last_frame[1]:

        #         red_tag  =  []
        #         for jj in range(rgp_green["coords"][kk_gr].shape[0]):
        #             red_tag.append(spts_red[rgp_green["coords"][kk_gr][jj][0], rgp_green["coords"][kk_gr][jj][1], rgp_green["coords"][kk_gr][jj][2]])
        #         red_tag  =  max(red_tag)

        #         bkg      =  bkg_list_a[tags_list_a.index(rgp_green["label"][kk_gr])]
        #         spt_int  =  np.sum(rgp_green["intensity_image"][kk_gr])
        #         sheet2_pair.write(sht2_idxs, 0, int(rgp_green["label"][kk_gr]))
        #         sheet2_pair.write(sht2_idxs, 1, int(spt_int))
        #         sheet2_pair.write(sht2_idxs, 2, bkg)
        #         sheet2_pair.write(sht2_idxs, 3, spt_int - bkg * rgp_green["area"][kk_gr])
        #         sheet2_pair.write(sht2_idxs, 4, spt_int / bkg)
        #         sheet2_pair.write(sht2_idxs, 5, int(rgp_green["area"][kk_gr]))
        #         sheet2_pair.write(sht2_idxs, 6, int(rgp_green["centroid-0"][kk_gr]))
        #         sheet2_pair.write(sht2_idxs, 7, int(rgp_green["centroid-1"][kk_gr]))
        #         sheet2_pair.write(sht2_idxs, 8, int(rgp_green["centroid-2"][kk_gr]))

        #         bkg      =  bkg_list_b[tags_list_b.index(red_tag)]
        #         spt_int  =  np.sum(rgp_red["intensity_image"][red_idx])
        #         sheet2_pair.write(sht2_idxs, 10, int(red_tag))
        #         sheet2_pair.write(sht2_idxs, 11, int(spt_int))
        #         sheet2_pair.write(sht2_idxs, 12, bkg)
        #         sheet2_pair.write(sht2_idxs, 13, spt_int - bkg * rgp_red["area"][red_idx])
        #         sheet2_pair.write(sht2_idxs, 14, spt_int / bkg)
        #         sheet2_pair.write(sht2_idxs, 15, int(rgp_red["area"][red_idx]))
        #         sheet2_pair.write(sht2_idxs, 16, int(rgp_red["centroid-0"][red_idx]))
        #         sheet2_pair.write(sht2_idxs, 17, int(rgp_red["centroid-1"][red_idx]))
        #         sheet2_pair.write(sht2_idxs, 18, int(rgp_red["centroid-2"][red_idx]))

        #         sht2_idxs  +=  1
 
        # book_overl_pair.save(analysis_folder + "/" + xls_a_name[:-22] + "_Green_Overlapping_volThr" + str(vol_thr) + "_clstrDist_" + str(clstr_dist) + ".xls")





        # def exp_model(x, a, b):
        #     return a * np.exp(- x / b)

        # bbss_vols                       =  np.linspace(min(min(vols_green_not_overl_above), min(vols_green_not_overl_below)), max(max(vols_green_not_overl_above), max(vols_green_not_overl_below)), 100)
        # hist_green_not_ovrl_above_vols  =  np.histogram(vols_green_not_overl_above, bins=bbss_vols)
        # hist_green_not_ovrl_below_vols  =  np.histogram(vols_green_not_overl_below, bins=bbss_vols)
        # w_vols  =  pg.plot(hist_green_not_ovrl_above_vols[1], hist_green_not_ovrl_above_vols[0], stepMode=True, fillLevel=0, brush=(0, 0, 255, 150), title="vols")
        # w_vols.plot(hist_green_not_ovrl_below_vols[1], hist_green_not_ovrl_below_vols[0], stepMode=True, pen='r')

        # params, cov                           =  curve_fit(exp_model, hist_green_not_ovrl_above_vols[1][:-1], hist_green_not_ovrl_above_vols[0], p0=[hist_green_not_ovrl_above_vols[1].max(), 1], maxfev=3000)
        # vols_green_not_overl_above_mean_fttg  =  params[1]
        # vols_green_not_overl_above_std_fttg   =  params[1]
        # vols_green_not_overl_above            =  np.sort(np.asarray(vols_green_not_overl_above))          
        # vols_green_not_overl_above_mean_surf  =  vols_green_not_overl_above[:int(vols_green_not_overl_above.size * 70 / 100)].mean()
        # vols_green_not_overl_above_std_surf   =  vols_green_not_overl_above[:int(vols_green_not_overl_above.size * 70 / 100)].std()

        # params, cov                           =  curve_fit(exp_model, hist_green_not_ovrl_below_vols[1][:-1], hist_green_not_ovrl_below_vols[0], p0=[hist_green_not_ovrl_below_vols[1].max(), 1], maxfev=3000)
        # vols_green_not_overl_below_mean_fttg  =  params[1]
        # vols_green_not_overl_below_std_fttg   =  params[1]
        # vols_green_not_overl_below            =  np.sort(np.asarray(vols_green_not_overl_below))          
        # vols_green_not_overl_below_mean_surf  =  vols_green_not_overl_above[:int(vols_green_not_overl_below.size * 70 / 100)].mean()
        # vols_green_not_overl_below_std_surf   =  vols_green_not_overl_above[:int(vols_green_not_overl_below.size * 70 / 100)].std()


        # bbss_ints                       =  np.linspace(min(min(ints_green_not_overl_above), min(ints_green_not_overl_below)), max(max(ints_green_not_overl_above), max(ints_green_not_overl_below)), 100)
        # hist_green_not_ovrl_above_ints  =  np.histogram(ints_green_not_overl_above, bins=bbss_ints)
        # hist_green_not_ovrl_below_ints  =  np.histogram(ints_green_not_overl_below, bins=bbss_ints)
        # w_ints                          =  pg.plot(hist_green_not_ovrl_above_ints[1], hist_green_not_ovrl_above_ints[0], stepMode=True, fillLevel=0, brush=(0, 0, 255, 150), title="ints")
        # w_ints.plot(hist_green_not_ovrl_below_ints[1], hist_green_not_ovrl_below_ints[0], stepMode=True, pen='r')

        # params, cov                           =  curve_fit(exp_model, hist_green_not_ovrl_above_ints[1][:-1], hist_green_not_ovrl_above_ints[0], p0=[hist_green_not_ovrl_above_ints[1].max(), 1], maxfev=3000)
        # ints_green_not_overl_above_mean_fttg  =  params[1]
        # ints_green_not_overl_above_std_fttg   =  params[1]
        # ints_green_not_overl_above            =  np.sort(np.asarray(ints_green_not_overl_above))          
        # ints_green_not_overl_above_mean_surf  =  ints_green_not_overl_above[:int(ints_green_not_overl_above.size * 70 / 100)].mean()
        # ints_green_not_overl_above_std_surf   =  ints_green_not_overl_above[:int(ints_green_not_overl_above.size * 70 / 100)].std()

        # params, cov                           =  curve_fit(exp_model, hist_green_not_ovrl_below_ints[1][:-1], hist_green_not_ovrl_below_ints[0], p0=[hist_green_not_ovrl_below_ints[1].max(), 1], maxfev=3000)
        # ints_green_not_overl_below_mean_fttg  =  params[1]
        # ints_green_not_overl_below_std_fttg   =  params[1]
        # ints_green_not_overl_below            =  np.sort(np.asarray(ints_green_not_overl_below))          
        # ints_green_not_overl_below_mean_surf  =  ints_green_not_overl_below[:int(ints_green_not_overl_below.size * 70 / 100)].mean()
        # ints_green_not_overl_below_std_surf   =  ints_green_not_overl_below[:int(ints_green_not_overl_below.size * 70 / 100)].std()


        # bbss_ints_mbkg                      =  np.linspace(min(min(ints_green_not_overl_abovembkg), min(ints_green_not_overl_belowmbkg)), max(max(ints_green_not_overl_abovembkg), max(ints_green_not_overl_belowmbkg)), 100)
        # hist_green_not_ovrl_above_intsmbkg  =  np.histogram(ints_green_not_overl_abovembkg, bins=bbss_ints_mbkg)
        # hist_green_not_ovrl_below_intsmbkg  =  np.histogram(ints_green_not_overl_belowmbkg, bins=bbss_ints_mbkg)
        # w_intsmbkg  =  pg.plot(hist_green_not_ovrl_above_intsmbkg[1], hist_green_not_ovrl_above_intsmbkg[0], stepMode=True, fillLevel=0, brush=(0, 0, 255, 150), title="mbkg")
        # w_intsmbkg.plot(hist_green_not_ovrl_below_intsmbkg[1], hist_green_not_ovrl_below_intsmbkg[0], stepMode=True, pen='r')

        # params, cov                               =  curve_fit(exp_model, hist_green_not_ovrl_above_intsmbkg[1][:-1], hist_green_not_ovrl_above_intsmbkg[0], p0=[hist_green_not_ovrl_above_intsmbkg[1].max(), 1], maxfev=3000)
        # ints_green_not_overl_abovembkg_mean_fttg  =  params[1]
        # ints_green_not_overl_abovembkg_std_fttg   =  params[1]
        # ints_green_not_overl_abovembkg            =  np.sort(np.asarray(ints_green_not_overl_abovembkg))          
        # ints_green_not_overl_abovembkg_mean_surf  =  ints_green_not_overl_abovembkg[:int(ints_green_not_overl_abovembkg.size * 70 / 100)].mean()
        # ints_green_not_overl_abovembkg_std_surf   =  ints_green_not_overl_abovembkg[:int(ints_green_not_overl_abovembkg.size * 70 / 100)].std()

        # params, cov                               =  curve_fit(exp_model, hist_green_not_ovrl_below_intsmbkg[1][:-1], hist_green_not_ovrl_below_intsmbkg[0], p0=[hist_green_not_ovrl_below_intsmbkg[1].max(), 1], maxfev=3000)
        # ints_green_not_overl_belowmbkg_mean_fttg  =  params[1]
        # ints_green_not_overl_belowmbkg_std_fttg   =  params[1]
        # ints_green_not_overl_belowmbkg            =  np.sort(np.asarray(ints_green_not_overl_belowmbkg))          
        # ints_green_not_overl_belowmbkg_mean_surf  =  ints_green_not_overl_belowmbkg[:int(ints_green_not_overl_belowmbkg.size * 70 / 100)].mean()
        # ints_green_not_overl_belowmbkg_std_surf   =  ints_green_not_overl_belowmbkg[:int(ints_green_not_overl_belowmbkg.size * 70 / 100)].std()

        # sht1_idxs  =  2
        # sht2_idxs  =  2
        # for kk_red in idxs_red_on_gr:
        #     # print(kk_red)
        #     if int(rgp_red["centroid-0"][kk_red]) <= first_last_frame[0]:

        #         line_list  =  []
        #         green_tag  =  []
        #         for jj in range(rgp_red["coords"][kk_red].shape[0]):
        #             green_tag.append(spts_green[rgp_red["coords"][kk_red][jj][0], rgp_red["coords"][kk_red][jj][1], rgp_red["coords"][kk_red][jj][2]])
        #         print(set(green_tag), "  ", kk_red)
        #         green_tag  =  max(green_tag)
        #         green_idx  =  np.where(rgp_green["label"] == green_tag)[0][0]
                
        #         bkg_r      =  bkg_list_b[tags_list_b.index(rgp_red["label"][kk_red])]
        #         spt_int_r  =  np.sum(rgp_red["intensity_image"][kk_red])
        #         bkg_g      =  bkg_list_a[tags_list_a.index(green_tag)]
        #         spt_int_g  =  np.sum(rgp_green["intensity_image"][green_idx])

        #         line_list.append([int(rgp_red["label"][kk_red]), int(spt_int_r), bkg_r, spt_int_r - bkg_r * rgp_red["area"][kk_red], spt_int_r / bkg_r, int(rgp_red["area"][kk_red]), int(rgp_red["centroid-0"][kk_red]), int(rgp_red["centroid-1"][kk_red]), int(rgp_red["centroid-2"][kk_red]), int(green_tag), int(spt_int_g), bkg_g, spt_int_g - bkg_g * rgp_green["area"][green_idx], spt_int_g / bkg_g, int(rgp_green["area"][green_idx]), int(rgp_green["centroid-0"][green_idx]), int(rgp_green["centroid-1"][green_idx]), int(rgp_green["centroid-2"][green_idx])])

        #         mtx_pair_above.append(line_list)
                
        #     if int(rgp_red["centroid-0"][kk_red]) >= first_last_frame[1]:

        #         line_list  =  []
        #         green_tag  =  []
        #         for jj in range(rgp_red["coords"][kk_red].shape[0]):
        #             green_tag.append(spts_green[rgp_red["coords"][kk_red][jj][0], rgp_red["coords"][kk_red][jj][1], rgp_red["coords"][kk_red][jj][2]])
        #         green_tag  =  max(green_tag)
        #         green_idx  =  np.where(rgp_green["label"] == green_tag)[0][0]
                
        #         bkg_r      =  bkg_list_b[tags_list_b.index(rgp_red["label"][kk_red])]
        #         spt_int_r  =  np.sum(rgp_red["intensity_image"][kk_red])
        #         bkg_g      =  bkg_list_a[tags_list_a.index(green_tag)]
        #         spt_int_g  =  np.sum(rgp_green["intensity_image"][green_idx])

        #         line_list.append([int(rgp_red["label"][kk_red]), int(spt_int_r), bkg_r, spt_int_r - bkg_r * rgp_red["area"][kk_red], spt_int_r / bkg_r, int(rgp_red["area"][kk_red]), int(rgp_red["centroid-0"][kk_red]), int(rgp_red["centroid-1"][kk_red]), int(rgp_red["centroid-2"][kk_red]), int(green_tag), int(spt_int_g), bkg_g, spt_int_g - bkg_g * rgp_green["area"][green_idx], spt_int_g / bkg_g, int(rgp_green["area"][green_idx]), int(rgp_green["centroid-0"][green_idx]), int(rgp_green["centroid-1"][green_idx]), int(rgp_green["centroid-2"][green_idx])])

        #         mtx_pair_below.append(line_list)


                        
        # hist_red_ints_above         =  np.histogram(red_ints_above, bins=100)
        # params_red_ints_above, cov  =  curve_fit(bimodal, hist_red_ints_above[1][:-1], hist_red_ints_above[0], p0=[hist_red_ints_above[1][np.argmax(hist_red_ints_above[0])], 1000, hist_red_ints_above[0].max(), 40000, 30, 100], maxfev=3000)
        # vv_red_ints_above           =  bimodal(hist_red_ints_above[1][1:], params_red_ints_above[0], params_red_ints_above[1], params_red_ints_above[2], params_red_ints_above[3], params_red_ints_above[4], params_red_ints_above[5])
        # w                           =  pg.plot(hist_red_ints_above[1], hist_red_ints_above[0], stepMode=True)
        # w.plot(hist_red_ints_above[1][1:], vv_red_ints_above, pen='r')

        # hist_red_ints_below         =  np.histogram(red_ints_below, bins=100)
        # params_red_ints_below, cov  =  curve_fit(bimodal, hist_red_ints_below[1][:-1], hist_red_ints_below[0], p0=[hist_red_ints_below[1][np.argmax(hist_red_ints_below[0])], 1000, hist_red_ints_below[0].max(), 40000, 30, 100], maxfev=3000)
        # vv_red_ints_below           =  bimodal(hist_red_ints_below[1][1:], params_red_ints_below[0], params_red_ints_below[1], params_red_ints_below[2], params_red_ints_below[3], params_red_ints_below[4], params_red_ints_below[5])
        # w                           =  pg.plot(hist_red_ints_below[1], hist_red_ints_below[0], stepMode=True)
        # w.plot(hist_red_ints_below[1][1:], vv_red_ints_below, pen='r')

        # hist_red_vols_above         =  np.histogram(red_vols_above, bins=100)
        # params_red_vols_above, cov  =  curve_fit(bimodal, hist_red_vols_above[1][:-1], hist_red_vols_above[0], p0=[hist_red_vols_above[1][np.argmax(hist_red_vols_above[0])], 1000, hist_red_vols_above[0].max(), 40, 30, 100], maxfev=3000)
        # vv_red_vols_above           =  bimodal(hist_red_vols_above[1][1:], params_red_vols_above[0], params_red_vols_above[1], params_red_vols_above[2], params_red_vols_above[3], params_red_vols_above[4], params_red_vols_above[5])
        # w                           =  pg.plot(hist_red_vols_above[1], hist_red_vols_above[0], stepMode=True)
        # w.plot(hist_red_vols_above[1][1:], vv_red_vols_above, pen='r')

        # hist_red_vols_below         =  np.histogram(red_vols_below, bins=100)
        # params_red_vols_below, cov  =  curve_fit(bimodal, hist_red_vols_below[1][:-1], hist_red_vols_below[0], p0=[hist_red_vols_below[1][np.argmax(hist_red_vols_below[0])], 1000, hist_red_vols_below[0].max(), 40, 30, 100], maxfev=3000)
        # vv_red_vols_below           =  bimodal(hist_red_vols_below[1][1:], params_red_vols_below[0], params_red_vols_below[1], params_red_vols_below[2], params_red_vols_below[3], params_red_vols_below[4], params_red_vols_below[5])
        # w                           =  pg.plot(hist_red_vols_below[1], hist_red_vols_below[0], stepMode=True)
        # w.plot(hist_red_vols_below[1][1:], vv_red_vols_below, pen='r')





    
        # hist4sm_red_vols_above      =  np.histogram(mtx_pair_above[:, 5], bins=100)
        # params_red_vols_above, cov  =  curve_fit(bimodal, hist4sm_red_vols_above[1][:-1], hist4sm_red_vols_above[0], p0=[hist4sm_red_vols_above[1][np.argmax(hist4sm_red_vols_above[0])], 10, hist4sm_red_vols_above[0].max(), 120, 30, 50], maxfev=3000)
        # vv_red_vols_above           =  bimodal(hist4sm_red_vols_above[1][1:], params_red_vols_above[0], params_red_vols_above[1], params_red_vols_above[2], params_red_vols_above[3], params_red_vols_above[4], params_red_vols_above[5])
        # w                           =  pg.plot(hist4sm_red_vols_above[1], hist4sm_red_vols_above[0], stepMode=True)
        # w.plot(hist4sm_red_vols_above[1][1:], vv_red_vols_above, pen='r')

        
        # hist4sm_red_vols_below      =  np.histogram(mtx_pair_below[:, 5], bins=100)
        # params_red_vols_below, cov  =  curve_fit(bimodal, hist4sm_red_vols_below[1][:-1], hist4sm_red_vols_below[0], p0=[hist4sm_red_vols_below[1][np.argmax(hist4sm_red_vols_below[0])], 10, hist4sm_red_vols_below[0].max(), 120, 30, 50], maxfev=3000)
        # vv_red_vols_below           =  bimodal(hist4sm_red_vols_below[1][1:], params_red_vols_below[0], params_red_vols_below[1], params_red_vols_below[2], params_red_vols_below[3], params_red_vols_below[4], params_red_vols_below[5])
        # w                           =  pg.plot(hist4sm_red_vols_below[1], hist4sm_red_vols_below[0], stepMode=True)
        # w.plot(hist4sm_red_vols_below[1][1:], vv_red_vols_below, pen='r')

        
        # hist4sm_red_ints_above      =  np.histogram(mtx_pair_above[:, 1], bins=100)
        # params_red_ints_above, cov  =  curve_fit(bimodal, hist4sm_red_ints_above[1][:-1], hist4sm_red_ints_above[0], p0=[hist4sm_red_ints_above[1][np.argmax(hist4sm_red_ints_above[0])], hist4sm_red_ints_above[1][np.argmax(hist4sm_red_ints_above[0])], hist4sm_red_ints_above[0].max(), 120, 30, 50], maxfev=3000)
        # vv_red_ints_above           =  bimodal(hist4sm_red_ints_above[1][1:], params_red_ints_above[0], params_red_ints_above[1], params_red_ints_above[2], params_red_ints_above[3], params_red_ints_above[4], params_red_ints_above[5])
        # w                           =  pg.plot(hist4sm_red_ints_above[1], hist4sm_red_ints_above[0], stepMode=True)
        # w.plot(hist4sm_red_ints_above[1][1:], vv_red_ints_above, pen='r')

        
        # hist4sm_red_ints_below      =  np.histogram(mtx_pair_below[:, 1], bins=100)
        # params_red_ints_below, cov  =  curve_fit(bimodal, hist4sm_red_ints_below[1][:-1], hist4sm_red_ints_below[0], p0=[hist4sm_red_ints_below[1][np.argmax(hist4sm_red_ints_below[0])], 10, hist4sm_red_ints_below[0].max(), 120, 30, 50], maxfev=3000)
        # vv_red_ints_below           =  bimodal(hist4sm_red_ints_below[1][1:], params_red_ints_below[0], params_red_ints_below[1], params_red_ints_below[2], params_red_ints_below[3], params_red_ints_below[4], params_red_ints_below[5])
        # w                           =  pg.plot(hist4sm_red_ints_below[1], hist4sm_red_ints_below[0], stepMode=True)
        # w.plot(hist4sm_red_ints_below[1][1:], vv_red_ints_below, pen='r')
