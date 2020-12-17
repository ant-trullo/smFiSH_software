"""This function write data of the spots detection."""


import datetime
import numpy as np
import xlwt
from skimage.measure import regionprops_table
from scipy.optimize import curve_fit
import pyqtgraph as pg
import pyqtgraph.exporters


class SaveSptsSmiFish:
    """Main class, does all the job"""
    def __init__(self, folder2write, soft_version, raw_data_fname, raw_spts, spts_thr_study, spts_lbls, bkg_cages, info_mature, info_nascent, thr_vals, int_thr_value, sld_thr_current, vol_thr_value, min_numb_z_value, sph_thr, segm_clstr_flag, a_b_flag):

        book    =  xlwt.Workbook(encoding='utf-8')
        sheet1  =  book.add_sheet("SM")
        sheet2  =  book.add_sheet("Clusters")
        sheet3  =  book.add_sheet("TS")
        sheet4  =  book.add_sheet("Total")

        sheet1.write(0, 13, "Filename")
        sheet1.write(0, 14, raw_data_fname)
        sheet1.write(1, 13, "Date")
        sheet1.write(1, 14, datetime.datetime.now().strftime("%d-%b-%Y"))
        sheet1.write(2, 13, "Software version")
        sheet1.write(2, 14, soft_version)

        sheet1.write(0, 0, "SingleMol")
        sheet1.write(2, 0, "Spot_ID")
        sheet1.write(2, 1, "Spot Vol")
        sheet1.write(2, 2, "Spot_Tot_Int")
        sheet1.write(2, 3, "Spot_Max_Int")
        sheet1.write(2, 4, "Spot_Background")
        sheet1.write(2, 5, "Spot_Int/Bkg")
        sheet1.write(2, 6, "Spot_IntMax/Bkg")
        sheet1.write(2, 7, "Spot_Int-Bkg")
        sheet1.write(2, 8, "Spot_Z_Pos")
        sheet1.write(2, 9, "Spot_X_Pos")
        sheet1.write(2, 10, "Spot_Y_Pos")

        sheet2.write(0, 0, "Cl")
        sheet2.write(2, 0, "Cl_ID")
        sheet2.write(2, 1, "Cl_Vol")
        sheet2.write(2, 2, "Cl_Tot_Int")
        sheet2.write(2, 3, "Cl_Max_Int")
        sheet2.write(2, 4, "Cl_Background")
        sheet2.write(2, 5, "Cl_Int/Bkg")
        sheet2.write(2, 6, "Cl_IntMax/Bkg")
        sheet2.write(2, 7, "Cl_Int-Bkg")
        sheet2.write(2, 8, "Cl_Z_Pos")
        sheet2.write(2, 9, "Cl_X_Pos")
        sheet2.write(2, 10, "Cl_Y_Pos")
        sheet2.write(2, 11, "Equiv SM numb")

        sheet3.write(0, 0, "TS")
        sheet3.write(2, 0, "TS_ID")
        sheet3.write(2, 1, "TS_Vol")
        sheet3.write(2, 2, "TS_Tot_Int")
        sheet3.write(2, 3, "TS_Max_Int")
        sheet3.write(2, 4, "TS_Background")
        sheet3.write(2, 5, "TS_Int/Bkg")
        sheet3.write(2, 6, "TS_IntMax/Bkg")
        sheet3.write(2, 7, "TS_Int-Bkg")
        sheet3.write(2, 8, "TS_Z_Pos")
        sheet3.write(2, 9, "TS_X_Pos")
        sheet3.write(2, 10, "TS_Y_Pos")

        sheet4.write(0, 0, "Total")
        sheet4.write(2, 0, "ID")
        sheet4.write(2, 1, "Vol")
        sheet4.write(2, 2, "Tot_Int")
        sheet4.write(2, 3, "Max_Int")
        sheet4.write(2, 4, "Background")
        sheet4.write(2, 5, "Int/Bkg")
        sheet4.write(2, 6, "IntMax/Bkg")
        sheet4.write(2, 7, "Int-Bkg")
        sheet4.write(2, 8, "Z_Pos")
        sheet4.write(2, 9, "X_Pos")
        sheet4.write(2, 10, "Y_Pos")


        def gauss(x, mu, sigma, A):
            """Gaussian function for fitting"""
            return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


#         def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
#             return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

        rgp_bkg   =  regionprops_table(bkg_cages, raw_spts, properties=["label", "area", "intensity_image"])
        rgp_spts  =  regionprops_table(spts_lbls, raw_spts, properties=["label", "area", "intensity_image", "centroid"])

        rgp_spts["centroid-0"]  =  np.round(rgp_spts["centroid-0"]).astype(int)
        rgp_spts["centroid-1"]  =  np.round(rgp_spts["centroid-1"]).astype(int)
        rgp_spts["centroid-2"]  =  np.round(rgp_spts["centroid-2"]).astype(int)

        for mm in range(len(rgp_spts["label"])):
            sheet4.write(3 + mm, 0, "Spt_" + str(rgp_spts["label"][mm]))
            sheet4.write(3 + mm, 1, int(rgp_spts["area"][mm]))
            sheet4.write(3 + mm, 2, int(np.sum(rgp_spts["intensity_image"][mm])))
            sheet4.write(3 + mm, 3, int(np.max(rgp_spts["intensity_image"][mm])))
            try:
                bkg_cage_idx  =  np.where(rgp_bkg["label"] == rgp_spts["label"][mm])[0][0]
                bkg_cage      =  rgp_bkg["intensity_image"][bkg_cage_idx].sum() / rgp_bkg["area"][bkg_cage_idx]
                sheet4.write(3 + mm, 4, bkg_cage)
                sheet4.write(3 + mm, 5, np.sum(rgp_spts["intensity_image"][mm] / bkg_cage))
                sheet4.write(3 + mm, 6, np.max(rgp_spts["intensity_image"][mm] / bkg_cage))
                sheet4.write(3 + mm, 7, np.sum(rgp_spts["intensity_image"][mm]) - bkg_cage * rgp_spts["area"][mm])
            except IndexError:
                pass

            sheet4.write(3 + mm, 8, int(rgp_spts["centroid-0"][mm]))
            sheet4.write(3 + mm, 9, int(rgp_spts["centroid-1"][mm]))
            sheet4.write(3 + mm, 10, int(rgp_spts["centroid-2"][mm]))

        spts_sm  =  np.sort(info_mature[2, :])
        hist_sm  =  np.histogram(spts_sm, bins=500)

        params, cov  =  curve_fit(gauss, hist_sm[1][:-1], hist_sm[0], p0=[hist_sm[1][np.argmax(hist_sm[0])], hist_sm[1][np.argmax(hist_sm[0])] / 2, hist_sm[0].max()], bounds=((0, 0, 0), (+np.inf, +np.inf, +np.inf)), maxfev=8000)
        # params, cov  =  curve_fit(bimodal, hist_sm[1][:-1], hist_sm[0], p0=[hist_sm[1][np.argmax(hist_sm[0])], 0, hist_sm[0].max(), 20000, 500, 200], bounds=((0, 0, 0, 0, 0, 0), (+np.inf, +np.inf, +np.inf, +np.inf, +np.inf, +np.inf)), maxfev=8000)
        sm_upbound   =  params[0] + 3 * params[1]

        sheet1.write(5, 13, "Mean")
        sheet1.write(5, 14, params[0])
        sheet1.write(6, 13, "STD")
        sheet1.write(6, 14, params[1])
        sheet1.write(7, 13, "Up bound")
        sheet1.write(7, 14, sm_upbound)
        sheet1.write(9, 13, "Numb SM")
        sheet2.write(9, 13, "Numb Cl")
        sheet2.write(10, 13, "Tot equiv SM numb")
        sheet3.write(9, 13, "Numb TS")
        sheet3.write(9, 14, info_nascent.shape[1])

        yy  =  gauss(hist_sm[1][:-1], params[0], params[1], params[2])     # , params[3], params[4], params[5])
        w   =  pg.plot(hist_sm[1], hist_sm[0], stepMode=True)
        w.plot(hist_sm[1][:-1], yy, pen='r')
        exporter  =  pg.exporters.ImageExporter(w.plotItem)
        exporter.export(folder2write + '/Hist_spots_fitting' + a_b_flag + '.png')

        book2       =  xlwt.Workbook(encoding='utf-8')
        sheet1_bis  =  book2.add_sheet("SM")

        sheet1_bis.write(0, 0, "Interval")
        sheet1_bis.write(0, 1, "Occurrence")
        sheet1_bis.write(0, 2, "Fit")

        for jj in range(len(hist_sm[0])):
            sheet1_bis.write(jj + 1, 0, int(hist_sm[1][jj]))
            sheet1_bis.write(jj + 1, 1, int(hist_sm[0][jj]))
            sheet1_bis.write(jj + 1, 2, float(yy[jj]))

        book2.save(folder2write + "/" + "DataForHistogramm" + a_b_flag + ".xls")

        idx_sm      =  0
        idx_cl      =  0
        eq_SM_numb  =  0
        for k in range(info_mature.shape[1]):
            if info_mature[2, k] <= sm_upbound:
                sheet1.write(3 + idx_sm, 0, "Spt_" + str(info_mature[0, k]))
                sheet1.write(3 + idx_sm, 1, int(info_mature[1, k]))
                sheet1.write(3 + idx_sm, 2, int(info_mature[2, k]))
                sheet1.write(3 + idx_sm, 3, int(info_mature[3, k]))

                try:
                    bkg_cage_idx  =  np.where(rgp_bkg["label"] == info_mature[0, k])[0][0]
                    bkg_cage      =  rgp_bkg["intensity_image"][bkg_cage_idx].sum() / rgp_bkg["area"][bkg_cage_idx]
                    sheet1.write(3 + idx_sm, 4, bkg_cage)
                    sheet1.write(3 + idx_sm, 5, info_mature[2, k] / bkg_cage)
                except IndexError:
                    pass

                sheet1.write(3 + idx_sm, 6, info_mature[3, k] / bkg_cage)
                sheet1.write(3 + idx_sm, 7, info_mature[2, k] - info_mature[1, k] * bkg_cage)
                sheet1.write(3 + idx_sm, 8, int(info_mature[4, k]))
                sheet1.write(3 + idx_sm, 9, int(info_mature[5, k]))
                sheet1.write(3 + idx_sm, 10, int(info_mature[6, k]))
                idx_sm  +=  1
            else:
                sheet2.write(3 + idx_cl, 0, "Spt_" + str(info_mature[0, k]))
                sheet2.write(3 + idx_cl, 1, int(info_mature[1, k]))
                sheet2.write(3 + idx_cl, 2, int(info_mature[2, k]))
                sheet2.write(3 + idx_cl, 3, int(info_mature[3, k]))

                try:
                    bkg_cage_idx  =  np.where(rgp_bkg["label"] == info_mature[0, k])[0][0]
                    bkg_cage      =  rgp_bkg["intensity_image"][bkg_cage_idx].sum() / rgp_bkg["area"][bkg_cage_idx]
                    sheet2.write(3 + idx_cl, 4, bkg_cage)
                    sheet2.write(3 + idx_cl, 5, info_mature[2, k] / bkg_cage)
                except IndexError:
                    pass

                sheet2.write(3 + idx_cl, 6, info_mature[3, k] / bkg_cage)
                sheet2.write(3 + idx_cl, 7, info_mature[2, k] - info_mature[1, k] * bkg_cage)
                sheet2.write(3 + idx_cl, 8, int(info_mature[4, k]))
                sheet2.write(3 + idx_cl, 9, int(info_mature[5, k]))
                sheet2.write(3 + idx_cl, 10, int(info_mature[6, k]))
                sheet2.write(3 + idx_cl, 11, info_mature[2, k] / params[0])
                eq_SM_numb  +=  int(np.round(info_mature[2, k] / params[0]))
                idx_cl      +=  1

        sheet1.write(9, 14, idx_sm)
        sheet2.write(9, 14, idx_cl)
        sheet2.write(10, 14, eq_SM_numb)

        for l in range(info_nascent.shape[1]):
            sheet3.write(3 + l, 0, "Spt_" + str(info_nascent[0, l]))
            sheet3.write(3 + l, 1, int(info_nascent[1, l]))
            sheet3.write(3 + l, 2, int(info_nascent[2, l]))
            sheet3.write(3 + l, 3, int(info_nascent[3, l]))

            try:
                bkg_cage_idx  =  np.where(rgp_bkg["label"] == info_nascent[0, l])[0][0]
                bkg_cage      =  rgp_bkg["intensity_image"][bkg_cage_idx].sum() / rgp_bkg["area"][bkg_cage_idx]
                sheet3.write(3 + l, 4, bkg_cage)
                sheet3.write(3 + l, 5, info_nascent[2, l] / bkg_cage)
            except IndexError:
                pass

            sheet3.write(3 + l, 6, info_nascent[3, l] / bkg_cage)
            sheet3.write(3 + l, 7, info_nascent[2, l] - info_nascent[1, l] * bkg_cage)
            sheet3.write(3 + l, 8, int(info_nascent[4, l]))
            sheet3.write(3 + l, 9, int(info_nascent[5, l]))
            sheet3.write(3 + l, 10, int(info_nascent[6, l]))

        book.save(folder2write + "/" + raw_data_fname[len(raw_data_fname) - raw_data_fname[::-1].find('/'): -4] + '_Spots_intensity' + a_b_flag + '.xls')

        np.save(folder2write + '/spts_lbls' + a_b_flag + '.npy', spts_lbls)
        np.save(folder2write + '/spts_g' + a_b_flag + '.npy', spts_thr_study.spts_g)
        np.save(folder2write + '/spts_other_vals' + a_b_flag + '.npy', np.append(np.array([spts_thr_study.mu, spts_thr_study.sigma]), spts_thr_study.spts_num))
        np.save(folder2write + '/bkg_cages' + a_b_flag + '.npy', bkg_cages)
        np.save(folder2write + '/info_mature' + a_b_flag + '.npy', info_mature)
        np.save(folder2write + '/info_nascent' + a_b_flag + '.npy', info_nascent)

        file  =  open(folder2write + '/params_text' + a_b_flag + '.txt', "w")
        file.write("Minimum Threshold = " + str(int(thr_vals[0])))
        file.write('\n' + "Maximum Threshold = " + str(int(thr_vals[-1])))
        file.write('\n' + "Number of Steps = " + str(int(thr_vals.size)))
        file.write('\n' + "Current Threshold Value = " + str(int(sld_thr_current)))
        file.write('\n' + "Average Intensity Threshold = " + str(int_thr_value))
        file.write('\n' + "Volume Threshold = " + str(int(vol_thr_value)))
        file.write('\n' + "Minimum number of z plains = " + str(int(min_numb_z_value)))
        file.write('\n' + "Sphericity Threshold  = " + str(int(sph_thr)))
        file.write('\n' + "Segment Cluster Flag  = " + segm_clstr_flag)
        file.close()




