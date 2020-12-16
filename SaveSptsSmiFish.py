"""This function write data of the spots detection."""


import datetime
import numpy as np
import xlwt
from skimage.measure import regionprops_table



class SaveSptsSmiFish:
    """Main class, does all the job"""
    def __init__(self, folder2write, soft_version, raw_data_fname, raw_spts, spts_thr_study, spts_lbls, bkg_cages, info_mature, info_nascent, thr_vals, int_thr_value, sld_thr_current, vol_thr_value, min_numb_z_value, sph_thr, segm_clstr_flag, a_b_flag):


        book   =  xlwt.Workbook(encoding='utf-8')
        sheet  =  book.add_sheet("Sheet 1")

        sheet.write(0, 24, "Filename")
        sheet.write(0, 25, raw_data_fname)
        sheet.write(1, 24, "Date")
        sheet.write(1, 25, datetime.datetime.now().strftime("%d-%b-%Y"))
        sheet.write(2, 24, "Software version")
        sheet.write(2, 25, soft_version)

        sheet.write(0, 0, "SingleMol")
        sheet.write(2, 0, "Spot_ID")
        sheet.write(2, 1, "Spot Vol")
        sheet.write(2, 2, "Spot_Tot_Int")
        sheet.write(2, 3, "Spot_Max_Int")
        sheet.write(2, 4, "Spot_Background")
        sheet.write(2, 5, "Spot_Int/Bkg")
        sheet.write(2, 6, "Spot_IntMax/Bkg")
        sheet.write(2, 7, "Spot_Int-Bkg")
        sheet.write(2, 8, "Spot_Z_Pos")
        sheet.write(2, 9, "Spot_X_Pos")
        sheet.write(2, 10, "Spot_Y_Pos")

        sheet.write(0, 12, "TS")
        sheet.write(2, 12, "TS_ID")
        sheet.write(2, 13, "TS_Vol")
        sheet.write(2, 14, "TS_Tot_Int")
        sheet.write(2, 15, "TS_Max_Int")
        sheet.write(2, 16, "TS_Background")
        sheet.write(2, 17, "TS_Int/Bkg")
        sheet.write(2, 18, "TS_IntMax/Bkg")
        sheet.write(2, 19, "TS_Int-Bkg")
        sheet.write(2, 20, "TS_Z_Pos")
        sheet.write(2, 21, "TS_X_Pos")
        sheet.write(2, 22, "TS_Y_Pos")

        rgp_bkg  =  regionprops_table(bkg_cages, raw_spts, properties=["label", "area", "intensity_image"])
        print(np.unique(rgp_bkg)[1:].size)

        for k in range(info_mature.shape[1]):
            print(k)
            sheet.write(3 + k, 0, "Spt_" + str(info_mature[0, k]))
            sheet.write(3 + k, 1, int(info_mature[1, k]))
            sheet.write(3 + k, 2, int(info_mature[2, k]))
            sheet.write(3 + k, 3, int(info_mature[3, k]))

            bkg_cage_idx  =  np.where(rgp_bkg["label"] == info_mature[0, k])[0][0]
            bkg_cage      =  rgp_bkg["intensity_image"][bkg_cage_idx].sum() / rgp_bkg["area"][bkg_cage_idx]
            sheet.write(3 + k, 4, bkg_cage)
            sheet.write(3 + k, 5, info_mature[2, k] / bkg_cage)

            sheet.write(3 + k, 6, info_mature[3, k] / bkg_cage)
            sheet.write(3 + k, 7, info_mature[2, k] - info_mature[1, k] * bkg_cage)
            sheet.write(3 + k, 8, int(info_mature[4, k]))
            sheet.write(3 + k, 9, int(info_mature[5, k]))
            sheet.write(3 + k, 10, int(info_mature[6, k]))

        for l in range(info_nascent.shape[1]):
            sheet.write(3 + l, 12, "Spt_" + str(info_nascent[0, l]))
            sheet.write(3 + l, 13, int(info_nascent[1, l]))
            sheet.write(3 + l, 14, int(info_nascent[2, l]))
            sheet.write(3 + l, 15, int(info_nascent[3, l]))

            bkg_cage_idx  =  np.where(rgp_bkg["label"] == info_nascent[0, l])[0][0]
            bkg_cage      =  rgp_bkg["intensity_image"][bkg_cage_idx].sum() / rgp_bkg["area"][bkg_cage_idx]
            sheet.write(3 + l, 4, bkg_cage)
            sheet.write(3 + l, 5, info_nascent[2, l] / bkg_cage)

            sheet.write(3 + l, 17, info_nascent[2, l] / bkg_cage)
            sheet.write(3 + l, 18, info_nascent[3, l] / bkg_cage)
            sheet.write(3 + l, 19, info_nascent[2, l] - info_nascent[1, l] * bkg_cage)
            sheet.write(3 + l, 20, int(info_nascent[4, l]))
            sheet.write(3 + l, 21, int(info_nascent[5, l]))
            sheet.write(3 + l, 22, int(info_nascent[6, l]))


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




