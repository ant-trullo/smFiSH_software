import numpy as np
from skimage.measure import regionprops_table, regionprops
from scipy.optimize import curve_fit
import xlwt

import RawDataLoader


def gauss(x, mu, sigma, A):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


# analysis_folder  =  '/home/atrullo/Dropbox/Suntag paper/Figure 3/JD134 CJD_cc14Zoom4_Airyscan Processing1/Cluster0_vol32/'
# raw_data_fname   =  '/home/atrullo/Dropbox/Suntag paper/Figure 3/JD134 CJD_cc14Zoom4_Airyscan Processing1/JD134 CJD_cc14Zoom4_Airyscan Processing1.czi'

raw_data  =  RawDataLoader.RawDataLoader(raw_data_fname, np.fromfile(analysis_folder + 'chs_spts_nucs.bin', 'uint16'))
spts_red  =  np.load(analysis_folder + 'spts_lbls_b.npy')


rgp_red  =  regionprops_table(spts_red, raw_data.spts_b, properties=["label", "coords", "intensity_image", "area", "centroid"])  # labels and coordinates of all the transcription spots

rgp_cc   =  regionprops(spts_red)  # regionprops to measure the centroids coordinates (the prvious gives approximation of it)
rgp_arr  =  []
for k in range(len(rgp_red["label"])):
    rgp_arr.append(
        [rgp_cc[k]["label"], int(np.round(rgp_cc[k]["centroid"][0])), int(np.round(rgp_cc[k]["centroid"][1])),
         int(np.round(rgp_cc[k]["centroid"][2]))])  # array with not approximate centroids-coordinates and label

rgp_arr  =  np.asarray(rgp_arr)  # convert to array

for k in range(rgp_red["label"].size):  # substitute the approximated centroids-coordinates with the non approximated ones
    idx2subst                                                                     =  np.where(rgp_arr[:, 0] == rgp_red["label"][k])[0][0]  # check the indexes of the label
    rgp_red["centroid-0"][k], rgp_red["centroid-1"][k], rgp_red["centroid-2"][k]  =  rgp_arr[idx2subst][1], rgp_arr[idx2subst][2], rgp_arr[idx2subst][3]  # substitute centroids-coordinates with the proper not approximated ones

first_last  =  np.load(analysis_folder + 'first_last_frame.npy')

red_above_1  =  []
red_above_2  =  []
red_above_3  =  []

for kk in range(rgp_red["label"].size):
    if rgp_red["centroid-0"][kk] > first_last[1] and rgp_red["centroid-0"][kk] <= first_last[1] + 10:
        red_above_1.append([rgp_red["label"][kk], rgp_red["area"][kk], np.sum(rgp_red["intensity_image"][kk]), rgp_red["centroid-0"][kk]])
    if rgp_red["centroid-0"][kk] > first_last[1] + 10 and rgp_red["centroid-0"][kk] <= first_last[1] + 20:
        red_above_2.append([rgp_red["label"][kk], rgp_red["area"][kk], np.sum(rgp_red["intensity_image"][kk]), rgp_red["centroid-0"][kk]])
    if rgp_red["centroid-0"][kk] > first_last[1] + 20:
        red_above_3.append([rgp_red["label"][kk], rgp_red["area"][kk], np.sum(rgp_red["intensity_image"][kk]), rgp_red["centroid-0"][kk]])

red_above_1  =   np.asarray(red_above_1)
red_above_2  =   np.asarray(red_above_2)
red_above_3  =   np.asarray(red_above_3)

hist_red_ints_1  =  np.histogram(red_above_1[:, 2], bins=100)
params_1, cov    =  curve_fit(gauss, hist_red_ints_1[1][:-1], hist_red_ints_1[0], p0=[hist_red_ints_1[1][np.argmax(hist_red_ints_1[0])], hist_red_ints_1[1][np.argmax(hist_red_ints_1[0])] / 10, hist_red_ints_1[0].max()], maxfev=3000)

hist_red_ints_2  =  np.histogram(red_above_2[:, 2], bins=100)
params_2, cov    =  curve_fit(gauss, hist_red_ints_2[1][:-1], hist_red_ints_2[0], p0=[hist_red_ints_2[1][np.argmax(hist_red_ints_2[0])], hist_red_ints_2[1][np.argmax(hist_red_ints_2[0])] / 10, hist_red_ints_2[0].max()], maxfev=3000)

hist_red_ints_3  =  np.histogram(red_above_3[:, 2], bins=100)
params_3, cov    =  curve_fit(gauss, hist_red_ints_3[1][:-1], hist_red_ints_3[0], p0=[hist_red_ints_3[1][np.argmax(hist_red_ints_3[0])], hist_red_ints_3[1][np.argmax(hist_red_ints_3[0])] / 10, hist_red_ints_3[0].max()], maxfev=3000)

book    =  xlwt.Workbook(encoding='utf-8')
sheet1  =  book.add_sheet("Sheet1")

sheet1.write(0, 0, raw_data_fname)

sheet1.write(1, 0, "z min")
sheet1.write(1, 1, "z max")
sheet1.write(1, 2, "numb os spts")
sheet1.write(1, 3, "mean")
sheet1.write(1, 4, "std")

sheet1.write(2, 0, int(first_last[1]))
sheet1.write(2, 1, int(first_last[1]) + 10)
sheet1.write(2, 2, red_above_1.size)
sheet1.write(2, 3, params_1[0])
sheet1.write(2, 4, params_1[1])

sheet1.write(3, 0, int(first_last[1]) + 10)
sheet1.write(3, 1, int(first_last[1]) +  20)
sheet1.write(3, 2, red_above_2.size)
sheet1.write(3, 3, params_2[0])
sheet1.write(3, 4, params_2[1])

sheet1.write(4, 0, int(first_last[1]) + 20)
sheet1.write(4, 1, spts_red.shape[0])
sheet1.write(4, 2, red_above_3.size)
sheet1.write(4, 3, params_3[0])
sheet1.write(4, 4, params_3[1])

xls2add  =  raw_data_fname[:-4]
book.save(analysis_folder + xls2add[len(xls2add) - xls2add[::-1].find('/'):] + '_SM_int_study_overZ.xls')