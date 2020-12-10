"""This function saves all the results of the smiFish software.

It pops-up the image of the dilated nuclei whit label numbers: you
must save it manually. It saves all the .bin files of the used images
and writes excel files with analysis info.
"""

import datetime
import numpy as np
import pyqtgraph as pg
import xlwt
from skimage.measure import regionprops
from skimage.color import label2rgb
from PyQt5 import QtGui, QtWidgets


class AnalysisSaver2:
    """Main function, does all the job"""
    def __init__(self, fwritename, soft_version, raw_spts, nucs_dil, nucs_3d_det, spts_segm, spts_mature, bkg_cages, raw_data_fname, int_thr_value):

        mycmap       =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3))
        nucs_dil_3c  =  label2rgb(nucs_dil, bg_label=0, bg_color=[0, 0, 0], colors=mycmap)

        w        =  pg.image(nucs_dil_3c)
        txt_pos  =  regionprops(nucs_dil)                                                       # regionprops of the nuclei: we use the centroids coordinates to print text on image and in the excel
        for t in range(len(txt_pos)):                                                           # map of dilated nuclei with tag-numbers on the top
            a  =  pg.TextItem(str(txt_pos[t]['label']), 'k')
            w.addItem(a)
            a.setPos(txt_pos[t]['centroid'][0] - 30, txt_pos[t]['centroid'][1] - 30)

        mtx_mature  =  np.zeros(nucs_3d_det.nucs_lbls.shape)                                                        # 3D matrix of all the spots center of mass
        steps       =  nucs_3d_det.nucs_lbls.shape[0]

        for k in range(spts_mature.shape[1]):                                                                       # populating the matrix
            mtx_mature[int(spts_mature[3, k]), int(spts_mature[4, k]), int(spts_mature[5, k])]  =  1

        mtx_mature_int  =  mtx_mature * np.sign(nucs_3d_det.nucs_lbls)
        mtx_mature_ext  =  mtx_mature * (1 - np.sign(nucs_3d_det.nucs_lbls))

        idxs  =  np.unique(nucs_dil)[1:]  # list of the nuclei tags

        book    =  xlwt.Workbook(encoding='utf-8')
        sheet1  =  book.add_sheet("Internal")
        sheet2  =  book.add_sheet("External")
        sheet3  =  book.add_sheet("Tot")
        print(a)

        sheet1.write(0, 0, "Nuc_id")
        sheet1.write(1, 0, "X coord")
        sheet1.write(2, 0, "Y coord")
        sheet1.write(3, 0, "Numb of Spts")
        sheet1.write(4, 0, "Region Volume")
        sheet1.write(5, 0, "Nucleus Volume")
        sheet1.write(6, 0, "Spots Intensity")

        sheet2.write(0, 0, "Nuc_id")
        sheet2.write(1, 0, "X coord")
        sheet2.write(2, 0, "Y coord")
        sheet2.write(3, 0, "Numb of Spts")
        sheet2.write(4, 0, "Region Volume")
        sheet2.write(5, 0, "Nucleus Volume")
        sheet2.write(6, 0, "Spots Intensity")

        sheet3.write(0, 0, "Nuc_id")
        sheet3.write(1, 0, "X coord")
        sheet3.write(2, 0, "Y coord")
        sheet3.write(3, 0, "Numb of Spts")
        sheet3.write(4, 0, "Region Volume")
        sheet3.write(5, 0, "Nucleus Volume")
        sheet3.write(6, 0, "Spots Intensity")

        book_b    =  xlwt.Workbook(encoding='utf-8')
        sheet1_b  =  book_b.add_sheet("Internal")
        sheet2_b  =  book_b.add_sheet("External")
        sheet3_b  =  book_b.add_sheet("Tot")

        sheet1_b.write(0, 0, "Nuc_id")
        sheet1_b.write(1, 0, "X coord")
        sheet1_b.write(2, 0, "Y coord")
        sheet1_b.write(3, 0, "Numb of Spts")
        sheet1_b.write(4, 0, "Region Volume")
        sheet1_b.write(5, 0, "Nucleus Volume")
        sheet1_b.write(6, 0, "Spots Intensity/Bkg")

        sheet2_b.write(0, 0, "Nuc_id")
        sheet2_b.write(1, 0, "X coord")
        sheet2_b.write(2, 0, "Y coord")
        sheet2_b.write(3, 0, "Numb of Spts")
        sheet2_b.write(4, 0, "Region Volume")
        sheet2_b.write(5, 0, "Nucleus Volume")
        sheet2_b.write(6, 0, "Spots Intensity/Bkg")

        sheet3_b.write(0, 0, "Nuc_id")
        sheet3_b.write(1, 0, "X coord")
        sheet3_b.write(2, 0, "Y coord")
        sheet3_b.write(3, 0, "Numb of Spts")
        sheet3_b.write(4, 0, "Region Volume")
        sheet3_b.write(5, 0, "Nucleus Volume")
        sheet3_b.write(6, 0, "Spots Intensity/Bkg")

        pbar1  =  ProgressBar(total1=idxs.size)
        pbar1.show()

        for q in range(idxs.size):
            # print(q)
            pbar1.update_progressbar(q)
            sheet1.write(0, q + 1, float(idxs[q]))                                                                      # writing the id
            sheet1.write(1, q + 1, txt_pos[q]['centroid'][0])                                                           # x coordinate of nuclei centroid
            sheet1.write(2, q + 1, txt_pos[q]['centroid'][1])                                                           # y coordinate of nuclei centroid
            sheet1.write(3, q + 1, np.sum((nucs_dil == idxs[q]) * mtx_mature))                                          # number of spots
            sheet1.write(4, q + 1, float((nucs_dil == idxs[q]).sum() * steps))                                          # region nuclei volume
            sheet1.write(5, q + 1, float(((nucs_dil == idxs[q]) * np.sign(nucs_3d_det.nucs_lbls)).sum()))               # nucleus volume

            sheet2.write(0, q + 1, float(idxs[q]))                                                                      # writing the id
            sheet2.write(1, q + 1, txt_pos[q]['centroid'][0])                                                           # x coordinate of nuclei centroid
            sheet2.write(2, q + 1, txt_pos[q]['centroid'][1])                                                           # y coordinate of nuclei centroid
            sheet2.write(3, q + 1, np.sum((nucs_dil == idxs[q]) * mtx_mature))                                          # number of spots
            sheet2.write(4, q + 1, float((nucs_dil == idxs[q]).sum() * steps))                                          # region nuclei volume
            sheet2.write(5, q + 1, float(((nucs_dil == idxs[q]) * np.sign(nucs_3d_det.nucs_lbls)).sum()))               # nucleus volume

            sheet3.write(0, q + 1, float(idxs[q]))                                                                      # writing the id
            sheet3.write(1, q + 1, txt_pos[q]['centroid'][0])                                                           # x coordinate of nuclei centroid
            sheet3.write(2, q + 1, txt_pos[q]['centroid'][1])                                                           # y coordinate of nuclei centroid
            sheet3.write(3, q + 1, np.sum((nucs_dil == idxs[q]) * mtx_mature))                                          # number of spots
            sheet3.write(4, q + 1, float((nucs_dil == idxs[q]).sum() * steps))                                          # region nuclei volume
            sheet3.write(5, q + 1, float(((nucs_dil == idxs[q]) * np.sign(nucs_3d_det.nucs_lbls)).sum()))               # nucleus volume

            sheet1_b.write(0, q + 1, float(idxs[q]))                                                                      # writing the id
            sheet1_b.write(1, q + 1, txt_pos[q]['centroid'][0])                                                           # x coordinate of nuclei centroid
            sheet1_b.write(2, q + 1, txt_pos[q]['centroid'][1])                                                           # y coordinate of nuclei centroid
            sheet1_b.write(3, q + 1, np.sum((nucs_dil == idxs[q]) * mtx_mature))                                          # number of spots
            sheet1_b.write(4, q + 1, float((nucs_dil == idxs[q]).sum() * steps))                                          # region nuclei volume
            sheet1_b.write(5, q + 1, float(((nucs_dil == idxs[q]) * np.sign(nucs_3d_det.nucs_lbls)).sum()))               # nucleus volume

            sheet2_b.write(0, q + 1, float(idxs[q]))                                                                      # writing the id
            sheet2_b.write(1, q + 1, txt_pos[q]['centroid'][0])                                                           # x coordinate of nuclei centroid
            sheet2_b.write(2, q + 1, txt_pos[q]['centroid'][1])                                                           # y coordinate of nuclei centroid
            sheet2_b.write(3, q + 1, np.sum((nucs_dil == idxs[q]) * mtx_mature))                                          # number of spots
            sheet2_b.write(4, q + 1, float((nucs_dil == idxs[q]).sum() * steps))                                          # region nuclei volume
            sheet2_b.write(5, q + 1, float(((nucs_dil == idxs[q]) * np.sign(nucs_3d_det.nucs_lbls)).sum()))               # nucleus volume

            sheet3_b.write(0, q + 1, float(idxs[q]))                                                                      # writing the id
            sheet3_b.write(1, q + 1, txt_pos[q]['centroid'][0])                                                           # x coordinate of nuclei centroid
            sheet3_b.write(2, q + 1, txt_pos[q]['centroid'][1])                                                           # y coordinate of nuclei centroid
            sheet3_b.write(3, q + 1, np.sum((nucs_dil == idxs[q]) * mtx_mature))                                          # number of spots
            sheet3_b.write(4, q + 1, float((nucs_dil == idxs[q]).sum() * steps))                                          # region nuclei volume
            sheet3_b.write(5, q + 1, float(((nucs_dil == idxs[q]) * np.sign(nucs_3d_det.nucs_lbls)).sum()))               # nucleus volume

        pbar1.close()

        pbar2  =  ProgressBar(total1=idxs.size)
        pbar2.show()

        for j in range(idxs.size):
            # print(j)
            pbar2.update_progressbar(j)
            idxs_ins  =  np.unique((nucs_dil == idxs[j]) * mtx_mature_int * spts_segm)[1:]                              # tags for the spots inside the dilated nucleus (nuclear region)

            jj  =  -1           # this is just for paging purposes
            for jj in range(idxs_ins.size):
                a_jj  =  np.where(spts_mature[0, :] == idxs_ins[jj])[0]
                if a_jj.size > 0:
                    sheet1.write(7 + jj, j + 1, float(spts_mature[1, a_jj[0]]))                                         # intensity value of the spots
                    sheet3.write(7 + jj, j + 1, float(spts_mature[1, a_jj[0]]))

                    bkg_bff  =  (bkg_cages == spts_mature[0, a_jj[0]]) * raw_spts
                    bkg_bff  =  bkg_bff[bkg_bff != 0].mean()
                    sheet1_b.write(7 + jj, j + 1, float(spts_mature[1, a_jj[0]]) / bkg_bff)             # intensity value of the spots
                    sheet3_b.write(7 + jj, j + 1, float(spts_mature[1, a_jj[0]]) / bkg_bff)

            # this part is just to write the names of the files and the data: if statement is just for paging purposes

            if j == 0:
                try:
                    sheet1.write(7 + jj + 2, 0, "File Name")
                    sheet1.write(7 + jj + 3, 0, raw_data_fname)
                    sheet1.write(7 + jj + 5, 0, "date")
                    sheet1.write(7 + jj + 6, 0, str(datetime.date.today().day) + '/' + str(datetime.date.today().month) + '/' + str(datetime.date.today().year))
                    sheet1.write(7 + jj + 8, 0, "Software Version")
                    sheet1.write(7 + jj + 9, 0, soft_version)

                    sheet1_b.write(7 + jj + 2, 0, "File Name")
                    sheet1_b.write(7 + jj + 3, 0, raw_data_fname)
                    sheet1_b.write(7 + jj + 5, 0, "date")
                    sheet1_b.write(7 + jj + 6, 0, str(datetime.date.today().day) + '/' + str(datetime.date.today().month) + '/' + str(datetime.date.today().year))
                    sheet1_b.write(7 + jj + 8, 0, "Software Version")
                    sheet1_b.write(7 + jj + 9, 0, soft_version)

                except UnboundLocalError:
                    sheet1.write(13, 0, "File Name")
                    sheet1.write(14, 0, raw_data_fname)
                    sheet1.write(16, 0, "date")
                    sheet1.write(17, 0, str(datetime.date.today().day) + '/' + str(datetime.date.today().month) + '/' + str(datetime.date.today().year))
                    sheet1.write(19, 0, "Software Version")
                    sheet1.write(20, 0, soft_version)

                    sheet1_b.write(13, 0, "File Name")
                    sheet1_b.write(14, 0, raw_data_fname)
                    sheet1_b.write(16, 0, "date")
                    sheet1_b.write(17, 0, str(datetime.date.today().day) + '/' + str(datetime.date.today().month) + '/' + str(datetime.date.today().year))
                    sheet1_b.write(19, 0, "Software Version")
                    sheet1_b.write(20, 0, soft_version)

            idxs_exs  =  np.unique((nucs_dil == idxs[j]) * mtx_mature_ext * spts_segm)[1:]                              # tags for the spots inside the dilated nucleus (nuclear region)

            for kk in range(idxs_exs.size):
                a_kk  =  np.where(spts_mature[0, :] == idxs_exs[kk])[0]
                if a_kk.size > 0:
                    sheet2.write(7 + kk, j + 1, float(spts_mature[1, a_kk[0]]))                                         # intensity value of the spots
                    sheet3.write(7 + jj + kk + 1, j + 1, float(spts_mature[1, a_kk[0]]))

                    bkg_bff  =  (bkg_cages == spts_mature[0, a_kk[0]]) * raw_spts
                    bkg_bff  =  bkg_bff[bkg_bff != 0].mean()
                    sheet2_b.write(7 + kk, j + 1, float(spts_mature[1, a_kk[0]]) / bkg_bff)                                         # intensity value of the spots
                    sheet3_b.write(7 + jj + kk + 1, j + 1, float(spts_mature[1, a_kk[0]]) / bkg_bff)


        book.save(fwritename + '/smiFish_journal.xls')
        book_b.save(fwritename + '/smiFish_background_journal.xls')

        steps, xsize, ysize  =  spts_segm.shape
        spts_segm            =  spts_segm.reshape(spts_segm.size)
        spts_segm            = np.append(steps, spts_segm)
        spts_segm            =  np.append(xsize, spts_segm)
        spts_segm            =  np.append(ysize, spts_segm)
        spts_segm.astype("uint16").tofile(str(fwritename) + "/spots_segmented.bin")

        nucs_lbls  =  nucs_3d_det.nucs_lbls
        nucs_lbls  =  nucs_lbls.reshape(nucs_3d_det.nucs_lbls.size)
        nucs_lbls  =  np.append(steps, nucs_lbls)
        nucs_lbls  =  np.append(xsize, nucs_lbls)
        nucs_lbls  =  np.append(ysize, nucs_lbls)
        nucs_lbls.astype("uint16").tofile(str(fwritename) + "/nucs_segmented.bin")

        nucs_dil  =  nucs_dil.reshape(nucs_dil.size)
        nucs_dil  =  np.append(xsize, nucs_dil)
        nucs_dil  =  np.append(ysize, nucs_dil)
        nucs_dil.astype("uint16").tofile(str(fwritename) + "/nucs_dil.bin")

        pbar2.close()



class ProgressBar(QtGui.QWidget):
    """ProgressBar"""
    def __init__(self, parent=None, total1=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1  =  QtGui.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtGui.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        """function to update the progress bar"""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()
