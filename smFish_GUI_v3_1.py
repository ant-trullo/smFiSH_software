"""This is the main window of the software to analyze smFiSH and smiFiSH. 
This is version 3.1. Since late october 2019
"""


import sys
import os.path
from importlib import reload
import traceback
import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets, QtCore
from skimage.morphology import label
from skimage.measure import regionprops
from skimage.color import label2rgb
from scipy.ndimage.morphology import binary_dilation

import RawDataLoader
import NucsSegmenter
import NucsPileUp
import AnalysisSaver2
import AnalysisLoader
import LabelsModify
import SpotsDetection3D_RecursiveDoG
import SpotsIntensityRecord
import Spots3D_Segmentation_beta
import SaveSptsSmiFish
import BackgroundSptsSmiFish
import NucsEllipsoidFitting
import SaveNucsAnalysis
import ColocalizeSpots
import MatureNascentUtility
import RemoveCentralSpots



class MainWindow(QtWidgets.QMainWindow):
    """Main windows: coordinates all the actions, algorithms, visualization tools and analysis tools"""
    def __init__(self, parent=None):

        QtWidgets.QMainWindow.__init__(self, parent)

        widget  =  QtWidgets.QWidget(self)
        self.setCentralWidget(widget)

        load_data_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/load-hi.png'), "&Load data", self)
        load_data_action.setShortcut("Ctrl+L")
        load_data_action.setStatusTip("Load lsm files and the xls output of FQ")
        load_data_action.triggered.connect(self.load_data)

        save_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/save-md.png'), "&Save", self)
        save_action.setShortcut("Ctrl+S")
        save_action.setStatusTip("Save the analysis in a folder")
        save_action.triggered.connect(self.save_analysis)

        save_spts_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/save-md.png'), "&Save Spots", self)
        save_spts_action.setShortcut("Ctrl+P")
        save_spts_action.setStatusTip("Save the analysis in a folder")
        save_spts_action.triggered.connect(self.save_spots)

        save_nucs_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/save-md.png'), "&Save Nuclei", self)
        save_nucs_action.setShortcut("Ctrl+N")
        save_nucs_action.setStatusTip("Save the analysis in a folder")
        save_nucs_action.triggered.connect(self.save_nucs)

        load_analysis_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/load-hi.png'), "&Load Analysis", self)
        load_analysis_action.setShortcut("Ctrl+A")
        load_analysis_action.setStatusTip("Load previously done analysis")
        load_analysis_action.triggered.connect(self.load_analysis)

        crop_stack_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/crop.png'), "&Crop Stack", self)
        crop_stack_action.setShortcut("Ctrl+C")
        crop_stack_action.setStatusTip("Crop the stack before the analysis")
        crop_stack_action.triggered.connect(self.crop_stack)

        exit_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/exit.png'), "&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit application")
        exit_action.triggered.connect(self.close)

        coloc_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/colocalization_icon.png'), "&Colocalize", self)
        coloc_action.setShortcut("Ctrl+B")
        coloc_action.setStatusTip("Exit application")
        coloc_action.triggered.connect(self.colocalize)

        chop_stack_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/scissors.png'), "&Chop", self)
        chop_stack_action.setShortcut("Ctrl+k")
        chop_stack_action.setStatusTip("Chop Stack")
        chop_stack_action.triggered.connect(self.chop_stack)

        # convert_analysis_action  =  QtWidgets.QAction("&CovertSaveAnalysis", self)
        # convert_analysis_action.triggered.connect(self.convert_analysis)

        menubar   =  self.menuBar()

        file_menu  =  menubar.addMenu("&File")
        file_menu.addAction(load_data_action)
        file_menu.addAction(save_spts_action)
        file_menu.addAction(save_nucs_action)
        file_menu.addAction(save_action)
        file_menu.addAction(load_analysis_action)
        file_menu.addAction(exit_action)

        modify_menu  =  menubar.addMenu("&Modify")
        modify_menu.addAction(crop_stack_action)
        modify_menu.addAction(chop_stack_action)

        post_processing_menu  =  menubar.addMenu("&PostProcessing")
        post_processing_menu.addAction(coloc_action)

        temporary_menu  =  menubar.addMenu("&Temporary")
        temporary_menu.addAction(coloc_action)

        busy_lbl  =  QtWidgets.QLabel("Ready")
        busy_lbl.setStyleSheet("color: green")

        pixsize_x_lbl  =  QtWidgets.QLabel("pix size XY =;")

        pixsize_z_lbl  =  QtWidgets.QLabel("Z step =")

        bottom_labels_box  =  QtWidgets.QHBoxLayout()
        bottom_labels_box.addWidget(busy_lbl)
        bottom_labels_box.addStretch()
        bottom_labels_box.addWidget(pixsize_x_lbl)
        bottom_labels_box.addWidget(pixsize_z_lbl)

        fname_raw_lbl  =  QtWidgets.QLabel("File: ", self)
        fname_raw_lbl.setToolTip("Name of the file you are working on")

        tabs_tot    =  QtWidgets.QTabWidget()
        tab_spts_a  =  QtWidgets.QWidget()
        tab_spts_b  =  QtWidgets.QWidget()
        tab_nucs    =  QtWidgets.QWidget()

        tabs_tot.addTab(tab_spts_a, "Spots A")
        tabs_tot.addTab(tab_spts_b, "Spots B")
        tabs_tot.addTab(tab_nucs, "Nuclei")

        tabs_nucs  =  QtWidgets.QTabWidget()
        tab1_nucs  =  QtWidgets.QWidget()
        tab2_nucs  =  QtWidgets.QWidget()
        tab3_nucs  =  QtWidgets.QWidget()
        tab4_nucs  =  QtWidgets.QWidget()


        # ~~~~~~~~~~~ NUCLEI START ~~~~~~~~~~~~ #


        frame_nucs_raw  =  pg.ImageView(self)
        frame_nucs_raw.ui.roiBtn.hide()
        frame_nucs_raw.ui.menuBtn.hide()
        frame_nucs_raw.view.setXLink("FrameNucsRaw")
        frame_nucs_raw.view.setYLink("FrameNucsRaw")
        frame_nucs_raw.timeLine.sigPositionChanged.connect(self.update_frames_from_raw)

        frame_nucs2D_segm  =  pg.ImageView(self, name="FrameNucsRaw")
        frame_nucs2D_segm.getImageItem().mouseClickEvent  =  self.click
        frame_nucs2D_segm.ui.roiBtn.hide()
        frame_nucs2D_segm.ui.menuBtn.hide()
        frame_nucs2D_segm.timeLine.sigPositionChanged.connect(self.update_frames_from_nucs2d)

        frame_nucs3D_segm  =  pg.ImageView(self)
        frame_nucs3D_segm.ui.roiBtn.hide()
        frame_nucs3D_segm.ui.menuBtn.hide()
        frame_nucs3D_segm.view.setXLink("FrameNucsRaw")
        frame_nucs3D_segm.view.setYLink("FrameNucsRaw")
        frame_nucs3D_segm.timeLine.sigPositionChanged.connect(self.update_frames_from_nucs3d)

        frame_nucs_ellips  =  pg.ImageView(self)
        frame_nucs_ellips.ui.roiBtn.hide()
        frame_nucs_ellips.ui.menuBtn.hide()
        frame_nucs_ellips.view.setXLink("FrameNucsRaw")
        frame_nucs_ellips.view.setYLink("FrameNucsRaw")
        frame_nucs_ellips.timeLine.sigPositionChanged.connect(self.update_frames_from_ellips)

        frame_nucs_raw_box  =  QtWidgets.QHBoxLayout()
        frame_nucs_raw_box.addWidget(frame_nucs_raw)

        frame_nucs2D_segm_box  =  QtWidgets.QHBoxLayout()
        frame_nucs2D_segm_box.addWidget(frame_nucs2D_segm)

        frame_nucs3D_segm_box  =  QtWidgets.QHBoxLayout()
        frame_nucs3D_segm_box.addWidget(frame_nucs3D_segm)

        frame_nucs_ellips_box  =  QtWidgets.QHBoxLayout()
        frame_nucs_ellips_box.addWidget(frame_nucs_ellips)

        tab1_nucs.setLayout(frame_nucs_raw_box)
        tab2_nucs.setLayout(frame_nucs2D_segm_box)
        tab3_nucs.setLayout(frame_nucs3D_segm_box)
        tab4_nucs.setLayout(frame_nucs_ellips_box)

        tabs_nucs.addTab(tab1_nucs, "Raw Data")
        tabs_nucs.addTab(tab2_nucs, "Segmented 2D")
        tabs_nucs.addTab(tab3_nucs, "Segmented 3D")
        tabs_nucs.addTab(tab4_nucs, "Ellipsoids")

        nucs_2dsegm_btn  =  QtWidgets.QPushButton("N Segm", self)
        nucs_2dsegm_btn.clicked.connect(self.nucs_2dsegm)
        nucs_2dsegm_btn.setToolTip("Segmentation of nuclei z by z")
        nucs_2dsegm_btn.setFixedSize(110, 25)

        nucs2D_pileup_btn  =  QtWidgets.QPushButton("Pile Up", self)
        nucs2D_pileup_btn.clicked.connect(self.nucs2D_pileup)
        nucs2D_pileup_btn.setToolTip("Pile up nuclei slices in all the z frames")
        nucs2D_pileup_btn.setFixedSize(110, 25)

        nucs_fit_visual_btn  =  QtWidgets.QPushButton("Visualize Fit", self)
        nucs_fit_visual_btn.clicked.connect(self.nucs_fit_visual)
        nucs_fit_visual_btn.setToolTip("Visually check the ellipsoidal fitting")
        nucs_fit_visual_btn.setFixedSize(110, 25)

        nucs_ellipsoids_fit_btn  =  QtWidgets.QPushButton("Ellipsod Fit", self)
        nucs_ellipsoids_fit_btn.clicked.connect(self.nucs_ellipsoids_fit)
        nucs_ellipsoids_fit_btn.setToolTip("Fit 3D nuclei with ellipsoids")
        nucs_ellipsoids_fit_btn.setFixedSize(110, 25)

        shuffle_clrs_nucs_btn  =  QtWidgets.QPushButton("Shuffle Clrs", self)
        shuffle_clrs_nucs_btn.clicked.connect(self.shuffle_clrs_nucs)
        shuffle_clrs_nucs_btn.setToolTip("Shuffle Colors")
        shuffle_clrs_nucs_btn.setFixedSize(110, 25)

        man_active_toggle = QtWidgets.QCheckBox("Hand Cut", self)
        man_active_toggle.setFixedSize(110, 25)
        man_active_toggle.stateChanged.connect(self.man_active)

        manual_cut_btn  =  QtWidgets.QPushButton("Modify", self)
        manual_cut_btn.clicked.connect(self.manual_cut)
        manual_cut_btn.setToolTip("Manual cutting of nuclei (Ctrl+Suppr)")
        manual_cut_btn.setFixedSize(110, 25)
        manual_cut_btn.setEnabled(True)
        
        frame_numb_nucs_lbl  =  QtWidgets.QLabel("Frame numb    ", self)
        frame_numb_nucs_lbl.setFixedSize(130, 25)

        hor_line_one  =  QtWidgets.QFrame()
        hor_line_one.setFrameStyle(QtWidgets.QFrame.HLine)

        hor_line_two  =  QtWidgets.QFrame()
        hor_line_two.setFrameStyle(QtWidgets.QFrame.HLine)

        nucs_fitting_combo  =  QtWidgets.QComboBox(self)
        nucs_fitting_combo.addItem("Ellips")
        nucs_fitting_combo.addItem("Semi-Ellip")
        nucs_fitting_combo.setFixedSize(90, 25)
        nucs_fitting_combo.activated[str].connect(self.nucs_fitting_switch)

        keys_nucs  =  QtWidgets.QVBoxLayout()
        keys_nucs.addWidget(nucs_2dsegm_btn)
        keys_nucs.addWidget(man_active_toggle)
        keys_nucs.addWidget(manual_cut_btn)
        keys_nucs.addWidget(nucs2D_pileup_btn)
        keys_nucs.addWidget(hor_line_one)
        keys_nucs.addWidget(nucs_fitting_combo)
        keys_nucs.addWidget(nucs_ellipsoids_fit_btn)
        keys_nucs.addWidget(hor_line_two)
        keys_nucs.addWidget(nucs_fit_visual_btn)
        keys_nucs.addStretch()
        keys_nucs.addWidget(shuffle_clrs_nucs_btn)
        keys_nucs.addWidget(frame_numb_nucs_lbl)

        layout_nucs  =  QtWidgets.QHBoxLayout()
        layout_nucs.addWidget(tabs_nucs)
        layout_nucs.addLayout(keys_nucs)

        # ~~~~~~~~~~~ NUCLEI END ~~~~~~~~~~~~ #




        # ~~~~~~~~~~~ SPOT-A START ~~~~~~~~~~~~ #

        frame_spts_raw_a  =  pg.ImageView(self, name="FrameSpts1")
        frame_spts_raw_a.ui.roiBtn.hide()
        frame_spts_raw_a.ui.menuBtn.hide()
        frame_spts_raw_a.timeLine.sigPositionChanged.connect(self.update_frame_spts_segm_a)

        frame_spts_segm_a  =  pg.ImageView(self)
        frame_spts_segm_a.ui.roiBtn.hide()
        frame_spts_segm_a.ui.menuBtn.hide()
        frame_spts_segm_a.view.setXLink("FrameSpts1")
        frame_spts_segm_a.view.setYLink("FrameSpts1")
        frame_spts_segm_a.timeLine.sigPositionChanged.connect(self.update_frame_spts_raw_a)

        frame_spts_plot_a  =  pg.PlotWidget(self)
        frame_spts_plot_a.setFixedSize(140, 100)

        tabs_spts_a  =  QtWidgets.QTabWidget()
        tab1_spts_a  =  QtWidgets.QWidget()
        tab2_spts_a  =  QtWidgets.QWidget()

        frame1_box_spts_a  =  QtWidgets.QHBoxLayout()
        frame1_box_spts_a.addWidget(frame_spts_raw_a)

        frame2_box_spts_a  =  QtWidgets.QHBoxLayout()
        frame2_box_spts_a.addWidget(frame_spts_segm_a)

        tab1_spts_a.setLayout(frame1_box_spts_a)
        tab2_spts_a.setLayout(frame2_box_spts_a)

        tabs_spts_a.addTab(tab1_spts_a, "Raw Data")
        tabs_spts_a.addTab(tab2_spts_a, "Segmented")

        thresh_lbl_a  =  QtWidgets.QLabel("Thresholds", self)
        thresh_lbl_a.setFixedSize(110, 25)

        min_thr_lbl_a  =  QtWidgets.QLabel("Min Thr", self)
        min_thr_lbl_a.setFixedSize(60, 25)

        min_thr_edt_a  =  QtWidgets.QLineEdit(self)
        min_thr_edt_a.textChanged[str].connect(self.min_thr_var_a)
        min_thr_edt_a.setToolTip("Sets the minimum value for the thresholding; suggested value 0")
        min_thr_edt_a.setFixedSize(35, 25)

        min_thr_lbl_edt_a  =  QtWidgets.QHBoxLayout()
        min_thr_lbl_edt_a.addWidget(min_thr_lbl_a)
        min_thr_lbl_edt_a.addWidget(min_thr_edt_a)

        max_thr_lbl_a  =  QtWidgets.QLabel("Max Thr", self)
        max_thr_lbl_a.setFixedSize(60, 25)

        max_thr_edt_a  =  QtWidgets.QLineEdit(self)
        max_thr_edt_a.textChanged[str].connect(self.max_thr_var_a)
        max_thr_edt_a.setToolTip("Sets the maximum value for the thresholding; suggested value 9")
        max_thr_edt_a.setFixedSize(35, 25)

        max_thr_lbl_edt_a  =  QtWidgets.QHBoxLayout()
        max_thr_lbl_edt_a.addWidget(max_thr_lbl_a)
        max_thr_lbl_edt_a.addWidget(max_thr_edt_a)

        steps_thr_lbl_a  =  QtWidgets.QLabel("Steps", self)
        steps_thr_lbl_a.setFixedSize(60, 25)

        steps_thr_edt_a  =  QtWidgets.QLineEdit(self)
        steps_thr_edt_a.textChanged[str].connect(self.steps_thr_var_a)
        steps_thr_edt_a.setToolTip("Sets the number of different thresholds between min and max to test")
        steps_thr_edt_a.setFixedSize(35, 25)

        steps_thr_lbl_edt_a  =  QtWidgets.QHBoxLayout()
        steps_thr_lbl_edt_a.addWidget(steps_thr_lbl_a)
        steps_thr_lbl_edt_a.addWidget(steps_thr_edt_a)

        spts_3d_detect_btn  =  QtWidgets.QPushButton("Thr Study", self)
        spts_3d_detect_btn.clicked.connect(self.spts_3d_detect)
        spts_3d_detect_btn.setToolTip("Segmentation study over several thresholds")
        spts_3d_detect_btn.setFixedSize(130, 25)

        select_thr_btn_a  =  QtWidgets.QPushButton("Select Thr", self)
        select_thr_btn_a.clicked.connect(self.select_thr)
        select_thr_btn_a.setToolTip("Select the threshold value pointed by the blue line")
        select_thr_btn_a.setFixedSize(130, 25)

        segm_cluster_combo_a  =  QtWidgets.QComboBox(self)
        segm_cluster_combo_a.addItem("Segm")
        segm_cluster_combo_a.addItem("Cluster")
        segm_cluster_combo_a.activated[str].connect(self.segm_cluster_switcher_a)
        segm_cluster_combo_a.setCurrentIndex(0)
        segm_cluster_combo_a.setFixedSize(65, 25)

        segment_spts_btn_a  =  QtWidgets.QPushButton("Segment", self)
        segment_spts_btn_a.clicked.connect(self.segment_spts)
        segment_spts_btn_a.setToolTip("Segment thresholded spots")
        segment_spts_btn_a.setFixedSize(130, 25)

        shuffle_clrs_spts_btn_a  =  QtWidgets.QPushButton("Shuffle Clrs", self)
        shuffle_clrs_spts_btn_a.clicked.connect(self.shuffle_clrs_spts)
        shuffle_clrs_spts_btn_a.setToolTip("Shuffle spots colors")
        shuffle_clrs_spts_btn_a.setFixedSize(130, 25)

        sld_thr_a  =  QtWidgets.QScrollBar(QtCore.Qt.Horizontal, self)
        sld_thr_a.valueChanged.connect(self.sld_thr_update_a)

        sld_val_lbl_a  =  QtWidgets.QLabel('Slider val 0', self)
        sld_val_lbl_a.setFixedSize(130, 25)

        numb_detected_spts_lbl_a  =  QtWidgets.QLabel("Detected ", self)
        numb_detected_spts_lbl_a.setFixedSize(130, 25)

        intensity_study_btn_a  =  QtWidgets.QPushButton("Spts Int", self)
        intensity_study_btn_a.clicked.connect(self.intensity_study)
        intensity_study_btn_a.setToolTip("Study Spots Intensity")
        intensity_study_btn_a.setFixedSize(130, 25)

        avint_thr_lbl_a  =  QtWidgets.QLabel("Av Int Thr", self)
        avint_thr_lbl_a.setFixedSize(60, 25)

        avint_thr_edt_a  =  QtWidgets.QLineEdit(self)
        avint_thr_edt_a.textChanged[str].connect(self.int_thr_var_a)
        avint_thr_edt_a.setToolTip("Average Intensity Threshold for Spots")
        avint_thr_edt_a.setFixedSize(40, 25)
        avint_thr_edt_a.setText(str(0))

        avint_box_a  =  QtWidgets.QHBoxLayout()
        avint_box_a.addWidget(avint_thr_lbl_a)
        avint_box_a.addWidget(avint_thr_edt_a)

        vol_thr_lbl_a  =  QtWidgets.QLabel("Vol Thr", self)
        vol_thr_lbl_a.setFixedSize(60, 25)

        vol_thr_edt_a  =  QtWidgets.QLineEdit(self)
        vol_thr_edt_a.textChanged[str].connect(self.vol_thr_var_a)
        vol_thr_edt_a.setToolTip("Volume Threshold for Spots")
        vol_thr_edt_a.setFixedSize(40, 25)
#        vol_thr_edt_a.setText(str(0))

        vol_box_a  =  QtWidgets.QHBoxLayout()
        vol_box_a.addWidget(vol_thr_lbl_a)
        vol_box_a.addWidget(vol_thr_edt_a)

        sph_thr_lbl_a  =  QtWidgets.QLabel("Sph Thr", self)
        sph_thr_lbl_a.setFixedSize(60, 25)

        sph_thr_edt_a  =  QtWidgets.QLineEdit(self)
        sph_thr_edt_a.textChanged[str].connect(self.sph_thr_var_a)
        sph_thr_edt_a.setToolTip("Sphericity Threshold for Spots (suggested value 9)")
        sph_thr_edt_a.setFixedSize(40, 25)
#        sph_thr_edt_a.setText(str(0))

        sph_box_a  =  QtWidgets.QHBoxLayout()
        sph_box_a.addWidget(sph_thr_lbl_a)
        sph_box_a.addWidget(sph_thr_edt_a)

        min_numb_z_lbl_a  =  QtWidgets.QLabel("Min # of z", self)
        min_numb_z_lbl_a.setFixedSize(60, 25)

        min_numb_z_edt_a  =  QtWidgets.QLineEdit(self)
        min_numb_z_edt_a.textChanged[str].connect(self.min_numb_z_var_a)
        min_numb_z_edt_a.setToolTip("Minimum number of z plains for a spot")
        min_numb_z_edt_a.setFixedSize(40, 25)
#        min_numb_z_edt_a.setText(str(0))

        min_numb_z_box_a  =  QtWidgets.QHBoxLayout()
        min_numb_z_box_a.addWidget(min_numb_z_lbl_a)
        min_numb_z_box_a.addWidget(min_numb_z_edt_a)

        popup_frame3_btn_a  =  QtWidgets.QPushButton("PopUp", self)
        popup_frame3_btn_a.clicked.connect(self.popup_frame3_a)
        popup_frame3_btn_a.setToolTip("Pop up the frame to enlarge it")
        popup_frame3_btn_a.setFixedSize(65, 20)
        
        frame_numb_spts_a_lbl  =  QtWidgets.QLabel("Frame numb    ", self)
        frame_numb_spts_a_lbl.setFixedSize(130, 25)

        keys_spts_a  =  QtWidgets.QVBoxLayout()
        keys_spts_a.addStretch()
        keys_spts_a.addWidget(thresh_lbl_a)
        keys_spts_a.addLayout(min_thr_lbl_edt_a)
        keys_spts_a.addLayout(max_thr_lbl_edt_a)
        keys_spts_a.addLayout(steps_thr_lbl_edt_a)
        keys_spts_a.addWidget(spts_3d_detect_btn)
        keys_spts_a.addStretch()
        keys_spts_a.addWidget(popup_frame3_btn_a)
        keys_spts_a.addWidget(frame_spts_plot_a)
        keys_spts_a.addWidget(sld_thr_a)
        keys_spts_a.addWidget(sld_val_lbl_a)
        keys_spts_a.addWidget(select_thr_btn_a)

        keys_spts_a.addLayout(avint_box_a)
        keys_spts_a.addLayout(vol_box_a)
        keys_spts_a.addLayout(min_numb_z_box_a)
        keys_spts_a.addWidget(segm_cluster_combo_a)
        keys_spts_a.addLayout(sph_box_a)
        keys_spts_a.addWidget(segment_spts_btn_a)
        keys_spts_a.addWidget(numb_detected_spts_lbl_a)
        keys_spts_a.addWidget(intensity_study_btn_a)
        keys_spts_a.addStretch()
        keys_spts_a.addWidget(shuffle_clrs_spts_btn_a)
        keys_spts_a.addWidget(frame_numb_spts_a_lbl)

        layout_spts_a  =  QtWidgets.QHBoxLayout()
        layout_spts_a.addWidget(tabs_spts_a)
        layout_spts_a.addLayout(keys_spts_a)


# ~~~~~~~~~~~ SPOT-A END ~~~~~~~~~~~~ #


# ~~~~~~~~~~~ SPOT-B START ~~~~~~~~~~~~ #

        frame_spts_raw_b  =  pg.ImageView(self, name='FrameSpts1')
        frame_spts_raw_b.ui.roiBtn.hide()
        frame_spts_raw_b.ui.menuBtn.hide()
        frame_spts_raw_b.timeLine.sigPositionChanged.connect(self.update_frame_spts_segm_b)

        frame_spts_segm_b  =  pg.ImageView(self)
        frame_spts_segm_b.ui.roiBtn.hide()
        frame_spts_segm_b.ui.menuBtn.hide()
        frame_spts_segm_b.view.setXLink('FrameSpts1')
        frame_spts_segm_b.view.setYLink('FrameSpts1')
        frame_spts_segm_b.timeLine.sigPositionChanged.connect(self.update_frame_spts_raw_b)

        frame_spts_plot_b  =  pg.PlotWidget(self)
        frame_spts_plot_b.setFixedSize(140, 100)

        tabs_spts_b  =  QtWidgets.QTabWidget()
        tab1_spts_b  =  QtWidgets.QWidget()
        tab2_spts_b  =  QtWidgets.QWidget()

        frame1_box_spts_b  =  QtWidgets.QHBoxLayout()
        frame1_box_spts_b.addWidget(frame_spts_raw_b)

        frame2_box_spts_b  =  QtWidgets.QHBoxLayout()
        frame2_box_spts_b.addWidget(frame_spts_segm_b)

        tab1_spts_b.setLayout(frame1_box_spts_b)
        tab2_spts_b.setLayout(frame2_box_spts_b)

        tabs_spts_b.addTab(tab1_spts_b, "Raw Data")
        tabs_spts_b.addTab(tab2_spts_b, "Segmented")

        thresh_lbl_b  =  QtWidgets.QLabel('Thresholds', self)
        thresh_lbl_b.setFixedSize(110, 25)

        min_thr_lbl_b  =  QtWidgets.QLabel('Min Thr', self)
        min_thr_lbl_b.setFixedSize(60, 25)

        min_thr_edt_b  =  QtWidgets.QLineEdit(self)
        min_thr_edt_b.textChanged[str].connect(self.min_thr_var_b)
        min_thr_edt_b.setToolTip('Sets the minimum value for the thresholding; suggested value 1')
        min_thr_edt_b.setFixedSize(35, 25)

        min_thr_lbl_edt_b  =  QtWidgets.QHBoxLayout()
        min_thr_lbl_edt_b.addWidget(min_thr_lbl_b)
        min_thr_lbl_edt_b.addWidget(min_thr_edt_b)

        max_thr_lbl_b  =  QtWidgets.QLabel('Max Thr', self)
        max_thr_lbl_b.setFixedSize(60, 25)

        max_thr_edt_b  =  QtWidgets.QLineEdit(self)
        max_thr_edt_b.textChanged[str].connect(self.max_thr_var_b)
        max_thr_edt_b.setToolTip('Sets the maximum value for the thresholding; suggested value 9')
        max_thr_edt_b.setFixedSize(35, 25)

        max_thr_lbl_edt_b  =  QtWidgets.QHBoxLayout()
        max_thr_lbl_edt_b.addWidget(max_thr_lbl_b)
        max_thr_lbl_edt_b.addWidget(max_thr_edt_b)

        steps_thr_lbl_b  =  QtWidgets.QLabel('Steps', self)
        steps_thr_lbl_b.setFixedSize(60, 25)

        steps_thr_edt_b  =  QtWidgets.QLineEdit(self)
        steps_thr_edt_b.textChanged[str].connect(self.steps_thr_var_b)
        steps_thr_edt_b.setToolTip('Sets the number of different thresholds between min and max to test')
        steps_thr_edt_b.setFixedSize(35, 25)

        steps_thr_lbl_edt_b  =  QtWidgets.QHBoxLayout()
        steps_thr_lbl_edt_b.addWidget(steps_thr_lbl_b)
        steps_thr_lbl_edt_b.addWidget(steps_thr_edt_b)

        spts_3d_detect_btn_b  =  QtWidgets.QPushButton("Thr Study", self)
        spts_3d_detect_btn_b.clicked.connect(self.spts_3d_detect)
        spts_3d_detect_btn_b.setToolTip('Segmentation study over several thresholds')
        spts_3d_detect_btn_b.setFixedSize(130, 25)

        select_thr_btn_b  =  QtWidgets.QPushButton("Select Thr", self)
        select_thr_btn_b.clicked.connect(self.select_thr)
        select_thr_btn_b.setToolTip('Select the threshold value pointed by the blue line')
        select_thr_btn_b.setFixedSize(130, 25)

        segm_cluster_combo_b  =  QtWidgets.QComboBox(self)
        segm_cluster_combo_b.addItem("Segm")
        segm_cluster_combo_b.addItem("Cluster")
        segm_cluster_combo_b.activated[str].connect(self.segm_cluster_switcher_b)
        segm_cluster_combo_b.setCurrentIndex(0)
        segm_cluster_combo_b.setFixedSize(65, 25)

        segment_spts_btn_b  =  QtWidgets.QPushButton("Segment", self)
        segment_spts_btn_b.clicked.connect(self.segment_spts)
        segment_spts_btn_b.setToolTip('Segment thresholded spots')
        segment_spts_btn_b.setFixedSize(130, 25)

        shuffle_clrs_spts_btn_b  =  QtWidgets.QPushButton("Shuffle Clrs", self)
        shuffle_clrs_spts_btn_b.clicked.connect(self.shuffle_clrs_spts)
        shuffle_clrs_spts_btn_b.setToolTip('Shuffle spots colors')
        shuffle_clrs_spts_btn_b.setFixedSize(130, 25)

        sld_thr_b  =  QtWidgets.QScrollBar(QtCore.Qt.Horizontal, self)
        sld_thr_b.valueChanged.connect(self.sld_thr_update_b)

        sld_val_lbl_b  =  QtWidgets.QLabel('Slider val 0', self)
        sld_val_lbl_b.setFixedSize(130, 25)

        numb_detected_spts_lbl_b  =  QtWidgets.QLabel('Detected ', self)
        numb_detected_spts_lbl_b.setFixedSize(130, 25)

        intensity_study_btn_b  =  QtWidgets.QPushButton("Spts Int", self)
        intensity_study_btn_b.clicked.connect(self.intensity_study)
        intensity_study_btn_b.setToolTip('Study Spots Intensity')
        intensity_study_btn_b.setFixedSize(130, 25)

        avint_thr_lbl_b  =  QtWidgets.QLabel('Av Int Thr', self)
        avint_thr_lbl_b.setFixedSize(60, 25)

        avint_thr_edt_b  =  QtWidgets.QLineEdit(self)
        avint_thr_edt_b.textChanged[str].connect(self.int_thr_var_b)
        avint_thr_edt_b.setToolTip('Average Intensity Threshold for Spots')
        avint_thr_edt_b.setFixedSize(40, 25)
        avint_thr_edt_b.setText(str(0))

        avint_box_b  =  QtWidgets.QHBoxLayout()
        avint_box_b.addWidget(avint_thr_lbl_b)
        avint_box_b.addWidget(avint_thr_edt_b)

        vol_thr_lbl_b  =  QtWidgets.QLabel('Vol Thr', self)
        vol_thr_lbl_b.setFixedSize(60, 25)

        vol_thr_edt_b  =  QtWidgets.QLineEdit(self)
        vol_thr_edt_b.textChanged[str].connect(self.vol_thr_var_b)
        vol_thr_edt_b.setToolTip('Volume Threshold for Spots')
        vol_thr_edt_b.setFixedSize(40, 25)
#        vol_thr_edt_b.setText(str(0))

        vol_box_b  =  QtWidgets.QHBoxLayout()
        vol_box_b.addWidget(vol_thr_lbl_b)
        vol_box_b.addWidget(vol_thr_edt_b)

        sph_thr_lbl_b  =  QtWidgets.QLabel("Sph Thr", self)
        sph_thr_lbl_b.setFixedSize(60, 25)

        sph_thr_edt_b  =  QtWidgets.QLineEdit(self)
        sph_thr_edt_b.textChanged[str].connect(self.sph_thr_var_b)
        sph_thr_edt_b.setToolTip("Sphericity Threshold for Spots (suggested value 9)")
        sph_thr_edt_b.setFixedSize(40, 25)
#        sph_thr_edt_b.setText(str(0))

        sph_box_b  =  QtWidgets.QHBoxLayout()
        sph_box_b.addWidget(sph_thr_lbl_b)
        sph_box_b.addWidget(sph_thr_edt_b)

        min_numb_z_lbl_b  =  QtWidgets.QLabel('Min # of z', self)
        min_numb_z_lbl_b.setFixedSize(60, 25)

        min_numb_z_edt_b  =  QtWidgets.QLineEdit(self)
        min_numb_z_edt_b.textChanged[str].connect(self.min_numb_z_var_b)
        min_numb_z_edt_b.setToolTip('Minimum number of z plains for a spot')
        min_numb_z_edt_b.setFixedSize(40, 25)
#        min_numb_z_edt_b.setText(str(0))

        min_numb_z_box_b  =  QtWidgets.QHBoxLayout()
        min_numb_z_box_b.addWidget(min_numb_z_lbl_b)
        min_numb_z_box_b.addWidget(min_numb_z_edt_b)

        popup_frame3_btn_b  =  QtWidgets.QPushButton("PopUp", self)
        popup_frame3_btn_b.clicked.connect(self.popup_frame3_b)
        popup_frame3_btn_b.setToolTip('Pop up the frame to enlarge it')
        popup_frame3_btn_b.setFixedSize(65, 20)
        
        frame_numb_spts_b_lbl  =  QtWidgets.QLabel("Frame numb    ", self)
        frame_numb_spts_b_lbl.setFixedSize(130, 25)

        keys_spts_b  =  QtWidgets.QVBoxLayout()
        keys_spts_b.addStretch()
        keys_spts_b.addWidget(thresh_lbl_b)
        keys_spts_b.addLayout(min_thr_lbl_edt_b)
        keys_spts_b.addLayout(max_thr_lbl_edt_b)
        keys_spts_b.addLayout(steps_thr_lbl_edt_b)
        keys_spts_b.addWidget(spts_3d_detect_btn_b)
        keys_spts_b.addStretch()
        keys_spts_b.addWidget(popup_frame3_btn_b)
        keys_spts_b.addWidget(frame_spts_plot_b)
        keys_spts_b.addWidget(sld_thr_b)
        keys_spts_b.addWidget(sld_val_lbl_b)
        keys_spts_b.addWidget(select_thr_btn_b)
        keys_spts_b.addLayout(avint_box_b)
        keys_spts_b.addLayout(vol_box_b)
        keys_spts_b.addLayout(min_numb_z_box_b)
        keys_spts_b.addWidget(segm_cluster_combo_b)
        keys_spts_b.addLayout(sph_box_b)
        keys_spts_b.addWidget(segment_spts_btn_b)
        keys_spts_b.addWidget(numb_detected_spts_lbl_b)
        keys_spts_b.addWidget(intensity_study_btn_b)
        keys_spts_b.addStretch()
        keys_spts_b.addWidget(shuffle_clrs_spts_btn_b)
        keys_spts_b.addWidget(frame_numb_spts_b_lbl)

        layout_spts_b  =  QtWidgets.QHBoxLayout()
        layout_spts_b.addWidget(tabs_spts_b)
        layout_spts_b.addLayout(keys_spts_b)

        # ~~~~~~~~~~~ SPOT-B END ~~~~~~~~~~~~ #


        tab_nucs.setLayout(layout_nucs)
        tab_spts_a.setLayout(layout_spts_a)
        tab_spts_b.setLayout(layout_spts_b)

        layout  =  QtWidgets.QVBoxLayout(widget)
        layout.addWidget(fname_raw_lbl)
        layout.addWidget(tabs_tot)
        layout.addLayout(bottom_labels_box)

        mycmap  =  np.fromfile("mycmap.bin", "uint16").reshape((10000, 3)) / 255.0
        self.colors4map  =  []
        for k in range(mycmap.shape[0]):
            self.colors4map.append(mycmap[k, :])
        self.colors4map     =  self.colors4map + self.colors4map + self.colors4map + self.colors4map + self.colors4map + self.colors4map
        self.colors4map[0]  =  np.array([0, 0, 0])
        self.colors4map[1]  =  np.array([1, 1, 1])

        self.frame_nucs2D_segm  =  frame_nucs2D_segm
        self.frame_nucs_raw     =  frame_nucs_raw
        self.frame_nucs3D_segm  =  frame_nucs3D_segm
        self.frame_nucs_ellips  =  frame_nucs_ellips
        self.fname_raw_lbl      =  fname_raw_lbl
        self.frame_spts_raw_a   =  frame_spts_raw_a
        self.frame_spts_segm_a  =  frame_spts_segm_a
        self.frame_spts_plot_a  =  frame_spts_plot_a
        self.frame_spts_raw_b   =  frame_spts_raw_b
        self.frame_spts_segm_b  =  frame_spts_segm_b
        self.frame_spts_plot_b  =  frame_spts_plot_b
        self.tabs_tot           =  tabs_tot
        self.sld_thr_b          =  sld_thr_b
        self.sld_thr_a          =  sld_thr_a
        self.min_thr_edt_a      =  min_thr_edt_a
        self.min_thr_edt_b      =  min_thr_edt_b
        self.max_thr_edt_a      =  max_thr_edt_a
        self.max_thr_edt_b      =  max_thr_edt_b
        self.steps_thr_edt_a    =  steps_thr_edt_a
        self.steps_thr_edt_b    =  steps_thr_edt_b
        self.avint_thr_edt_a    =  avint_thr_edt_a
        self.avint_thr_edt_b    =  avint_thr_edt_b
        self.vol_thr_edt_a      =  vol_thr_edt_a
        self.vol_thr_edt_b      =  vol_thr_edt_b
        self.min_numb_z_edt_a   =  min_numb_z_edt_a
        self.min_numb_z_edt_b   =  min_numb_z_edt_b
        self.sph_thr_edt_a      =  sph_thr_edt_a
        self.sph_thr_edt_b      =  sph_thr_edt_b
        self.man_active_flag    =  0
        self.c_count            =  0
        self.frame_nucs_cntr    =  0
        self.int_thr_value_a    =  0
        self.int_thr_value_b    =  0
        self.rmv_central_flag   =  0
        self.mycmap             =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3))
        self.sph_thr_lbl_a      =  sph_thr_lbl_a
        self.sph_thr_lbl_b      =  sph_thr_lbl_b
        self.soft_version       =  "smFish_v3.1"
        self.nucs_fitting_flag  =  "Ellips"
        self.busy_lbl           =  busy_lbl
        self.pixsize_x_lbl      =  pixsize_x_lbl
        self.pixsize_z_lbl      =  pixsize_z_lbl

        self.nucs_fit_visual_btn       =  nucs_fit_visual_btn
        self.segment_spts_btn_a        =  segment_spts_btn_a
        self.segment_spts_btn_b        =  segment_spts_btn_b
        self.nucs_2dsegm_btn           =  nucs_2dsegm_btn
        self.manual_cut_btn            =  manual_cut_btn
        self.numb_detected_spts_lbl_a  =  numb_detected_spts_lbl_a
        self.numb_detected_spts_lbl_b  =  numb_detected_spts_lbl_b
        self.sld_val_lbl_a             =  sld_val_lbl_a
        self.sld_val_lbl_b             =  sld_val_lbl_b
        self.first_last_frame          =  np.array([0, 0])
        self.frame_numb_spts_a_lbl     =  frame_numb_spts_a_lbl
        self.frame_numb_spts_b_lbl     =  frame_numb_spts_b_lbl
        self.frame_numb_nucs_lbl       =  frame_numb_nucs_lbl
        self.nucs_fitting_combo        =  nucs_fitting_combo
        self.segm_cluster_combo_a      =  segm_cluster_combo_a
        self.segm_cluster_combo_b      =  segm_cluster_combo_b
        self.segm_clstr_flag_a         =  "s"
        self.segm_clstr_flag_b         =  "s"

        self.setGeometry(800, 100, 900, 800)
        self.setWindowTitle(self.soft_version)
        self.setWindowIcon(QtGui.QIcon('Icons/DrosophilaIcon.png'))
        self.show()


    def closeEvent(self, event):
        "Close the GUI, asking confirmation"
        quit_msg  =  "Are you sure you want to exit the program?"
        reply     =  QtWidgets.QMessageBox.question(self, 'Message', quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


    def busy_indicator(self):
        """Write a red text (BUSY) as a label on the GUI (bottom left)"""
        self.busy_lbl.setText("Busy")
        self.busy_lbl.setStyleSheet('color: red')


    def ready_indicator(self):
        """Write a green text (READY) as a label on the GUI (bottom left)"""
        self.busy_lbl.setText("Ready")
        self.busy_lbl.setStyleSheet('color: green')


    def update_frames_from_raw(self):
        """Index changes of the raw frame activates all the other frames in which an image is set"""
        self.frame_numb_nucs_lbl.setText("Frame numb " + str(self.frame_nucs_raw.currentIndex))
        try:
            self.frame_nucs2D_segm.setCurrentIndex(self.frame_nucs_raw.currentIndex)
        except AttributeError:
            pass

        try:
            self.frame_nucs3D_segm.setCurrentIndex(self.frame_nucs_raw.currentIndex)
        except AttributeError:
            pass

        try:
            self.frame_nucs_ellips.setCurrentIndex(self.frame_nucs_raw.currentIndex)
        except AttributeError:
            pass


    def update_frames_from_nucs2d(self):
        """Frame index changes of the nucs2d activates 'self.update_frames_from_raw' """
        self.frame_nucs_raw.setCurrentIndex(self.frame_nucs2D_segm.currentIndex)


    def update_frames_from_nucs3d(self):
        """Frame index changes of the nucs3d activates 'self.update_frames_from_raw' """
        self.frame_nucs_raw.setCurrentIndex(self.frame_nucs3D_segm.currentIndex)


    def update_frames_from_ellips(self):
        """Frame index changes of the ellips activates 'self.update_frames_from_raw' """
        self.frame_nucs_raw.setCurrentIndex(self.frame_nucs_ellips.currentIndex)


    def update_frame_spts_segm_a(self):
        """Frame index changes in raw spots a activates changes in segmented spots a and updates frame number edit_box"""
        self.frame_numb_spts_a_lbl.setText("Frame numb " + str(self.frame_spts_raw_a.currentIndex))
        try:
            self.frame_spts_segm_a.setCurrentIndex(self.frame_spts_raw_a.currentIndex)
        except AttributeError:
            pass


    def update_frame_spts_raw_a(self):
        """Frame index changes in segmented spots a activates changes in raw spots a"""
        self.frame_spts_raw_a.setCurrentIndex(self.frame_spts_segm_a.currentIndex)


    def update_frame_spts_segm_b(self):
        """Frame index changes in raw spots b activates changes in segmented spots b and updates frame number edit_box"""
        self.frame_numb_spts_b_lbl.setText("Frame numb " + str(self.frame_spts_raw_b.currentIndex))
        try:
            self.frame_spts_segm_b.setCurrentIndex(self.frame_spts_raw_b.currentIndex)
        except AttributeError:
            pass


    def update_frame_spts_raw_b(self):
        """Frame index changes in segmented spots b activates changes in raw spots b"""
        self.frame_spts_raw_b.setCurrentIndex(self.frame_spts_segm_b.currentIndex)


    def segm_cluster_switcher_a(self, text):
        """Switch between clustering and segmentation for nuclei, channel a"""
        if text == "Segm":
            self.sph_thr_lbl_a.setText("Sph Thr")
            self.segment_spts_btn_a.setText("Segment")
            self.segm_clstr_flag_a  =  "s"
        else:
            self.sph_thr_lbl_a.setText("Dist Clst")
            self.segment_spts_btn_a.setText("Cluster")
            self.segm_clstr_flag_a  =  "c"


    def segm_cluster_switcher_b(self, text):
        """Switch between clustering and segmentation for nuclei, channel b"""
        if text == "Segm":
            self.sph_thr_lbl_b.setText("Sph Thr")
            self.segment_spts_btn_b.setText("Segment")
            self.segm_clstr_flag_b  =  "s"
        else:
            self.sph_thr_lbl_b.setText("Dist Clst")
            self.segment_spts_btn_b.setText("Cluster")
            self.segm_clstr_flag_b  =  "c"


    def crop_stack(self):
        """Call CroppingTool Tool"""
        if self.tabs_tot.currentIndex() == 0:
            self.mpp1  =  CroppingTool(self.raw_data.spts_a)

        if self.tabs_tot.currentIndex() == 1:
            self.mpp1  =  CroppingTool(self.raw_data.spts_b)

        self.mpp1.show()
        self.mpp1.procStart.connect(self.crop_tool_sgnl)


    def crop_tool_sgnl(self):
        """Crop raw data stack following the prescription of the Cropping Tool"""
        pts  =  self.mpp1.roi.parentBounds()
        x0   =  np.round(np.max([0, pts.x()])).astype(np.int)
        y0   =  np.round(np.max([0, pts.y()])).astype(np.int)
        x1   =  np.round(np.min([pts.x() + pts.width(), self.raw_data.spts_a.shape[1]])).astype(np.int)
        y1   =  np.round(np.min([pts.y() + pts.height(), self.raw_data.spts_a.shape[2]])).astype(np.int)

        self.raw_data.spts_a  =  self.raw_data.spts_a[:, x0:x1, y0:y1]
        if self.chs_spts_nucs.size > 2:
            self.raw_data.spts_b  =  self.raw_data.spts_b[:, x0:x1, y0:y1]
        self.raw_data.nucs    =  self.raw_data.nucs[:, x0:x1, y0:y1]
        self.frame_spts_raw_a.setImage(self.raw_data.spts_a)
        self.frame_nucs_raw.setImage(self.raw_data.nucs)
        if self.chs_spts_nucs.size > 2:
            self.frame_spts_raw_b.setImage(self.raw_data.spts_b)

        self.mpp1.close()


    def nucs_fitting_switch(self, text):
        """Switch between different nuclei fitting methods"""
        self.nucs_fitting_flag  =  text


    def load_data(self):
        """Load and visualize raw data to start the analysis"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        self.chs_spts_nucs  =  SetColorChannel.getChannels()

        try:
            self.raw_data_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data file to analyze", filter='*.czi')[0])
            # self.raw_data_fname  =  '/home/atrullo/Dropbox/smiFish/JD_data/JD129(m142het_GFP_Sun_Armadillo)_cc14_4_Airyscan Processing.czi'
            self.fname_raw_lbl.setText(self.raw_data_fname)

            self.raw_data  =  RawDataLoader.RawDataLoader(self.raw_data_fname, self.chs_spts_nucs)
            self.pixsize_x_lbl.setText("pix size XY = " + str(self.raw_data.pix_sizeX) + "µm;")
            self.pixsize_z_lbl.setText("Z step = " + str(self.raw_data.pix_sizeZ) + "µm")

            try:
                self.frame_spts_raw_a.setImage(self.raw_data.spts_a)
            except AttributeError:
                pass
            try:
                self.frame_nucs_raw.setImage(self.raw_data.nucs)
            except AttributeError:
                pass
            try:
                self.frame_spts_raw_b.setImage(self.raw_data.spts_b)
            except AttributeError:
                pass

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def load_analysis(self):
        """Load an already done analysis"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        self.analysis_folder  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select the Directory with the Analysis"))
        # self.analysis_folder  =  '/home/atrullo/Desktop/DeleteSmiFishTest'
        self.raw_data_fname   =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data file of the analysis", self.analysis_folder, filter='*.tif *.czi *.lsm')[0])
        # self.raw_data_fname   =  '/home/atrullo/Dropbox/smiFish/smFISH_Matthieu/C10DsRedlessxYw_emb2_Border_Out.czi'
        self.fname_raw_lbl.setText("File: " + self.raw_data_fname)

        try:
            self.chs_spts_nucs     =  np.fromfile(self.analysis_folder + '/chs_spts_nucs.bin', 'uint16')
            self.raw_data          =  RawDataLoader.RawDataLoader(self.raw_data_fname, self.chs_spts_nucs)
            self.first_last_frame  =  np.load(self.analysis_folder + '/first_last_frame.npy')
            if self.first_last_frame.size == 3:

                if self.first_last_frame[2] == 0:
                    self.raw_data.nucs    =  self.raw_data.nucs[self.first_last_frame[0]:self.first_last_frame[1]]
                    self.raw_data.spts_a  =  self.raw_data.spts_a[self.first_last_frame[0]:self.first_last_frame[1]]
                    self.raw_data.spts_b  =  self.raw_data.spts_b[self.first_last_frame[0]:self.first_last_frame[1]]

                if self.first_last_frame[2] == 1:
                    self.rmv_central_flag  =  1

            self.frame_nucs_raw.setImage(self.raw_data.nucs)
            self.frame_spts_raw_a.setImage(self.raw_data.spts_a)
            self.frame_spts_raw_b.setImage(self.raw_data.spts_b)
            self.pixsize_x_lbl.setText("pix size XY = " + str(self.raw_data.pix_sizeX) + "µm;")
            self.pixsize_z_lbl.setText("Z step = " + str(self.raw_data.pix_sizeZ) + "µm")


            if os.path.isfile(self.analysis_folder + '/bkg_cages_a.npy'):
                self.spts_thr_study_a   =  AnalysisLoader.AnalysisLoaderSptsStudy(self.analysis_folder, "_a")
                self.spts_segm_a        =  AnalysisLoader.AnalysisLoaderSptsSegm(self.analysis_folder, "_a").spts_lbls.astype(np.int32)
                self.spts_detect_a      =  AnalysisLoader.AnalysisLoaderSptsDetected(self.analysis_folder, "_a").spts_detected
                self.frame_nucs_cntr_a  =  self.frame_spts_segm_a.currentIndex
                self.frame_spts_segm_a.setImage(self.spts_segm_a, levels=(0, self.spts_segm_a.max()))
                self.spts_segm2show_a    =  label2rgb(self.spts_segm_a, bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
                self.frame_spts_segm_a.setImage(self.spts_segm2show_a)
                params                   =  AnalysisLoader.AnalysisLoaderParams(self.analysis_folder, "_a")
                self.min_thr_value_a     =  params.min_thr_value
                self.max_thr_value_a     =  params.max_thr_value
                self.steps_thr_value_a   =  params.steps_thr_value
                self.sld_thr_a_current   =  params.sld_thr_current
                self.av_int_thr_a        =  params.av_int_thr
                self.vol_thr_value_a     =  params.vol_thr
                self.min_numb_z_value_a  =  params.min_numb_z_value
                self.sph_thr_value_a     =  params.sph_thr
                self.segm_clstr_flag_a   =  params.segm_clstr_flag
                if self.segm_clstr_flag_a == "c":
                    self.segm_cluster_combo_a.setCurrentIndex(1)
                    self.segm_cluster_switcher_a("c")
                self.thr_vals_a          =  np.linspace(self.min_thr_value_a, self.max_thr_value_a, self.steps_thr_value_a)
                self.min_thr_edt_a.setText(str(self.min_thr_value_a))
                self.max_thr_edt_a.setText(str(self.max_thr_value_a))
                self.steps_thr_edt_a.setText(str(self.steps_thr_value_a))
                self.avint_thr_edt_a.setText(str(self.av_int_thr_a))
                self.vol_thr_edt_a.setText(str(self.vol_thr_value_a))
                self.min_numb_z_edt_a.setText(str(self.min_numb_z_value_a))
                self.sph_thr_edt_a.setText(str(self.sph_thr_value_a))
                self.frame_spts_plot_a.clear()

                try:
                    self.frame_spts_plot_a.removeItem(self.roi_sld_a)
                except AttributeError:
                    pass

                x_dict_lbls  =  {0: ''}
                for k in range(1, self.thr_vals_a.size):
                    x_dict_lbls.update({k: ''})
                y_dict_lbls  =  x_dict_lbls
                self.frame_spts_plot_a.getPlotItem().axes['bottom']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.frame_spts_plot_a.getPlotItem().axes['left']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.frame_spts_plot_a.plot(self.thr_vals_a, self.spts_thr_study_a.spts_num, symbol='x', pen='r')
                self.roi_sld_a  =  pg.LineSegmentROI([[self.thr_vals_a[0], 0], [self.thr_vals_a[0], self.spts_thr_study_a.spts_num.max()]], pen='b')
                self.frame_spts_plot_a.addItem(self.roi_sld_a)
                self.sld_thr_a.setMaximum(self.steps_thr_value_a - 1)
                self.sld_thr_a.setValue(self.sld_thr_a_current)
                self.numb_detected_spts_lbl_a.setText("Detected " + str(len(regionprops(self.spts_segm_a))))
                self.info_mature_a   =  np.load(self.analysis_folder + '/info_mature_a.npy')
                self.info_nascent_a  =  np.load(self.analysis_folder + '/info_nascent_a.npy')


            if os.path.isfile(self.analysis_folder + '/bkg_cages_b.npy'):
                self.spts_thr_study_b   =  AnalysisLoader.AnalysisLoaderSptsStudy(self.analysis_folder, "_b")
                self.spts_segm_b        =  AnalysisLoader.AnalysisLoaderSptsSegm(self.analysis_folder, "_b").spts_lbls
                self.spts_detect_b      =  AnalysisLoader.AnalysisLoaderSptsDetected(self.analysis_folder, "_b").spts_detected
                self.frame_nucs_cntr_b  =  self.frame_spts_segm_b.currentIndex
                self.frame_spts_segm_b.setImage(self.spts_segm_b, levels=(0, self.spts_segm_b.max()))
                self.spts_segm2show_b  =  label2rgb(self.spts_segm_b, bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
                self.frame_spts_segm_b.setImage(self.spts_segm2show_b)
                params                   =  AnalysisLoader.AnalysisLoaderParams(self.analysis_folder, "_b")
                self.min_thr_value_b     =  params.min_thr_value
                self.max_thr_value_b     =  params.max_thr_value
                self.steps_thr_value_b   =  params.steps_thr_value
                self.sld_thr_b_current   =  params.sld_thr_current
                self.av_int_thr_b        =  params.av_int_thr
                self.vol_thr_value_b     =  params.vol_thr
                self.min_numb_z_value_b  =  params.min_numb_z_value
                self.sph_thr_value_b     =  params.sph_thr
                self.thr_vals_b          =  np.linspace(self.min_thr_value_b, self.max_thr_value_b, self.steps_thr_value_b)
                self.segm_clstr_flag_b   =  params.segm_clstr_flag
                if self.segm_clstr_flag_b == "c":
                    self.segm_cluster_combo_b.setCurrentIndex(1)
                    self.segm_cluster_switcher_b("c")                    
                self.min_thr_edt_b.setText(str(self.min_thr_value_b))
                self.max_thr_edt_b.setText(str(self.max_thr_value_b))
                self.steps_thr_edt_b.setText(str(self.steps_thr_value_b))
                self.avint_thr_edt_b.setText(str(self.av_int_thr_b))
                self.vol_thr_edt_b.setText(str(self.vol_thr_value_b))
                self.min_numb_z_edt_b.setText(str(self.min_numb_z_value_b))
                self.sph_thr_edt_b.setText(str(self.sph_thr_value_b))
                self.frame_spts_plot_b.clear()

                try:
                    self.frame_spts_plot_b.removeItem(self.roi_sld_b)
                except AttributeError:
                    pass

                x_dict_lbls  =  {0: ''}
                for k in range(1, self.thr_vals_b.size):
                    x_dict_lbls.update({k: ''})
                y_dict_lbls  =  x_dict_lbls
                self.frame_spts_plot_b.getPlotItem().axes['bottom']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.frame_spts_plot_b.getPlotItem().axes['left']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.frame_spts_plot_b.plot(self.thr_vals_b, self.spts_thr_study_b.spts_num, symbol='x', pen='r')
                self.roi_sld_b  =  pg.LineSegmentROI([[self.thr_vals_b[0], 0], [self.thr_vals_b[0], self.spts_thr_study_b.spts_num.max()]], pen='b')
                self.frame_spts_plot_b.addItem(self.roi_sld_b)
                self.sld_thr_b.setMaximum(self.steps_thr_value_b - 1)
                self.sld_thr_b.setValue(self.sld_thr_b_current)
                self.numb_detected_spts_lbl_b.setText("Detected " + str(len(regionprops(self.spts_segm_b))))
                self.info_mature_b   =  np.load(self.analysis_folder + '/info_mature_b.npy')
                self.info_nascent_b  =  np.load(self.analysis_folder + '/info_nascent_b.npy')

            if os.path.isfile(self.analysis_folder + '/nucs_ellips.npy'):
                nucs_bff  =  AnalysisLoader.AnalysisLoaderNucs(self.analysis_folder)
                self.nucs_2d_det        =  nucs_bff.nucs_2d_det
                self.nucs_3d_det        =  nucs_bff.nucs_3d_det
                self.nucs_ellips        =  nucs_bff.nucs_ellips
                self.nucs_fitting_flag  =  str(np.loadtxt(self.analysis_folder + '/NucsFittingInfo.txt', dtype=str))
                if self.nucs_fitting_flag == "Ellips":
                    self.nucs_fitting_combo.setCurrentIndex(0)
                if self.nucs_fitting_flag == "Semi-Ellip":
                    self.nucs_fitting_combo.setCurrentIndex(1)

                self.frame_nucs2D_segm.setImage(self.nucs_2d_det)
                self.rnd_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_2d_det.max()), color=self.colors4map)
                self.frame_nucs2D_segm.setColorMap(self.rnd_cmap)

                self.frame_nucs3D_segm.setImage(self.nucs_3d_det)
                self.rnd_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_3d_det.max()), color=self.colors4map)
                self.frame_nucs3D_segm.setColorMap(self.rnd_cmap)

                self.frame_nucs_ellips.setImage(self.nucs_ellips)
                self.rnd_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_ellips.max()), color=self.colors4map)
                self.frame_nucs_ellips.setColorMap(self.rnd_cmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def man_active(self, state):
        """Activate manual corrections for segmented nuclei"""
        if state == QtCore.Qt.Checked:
            self.man_active_flag  =  1
            self.manual_cut_btn.setEnabled(True)
            self.nucs_2dsegm_btn.setEnabled(False)
        else:
            self.man_active_flag  =  0
            self.manual_cut_btn.setEnabled(False)
            self.nucs_2dsegm_btn.setEnabled(True)


    def click(self, event):
        """Put a segment on the nuclei image to cut or merge"""
        event.accept()
        pos        =  event.pos()
        modifiers  =  QtWidgets.QApplication.keyboardModifiers()

        if self.man_active_flag == 1:
            if modifiers  ==  QtCore.Qt.ShiftModifier:
                if self.c_count - 2 * int(self.c_count / 2) == 0:
                    self.pos1  =  pos
                else:
                    self.roi  =  pg.LineSegmentROI([self.pos1, pos], pen='r')
                    self.frame_nucs2D_segm.addItem(self.roi)

                self.c_count  +=  1


    def manual_cut(self):
        """Manual nuclei cutting"""
        cif           =  self.frame_nucs2D_segm.currentIndex
        hh            =  self.frame_nucs2D_segm.view.viewRange()

        pp       =  self.roi.getHandles()
        pp       =  [self.roi.mapToItem(self.frame_nucs2D_segm.imageItem, p.pos()) for p in pp]
        end_pts  =  np.array([[int(pp[0].x()), int(pp[0].y())], [int(pp[1].x()), int(pp[1].y())]])

        self.bufframe                   =  np.copy(self.nucs_2d_det[cif, :, :])
        self.nucs_2d_det[cif, :, :]     =  LabelsModify.LabelsModify(self.nucs_2d_det[cif, :, :], end_pts).labels_fin
        self.frame_nucs2D_segm.setImage(self.nucs_2d_det)
        self.frame_nucs2D_segm.setCurrentIndex(cif)
        self.frame_nucs2D_segm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_nucs2D_segm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
        self.frame_nucs2D_segm.removeItem(self.roi)


    def keyPressEvent(self, event):
        """Undo command for nuclei manual modification"""
        if event.key() == (QtCore.Qt.ControlModifier and Qt.Key_Z):
            cif                                       =  self.frame_nucs2D_segm.currentIndex
            hh                                        =  self.frame_nucs2D_segm.view.viewRange()
            self.nucs_2d_det[cif, :, :]     =  self.bufframe
            self.frame_nucs2D_segm.setImage(self.nucs_2d_det)
            self.frame_nucs2D_segm.setCurrentIndex(cif)
            self.frame_nucs2D_segm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
            self.frame_nucs2D_segm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

        if event.key() == (QtCore.Qt.ShiftModifier and Qt.Key_Delete):
            self.manual_cut()


    def shuffle_clrs_nucs(self):
        """Shuffle nuclei colors"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:

            colors_bff  =  self.colors4map[1:]
            np.random.shuffle(colors_bff)
            self.colors4map[1:]  =  colors_bff
            mycmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_2d_det.max()), color=self.colors4map)
            self.frame_nucs2D_segm.setColorMap(mycmap)
            self.frame_nucs3D_segm.setColorMap(mycmap)
            self.frame_nucs_ellips.setColorMap(mycmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def nucs_2dsegm(self):
        """Nuclei segmentation frame by frame"""
        reload(NucsSegmenter)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:

            self.nucs_2d_det  =  NucsSegmenter.NucsSegmenter(self.raw_data.nucs).nucs_lbls
            self.frame_nucs2D_segm.setImage(self.nucs_2d_det)
            self.rnd_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_2d_det.max()), color=self.colors4map)
            self.frame_nucs2D_segm.setColorMap(self.rnd_cmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def nucs2D_pileup(self):
        """Pile up segmented 2D nuclei"""
        reload(NucsPileUp)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:

            self.nucs_3d_det  =  NucsPileUp.NucsPileUp(self.nucs_2d_det).nucs_lbls_piled
            self.frame_nucs3D_segm.setImage(self.nucs_3d_det)
            self.rnd_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_2d_det.max()), color=self.colors4map)
            self.frame_nucs3D_segm.setColorMap(self.rnd_cmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def nucs_ellipsoids_fit(self):
        """Pile up segmented 2D nuclei"""
        reload(NucsEllipsoidFitting)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            if self.nucs_fitting_flag == "Semi-Ellip":
                self.nucs_ellips  =  NucsEllipsoidFitting.NucsEllipsoidFitting(np.copy(self.nucs_3d_det), "Semi-Ellip").nucs_elips
                print("Semi")
            else:
                self.nucs_ellips  =  NucsEllipsoidFitting.NucsEllipsoidFitting(np.copy(self.nucs_3d_det)).nucs_elips
                print("Ellips")
            self.frame_nucs_ellips.setImage(self.nucs_ellips)
            self.rnd_cmap  =  pg.ColorMap(np.linspace(0, 1, self.nucs_2d_det.max()), color=self.colors4map)
            self.frame_nucs_ellips.setColorMap(self.rnd_cmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def nucs_fit_visual(self):
        """Popup tool to visualize nuclei ellipsoidal fitting"""
        self.mpp5  =  NucsFittingVisual(self.raw_data.nucs / self.raw_data.nucs.max(), np.sign(self.nucs_ellips))
        self.mpp5.show()


    def save_analysis(self):
        """Save the analysis"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:

            folder  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select Directory"))
            reload(AnalysisSaver2)
            AnalysisSaver2.AnalysisSaver2(folder, self.soft_version, self.raw_data.spts_a, self.nucs_area.nucs_dil, self.nucs_3d_det, self.spts_segm.spts_lbls, self.spts_mature, self.bkg_cages, self.raw_data_fname, self.int_thr_value_a)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def min_thr_var_a(self, text):
        """Minimum value of the threshold for A spots"""
        self.min_thr_value_a  =  np.float(text)


    def max_thr_var_a(self, text):
        """Maximum value of the threshold for A spots"""
        self.max_thr_value_a  =  np.float(text)


    def int_thr_var_a(self, text):
        """Intensity threshold for spots A"""
        self.int_thr_value_a  =  np.float(text)


    def vol_thr_var_a(self, text):
        """Intensity threshold for spots A"""
        self.vol_thr_value_a  =  np.int(text)


    def steps_thr_var_a(self, text):
        """Number of steps for the preliminary threshold study for spots A"""
        self.steps_thr_value_a  =  np.int(text)


    def min_thr_var_b(self, text):
        """Minimum value of the threshold for B spots"""
        self.min_thr_value_b  =  np.float(text)


    def max_thr_var_b(self, text):
        """Maximum value of the threshold for B spots"""
        self.max_thr_value_b  =  np.float(text)


    def int_thr_var_b(self, text):
        """Intensity threshold for spots B"""
        self.int_thr_value_b  =  np.float(text)


    def vol_thr_var_b(self, text):
        """Volume threshold for spots A"""
        self.vol_thr_value_b  =  np.int(text)


    def sph_thr_var_a(self, text):
        """Sphericity threshold for spots A"""
        self.sph_thr_value_a  =  np.int(text)


    def sph_thr_var_b(self, text):
        """Sphericity threshold for spots A"""
        self.sph_thr_value_b  =  np.int(text)


    def steps_thr_var_b(self, text):
        """Number of steps for the preliminary threshold study for spots B"""
        self.steps_thr_value_b  =  np.int(text)


    def min_numb_z_var_a(self, text):
        """Number of z plains the spot must be present in"""
        self.min_numb_z_value_a  =  np.int(text)


    def min_numb_z_var_b(self, text):
        """Number of z plains the spot must be present in"""
        self.min_numb_z_value_b  =  np.int(text)


    def spts_3d_detect(self):
        """3D spots detection"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        if self.tabs_tot.currentIndex() == 0:
            try:

                self.thr_vals_a  =  np.linspace(self.min_thr_value_a, self.max_thr_value_a, int(self.steps_thr_value_a))
                self.sld_thr_a.setMaximum(self.steps_thr_value_a - 1)
                reload(SpotsDetection3D_RecursiveDoG)
                # self.spts_thr_study  =  SpotsDetection3D_Recursive3.SpotsDetection3D_Recursive3(self.raw_data.spts, self.kern_gauss_value, self.thr_vals)
                self.spts_thr_study_a  =  SpotsDetection3D_RecursiveDoG.SpotsDetection3D_RecursiveDoG(self.raw_data.spts_a, self.thr_vals_a)

                self.frame_spts_plot_a.clear()
                try:
                    self.frame_spts_plot_a.removeItem(self.roi_sld_a)
                except AttributeError:
                    pass

                self.frame_spts_plot_a.plot(self.thr_vals_a, self.spts_thr_study_a.spts_num, symbol='x', pen='r')
                x_dict_lbls  =  {0: ''}
                for k in range(1, self.thr_vals_a.size):
                    x_dict_lbls.update({k: ''})
                y_dict_lbls  =  x_dict_lbls
                self.frame_spts_plot_a.getPlotItem().axes['bottom']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.frame_spts_plot_a.getPlotItem().axes['left']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.roi_sld_a  =  pg.LineSegmentROI([[self.thr_vals_a[0], 0], [self.thr_vals_a[0], self.spts_thr_study_a.spts_num.max()]], pen='b')
                self.frame_spts_plot_a.addItem(self.roi_sld_a)
                self.sld_thr_a.setValue(0)

            except Exception:
                traceback.print_exc()

        if self.tabs_tot.currentIndex() == 1:

            try:

                self.thr_vals_b  =  np.linspace(self.min_thr_value_b, self.max_thr_value_b, int(self.steps_thr_value_b))
                self.sld_thr_b.setMaximum(self.steps_thr_value_b - 1)
                reload(SpotsDetection3D_RecursiveDoG)
                # self.spts_thr_study  =  SpotsDetection3D_Recursive3.SpotsDetection3D_Recursive3(self.raw_data.spts, self.kern_gauss_value, self.thr_vals)
                self.spts_thr_study_b  =  SpotsDetection3D_RecursiveDoG.SpotsDetection3D_RecursiveDoG(self.raw_data.spts_b, self.thr_vals_b)

                self.frame_spts_plot_b.clear()
                try:
                    self.frame_spts_plot_b.removeItem(self.roi_sld_b)
                except AttributeError:
                    pass

                self.frame_spts_plot_b.plot(self.thr_vals_b, self.spts_thr_study_b.spts_num, symbol='x', pen='r')
                x_dict_lbls  =  {0: ''}
                for k in range(1, self.thr_vals_b.size):
                    x_dict_lbls.update({k: ''})
                y_dict_lbls  =  x_dict_lbls
                self.frame_spts_plot_b.getPlotItem().axes['bottom']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.frame_spts_plot_b.getPlotItem().axes['left']['item'].setTicks([list(x_dict_lbls.items()), list(y_dict_lbls.items())])
                self.roi_sld_b  =  pg.LineSegmentROI([[self.thr_vals_b[0], 0], [self.thr_vals_b[0], self.spts_thr_study_b.spts_num.max()]], pen='b')
                self.frame_spts_plot_b.addItem(self.roi_sld_b)
                self.sld_thr_b.setValue(0)

            except Exception:
                traceback.print_exc()

        self.ready_indicator()


    def select_thr(self):
        """Select the threshold to use for the analysis"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        if self.tabs_tot.currentIndex() == 0:

            try:

                self.frame_nucs_cntr_a  =  self.frame_spts_segm_a.currentIndex
                hh                      =  self.frame_spts_segm_a.view.viewRange()
                reload(SpotsDetection3D_RecursiveDoG)
                self.spts_detect_a      =  SpotsDetection3D_RecursiveDoG.SpotsDetection3D_SingleThr(self.spts_thr_study_a.spts_g, self.spts_thr_study_a.mu, self.spts_thr_study_a.sigma, self.thr_vals_a[self.sld_thr_a.value()]).spts_detect
                if self.rmv_central_flag  == 1:
                    self.spts_detect_a  =  RemoveCentralSpots.RemoveCentralSpots(self.spts_detect_a, self.first_last_frame).spts_selected  
                self.frame_spts_segm_a.setImage(np.sign(self.spts_detect_a))
                self.frame_spts_segm_a.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
                self.frame_spts_segm_a.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
                self.frame_spts_segm_a.setCurrentIndex(self.frame_nucs_cntr_a)

            except Exception:
                traceback.print_exc()

        if self.tabs_tot.currentIndex() == 1:

            try:

                self.frame_nucs_cntr_b  =  self.frame_spts_segm_b.currentIndex
                hh                      =  self.frame_spts_segm_b.view.viewRange()
                reload(SpotsDetection3D_RecursiveDoG)
                self.spts_detect_b      =  SpotsDetection3D_RecursiveDoG.SpotsDetection3D_SingleThr(self.spts_thr_study_b.spts_g, self.spts_thr_study_b.mu, self.spts_thr_study_b.sigma, self.thr_vals_b[self.sld_thr_b.value()]).spts_detect
                if self.rmv_central_flag  == 1:
                    self.spts_detect_b  =  RemoveCentralSpots.RemoveCentralSpots(self.spts_detect_b, self.first_last_frame).spts_selected  
                self.frame_spts_segm_b.setImage(np.sign(self.spts_detect_b))
                self.frame_spts_segm_b.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
                self.frame_spts_segm_b.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
                self.frame_spts_segm_b.setCurrentIndex(self.frame_nucs_cntr_b)

            except Exception:
                traceback.print_exc()

        self.ready_indicator()


    def segment_spts(self):
        """Spots 3D segmentation"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        reload(Spots3D_Segmentation_beta)

        if self.tabs_tot.currentIndex() == 0:

            try:

                hh                    =  self.frame_spts_segm_a.view.viewRange()
                self.frame_nucs_cntr  =  self.frame_spts_segm_a.currentIndex
                # self.spts_segm        =  Spots3D_Segmentation.Spots3D_Segmentation(self.spts_detect.spts_detect)   # need to define int_thr from the backgroung image
                self.spts_segm_a      =  Spots3D_Segmentation_beta.Spots3D_Segmentation(self.spts_detect_a, self.raw_data.spts_a, self.int_thr_value_a, self.vol_thr_value_a, self.min_numb_z_value_a, self.sph_thr_value_a, self.segm_clstr_flag_a).spts_lbls   # need to define int_thr from the backgroung image
                if self.rmv_central_flag  == 1:
                    self.spts_segm_a  =  RemoveCentralSpots.RemoveCentralSegmentedSpots(self.spts_segm_a, self.first_last_frame).spts_selected  
                self.frame_spts_segm_a.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
                self.frame_spts_segm_a.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
                self.spts_segm2show_a  =  label2rgb(self.spts_segm_a, bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
                self.frame_spts_segm_a.setImage(self.spts_segm2show_a)

                self.frame_spts_segm_a.setCurrentIndex(self.frame_nucs_cntr)
                self.numb_detected_spts_lbl_a.setText("Detected " + str(len(regionprops(self.spts_segm_a))))

            except Exception:
                traceback.print_exc()

        if self.tabs_tot.currentIndex() == 1:

            try:

                hh                      =  self.frame_spts_segm_b.view.viewRange()
                self.frame_nucs_cntr_b  =  self.frame_spts_segm_b.currentIndex
                self.spts_segm_b        =  Spots3D_Segmentation_beta.Spots3D_Segmentation(self.spts_detect_b, self.int_thr_value_b, self.vol_thr_value_b, self.min_numb_z_value_b, self.sph_thr_value_b, self.segm_clstr_flag_b).spts_lbls   # need to define int_thr from the backgroung image
                # self.spts_segm        =  Spots3D_Segmentation.Spots3D_Segmentation(self.spts_detect.spts_detect)   # need to define int_thr from the backgroung image
                if self.rmv_central_flag  == 1:
                    self.spts_segm_b  =  RemoveCentralSpots.RemoveCentralSegmentedSpots(self.spts_segm_b, self.first_last_frame).spts_selected  
                self.frame_spts_segm_b.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
                self.frame_spts_segm_b.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
                self.spts_segm2show_b  =  label2rgb(self.spts_segm_b, bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
                self.frame_spts_segm_b.setImage(self.spts_segm2show_b)
                self.numb_detected_spts_lbl_b.setText("Detected " + str(len(regionprops(self.spts_segm_b))))

            except Exception:
                traceback.print_exc()

        self.ready_indicator()


    def sld_thr_update_a(self):
        """Update the roi threshold position in the plot and the threshold label for a spots"""
        self.roi_sld_a.setPos([self.thr_vals_a[self.sld_thr_a.value()] - self.thr_vals_a[0], 0], update=True)
        self.sld_val_lbl_a.setText("Slider val " + str(self.sld_thr_a.value()))


    def sld_thr_update_b(self):
        """Update the roi threshold position in the plot and the threshold label for b spots"""
        self.roi_sld_b.setPos([self.thr_vals_b[self.sld_thr_b.value()] - self.thr_vals_b[0], 0], update=True)
        self.sld_val_lbl_b.setText("Slider val " + str(self.sld_thr_b.value()))


    def save_spots(self):
        """Save spots analysis (nuclei are not included)"""
        folder2write  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select or Define a Folder to Store the analysis results"))
        self.chs_spts_nucs.astype('uint16').tofile(folder2write + '/chs_spts_nucs.bin')
        np.save(folder2write + '/first_last_frame.npy', self.first_last_frame)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        reload(BackgroundSptsSmiFish)
        reload(SaveSptsSmiFish)

        if self.tabs_tot.currentIndex() == 0:
            try:
                bkg_cages  =  BackgroundSptsSmiFish.BackgroundSptsSmiFish(self.spts_segm_a).cages
                SaveSptsSmiFish.SaveSptsSmiFish(folder2write, self.soft_version, self.raw_data_fname, self.raw_data.spts_a, self.spts_thr_study_a, self.spts_segm_a, bkg_cages, self.info_mature_a, self.info_nascent_a, self.thr_vals_a, self.int_thr_value_a, self.sld_thr_a.value(), self.vol_thr_value_a, self.min_numb_z_value_a, self.sph_thr_value_a, self.segm_clstr_flag_a, "_a")

            except Exception:
                traceback.print_exc()

        if self.tabs_tot.currentIndex() == 1:
            try:
                bkg_cages  =  BackgroundSptsSmiFish.BackgroundSptsSmiFish(self.spts_segm_b).cages
                SaveSptsSmiFish.SaveSptsSmiFish(folder2write, self.soft_version, self.raw_data_fname, self.raw_data.spts_b, self.spts_thr_study_b, self.spts_segm_b, bkg_cages, self.info_mature_b, self.info_nascent_b, self.thr_vals_b, self.int_thr_value_b, self.sld_thr_b.value(), self.vol_thr_value_b, self.min_numb_z_value_b, self.sph_thr_value_b, self.segm_clstr_flag_b, "_b")

            except Exception:
                traceback.print_exc()

        self.ready_indicator()


    def save_nucs(self):
        """Save nuclei analysis (spots are not included)"""
        reload(SaveNucsAnalysis)
        folder2write  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select or Define a Folder to Store the analysis results"))
        self.chs_spts_nucs.astype('uint16').tofile(folder2write + '/chs_spts_nucs.bin')
        np.save(folder2write + '/first_last_frame.npy', self.first_last_frame)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:
            SaveNucsAnalysis.SaveNucsAnalysis(folder2write, self.nucs_2d_det, self.nucs_3d_det, self.nucs_ellips)
            np.savetxt(folder2write + '/NucsFittingInfo.txt', [self.nucs_fitting_flag], fmt='%s')

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def shuffle_clrs_spts(self):
        """Shuffle colors of the spots"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:
            np.random.shuffle(self.mycmap)
            if self.tabs_tot.currentIndex() == 0:

                hh                      =  self.frame_spts_segm_a.view.viewRange()
                self.frame_nucs_cntr_a  =  self.frame_spts_segm_a.currentIndex
                self.spts_segm2show_a   =  label2rgb(self.spts_segm_a, bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
                self.frame_spts_segm_a.setImage(self.spts_segm2show_a)
                self.frame_spts_segm_a.setCurrentIndex(self.frame_nucs_cntr_a)
                self.frame_spts_segm_a.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
                self.frame_spts_segm_a.view.setYRange(hh[1][0], hh[1][1], padding=.0002)


            if self.tabs_tot.currentIndex() == 1:

                hh                      =  self.frame_spts_segm_b.view.viewRange()
                self.frame_nucs_cntr_b  =  self.frame_spts_segm_b.currentIndex
                self.spts_segm2show_b   =  label2rgb(self.spts_segm_b, bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
                self.frame_spts_segm_b.setImage(self.spts_segm2show_b)
                self.frame_spts_segm_b.setCurrentIndex(self.frame_nucs_cntr_b)
                self.frame_spts_segm_b.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
                self.frame_spts_segm_b.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def intensity_study(self):
        """Collect and organize all the info of the spots"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()     # slkdfjslkd
        reload(SpotsIntensityRecord)

        if self.tabs_tot.currentIndex() == 0:
            try:

                self.spts_ints_ctrs_a  =  SpotsIntensityRecord.SpotsIntensityRecord(self.spts_segm_a, self.raw_data.spts_a).spts_ints_ctrs
                self.mpp1_a            =  MatureNascent(self.spts_ints_ctrs_a, self.spts_segm_a, self.raw_data.nucs_mip)
                self.mpp1_a.show()
                self.mpp1_a.procStart.connect(self.maturenascent_thr_a)
                
            except Exception:
                traceback.print_exc()

        if self.tabs_tot.currentIndex() == 1:
            try:

                self.spts_ints_ctrs_b  =  SpotsIntensityRecord.SpotsIntensityRecord(self.spts_segm_b, self.raw_data.spts_b).spts_ints_ctrs
                self.mpp1_b            =  MatureNascent(self.spts_ints_ctrs_b, self.spts_segm_b, self.raw_data.nucs_mip)
                self.mpp1_b.show()
                self.mpp1_b.procStart.connect(self.maturenascent_thr_b)

            except Exception:
                traceback.print_exc()

        self.ready_indicator()


    def maturenascent_thr_a(self, message):
        """Input the threshold mature/nascent for the detected objects"""
        print(message)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:
            j_mt  =  np.where(self.spts_ints_ctrs_a[2, :] <  message)
            j_ns  =  np.where(self.spts_ints_ctrs_a[2, :] >= message)

            self.info_mature_a   =  np.delete(self.spts_ints_ctrs_a, j_ns, axis=1)
            self.info_nascent_a  =  np.delete(self.spts_ints_ctrs_a, j_mt, axis=1)

            spts_mature_a   =  np.copy(self.spts_segm_a)
            spts_nascent_a  =  np.zeros((self.spts_segm_a.shape), dtype=np.int)

            spts_nas  =  np.zeros(self.spts_segm_a.shape)
            for k in self.info_nascent_a[0, :]:
                spts_nas  +=  k * (self.spts_segm_a == k)

            spts_mature_a  *=  (1 - np.sign(spts_nas)).astype(bool)

            for kk in self.info_nascent_a[0, :]:
                sing_bff        =   binary_dilation((spts_nas == kk), iterations=1) * (1 - np.sign(spts_nascent_a))
                spts_nascent_a  +=  kk * sing_bff

            self.info_nascent_a  =  SpotsIntensityRecord.SpotsIntensityRecord(spts_nascent_a, self.raw_data.spts_a).spts_ints_ctrs

                
            self.frame_spts_segm_a.setImage(np.sign(spts_mature_a) + 2 * np.sign(spts_nascent_a))

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def maturenascent_thr_b(self, message):
        """Input the threshold mature/nascent for the detected objects"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        reload(SpotsIntensityRecord)
        
        try:

            j_mt  =  np.where(self.spts_ints_ctrs_b[2, :] <  message)
            j_ns  =  np.where(self.spts_ints_ctrs_b[2, :] >= message)

            self.info_mature_b   =  np.delete(self.spts_ints_ctrs_b, j_ns, axis=1)
            self.info_nascent_b  =  np.delete(self.spts_ints_ctrs_b, j_mt, axis=1)

            spts_mature_b   =  np.copy(self.spts_segm_b)
            spts_nascent_b  =  np.zeros((self.spts_segm_b.shape))

            spts_nas  =  np.zeros(self.spts_segm_b.shape)
            for k in self.info_nascent_b[0, :]:
                spts_nas  +=  k * (self.spts_segm_b == k)

            spts_mature_b  *=  (1 - np.sign(spts_nas)).astype(bool)

            for kk in self.info_nascent_b[0, :]:
                sing_bff         =  binary_dilation((spts_nas == kk), iterations=1) * (1 - np.sign(spts_nascent_b))
                spts_nascent_b  +=  kk * sing_bff

            self.info_nascent_b  =  SpotsIntensityRecord.SpotsIntensityRecord(spts_nascent_b, self.raw_data.spts_b).spts_ints_ctrs

            self.frame_spts_segm_b.setImage(np.sign(spts_mature_b) + 2 * np.sign(spts_nascent_b))

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def popup_frame3_a(self):
        """Call popup tool for better visualization of the spots curve"""
        self.mpp2_a  =  PopUpPlot(self.thr_vals_a, self.spts_thr_study_a.spts_num)
        self.mpp2_a.show()
        self.mpp2_a.procStart.connect(self.drive_sld_thr_a)


    def drive_sld_thr_a(self, message):
        """Update the spots curve in the main Window"""
        self.sld_thr_a.setValue(message)
        self.mpp2_a.close()
        del self.mpp2_a


    def popup_frame3_b(self):
        """Call popup tool for better visualization of the spots curve"""
        self.mpp2_b  =  PopUpPlot(self.thr_vals_b, self.spts_thr_study_b.spts_num)
        self.mpp2_b.show()
        self.mpp2_b.procStart.connect(self.drive_sld_thr_b)


    def drive_sld_thr_b(self, message):
        """Update the spots curve in the main Window"""
        self.sld_thr_b.setValue(message)
        self.mpp2_b.close()
        del self.mpp2_b


    def colocalize(self):
        """Activate tool for colocalization"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()

        try:

            post_proc_folder  =  str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select the Directory with the Analysis"))
            raw_data_fname    =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data file of the analyze", filter='*.czi')[0])
            # self.post_proc_folder  =  '/home/atrullo/Dropbox/smiFish//JD134_CJD_cc11Zoom4/'
            params            =  AnalysisLoader.AnalysisLoaderParams(post_proc_folder, "_b")
            self.mpp3         =  ColocalizeTool(post_proc_folder, raw_data_fname, self.soft_version, params.vol_thr, params.sph_thr)
            self.mpp3.show()

        except Exception:
            traceback.print_exc()

        self.ready_indicator()


    def chop_stack(self):
        """Call the popup tool to chop the stack"""
        self.mpp4  =  ChopStack(self.raw_data)
        self.mpp4.show()
        self.mpp4.procStart.connect(self.chop_stack_sgnl)


    def chop_stack_sgnl(self, message):
        """Chop raw data stack following the prescription of the ChopStack tool"""
        if message == 0:
            self.raw_data.nucs    =  self.raw_data.nucs[self.mpp4.first_last_frame[0]:self.mpp4.first_last_frame[1]]
            self.raw_data.spts_a  =  self.raw_data.spts_a[self.mpp4.first_last_frame[0]:self.mpp4.first_last_frame[1]]
            self.raw_data.spts_b  =  self.raw_data.spts_b[self.mpp4.first_last_frame[0]:self.mpp4.first_last_frame[1]]
            self.frame_spts_raw_a.setImage(self.raw_data.spts_a)
            self.frame_nucs_raw.setImage(self.raw_data.nucs)
            if self.chs_spts_nucs.size > 2:
                self.frame_spts_raw_b.setImage(self.raw_data.spts_b)
            self.first_last_frame  =  np.append(self.mpp4.first_last_frame, 0)

        else:
            self.rmv_central_flag  =  1
            self.first_last_frame  =  np.append(self.mpp4.first_last_frame, 1)


        self.mpp4.close()



class MatureNascent(QtWidgets.QWidget):
    """Tool to discriminate between TS and mature mRNA in terms of intensity"""
    procStart  =  QtCore.pyqtSignal(int)

    def __init__(self, spts_ints_ctrs, spts_lbls, nucs_mip):
        QtWidgets.QWidget.__init__(self)

        # if spts_ints_ctrs.shape[1] < 260:
        #     # mtx2thr  =  MatureNascentUtility.MatureNascentUtility(spts_lbls.astype(np.int32), spts_ints_ctrs[(0, 2), :].astype(np.int32))
        #     mtx2thr  =  MatureNascentUtilityPY.MatureNascent(spts_lbls.astype(np.int32), spts_ints_ctrs[(0, 2), :].astype(np.int32)).mtx2thr
        # else:
        #     cpu_ow   =  multiprocessing.cpu_count()
        #     t_chops  =  1 + spts_ints_ctrs.shape[1] // cpu_ow

        #     a  =  []
        #     for t in range(cpu_ow - 1):
        #         a.append(spts_ints_ctrs[:, t * t_chops:(t + 1) * t_chops])           # in the multiprocessing pool each core will work on a certain number of frames: here we chop the frames

        #     a.append(spts_ints_ctrs[:, (t + 1) * t_chops:])
        #     job_args  =  []
        #     for k in range(cpu_ow):
        #         job_args.append([spts_lbls.astype(np.int32), a[k][(0, 2), :].astype(np.int32)])
    
        #     pool     =  multiprocessing.Pool()
        #     results  =  pool.map(func_patch, job_args)
        #     pool.close()
        #     mtx2thr  =  results[0].mtx2thr
        #     for p in range(1, cpu_ow):
        #         mtx2thr  +=  results[p].mtx2thr

        mtx2thr  =  MatureNascentUtility.MatureNascentUtility(spts_lbls.astype(np.int32), spts_ints_ctrs[(0, 2), :].astype(np.int32))
        
        lbls2show           =  np.zeros(spts_lbls.shape[1:] + (3,))
        lbls2show[:, :, 0]  =  np.sign(mtx2thr.sum(0))
        lbls2show[:, :, 2]  =  nucs_mip / nucs_mip.max()

        roi  =  pg.LinearRegionItem(orientation=False)
        roi.sigRegionChanged.connect(self.update_ts_sm_map)

        hh_ints  =  np.histogram(spts_ints_ctrs[2, :], bins=100)

        frame1_plot  =  pg.PlotWidget(self)
        frame1_plot.plot(hh_ints[1], hh_ints[0], stepMode=True, fillLevel=0, brush=(0, 0, 255, 150))
        frame1_plot.addItem(roi)

        frame2_img  =  pg.ImageView(self)
        frame2_img.ui.menuBtn.hide()
        frame2_img.ui.roiBtn.hide()
        frame2_img.setFixedSize(550, 550)
        frame2_img.setImage(lbls2show)

        ts_number_label  =  QtWidgets.QLabel()
        ts_number_label.setText("TS Numb       ")

        select_update_btn  =  QtWidgets.QPushButton("Select", self)
        select_update_btn.clicked.connect(self.select_update)
        select_update_btn.setToolTip('Split nascent from mature depending on the threshold')
        select_update_btn.setFixedSize(110, 25)

        zoomon_hist_btn  =  QtWidgets.QPushButton("Zoom +", self)
        zoomon_hist_btn.clicked.connect(self.zoomon_hist)

        zoomoff_hist_btn  =  QtWidgets.QPushButton("Zoom -", self)
        zoomoff_hist_btn.clicked.connect(self.zoomoff_hist)

        frames_box  =  QtWidgets.QHBoxLayout()
        frames_box.addWidget(frame1_plot)
        frames_box.addWidget(frame2_img)

        btns_box  =  QtWidgets.QHBoxLayout()
        btns_box.addWidget(zoomon_hist_btn)
        btns_box.addWidget(zoomoff_hist_btn)
        btns_box.addStretch()
        btns_box.addWidget(ts_number_label)
        btns_box.addWidget(select_update_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(frames_box)
        layout.addLayout(btns_box)

        self.roi              =  roi
        self.spts_ints_ctrs   =  spts_ints_ctrs
        self.spts_lbls        =  spts_lbls
        self.nucs_mip         =  nucs_mip
        self.frame2_img       =  frame2_img
        self.mtx2thr          =  mtx2thr
        self.ts_number_label  =  ts_number_label
        self.frame1_plot      =  frame1_plot
        self.y_hist_size      =  hh_ints[0].max()

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("MatureNascent")


    def update_ts_sm_map(self):
        """Updates the blue and red map to distinguish betweeb TS and sm"""
        lbls2show           =  np.zeros(self.spts_lbls.shape[1:] + (3,))
        lbls2show[:, :, 2]  =  self.nucs_mip / self.nucs_mip.max()
        lbls2show[:, :, 0]  =  np.sign((self.mtx2thr > self.roi.getRegion()[1]).sum(0))
#        lbls2show[:, :, 2]  =  np.sign(((self.mtx2thr < self.roi.getRegion()[1]) * (self.mtx2thr != 0)).sum(0))
        self.frame2_img.setImage(lbls2show)
        self.ts_number_label.setText("TS Numb = " + str(label(lbls2show[:, :, 0]).max()))


    def zoomon_hist(self):
        """Zoom in the histogram"""
        self.y_hist_size  /=  2
        self.frame1_plot.setYRange(0, self.y_hist_size)


    def zoomoff_hist(self):
        """Zoom out the histogram"""
        self.y_hist_size  *=  2
        self.frame1_plot.setYRange(0, self.y_hist_size)


    @QtCore.pyqtSlot()
    def select_update(self):
        """Select threshold value and send to main GUI"""
        thresh  =  self.roi.getRegion()[1]
        self.procStart.emit(thresh)


    # def closeEvent(self, event):
    #     """Confirmation dialog to close the popup tool"""
    #     quit_msg  =  "Are you sure you want to exit the program?"
    #     reply     =  QtWidgets.QMessageBox.question(self, 'Message', quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

    #     if reply == QtWidgets.QMessageBox.Yes:
    #         event.accept()
    #     else:
    #         event.ignore()



class SetColorChannel(QtWidgets.QDialog):
    """Set the color channels of the raw data to put in the gui"""
    def __init__(self, parent=None):
        super(SetColorChannel, self).__init__(parent)

        nuclei_channel_lbl  =  QtWidgets.QLabel("Nuclei Channel", self)
        nuclei_channel_lbl.setFixedSize(120, 22)

        spots_a_channel_lbl  =  QtWidgets.QLabel("Spots A Channel", self)
        spots_a_channel_lbl.setFixedSize(120, 22)

        spots_b_channel_lbl  =  QtWidgets.QLabel("Spots B Channel", self)
        spots_b_channel_lbl.setFixedSize(120, 22)

        nuclei_channel_combo  =  QtWidgets.QComboBox(self)
        nuclei_channel_combo.addItem("1")
        nuclei_channel_combo.addItem("2")
        nuclei_channel_combo.addItem("3")
        nuclei_channel_combo.addItem("4")
        nuclei_channel_combo.addItem("None")
        nuclei_channel_combo.activated[str].connect(self.nuclei_channel_switch)
        nuclei_channel_combo.setCurrentIndex(2)
        nuclei_channel_combo.setFixedSize(65, 25)

        spots_a_channel_combo  =  QtWidgets.QComboBox(self)
        spots_a_channel_combo.addItem("1")
        spots_a_channel_combo.addItem("2")
        spots_a_channel_combo.addItem("3")
        spots_a_channel_combo.addItem("4")
        spots_a_channel_combo.addItem("None")
        spots_a_channel_combo.activated[str].connect(self.spots_a_channel_switch)
        spots_a_channel_combo.setCurrentIndex(0)
        spots_a_channel_combo.setFixedSize(65, 25)

        spots_b_channel_combo  =  QtWidgets.QComboBox(self)
        spots_b_channel_combo.addItem("1")
        spots_b_channel_combo.addItem("2")
        spots_b_channel_combo.addItem("3")
        spots_b_channel_combo.addItem("4")
        spots_b_channel_combo.addItem("None")
        spots_b_channel_combo.activated[str].connect(self.spots_b_channel_switch)
        spots_b_channel_combo.setCurrentIndex(1)
        spots_b_channel_combo.setFixedSize(65, 25)

        enter_values_btn  =  QtWidgets.QPushButton("OK", self)
        enter_values_btn.setToolTip("Set Channels Number")
        enter_values_btn.setFixedSize(60, 25)
        enter_values_btn.clicked.connect(self.enter_values)

        nuclei_box  =  QtWidgets.QHBoxLayout()
        nuclei_box.addWidget(nuclei_channel_lbl)
        nuclei_box.addWidget(nuclei_channel_combo)

        spots_a_box  =  QtWidgets.QHBoxLayout()
        spots_a_box.addWidget(spots_a_channel_lbl)
        spots_a_box.addWidget(spots_a_channel_combo)

        spots_b_box  =  QtWidgets.QHBoxLayout()
        spots_b_box.addWidget(spots_b_channel_lbl)
        spots_b_box.addWidget(spots_b_channel_combo)

        enter_box  =  QtWidgets.QHBoxLayout()
        enter_box.addStretch()
        enter_box.addWidget(enter_values_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addLayout(spots_a_box)
        layout.addLayout(spots_b_box)
        layout.addLayout(nuclei_box)
        layout.addLayout(enter_box)

        self.nuclei_channel   =  3
        self.spots_a_channel  =  1
        self.spots_b_channel  =  2

        self.setWindowModality(Qt.ApplicationModal)
        self.setLayout(layout)
        self.setGeometry(300, 300, 250, 150)
        self.setWindowTitle("Set Channels")


    def nuclei_channel_switch(self, text):
        """Set nuclei channel"""
        self.nuclei_channel  =  text


    def spots_b_channel_switch(self, text):
        """Set spots b channel"""
        self.spots_b_channel  =  text


    def spots_a_channel_switch(self, text):
        """Set spots a channel"""
        self.spots_a_channel  =  text


    def enter_values(self):
        """Organizing channels info for the output"""
        self.chs_spts_nucs  =  []
        if self.spots_a_channel != "None":
            self.chs_spts_nucs.append(int(self.spots_a_channel))
        else:
            self.chs_spts_nucs.append(0)

        if self.spots_b_channel != "None":
            self.chs_spts_nucs.append(int(self.spots_b_channel))
        else:
            self.chs_spts_nucs.append(0)

        if self.nuclei_channel != "None":
            self.chs_spts_nucs.append(int(self.nuclei_channel))
        else:
            self.chs_spts_nucs.append(0)
#        self.chs_spts_nucs.append(int(self.nuclei_channel))

        self.chs_spts_nucs  =  np.asarray(self.chs_spts_nucs) - 1

        self.close()


    def params(self):
        """Function to send results"""
        return self.chs_spts_nucs


    @staticmethod
    def getChannels(parent=None):
        """Send results"""
        dialog  =  SetColorChannel(parent)
        result  =  dialog.exec_()
        flag    =  dialog.params()
        return flag



class CroppingTool(QtWidgets.QWidget):
    """Popup tool to crop raw data"""
    procStart  =  QtCore.pyqtSignal(int)

    def __init__(self, stack):
        QtWidgets.QWidget.__init__(self)

        framepp1  =  pg.ImageView(self)
        framepp1.ui.roiBtn.hide()
        framepp1.ui.menuBtn.hide()
        framepp1.setImage(stack)

        roi  =  pg.RectROI([80, 80], [80, 80], pen='r')
        framepp1.addItem(roi)

        send_crop_btn  =  QtWidgets.QPushButton("Crop", self)
        send_crop_btn.setFixedSize(120, 25)
        send_crop_btn.clicked.connect(self.crop_to_mainwindows)

        keys  =  QtWidgets.QHBoxLayout()
        keys.addStretch()
        keys.addWidget(send_crop_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(framepp1)
        layout.addLayout(keys)

        self.framepp1  =  framepp1
        self.roi       =  roi

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("Modifier Tool")


    @QtCore.pyqtSlot()
    def crop_to_mainwindows(self):
        val  =  self.framepp1.currentIndex
        self.procStart.emit(val)



class PopUpPlot(QtWidgets.QWidget):
    """Popup the plot in the GUI for better visualization"""
    procStart  =  QtCore.pyqtSignal(int)

    def __init__(self, thr_vals, spts_num):
        QtWidgets.QWidget.__init__(self)

        frame3_spts  =  pg.PlotWidget(self)
        frame3_spts.plot(thr_vals, spts_num, symbol='x', pen='r')

        sld_thr  =  QtWidgets.QScrollBar(QtCore.Qt.Horizontal, self)
        sld_thr.valueChanged.connect(self.sld_thr_update)
        sld_thr.setMaximum(thr_vals.size - 1)

        sld_val_lbl  =  QtWidgets.QLabel('Slider val 0', self)

        roi_sld  =  pg.LineSegmentROI([[thr_vals[0], 0], [thr_vals[0], spts_num.max()]], pen='b')
        frame3_spts.addItem(roi_sld)

        select_thr_btn  =  QtWidgets.QPushButton("Send", self)
        select_thr_btn.setToolTip("Select value and send to the main Window")
        select_thr_btn.setFixedSize(120, 25)
        select_thr_btn.clicked.connect(self.send_thr_val)

        keys_box  =  QtWidgets.QHBoxLayout()
        keys_box.addStretch()
        keys_box.addWidget(select_thr_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(frame3_spts)
        layout.addWidget(sld_thr)
        layout.addWidget(sld_val_lbl)
        layout.addLayout(keys_box)

        self.roi_sld      =  roi_sld
        self.sld_val_lbl  =  sld_val_lbl
        self.thr_vals     =  thr_vals
        self.sld_thr      =  sld_thr

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("PopUp Thresholding")


    def sld_thr_update(self):
        self.roi_sld.setPos([self.thr_vals[self.sld_thr.value()] - self.thr_vals[0], 0], update=True)
        self.sld_val_lbl.setText("Slider val " + str(self.sld_thr.value()))


    @QtCore.pyqtSlot()
    def send_thr_val(self):
        val  =  self.sld_thr.value()
        self.procStart.emit(val)



class ColocalizeTool(QtWidgets.QWidget):
    """Popup tool for colocalization study"""
    def __init__(self, post_proc_folder, raw_data_fname, soft_version, vol_thr, clust_dist):
        QtWidgets.QWidget.__init__(self)

        folder_name_lbl  =  QtWidgets.QLabel(post_proc_folder, self)
        folder_name_lbl.setToolTip('Name of the file you are working on')

        coloc_vals      =  ColocalizeSpots.ColocalizeSpots(post_proc_folder, raw_data_fname, soft_version, vol_thr, clust_dist)
        self.spts_clcz  =  coloc_vals.coloc_mtx
        
        perc_mrna_lbl  =  QtWidgets.QLabel("% mRNA in Transl", self)
        # perc_mrna_lbl.setToolTip('Name of the file you are working on')
        perc_mrna_lbl.setFixedSize(120, 25)

        perc_mrna_value_lbl  =  QtWidgets.QLabel(str(coloc_vals.mrna_in_transl), self)
        perc_mrna_value_lbl.setFixedSize(120, 25)
        # perc_mrna_value_lbl.setToolTip('Name of the file you are working on')

        frame1  =  pg.ImageView(self)
        frame1.ui.roiBtn.hide()
        frame1.ui.menuBtn.hide()
        frame1.setImage(self.spts_clcz)

        cmmd_box  =  QtWidgets.QVBoxLayout()
        cmmd_box.addWidget(perc_mrna_lbl)
        cmmd_box.addWidget(perc_mrna_value_lbl)
        cmmd_box.addStretch()

        frame_cmmd_box  =  QtWidgets.QHBoxLayout()
        frame_cmmd_box.addWidget(frame1)
        frame_cmmd_box.addLayout(cmmd_box)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(folder_name_lbl)
        layout.addLayout(frame_cmmd_box)

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("Colocalize Tool")



class ChopStack(QtWidgets.QWidget):
    """Popup tool to remove frames from the raw data stack"""
    procStart  =  QtCore.pyqtSignal(int)
    
    def __init__(self, raw_data):
        QtWidgets.QWidget.__init__(self)

        frame_spts_raw_a  =  pg.ImageView(self, name='FrameSpts1')
        frame_spts_raw_a.ui.roiBtn.hide()
        frame_spts_raw_a.ui.menuBtn.hide()
        frame_spts_raw_a.timeLine.sigPositionChanged.connect(self.update_frame_from_spts_a)
        frame_spts_raw_a.setImage(raw_data.spts_a)

        frame_spts_raw_b  =  pg.ImageView(self)
        frame_spts_raw_b.ui.roiBtn.hide()
        frame_spts_raw_b.ui.menuBtn.hide()
        frame_spts_raw_b.view.setXLink('FrameSpts1')
        frame_spts_raw_b.view.setYLink('FrameSpts1')
        frame_spts_raw_b.timeLine.sigPositionChanged.connect(self.update_frame_from_spts_b)
        frame_spts_raw_b.setImage(raw_data.spts_b)

        frame_nucs_raw  =  pg.ImageView(self)
        frame_nucs_raw.ui.roiBtn.hide()
        frame_nucs_raw.ui.menuBtn.hide()
        frame_nucs_raw.view.setXLink('FrameSpts1')
        frame_nucs_raw.view.setYLink('FrameSpts1')
        frame_nucs_raw.timeLine.sigPositionChanged.connect(self.update_frame_from_nucs)
        frame_nucs_raw.setImage(raw_data.nucs)

        tabs_tot    =  QtWidgets.QTabWidget()
        tab_spts_a  =  QtWidgets.QWidget()
        tab_spts_b  =  QtWidgets.QWidget()
        tab_nucs    =  QtWidgets.QWidget()

        frame_box_spts_a  =  QtWidgets.QHBoxLayout()
        frame_box_spts_a.addWidget(frame_spts_raw_a)

        frame_box_spts_b  =  QtWidgets.QHBoxLayout()
        frame_box_spts_b.addWidget(frame_spts_raw_b)

        frame_box_nucs  =  QtWidgets.QHBoxLayout()
        frame_box_nucs.addWidget(frame_nucs_raw)

        tab_spts_a.setLayout(frame_box_spts_a)
        tab_spts_b.setLayout(frame_box_spts_b)
        tab_nucs.setLayout(frame_box_nucs)

        tabs_tot.addTab(tab_spts_a, "Spots A")
        tabs_tot.addTab(tab_spts_b, "Spots B")
        tabs_tot.addTab(tab_nucs, "Nucs")

        current_frame_lbl  =  QtWidgets.QLabel("Current Frame 0", self)
        current_frame_lbl.setFixedSize(120, 25)

        first_frame_lbl  =  QtWidgets.QLabel("First Frame 0", self)
        first_frame_lbl.setFixedSize(120, 25)

        last_frame_lbl  =  QtWidgets.QLabel("Last Frame " + str(raw_data.nucs.shape[0]), self)
        last_frame_lbl.setFixedSize(120, 25)

        set_first_frame_btn  =  QtWidgets.QPushButton("Set First", self)
        set_first_frame_btn.clicked.connect(self.set_first_frame)
        set_first_frame_btn.setToolTip('Set current frame as first analysis frame')
        set_first_frame_btn.setFixedSize(130, 25)

        set_last_frame_btn  =  QtWidgets.QPushButton("Set Last", self)
        set_last_frame_btn.clicked.connect(self.set_last_frame)
        set_last_frame_btn.setToolTip('Set current frame as last analysis frame')
        set_last_frame_btn.setFixedSize(130, 25)

        # close_send_btn  =  QtWidgets.QPushButton("Close && Insert", self)
        # close_send_btn.clicked.connect(self.close_insert)
        # close_send_btn.setToolTip('Apply changes in main GUI')
        # close_send_btn.setFixedSize(130, 25)

        rmv_selected_btn  =  QtWidgets.QPushButton("Remove Non Selected", self)
        rmv_selected_btn.clicked.connect(self.close_insert_rmv)
        rmv_selected_btn.setToolTip('Remove selected frames and apply changes in main GUI')
        rmv_selected_btn.setFixedSize(140, 25)

        keep_selected_btn  =  QtWidgets.QPushButton("Zero On Selected", self)
        keep_selected_btn.clicked.connect(self.close_insert_keep)
        keep_selected_btn.setToolTip('Keep selected frames and apply changes in main GUI')
        keep_selected_btn.setFixedSize(140, 25)

        rmv_keep_btn_box  =  QtWidgets.QVBoxLayout()
        rmv_keep_btn_box.addWidget(rmv_selected_btn)
        rmv_keep_btn_box.addWidget(keep_selected_btn)

        first_lbl_btn_box  =  QtWidgets.QHBoxLayout()
        first_lbl_btn_box.addWidget(set_first_frame_btn)
        first_lbl_btn_box.addWidget(first_frame_lbl)

        last_lbl_btn_box  =  QtWidgets.QHBoxLayout()
        last_lbl_btn_box.addWidget(set_last_frame_btn)
        last_lbl_btn_box.addWidget(last_frame_lbl)

        first_last_box  =  QtWidgets.QVBoxLayout()
        first_last_box.addLayout(first_lbl_btn_box)
        first_last_box.addLayout(last_lbl_btn_box)

        lbls_buttons_box  =  QtWidgets.QHBoxLayout()
        lbls_buttons_box.addWidget(current_frame_lbl)
        lbls_buttons_box.addStretch()
        lbls_buttons_box.addLayout(first_last_box)
        lbls_buttons_box.addLayout(rmv_keep_btn_box)
        # lbls_buttons_box.addWidget(close_send_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(tabs_tot)
        layout.addLayout(lbls_buttons_box)

        self.frame_spts_raw_a   =  frame_spts_raw_a
        self.frame_spts_raw_b   =  frame_spts_raw_b
        self.frame_nucs_raw     =  frame_nucs_raw
        self.current_frame_lbl  =  current_frame_lbl
        self.first_frame_lbl    =  first_frame_lbl
        self.last_frame_lbl     =  last_frame_lbl
        self.first_last_frame   =  np.array([0, raw_data.nucs.shape[0]])

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("ChopStack Tool")


    def update_frame_from_spts_a(self):
        """Update frame index from spots raw a to all the other"""
        self.frame_spts_raw_b.setCurrentIndex(self.frame_spts_raw_a.currentIndex)
        self.frame_nucs_raw.setCurrentIndex(self.frame_spts_raw_a.currentIndex)
        self.current_frame_lbl.setText("Current Frame " + str(self.frame_spts_raw_a.currentIndex))


    def update_frame_from_spts_b(self):
        """Update frame index from raw spots b to raw spots a"""
        self.frame_spts_raw_a.setCurrentIndex(self.frame_spts_raw_b.currentIndex)


    def update_frame_from_nucs(self):
        """Update frame index from raw nucs to raw spots a"""
        self.frame_spts_raw_a.setCurrentIndex(self.frame_nucs_raw.currentIndex)


    def set_first_frame(self):
        """Set first frame choise"""
        self.first_last_frame[0]  =  self.frame_spts_raw_a.currentIndex
        self.first_frame_lbl.setText("First Frame " + str(self.first_last_frame[0]))


    def set_last_frame(self):
        """Set last frame choise"""
        self.first_last_frame[1]  =  self.frame_spts_raw_a.currentIndex
        self.last_frame_lbl.setText("Last Frame " + str(self.first_last_frame[1]))


    # @QtCore.pyqtSlot()
    # def close_insert(self):
    #     """Send message to main GUI"""
    #     self.procStart.emit()


    @QtCore.pyqtSlot()
    def close_insert_rmv(self):
        """Send message to main GUI"""
        self.procStart.emit(0)


    @QtCore.pyqtSlot()
    def close_insert_keep(self):
        """Send message to main GUI"""
        self.procStart.emit(1)



class NucsFittingVisual(QtWidgets.QWidget):
    """Tool to visualize and verify the nuclei ellipsoidal fitting"""
    def __init__(self, raw_nucs_norm, nucs_ellips_sign):
        QtWidgets.QWidget.__init__(self)
        
        img2show              =  np.zeros(nucs_ellips_sign.shape + (3, ))
        img2show[:, :, :, 2]  =  nucs_ellips_sign
        img2show[:, :, :, 0]  =  raw_nucs_norm

        frame_nucs_ellips  =  pg.ImageView(self)
        frame_nucs_ellips.ui.menuBtn.hide()
        frame_nucs_ellips.ui.roiBtn.hide()
        frame_nucs_ellips.setImage(img2show)

        raw_max_sld  =  QtWidgets.QScrollBar(QtCore.Qt.Vertical, self)
        raw_max_sld.setMaximum(40)
        raw_max_sld.setMinimum(1)
        raw_max_sld.valueChanged.connect(self.raw_max_sld_update)

        layout  =  QtWidgets.QHBoxLayout()
        layout.addWidget(frame_nucs_ellips)
        layout.addWidget(raw_max_sld)

        self.raw_max_sld        =  raw_max_sld
        self.img2show           =  img2show   
        self.frame_nucs_ellips  =  frame_nucs_ellips  
        self.nucs_ellips_sign   =  nucs_ellips_sign

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("Visualize Fitting Tool")


    def raw_max_sld_update(self):
        self.img2show[:, :, :, 2]  =  (1 / self.raw_max_sld.value()) * self.nucs_ellips_sign
        self.frame_nucs_ellips.updateImage()
#        self.frame_nucs_ellips.autoLevels()



# def except_hook(cls, exception, traceback):
#     sys.__excepthook__(cls, exception, traceback)



#if __name__ == "__main__":
#
#    app         =  QtWidgets.QApplication(sys.argv)
#    splash_pix  =  QtWidgets.QPixmap('Icons/DrosophilaIcon.png')                     #
#    splash      =  QtWidgets.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)      #
#    splash.setMask(splash_pix.mask())                                            #
#    splash.show()                                                                #
#    app.processEvents()                                                          #
#    window      =  MainWindow()
#
#    splash.finish(window)                                                        #
#    window.show()
#    sys.exit(app.exec_())
#


#def main():
#    app  =  QtWidgets.QApplication(sys.argv)
#    ex   =  MainWindow()
#    sys.exit(app.exec_())
#
#
#
#def except_hook(cls, exception, traceback):
#    sys.__excepthook__(cls, exception, traceback)
#
#
#if __name__ == '__main__':
#
#    sys.excepthook = except_hook
#    main()



def main():
    app         =  QtWidgets.QApplication(sys.argv)
    splash_pix  =  QtGui.QPixmap('Icons/DrosophilaIcon.png')                     #
    splash      =  QtWidgets.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)      #
    splash.setMask(splash_pix.mask())                                            #
    splash.show()                                                                #
    app.processEvents()                                                          #
    ex   =  MainWindow()
    splash.finish(ex)                                                            #
    sys.exit(app.exec_())



def except_hook(cls, exception, traceback):
    sys.__excepthook__(cls, exception, traceback)


if __name__ == '__main__':

    sys.excepthook = except_hook
    main()

