import numpy as np
from PyQt5 import QtWidgets
from skimage.measure import regionprops_table
import tifffile

xls_fname        =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data file of the analyze", filter='*.xls')[0])
analysis_folder  =  xls_fname[: len(xls_fname) - xls_fname[::-1].find('/')]

spts_red    =  np.load(analysis_folder + 'spts_lbls_b.npy')
spts_green  =  np.load(analysis_folder + 'spts_lbls_a.npy')

# spts_green  *=  np.sign(spts_red)


rgp         =  regionprops_table(spts_red, spts_green, properties=["label", "coords"])
tags_green  =  []
for index, k in enumerate(rgp["coords"]):
    print(str(index) + "/" + str(rgp["coords"].shape[0]))
    bff_tags  =  []
    for p in range(k.shape[0]):
        bff_tags.append(spts_green[k[p, 0], k[p, 1], k[p, 2]])
    bff_tags  =  np.asarray(bff_tags)
    bff_tags  =  bff_tags[bff_tags != 0]
    if bff_tags.shape[0] != 0:
        if np.diff(bff_tags).sum() > 0:
            # break
            tags_green = tags_green +  list(np.unique(bff_tags))

ant              =  np.zeros(spts_green.shape + (3,))
ant[:, :, :, 0]  =  np.sign(spts_red)
idx  =  0
for tag in tags_green:
    ant[:, :, :, 1] += spts_green == tag
    print(idx)
    idx += 1


tifffile.imwrite(analysis_folder + xls_fname[len(xls_fname) - xls_fname[::-1].find('/'):xls_fname.find('Green')] + 'RedOnSeveralGreen.tif', np.rot90(ant[:, :, ::-1, :], axes=(1, 2)).astype("uint8"))
