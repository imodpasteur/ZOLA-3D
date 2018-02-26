# Running ZOLA PSF localization on test dataset
We provide instructions on how to run 3D localizations on the test data from the ZOLA package.

## Test data

Single-molecule localilization data can be found in the TEST_DATA directory of the project.

Mitochondrial protein TOM22 immunolabelled with Alexa 647 with corresponding bead calibration stack is imaged using tetrapod PSF.

## Calibration

1. Open ZOLA-3D/TEST_DATA/Mitochondria-tetrapod/data/Mito-sp4-SR-data-crop-small.tif
2. Select `plugins` -> `ZOLA` -> `Localization` and enter optical and computational parameters. Click OK to start localization.

![localization dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_loc_mito_screenshot.png)

3. You will soon see automatically updated color-coded histogram.

![automatic histogram](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_loc_mito_2D_hist20nm.png)

4. Once localization is done (elapsed time ~4 minutes on our Tesla GPU), you can filter the localization table. 
Select `ZOLA -> Filtering`.
You will see three new windows: disclaimer, scatter plot which can be used to select ROI, and brightness and contrast dialog, which represents now histogram on z localizations.
We select a broad peak coresponding to the useful data and avoinig artefacts in the lower part of the axial range.
Click OK, enter maximum value of Chi2 = 3 and click OK once again. 
Now our table is filtered.

![filtering dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_loc_filter.png)

5. In order to an axial projection we need to generate 3D histogram `ZOLA -> Visualization -> 2D/3D historgam`

![histogram dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_loc_hist_dialog.png)

6. Next, run reslice (keyboard shortcut '/') form the left to see YZ stack. Mitochondrial organelles' hollow structure is readily apparent.

![mito ZY projection](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_loc_mito_reslice.gif)

7. Don't forget to save updated localization table `ZOLA -> Import/Export -> Export table`


