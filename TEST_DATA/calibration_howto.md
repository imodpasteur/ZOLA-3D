# Running ZOLA PSF calibration on test dataset
We provide instructions to run calibration on the test data in the package.

## Test data

Located in the TEST_DATA directory of the project.

Mitochondria with corresponding calibration stack is imaged using tetrapod PSF.

## Calibration

1. Open ZOLA-3D/TEST_DATA/Mitochondria-tetrapod/calibration/cal-sp4-50nm_crop.tif
2. Select a single bead

![calibration bead selection](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_cal_bead_screenshot.png)

3. Select `plugins` -> `ZOLA` -> `Calibration: PSF modeling` and enter optical and computational parameters. Data acquired in photon counts, so we dont' need to enter camera sensitivity values. If you want to to open save file dialot, double click on the path below. Click OK to start calibration.

![calibration dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_cal_dialog_screenshot.png)

4. Check the resulting stack with data/model comparison.

![calibration output](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_cal_bead_output.gif)
