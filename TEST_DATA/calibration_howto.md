# Running ZOLA PSF calibration on test dataset
We provide instructions to run calibration on the test data in the package.

## Test data

Located in the TEST_DATA directory of the project.

Mitochondria with corresponding calibration stack is imaged using tetrapod PSF.

## Calibration

1. Open ZOLA-3D/TEST_DATA/Mitochondria-tetrapod/calibration/cal-sp4-50nm_crop.tif
2. Select a single bead

![calibration bead selection](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_cal_bead_screenshot.png)

3. Setup the camera first. Select `plugins` -> `ZOLA` -> `Camera setup` -> `EMCCD`. Data was acquired in photon counts, so we don't need to enter camera sensitivity values.

![EMCCD setup dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_camera_setup_EMCCD.png)

4. Select `plugins` -> `ZOLA` -> `Calibration: PSF modeling` and enter optical and computational parameters.

Pixel size and z steps should be adjusted, otherwise default values of fitting parameters should be fine.

If you want to open save file dialog, double click on the path below. Click OK to start calibration.

![calibration dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_cal_dialog_screenshot.png)

5. At the end of calibration process you will see the combined stack with data/model comparison.

![calibration output](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/img/ZOLA_cal_bead_output.gif)
