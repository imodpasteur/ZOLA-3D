# Running ZOLA PSF calibration on test dataset
We provide instructions to run calibration on the test data in the package.

## Test data

Download and unzip [TEST-DATA.zip](https://github.com/imodpasteur/ZOLA-3D/releases/download/v0.1.9/TEST-DATA.zip) to get Mitochondria test data.

Mitochondria with corresponding calibration stack is imaged using tetrapod PSF.

## Calibration

1. Open Mitochondria-tetrapod/calibration/cal-sp4-50nm_crop.tif
2. Select a single bead

![calibration bead selection](https://github.com/imodpasteur/ZOLA-3D/blob/master/img/ZOLA_cal_bead_screenshot.png)

3. Setup the camera first. Select `Plugins` -> `ZOLA` -> `Camera setup` -> `EMCCD`. Data was acquired in photon counts, so we don't need to enter camera sensitivity values.

![EMCCD setup dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/img/ZOLA_camera_setup_EMCCD.png)

4. Select `Plugins` -> `ZOLA` -> `Calibration: PSF modeling` and enter optical and computational parameters.

Pixel size and z steps should be adjusted, otherwise default values of fitting parameters should be fine.

Choosing the crop size try to avoid multiple beads per crop.

If you want to open save file dialog, double click on the path below. Click OK to start calibration.

![calibration dialog](https://github.com/imodpasteur/ZOLA-3D/blob/master/img/ZOLA_cal_dialog_screenshot.png)

5. At the end of calibration process you will see the combined stack with data/model comparison.

![calibration output](https://github.com/imodpasteur/ZOLA-3D/blob/master/img/ZOLA_cal_bead_output.gif)
