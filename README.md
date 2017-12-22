# ZOLA-3D
Zernike Optimized Localization Algorithm for 3D single molecule localizations

Here we provide a CUDA-enabled 3D localization algorithm for ImageJ/Fiji distribution. 

## Requirements

Nvidia CUDA 8 must be installed and included in the path. 

Follow instructions on https://developer.nvidia.com/cuda-80-ga2-download-archive. 

Please make sure `nvcc --version` is operational from the terminal.


The plugin is tested on Ubuntu and Windows. MacOS libraries are provided but CUDA ubstallation can be tricky.

## Installation

The easy way is to use Fiji update system. 
1. Launch Fiji
2. Select `Help` -> `Update`
3. Select `Manage update sites`
4. Select `IMOD-ZOLA` from the list
5. Close and Apply changes
6. Restart Fiji
7. ZOLA will appear in Plugins dropdown menu
![Fiji update site](https://github.com/imodpasteur/ZOLA-3D/blob/master/images/fiji%20update.png)

## Usage

The full manual and test data can be found [HERE](https://www.dropbox.com/sh/5h4kz7ruuv3iw0b/AAD4JSNIT-L17mVr1EqMi2WRa?dl=0)

1. For optimal performance your images should obey poisson noise statistics, equivalent of photon counting mode of EMCCD camera. If your camera does not support this mode of acquisition, photon conversion step is needed:
(a) Acquire 50 frames with static intensity distribution (beads or gradient illumination). 
(b) Select ~100x100 px region with highest intensity range (proximity of the beads, avoiding saturation). 
(c) Run `Zola -> Additional tools -> Photon conversion`. This will autimatically find offset and gain to get poison noise statistics out of raw counts. The optimal parameters of count conversion are recorded automatically and propagate throughout the ZOLA interface.

2. Calibrate your PSF from a stack of fluorescent beads using Zernike polynomials (ideal for cylyndrical lens or deformable mirror). 50-100 nm z step is optimal.
Open your stack, select central position of one or few beads with `Point selector` and run `ZOLA -> Calibration: PSF modeling`

3. Next, this calibrations is used for localization of single molecules in 3D. 
Open your stack with single molecule blinking images, select region of interest and run `ZOLA -> Localization`

![input](https://github.com/imodpasteur/ZOLA-3D/blob/master/images/frames20130%2B50.gif) ![output](https://github.com/imodpasteur/ZOLA-3D/blob/master/images/anim-slow.gif)


## Distinctive features:

* User-fliendly interface implemented in Fiji.
* We incorporate 3D predetection (instead of just finding a peak) into the algorithm, so the PSF can be spatially extended, such as tetrapod PSF.
* We use integrated model of PSF encoded into pupil function. This makes PSF calibration extremely robust to the noise. As a result, only one or two beads are needed to calibrate PSF.
* Index mismatch is included. Just select the correct immersion and mounting medium parameters.
* All calculations are done on GPU.


