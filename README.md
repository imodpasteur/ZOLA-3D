# ZOLA-3D
ZOLA-3D (Zernike Optimized Localization Approach) is a full software package to reconstruct 3D single molecule localization images for a wide range of point spread functions.
Here we provide a CUDA-enabled 3D localization plugin for ImageJ/Fiji software. 

Please, open the  `User manual - HOW TO USE ZOLA.pdf` file for an extended explanation of how to use ZOLA.

## Requirements


To enable fast processing on GPU, Nvidia CUDA must be installed and included in the path. CUDA 8 is  compatible with our tools.
Check your version by executing `nvcc --version` in the terminal.

Follow instructions on https://developer.nvidia.com/cuda-toolkit to install CUDA

The plugin was tested on Ubuntu, Windows and Mac with Cuda 8.



## Installation

The easiest way to install ZOLA-3D is to use the Fiji update system.

1. Launch Fiji
2. Select `Help` -> `Update`
3. Select `Manage update sites`
4. Select `IMOD-ZOLA` from the list
5. Click Close
6. Click Apply changes
7. Restart Fiji
7. ZOLA should now appear in the Fiji->Plugins dropdown menu
![Fiji update site](https://github.com/imodpasteur/ZOLA-3D/blob/master/images/fiji%20update.png)

If you prefer to use ImageJ instead, you can download `ZOLA_-0.1.X-SNAPSHOT.jar` plugin and corresponding libraries from here.

If you want to benefit from GPU acceleration, `lib/*` and `jar/*` folders containing jcuda bindings should be copied into Fiji/ImageJ folder. 

## Usage

The full manual and test data can be found [HERE](https://github.com/imodpasteur/ZOLA-3D/tree/master/TEST_DATA)

1. For optimal performance, your images should be in photon counts, as usually the case with recent EMCCD cameras. images are not in photon counts, you can convert them as explained in [full manual](https://github.com/imodpasteur/ZOLA-3D/blob/master/User%20manual%20-%20HOW%20TO%20USE%20ZOLA.pdf) .

2. Calibrate your PSF from a z-stack of fluorescent beads using Zernike polynomials (this is ideal for PSFs shaped using a cylindrical lens or deformable mirror). We recommend using z-steps of 50 or 100 nm.
Open your stack, select central position of one or few beads with the Fiji `Point selector tool` and run `ZOLA -> Calibration: PSF modeling`
Choose 28-45 Zernike modes.

3. Now that the PSF is calibrated, you can use ZOLA for 3D single molecule localization based super-resolution microscopy. 
Open your stack with single molecule blinking images, select region of interest and run `ZOLA -> Localization`

4. The single molecule Localization table is saved in csv, Thunderstorm format. You can correct the drift in 3D and render 3D color-stack using ZOLA->Drift correction->3D Drift correction  or any other software. 

### Example of 3D reconstruction of immunolabeled mitochondrial protein TOM22:

| ![input](https://github.com/imodpasteur/ZOLA-3D/blob/master/images/frames20130%2B50.gif)        | ![output](https://github.com/imodpasteur/ZOLA-3D/blob/master/images/anim-slow.gif) | 
| ------------- |:-------------:|
| *raw image sequence with blinking tetrapod PSFs*| *3D scatter-plot reconstruction (Visp rendering)* |

For detailed instructions please follow step by step guide on how to run [PSF calibration](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/calibration_howto.md) and [localization](https://github.com/imodpasteur/ZOLA-3D/blob/master/TEST_DATA/localization_howto.md)



## Distinctive features:

* User-friendly interface implemented in Fiji.
* Ability to reconstruct 3D single molecule images with variable axial range (e.g. 1-5 micrometers).
* Ability to model a wide range of spatially extended PSFs (including astigmatism, saddle point and tetrapod) and to detect and localize single molecules using cross-correlation and maximum likelihood algorithms.
* Realistic and robust calibration of the PSF using only one or a few beads thanks to a Zernike-based PSF model. Depth-dependent aberrations due to refractive index mismatch  are automatically handled using beads on the coverslip only. 
* All calculations can run on GPU. Main features also run on CPU, but processing on GPUs is much faster.


