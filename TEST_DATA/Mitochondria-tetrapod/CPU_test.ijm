
//open("cal-sp4-50nm_crop.tif");
//makePoint(127, 92);
//run(" Calibration: PSF modeling", "run_on_cpu camera_adu=1.0000 camera_gain=1.0000 camera_offset=0.0000 pixel_size=106 z_step=50 bead_moving=[far -> close to objective] numerical_aperture=1.200 immersion_refractive=1.330 wavelength=0.640 patch_size=36 zernike_coefficient=45 iteration=10 result_calibration_file=ZOLA_calibration_PSF.csv");

//run("Brightness/Contrast...");
//run("Enhance Contrast", "saturated=0.35");
//selectWindow("cal-sp4-50nm_1_MMStack_Pos0.ome.tif");
open("Mito-sp4-SR-data-crop-small.tif");
run(" Localization", "run_on_gpu camera_emccd_adu=1.0000 camera_emccd_gain=1.0000 camera_emccd_offset=0.0000 calibration_file=/home/andrey/atlas/Andrey/data/20170824-Wei-Mito-Tom22-best/s2-10beta-10Cot-20GL/ZOLA_calibration_PSF_crop.csv mounting_medium_refractive=1.330 distance_focus_to_coverslip=0.5 patch_size=16 expected_axial_range=4 min_number_of_photons=500 localization_table=ZOLA_calibration_PSF.csv");
run("2D/3D histogram");
selectWindow("3D histogram 20.0nm/px");
//run("Brightness/Contrast...");
/*
run("Enhance Contrast", "saturated=0.35");
run("Reslice [/]...", "output=20.000 start=Left flip");
run("Out [-]");
run("Out [-]");
run("Out [-]");
run("Out [-]");
selectWindow("c1-sp4_1_MMStack_Pos0.ome-1.tif");
run("Save");
selectWindow("Live rendering 20nm/pixel");
selectWindow("c1-sp4_1_MMStack_Pos0.ome-1.tif");
selectWindow("Live rendering 20nm/pixel");
run("3D Drift correction", "run_on_cpu cross-correlation_pixel_size=50 number=5 maximum_drift=6 localization_table_attached=[]");
run("3D Drift correction");
selectWindow("Reslice of 3D");
selectWindow("c1-sp4_1_MMStack_Pos0.ome.tif");
selectWindow("Mito-sp4-SR-data-crop.tif");
selectWindow("Live rendering 20nm/pixel");
*/
