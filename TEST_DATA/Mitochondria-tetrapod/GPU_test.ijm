CAL_DATA="calibration"+File.separator+"cal-sp4-50nm_crop.tif";
CAL_RES = "calibration"+File.separator+"ZOLA_calibration.csv";
DATA="data"+File.separator+"Mito-sp4-SR-data-crop-small.tif";
LOC="data"+File.separator+"ZOLA_localization_table";

open(CAL_DATA);
makePoint(127, 92);
run(" Calibration: PSF modeling", "run_on_gpu camera_adu=1.0000 camera_gain=1.0000 camera_offset=0.0000 pixel_size=106 z_step=50 bead_moving=[far -> close to objective] numerical_aperture=1.200 immersion_refractive=1.330 wavelength=0.640 patch_size=36 zernike_coefficient=45 iteration=10 result_calibration_file="+CAL_RES);

//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
//selectWindow("cal-sp4-50nm_1_MMStack_Pos0.ome.tif");
open(DATA);
run(" Localization", "run_on_gpu camera_emccd_adu=1.0000 camera_emccd_gain=1.0000 camera_emccd_offset=0.0000 calibration_file="+CAL+" mounting_medium_refractive=1.330 distance_focus_to_coverslip=0.5 patch_size=16 expected_axial_range=1.5 min_number_of_photons=1000 localization_table="+LOC);
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
run("3D Drift correction", "run_on_gpu cross-correlation_pixel_size=50 number=5 maximum_drift=6 localization_table_attached=[]");
run("3D Drift correction");
selectWindow("Reslice of 3D");
selectWindow("c1-sp4_1_MMStack_Pos0.ome.tif");
selectWindow("Mito-sp4-SR-data-crop.tif");
selectWindow("Live rendering 20nm/pixel");
*/
