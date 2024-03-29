/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.configuration;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.WaitForUserDialog;
import java.awt.Font;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.io.InputStream;
import java.util.jar.JarOutputStream;
import java.util.jar.JarEntry;
import java.io.FileOutputStream;
import java.io.OutputStream;
/**
 *
 * @author benoit
 */
public class ConfigFile {
    
    
    
    
    static String [] option = {
        
        "Plugins>ZOLA>  Camera setup>, \"EMCCD calibration photon count\", org.pasteur.imagej.ZOLA(\"zola_photonCountEstimation\")\n"+
        
        "Plugins>ZOLA>  Camera setup>SCMOS calibration photon count, \"compute offset&variance maps (1/2)\", org.pasteur.imagej.ZOLA(\"zola_SCMOS_offset_variance_map\")\n"+
        "Plugins>ZOLA>  Camera setup>SCMOS calibration photon count, \"compute gain map (2/2)\", org.pasteur.imagej.ZOLA(\"zola_SCMOS_gain_map\")\n"+
        "Plugins>ZOLA>  Camera setup>, \"test photon count conversion\", org.pasteur.imagej.ZOLA(\"zola_convertImageToPhotonCount\")\n",
        
        "Plugins>ZOLA>additional tools, \"Best PSF model optimization\", org.pasteur.imagej.ZOLA(\"zola_phaseoptim\")\n",
        
        
                
        
        
        "Plugins>ZOLA>additional tools>Simulation, \"imageFromTable\", org.pasteur.imagej.ZOLA(\"zola_imageFromTable\")\n" +
        "Plugins>ZOLA>additional tools>Simulation>PSF generator, \"Initialisation\", org.pasteur.imagej.ZOLA(\"zola_getModels_GPUalloc\")\n" +
"Plugins>ZOLA>additional tools>Simulation>PSF generator, \"computation\", org.pasteur.imagej.ZOLA(\"zola_getModels_GPUcomput\")\n" +
"Plugins>ZOLA>additional tools>Simulation>PSF generator, \"freeMemory\", org.pasteur.imagej.ZOLA(\"zola_getModels_GPUfree\")\n",
        
        "Plugins>ZOLA>Dual camera, \"Simple Registration\", org.pasteur.imagej.ZOLA(\"zola_dual_Cam_registration\")\n" +
"Plugins>ZOLA>Dual camera, \"Fusion\", org.pasteur.imagej.ZOLA(\"zola_dual_Cam_localization\")\n" +
"Plugins>ZOLA>Dual camera, \"CRLB_\", org.pasteur.imagej.ZOLA(\"zola_dualCamCRLB\")\n",
        
        
        
        
        "Plugins>ZOLA>Dual color, \"Compute chromatic aberrations\", org.pasteur.imagej.ZOLA(\"zola_computeChromaticAberrations\")\n" +
"Plugins>ZOLA>Dual color, \"Compute SR 2 LR registration\", org.pasteur.imagej.ZOLA(\"zola_computeSR2LRregistration\")\n" +
"Plugins>ZOLA>Dual color, \"Color Registration\", org.pasteur.imagej.ZOLA(\"zola_colorRegistration\")\n",
    
        
        "Plugins>ZOLA>Drift correction, \"3D Fiducial-markers\", org.pasteur.imagej.ZOLA(\"zola_driftCorrectionFiducialMarker\")\n"+
        "Plugins>ZOLA>Drift correction, \"set Drift\", org.pasteur.imagej.ZOLA(\"zola_setDrift\")\n",
        
        "Plugins>ZOLA>additional tools>statistic, \"localization histogram\", org.pasteur.imagej.ZOLA(\"zola_stathisto\")\n"+
        "Plugins>ZOLA>additional tools>statistic, \"measure residence time\", org.pasteur.imagej.ZOLA(\"zola_statAttachTime1\")\n",
        
        "Plugins>ZOLA>additional tools>wobble correction, \"calibration\", org.pasteur.imagej.ZOLA(\"zola_wobblecalibration\")\n"+
        "Plugins>ZOLA>additional tools>wobble correction, \"correction\", org.pasteur.imagej.ZOLA(\"zola_wobblecorrection\")\n",
        
        "Plugins>ZOLA>Visualization, \"frame color coding\", org.pasteur.imagej.ZOLA(\"zola_timeRendering\")\n"+
        "Plugins>ZOLA>Visualization, \"Multiple_ROI\", org.pasteur.imagej.ZOLA(\"zola_multipleROIRendering\")\n"+
        "Plugins>ZOLA>Visualization, \"localMaxima\", org.pasteur.imagej.ZOLA(\"zola_localMaxima\")\n"+
            "Plugins>ZOLA>Visualization, \"split_4_FRC-FSC\", org.pasteur.imagej.ZOLA(\"zola_plot4FRC_FSC\")\n"+
        "Plugins>ZOLA>Visualization, \"scatter_plot\", org.pasteur.imagej.ZOLA(\"zola_scatterPlot\")\n",
        
    "Plugins>ZOLA>dev, \"Multi-emitter localization\", org.pasteur.imagej.ZOLA(\"zola_multi-emitter localization\")\n" +
"Plugins>ZOLA>dev, \"merge frames\", org.pasteur.imagej.ZOLA(\"zola_mergeframe\")\n" +
"Plugins>ZOLA>dev, \"test\", org.pasteur.imagej.ZOLA(\"zola_test\")\n" +
"Plugins>ZOLA>dev, \"test2\", org.pasteur.imagej.ZOLA(\"zola_test2\")\n" +
"Plugins>ZOLA>dev, \"test3\", org.pasteur.imagej.ZOLA(\"zola_test3\")\n" +
"Plugins>ZOLA>dev, \"test4\", org.pasteur.imagej.ZOLA(\"zola_test4\")\n" +
"Plugins>ZOLA>dev>Simulation, \"Simulation\", org.pasteur.imagej.ZOLA(\"zola_simulation\")\n" +
"Plugins>ZOLA>dev>Simulation, \"Simulation2beads\", org.pasteur.imagej.ZOLA(\"zola_simulation2beads\")\n" +
"Plugins>ZOLA>dev>Simulation, \"SimulationOverlappingBeads\", org.pasteur.imagej.ZOLA(\"zola_simulationOverlappingBeads\")\n" +
"Plugins>ZOLA>dev>Simulation, \"SimulationOneBead\", org.pasteur.imagej.ZOLA(\"zola_simulationOneBead\")\n" +
"Plugins>ZOLA>dev>Simulation, \"addPoissonNoise\", org.pasteur.imagej.ZOLA(\"zola_addPoissonNoise\")\n" +
"Plugins>ZOLA>dev, \"CRLBfromFile\", org.pasteur.imagej.ZOLA(\"zola_crlb_FromFile\")\n" +
"Plugins>ZOLA>dev, \"CRLBfromFileDualObj\", org.pasteur.imagej.ZOLA(\"zola_crlb_FromFileDualObj\")\n" +
"Plugins>ZOLA>dev, \"Photon conversion wrong\", org.pasteur.imagej.ZOLA(\"zola_photonConversion\")\n"+
"Plugins>ZOLA>dev>Resolution, \"points Measurment\", org.pasteur.imagej.ZOLA(\"zola_pointResolution\")\n" +
"Plugins>ZOLA>dev>Resolution, \"filament Measurment\", org.pasteur.imagej.ZOLA(\"zola_filamentResolution\")\n" +
"Plugins>ZOLA>dev>, \"RemoveLocalizations\", org.pasteur.imagej.ZOLA(\"zola_removeLocalizations\")\n" +
"Plugins>ZOLA>dev>, \"NPC detection\", org.pasteur.imagej.ZOLA(\"zola_NPCdetection\")\n" +
"Plugins>ZOLA>dev>, \"CRLB_3D\", org.pasteur.imagej.ZOLA(\"zola_crlb_3D\")\n" +
"Plugins>ZOLA>dev>, \"PSF_Generator\", org.pasteur.imagej.ZOLA(\"zola_PSF_Generator\")\n" +
"Plugins>ZOLA>dev>, \"Dec\", org.pasteur.imagej.ZOLA(\"zola_deconvolution\") "};
    
    
    static String [] message = {
        "Photon_count_estimation :",
        "Phase_optimization :",
        "PSF_simulation :",
        "Dual_cam_localization (beta) :",
        "Dual_color_registration (beta) :",
        "Drift correction using fiducial markers (alpha) :",
        "Statistics: histograms/mean particle attach time :",
        "Wobble correction :",
        "Visualization: additional color rendering :",
        "developpement (prototype: not recommended) :",
        };
                
    
    static String fileName="ZOLA_AddOn.jar";
    
    public static int configure(int actualValue){
        
        if (message.length!=option.length){
            IJ.log("ERROR: configuration, option and message size are different");
        }
        
        
        String path=IJ.getDirectory("plugins");
        String jarpath =path+fileName;
        File file = new File(jarpath);
        if (!file.exists()){
            actualValue=0;
        }
        
        
        int numBerOfOption=option.length;
        
        int newValue=actualValue;
        
        
        GenericDialog gd = new GenericDialog("ZOLA: configuration");
        
        
        
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
        try{
            fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
        }catch(Exception e){}
        
        boolean [] bool = new boolean[numBerOfOption];
        
        
        for (int i=0;i<numBerOfOption;i++){
            
            if ((actualValue&(int)Math.pow(2, i))==(int)Math.pow(2, i)){
                bool[i]=true;
            }
        }
        gd.addMessage("Choose the tools you want to add to ZOLA",fontBold);
        for (int i=0;i<numBerOfOption;i++){
            gd.addCheckbox(message[i], bool[i]);
        }
        
        
        gd.showDialog();
        
        boolean [] boolResult = new boolean[numBerOfOption];
        
        if (!gd.wasCanceled()){
            newValue=0;
            for (int i=0;i<numBerOfOption;i++){
                boolResult[i] = (boolean)gd.getNextBoolean();
                
                if (boolResult[i]){
                    newValue+=Math.pow(2, i);
                    
                }
            }
            
        }
        
        modifPulginConfigFile(file,boolResult,actualValue,newValue);
        
            
            
        
        return newValue;
    }
    
    
    
    static void modifPulginConfigFile(File file,boolean [] bool,int actualValue, int newValue){
        
        

        
        
        
        
        
        if (newValue==0){
            if (file.exists()){
                file.delete();
                if (newValue!=actualValue){
                    new WaitForUserDialog("ZOLA configuration", "Please restart ImageJ or refresh menu (help -> refresh menu)").show();
                }
            }
        }    
        else if ((newValue!=actualValue)||(!file.exists())){
            JarOutputStream jos = null;
            try {
               OutputStream os = new FileOutputStream(file);
               jos = new JarOutputStream(os);
            } catch (IOException e) {
                  e.printStackTrace();
            }

            int len =0;
            byte[][] buffer = new byte[option.length][];
            for (int i=0;i<buffer.length;i++){
                if (bool[i]){
                    buffer[i]=option[i].getBytes();
                }
            }
            try {
                JarEntry je = new JarEntry("plugins.config");
                jos.putNextEntry(je);
                for (int i=0;i<buffer.length;i++){
                    if (bool[i]){
                        jos.write(buffer[i]);
                    }
                }
                jos.closeEntry();
                jos.close();
            } catch (IOException e) {
                  e.printStackTrace();
            }
            
            new WaitForUserDialog("ZOLA configuration", "Please restart ImageJ or refresh menu (help -> refresh menu)").show();
        }
        
        
    }
    
    
}
