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
            
        
        "Plugins>ZOLA>additional tools, \"Best PSF model optimization\", org.pasteur.imagej.ZOLA(\"zola_phaseoptim\")\n",
        
        "Plugins>ZOLA>additional tools>Simulation>PSF generator, \"Initialisation\", org.pasteur.imagej.ZOLA(\"zola_getModels_GPUalloc\")\n" +
"Plugins>ZOLA>additional tools>Simulation>PSF generator, \"computation\", org.pasteur.imagej.ZOLA(\"zola_getModels_GPUcomput\")\n" +
"Plugins>ZOLA>additional tools>Simulation>PSF generator, \"freeMemory\", org.pasteur.imagej.ZOLA(\"zola_getModels_GPUfree\")\n",
        
        "Plugins>ZOLA>Dual camera, \"Simple Registration\", org.pasteur.imagej.ZOLA(\"zola_dual_Cam_registration\")\n" +
"Plugins>ZOLA>Dual camera, \"Fusion\", org.pasteur.imagej.ZOLA(\"zola_dual_Cam_localization\")\n" +
"Plugins>ZOLA>Dual camera, \"CRLB_\", org.pasteur.imagej.ZOLA(\"zola_dualCamCRLB\")\n",
        
        
        
        
        "Plugins>ZOLA>Dual color, \"Compute chromatic aberrations\", org.pasteur.imagej.ZOLA(\"zola_computeChromaticAberrations\")\n" +
"Plugins>ZOLA>Dual color, \"Compute SR 2 LR registration\", org.pasteur.imagej.ZOLA(\"zola_computeSR2LRregistration\")\n" +
"Plugins>ZOLA>Dual color, \"Color Registration\", org.pasteur.imagej.ZOLA(\"zola_colorRegistration\")\n",
    
        
        "Plugins>ZOLA>Drift correction, \"3D Fiducial-markers\", org.pasteur.imagej.ZOLA(\"zola_driftCorrectionFiducial\")\n",
        
        "Plugins>ZOLA>additional tools>statistic, \"localization histogram\", org.pasteur.imagej.ZOLA(\"zola_stathisto\")\n"+
        "Plugins>ZOLA>additional tools>statistic, \"measure mean attach time\", org.pasteur.imagej.ZOLA(\"zola_statAttachTime1\")\n",
        
        "Plugins>ZOLA>Visualization, \"frame color coding\", org.pasteur.imagej.ZOLA(\"zola_timeRendering\")\n",
        
    "Plugins>ZOLA>dev, \"Multi-emitter localization\", org.pasteur.imagej.ZOLA(\"zola_multi-emitter localization\")\n" +
"Plugins>ZOLA>dev, \"test\", org.pasteur.imagej.ZOLA(\"zola_test\")\n" +
"Plugins>ZOLA>dev>Simulation, \"Simulation\", org.pasteur.imagej.ZOLA(\"zola_simulation\")\n" +
"Plugins>ZOLA>dev>Simulation, \"Simulation2beads\", org.pasteur.imagej.ZOLA(\"zola_simulation2beads\")\n" +
"Plugins>ZOLA>dev, \"CRLBfromFile\", org.pasteur.imagej.ZOLA(\"zola_crlb_FromFile\")\n" +
"Plugins>ZOLA>dev, \"CRLBfromFileDualObj\", org.pasteur.imagej.ZOLA(\"zola_crlb_FromFileDualObj\")\n" +
"Plugins>ZOLA>dev, \"Photon conversion wrong\", org.pasteur.imagej.ZOLA(\"zola_photonConversion\")\n"+
"Plugins>ZOLA>dev>Resolution, \"points Measurment\", org.pasteur.imagej.ZOLA(\"zola_pointResolution\")\n" +
"Plugins>ZOLA>dev>Resolution, \"filament Measurment\", org.pasteur.imagej.ZOLA(\"zola_filamentResolution\")"};
    
    
    static String [] message = {
        "Phase_optimization :",
        "PSF_simulation :",
        "Dual_cam_localization (beta) :",
        "Dual_color_registration (beta) :",
        "Drift correction using fiducial markers (alpha) :",
        "Statistics: histograms/mean particle attach time :",
        "Visualization: frame color coding :",
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
            
            new WaitForUserDialog("ZOLA configuration", "Please restart ImageJ").show();
        }
        
        
    }
    
    
}
