/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

import ij.IJ;
import java.io.File;

/**
 *
 * @author benoit
 */
public class Protocol_localization {
    
    public static String pathCalib="";
    public static String pathLoc;
    public static double time_minut;
    public static double focus;
    public static double refIndexSampMedium;
    public static int frameNumber;
    public static int frameNumberReconstructed;
    public static boolean abortedReconstruction;
    public static int localizationNumber;
    public static int patchSize;
    public static boolean isSCMOS;
    public static double adu;
    public static double gain;
    public static double offset;
    public static String path_SCMOSvariance;
    public static String path_SCMOSgain; 
    public static String path_SCMOSoffset;
    
    public static void saveProtocolJSON(String path){
        Protocol_localization.saveProtocolJSON(path,pathCalib, pathLoc,time_minut,focus,refIndexSampMedium,frameNumber,frameNumberReconstructed,abortedReconstruction,localizationNumber,patchSize,isSCMOS, adu, gain, offset, path_SCMOSvariance, path_SCMOSgain, path_SCMOSoffset);
    }
    
    
    public static void saveProtocolJSON(String path,String pathCalib, String pathLoc,double time_minut, double focus,double refIndexSampMedium,int frameNumber,int frameNumberReconstructed,boolean abortedReconstruction,int localizationNumber,int patchSize,boolean isSCMOS, double adu, double gain, double offset, String path_SCMOSvariance, String path_SCMOSgain, String path_SCMOSoffset){
        if (path!=null){
            
            

            JSONformat obj = new JSONformat();
            //JSONObject obj = new JSONObject();
            obj.put("path_calibration",pathCalib);
            obj.put("path_localization",pathLoc);
            obj.put("reconstruction_time(min)",time_minut);
            obj.put("focus",focus);
            obj.put("refractive_index_sample_medium",refIndexSampMedium);
            obj.put("frame_number_total",frameNumber);
            obj.put("frame_number_reconstructed",frameNumberReconstructed);
            obj.put("aborted_reconstruction",abortedReconstruction);
            obj.put("localization_number",localizationNumber);
            obj.put("patch_size",patchSize);
            obj.put("is_SCMOS",isSCMOS);
            if (!isSCMOS){
                obj.put("camera_adu",adu);
                obj.put("camera_gain",gain);
                obj.put("camera_offset",offset);
                
            }
            else{
                obj.put("camera_SCMOS_path_gain",path_SCMOSgain);
                obj.put("camera_SCMOS_path_offset",path_SCMOSoffset);
                obj.put("camera_SCMOS_path_variance",path_SCMOSvariance);
            }
            

            obj.save(path);

        }
    }
    
    
    
    
    
    public static void loadProtocolJSON(String path){
        if (path!=null){
            
            try{
                File f = new File(path);
                
                if (!f.exists()){
                    IJ.log("ERROR: file does not exists: "+path);
                }
            }catch(Exception eee){
                IJ.log("ERROR: file may not exists: "+path);
            }
            
            
            JSONformat obj= new JSONformat();

            obj.load(path);
            

            
            pathCalib=obj.get("path_calibration");
            pathLoc=obj.get("path_localization");
            focus = Double.parseDouble(obj.get("focus").toString());
            refIndexSampMedium = Double.parseDouble(obj.get("refractive_index_sample_medium").toString());
            frameNumber = Integer.parseInt(obj.get("frame_number_total").toString());
            frameNumberReconstructed = Integer.parseInt(obj.get("frame_number_reconstructed").toString());
            time_minut = Double.parseDouble(obj.get("reconstruction_time(min)").toString());
            abortedReconstruction = Boolean.parseBoolean(obj.get("aborted_reconstruction").toString());
            localizationNumber = Integer.parseInt(obj.get("localization_number").toString());
            patchSize = Integer.parseInt(obj.get("patch_size").toString());
            isSCMOS = Boolean.parseBoolean(obj.get("is_SCMOS").toString());
            
            if (!isSCMOS){
                adu = Double.parseDouble(obj.get("camera_adu").toString());
                gain = Double.parseDouble(obj.get("camera_gain").toString());
                offset = Double.parseDouble(obj.get("camera_offset").toString());
                
                
            }
            else{
                path_SCMOSgain=obj.get("camera_SCMOS_path_gain");
                path_SCMOSoffset=obj.get("camera_SCMOS_path_offset");
                path_SCMOSvariance=obj.get("camera_SCMOS_path_variance");
                
            }
            


        }
    }
    
    
    
    
    
    
    
}
