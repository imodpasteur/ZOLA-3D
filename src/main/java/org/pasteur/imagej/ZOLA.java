/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej;



import org.pasteur.imagej.process.*;
import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.postprocess.*;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.cuda.*;



import java.awt.Checkbox;
import ij.IJ;
import ij.WindowManager;
import ij.Prefs;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Vector;
import java.util.ArrayList;
import java.awt.TextField; 
import java.awt.event.MouseEvent; 
import java.awt.event.MouseListener; 
import java.awt.event.TextEvent; 
import java.awt.event.TextListener; 
import javax.swing.JFileChooser;
import ij.process.FloatPolygon;
import javax.swing.JOptionPane;
import ij.gui.WaitForUserDialog;
import ij.plugin.frame.RoiManager;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaError;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.YesNoCancelDialog;
import ij.process.FloatProcessor;
import java.awt.Button;
import java.awt.Choice;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;


/**
 *
 * @author benoit
 */
public class ZOLA implements PlugIn  {
    
    boolean isSCMOS=false;
    
    int concecutiveFrameThreshold=1000;
    
    int offFrame=30;
    
    int smoothingFrameNumber=100;
            
            
    int maxFFT=140;
    double dualCamMaxDistanceMergingPlan=200;
    
    
    double maxDistanceMergingPlan=100;
    
    double maxDistanceMergingZ=200;
    
    
    int sizeTextString=45;
    
    Prefs prefs = new Prefs(); 
    
    
    
    String path_SCMOSvariance="";
    String path_SCMOSoffset="";
    String path_SCMOSgain="";
    
    String path_wobble="";
    
    String path_localization="";
    String path_result="";
    String path_localization2="";
    String path_calibration="";
    
    String path_registration="";
    
    
    int model_number=100;
    int stream_number=1;
            
    int shiftrendering;
    boolean showLUT;
    boolean is3Drendering;
    boolean isCOLORrendering;
    double binSize;
    int sizePatch;
    int sizePatchPhaseRet;
    int nbThread;
    double minZ;
    double maxZ;
    double axialRange;
    double stepZloc;
    int photonThreshold;
    int nbStream;
    double sigma=.9;
    
    double parameterFit1=1;
    double parameterFit2=1;
    
    double xystep;
    double zstep;
    
    double na;
    double noil;
    double nwat;
    double zfocus;
    double wavelength;
    
    //double rescale_slope;
    //double rescale_intercept;
    
    double adu;//Analog2DigitalUnit
    double gain;//Analog2DigitalUnit
    double offset;//Analog2DigitalUnit
    
    double radius=5;
    int axialside=0;
    int zernikeCoef=1;
    int orderFit=1;
    int iterationNumber=15;
    
    int frameSplit=100000;
    
    int idVariable=0;
    int binNumberVariable=100;
    
    double drift_sub_image_sizeUm=6;
    double drift_pixelsizeNM=100;
    int drift_bin=3;
    
    boolean GPU_computation;
       
    boolean withApoFactor;
    
    double sizeRendering=20;
            
    double maximumDistance;
    
    static org.pasteur.imagej.process.gpu.Model_generator_ mg;
    
    static StackLocalization sl;
    static Boolean processing;
    
    
    
    double [][] scmosvargain=null;
    double [][] scmosvariance=null;
    double [][] scmosgain=null;
    double [][] scmosoffset=null;
    
    int configure=0;
    
    public ZOLA(){
        
        
    }
    
    
    
    
    public void run(String command){
        
        File dir = new File(IJ.getDirectory("plugins"));
        File[] filesList = dir.listFiles();
        for (File file : filesList) {
            if (file.isFile()) {
                String start="ZOLA_-";
                if (file.getName().startsWith("ZOLA_-")){
                    int ends=file.getName().indexOf(".jar");
                    IJ.log("ZOLA version: "+file.getName().substring(start.length(), ends));
                }
            }
        }

         
         
        
        
        configure=(int)prefs.get("Zola.configure", configure);
        
        path_result=prefs.get("Zola.path_result", IJ.getDirectory("image")+"path_result.csv");
        
        path_wobble=prefs.get("Zola.pathwobble", IJ.getDirectory("image")+"calibration_wobble.csv");
        
        path_localization=prefs.get("Zola.pathlocalization", IJ.getDirectory("image")+"localization_table.csv");
        path_localization2=prefs.get("Zola.pathlocalization2", IJ.getDirectory("image")+"localization_table.csv");
        path_calibration=prefs.get("Zola.pathcalibration", IJ.getDirectory("image")+"calibration_table.json");
        path_registration=prefs.get("Zola.path_registration", IJ.getDirectory("image")+"calibration_table.json");
        
        path_SCMOSvariance=prefs.get("Zola.path_SCMOSvariance", "");
        path_SCMOSoffset=prefs.get("Zola.path_SCMOSoffset", "");
        path_SCMOSgain=prefs.get("Zola.path_SCMOSgain", "");
        
        //rescale_slope=prefs.get("Zola.rescaleslope", 1);
        //rescale_intercept=prefs.get("Zola.rescaleintercept", 0);
        
        adu=prefs.get("Zola.adu", 1);
        gain=prefs.get("Zola.gain", 1);
        offset=prefs.get("Zola.offset", 0);
        model_number=(int)prefs.get("Zola.model_number", 100);
        stream_number=(int)prefs.get("Zola.stream_number", 1);
        
        frameSplit=(int)prefs.get("Zola.frameSplit", 1);
        
        shiftrendering=(int)prefs.get("Zola.shiftrendering", 1);
        sizePatch=(int)prefs.get("Zola.sizePatchLoc", 32);
        sizePatchPhaseRet=(int)prefs.get("Zola.sizePatchPhaseRet", 64);
        concecutiveFrameThreshold=(int)prefs.get("Zola.concecutiveFrameThreshold", 1000);
        smoothingFrameNumber=(int)prefs.get("Zola.smoothingFrameNumber", 100);
        nbThread=(int)prefs.get("Zola.nbThread", 100);
        axialRange=prefs.get("Zola.axialRange", 4);
        minZ=prefs.get("Zola.minZ", -2);
        maxZ=prefs.get("Zola.maxZ", 2);
        stepZloc=prefs.get("Zola.stepZloc", .05);
        photonThreshold=(int)prefs.get("Zola.photonThreshold", 500);
        nbStream=(int)prefs.get("Zola.nbStream", 4);
        
        binSize=prefs.get("Zola.binSize", 10);
        
        radius=prefs.get("Zola.radius", 5);
        
        isSCMOS=(boolean)((prefs.get("Zola.isSCMOS", false)));
        
        is3Drendering=(boolean)((prefs.get("Zola.is3Drendering", true)));
        isCOLORrendering=(boolean)((prefs.get("Zola.isCOLORrendering", true)));
        
        showLUT=(boolean)((prefs.get("Zola.showLUT", true)));
        
        
        maximumDistance=prefs.get("Zola.maximumDistance", 200);
        
        xystep=prefs.get("Zola.xystep", .11);
        zstep=prefs.get("Zola.zstep", .05);
        
        na=prefs.get("Zola.na", 1.2);
        noil=prefs.get("Zola.noil", 1.33);
        nwat=prefs.get("Zola.nwat", 1.33);
        zfocus=prefs.get("Zola.zfocus", 0);
        sigma=prefs.get("Zola.sigma", 0.9);
        parameterFit1=prefs.get("Zola.parameterFit1", 1);
        parameterFit2=prefs.get("Zola.parameterFit2", 1);
        wavelength=prefs.get("Zola.wavelength", .64);
        axialside=(int)prefs.get("Zola.axialside", 0);
        zernikeCoef=(int)prefs.get("Zola.zernikeCoef", 1);
        
        idVariable=(int)prefs.get("Zola.idVariable", 0);
        binNumberVariable=(int)prefs.get("Zola.binNumberVariable", 100);
        
        
        
    
    
        iterationNumber=(int)prefs.get("Zola.iterationNumber", 30);
        
        drift_sub_image_sizeUm=prefs.get("Zola.driftsubimagesizeUm", 6);
        drift_pixelsizeNM=prefs.get("Zola.driftpixelsizeNM", 50);
        drift_bin=(int)prefs.get("Zola.driftbin", 3);
        
        
        GPU_computation=(boolean)prefs.get("Zola.GPU_computation", false);
        withApoFactor=(boolean)prefs.get("Zola.withApoFactor", false);
        sizeRendering=prefs.get("Zola.sizeRendering", 20);
        
        dualCamMaxDistanceMergingPlan=prefs.get("Zola.dualCamMaxDistanceMergingPlan", 250);
        
        
        maxDistanceMergingPlan=prefs.get("Zola.maxDistanceMergingPlan", 100);
        
        maxDistanceMergingZ=prefs.get("Zola.maxDistanceMergingZ", 200);
        
        
        
        
        
//        IJ.log("library path: "+System.getProperty("java.library.path"));
//        IJ.log("should we add... : "+IJ.getDirectory("startup")+""+"lib");
        
        //IJ.log("System  :"++"       "+System.getProperty("user.dir")+"/lib");
        
        
        /*
        String lib=System.getProperty("java.library.path")+":"+System.getProperty("user.dir")+"/lib";
        
        IJ.log("system before "+System.getProperty("java.library.path"));
        
        IJ.log("lib:"+lib);
        System.setProperty("java.library.path",lib);
        IJ.log("system after "+System.getProperty("java.library.path"));*/
        
        /*System.load(System.getProperty("user.dir")+"/lib/JCudaRuntime");*/
        if (path_wobble==null){
            path_wobble="";
        }
        if (path_localization==null){
            path_localization="";
        }
        if (path_localization2==null){
            path_localization2="";
        }
        if (path_calibration==null){
            path_calibration="";
        }
        if (path_registration==null){
            path_registration="";
        }
        if (path_result==null){
            path_result="";
        }
        
        
        if (command.startsWith("zola_calibration")){
            
            calibration();
            
        }
        
        if (command.startsWith("zola_wobblecalibration")){
            
            wobblecalibration();
            
        }
        if (command.startsWith("zola_wobblecorrection")){
            
            wobblecorrection();
            
        }
        
        else if (command.startsWith("zola_filteringHistogram")){
            
            stathistoFilter();
            
            
        }
        
        else if (command.startsWith("zola_filtering")){
            
            filtering();
            
        }
        else if (command.startsWith("zola_removeLocalizations")){
            
            inversefiltering();
            
        }
        else if (command.startsWith("zola_deconvolution")){
            
            deconvolution();
            
        }
        else if (command.startsWith("zola_crlb_FromFile")){
            
            crlbfromfile();
            
            
        }
        else if (command.startsWith("zola_crlb")){
            
            crlb();
            
            
        }
        
        else if (command.startsWith("zola_stathisto")){
            
            stathisto();
            
            
        }
        
        else if (command.startsWith("zola_setcameraEMCCD")){
            
            setCameraEMCCD();
            
            
        }
        else if (command.startsWith("zola_setcameraSCMOS")){
            
            setCameraSCMOS();
            
            
        }
        else if (command.startsWith("zola_statmean")){
            
            statmean();
            
            
        }
        
        
        else if (command.startsWith("zola_condensation")){
            
            condensation();
            
            
        }
        
        else if (command.startsWith("zola_statAttachTime1")){
            
            statAttachTime1();
            
            
        }
        
        
        else if (command.startsWith("zola_test")){
            
            test();
            
            
        }
        
        else if (command.startsWith("zola_driftCorrectionFiducialMarker")){
            
            driftCorrectionFiducialMarker();
            
            
        }
        
        else if (command.startsWith("zola_driftCorrection")){
            
            driftCorrection();
            
            
        }
        
        
        
        else if (command.startsWith("zola_resetDrift")){
            
            resetDriftCorrection();
            
            
        }
        
        else if (command.startsWith("zola_setDrift")){
            
            setDriftCorrection();
            
            
        }
        
        
        else if (command.startsWith("zola_localization_test")){
            
            localizationtest();
            
            
            
        }
        
        else if (command.startsWith("zola_localization")){
            
            localization();
            
            
            
        }
        else if (command.startsWith("zola_configure")){
            
            configure();
            
            
            
        }
        
        
        
        else if (command.startsWith("zola_colorHist")){
            
            colorHist();
            
        }
        
        else if (command.startsWith("zola_hist")){
            
            hist();
            
        }
        
        else if (command.startsWith("zola_timeRendering")){
            
            timeRender();
            
        }
        
           
                
        else if (command.startsWith("zola_scatterPlot")){
            
            scatterPlot();
            
        }
        else if (command.startsWith("zola_colorize")){
            
            colorizeHist();
            
        }
        
        
        
        else if (command.startsWith("zola_import")){
            importData();
            
        }
        else if (command.startsWith("zola_append")){
            appendData();
            
        }
        
        else if (command.startsWith("zola_export")){
            exportData();
            
        }
        
        
        else if (command.startsWith("zola_phaseoptim")){
            phaseoptim();
            
        }
        
        
        else if (command.startsWith("zola_getModels_GPUalloc")){
            getModels_GPUalloc();
        }
        else if (command.startsWith("zola_getModels_GPUcomput")){
            getModels_GPUcomput();
        }
        else if (command.startsWith("zola_getModels_GPUfree")){
            getModels_GPUfree();
        }
        
        else if (command.startsWith("zola_dual_Cam_localization")){
            
            dualcamlocalization();
        }
        
        else if (command.startsWith("zola_dual_Cam_registration")){
            
            dualcamregistration();
        }
        
        else if (command.startsWith("zola_multi-emitter localization")){
            
            multiEmitterLocalization();
        }
        else if (command.startsWith("zola_mergeframe")){
            
            mergeframe();
        }
        
        else if (command.startsWith("zola_test")){
            
            test();
        }
        else if (command.startsWith("zola_addPoissonNoise")){
            
            addPoissonNoise();
        }
        else if (command.startsWith("zola_simulation2beads")){
            
            simulation2beads();
        }
        else if (command.startsWith("zola_simulationOverlappingBeads")){
            
            simulationOverlappingBeads();
        }
        else if (command.startsWith("zola_simulation")){
            
            simulation();
        }
        
        else if (command.startsWith("zola_crlb_FromFile")){
            
            crlbfromfile();
        }
        else if (command.startsWith("zola_crlb_FromFileDualObj")){
            
            crlbfromfileDualObj();
        }
        else if (command.startsWith("zola_photonConversion")){
            
            IJ.log("not implemented");
        }
        else if (command.startsWith("zola_pointResolution")){
            
            pointResolution();
        }
        else if (command.startsWith("zola_filamentResolution")){
            
            filamentResolution();
        }
        
        else if (command.startsWith("zola_computeChromaticAberrations")){
            
            computeChromaticAberrations();
        }
        else if (command.startsWith("zola_computeSR2LRregistration")){
            
            computeSR2LRregistration();
        }
        else if (command.startsWith("zola_colorRegistration")){
            
            colorRegistration();
        }
        
        
        
        
        
        
        
        /*"Plugins>ZOLA>dev, \"Multi-emitter localization\", org.pasteur.imagej.ZOLA(\"zola_multi-emitter localization\")\n" +
"Plugins>ZOLA>dev, \"test\", PasteurLocalization.ZOLA(\"zola_test\")\n" +
"Plugins>ZOLA>dev>Simulation, \"Simulation\", PasteurLocalization.ZOLA(\"zola_simulation\")\n" +
"Plugins>ZOLA>dev>Simulation, \"Simulation2beads\", PasteurLocalization.ZOLA(\"zola_simulation2beads\")\n" +
"Plugins>ZOLA>dev, \"CRLBfromFile\", PasteurLocalization.ZOLA(\"zola_crlb_FromFile\")\n" +
"Plugins>ZOLA>dev, \"CRLBfromFileDualObj\", PasteurLocalization.ZOLA(\"zola_crlb_FromFileDualObj\")\n" +
"Plugins>ZOLA>dev, \"Photon conversion wrong\", PasteurLocalization.ZOLA(\"zola_photonConversion\")\n"+
"Plugins>ZOLA>dev>Resolution, \"points Measurment\", PasteurLocalization.ZOLA(\"zola_pointResolution\")\n" +
"Plugins>ZOLA>dev>Resolution, \"filament Measurment\", PasteurLocalization.ZOLA(\"zola_filamentResolution\")"};*/
        
        prefs.set("Zola.configure", configure);
        prefs.set("Zola.path_SCMOSvariance", path_SCMOSvariance);
        prefs.set("Zola.path_SCMOSgain", path_SCMOSgain);
        prefs.set("Zola.path_SCMOSoffset", path_SCMOSoffset);
        prefs.set("Zola.GPU_computation", GPU_computation);
        prefs.set("Zola.withApoFactor", withApoFactor);
        prefs.set("Zola.pathwobble", path_wobble);
        prefs.set("Zola.pathlocalization", path_localization);
        prefs.set("Zola.pathlocalization2", path_localization2);
        prefs.set("Zola.path_registration", path_registration);
        if (path_calibration.length()>2){
            prefs.set("Zola.pathcalibration", path_calibration);
        }
        prefs.set("Zola.frameSplit", frameSplit);
        prefs.set("Zola.model_number", model_number);
        prefs.set("Zola.stream_number", stream_number);
        prefs.set("Zola.shiftrendering", shiftrendering);
        prefs.set("Zola.sizePatchLoc", sizePatch);
        prefs.set("Zola.sizePatchPhaseRet", sizePatchPhaseRet);
        prefs.set("Zola.concecutiveFrameThreshold", concecutiveFrameThreshold);
        prefs.set("Zola.smoothingFrameNumber", smoothingFrameNumber);
        prefs.set("Zola.nbThread", nbThread);
        prefs.set("Zola.minZ", minZ);
        prefs.set("Zola.maxZ", maxZ);
        prefs.set("Zola.axialRange", axialRange);
        
        prefs.set("Zola.stepZloc", stepZloc);
        prefs.set("Zola.photonThreshold", photonThreshold);
        prefs.set("Zola.nbStream", nbStream);
        
        prefs.set("Zola.binSize", binSize);
        
        
        
        prefs.set("Zola.xystep", xystep);
        prefs.set("Zola.zstep", zstep);
        
        prefs.set("Zola.radius", radius);
        
        prefs.set("Zola.adu", adu);
        prefs.set("Zola.gain", gain);
        prefs.set("Zola.offset", offset);
        
        prefs.set("Zola.isSCMOS", isSCMOS);
        
        
        prefs.set("Zola.is3Drendering", is3Drendering);
        prefs.set("Zola.showLUT", showLUT);
        prefs.set("Zola.isCOLORrendering", isCOLORrendering);
        prefs.set("Zola.na", na);
        prefs.set("Zola.noil", noil);
        prefs.set("Zola.nwat", nwat);
        prefs.set("Zola.sigma", sigma);
        prefs.set("Zola.parameterFit1", parameterFit1);
        prefs.set("Zola.parameterFit2", parameterFit2);
        prefs.set("Zola.zfocus", zfocus);
        prefs.set("Zola.wavelength", wavelength);
        
        prefs.set("Zola.maximumDistance", maximumDistance);
        prefs.set("Zola.axialside", axialside);
        prefs.set("Zola.zernikeCoef", zernikeCoef);
        prefs.set("Zola.iterationNumber", iterationNumber);
        prefs.set("Zola.driftsubimagesizeUm", drift_sub_image_sizeUm);
        prefs.set("Zola.driftpixelsizeNM", drift_pixelsizeNM);
        prefs.set("Zola.driftbin", drift_bin);
        
        prefs.set("Zola.idVariable", idVariable);
        prefs.set("Zola.binNumberVariable", binNumberVariable);
        
        prefs.set("Zola.sizeRendering", sizeRendering);
        
        
        prefs.set("Zola.dualCamMaxDistanceMergingPlan", dualCamMaxDistanceMergingPlan);
        
        
        prefs.set("Zola.maxDistanceMergingPlan", maxDistanceMergingPlan);
        
        prefs.set("Zola.maxDistanceMergingZ", maxDistanceMergingZ);      
                
        prefs.savePreferences();
    }
    
    
    
    
    
    
    
    
    void configure(){
        
        configure=org.pasteur.imagej.configuration.ConfigFile.configure(configure);
        
        
    }
    
    
    
    
    
    
    
    void localization(){
        
        
        nbStream=5;
        nbThread=100;  
            
        ImagePlus imp=IJ.getImage();


        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
            processing=new Boolean(true);
            
            
            if (imp==null){
                IJ.log("no image opened");
                processing=new Boolean(false);
            }
            /*else if ((imp.getWidth()>180)||(imp.getHeight()>180)){
                IJ.log("image should be < 180 pixels for the moment");
            }*/
            else{
                
                
                Roi r=imp.getRoi();
                int width=imp.getWidth();
                int height=imp.getHeight();
                
                boolean scmosLoad=false;
                if (isSCMOS){
                    scmosLoad=loadSCMOScameraFiles(width, height, r);
                }
                
                
                
                String pathdir=IJ.getDirectory("image");
                if (pathdir==null){
                    pathdir=IJ.getDirectory("current");
                }
                String path=pathdir;
                
                path=path+"ZOLA_localization_table.csv";
                
                
                
                
                GenericDialog gd=null;
                
                //since NonBlockingGenericDialog incompatible with headless mode:
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())){
                    gd = new NonBlockingGenericDialog("ZOLA: Localization");
                }
                else{
                    gd = new GenericDialog("ZOLA: Localization");
                }

                    Font font = gd.getFont();
                    Font fontBold= gd.getFont();
                try{
                    fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
                }catch(Exception e){}

                    gd.addCheckbox("Run_on_GPU", GPU_computation);






                    if (!isSCMOS){
                        gd.addMessage("Camera is EMCCD",fontBold);
                        gd.addMessage("ADU = "+adu);
                        gd.addNumericField("Gain = ", gain, 1,6,"(1 if camera provides photon counts)");
                        gd.addMessage("Offset = "+offset);

                    }
                    else{
                        if (scmosLoad){
                            gd.addMessage("Camera is SCMOS",fontBold);
                        }
                        else{
                            gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                        }
    //                    String sep=File.separator;
    //                    String [] var=this.path_SCMOSvariance.split(sep);
    //                    String [] off=this.path_SCMOSoffset.split(sep);
    //                    String [] gain=this.path_SCMOSgain.split(sep);
    //                    gd.addMessage("Offset file = "+off[var.length-1]);
    //                    gd.addMessage("Variance file = "+var[var.length-1]);
    //                    gd.addMessage("Gain file = "+gain[var.length-1]);

                        gd.addMessage("Offset file = "+path_SCMOSoffset);
                        gd.addMessage("Variance file = "+path_SCMOSvariance);
                        gd.addMessage("Gain file = "+path_SCMOSgain);

                    }



                    gd.addStringField("PSF_calibration_file:",path_calibration,sizeTextString); 



                    gd.addMessage("Sample parameters", fontBold);

                    gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");

                    gd.addNumericField("Distance_focus_to_coverslip", zfocus, 3,6,"(µm)");



                    gd.addMessage("Advanced parameter", fontBold);


                    gd.addNumericField("Patch_size: ", sizePatch,0,6,"(pixels)");

                    //gd.addNumericField("number_of_streams: ", nbStream,0);
                    //gd.addNumericField("number_of_threads: ", nbThread,0);





                    gd.addNumericField("Expected_axial_range: ", axialRange,2,6,"(µm)");

                    gd.addNumericField("Min_number_of_photons: ", photonThreshold,0);


                    gd.addStringField("Localization_table:",path,sizeTextString); 

                    Choice choicefield=null;



                    TextField textField; 
                    TextField textField2; 
                    // Add a mouse listener to the config file field 
                    if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
                    {
                        int t=0;
                        //textField
                        Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

                        textField = texts.get(t++); 
                        MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
                        textField.addMouseListener(mol); 

                        textField2 = texts.get(t++); 
                        MouseOptionSave mos=new MouseOptionSave(textField2,pathdir,"Export localization result table");
                        textField2.addMouseListener(mos); 



                        //NumericField
                        Vector<TextField> nums = (Vector<TextField>) gd.getNumericFields(); 



                        PreviewButton but = new PreviewButton(imp,texts,nums); 

                        gd.add(but); 

                        AutoFocusButton afbut = new AutoFocusButton(imp,texts,nums); 

                        gd.add(afbut); 

                    }


                    gd.showDialog();





                    if (!gd.wasCanceled()){
                        GPU_computation = (boolean)gd.getNextBoolean();

                        if (!isSCMOS){
                            gain = (double)gd.getNextNumber();
                        }

                        path_calibration = gd.getNextString(); 

                        nwat= (double)gd.getNextNumber();

                        zfocus= (double)gd.getNextNumber();

                        sizePatch = (int)gd.getNextNumber();
                        //nbStream = (int)gd.getNextNumber();
                        //nbThread = (int)gd.getNextNumber();

                        axialRange= (double)gd.getNextNumber();




                        photonThreshold = (int)gd.getNextNumber();

                        path_localization = gd.getNextString(); 
                    }
                    else{
                        processing=new Boolean(false);
                        return;
                    }
                
                
                
                int nnn=(path_localization.lastIndexOf("."));
                if ((nnn>=0)&&(nnn>path_localization.length()-6)){
                    path_localization=path_localization.substring(0,nnn);
                }
                path_localization+=".csv";
            
                if (path_localization.equals(path_calibration)){
                    IJ.log("Process stopped: the localization path has to be different than calibration path");
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                if (sizePatch%2==1){
                    sizePatch++;
                }
                
                
                
                
                
                
                
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                }
                if ((width<sizePatch)||(height<sizePatch)){
                    IJ.log("Process stopped: Localization impossible. The image, or selected ROI should be larger than patch size");
                    IJ.log("Please, select a large rectangle in the image before running the localization");
                    IJ.log("image width:"+width+" ; height:"+height+" ; sizePatch:"+sizePatch);
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                
                
                
                if (gain==0||adu==0){
                    new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
                }
                else{
                    
                    
                    
                    
                    int sizeFFT=128;
                    sizeFFT=Math.max(sizeFFT,sizePatch*2);
                    sizeFFT=Math.min(sizeFFT,maxFFT);
                    
                    
                    if (GPU_computation){
                        MyCudaStream.init(nbStream+1);
                        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(nwat);
                        dp.param.Zfocus=zfocus;
                        //dp.param.weightZ=.9;
                        IJ.log("Localization started");


                        org.pasteur.imagej.process.gpu.LocalizationPipeline_ lp;
                        if (isSCMOS){
                            lp=new org.pasteur.imagej.process.gpu.LocalizationPipeline_(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain,false);
                        }
                        else{
                            lp=new org.pasteur.imagej.process.gpu.LocalizationPipeline_(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,adu, gain,offset,false);
                        }
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);


                        /*classes.testClass lp =new classes.testClass(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,rescale_slope, rescale_intercept);
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);*/



                        dp.free();
                    
                    
                        MyCudaStream.destroy();
                    }
                    else{
                        org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(nwat);
                        dp.param.Zfocus=zfocus;
                        
                        IJ.log("Localization started");
                        
                        org.pasteur.imagej.process.cpu.LocalizationPipeline lp;
                        if (isSCMOS){
                            
                            lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain,false);
                        }
                        else{
                            lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,adu, gain,offset,false);
                        }
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);
                        
                        
                    }
                    
                }
            }
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    void localizationtest(){
        
            
        ImagePlus imp=IJ.getImage();
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
            processing=new Boolean(true);
            
            
            if (imp==null){
                IJ.log("no image opened");
                processing=new Boolean(false);
            }
            /*else if ((imp.getWidth()>180)||(imp.getHeight()>180)){
                IJ.log("image should be < 180 pixels for the moment");
            }*/
            else{
                
                int select=imp.getCurrentSlice();
                IJ.log("localization is applied only on frame "+select);
                ImageProcessor ip = imp.getProcessor();
                ImagePlus impTest=new ImagePlus("test",ip);
                
                Roi r=imp.getRoi();
                int width=imp.getWidth();
                int height=imp.getHeight();
                
                boolean scmosLoad=false;
                if (isSCMOS){
                    scmosLoad=loadSCMOScameraFiles(width, height, r);
                }
                
                
                
                
                GenericDialog gd = new GenericDialog("ZOLA: Localization test");
                
                
                Font font = gd.getFont();
                Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
                
                
                
                
                
                
                
                
                if (!isSCMOS){
                    gd.addMessage("Camera is EMCCD",fontBold);
                    gd.addMessage("ADU = "+adu);
                    gd.addNumericField("Gain = ", gain, 1,6,"(1 if camera provides photon counts)");
                    gd.addMessage("Offset = "+offset);
                    
                }
                else{
                    if (scmosLoad){
                        gd.addMessage("Camera is SCMOS",fontBold);
                    }
                    else{
                        gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                    }
//                    String sep=File.separator;
//                    String [] var=this.path_SCMOSvariance.split(sep);
//                    String [] off=this.path_SCMOSoffset.split(sep);
//                    String [] gain=this.path_SCMOSgain.split(sep);
//                    gd.addMessage("Offset file = "+off[var.length-1]);
//                    gd.addMessage("Variance file = "+var[var.length-1]);
//                    gd.addMessage("Gain file = "+gain[var.length-1]);
                    
                    gd.addMessage("Offset file = "+path_SCMOSoffset);
                    gd.addMessage("Variance file = "+path_SCMOSvariance);
                    gd.addMessage("Gain file = "+path_SCMOSgain);
                    
                }
                
                
                
                gd.addStringField("PSF_calibration_file:",path_calibration,sizeTextString); 
                
                
                
                gd.addMessage("Sample parameters", fontBold);
                
                gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        
                gd.addNumericField("Distance_focus_to_coverslip", zfocus, 3,6,"(µm)");
                
                
                
                gd.addMessage("Advanced parameter", fontBold);
                
                
                gd.addNumericField("Patch_size: ", sizePatch,0,6,"(pixels)");
                
                //gd.addNumericField("number_of_streams: ", nbStream,0);
                //gd.addNumericField("number_of_threads: ", nbThread,0);
                
                
                
                
                
                gd.addNumericField("Expected_axial_range: ", axialRange,2,6,"(µm)");
                
                gd.addNumericField("Min_number_of_photons: ", photonThreshold,0);
                
                
                
        
                TextField textField; 
            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;
                
                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
                
                
                textField = texts.get(t++); 
                MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
                textField.addMouseListener(mol); 
                
                


            }
        
        
                gd.showDialog();
                
                
                
                
                
                if (!gd.wasCanceled()){
                    
                    
                    if (!isSCMOS){
                        gain = (double)gd.getNextNumber();
                    }
                    
                    path_calibration = gd.getNextString(); 
                    
                    nwat= (double)gd.getNextNumber();
                    
                    zfocus= (double)gd.getNextNumber();
                    
                    sizePatch = (int)gd.getNextNumber();
                    //nbStream = (int)gd.getNextNumber();
                    //nbThread = (int)gd.getNextNumber();
                    
                    axialRange= (double)gd.getNextNumber();
                    
                    
                            
                    
                    photonThreshold = (int)gd.getNextNumber();
                    
                }
                else{
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                
                
                
                
                if (sizePatch%2==1){
                    sizePatch++;
                }
                
                
                
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                }
                if ((width<sizePatch)||(height<sizePatch)){
                    IJ.log("Process stopped: Localization impossible. The image, or selected ROI should be larger than patch size");
                    IJ.log("Please, select a large rectangle in the image before running the localization");
                    IJ.log("image width:"+width+" ; height:"+height+" ; sizePatch:"+sizePatch);
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                
                
                
                if (gain==0||adu==0){
                    new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
                }
                else{
                    
                    
                    
                    
                    int sizeFFT=128;
                    sizeFFT=Math.max(sizeFFT,sizePatch*2);
                    sizeFFT=Math.min(sizeFFT,maxFFT);
                    
                    
                    
                    org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
                    if (!dp.loading){
                        IJ.log("impossible to load "+path_calibration);
                        processing=new Boolean(false);
                        return;
                    }
                    dp.setSizeoutput(sizePatch);
                    dp.setNwat(nwat);
                    dp.param.Zfocus=zfocus;
                    IJ.log("Localization started");

                    org.pasteur.imagej.process.cpu.LocalizationPipeline lp;
                    if (isSCMOS){

                        lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,"",scmosvariance,scmosgain,scmosoffset,scmosvargain,true);
                    }
                    else{
                        lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,"",adu, gain,offset,true);
                    }
                    sl=new StackLocalization();
                    sl=lp.detectParticles(impTest,sl);
                    
                    
                }
            }
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    void crlb(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
        
        
        double maxValuePlot=-1;
        
        String pathResult="";
        
        double photonNumber=3000;
        double background=50;
        
        String [] fieldSimul = {"molecule moves","objective moves"};
        
        GenericDialog gd = new GenericDialog("ZOLA: Fundamental limit (CRLB)");
        
        
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}        
                
        gd.addCheckbox("Run_on_GPU", GPU_computation);
        
        //String path_calibration2="";
        
        gd.addStringField("Calibration_file:",path_calibration,sizeTextString); 
        
        //gd.addStringField("Calibration_file_cam2 [optional]:",path_calibration2,sizeTextString); 
        
        
        
        
        gd.addMessage("Parameters",fontBold);

        gd.addNumericField("patch_size: ", sizePatchPhaseRet,0,6,"(pixels)");
        
        
        
        
        gd.addNumericField("Axial_range: ", axialRange,2,6,"(µm)");
        
        gd.addNumericField("Z_step: ", stepZloc*1000,1,6,"(nm)");
        
        gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        
        gd.addNumericField("Distance_focus_to_coverslip (µm)", zfocus, 3);
        
        
        gd.addNumericField("Photon_number: ", photonNumber,0);
        gd.addNumericField("Background_intensity: ", background,0);
        
        //gd.addStringField("Plot_max_value (µm) [optional]: ", "",6);
        
        //gd.addStringField("Result_table [optional]:",pathResult,sizeTextString); 
        
        TextField textField; 
        TextField textField2; 
        TextField textField3; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
            //textField3 = texts.get(t++); 
            //mol=new MouseOptionLoad(textField3,path_calibration,"Import cal. file cam 2 (optional)");
            //textField3.addMouseListener(mol); 
//            textField2 = texts.get(t++); 
//            textField3 = texts.get(t++); 
//            MouseOptionSave mos=new MouseOptionSave(textField3,path_calibration,"Export CRLB result file");
//            textField3.addMouseListener(mos); 
            
            
        }
        
        
        gd.showDialog();



        if (!gd.wasCanceled()){
            
            GPU_computation = (boolean)gd.getNextBoolean();
            
            path_calibration = gd.getNextString(); 
            //path_calibration2 = gd.getNextString(); 
            
            sizePatchPhaseRet = (int)gd.getNextNumber();
            
            
            axialRange= (double)gd.getNextNumber();
            stepZloc=((double)gd.getNextNumber())/1000.;
            nwat=(double)gd.getNextNumber();
            zfocus=(double)gd.getNextNumber();
            
            
            photonNumber=(double)gd.getNextNumber();
            background=(double)gd.getNextNumber();
            
//            try{
//                String str_=gd.getNextString();
//                if (str_.length()>0){
//                    maxValuePlot=(double)Double.parseDouble(str_);
//                }
//                else{
//                    maxValuePlot=-1;
//                }
//            }catch(Exception ee){maxValuePlot=-1;}
//            pathResult = gd.getNextString(); 
            

        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        if (pathResult.equals(path_calibration)){
            IJ.log("Process stopped: the result crlb path has to be different than calibration path");
            processing=new Boolean(false);
            return;
        }
        
//        if (pathResult.equals(path_calibration2)){
//            if (pathResult.length()>0){//only if path not null
//                IJ.log("Process stopped: the result crlb path has to be different than calibration path");
//                processing=new Boolean(false);
//                return;
//            }
//        }
        
        
        
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        if (GPU_computation){
            MyCudaStream.init(1);
            org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
            if (!dp.loading){
                IJ.log("impossible to load "+path_calibration);
                processing=new Boolean(false);
                return;
            }

            if (dp.param!=null){


                //dp.
                dp.setSizeoutput(sizePatchPhaseRet);
                dp.setNwat(nwat);
                dp.param.Zfocus=zfocus;
                org.pasteur.imagej.process.gpu.CRLB_ crlb;
    //            if (dp2!=null){
    //                dp2.setSizeoutput(sizePatchPhaseRet);
    //                crlb =new CRLB(dp,dp2,-axialRange/2.,axialRange/2.,stepZloc,photonNumber,background,maxValuePlot); 
    //            }
    //            else{
                    crlb =new org.pasteur.imagej.process.gpu.CRLB_(dp,axialRange,stepZloc,photonNumber,background,-1); 
    //            }


                crlb.run(null);

            }

            MyCudaStream.destroy();
        }
        else{
            org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
            if (!dp.loading){
                IJ.log("impossible to load "+path_calibration);
                processing=new Boolean(false);
                return;
            }

            if (dp.param!=null){


                //dp.
                dp.setSizeoutput(sizePatchPhaseRet);
                dp.setNwat(nwat);
                dp.param.Zfocus=zfocus;
                org.pasteur.imagej.process.cpu.CRLB crlb;
    //            if (dp2!=null){
    //                dp2.setSizeoutput(sizePatchPhaseRet);
    //                crlb =new CRLB(dp,dp2,-axialRange/2.,axialRange/2.,stepZloc,photonNumber,background,maxValuePlot); 
    //            }
    //            else{
                    crlb =new org.pasteur.imagej.process.cpu.CRLB(dp,axialRange,stepZloc,photonNumber,background,-1); 
    //            }


                crlb.run(null);

            }

        }
        
//        classes.cudaProcess.DataPhase dp2=null;
//        if (path_calibration2.length()>3){
//            dp2 = new classes.cudaProcess.DataPhase(sizeFFT,path_calibration2);
//        }
        

        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    void phaseoptim(){
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        double photonNumber=3000;
        double background=50;
        
        
        String pathResult="";
        
        
        if (false){
        }
        else{
            
            String path=IJ.getDirectory("image");
            if (path==null){
                path=IJ.getDirectory("current");
            }
            
            
            
            //String [] fieldZernike = {"pixel-based (prototype)","Zernike 15 coefs","Zernike 28 coefs","Zernike 45 coefs"};
            String [] fieldZernike = {"pixel-based (long)","Zernike 15 coefs","Zernike 28 coefs","Zernike 45 coefs"};
            int [] correspondingFieldZernike = {0,15,28,45};
            if (zernikeCoef>=fieldZernike.length){
                zernikeCoef=fieldZernike.length-1;
            }
            if (zernikeCoef<0){
                zernikeCoef=0;
            }
            
//            int nn=(path.lastIndexOf("."));
//            if ((nn>path.length()-6)&&(nn>=0)){
//                path=path.substring(0,nn);
//            }
            path=path+"phase_optimized.json";
            

            GenericDialog gd = new GenericDialog("Phase optimization!");

            gd.addMessage("This method is works only with GPU");
            
            gd.addMessage("Image parameters");
            
            gd.addNumericField("pixel size (µm): ", xystep, 3);
            
            
            gd.addNumericField("Numerical aperture: ", na, 3);
            gd.addNumericField("Immersion_refractive index: ", noil, 3,6,"~1.518 for oil, ~1.33 for water");
            gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
            
            gd.addNumericField("distance_focus_to_coverslip (µm): ", zfocus, 3);
            gd.addNumericField("wavelength: ", wavelength, 3);

            gd.addMessage("Optimization parameters");
            
            
            gd.addNumericField("axial_range (µm): ", axialRange,2);
            gd.addNumericField("z_step (µm): ", stepZloc,2);
            
            gd.addNumericField("Photon number: ", photonNumber,0);
            gd.addNumericField("Background intensity: ", background,0);
            gd.addCheckbox("Use_pupil_apodization factor", withApoFactor);
            gd.addChoice("Phase-based method", fieldZernike,  fieldZernike[zernikeCoef]);
            
            gd.addNumericField("Iteration #: ", iterationNumber, 0);

            
            
            gd.addStringField("Result optimized phase file:",pathResult,sizeTextString); 
        

            TextField textField2; 
            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;

                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

                textField2 = texts.get(t++); 
                MouseOptionSave mos=new MouseOptionSave(textField2,pathResult,"Export calibration file");
                textField2.addMouseListener(mos); 


            }
            
            
            gd.showDialog();
            
            int orderFit=1;
            
            if (!gd.wasCanceled()){
                
                
                xystep = (double)gd.getNextNumber();
                na = (double)gd.getNextNumber();
                noil = (double)gd.getNextNumber();
                nwat = (double)gd.getNextNumber();
                zfocus = (double)gd.getNextNumber();
                wavelength = (double)gd.getNextNumber();

                axialRange= (double)gd.getNextNumber();
                stepZloc=(double)gd.getNextNumber();

                photonNumber=(double)gd.getNextNumber();
                background=(double)gd.getNextNumber();

                
                withApoFactor=gd.getNextBoolean();
                
                zernikeCoef=gd.getNextChoiceIndex();
                iterationNumber = (int)gd.getNextNumber();
                
                pathResult = gd.getNextString(); 
            }
            else{
                processing=new Boolean(false);
                return;
            }
            IJ.log("zernikeCoef "+zernikeCoef);
            
            int sizeFFT=128;

            MyCudaStream.init(1);
            
            
            if (zernikeCoef==0){
                        
                org.pasteur.imagej.process.gpu.GenericPhaseOptimization_test po = new org.pasteur.imagej.process.gpu.GenericPhaseOptimization_test(sizeFFT,sizeFFT,xystep,axialRange, stepZloc,photonNumber, background,wavelength,noil,nwat,zfocus,na,pathResult,withApoFactor);
                po.run(iterationNumber);
                po.free();
            }
            else{
                org.pasteur.imagej.process.gpu.ZernikePhaseOptimization_ po = new org.pasteur.imagej.process.gpu.ZernikePhaseOptimization_(sizeFFT,sizeFFT,xystep,axialRange, stepZloc,photonNumber, background,wavelength,noil,nwat,zfocus,na,(correspondingFieldZernike[zernikeCoef]),pathResult,withApoFactor);
                po.run(iterationNumber);
                po.free();
            }
            
            
            
            
            
            
            MyCudaStream.destroy();

            
            
            
        }
        
        IJ.log("well done");
        processing=new Boolean(false);
        
        
        
    }
    
    
    
    
    
    
        
    void calibration(){
        
        
        
        ImagePlus imp = IJ.getImage();
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        
        if (imp==null){
            IJ.log("no image opened");
        }
        else{
            
            String pathdir=IJ.getDirectory("image");
            String path=pathdir;
            if (path==null){
                path=IJ.getDirectory("current");
            }
//            int nn=(path.lastIndexOf("."));
//            if ((nn>path.length()-6)&&(nn>=0)){
//                path=path.substring(0,nn);
//            }
            path=path+"ZOLA_calibration_PSF.json";
            
            boolean stacked=true;
            int nbImage=imp.getNSlices();
            if (nbImage==1){
                nbImage=imp.getNFrames();
                stacked=false;
            }
            if (nbImage==1){
                IJ.log("WARNING, you should provide a stack with more than 1 image");
            }
            ImageStack ims= imp.getStack();
            
            int width=imp.getWidth();
            int height=imp.getHeight();
            
            boolean scmosLoad=false;
            if (isSCMOS){
                scmosLoad=loadSCMOScameraFiles(width, height, null);
            }
                
            
            
            float [] xp;
            float [] yp;
            
            try{
                xp=imp.getRoi().getFloatPolygon().xpoints;
                yp=imp.getRoi().getFloatPolygon().ypoints ;
            }
            catch (Exception ee){
                IJ.error("Please, click on at least one bead using ImageJ point selection tool before doing calibration.");
                processing=new Boolean(false);
                return;
            }
            
            sigma=0.9;
            
            double [][][] image =new double [nbImage][width][height];
            
            for (int i=0;i<nbImage;i++){
                ImageProcessor ip;
                if (stacked){
                    ip = ims.getProcessor(i+1);
                }
                else{
                    imp.setSlice(i+1);
                    ip = imp.getProcessor();
                }
                for (int ii=0,jj=0;ii<width;ii++,jj++){
                    for (int iii=0,jjj=0;iii<height;iii++,jjj++){
                        image[i][ii][iii]=ip.getPixelValue(jj, jjj);
                    }
                }
            }
            
            
            //classes.ImageShow.imshow(image,"image input");
            //String [] fieldZernike = {"6","10","15","21","28","36","45"};
            //String [] fieldZernike = {"pixel-based (prototype)","Zernike 15 coefs","Zernike 28 coefs","Zernike 45 coefs","Zernike 66 coefs"};
            String [] fieldZernike = {"double-helix (long)","Zernike 15 coefs","Zernike 28 coefs","Zernike 45 coefs","Zernike 66 coefs"};
            int [] correspondingFieldZernike = {0,15,28,45,66};
            if (zernikeCoef>=fieldZernike.length){
                zernikeCoef=fieldZernike.length-1;
            }
            if (zernikeCoef<0){
                zernikeCoef=0;
            }
            String [] fieldSide = {"far -> close to objective","close -> far from objective"};
            if (axialside>=fieldSide.length){
                axialside=fieldSide.length-1;
            }
            if (axialside<0){
                axialside=0;
            }
            GenericDialog gd = new GenericDialog("ZOLA: Calibration - PSF modeling");

            
            
                
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            gd.addCheckbox("Run_on_GPU", GPU_computation);
            
            
                if (!isSCMOS){
                    gd.addMessage("Camera is EMCCD",fontBold);
                    gd.addMessage("ADU = "+adu);
                    gd.addNumericField("Gain = ", gain, 1,6,"(1 if camera provides photon counts)");
                    gd.addMessage("Offset = "+offset);
                    
                }
                else{
                    
                    if (scmosLoad){
                        gd.addMessage("Camera is SCMOS",fontBold);
                    }
                    else{
                        gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                    }
//                    String sep=File.separator;
//                    String [] var=this.path_SCMOSvariance.split(sep);
//                    String [] off=this.path_SCMOSoffset.split(sep);
//                    String [] gain=this.path_SCMOSgain.split(sep);
//                    gd.addMessage("Offset file = "+off[var.length-1]);
//                    gd.addMessage("Variance file = "+var[var.length-1]);
//                    gd.addMessage("Gain file = "+gain[var.length-1]);
                    
                    gd.addMessage("Offset file = "+path_SCMOSoffset);
                    gd.addMessage("Variance file = "+path_SCMOSvariance);
                    gd.addMessage("Gain file = "+path_SCMOSgain);
                }
                
            
            gd.addNumericField("Pixel_size:", xystep*1000, 1,6,"(nm)");
            
            
            
            gd.addMessage("Sample parameters",fontBold);
            
            gd.addNumericField("Z_step: ", zstep*1000, 1,6,"(nm)");
            gd.addChoice("Bead_moving direction", fieldSide,  fieldSide[axialside]);
            
            gd.addNumericField("Numerical_aperture: ", na, 3);
            gd.addNumericField("immersion_refractive index: ", noil, 3,6,"~1.518 for oil, ~1.33 for water");
            gd.addNumericField("wavelength: ", wavelength, 3,6,"(µm)");

            gd.addMessage("Advanced parameters",fontBold);
            
            
            gd.addNumericField("Patch_size: ", sizePatchPhaseRet, 0,6,"(pixels)");
            //gd.addSlider("Zernike_coefficient #", 0, 2*((int)Math.sqrt(nbImage)), (int)Math.sqrt(nbImage));
            gd.addCheckbox("Use_pupil_apodization factor", withApoFactor);
            gd.addChoice("Zernike_coefficient number", fieldZernike,  fieldZernike[zernikeCoef]);
            //gd.addNumericField("Sigma_gaussian_blur", sigma,  2,6,"(pixels)");
            //gd.addNumericField("Zernike_coefficient #: ", zernikeCoef, 0);
            gd.addNumericField("Iteration number: ", iterationNumber, 0);

            gd.addStringField("Result_calibration_file:",path,sizeTextString); 
            
        TextField textField2; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            
            
            textField2 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField2,pathdir,"Export calibration file");
            textField2.addMouseListener(mos); 
            
            
        }
            
            gd.showDialog();
            
            
            if (!gd.wasCanceled()){
                
                GPU_computation = (boolean)gd.getNextBoolean();
                if (!isSCMOS){
                    gain = (double)gd.getNextNumber();
                }
                xystep = (double)gd.getNextNumber()/1000;
                zstep = (double)gd.getNextNumber()/1000;
                
                axialside=gd.getNextChoiceIndex();
                
                na = (double)gd.getNextNumber();
                noil = (double)gd.getNextNumber();
                wavelength = (double)gd.getNextNumber();
                
                sizePatchPhaseRet = (int)gd.getNextNumber();
                withApoFactor=(boolean)gd.getNextBoolean();
                zernikeCoef=gd.getNextChoiceIndex();
                
                //sigma = (double)gd.getNextNumber();
                iterationNumber = (int)gd.getNextNumber();
                
                
                path_calibration = gd.getNextString(); 
            }
            else{
                processing=new Boolean(false);
                return;
            }
            
            int nnn=(path_calibration.lastIndexOf("."));
            //IJ.log("path_calibration "+path_calibration+"  "+nnn+"  "+path_calibration.length());
            if ((nnn>=0)&&(nnn>path_calibration.length()-6)){
                path_calibration=path_calibration.substring(0,nnn);
            }
            
            path_calibration+=".json";
            //IJ.log("path_calibration "+path_calibration);
            if (sizePatchPhaseRet%2==1){
                sizePatchPhaseRet++;
            }
            
            int sizeFFT=180;
            sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
            sizeFFT=Math.min(sizeFFT,maxFFT);
            
            
            
            
            if (gain==0||adu==0){
                new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
            }
            else{
                IJ.log("ZOLA PSF modeling started");
                
                InitBackgroundAndPhotonNumber paramImage ;
                if (isSCMOS){
                    paramImage = new InitBackgroundAndPhotonNumber(image,xp,yp,sizePatchPhaseRet,this.scmosvariance,this.scmosgain,this.scmosoffset,this.scmosvargain,zstep);

                }
                else{
                    paramImage = new InitBackgroundAndPhotonNumber(image,xp,yp,sizePatchPhaseRet,adu,gain,offset,zstep);

                }
                long t=System.currentTimeMillis();
                
                if (GPU_computation){
                    //maybe test if image is a square
                    //IJ.log("init...");
                    MyCudaStream.init(1);
                    //IJ.log("init cuda ok "+xp+"  "+yp+"  "+sizePatchPhaseRet+"  "+rescale_slope+"  "+rescale_intercept+"  "+zstep+"  "+image);
                    
                    //IJ.log("param ok");
                    if (zernikeCoef==0){
                        
                        org.pasteur.imagej.process.gpu.GenericPhaseRetrieval_ gp = new org.pasteur.imagej.process.gpu.GenericPhaseRetrieval_(sizeFFT,xystep,zstep,wavelength,noil,na,paramImage,path_calibration,sigma,axialside,withApoFactor);
                        gp.run(iterationNumber);
                        gp.free();
                    }
                    else{
                        org.pasteur.imagej.process.gpu.ZernikePhaseRetrieval_ prp = new org.pasteur.imagej.process.gpu.ZernikePhaseRetrieval_(sizeFFT,xystep,zstep,wavelength,noil,na,(correspondingFieldZernike[zernikeCoef]),paramImage,path_calibration,sigma,axialside,withApoFactor);
                        //IJ.log("process many ok");
                        prp.run(iterationNumber);
                        //IJ.log("run ok");
                        prp.free();
                    }
                    MyCudaStream.destroy();
                }
                else{
                    
                    if (zernikeCoef==0){
                        
                        org.pasteur.imagej.process.cpu.GenericPhaseRetrieval prp = new org.pasteur.imagej.process.cpu.GenericPhaseRetrieval(sizeFFT,xystep,zstep,wavelength,noil,na,paramImage,path_calibration,sigma,axialside,withApoFactor);
                        prp.run(iterationNumber);
                    }
                    else{
                        org.pasteur.imagej.process.cpu.ZernikePhaseRetrieval prp = new org.pasteur.imagej.process.cpu.ZernikePhaseRetrieval(sizeFFT,xystep,zstep,wavelength,noil,na,correspondingFieldZernike[zernikeCoef],paramImage,path_calibration,sigma,axialside,withApoFactor);
                        prp.run(iterationNumber);
                    }
                    
                    
                }
                long t2=System.currentTimeMillis();
                IJ.log("elapsed time "+Math.floor((double)(t2-t)/(double)(1000*60))+" min "+Math.floor(60*((((double)(t2-t)/(double)(1000*60)))-Math.floor((double)(t2-t)/(double)(1000*60))))+" sec");
            }

            
            
        }
        
        processing=new Boolean(false);
        
        
        
        
        
    }
    
    
    
    
    
    void wobblecorrection(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        if (sl!=null){
            
            String path=path_localization;
            int lllll=path.lastIndexOf(File.separator);
            path=path.substring(0,lllll);
            
            
            
            
            GenericDialog gd = new GenericDialog("Wobble correction");
            
            
            
            gd.addStringField("Wobble calibration file: ",path_wobble,sizeTextString); 

            
            TextField textField;
            
            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;

                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

                textField = texts.get(t++); 
                MouseOptionLoad mol=new MouseOptionLoad(textField,path_wobble,"Import calibration file wobble");
                textField.addMouseListener(mol); 
                
                


            }
            
            gd.showDialog();
            
            
            if (!gd.wasCanceled()){

                
                path_wobble  = gd.getNextString();
                
                
                sl=WobbleCorrection.correction(path_wobble,sl);
                
                ZRendering.colorRendering(null,sl, 20,1,true);
                IJ.log("wobble effect corrected. Please, export the localization table to save the correction");
                
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
        
    void wobblecalibration(){
        
        
        nbStream=5;
        nbThread=100;  
            
        ImagePlus imp=IJ.getImage();


        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
            processing=new Boolean(true);
            
            
            if (imp==null){
                IJ.log("no image opened");
                processing=new Boolean(false);
            }
            /*else if ((imp.getWidth()>180)||(imp.getHeight()>180)){
                IJ.log("image should be < 180 pixels for the moment");
            }*/
            else{
                
                
                Roi r=imp.getRoi();
                int width=imp.getWidth();
                int height=imp.getHeight();
                
                boolean scmosLoad=false;
                if (isSCMOS){
                    scmosLoad=loadSCMOScameraFiles(width, height, r);
                }
                
                
                
                String pathdir=IJ.getDirectory("image");
                if (pathdir==null){
                    pathdir=IJ.getDirectory("current");
                }
                String path=pathdir;
                
                path=path+"ZOLA_wobble_calibration.csv";
                
                GenericDialog gd=null;
                
                //since NonBlockingGenericDialog incompatible with headless mode:
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())){
                    gd = new NonBlockingGenericDialog("ZOLA: wobble correction");
                }
                else{
                    gd = new GenericDialog("ZOLA: wobble correction");
                }
                
                
                
                Font font = gd.getFont();
                Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
                
                gd.addCheckbox("Run_on_GPU", GPU_computation);
                
                
                
                
                
                
                if (!isSCMOS){
                    gd.addMessage("Camera is EMCCD",fontBold);
                    gd.addMessage("ADU = "+adu);
                    gd.addNumericField("Gain = ", gain, 1,6,"(1 if camera provides photon counts)");
                    gd.addMessage("Offset = "+offset);
                    
                }
                else{
                    if (scmosLoad){
                        gd.addMessage("Camera is SCMOS",fontBold);
                    }
                    else{
                        gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                    }
//                    String sep=File.separator;
//                    String [] var=this.path_SCMOSvariance.split(sep);
//                    String [] off=this.path_SCMOSoffset.split(sep);
//                    String [] gain=this.path_SCMOSgain.split(sep);
//                    gd.addMessage("Offset file = "+off[var.length-1]);
//                    gd.addMessage("Variance file = "+var[var.length-1]);
//                    gd.addMessage("Gain file = "+gain[var.length-1]);
                    
                    gd.addMessage("Offset file = "+path_SCMOSoffset);
                    gd.addMessage("Variance file = "+path_SCMOSvariance);
                    gd.addMessage("Gain file = "+path_SCMOSgain);
                    
                }
                
                
                
                gd.addStringField("PSF_calibration_file:",path_calibration,sizeTextString); 
                
                
                
                gd.addMessage("Sample parameters", fontBold);
                
                gd.addNumericField("Z_step: ", zstep*1000, 1,6,"(nm)");
                
                
                
                
                
                gd.addMessage("Advanced parameter", fontBold);
                
                
                gd.addNumericField("Patch_size: ", sizePatch,0,6,"(pixels)");
                
                //gd.addNumericField("number_of_streams: ", nbStream,0);
                //gd.addNumericField("number_of_threads: ", nbThread,0);
                
                
                
                
                
                gd.addNumericField("Expected_axial_range: ", axialRange,2,6,"(µm)");
                
                gd.addNumericField("Min_number_of_photons: ", photonThreshold,0);
                
                
                gd.addStringField("wobble_calibration file:",path,sizeTextString); 
                
                Choice choicefield=null;
                
                
                
                TextField textField; 
                TextField textField2; 
                // Add a mouse listener to the config file field 
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
                {
                    int t=0;
                    //textField
                    Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

                    textField = texts.get(t++); 
                    MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
                    textField.addMouseListener(mol); 

                    textField2 = texts.get(t++); 
                    MouseOptionSave mos=new MouseOptionSave(textField2,pathdir,"Export wobble_calibration file");
                    textField2.addMouseListener(mos); 
                    
                    
                    

                }
        
        
                gd.showDialog();
                
                
                
                
                
                if (!gd.wasCanceled()){
                    GPU_computation = (boolean)gd.getNextBoolean();
                    
                    if (!isSCMOS){
                        gain = (double)gd.getNextNumber();
                    }
                    
                    path_calibration = gd.getNextString(); 
                    
                    zstep = (double)gd.getNextNumber()/1000;
                    
                    sizePatch = (int)gd.getNextNumber();
                    //nbStream = (int)gd.getNextNumber();
                    //nbThread = (int)gd.getNextNumber();
                    
                    axialRange= (double)gd.getNextNumber();
                    
                    
                            
                    
                    photonThreshold = (int)gd.getNextNumber();
                    
                    path = gd.getNextString(); 
                }
                else{
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                int nnn=(path.lastIndexOf("."));
                if ((nnn>=0)&&(nnn>path.length()-6)){
                    path=path.substring(0,nnn);
                }
                path+=".csv";
            
                if (path_localization.equals(path_calibration)){
                    IJ.log("Process stopped: the localization path has to be different than calibration path");
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                if (sizePatch%2==1){
                    sizePatch++;
                }
                
                
                
                
                
                
                
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                }
                if ((width<sizePatch)||(height<sizePatch)){
                    IJ.log("Process stopped: Localization impossible. The image, or selected ROI should be larger than patch size");
                    IJ.log("Please, select a large rectangle in the image before running the localization");
                    IJ.log("image width:"+width+" ; height:"+height+" ; sizePatch:"+sizePatch);
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                
                
                
                if (gain==0||adu==0){
                    new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
                }
                else{
                    
                    
                    
                    
                    int sizeFFT=128;
                    sizeFFT=Math.max(sizeFFT,sizePatch*2);
                    sizeFFT=Math.min(sizeFFT,maxFFT);
                    
                    
                    if (GPU_computation){
                        MyCudaStream.init(nbStream+1);
                        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(dp.param.noil);
                        dp.param.Zfocus=0;
                        IJ.log("ref index = "+dp.param.noil+"    "+dp.param.nwat+"  "+noil+"   na"+dp.param.na);
                        xystep=dp.param.xystep;
                        IJ.log("Localization started");


                        org.pasteur.imagej.process.gpu.LocalizationPipeline_ lp;
                        if (isSCMOS){
                            lp=new org.pasteur.imagej.process.gpu.LocalizationPipeline_(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain,false);
                        }
                        else{
                            lp=new org.pasteur.imagej.process.gpu.LocalizationPipeline_(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,adu, gain,offset,false);
                        }
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);


                        /*classes.testClass lp =new classes.testClass(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,rescale_slope, rescale_intercept);
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);*/



                        dp.free();
                    
                    
                        MyCudaStream.destroy();
                    }
                    else{
                        org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(dp.param.noil);
                        
                        dp.param.Zfocus=0;
                        xystep=dp.param.xystep;
                        IJ.log("Localization started");
                        
                        org.pasteur.imagej.process.cpu.LocalizationPipeline lp;
                        if (isSCMOS){
                            
                            lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain,false);
                        }
                        else{
                            lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,adu, gain,offset,false);
                        }
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);
                        
                        
                        
                    }
                    
                    
                    
                    WobbleCorrection.calibration(path, sl, zstep*1000, xystep*1000);
                    
                    this.path_wobble=path;
                    
                }
            }
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    void setCameraEMCCD(){
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            
            return;
        }
        
        processing=new Boolean(true);


        GenericDialog gd = new GenericDialog("ZOLA: EMCCD Camera parameters");

        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}

        gd.addNumericField("EMCCD_Camera_ADU: ", adu, 4,6,"(1 if camera provides photon counts)");
        gd.addNumericField("EMCCD_Camera_Gain: ", gain, 4,6,"(1 if camera provides photon counts)");
        gd.addNumericField("EMCCD_Camera_Offset: ", offset, 4,6,"(0 if camera provides photon counts)");



        gd.showDialog();
        
        

        if (!gd.wasCanceled()){


            adu=(double)gd.getNextNumber();
            gain=(double)gd.getNextNumber();
            offset=(double)gd.getNextNumber();
            this.isSCMOS=false;
        }
        else{
            processing=new Boolean(false);
            return;
        }




        
        processing=new Boolean(false);
        
        
        
        
        
    }
    
    
    
    
    
    
    
    
    
    void setCameraSCMOS(){
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            
            return;
        }
        
        processing=new Boolean(true);


        GenericDialog gd = new GenericDialog("ZOLA: SCMOS Camera parameters");

        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}    
        String path_SCMOSoffset=this.path_SCMOSoffset;
        String path_SCMOSvariance=this.path_SCMOSvariance;
        String path_SCMOSgain=this.path_SCMOSgain;

        gd.addMessage("SCMOS Camera have pixel-based readout noise",fontBold);

        gd.addMessage("");
        gd.addStringField("TIF image path of offset:",path_SCMOSoffset,sizeTextString); 
        gd.addStringField("TIF image path of variance:",path_SCMOSvariance,sizeTextString); 
        gd.addStringField("TIF image path of gain:",path_SCMOSgain,sizeTextString); 



        TextField textField1,textField2,textField3; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            
            
            textField1 = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField1,path_SCMOSoffset,"Import SCMOS offset calibration image");
            textField1.addMouseListener(mol); 
            
            textField2 = texts.get(t++); 
            mol=new MouseOptionLoad(textField2,path_SCMOSvariance,"Import SCMOS variance calibration image");
            textField2.addMouseListener(mol); 
            
            textField3 = texts.get(t++); 
            mol=new MouseOptionLoad(textField3,path_SCMOSgain,"Import SCMOS gain calibration image");
            textField3.addMouseListener(mol); 
            
            
        }
        
        
        gd.showDialog();
        
        if (!gd.wasCanceled()){


            path_SCMOSoffset = gd.getNextString(); 
            path_SCMOSvariance = gd.getNextString(); 
            path_SCMOSgain = gd.getNextString(); 
            boolean isOK=false;
            if ((path_SCMOSoffset.length()>3)&&(path_SCMOSvariance.length()>3)&&(path_SCMOSgain.length()>3)){
                try{
                    ImagePlus imp1=new ImagePlus(path_SCMOSoffset);
                    ImagePlus imp2=new ImagePlus(path_SCMOSvariance);
                    ImagePlus imp3=new ImagePlus(path_SCMOSgain);
                    int w1=imp1.getWidth();
                    int w2=imp2.getWidth();
                    int w3=imp3.getWidth();
                    int h1=imp1.getHeight();
                    int h2=imp2.getHeight();
                    int h3=imp3.getHeight();
                    if ((w1>0)&&(h1>0)&&(w2>0)&&(h2>0)&&(w2>0)&&(h2>0)){
                        if ((w1!=w2)||(w1!=w3)||(w2!=w3)||(h1!=h2)||(h1!=h3)||(h2!=h3)){
                            IJ.log("ERROR:images should have the same size");
                            IJ.log("widths "+w1+"  "+w2+"  "+w3);
                            IJ.log("heights "+h1+"  "+h2+"  "+h3);
                        }
                        else{
                            isOK=true;
                        }
                    }
                    else{
                        IJ.log("ERROR: something wrong with one of the image: size=0");
                    }

                }catch(Exception e){
                    IJ.log("ERROR occured while oppening SCMOS images. Please, set again the 3 image paths to configure SCMOS camera noise "+e);
                }
            }
            else{
                isOK=false;
                IJ.log("Please, select variance, gain and offset files to set SCMOS camera parameters ");
            }
            
            
            if (isOK){
                this.path_SCMOSoffset=path_SCMOSoffset;
                this.path_SCMOSvariance=path_SCMOSvariance;
                this.path_SCMOSgain=path_SCMOSgain;
                this.isSCMOS=true;
                IJ.log("SCMOS camera parameters are set.");
            }
            else{
                IJ.log("ERROR: something wrong happens");
            }
            
            
            
            
        }

        


        
        processing=new Boolean(false);
        
        
        
        
        
    }
    
    
    
    
    
    
    
        
    
    
    
    
    void resetDriftCorrection(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        for (int i=0;i<sl.fl.size();i++){
            for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                sl.fl.get(i).loc.get(ii).removeDrift();
            }
        }
        
        ZRendering.colorRendering(null,sl, 20,0,true);
        processing=new Boolean(false);
        IJ.log("Drift correction removed");
    }
    
    
    
    
    
    
    void setDriftCorrection(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        String path=path_localization;
        int lllll=path.lastIndexOf(File.separator);
        path=path.substring(0,lllll);

        GenericDialog gd = new GenericDialog("ZOLA: 3D Drift Correction");

        Font font = gd.getFont();
        Font fontBold= gd.getFont();
        
        gd.addStringField("4 columns drift file: [frame,Xdrift,Ydrift,Zdrift]: ","",sizeTextString);


        TextField textField;
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;

            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path,"Import drift table");
            textField.addMouseListener(mol); 



        }


        gd.showDialog();

        if (!gd.wasCanceled()){
            path = gd.getNextString();
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        processing=new Boolean(true);
        
        double [][] drift=FileVectorLoader.getTableFromFile(path, ",");
        Arrays.sort(drift, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o2[0]).compareTo(o1[0]);
            }
        });
        int maxFrame=0;
        for (int i=0;i<drift.length;i++){
            if (maxFrame<(int)drift[i][0]){
                maxFrame=(int)drift[i][0];
            }
        }
        IJ.log("max frame "+maxFrame);
        int [] link = new int[maxFrame+1];
        for (int i=0;i<maxFrame;i++){
            link[i]=-1;
        }
        for (int i=0;i<drift.length;i++){
            link[(int)drift[i][0]]=i;
        }
        
        StackLocalization slb = new StackLocalization();
        
        
        for (int i=0;i<sl.fl.size();i++){
            int frame=sl.fl.get(i).numFrame;
            if (frame<maxFrame){
                if (link[frame]>=0){
                    //drift found in drift file for current frame
                    FrameLocalization flb=new FrameLocalization(sl.fl.get(i).numFrame);
                    for (int j=0;j<sl.fl.get(i).loc.size();j++){
                        PLocalization p = sl.fl.get(i).loc.get(j);
                        PLocalization pb= new PLocalization();
                        pb=p.copy();
                        pb.setDrift(drift[link[frame]][1], drift[link[frame]][2], drift[link[frame]][3]);
                        flb.loc.add(pb);
                    }
                    slb.fl.add(flb);
                }
            }
            
        }
        sl=slb;
        ZRendering.colorRendering(null,sl, 20,0,true);
        processing=new Boolean(false);
        IJ.log("Drift correction removed");
    }
    
    
    
    
    
    
    void driftCorrection(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        if (sl!=null){
            
            String path=path_localization;
            int lllll=path.lastIndexOf(File.separator);
            path=path.substring(0,lllll);
            
            GenericDialog gd = new GenericDialog("ZOLA: 3D Drift Correction");

            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            gd.addCheckbox("Run_on_GPU", GPU_computation);
            
            String [] binNumber = {"2","3","4","5","6","7","8","9","10"};
            if (drift_bin>=binNumber.length){
                drift_bin=binNumber.length-1;
            }
            if (drift_bin<0){
                drift_bin=0;
            }
            
            
            gd.addNumericField("Cross-correlation_pixel_size: ", drift_pixelsizeNM,1,6,"(nm)");

            gd.addChoice("Number of time points:", binNumber,  binNumber[drift_bin]);
            
            gd.addNumericField("Maximum_drift: ", drift_sub_image_sizeUm,0,6,"(µm)");
            
            gd.addMessage("If non zero initial drift -> attach a localization table previously reconstructed and drift corrected:");
            gd.addStringField("Localization_table_attached [optional]: ","",sizeTextString);

            
            TextField textField;
            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;

                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

                textField = texts.get(t++); 
                MouseOptionLoad mol=new MouseOptionLoad(textField,path,"Import localization table");
                textField.addMouseListener(mol); 
                


            }
        

            gd.showDialog();

            if (!gd.wasCanceled()){
                GPU_computation = (boolean)gd.getNextBoolean();
                drift_pixelsizeNM = (double)gd.getNextNumber();
                drift_bin=gd.getNextChoiceIndex();
                drift_sub_image_sizeUm = (double)gd.getNextNumber();
                path = gd.getNextString();
            }
            else{
                processing=new Boolean(false);
                return;
            }
            
            
            if (GPU_computation){
                MyCudaStream.init(1);
                new org.pasteur.imagej.process.gpu.DriftCorrection_(sl,drift_pixelsizeNM, Integer.parseInt(binNumber[drift_bin]),drift_sub_image_sizeUm,path);
                MyCudaStream.destroy();
            }
            else{
                new org.pasteur.imagej.process.cpu.DriftCorrection(sl,drift_pixelsizeNM, Integer.parseInt(binNumber[drift_bin]),drift_sub_image_sizeUm);
            }
            ZRendering.colorRendering(null,sl, 20,0,true);
            
            IJ.log("Please, export the localization table to save drift correction");
            
            processing=new Boolean(false);
        }
        else{
            IJ.log("Please, import first a localization table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    void driftCorrectionFiducialMarker(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        if (sl!=null){
            
            String path=path_localization;
            int lllll=path.lastIndexOf(File.separator);
            path=path.substring(0,lllll);
            
            
            
            
            GenericDialog gd = new GenericDialog("Drift correction Fiducial");
            
            
            gd.addMessage("Conditions to merge consecutives beads");
            
            gd.addNumericField("max._X/Y distance between consecutive localizations: ", maxDistanceMergingPlan,1,6,"(nm)");
            gd.addNumericField("max._Z distance between consecutive localizations: ", maxDistanceMergingZ,1,6,"(nm)");
            
            gd.addNumericField("On_time: ", concecutiveFrameThreshold,0,6,"(~100 by default)");
            gd.addNumericField("Off_time: ",offFrame,0,6,"(~30 by default)");
            gd.addNumericField("frame_smoothing_factor: ", smoothingFrameNumber,0,6,"(~50 by default)");
            
            
            gd.addMessage("if non zero initial drift -> attach a localization table previously reconstructed and drift corrected:");
            gd.addStringField("localization_table_attached [optional]: ","",sizeTextString);
            
            
            gd.addStringField("Save_drift_table [optional]: ","",sizeTextString);

            
            TextField textField;
            TextField textField2;
            String path2=path;
            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;

                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 

                textField = texts.get(t++); 
                MouseOptionLoad mol=new MouseOptionLoad(textField,path,"Import localization table");
                textField.addMouseListener(mol); 
                
                textField2 = texts.get(t++); 
                MouseOptionSave mos=new MouseOptionSave(textField2,path2,"Export drift table");
                textField2.addMouseListener(mos); 
                
                


            }
            
            gd.showDialog();
            
            
            if (!gd.wasCanceled()){


                maxDistanceMergingPlan=gd.getNextNumber();

                maxDistanceMergingZ=gd.getNextNumber();
                
                concecutiveFrameThreshold=(int)gd.getNextNumber();
                offFrame = (int)gd.getNextNumber();
                smoothingFrameNumber=(int)gd.getNextNumber();
                
                path = gd.getNextString();
                String path_DriftTable = gd.getNextString();
                
                
                
                DriftCorrectionFiducial dcb=new DriftCorrectionFiducial(sl,maxDistanceMergingPlan, maxDistanceMergingZ,concecutiveFrameThreshold,50,smoothingFrameNumber,path,path_DriftTable,offFrame);
                sl=dcb.run();
                
                ZRendering.colorRendering(null,sl, 20,0,true);

                
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    void condensation(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        if (sl!=null){
            
            int offFrame=0;
            GenericDialog gd = new GenericDialog("ZOLA: Merging");

            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            gd.addNumericField("Max._lateral distance between consecutive localizations: ", maxDistanceMergingPlan,1,6,"(nm)");
            gd.addNumericField("Max._axial distance between consecutive localizations: ", maxDistanceMergingZ,1,6,"(nm)");
            
            gd.addNumericField("Off_frame: ", offFrame,0);


            

            gd.showDialog();
            
            
            if (!gd.wasCanceled()){


                maxDistanceMergingPlan=gd.getNextNumber();

                maxDistanceMergingZ=gd.getNextNumber();
                
                offFrame=(int)gd.getNextNumber();
            
            
                Condensation cond=new Condensation(sl,maxDistanceMergingPlan, maxDistanceMergingZ,offFrame);
                sl=cond.run();
                
                ZRendering.colorRendering(null,sl, 20,0,true);

                
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    void statAttachTime1(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            double pixSize=20;
            ImagePlus imp=ZRendering.scatter2D(sl, pixSize,true);
            IJ.run("Brightness/Contrast...");
            double minrange=imp.getDisplayRangeMin();
            double maxrange=imp.getDisplayRangeMax();
            
            int xinit=0;
            int yinit=0;
            int width=imp.getWidth();
            int height=imp.getHeight();
            new WaitForUserDialog("ROI selection", "Select a rectangle in the image to control X/Y range, and adjust brightness/Cont. to control Z range, Then press OK.").show();
            if (imp.isVisible()){
                Roi r=imp.getRoi();
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    xinit=bound.x;
                    yinit=bound.y;
                    width=bound.width;
                    height=bound.height;
                }
                minrange=imp.getDisplayRangeMin();
                maxrange=imp.getDisplayRangeMax();
            }
            
            
            
            double minX=Double.POSITIVE_INFINITY,maxX=Double.NEGATIVE_INFINITY,minY=Double.POSITIVE_INFINITY,maxY=Double.NEGATIVE_INFINITY,minZ=Double.POSITIVE_INFINITY,maxZ=Double.NEGATIVE_INFINITY;
            
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    PLocalization p = sl.fl.get(i).loc.get(j);
                    
                    
                    if (p.X<minX){
                        minX=p.X;
                    }
                    if (p.Y<minY){
                        minY=p.Y;
                    }
                    if (p.Z<minZ){
                        minZ=p.Z;
                    }
                    if (p.X>maxX){
                        maxX=p.X;
                    }
                    if (p.Y>maxY){
                        maxY=p.Y;
                    }
                    if (p.Z>maxZ){
                        maxZ=p.Z;
                    }
                    
                    
                    
                }
            }
            
            double minX2=xinit*pixSize+minX;
            double maxX2=xinit*pixSize+width*pixSize+minX;
            double minY2=yinit*pixSize+minY;
            double maxY2=yinit*pixSize+height*pixSize+minY;
            
            int offFrame=0;
            int start=0;
            int end=100;
            GenericDialog gd = new GenericDialog("ZOLA: stat residence time curve");

            gd.addNumericField("max._lateral distance between consecutive localizations: ", maxDistanceMergingPlan,1,6,"(nm)");
            gd.addNumericField("max._axial distance between consecutive localizations: ", maxDistanceMergingZ,1,6,"(nm)");
            
            gd.addNumericField("off_frame: ", offFrame,0);
            gd.addNumericField("start_histogram: ", start,0);
            gd.addNumericField("end_histogram: ", end,0);

            

            gd.showDialog();
            
            
            if (!gd.wasCanceled()){


                /*minX=gd.getNextNumber();
                maxX=gd.getNextNumber();
                
                minY=gd.getNextNumber();
                maxY=gd.getNextNumber();
                
                minZ=gd.getNextNumber();
                maxZ=gd.getNextNumber();*/
                
                
                StackLocalization slb = new StackLocalization();
                int nbtot=0;
                int nbafter=0;
                for (int i=0;i<sl.fl.size();i++){
                    int k=0;
                    FrameLocalization flb=new FrameLocalization(sl.fl.get(i).numFrame);
                    for (int j=0;j<sl.fl.get(i).loc.size();j++){
                        PLocalization p = sl.fl.get(i).loc.get(j);
                        if ((p.X>=minX2)&&(p.X<=maxX2)&&(p.Y>=minY2)&&(p.Y<=maxY2)&&(p.Z>=minrange+minZ)&&(p.Z<=maxrange+minZ)){
                            PLocalization pb= new PLocalization();
                            pb=p.copy();
                            flb.loc.add(pb);
                            nbafter++;
                        }
                        else{
                            nbtot++;
                        }
                    }
                    if (flb.loc.size()>0){
                        slb.fl.add(flb);
                    }
                }
                nbtot+=nbafter;
                
                IJ.log("Localization number after filtering: "+nbafter+"/"+nbtot);
                ZRendering.colorRendering(null,slb, 20,0,true);
                
                maxDistanceMergingPlan=gd.getNextNumber();

                maxDistanceMergingZ=gd.getNextNumber();
                
                offFrame=(int)gd.getNextNumber();
                start=(int)gd.getNextNumber();
                end=(int)gd.getNextNumber();
                
                
                IJ.log("");
                IJ.log("**********ON TIME***********");
                Condensation cond=new Condensation(slb,maxDistanceMergingPlan, maxDistanceMergingZ,offFrame);
                StackLocalization smerged=cond.run();
                
                
                double [] p=cond.fitExponential1(start,end);
                
                parameterFit1=p[0];
                parameterFit2=p[1];
                IJ.log("Mean residence time="+p[1]+" frames");
                IJ.log("");
                
                
                
                
                //IJ.log("**********OFF TIME***********");
                //IJ.log("note that here we use condensation with off_time=infinity");
                //classes.PasteurLocalization.Condensation cond2=new classes.PasteurLocalization.Condensation(slb,maxDistanceMergingPlan, maxDistanceMergingZ,Integer.MAX_VALUE);
                //StackLocalization smerged2=cond2.run();
                //cond.computeOffTimes(smerged2);
                //cond.fitExponential1(start);
                //IJ.log("");
                
                
                
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        
        
        
        processing=new Boolean(false);
    }
    
    
    
    
    void filtering(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            double pixSize=20;
            ImagePlus imp=ZRendering.scatter2D(sl, pixSize,true);
            IJ.run("Brightness/Contrast...");
            double minrange=imp.getDisplayRangeMin();
            double maxrange=imp.getDisplayRangeMax();
            
            int xinit=0;
            int yinit=0;
            int width=imp.getWidth();
            int height=imp.getHeight();
            new WaitForUserDialog("ROI selection", "Select a rectangle in the image to control X/Y range, and adjust brightness/Cont. to control Z range, Then press OK.").show();
            if (imp.isVisible()){
                Roi r=imp.getRoi();
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    xinit=bound.x;
                    yinit=bound.y;
                    width=bound.width;
                    height=bound.height;
                }
                minrange=imp.getDisplayRangeMin();
                maxrange=imp.getDisplayRangeMax();
            }
            
            
            int offFrame=0;
            GenericDialog gd = new GenericDialog("ZOLA: Filtering");
            double minX=Double.POSITIVE_INFINITY,maxX=Double.NEGATIVE_INFINITY,minY=Double.POSITIVE_INFINITY,maxY=Double.NEGATIVE_INFINITY,minZ=Double.POSITIVE_INFINITY,maxZ=Double.NEGATIVE_INFINITY;
            double maxCRLBX=Double.NEGATIVE_INFINITY,maxCRLBY=Double.NEGATIVE_INFINITY,maxCRLBZ=Double.NEGATIVE_INFINITY;
            double maxScore=Double.NEGATIVE_INFINITY;
            
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    PLocalization p = sl.fl.get(i).loc.get(j);
                    
                    //if (p.score<minScore){
                    //    minScore=p.score;
                    //}
                    if (p.score>maxScore){
                        maxScore=p.score;
                    }
                    if (p.X<minX){
                        minX=p.X;
                    }
                    if (p.Y<minY){
                        minY=p.Y;
                    }
                    if (p.Z<minZ){
                        minZ=p.Z;
                    }
                    if (p.X>maxX){
                        maxX=p.X;
                    }
                    if (p.Y>maxY){
                        maxY=p.Y;
                    }
                    if (p.Z>maxZ){
                        maxZ=p.Z;
                    }
                    
                    
                    if (p.crlb_X>maxCRLBX){
                        maxCRLBX=p.crlb_X;
                    }
                    if (p.crlb_Y>maxCRLBY){
                        maxCRLBY=p.crlb_Y;
                    }
                    if (p.crlb_Z>maxCRLBZ){
                        maxCRLBZ=p.crlb_Z;
                    }
                }
            }
            
            double minX2=xinit*pixSize+minX;
            double maxX2=xinit*pixSize+width*pixSize+minX;
            double minY2=yinit*pixSize+minY;
            double maxY2=yinit*pixSize+height*pixSize+minY;
            
            /*gd.addNumericField("min_X (nm): ", minX,1);
            gd.addNumericField("max_X (nm): ", maxX,1);
            
            gd.addNumericField("min_Y (nm): ", minY,1);
            gd.addNumericField("max_Y (nm): ", maxY,1);
            
            gd.addNumericField("min_Z (nm): ", minZ,1);
            gd.addNumericField("max_Z (nm): ", maxZ,1);*/
            
            
            
            //gd.addNumericField("min_chi2: ", minScore,2,6,"<1");
            gd.addNumericField("max_chi2: ", maxScore,2,6,">1");
            
            

            gd.addNumericField("max_CRLB_X (nm): ", maxCRLBX,1,6,"(nm)");
            
            gd.addNumericField("max_CRLB_Y (nm): ", maxCRLBY,1,6,"(nm)");
            
            gd.addNumericField("max_CRLB_Z (nm): ", maxCRLBZ,1,6,"(nm)");
            
            
            

            gd.showDialog();
            
            
            if (!gd.wasCanceled()){


                /*minX=gd.getNextNumber();
                maxX=gd.getNextNumber();
                
                minY=gd.getNextNumber();
                maxY=gd.getNextNumber();
                
                minZ=gd.getNextNumber();
                maxZ=gd.getNextNumber();*/
                
                //minScore=gd.getNextNumber();
                maxScore=gd.getNextNumber();
                
                maxCRLBX=gd.getNextNumber();
                maxCRLBY=gd.getNextNumber();
                maxCRLBZ=gd.getNextNumber();
                StackLocalization slb = new StackLocalization();
                int nbtot=0;
                int nbafter=0;
                for (int i=0;i<sl.fl.size();i++){
                    int k=0;
                    FrameLocalization flb=new FrameLocalization(sl.fl.get(i).numFrame);
                    for (int j=0;j<sl.fl.get(i).loc.size();j++){
                        PLocalization p = sl.fl.get(i).loc.get(j);
                        if ((p.X>=minX2)&&(p.X<=maxX2)&&(p.Y>=minY2)&&(p.Y<=maxY2)&&(p.Z>=minrange+minZ)&&(p.Z<=maxrange+minZ)&&((p.score<=maxScore)||(p.score<0))&&((p.crlb_X<=maxCRLBX)||(p.crlb_X<0))&&((p.crlb_Y<=maxCRLBY)||(p.crlb_Y<0))&&((p.crlb_Z<=maxCRLBZ)||(p.crlb_Z<0))){
                            PLocalization pb= new PLocalization();
                            pb=p.copy();
                            flb.loc.add(pb);
                            nbafter++;
                        }
                        else{
                            nbtot++;
                        }
                    }
                    if (flb.loc.size()>0){
                        slb.fl.add(flb);
                    }
                }
                nbtot+=nbafter;
                sl=slb;
                IJ.log("Localization number after filtering: "+nbafter+"/"+nbtot);
                ZRendering.colorRendering(null,sl, 20,0,true);

                
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    void inversefiltering(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            double pixSize=20;
            ImagePlus imp=ZRendering.scatter2D(sl, pixSize,true);
            IJ.run("Brightness/Contrast...");
            double minrange=imp.getDisplayRangeMin();
            double maxrange=imp.getDisplayRangeMax();
            
            int xinit=0;
            int yinit=0;
            int width=imp.getWidth();
            int height=imp.getHeight();
            new WaitForUserDialog("ROI selection", "Select a rectangle in the image to control X/Y range, and adjust brightness/Cont. to control Z range, Then press OK.").show();
            if (imp.isVisible()){
                Roi r=imp.getRoi();
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    xinit=bound.x;
                    yinit=bound.y;
                    width=bound.width;
                    height=bound.height;
                }
                minrange=imp.getDisplayRangeMin();
                maxrange=imp.getDisplayRangeMax();
            }
            
            
            int offFrame=0;
            GenericDialog gd = new GenericDialog("ZOLA: Filtering");
            double minX=Double.POSITIVE_INFINITY,maxX=Double.NEGATIVE_INFINITY,minY=Double.POSITIVE_INFINITY,maxY=Double.NEGATIVE_INFINITY,minZ=Double.POSITIVE_INFINITY,maxZ=Double.NEGATIVE_INFINITY;
            double maxCRLBX=Double.NEGATIVE_INFINITY,maxCRLBY=Double.NEGATIVE_INFINITY,maxCRLBZ=Double.NEGATIVE_INFINITY;
            double maxScore=Double.NEGATIVE_INFINITY;
            
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    PLocalization p = sl.fl.get(i).loc.get(j);
                    
                    //if (p.score<minScore){
                    //    minScore=p.score;
                    //}
                    if (p.score>maxScore){
                        maxScore=p.score;
                    }
                    if (p.X<minX){
                        minX=p.X;
                    }
                    if (p.Y<minY){
                        minY=p.Y;
                    }
                    if (p.Z<minZ){
                        minZ=p.Z;
                    }
                    if (p.X>maxX){
                        maxX=p.X;
                    }
                    if (p.Y>maxY){
                        maxY=p.Y;
                    }
                    if (p.Z>maxZ){
                        maxZ=p.Z;
                    }
                    
                    
                    if (p.crlb_X>maxCRLBX){
                        maxCRLBX=p.crlb_X;
                    }
                    if (p.crlb_Y>maxCRLBY){
                        maxCRLBY=p.crlb_Y;
                    }
                    if (p.crlb_Z>maxCRLBZ){
                        maxCRLBZ=p.crlb_Z;
                    }
                }
            }
            
            double minX2=xinit*pixSize+minX;
            double maxX2=xinit*pixSize+width*pixSize+minX;
            double minY2=yinit*pixSize+minY;
            double maxY2=yinit*pixSize+height*pixSize+minY;
            
            
            
            if (true){


                /*minX=gd.getNextNumber();
                maxX=gd.getNextNumber();
                
                minY=gd.getNextNumber();
                maxY=gd.getNextNumber();
                
                minZ=gd.getNextNumber();
                maxZ=gd.getNextNumber();*/
                
                //minScore=gd.getNextNumber();
                maxScore=gd.getNextNumber();
                
                maxCRLBX=gd.getNextNumber();
                maxCRLBY=gd.getNextNumber();
                maxCRLBZ=gd.getNextNumber();
                StackLocalization slb = new StackLocalization();
                int nbtot=0;
                int nbafter=0;
                for (int i=0;i<sl.fl.size();i++){
                    int k=0;
                    FrameLocalization flb=new FrameLocalization(sl.fl.get(i).numFrame);
                    for (int j=0;j<sl.fl.get(i).loc.size();j++){
                        PLocalization p = sl.fl.get(i).loc.get(j);
                        if ((p.X>=minX2)&&(p.X<=maxX2)&&(p.Y>=minY2)&&(p.Y<=maxY2)&&(p.Z>=minrange+minZ)&&(p.Z<=maxrange+minZ)){
                            nbtot++;
                        }
                        else{
                            PLocalization pb= new PLocalization();
                            pb=p.copy();
                            flb.loc.add(pb);
                            nbafter++;
                            
                        }
                    }
                    if (flb.loc.size()>0){
                        slb.fl.add(flb);
                    }
                }
                nbtot+=nbafter;
                sl=slb;
                IJ.log("Localization number after filtering: "+nbafter+"/"+nbtot);
                ZRendering.colorRendering(null,sl, 20,0,true);

                
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    void deconvolution(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            
            ZRendering.hist2D(sl, 10, 0);
            Deconvolution dec = new Deconvolution(sl);
            sl=dec.runMulti(10);
            ZRendering.hist2D(sl, 10, 0);
            
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    void stathistoFilter(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            
            
            GenericDialog gd=null;
                
                //since NonBlockingGenericDialog incompatible with headless mode:
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())){
                    gd = new NonBlockingGenericDialog("ZOLA: Filtering");
                }
                else{
                    gd = new GenericDialog("ZOLA: Filtering");
                }
            PLocalization p=null;
            
            looper:for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    p=sl.fl.get(i).loc.get(j);
                    break looper;
                }
            }
            if (p!=null){
                int variableNumber=p.getNumberVariable();
                

                String [] fieldVariable =new String[variableNumber];
                for (int i=0;i<variableNumber;i++){
                    fieldVariable[i]=p.getLabel(i);
                }
                if (idVariable>=fieldVariable.length){
                    idVariable=fieldVariable.length-1;
                }
                if (idVariable<0){
                    idVariable=0;
                }
                
                
                double mini=Double.MAX_VALUE;
                double maxi=Double.MIN_VALUE;
                
                
                for (int i=0;i<sl.fl.size();i++){
                    for (int j=0;j<sl.fl.get(i).loc.size();j++){
                        PLocalization pp = sl.fl.get(i).loc.get(j);
                        double val=pp.getValueVariable(idVariable);
                        if (val<mini){
                            mini=val;
                        }
                        if (val>maxi){
                            maxi=val;
                        }
                        
                    }
                }
                
                
                
                gd.addChoice("Variable: ", fieldVariable,  fieldVariable[idVariable]);
                
                gd.addNumericField("Minimum: ", mini,2,12,"");
                gd.addNumericField("Maximum: ", maxi,2,12,"");
                
                
                TextField tf1,tf2;
                Choice choicefield=null;
                // Add a mouse listener to the config file field 
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
                {
                    
                    Vector<Choice> chois = (Vector<Choice>) gd.getChoices();

                    choicefield = chois.get(0); 
                    //MouseOptionSelect mosel=new MouseOptionSelect(choicefield,fieldVariable);
                    //choicefield.addMouseListener(mosel); 
                    
                    Vector<TextField> texts = (Vector<TextField>) gd.getNumericFields(); 
                    
                    tf1 = texts.get(0); 
                    tf2 = texts.get(1); 
                    
                    ShowHistogramButton but = new ShowHistogramButton(choicefield,tf1,tf2); 
                    gd.add(but); 
                    
                    
                    

                }
                
                
                
                
                
                
                
                
                gd.showDialog();


                if (!gd.wasCanceled()){
                    idVariable=gd.getNextChoiceIndex();
                    mini=gd.getNextNumber();
                    maxi=gd.getNextNumber();
                }

                IJ.log("mini "+mini+"  "+maxi);


                if (!gd.wasCanceled()){
                    
                    
                    StackLocalization slb = new StackLocalization();
                    int nbtot=0;
                    int nbafter=0;
                    for (int i=0;i<sl.fl.size();i++){
                        int k=0;
                        FrameLocalization flb=new FrameLocalization(sl.fl.get(i).numFrame);
                        for (int j=0;j<sl.fl.get(i).loc.size();j++){
                            PLocalization pp = sl.fl.get(i).loc.get(j);
                            double val=pp.getValueVariable(idVariable);
                            if ((val>=mini)&&(val<=maxi)){
                                PLocalization pb= new PLocalization();
                                pb=pp.copy();
                                flb.loc.add(pb);
                                nbafter++;
                            }
                            else{
                                nbtot++;
                            }
                        }
                        if (flb.loc.size()>0){
                            slb.fl.add(flb);
                        }
                    }
                    nbtot+=nbafter;
                    sl=slb;
                    IJ.log("Localization number after filtering: "+nbafter+"/"+nbtot);
                    ZRendering.colorRendering(null,sl, 20,0,true);
                    

                }
            }
            else{
                IJ.log("computation not possible: No localization in the table.");
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    
    void stathisto(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            
            GenericDialog gd = new GenericDialog("ZOLA: Statistic");
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            PLocalization p=null;
            
            looper:for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    p=sl.fl.get(i).loc.get(j);
                    break looper;
                }
            }
            if (p!=null){
                int variableNumber=p.getNumberVariable();
                

                String [] fieldVariable =new String[variableNumber];
                for (int i=0;i<variableNumber;i++){
                    fieldVariable[i]=p.getLabel(i);
                }
                if (idVariable>=fieldVariable.length){
                    idVariable=fieldVariable.length-1;
                }
                if (idVariable<0){
                    idVariable=0;
                }



                gd.addChoice("Variable histogram", fieldVariable,  fieldVariable[idVariable]);
                gd.addNumericField("Bin_number: ", this.binNumberVariable,2,6,"");


                gd.showDialog();


                if (!gd.wasCanceled()){
                    idVariable=gd.getNextChoiceIndex();
                    binNumberVariable=(int)gd.getNextNumber();
                }




                if (!gd.wasCanceled()){
                    
                    Statistic.computeHist(sl,idVariable,binNumberVariable);
                    //ZRendering.smartColorRendering(null,sl, 20,1,true);


                }
            }
            else{
                IJ.log("computation not possible: No localization in the table.");
            }
        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    
    
    void statmean(){
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        
        
        if (sl!=null){
            
                    
                Statistic.computeMean(sl);
                //ZRendering.smartColorRendering(null,sl, 20,1,true);

        }
        else{
            IJ.log("Please, import first a table");
        }
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
    
    
    void scatterPlot(){
        
        if (sl!=null){
            
            double minX=Double.POSITIVE_INFINITY;
            double maxX=Double.NEGATIVE_INFINITY;

            double minY=Double.POSITIVE_INFINITY;
            double maxY=Double.NEGATIVE_INFINITY;

            double minZ=Double.POSITIVE_INFINITY;
            double maxZ=Double.NEGATIVE_INFINITY;
            double x;
            double y;
            double z;
            
            
            
                
                
                
                
                
                

            //maybe use Arrays.sort to be fast
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        x=sl.fl.get(i).loc.get(j).X;
                        y=sl.fl.get(i).loc.get(j).Y;
                        z=sl.fl.get(i).loc.get(j).Z;
                        if (x<minX){
                            minX=x;
                        }
                        if (y<minY){
                            minY=y;
                        }
                        if (x>maxX){
                            maxX=x;
                        }
                        if (y>maxY){
                            maxY=y;
                        }
                        if (z<minZ){
                            minZ=z;
                        }
                        if (z>maxZ){
                            maxZ=z;
                        }
                    }
                }
            }
            
            PLocalization p=null;
            looper:for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    p=sl.fl.get(i).loc.get(j);
                    break looper;
                }
            }
            if (p!=null){
                int variableNumber=p.getNumberVariable();
                String [] fieldVariable =new String[variableNumber];
                
                for (int i=0;i<variableNumber;i++){
                    fieldVariable[i]=p.getLabel(i);
                }
                if (idVariable>=fieldVariable.length){
                    idVariable=fieldVariable.length-1;
                }
                if (idVariable<0){
                    idVariable=0;
                }

                GenericDialog gd = new GenericDialog("ZOLA: Render scatter plot");
                Font font = gd.getFont();
                Font fontBold= gd.getFont();
                try{
                    fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
                }catch(Exception e){}

                gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");

                gd.addChoice("Variable: ", fieldVariable,  fieldVariable[idVariable]);
                
                gd.addNumericField("Min_X: ", Math.floor(minX),0,6,"(nm)");
                gd.addNumericField("Max_X: ", Math.ceil(maxX),0,6,"(nm)");
                gd.addNumericField("Min_Y: ", Math.floor(minY),0,6,"(nm)");
                gd.addNumericField("Max_Y: ", Math.ceil(maxY),0,6,"(nm)");
                gd.addNumericField("Min_Z: ", Math.floor(minZ),0,6,"(nm)");
                gd.addNumericField("Max_Z: ", Math.ceil(maxZ),0,6,"(nm)");
                
                

                gd.showDialog();

                if (!gd.wasCanceled()){

                    sizeRendering = (double)gd.getNextNumber();
                    idVariable=gd.getNextChoiceIndex();
                    minX = (double)gd.getNextNumber();
                    maxX = (double)gd.getNextNumber();
                    minY = (double)gd.getNextNumber();
                    maxY = (double)gd.getNextNumber();
                    minZ = (double)gd.getNextNumber();
                    maxZ = (double)gd.getNextNumber();
                }
                else{
                    return;
                }

                ZRendering.scatter3D(sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,idVariable);
                
            }
            else{
                IJ.log("please import a localization table first");
            }
        }
        
    }
    
    void hist(){
        
        if (sl!=null){
            
            double sizeRenderingZ=sizeRendering;
            double minX=Double.POSITIVE_INFINITY;
            double maxX=Double.NEGATIVE_INFINITY;

            double minY=Double.POSITIVE_INFINITY;
            double maxY=Double.NEGATIVE_INFINITY;

            double minZ=Double.POSITIVE_INFINITY;
            double maxZ=Double.NEGATIVE_INFINITY;
            double x;
            double y;
            double z;

            //maybe use Arrays.sort to be fast
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        x=sl.fl.get(i).loc.get(j).X;
                        y=sl.fl.get(i).loc.get(j).Y;
                        z=sl.fl.get(i).loc.get(j).Z;
                        if (x<minX){
                            minX=x;
                        }
                        if (y<minY){
                            minY=y;
                        }
                        if (x>maxX){
                            maxX=x;
                        }
                        if (y>maxY){
                            maxY=y;
                        }
                        if (z<minZ){
                            minZ=z;
                        }
                        if (z>maxZ){
                            maxZ=z;
                        }
                    }
                }
            }
            
            GenericDialog gd = new GenericDialog("ZOLA: Render histogram");
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            gd.addNumericField("Pixel_size: ", sizeRendering,1,6,"(nm)");
            gd.addNumericField("Axial_pixel_size: ", sizeRenderingZ,1,6,"(nm)");
            gd.addCheckbox("3D_rendering :", is3Drendering);
            
            gd.addNumericField("Shift_histogram (px): ", shiftrendering,0);
            
            gd.addMessage("Optional parameters");
            
            gd.addNumericField("Min_X: ", Math.floor(minX),0,6,"(nm)");
            gd.addNumericField("Max_X: ", Math.ceil(maxX),0,6,"(nm)");
            gd.addNumericField("Min_Y: ", Math.floor(minY),0,6,"(nm)");
            gd.addNumericField("Max_Y: ", Math.ceil(maxY),0,6,"(nm)");
            gd.addNumericField("Min_Z: ", Math.floor(minZ),0,6,"(nm)");
            gd.addNumericField("Max_Z: ", Math.ceil(maxZ),0,6,"(nm)");
            
            
            gd.showDialog();

            if (!gd.wasCanceled()){

                sizeRendering = (double)gd.getNextNumber();
                sizeRenderingZ = (double)gd.getNextNumber();
                is3Drendering = (boolean)gd.getNextBoolean();
                shiftrendering = (int)gd.getNextNumber();
                
                minX = (double)gd.getNextNumber();
                maxX = (double)gd.getNextNumber();
                minY = (double)gd.getNextNumber();
                maxY = (double)gd.getNextNumber();
                minZ = (double)gd.getNextNumber();
                maxZ = (double)gd.getNextNumber();
                
            }
            else{
                return;
            }
            if (shiftrendering<0){
                shiftrendering=0;
                IJ.log("render shift should be at least 0");
            }
            if (is3Drendering){
                
                ZRendering.hist3D(sl,sizeRendering,sizeRenderingZ,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering);
            }
            else{
                ZRendering.hist2D(sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering);
            }
            
        }
        else{
            IJ.log("please import a localization table first");
        }
    }
    
    
    
    
    
    void timeRender(){
        
        if (sl!=null){
            
            
            double minX=Double.POSITIVE_INFINITY;
            double maxX=Double.NEGATIVE_INFINITY;

            double minY=Double.POSITIVE_INFINITY;
            double maxY=Double.NEGATIVE_INFINITY;

            double minZ=Double.POSITIVE_INFINITY;
            double maxZ=Double.NEGATIVE_INFINITY;
            double x;
            double y;
            double z;

            //maybe use Arrays.sort to be fast
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        x=sl.fl.get(i).loc.get(j).X;
                        y=sl.fl.get(i).loc.get(j).Y;
                        z=sl.fl.get(i).loc.get(j).Z;
                        if (x<minX){
                            minX=x;
                        }
                        if (y<minY){
                            minY=y;
                        }
                        if (x>maxX){
                            maxX=x;
                        }
                        if (y>maxY){
                            maxY=y;
                        }
                        if (z<minZ){
                            minZ=z;
                        }
                        if (z>maxZ){
                            maxZ=z;
                        }
                    }
                }
            }
            
            GenericDialog gd = new GenericDialog("ZOLA: Render color image");


            gd.addNumericField("pixel size: ", sizeRendering,1,6,"(nm)");
            
            gd.addCheckbox("Show_calibration_bar :", showLUT);
            
            
            gd.addNumericField("frame_number_per_split: ", frameSplit,0,6,"(frames)");
            
            gd.addMessage("Optional parameters");
            gd.addNumericField("min_X: ", Math.floor(minX),0,6,"(nm)");
            gd.addNumericField("max_X: ", Math.ceil(maxX),0,6,"(nm)");
            gd.addNumericField("min_Y: ", Math.floor(minY),0,6,"(nm)");
            gd.addNumericField("max_Y: ", Math.ceil(maxY),0,6,"(nm)");
            gd.addNumericField("min_Z: ", Math.floor(minZ),0,6,"(nm)");
            gd.addNumericField("max_Z: ", Math.ceil(maxZ),0,6,"(nm)");
            
            
            gd.showDialog();

            if (!gd.wasCanceled()){

                sizeRendering = (double)gd.getNextNumber();
                showLUT = (boolean)gd.getNextBoolean();
                
                frameSplit = (int)gd.getNextNumber();
                minX = (double)gd.getNextNumber();
                maxX = (double)gd.getNextNumber();
                minY = (double)gd.getNextNumber();
                maxY = (double)gd.getNextNumber();
                minZ = (double)gd.getNextNumber();
                maxZ = (double)gd.getNextNumber();
                
            }
            else{
                return;
            }
            
            
           
            ZRendering.timeRendering2D(sl, sizeRendering, minX, maxX, minY, maxY, minZ, maxZ, frameSplit);
            
            
            
        }
        else{
            IJ.log("please import a localization table first");
        }
    }
    
    
    
//    void smartColorHist(){
//        
//        if (sl!=null){
//            
//            
//            double minX=Double.POSITIVE_INFINITY;
//            double maxX=Double.NEGATIVE_INFINITY;
//
//            double minY=Double.POSITIVE_INFINITY;
//            double maxY=Double.NEGATIVE_INFINITY;
//
//            double minZ=Double.POSITIVE_INFINITY;
//            double maxZ=Double.NEGATIVE_INFINITY;
//            double x;
//            double y;
//            double z;
//
//            //maybe use Arrays.sort to be fast
//            for (int i=0;i<sl.fl.size();i++){
//                for (int j=0;j<sl.fl.get(i).loc.size();j++){
//                    if (sl.fl.get(i).loc.get(j).exists){
//                        x=sl.fl.get(i).loc.get(j).X;
//                        y=sl.fl.get(i).loc.get(j).Y;
//                        z=sl.fl.get(i).loc.get(j).Z;
//                        if (x<minX){
//                            minX=x;
//                        }
//                        if (y<minY){
//                            minY=y;
//                        }
//                        if (x>maxX){
//                            maxX=x;
//                        }
//                        if (y>maxY){
//                            maxY=y;
//                        }
//                        if (z<minZ){
//                            minZ=z;
//                        }
//                        if (z>maxZ){
//                            maxZ=z;
//                        }
//                    }
//                }
//            }
//            
//            GenericDialog gd = new GenericDialog("ZOLA: Render color image");
//            Font font = gd.getFont();
//            Font fontBold= gd.getFont();
//            try{
//                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
//            }catch(Exception e){}
//            
//            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
//            gd.addNumericField("Shift_histogram: ", shiftrendering,0,6,"(pixels)");
//            gd.addCheckbox("Show_calibration_bar:", showLUT);
//            gd.showDialog();
//
//            if (!gd.wasCanceled()){
//
//                sizeRendering = (double)gd.getNextNumber();
//                shiftrendering = (int)gd.getNextNumber();
//                showLUT = (boolean)gd.getNextBoolean();
//            }
//            else{
//                return;
//            }
//            if (shiftrendering<0){
//                shiftrendering=0;
//                IJ.log("render shift should be at least 0");
//            }
//            
//            
//            
//            ZRendering.colorRendering(null,sl,sizeRendering,shiftrendering,showLUT);
//            
//            
//        }
//        else{
//            IJ.log("please import a localization table first");
//        }
//    }
    
    
    
    
    
    
    void colorHist(){
        
        if (sl!=null){
            
            String [] lutstmp=ZRendering.getFijiLUTs();
            
            String [] luts = new String[lutstmp.length+1];
            for (int i=0;i<lutstmp.length;i++){
                luts[i+1]=lutstmp[i];
            }
            luts[0]=ZRendering.lut;
            
            
            double minX=Double.POSITIVE_INFINITY;
            double maxX=Double.NEGATIVE_INFINITY;

            double minY=Double.POSITIVE_INFINITY;
            double maxY=Double.NEGATIVE_INFINITY;

            double minZ=Double.POSITIVE_INFINITY;
            double maxZ=Double.NEGATIVE_INFINITY;
            double x;
            double y;
            double z;

            //maybe use Arrays.sort to be fast
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        x=sl.fl.get(i).loc.get(j).X;
                        y=sl.fl.get(i).loc.get(j).Y;
                        z=sl.fl.get(i).loc.get(j).Z;
                        if (x<minX){
                            minX=x;
                        }
                        if (y<minY){
                            minY=y;
                        }
                        if (x>maxX){
                            maxX=x;
                        }
                        if (y>maxY){
                            maxY=y;
                        }
                        if (z<minZ){
                            minZ=z;
                        }
                        if (z>maxZ){
                            maxZ=z;
                        }
                    }
                }
            }
            
            
            //IJ.log("WARNING: code line to delete 3000");
            //maxZ=minZ+3000;
            
            GenericDialog gd = new GenericDialog("ZOLA: Render color image");
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
            
            gd.addChoice("Lookup_table", luts,  luts[0]);
            
            gd.addCheckbox("Show_calibration_bar :", showLUT);
            
            gd.addNumericField("Shift_histogram: ", shiftrendering,0,6,"(pixels)");
            
            gd.addMessage("Optional parameters");
            gd.addNumericField("Min_X: ", Math.floor(minX),0,6,"(nm)");
            gd.addNumericField("Max_X: ", Math.ceil(maxX),0,6,"(nm)");
            gd.addNumericField("Min_Y: ", Math.floor(minY),0,6,"(nm)");
            gd.addNumericField("Max_Y: ", Math.ceil(maxY),0,6,"(nm)");
            gd.addNumericField("Min_Z: ", Math.floor(minZ),0,6,"(nm)");
            gd.addNumericField("Max_Z: ", Math.ceil(maxZ),0,6,"(nm)");
            
            
            gd.showDialog();
            String lut=ZRendering.lut;
            if (!gd.wasCanceled()){

                sizeRendering = (double)gd.getNextNumber();
                lut=gd.getNextChoice();
                showLUT = (boolean)gd.getNextBoolean();
                shiftrendering = (int)gd.getNextNumber();
                minX = (double)gd.getNextNumber();
                maxX = (double)gd.getNextNumber();
                minY = (double)gd.getNextNumber();
                maxY = (double)gd.getNextNumber();
                minZ = (double)gd.getNextNumber();
                maxZ = (double)gd.getNextNumber();
                
            }
            else{
                return;
            }
            if (shiftrendering<0){
                shiftrendering=0;
                IJ.log("render shift should be at least 0");
            }
            
            
            
            
            ZRendering.colorRendering(null,sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering,showLUT,lut);
            
            
        }
        else{
            IJ.log("please import a localization table first");
        }
    }
    
    
    
//    
//    void scatter(){
//        
//        if (sl!=null){
//            
//            String [] lutstmp=ZRendering.getFijiLUTs();
//            
//            String [] luts = new String[lutstmp.length+1];
//            for (int i=0;i<lutstmp.length;i++){
//                luts[i+1]=lutstmp[i];
//            }
//            luts[0]=ZRendering.lut;
//            
//            
//            double minX=Double.POSITIVE_INFINITY;
//            double maxX=Double.NEGATIVE_INFINITY;
//
//            double minY=Double.POSITIVE_INFINITY;
//            double maxY=Double.NEGATIVE_INFINITY;
//
//            double minZ=Double.POSITIVE_INFINITY;
//            double maxZ=Double.NEGATIVE_INFINITY;
//            double x;
//            double y;
//            double z;
//
//            //maybe use Arrays.sort to be fast
//            for (int i=0;i<sl.fl.size();i++){
//                for (int j=0;j<sl.fl.get(i).loc.size();j++){
//                    if (sl.fl.get(i).loc.get(j).exists){
//                        x=sl.fl.get(i).loc.get(j).X;
//                        y=sl.fl.get(i).loc.get(j).Y;
//                        z=sl.fl.get(i).loc.get(j).Z;
//                        if (x<minX){
//                            minX=x;
//                        }
//                        if (y<minY){
//                            minY=y;
//                        }
//                        if (x>maxX){
//                            maxX=x;
//                        }
//                        if (y>maxY){
//                            maxY=y;
//                        }
//                        if (z<minZ){
//                            minZ=z;
//                        }
//                        if (z>maxZ){
//                            maxZ=z;
//                        }
//                    }
//                }
//            }
//            
//            GenericDialog gd = new GenericDialog("ZOLA: Render color image");
//            Font font = gd.getFont();
//            Font fontBold= gd.getFont();
//            try{
//                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
//            }catch(Exception e){}
//            
//            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
//            
//            gd.addChoice("Lookup_table", luts,  luts[0]);
//            
//            gd.addCheckbox("Show_calibration_bar :", showLUT);
//            
//            gd.addNumericField("Shift_histogram: ", shiftrendering,0,6,"(pixels)");
//            
//            gd.addMessage("Optional parameters");
//            gd.addNumericField("Min_X: ", Math.floor(minX),0,6,"(nm)");
//            gd.addNumericField("Max_X: ", Math.ceil(maxX),0,6,"(nm)");
//            gd.addNumericField("Min_Y: ", Math.floor(minY),0,6,"(nm)");
//            gd.addNumericField("Max_Y: ", Math.ceil(maxY),0,6,"(nm)");
//            gd.addNumericField("Min_Z: ", Math.floor(minZ),0,6,"(nm)");
//            gd.addNumericField("Max_Z: ", Math.ceil(maxZ),0,6,"(nm)");
//            
//            
//            gd.showDialog();
//            String lut=ZRendering.lut;
//            if (!gd.wasCanceled()){
//
//                sizeRendering = (double)gd.getNextNumber();
//                lut=gd.getNextChoice();
//                showLUT = (boolean)gd.getNextBoolean();
//                shiftrendering = (int)gd.getNextNumber();
//                minX = (double)gd.getNextNumber();
//                maxX = (double)gd.getNextNumber();
//                minY = (double)gd.getNextNumber();
//                maxY = (double)gd.getNextNumber();
//                minZ = (double)gd.getNextNumber();
//                maxZ = (double)gd.getNextNumber();
//                
//            }
//            else{
//                return;
//            }
//            if (shiftrendering<0){
//                shiftrendering=0;
//                IJ.log("render shift should be at least 0");
//            }
//            
//            
//            
//            
//            ZRendering.colorRendering(null,sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering,showLUT,lut);
//            
//            
//        }
//        else{
//            IJ.log("please import a localization table first");
//        }
//    }
    
    
    
    void colorizeHist(){
        
        ImagePlus imp=IJ.getImage();
        String title = imp.getTitle();
        if (title.startsWith(ZRendering.nameHistPlot) ||title.startsWith(ZRendering.nameScatterPlot)){

            
            
            ZRendering.colorizeHist(imp);
            
        }
        else{
            //IJ.log("WARNING, this colorize function should be used from 3D histogram images");
            ZRendering.colorizeHist(imp);
        }
        
    }
    
    
    void getModels_GPUalloc(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            IJ.log("Sorry but it is not possible to import data during processing");
            return;
        }
        
        processing=new Boolean(true);
        
        //process
        
        
        TextField textField; 
        int sizeImage=sizePatch;
        
        GenericDialog gd = new GenericDialog("ZOLA: GPU memory allocation");
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
        gd.addStringField("File_path:",path_calibration ,sizeTextString); 
        gd.addNumericField("Model_number: ", model_number,0,6,"");
        
        gd.addNumericField("patch_size: ", sizePatch,0,6,"(pixels)");
        gd.addNumericField("image_size: ", sizeImage,0,6,"(pixels)");
        
        gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        
        
        
        
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
        }
        
        gd.showDialog();
        
        if (!gd.wasCanceled()){
            path_calibration = gd.getNextString(); 
            model_number = (int)gd.getNextNumber();
            sizePatch = (int)gd.getNextNumber();
            sizeImage = (int)gd.getNextNumber();
            nwat = gd.getNextNumber();
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        
        MyCudaStream.init(1);
        
        
        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
        if (!dp.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        dp.setSizeoutput(sizePatch);
        dp.setNwat(nwat);
        dp.setMany(model_number);
        
        
        
        
        
        
        mg=new org.pasteur.imagej.process.gpu.Model_generator_(dp,sizeImage);
        
        
        
        processing=new Boolean(false);
        
    }
    void getModels_GPUcomput(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            IJ.log("Sorry but it is not possible to import data during processing");
            return;
        }
        
        processing=new Boolean(true);
        
        //process
        
        
        TextField textField; 
        
        String stringvalue="";
        
        GenericDialog gd = new GenericDialog("ZOLA: GPU memory allocation");
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
        gd.addStringField("list_of_positions {x1,y1,z1,f1;x2,y2,z2,f2;...} (nm):",stringvalue ,sizeTextString); 
        
        gd.addStringField("Result_image_path:",path_result ,sizeTextString); 
        
        
        gd.showDialog();
        
        if (!gd.wasCanceled()){
            stringvalue = gd.getNextString(); 
            path_result = gd.getNextString(); 
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        String [] stringvalue2=stringvalue.split(";");
        String [] stringvalue3;
        if (stringvalue2.length!=0){
            
            int count=0;
            for (int i=0;i<stringvalue2.length;i++){
                stringvalue3=stringvalue2[i].split(",");
                if (stringvalue3.length!=4){
                    IJ.log("WARNING, only "+stringvalue3.length+" values instead of 4 are detected line "+i);
                }
                else{
                    count++;
                }
            }
            double [][] value = new double [count][4];
            count=0;
            for (int i=0;i<stringvalue2.length;i++){
                stringvalue3=stringvalue2[i].split(",");
                if (stringvalue3.length==4){
                    value[count][0]=Double.parseDouble(stringvalue3[0]);
                    value[count][1]=Double.parseDouble(stringvalue3[1]);
                    value[count][2]=Double.parseDouble(stringvalue3[2]);
                    value[count][3]=Double.parseDouble(stringvalue3[3]);
                    count++;
                }
            }
            mg.computePSF(value,path_result);
        }
        else{
            IJ.log("WARNING, 0 position detected");
        }
        
        processing=new Boolean(false);
    }
    void getModels_GPUfree(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            IJ.log("Sorry but it is not possible to import data during processing");
            return;
        }
        
        processing=new Boolean(true);
        
        //process
        mg.free();
        
        MyCudaStream.destroy();
        
        processing=new Boolean(false);
    }
    
    
    

    void importData(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            IJ.log("Sorry but it is not possible to import data during processing");
            return;
        }
        
        processing=new Boolean(true);
        
        TextField textField; 
        
        
        GenericDialog gd = new GenericDialog("ZOLA: Import localization table");

        
        gd.addStringField("File_path:",path_localization ,sizeTextString); 
        
        
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_localization,"Import localization table");
            textField.addMouseListener(mol); 
            
        }
        
        gd.showDialog();
        
        if (!gd.wasCanceled()){
            path_localization = gd.getNextString(); 
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        
        if (path_localization!=null){
            File f = new File(path_localization);
            if (f.exists()){
                sl=new StackLocalization(path_localization);
                
                
                
                /*for (int i=0;i<sl.fl.size();i++){
                    for (int j=0;j<sl.fl.get(i).loc.size();j++){
                        sl.fl.get(i).loc.get(j).X+=35.74;
                        sl.fl.get(i).loc.get(j).Y+=103.25;
                        sl.fl.get(i).loc.get(j).Z+=139.55;
                    }
                }*/
                
                
                ZRendering.colorRendering(null,sl, 20,0,true);


                IJ.log("localization table loaded");
            }
            else{
                IJ.log("Import impossible: File not found: "+path_localization);
            }
            
        }
        
        
        
        processing=new Boolean(false);
    }
    
    
    
    
    

    void appendData(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            IJ.log("Sorry but it is not possible to import data during processing");
            return;
        }
        
        processing=new Boolean(true);
        
        TextField textField; 
        
        
        GenericDialog gd = new GenericDialog("ZOLA: Append localization table");

        
        gd.addStringField("File_path:",path_localization ,sizeTextString); 
        
        
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_localization,"Import localization table");
            textField.addMouseListener(mol); 
            
        }
        
        gd.showDialog();
        
        if (!gd.wasCanceled()){
            path_localization = gd.getNextString(); 
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        
        if (path_localization!=null){
            File f = new File(path_localization);
            if (f.exists()){
                if (sl==null){
                    sl=new StackLocalization(path_localization);
                }
                else{
                    sl.append(path_localization);
                }
                ZRendering.colorRendering(null,sl, 20,0,true);
        
                IJ.log("localization table appended");
            }
            else{
                IJ.log("Import impossible: File not found: "+path_localization);
            }
        }
        
        
        
        processing=new Boolean(false);
    }
    
    
    
    void fuseData(){
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            IJ.log("Sorry but it is not possible to import data during processing");
            return;
        }
        
        processing=new Boolean(true);
        
        TextField textField; 
        
        
        GenericDialog gd = new GenericDialog("ZOLA: Fuse localization table");

        
        gd.addStringField("File_path:",path_localization ,sizeTextString); 
        
        
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_localization,"Import localization table");
            textField.addMouseListener(mol); 
            
        }
        
        gd.showDialog();
        
        if (!gd.wasCanceled()){
            path_localization = gd.getNextString(); 
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        
        if (path_localization!=null){
            File f = new File(path_localization);
            if (f.exists()){
                if (sl==null){
                    sl=new StackLocalization(path_localization);
                }
                else{
                    sl.fuse(path_localization);
                }
                ZRendering.colorRendering(null,sl, 20,0,true);
        
                IJ.log("localization table fused");
            }
            else{
                IJ.log("Import impossible: File not found: "+path_localization);
            }
        }
        
        
        
        processing=new Boolean(false);
    }
    
    
    
    void exportData(){
        
        boolean nopressed=true;
        
        while (nopressed){
            TextField textField; 


            GenericDialog gd = new GenericDialog("ZOLA: Export localization table");


            gd.addStringField("File_path:",path_localization,sizeTextString); 


            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;

                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
                textField = texts.get(t++); 
                MouseOptionSave mos=new MouseOptionSave(textField,path_localization,"Export localization table");
                textField.addMouseListener(mos); 

            }



            gd.showDialog();

            if (!gd.wasCanceled()){
                path_localization = gd.getNextString(); 
            }
            else{
                return;
            }




            if (path_localization!=null){
                boolean ok=true;
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
                {
                    File f = new File(path_localization);
                    if (f.exists()){

                        YesNoCancelDialog yn = new YesNoCancelDialog(IJ.getInstance(),"WARNING", "file already exists. Do you want to overwrite the file ?");
                        
                        if (yn.cancelPressed()){
                            //IJ.log("cancel pressed");
                            nopressed=false;
                            ok=false;
                        }
                        else if (yn.yesPressed()){
                            //IJ.log("yes pressed");
                        }
                        else{
                            //IJ.log("no pressed");
                            ok=false;
                        }
                    }
                }
                if (ok){
                    if (sl!=null){
                        sl.save(path_localization);
                    }
                    IJ.log("localization table saved");
                    nopressed=false;
                }
                else if (nopressed==false){
                    IJ.log("localization table NOT saved");
                }

            }
        }
        
        
        
        
    }
    
    
    
    
    
    void dualcamregistration(){
        
        
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        
        
        String pathRes1=path_localization; 
        String pathRes2=path_localization2; 
        String pathRes1_=path_localization; 
        
        int nn=(pathRes1.lastIndexOf("."));
        if ((nn>pathRes1.length()-6)&&(nn>=0)){
            pathRes1=pathRes1.substring(0,nn);
        }
        pathRes1=pathRes1+"_registredResult.csv";
        
        
        nn=(pathRes2.lastIndexOf("."));
        if ((nn>pathRes2.length()-6)&&(nn>=0)){
            pathRes2=pathRes2.substring(0,nn);
        }
        pathRes2=pathRes2+"_registredResult.csv";
        
        nn=(pathRes1_.lastIndexOf("."));
        if ((nn>pathRes1_.length()-6)&&(nn>=0)){
            pathRes1_=pathRes1_.substring(0,nn);
        }
        pathRes1_=pathRes1_+"_aggregatedResult.csv";
        
        
        processing=new Boolean(true);


        /*String [] imTitle=WindowManager.getImageTitles();

        if (imTitle.length!=2){
            IJ.log("Please, open 2 images first");
            processing=new Boolean(false);
        }*/


        
        



        GenericDialog gd = new GenericDialog(" 2 cameras registration");
        
        gd.addMessage("At this stage, we assume the drift has been corrected");

        gd.addMessage("Input localization tables");

        gd.addStringField("Camera_1:",path_localization,sizeTextString); 
        
        

        gd.addStringField("Camera_2:",path_localization2,sizeTextString); 

        
        
        gd.addMessage("Fusion parameters:");
        
        gd.addNumericField("max._X/Y distance between cam1 and cam2 localizations (nm): ", dualCamMaxDistanceMergingPlan,1);
        
        
        
        gd.addMessage("Output localization tables");

        gd.addStringField("Camera_1_(registred):",pathRes1,sizeTextString); 
        

        gd.addStringField("Camera_2_(registred):",pathRes2,sizeTextString); 
        
        gd.addStringField("Camera_1_+_2_(fusion):",pathRes1_,sizeTextString); 
        
        
        
        TextField textField; 

        TextField textField2; 
        TextField textField3; 
        TextField textField4; 
        TextField textField5; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;

            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_localization,"Import localization table (cam 1)");
            textField.addMouseListener(mol); 

            textField2 = texts.get(t++); 
            mol=new MouseOptionLoad(textField2,path_localization,"Import localization table (cam 2)");
            textField2.addMouseListener(mol); 
            
            
            textField3 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField3,path_localization,"Save unfused localization table (cam 1)");
            textField3.addMouseListener(mos); 
            
            
            textField4 = texts.get(t++); 
            mos=new MouseOptionSave(textField4,path_localization,"Save unfused localization table (cam 2)");
            textField4.addMouseListener(mos); 
            
            textField5 = texts.get(t++); 
            mos=new MouseOptionSave(textField5,path_localization,"Save fused localization table (cam 1+2)");
            textField5.addMouseListener(mos); 

            

        }


        gd.showDialog();



        if (!gd.wasCanceled()){


            path_localization = gd.getNextString(); 
            
            path_localization2 = gd.getNextString(); 

            //imSize=(int)gd.getNextNumber();
            //xystep=gd.getNextNumber();
            
            dualCamMaxDistanceMergingPlan=gd.getNextNumber();
            
            
            
            pathRes1 = gd.getNextString(); 
            
            pathRes2 = gd.getNextString(); 
            
            pathRes1_ = gd.getNextString(); 

        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        if ((path_localization.equals(path_localization2))||(path_localization.equals(pathRes1))||(path_localization.equals(pathRes2))||(path_localization.equals(pathRes1_))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((path_localization2.equals(pathRes1))||(path_localization2.equals(pathRes2))||(path_localization2.equals(pathRes1_))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((pathRes1.equals(pathRes2))||(pathRes1.equals(pathRes1_))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((pathRes2.equals(pathRes1_))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        
        
        StackLocalization sl1=new StackLocalization(path_localization);
        StackLocalization sl2=new StackLocalization(path_localization2);
        
        //double imSizeNM=((double)imSize)*xystep*1000;
        
        
        org.pasteur.imagej.postprocess.DualCamRegistration dcf=new org.pasteur.imagej.postprocess.DualCamRegistration(sl1,sl2,pathRes1,pathRes2,pathRes1_,dualCamMaxDistanceMergingPlan);
        dcf.run();
        IJ.log("well done");
        processing=new Boolean(false);
    }
    
    
    
    
    
    void dualcamlocalization(){
        
        
        nbStream=5;
        nbThread=100;  
        
        String [] images=WindowManager.getImageTitles();
        
        
        //String[] images = { "Try all",  "Huang",  "coucou" , "Huang"};
        
        if (images.length<2){
            IJ.log("Sorry but you should import at lease 2 stacks of images");
            return;
        }
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
        
        String pathRes=path_localization; 
        
        int nn=(pathRes.lastIndexOf("."));
        
        nn=(pathRes.lastIndexOf("."));
        if ((nn>pathRes.length()-6)&&(nn>=0)){
            pathRes=pathRes.substring(0,nn);
        }
        pathRes=pathRes+"_fusionResult.csv";
        
        
        


        /*String [] imTitle=WindowManager.getImageTitles();

        if (imTitle.length!=2){
            IJ.log("Please, open 2 images first");
            processing=new Boolean(false);
        }*/


        
        double gain2=gain;
        
        
        
        double zfocus1=0;
        double zfocus2=0;
        
        String path_calibration2=path_calibration;


        GenericDialog gd = new GenericDialog(" 2 cameras localization");
        
        gd.addMessage("At this stage, we assume the drift has been corrected");
        
        gd.addMessage("Camera 1");
        
        
        gd.addNumericField("Camera_Gain_1: ", gain, 4,6,"(1 if camera provides photon counts)");

        gd.addStringField("Calibration_file_1:",path_calibration,sizeTextString); 
                
        gd.addChoice("Image_cam_1:", images, images[0]);
        
        gd.addStringField("Localization_cam_1:",path_localization,sizeTextString); 
        
                

        gd.addNumericField("distance_focus_to_coverslip_1", zfocus1, 3,6,"(µm)");
                
                
                
        
        
        gd.addMessage("Camera 2");
        
            gd.addNumericField("Camera_Gain_2: ", gain2, 4,6,"(1 if camera provides photon counts)");
        
        gd.addStringField("Calibration_file_2:",path_calibration,sizeTextString); 
        
        gd.addChoice("Image_cam_2:", images, images[1]);
        
        gd.addStringField("Localization_cam_2:",path_localization,sizeTextString); 
        

        gd.addNumericField("distance_focus_to_coverslip_2", zfocus2, 3,6,"(µm)");

        
        
        gd.addMessage("Additional parameters");
        
        gd.addMessage("Camera_ADU (1&2): "+adu);
        gd.addMessage("Camera_Offset (1&2): "+ offset);
        gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        
        
        
        gd.addNumericField("max._X/Y distance between cam1 and cam2 localizations (nm): ", dualCamMaxDistanceMergingPlan,1);
        
        gd.addNumericField("patch size: ", sizePatch,0);
        
        
        gd.addMessage("Output localization table");
        
        gd.addStringField("fusion:",pathRes,sizeTextString); 
        
        
        
        TextField textField; 

        TextField textField2; 
        TextField textField3; 
        TextField textField4; 
        TextField textField5; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file (cam 1)");
            textField.addMouseListener(mol); 
                
                
            textField2 = texts.get(t++); 
            mol=new MouseOptionLoad(textField2,path_localization,"Import localization table (cam 1)");
            textField2.addMouseListener(mol); 
            
            
            textField3 = texts.get(t++); 
            mol=new MouseOptionLoad(textField3,path_calibration,"Import calibration file (cam 2)");
            textField3.addMouseListener(mol); 
            
            
            textField4 = texts.get(t++); 
            mol=new MouseOptionLoad(textField4,path_localization,"Import localization table (cam 2)");
            textField4.addMouseListener(mol); 
            
            
            textField5 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField5,pathRes,"Save fused localization table (cam 1+2)");
            textField5.addMouseListener(mos); 

            

        }


        gd.showDialog();

        String imageCam1;
        String imageCam2;

        if (!gd.wasCanceled()){
            
            
            gain=(double)gd.getNextNumber();
            
            
            path_calibration = gd.getNextString(); 
            imageCam1 = gd.getNextChoice();
        
            path_localization = gd.getNextString(); 
            
            
            zfocus1=(double)gd.getNextNumber();
            
            gain2=(double)gd.getNextNumber();
            
            path_calibration2 = gd.getNextString(); 
            imageCam2 = gd.getNextChoice();
            path_localization2 = gd.getNextString(); 
            zfocus2=(double)gd.getNextNumber();
            
            
            
            nwat=(double)gd.getNextNumber();
            
            
            dualCamMaxDistanceMergingPlan=gd.getNextNumber();
            
            
            sizePatch = (int)gd.getNextNumber();
            
            
            pathRes = gd.getNextString(); 

        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        if ((path_localization.equals(path_localization2))||(path_localization.equals(pathRes))||(path_localization.equals(path_calibration))||(path_localization.equals(path_calibration2))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((path_localization2.equals(pathRes))||(path_localization2.equals(path_calibration))||(path_localization2.equals(path_calibration2))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((pathRes.equals(path_calibration))||(pathRes.equals(path_calibration2))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((path_calibration.equals(path_calibration2))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        
        IJ.log("calib1 "+path_calibration);
        
        IJ.log("calib2 "+path_calibration2);
        
        IJ.log("loc1 "+path_localization);
        
        IJ.log("loc2 "+path_localization2);
        
        IJ.log("imageCam1 "+imageCam1);
        
        IJ.log("imageCam2 "+imageCam2);
        
        
        
        
        StackLocalization sl1=new StackLocalization(path_localization);
        StackLocalization sl2=new StackLocalization(path_localization2);
        
        MyCudaStream.init(nbStream*2);//nbStream*nbCam
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatch*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        org.pasteur.imagej.process.gpu.DataPhase_ dp1 = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
        if (!dp1.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        dp1.setSizeoutput(sizePatch);
        dp1.setNwat(nwat);
        dp1.param.Zfocus=zfocus1;
        
        org.pasteur.imagej.process.gpu.DataPhase_ dp2 = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration2);
        if (!dp2.loading){
            IJ.log("impossible to load "+path_calibration2);
            processing=new Boolean(false);
            return;
        }
        dp2.setSizeoutput(sizePatch);
        dp2.setNwat(nwat);
        dp2.param.Zfocus=zfocus2;
        //double imSizeNM=((double)imSize)*xystep*1000;
        
        
        org.pasteur.imagej.process.gpu.DualCamLocalizationPipeline_ dcf=new org.pasteur.imagej.process.gpu.DualCamLocalizationPipeline_(imageCam1,sl1,dp1,imageCam2,sl2,dp2,pathRes,dualCamMaxDistanceMergingPlan,nbStream, nbThread,adu,gain,offset,gain2);
        
        sl=new StackLocalization();
        
        dcf.run(sl);
        IJ.log("well done");
        org.pasteur.imagej.cuda.MyCudaStream.destroy();
        processing=new Boolean(false);
    }
    
    
    
    
    
    
    
    
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************           DEVELOPMENT                   ***************
//    ************             METHODS                     ***************
//    ************               ...                       ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
    
    
    
    
    void filamentResolution(){
        
        
        
        
        double sizePix=20;
        
        
        
        String pathResinit=path_localization; 
        
        
        int nn=(pathResinit.lastIndexOf("."));
        if ((nn>pathResinit.length()-6)&&(nn>=0)){
            pathResinit=pathResinit.substring(0,nn);
        }
        int num=0;
        String pathRes=pathResinit+"_withFilament_"+num+".csv";
        String pathResex=pathResinit+"_withoutFilament_"+num+".csv";
        File f=new File(pathRes);
        File ff=new File(pathResex);
        while (f.exists()||ff.exists()){
            num++;
            pathRes=pathResinit+"_withFilament_"+num+".csv";
            f=new File(pathRes);
            pathResex=pathResinit+"_withoutFilament_"+num+".csv";
            ff=new File(pathResex);
        }
        
        
        
        if (sl!=null){
            ZRendering.hist3D(sl, sizePix,0);
        }
        else{
            IJ.log("please import a localization table first");
            return;
        }
        
        
        new WaitForUserDialog("Line selection", "Select a line and add it in RoiManager (ctrl+t) / Then press OK.").show();
        
        

        GenericDialog gd = new GenericDialog("ZOLA: Resolution computation");

        
        gd.addNumericField("maximum filament radius (nm): ", maximumDistance,1);
        gd.addNumericField("filament polynomiaml fit order: ", orderFit,0);
        gd.addNumericField("histogram bin size: ", binSize,1);
        
        gd.addStringField("Result localization table of filament :","",sizeTextString); 
        
        gd.addStringField("Result localization table without filament :","",sizeTextString); 
        
            
        TextField textField;
        TextField textField2;
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            
            textField = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField,pathRes,"Export localization table");
            textField.addMouseListener(mos); 
            
            textField2 = texts.get(t++); 
            mos=new MouseOptionSave(textField2,pathResex,"Export localization table");
            textField2.addMouseListener(mos); 
            
            
        }
        
        gd.showDialog();

        if (!gd.wasCanceled()){
            
            maximumDistance = (double)gd.getNextNumber();
            orderFit = (int)gd.getNextNumber();
            binSize = (double)gd.getNextNumber();
            pathRes = gd.getNextString(); 
            pathResex = gd.getNextString(); 
        }
        else{
            return;
        }
        
        double [] xp;
        double [] yp;
        double [] zp;
        ImagePlus imp = IJ.getImage();
        RoiManager rm = RoiManager.getInstance();
        Roi [] roi;
        if (rm!=null){
            try{
                roi = rm.getRoisAsArray();
            }
            catch(Exception ee){
                IJ.log("are you sure you added points in RoiManager tool");
                imp.close();return;
            }
        }
        else{
            IJ.log("are you sure you added points in RoiManager tool");
            imp.close();return;
        }
        
        
        try{
            if (roi.length<=0){
                IJ.log("please.. dont forget to select seeds");
                imp.close();
                return;
            }
            int nbPoints=0;
            for (int j=0; j<roi.length; j++) {
                FloatPolygon p = roi[j].getInterpolatedPolygon();
                nbPoints+=p.npoints;
            }
            IJ.log("nbPOints "+nbPoints);
            xp =new double[nbPoints];
            yp =new double[nbPoints];
            zp =new double[nbPoints];
            
            int t=0;
            for (int j=0; j<roi.length; j++) {
               
               
               FloatPolygon p = roi[j].getInterpolatedPolygon();
               for (int i=0;i<p.npoints;i++){
                   zp[t] = (double)roi[j].getPosition()-1;
                   xp[t] = (double)p.xpoints[i];
                   yp[t] = (double)p.ypoints[i];
                   //IJ.log(xp[t]+", "+yp[t]+", "+zp[t]);
                   t++;
               }
               
               
            } 
            //yp=imp.getRoi().getFloatPolygon().
            
        }
        catch (Exception ee){
            IJ.log("No measurment. You should add seeds on beads");
            imp.close();
            return;
        }
        imp.close();
        org.pasteur.imagej.development.FilamentResolution fr=new org.pasteur.imagej.development.FilamentResolution(sl,sizePix,xp,yp,zp,orderFit,maximumDistance,binSize,pathRes,pathResex);
        
            
    }
    
    
    
    
    void pointResolution(){
        
        
        
        
        double sizePix=20;
        
        
        
        
        String pathResinit=path_localization; 
        
        int nn=(pathResinit.lastIndexOf("."));
        if ((nn>pathResinit.length()-6)&&(nn>=0)){
            pathResinit=pathResinit.substring(0,nn);
        }
        int num=0;
        String pathRes=pathResinit+"_withFilament_"+num+".csv";
        String pathResex=pathResinit+"_withoutFilament_"+num+".csv";
        File f=new File(pathRes);
        File ff=new File(pathResex);
        while (f.exists()||ff.exists()){
            num++;
            pathRes=pathResinit+"_withFilament_"+num+".csv";
            f=new File(pathRes);
            pathResex=pathResinit+"_withoutFilament_"+num+".csv";
            ff=new File(pathResex);
        }
        
        
        
        if (sl!=null){
            ZRendering.hist3D(sl, sizePix,0);
            
        }
        else{
            IJ.log("please import a localization table first");
            return;
        }
        
        new WaitForUserDialog("Point selection", "Select a point and add it in RoiManager (ctrl+t) / Then press OK.").show();
        
        
        

        GenericDialog gd = new GenericDialog("ZOLA: Resolution computation");

        
        gd.addNumericField("maximum point radius (nm): ", maximumDistance,1);
        gd.addNumericField("histogram bin size: ", binSize,1);
        
        gd.addStringField("Result localization sub region :",pathRes,sizeTextString); 
        gd.addStringField("Result localization table without filament :",pathResex,sizeTextString); 
            
        TextField textField; 
        TextField textField2; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            
            
            textField = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField,pathRes,"Export localization table");
            textField.addMouseListener(mos); 
            
            textField2 = texts.get(t++); 
            mos=new MouseOptionSave(textField2,pathResex,"Export localization table");
            textField2.addMouseListener(mos); 
            
            
        }
            
        
        gd.showDialog();

        if (!gd.wasCanceled()){
            
            maximumDistance = (double)gd.getNextNumber();
            
            binSize = (double)gd.getNextNumber();
            pathRes = gd.getNextString(); 
            pathResex = gd.getNextString(); 
        }
        else{
            return;
        }
        
        double [] xp;
        double [] yp;
        double [] zp;
        ImagePlus imp = IJ.getImage();
        RoiManager rm = RoiManager.getInstance();
        Roi [] roi;
        if (rm!=null){
            try{
                roi = rm.getRoisAsArray();
            }
            catch(Exception ee){
                IJ.log("are you sure you added points in RoiManager tool");
                imp.close();return;
            }
        }
        else{
            IJ.log("are you sure you added points in RoiManager tool");
            imp.close();return;
        }
        
        
        try{
            if (roi.length<=0){
                IJ.log("please.. dont forget to select seeds");
                imp.close();
                return;
            }
            int nbPoints=0;
            for (int j=0; j<roi.length; j++) {
                FloatPolygon p = roi[j].getFloatPolygon();
                nbPoints+=p.npoints;
            }
            IJ.log("nbPOints "+nbPoints);
            xp =new double[nbPoints];
            yp =new double[nbPoints];
            zp =new double[nbPoints];
            
            int t=0;
            for (int j=0; j<roi.length; j++) {
               
               
               FloatPolygon p = roi[j].getFloatPolygon();
               for (int i=0;i<p.npoints;i++){
                   zp[t] = (double)roi[j].getPosition()-1;
                   xp[t] = (double)p.xpoints[i];
                   yp[t] = (double)p.ypoints[i];
                   //IJ.log(xp[t]+", "+yp[t]+", "+zp[t]);
                   t++;
               }
               
               
            } 
            //yp=imp.getRoi().getFloatPolygon().
            
        }
        catch (Exception ee){
            IJ.log("No measurment. You should add seeds on beads");
            imp.close();
            return;
        }
        imp.close();
        
        org.pasteur.imagej.development.PointResolution pr=new org.pasteur.imagej.development.PointResolution(sl,sizePix,xp,yp,zp,maximumDistance,binSize,pathRes,pathResex);
        
            
    }
    
    
    
    
    
    
    void addPoissonNoise(){
        
        ImagePlus imp=IJ.getImage();


        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }

        processing=new Boolean(true);


        if (imp==null){
            IJ.log("no image opened");
            processing=new Boolean(false);
        }
        /*else if ((imp.getWidth()>180)||(imp.getHeight()>180)){
            IJ.log("image should be < 180 pixels for the moment");
        }*/
        else{


            Roi r=imp.getRoi();
            int width=imp.getWidth();
            int height=imp.getHeight();
            int startX=0;
            int startY=0;
            if (r!=null){
                Rectangle bound=imp.getRoi().getBounds();
                startX=bound.x;
                startY=bound.y;
                width=bound.width;
                height=bound.height;
            }
            
            ImageStack ims = new ImageStack(width,height);
            for (int i=0;i<imp.getNSlices();i++){
                for (int ii=0;ii<imp.getNFrames();ii++){
                    imp.setSlice(i*imp.getNFrames()+ii+1);
                    ImageProcessor p=imp.getProcessor();
                    ImageProcessor p2 = new FloatProcessor(width,height);
                    for (int u=startX;u<width;u++){
                        for (int uu=startY;uu<height;uu++){
                            double pixVal = p.getPixelValue(u, uu);
                            java.util.Random rand = new java.util.Random();
                            double L = (-(pixVal));
                            int k = 0;
                            double pp = 0;
                            do {
                               k++;
                               // Generate uniform random number u in [0,1] and let p ← p × u.
                               pp += Math.log(rand.nextDouble());
                            } while (pp >= L);
                            p2.putPixelValue(u, uu, k-1);
                        }
                    }
                    ims.addSlice(p2);
                }
            }
            ImagePlus imp2 = new ImagePlus("Poisson noise image",ims);
            imp2.show();

        }
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    void simulation(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
        double backSecOrder=-0.001;
        double photonNumber=2000;
        double background=30;
        int copies=1;
        double stepXloc=0;
        double stepYloc=0;
        
        
        
        GenericDialog gd = new GenericDialog("Simulation");
        
        gd.addStringField("Calibration file:",path_calibration,sizeTextString); 
        
        
        double refIndex=1.33;
        IJ.log("WARNING: REFRACTIVE INDEX SIMULATION IMMERSION MEDIUM="+refIndex);
        
        gd.addMessage("Computation parameter");

        
        gd.addNumericField("patch size: ", sizePatchPhaseRet,0);

        gd.addMessage("Simulation parameters");
        
        gd.addNumericField("distance_focus_to_coverslip (µm): ", zfocus, 3);
        
        
        gd.addNumericField("Axial_range: ", axialRange,2,6,"(µm)");
        
        gd.addNumericField("z_step (µm): ", stepZloc,4);
        gd.addNumericField("x_step (µm): ", stepXloc,4);
        gd.addNumericField("y_step (µm): ", stepYloc,4);
        
        //gd.addNumericField("sample ref index: ", refIndex,2);

        gd.addNumericField("Photon_number: ", photonNumber,0);
        gd.addNumericField("Background_intensity: ", background,0);
        
        gd.addNumericField("Background_slope: ", backSecOrder,3);
        
        
        gd.addNumericField("copy number at each position: ", copies,0);
        
        
        TextField textField;  
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
            
            
            
        }
        
        
        gd.showDialog();



        if (!gd.wasCanceled()){
            
            path_calibration = gd.getNextString(); 
            
            sizePatchPhaseRet = (int)gd.getNextNumber();
            
            zfocus = (double)gd.getNextNumber();
            
            axialRange= (double)gd.getNextNumber();
            
            stepZloc=(double)gd.getNextNumber();
            
            stepXloc=(double)gd.getNextNumber();
            stepYloc=(double)gd.getNextNumber();
            
            //refIndex=(double)gd.getNextNumber();
            
            
            photonNumber=(double)gd.getNextNumber();
            background=(double)gd.getNextNumber();
            
            backSecOrder=(double)gd.getNextNumber();
            IJ.log(" "+backSecOrder);
            copies=(int)gd.getNextNumber();
            
            

        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        if (copies<1){
            IJ.log("Process stopped: the number of copies should be >= 1");
            processing=new Boolean(false);
            return;
        }
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        MyCudaStream.init(1);

        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
        if (!dp.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        if (dp.param!=null){
            dp.setSizeoutput(sizePatchPhaseRet);
            dp.setNwat(refIndex);
            dp.param.Zfocus=zfocus;
            //warning: we should use Zfocus and use SearchPSFcenter class
            org.pasteur.imagej.development.Simulation sim =new org.pasteur.imagej.development.Simulation(copies, dp,axialRange,stepZloc,stepXloc,stepYloc,photonNumber,background,backSecOrder,sizePatchPhaseRet); 
            
            sim.run();
            
            dp.free();
        }
        
        MyCudaStream.destroy();
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    void simulation2beads(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
        double backSecOrder=-0.001;
        double photonNumber=2000;
        double background=30;
        double stepXloc=0;
        double stepYloc=0;
        
        
        
        
        GenericDialog gd = new GenericDialog("Simulation 2");
        
        gd.addStringField("Calibration file:",path_calibration,sizeTextString); 
        
        
        double refIndex=1.33;
        IJ.log("WARNING: REFRACTIVE INDEX SIMULATION IMMERSION MEDIUM="+refIndex);
        
        gd.addMessage("Computation parameter");

        
        gd.addNumericField("patch size: ", sizePatchPhaseRet,0);

        gd.addMessage("Simulation parameters");
        
        gd.addNumericField("distance_focus_to_coverslip (µm): ", zfocus, 3);
        
        
        gd.addNumericField("Axial_range: ", axialRange,2,6,"(µm)");
        
        
        //gd.addNumericField("sample ref index: ", refIndex,2);

        gd.addNumericField("Photon_number: ", photonNumber,0);
        gd.addNumericField("Background_intensity: ", background,0);
        
        
        
        
        
        TextField textField;  
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
            
            
            
        }
        
        
        gd.showDialog();



        if (!gd.wasCanceled()){
            
            path_calibration = gd.getNextString(); 
            
            sizePatchPhaseRet = (int)gd.getNextNumber();
            
            zfocus = (double)gd.getNextNumber();
            
            
            axialRange= (double)gd.getNextNumber();
            
            
            //refIndex=(double)gd.getNextNumber();
            
            
            photonNumber=(double)gd.getNextNumber();
            background=(double)gd.getNextNumber();
            
            
            
            
            

        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        MyCudaStream.init(1);

        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
        if (!dp.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        if (dp.param!=null){
            dp.setSizeoutput(sizePatchPhaseRet);
            dp.setNwat(refIndex);
            dp.param.Zfocus=zfocus;
            //warning: we should use Zfocus and use SearchPSFcenter class
            org.pasteur.imagej.development.Simulation2Obj sim =new org.pasteur.imagej.development.Simulation2Obj( dp,axialRange,photonNumber,background,sizePatchPhaseRet); 
            
            sim.run();
            
            dp.free();
        }
        
        MyCudaStream.destroy();
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    void simulationOverlappingBeads(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
        double backSecOrder=-0.001;
        double photonNumber=2000;
        double background=30;
        double stepXloc=0;
        double stepYloc=0;
        
        double distance=1;//µm
        
        
        GenericDialog gd = new GenericDialog("Simulation Overlapping beads");
        
        gd.addStringField("Calibration file:",path_calibration,sizeTextString); 
        
        
        double refIndex=1.33;
        IJ.log("WARNING: REFRACTIVE INDEX SIMULATION IMMERSION MEDIUM="+refIndex);
        
        gd.addMessage("Computation parameter");

        
        gd.addNumericField("patch size: ", sizePatchPhaseRet,0);

        gd.addMessage("Simulation parameters");
        
        gd.addNumericField("distance_focus_to_coverslip (µm): ", zfocus, 3);
        
        
        gd.addNumericField("Axial_range: ", axialRange,2,6,"(µm)");
        
        
        //gd.addNumericField("sample ref index: ", refIndex,2);

        gd.addNumericField("Photon_number: ", photonNumber,0);
        gd.addNumericField("Background_intensity: ", background,0);
        
        
        gd.addNumericField("max_distance_overlapping: ", distance,2,6,"(µm)");
        
        
        TextField textField;  
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
            
            
            
        }
        
        
        gd.showDialog();



        if (!gd.wasCanceled()){
            
            path_calibration = gd.getNextString(); 
            
            sizePatchPhaseRet = (int)gd.getNextNumber();
            
            zfocus = (double)gd.getNextNumber();
            
            
            axialRange= (double)gd.getNextNumber();
            
            
            //refIndex=(double)gd.getNextNumber();
            
            
            photonNumber=(double)gd.getNextNumber();
            background=(double)gd.getNextNumber();
            
            
            distance=(double)gd.getNextNumber();
            
            

        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        MyCudaStream.init(1);

        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
        if (!dp.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        if (dp.param!=null){
            dp.setSizeoutput(sizePatchPhaseRet);
            dp.setNwat(refIndex);
            dp.param.Zfocus=zfocus;
            //warning: we should use Zfocus and use SearchPSFcenter class
            org.pasteur.imagej.development.SimulationOverlappingBeads sim =new org.pasteur.imagej.development.SimulationOverlappingBeads( dp,axialRange,photonNumber,background,sizePatchPhaseRet,distance); 
            
            sim.run();
            
            dp.free();
        }
        
        MyCudaStream.destroy();
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    void crlbfromfile(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
            

        GenericDialog gd = new GenericDialog("Fundamental limit (CRLB)");
        
        String pathdir=IJ.getDirectory("image");
        if (pathdir==null){
            pathdir=IJ.getDirectory("current");
        }
        
        String pathResult=pathdir+"crlbtable.csv";
        
        
        
        String path_calibration2="";
        
        gd.addStringField("Calibration file:",path_calibration,sizeTextString); 
        
        gd.addStringField("Calibration_file_cam2 [optional]:",path_calibration2,sizeTextString); 
        
        
        gd.addMessage("Computation parameter");

        
        gd.addNumericField("size of patch: ", sizePatchPhaseRet,0);
        
        gd.addNumericField("axial_range (µm): ", axialRange,2);
        gd.addNumericField("distance_focus_to_coverslip (µm): ", zfocus,2);
        gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        gd.addMessage("CRLB parameters");
        
        
        gd.addStringField("Result table:",pathResult,sizeTextString); 
        
        TextField textField; 
        TextField textField2,textField3; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
            textField3 = texts.get(t++); 
            mol=new MouseOptionLoad(textField3,path_calibration,"Import cal. file cam 2 (optional)");
            textField3.addMouseListener(mol); 
            
            textField2 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField2,path_calibration,"Export CRLB result file");
            textField2.addMouseListener(mos); 
            
            
        }
        

        gd.showDialog();



        if (!gd.wasCanceled()){
            
            
            
            path_calibration = gd.getNextString(); 
            path_calibration2 = gd.getNextString(); 

            sizePatchPhaseRet = (int)gd.getNextNumber();
            
            axialRange = (double)gd.getNextNumber();
            zfocus = (double)gd.getNextNumber();
            nwat = (double)gd.getNextNumber();
            
            pathResult = gd.getNextString(); 

        }
        else{
             processing=new Boolean(false);
            return;
        }
        
        
        if (pathResult.equals(path_calibration)){
            IJ.log("Process stopped: the result crlb path has to be different than calibration path");
            processing=new Boolean(false);
            return;
        }
        
        
        if (pathResult.equals(path_calibration2)){
            if (pathResult.length()>0){//only if path not null
                IJ.log("Process stopped: the result crlb path has to be different than calibration path");
                processing=new Boolean(false);
                return;
            }
        }
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        
        org.pasteur.imagej.cuda.MyCudaStream.init(1);

        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
        if (!dp.loading){
            
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        org.pasteur.imagej.process.gpu.DataPhase_ dp2=null;
        if (path_calibration2.length()>3){
            dp2 = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration2);
            if (!dp2.loading){
                IJ.log("impossible to load "+path_calibration2);
                processing=new Boolean(false);
                return;
            }
        }
        
        dp.setSizeoutput(sizePatchPhaseRet);
        dp.param.Zfocus=zfocus;
        dp.setNwat(nwat);
        org.pasteur.imagej.process.gpu.CRLB_ crlb; 
        String path=FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B");
//        if (dp2!=null){
//            String path2=classes.FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B for cam 2");
//            dp2.setSizeoutput(sizePatchPhaseRet);
//            crlb =new CRLB(dp,dp2,path,path2); 
//        }
//        else{
            crlb =new org.pasteur.imagej.process.gpu.CRLB_(dp,path,axialRange); 
//        }
        
        
        
        
        
        
        crlb.run(pathResult);
        
        dp.free();
        
        org.pasteur.imagej.cuda.MyCudaStream.destroy();
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    void crlbfromfileDualObj(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
            
        


        GenericDialog gd = new GenericDialog("Fundamental limit (CRLB)");
        
        String pathdir=IJ.getDirectory("image");
        if (pathdir==null){
            pathdir=IJ.getDirectory("current");
        }
        
        String pathResult=pathdir+"crlbtable.csv";
        
        double zfocus2=zfocus;
        
        String path_calibration2="";
        
        
        gd.addMessage("-------------Camera 1-------------");
        
        gd.addStringField("Calibration_file_A1:",path_calibration,sizeTextString); 
        
        //gd.addStringField("Calibration_file_cam2 [optional]:",path_calibration2,sizeTextString); 
        
        gd.addNumericField("distance_focus_to_coverslip_A1 (µm)", zfocus, 3);
        
        
        gd.addMessage("-------------Camera 2-------------");
        
        gd.addStringField("Calibration_file_A2:",path_calibration,sizeTextString); 
        
        //gd.addStringField("Calibration_file_cam2 [optional]:",path_calibration2,sizeTextString); 
        
        gd.addNumericField("distance_focus_to_coverslip_A2 (µm)", zfocus2, 3);
        
        
        
        gd.addMessage("Computation parameter");

        
        gd.addNumericField("size of patch: ", sizePatchPhaseRet,0);
        
        gd.addNumericField("axial_range (µm): ", axialRange,2);
        gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        gd.addMessage("CRLB parameters");
        
        
        gd.addStringField("Result table:",pathResult,sizeTextString); 
        
        TextField textField; 
        TextField textField2,textField3; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;
            
            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
            textField.addMouseListener(mol); 
            
            textField3 = texts.get(t++); 
            mol=new MouseOptionLoad(textField3,path_calibration,"Import cal. file cam 2 (optional)");
            textField3.addMouseListener(mol); 
            
            textField2 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField2,path_calibration,"Export CRLB result file");
            textField2.addMouseListener(mos); 
            
            
        }
        

        gd.showDialog();



        if (!gd.wasCanceled()){
            
            
            
            path_calibration = gd.getNextString(); 
            zfocus = (double)gd.getNextNumber();
            
            path_calibration2 = gd.getNextString(); 
            zfocus2 = (double)gd.getNextNumber();
            
            sizePatchPhaseRet = (int)gd.getNextNumber();
            
            axialRange = (double)gd.getNextNumber();
            
            nwat = (double)gd.getNextNumber();
            
            pathResult = gd.getNextString(); 

        }
        else{
             processing=new Boolean(false);
            return;
        }
        
        
        if (pathResult.equals(path_calibration)){
            IJ.log("Process stopped: the result crlb path has to be different than calibration path");
            processing=new Boolean(false);
            return;
        }
        
        
        if (pathResult.equals(path_calibration2)){
            if (pathResult.length()>0){//only if path not null
                IJ.log("Process stopped: the result crlb path has to be different than calibration path");
                processing=new Boolean(false);
                return;
            }
        }
        
        int sizeFFT=128;
        sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
        sizeFFT=Math.min(sizeFFT,maxFFT);
        

        org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
        if (!dp.loading){
            
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        org.pasteur.imagej.process.cpu.DataPhase dp2=null;
        if (path_calibration2.length()>3){
            dp2 = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration2);
            if (!dp2.loading){
                IJ.log("impossible to load "+path_calibration2);
                processing=new Boolean(false);
                return;
            }
        }
        
        org.pasteur.imagej.process.cpu.CRLBdualCam crlb; 
        
        
        dp.setSizeoutput(sizePatchPhaseRet);
        dp.param.Zfocus=zfocus;
        dp.setNwat(nwat);
        String path1=FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B (cam 1)");
        
        
        
        dp2.setSizeoutput(sizePatchPhaseRet);
        dp2.param.Zfocus=zfocus2;
        dp2.setNwat(nwat);
        String path2=FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B (cam 2)");
        
        
//        if (dp2!=null){
//            String path2=classes.FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B for cam 2");
//            dp2.setSizeoutput(sizePatchPhaseRet);
//            crlb =new CRLB(dp,dp2,path,path2); 
//        }
//        else{
            crlb =new org.pasteur.imagej.process.cpu.CRLBdualCam(dp,dp2,path1,path2,axialRange); 
//        }
        
        
        
        
        
        
        crlb.run(pathResult);
        
        
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    void multiEmitterLocalization(){
        
        
        nbStream=1;
            
        ImagePlus imp=IJ.getImage();


        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
            processing=new Boolean(true);
            
            
            if (imp==null){
                IJ.log("no image opened");
                processing=new Boolean(false);
            }
            /*else if ((imp.getWidth()>180)||(imp.getHeight()>180)){
                IJ.log("image should be < 180 pixels for the moment");
            }*/
            else{
                
                
                
                Roi r=imp.getRoi();
                int width=imp.getWidth();
                int height=imp.getHeight();
                
                boolean scmosLoad=false;
                if (isSCMOS){
                    scmosLoad=loadSCMOScameraFiles(width, height, r);
                }
                
                
                
                String pathdir=IJ.getDirectory("image");
                if (pathdir==null){
                    pathdir=IJ.getDirectory("current");
                }
                String path=pathdir;
                
                path=path+"ZOLA_localization_table.csv";
                
                
                GenericDialog gd = new GenericDialog("ZOLA: Localization");
                
                
                Font font = gd.getFont();
                Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
                
                gd.addCheckbox("Run_on_GPU", GPU_computation);
                
                
                
                
                
                
                if (!isSCMOS){
                    gd.addMessage("Camera is EMCCD",fontBold);
                    gd.addMessage("ADU = "+adu);
                    gd.addNumericField("Gain = ", gain, 1,6,"(1 if camera provides photon counts)");
                    gd.addMessage("Offset = "+offset);
                    
                }
                else{
                    if (scmosLoad){
                        gd.addMessage("Camera is SCMOS",fontBold);
                    }
                    else{
                        gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                    }
//                    String sep=File.separator;
//                    String [] var=this.path_SCMOSvariance.split(sep);
//                    String [] off=this.path_SCMOSoffset.split(sep);
//                    String [] gain=this.path_SCMOSgain.split(sep);
//                    gd.addMessage("Offset file = "+off[var.length-1]);
//                    gd.addMessage("Variance file = "+var[var.length-1]);
//                    gd.addMessage("Gain file = "+gain[var.length-1]);
                    
                    gd.addMessage("Offset file = "+path_SCMOSoffset);
                    gd.addMessage("Variance file = "+path_SCMOSvariance);
                    gd.addMessage("Gain file = "+path_SCMOSgain);
                    
                }
                
                
                
                gd.addStringField("PSF_calibration_file:",path_calibration,sizeTextString); 
                
                
                
                gd.addMessage("Sample parameters", fontBold);
                
                gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
        
                gd.addNumericField("Distance_focus_to_coverslip", zfocus, 3,6,"(µm)");
                
                
                
                gd.addMessage("Advanced parameter", fontBold);
                
                
                gd.addNumericField("Patch_size: ", sizePatch,0,6,"(pixels)");
                
                //gd.addNumericField("number_of_streams: ", nbStream,0);
                //gd.addNumericField("number_of_threads: ", nbThread,0);
                
                
                
                
                
                gd.addNumericField("Expected_axial_range: ", axialRange,2,6,"(µm)");
                        
                
                gd.addNumericField("Min_number_of_photons: ", photonThreshold,0);
                
                
                gd.addStringField("Localization_table:",path,sizeTextString); 
                
        
                TextField textField; 
            TextField textField2; 
            // Add a mouse listener to the config file field 
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
            {
                int t=0;
                
                Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
                
                
                textField = texts.get(t++); 
                MouseOptionLoad mol=new MouseOptionLoad(textField,path_calibration,"Import calibration file");
                textField.addMouseListener(mol); 
                
                

                textField2 = texts.get(t++); 
                MouseOptionSave mos=new MouseOptionSave(textField2,pathdir,"Export localization result table");
                textField2.addMouseListener(mos); 
                
                


            }
        
        
                gd.showDialog();
                
                
                
                
                
                if (!gd.wasCanceled()){
                    GPU_computation = (boolean)gd.getNextBoolean();

                    gain = (double)gd.getNextNumber(); 
                    
                    path_calibration = gd.getNextString(); 
                    
                    nwat= (double)gd.getNextNumber();
                    
                    zfocus= (double)gd.getNextNumber();
                    
                    sizePatch = (int)gd.getNextNumber();
                    //nbStream = (int)gd.getNextNumber();
                    //nbThread = (int)gd.getNextNumber();
                    
                    axialRange= (double)gd.getNextNumber();
                    
                    
                            
                    
                    photonThreshold = (int)gd.getNextNumber();
                    
                    path_localization = gd.getNextString(); 
                }
                else{
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                int nnn=(path_localization.lastIndexOf("."));
                if ((nnn>=0)&&(nnn>path_localization.length()-6)){
                    path_localization=path_localization.substring(0,nnn);
                }
                path_localization+=".csv";
            
                if (path_localization.equals(path_calibration)){
                    IJ.log("Process stopped: the localization path has to be different than calibration path");
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                if (sizePatch%2==1){
                    sizePatch++;
                }
                
                
                
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                }
                if ((width<sizePatch)||(height<sizePatch)){
                    IJ.log("Process stopped: Localization impossible. The image, or selected ROI should be larger than patch size");
                    IJ.log("Please, select a large rectangle in the image before running the localization");
                    IJ.log("image width:"+width+" ; height:"+height+" ; sizePatch:"+sizePatch);
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                
                
                
                if (gain==0||adu==0){
                    new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
                }
                else{
                    
                    
                    
                    
                    int sizeFFT=128;
                    sizeFFT=Math.max(sizeFFT,sizePatch*2);
                    sizeFFT=Math.min(sizeFFT,maxFFT);
                    
                    
                    if (GPU_computation){
                        MyCudaStream.init(nbStream+1);
                        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(nwat);
                        dp.param.Zfocus=zfocus;
                        IJ.log("Localization started");
                        
                        
                        org.pasteur.imagej.process.gpu.LocalizationMotionBlurPipeline_ lmbp=null;
                        if (isSCMOS){
                            IJ.log("SCMOS not yet implemented");
                        }
                        else{
                            lmbp=new org.pasteur.imagej.process.gpu.LocalizationMotionBlurPipeline_(dp,axialRange,photonThreshold,nbStream,sizePatch,path_localization,adu, gain,offset,true);
                        }
                        
                        sl=new StackLocalization();
                        sl=lmbp.detectParticles(imp,sl);


                        
                        
                        dp.free();
                        
                        
                        MyCudaStream.destroy();
                    }
                    else{
                        
                        IJ.log("CPU not yet implemented");
                        
                    }
                    
                }
            }
        
        processing=new Boolean(false);
        
    }
    
    
    
    
    void computeChromaticAberrations(){
        
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        File f = new File(path_localization);
        String path=f.getParent();
        int nn=(path.lastIndexOf("."));
        if ((nn>path.length()-6)&&(nn>=0)){
            path=path.substring(0,nn);
        }
        
        path=path+File.separator+"ZOLA_registrationParameters.csv";
        
            
        
        
        
        int orderFit=2;
        int minFrameBead=10000;
        
        processing=new Boolean(true);


        /*String [] imTitle=WindowManager.getImageTitles();

        if (imTitle.length!=2){
            IJ.log("Please, open 2 images first");
            processing=new Boolean(false);
        }*/


        
        



        GenericDialog gd = new GenericDialog("2 colors registration");
        
        gd.addMessage("At this stage, we assume the drift has been corrected");

        gd.addMessage("Input localization tables");

        gd.addStringField("Color_1 (reference):",path_localization,sizeTextString); 
        
        

        gd.addStringField("Color_2:",path_localization2,sizeTextString); 

        
        //gd.addNumericField("if flip -> pixel number (X):", imSize,0);
        //gd.addNumericField("if flip -> pixel size (µm):", this.xystep,3);
        
        gd.addMessage("Fusion parameters:");
        
        gd.addNumericField("max._X/Y distance between cam1 and cam2 localizations (nm): ", dualCamMaxDistanceMergingPlan,1);
        
        
        
        gd.addNumericField("polynomial fit order: ", orderFit,1);
        
        gd.addMessage("Output registration parameters");

        gd.addStringField("path_registration_parameters:",path,sizeTextString); 
        
gd.addMessage("WARNING !! if order>0 -> size of images used to compute registration parameters");
gd.addMessage("should have the same size as image to register !!");
        
        
        
        TextField textField; 

        TextField textField2; 
        TextField textField3; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;

            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_localization,"Import localization table (cam 1)");
            textField.addMouseListener(mol); 

            textField2 = texts.get(t++); 
            mol=new MouseOptionLoad(textField2,path_localization2,"Import localization table (cam 2)");
            textField2.addMouseListener(mol); 
            
            
            textField3 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField3,path,"Save registred localization table (cam 2)");
            textField3.addMouseListener(mos); 
            

            

        }


        gd.showDialog();



        if (!gd.wasCanceled()){


            path_localization = gd.getNextString(); 
            
            path_localization2 = gd.getNextString(); 

            //imSize=(int)gd.getNextNumber();
            //xystep=gd.getNextNumber();
            
            dualCamMaxDistanceMergingPlan=gd.getNextNumber();
            
            
            
            
            orderFit=(int)gd.getNextNumber();
            
            
            path = gd.getNextString(); 
            
        }
        else{
            processing=new Boolean(false);
            return;
        }
        path_registration=path;
        
        
        if ((path_localization.equals(path_localization2))||(path_localization.equals(path))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((path_localization2.equals(path))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        
        
        
        StackLocalization sl1=new StackLocalization(path_localization);
        StackLocalization sl2=new StackLocalization(path_localization2);
        
        //double imSizeNM=((double)imSize)*xystep*1000;
        
        
        
        
        DualColorFusion dcf=new DualColorFusion(sl1,sl2,path,dualCamMaxDistanceMergingPlan,orderFit);
        dcf.run();
        IJ.log("Registration parameters saved");
        processing=new Boolean(false);
        
        
        
    }
    
    
    
    
    
    
    
    void computeSR2LRregistration(){
        
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        File f = new File(path_localization);
        String path=f.getParent();
        int nn=(path.lastIndexOf("."));
        if ((nn>path.length()-6)&&(nn>=0)){
            path=path.substring(0,nn);
        }
        
        path=path+File.separator+"ZOLA_localization_table_registered.csv";
        
            
        String path_widefield2="";
        String path_widefield3="";
        int magnification=1;
        
        double sigmaNMxy=100;
        double sigmaNMz=200;
        
        processing=new Boolean(true);
        
        
        String path_widefield=IJ.getDirectory("current");



        GenericDialog gd = new GenericDialog("2 colors registration");
        
        gd.addMessage("At this stage, we assume the drift has been corrected");


        gd.addStringField("Widefield_Image (reference):",path_widefield,sizeTextString); 
        gd.addStringField("Widefield_Image_green (optional):",path_widefield2,sizeTextString); 
        gd.addStringField("Widefield_Image_blue (optional):",path_widefield3,sizeTextString); 
        
        gd.addNumericField("Widefield_pixel_size:", xystep*1000, 1,6,"(nm)");
        gd.addNumericField("Widefield_Z-step:", zstep*1000, 1,6,"(nm)");

        gd.addStringField("Localization table:",path_localization2,sizeTextString); 
        
        gd.addNumericField("Magnification:", magnification, 0,6,"");
        
        gd.addNumericField("Sigma_XY:", sigmaNMxy, 1,6,"(nm)");
        gd.addNumericField("Sigma_Z:", sigmaNMz, 1,6,"(nm)");
        
        
        gd.addMessage("Output registration parameters");

        gd.addStringField("path_registration_parameters:",path,sizeTextString); 
        
        
        
        TextField textField; 

        TextField textField1; 
        TextField textField2; 
        TextField textField3;
        TextField textField4; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;

            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_widefield,"Import widefield image 1");
            textField.addMouseListener(mol); 
            
            textField1 = texts.get(t++); 
            mol=new MouseOptionLoad(textField1,path_widefield2,"Import widefield image 2");
            textField1.addMouseListener(mol); 
            
            textField4 = texts.get(t++); 
            mol=new MouseOptionLoad(textField4,path_widefield3,"Import widefield image 3");
            textField4.addMouseListener(mol); 
            
            textField2 = texts.get(t++); 
            mol=new MouseOptionLoad(textField2,path_localization,"Import localization table");
            textField2.addMouseListener(mol); 
            
            
            textField3 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField3,path,"Save registered localization table");
            textField3.addMouseListener(mos); 
            

            

        }


        gd.showDialog();


        IJ.log("zstep"+zstep);
        if (!gd.wasCanceled()){


            path_widefield = gd.getNextString(); 
            path_widefield2 = gd.getNextString(); 
            path_widefield3 = gd.getNextString();
            
            xystep = (double)gd.getNextNumber()/1000;
            zstep = (double)gd.getNextNumber()/1000;
            
            
            path_localization = gd.getNextString(); 

            
            magnification = (int)gd.getNextNumber();
            sigmaNMxy = (double)gd.getNextNumber();
            sigmaNMz = (double)gd.getNextNumber();
            path = gd.getNextString(); 
            
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        
        if ((path_localization.equals(path))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        
        
        
        StackLocalization sl1=new StackLocalization(path_localization);
        
        //double imSizeNM=((double)imSize)*xystep*1000;
        
        ImagePlus imp = new ImagePlus(path_widefield);
        ImagePlus imp2=null;
        if (path_widefield2.length()>3){
            try{
                imp2 = new ImagePlus(path_widefield2);
            }catch(Exception e){}
        }
        ImagePlus imp3=null;
        if (path_widefield3.length()>3){
            try{
                imp3 = new ImagePlus(path_widefield3);
            }catch(Exception e){}
        }
        DualColorFusion dcf=new DualColorFusion(imp,sl1,path,xystep*1000,zstep*1000,magnification,sigmaNMxy,sigmaNMz,imp2,imp3);
        dcf.run();
        IJ.log("Registration parameters saved");
        processing=new Boolean(false);
        
        
        
    }
    
    
    
    
    
    
    
    void colorRegistration(){
        
        
        
        
        
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        
        String path=path_localization2;
        
        int nn=(path.lastIndexOf("."));
        
        nn=(path.lastIndexOf("."));
        if ((nn>path.length()-6)&&(nn>=0)){
            path=path.substring(0,nn);
        }
        path=path+"_registered.csv";


        GenericDialog gd = new GenericDialog("Second color registration");
        

        gd.addMessage("Input localization tables");

        gd.addStringField("Localization_Table_to_register:",path_localization2,sizeTextString); 
        
        gd.addMessage("Input registration parameters");

        gd.addStringField("registration_parameter_file:",path_registration,sizeTextString); 

        gd.addMessage("Output registered localization table");
        
        gd.addStringField("Localization_Table_registered:",path,sizeTextString); 
        
        TextField textField; 

        TextField textField2; 
        TextField textField3; 
        // Add a mouse listener to the config file field 
        if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())) 
        {
            int t=0;

            Vector<TextField> texts = (Vector<TextField>) gd.getStringFields(); 
            textField = texts.get(t++); 
            MouseOptionLoad mol=new MouseOptionLoad(textField,path_localization2,"Import localization table");
            textField.addMouseListener(mol); 

            textField2 = texts.get(t++); 
            mol=new MouseOptionLoad(textField2,path_registration,"Import registration parameters");
            textField2.addMouseListener(mol); 
            
            
            textField3 = texts.get(t++); 
            MouseOptionSave mos=new MouseOptionSave(textField3,path,"Save registred localization table");
            textField3.addMouseListener(mos); 
            

            

        }


        gd.showDialog();



        if (!gd.wasCanceled()){


            path_localization2 = gd.getNextString(); 
            
            path_registration = gd.getNextString(); 

            
            path=gd.getNextString();
            
            
            
        }
        else{
            processing=new Boolean(false);
            return;
        }
        
        
        if ((path_localization2.equals(path_registration))||(path_localization2.equals(path))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        else if ((path_registration.equals(path))){
            IJ.log("Process stopped: at least 2 paths have the same name");
            processing=new Boolean(false);
            return;
        }
        
        
        
        StackLocalization sl1=new StackLocalization(path_localization2);
        
        //double imSizeNM=((double)imSize)*xystep*1000;
        
        
        
        
        DualColorFusion dcf=new DualColorFusion(sl1,path_registration,path);
        dcf.run();
        
        IJ.log("Registered table saved");
        processing=new Boolean(false);
        
        
        
    }
    
    
    
    
    void mergeframe(){
        
        
        
        
            
        ImagePlus imp=IJ.getImage();


        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
            processing=new Boolean(true);
            
            
            if (imp==null){
                IJ.log("no image opened");
                processing=new Boolean(false);
            }
            /*else if ((imp.getWidth()>180)||(imp.getHeight()>180)){
                IJ.log("image should be < 180 pixels for the moment");
            }*/
            else{
                
                
                Roi r=imp.getRoi();
                int width=imp.getWidth();
                int height=imp.getHeight();
                
                
                
                
                
                GenericDialog gd=null;
                
                //since NonBlockingGenericDialog incompatible with headless mode:
                if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())){
                    gd = new NonBlockingGenericDialog("ZOLA: Merge frames");
                }
                else{
                    gd = new GenericDialog("ZOLA: Merge frames");
                }

                    Font font = gd.getFont();
                    Font fontBold= gd.getFont();
                    try{
                        fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
                    }catch(Exception e){}

                    int merge=10;
                    
                    gd.addNumericField("merged frame number: ", merge, 3,6,"");
                    
                    gd.showDialog();





                    if (!gd.wasCanceled()){
                        
                        merge = (int)gd.getNextNumber();
                        
                    }
                    else{
                        processing=new Boolean(false);
                        return;
                    }
                
                int xstart=0;
                int ystart=0;
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                    xstart=bound.x;
                    ystart=bound.y;
                }
                
                
                
                
                
                boolean stacked=true;
                int nbImage=imp.getNSlices();
                if (imp.getNSlices()==1){
                    nbImage=imp.getNFrames();
                    stacked=false;
                }
                ImageStack [] ims=null;
                if (stacked){
                    ims = new ImageStack[1];
                    ims[0]=imp.getStack();

                }
            
            
            
                
                ImageProcessor ip;
                
                ArrayList<float [][]> al = new ArrayList<float [][]>();
                
                float[][] summed = new float[width][height];
                ImageStack imsres = new ImageStack(width,height);
                int index=0;
                for (int i=0;i<nbImage;i++){
                    IJ.showProgress((double)i/(double)(nbImage));
                    if (stacked){
                        ip = ims[0].getProcessor(i+1);
                    }
                    else{
                        imp.setSlice(i+1);
                        ip = imp.getProcessor();
                    }
                    
                    //fill tmp stack at the beginning
                    float[][] tmp = new float[width][height];
                    for (int u=0;u<width;u++){
                        for (int uu=0;uu<height;uu++){
                            float v=ip.getPixelValue(u+xstart, uu+ystart);
                            tmp[u][uu]=v;
                            summed[u][uu]+=v;
                        }
                    }
                    al.add(tmp);
                        
                    
                    if (i>=merge/2&&i<=merge){
                        FloatProcessor ipr = new FloatProcessor(width,height);
                        for (int u=0;u<width;u++){
                            for (int uu=0;uu<height;uu++){
                                ipr.putPixelValue(u, uu,summed[u][uu]);
                            }
                        }
                        imsres.addSlice(ipr);
                    }
                    else if (i>merge){
                        
                        FloatProcessor ipr = new FloatProcessor(width,height);
                        for (int u=0;u<width;u++){
                            for (int uu=0;uu<height;uu++){
                                summed[u][uu]-=al.get(0)[u][uu];
                                ipr.putPixelValue(u, uu,summed[u][uu]);
                            }
                        }
                        imsres.addSlice(ipr);
                        al.remove(0);
                        
                    }
                    
                    
                }
                
                
                for (int i=0;i<merge/2;i++){
                    FloatProcessor ipr = new FloatProcessor(width,height);
                    for (int u=0;u<width;u++){
                        for (int uu=0;uu<height;uu++){
                            summed[u][uu]-=al.get(0)[u][uu];
                            ipr.putPixelValue(u, uu,summed[u][uu]);
                        }
                    }
                    imsres.addSlice(ipr);
                    al.remove(0);
                }
                
                ImagePlus impres = new ImagePlus("result",imsres);
                impres.show();
                
                IJ.showProgress(1);
                
                
            }
        
        processing=new Boolean(false);
        
        
    }
    
    
    
    
    
    
    void test(){
        
        
        Zernike z = new Zernike(128, 32,1);
        
        double [][] image=z.getGL(6, true);
        
        ImageShow.imshow(image);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************             PRIVATE                     ***************
//    ************             CLASSES                     ***************
//    ************               ...                       ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ************                                         ***************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
//    ********************************************************************
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class MouseOptionSave implements  MouseListener{
        
        String actualPath;
        String message;
        TextField textField;
        MouseOptionSave(TextField textField,String actualPath,String message){
            this.message=message;
            this.actualPath=actualPath;
            this.textField=textField;
        }
        
        



    public void mouseClicked(MouseEvent e) 
    { 
     if (e.getClickCount() > 1) // Double-click 
     { 
      if (e.getSource() == textField) 
      { 
       String path=null;
       
        path = chooseSavingPath(actualPath,message); 
       if (path!=null){
            textField.setText(path); 
       }
       
      } 
     } 
    } 

    public void mousePressed(MouseEvent e) 
    { 

    } 

    public void mouseReleased(MouseEvent e) 
    { 

    } 

    public void mouseEntered(MouseEvent e) 
    { 

    } 

    public void mouseExited(MouseEvent e) 
    { 

    } 
    
 
 
    String chooseSavingPath(String actualPath,String message){
        
        try {
            JFileChooser dialogue ;

            dialogue = new JFileChooser(new File(actualPath));

            File fichier;
            dialogue.setDialogTitle(message);
            int x=dialogue.showSaveDialog(null);
            if (x==JFileChooser.APPROVE_OPTION) {
                fichier = dialogue.getSelectedFile();
                OpenDialog.setDefaultDirectory(fichier.getParent());
                return fichier.getPath();
                
            }
            else if (x==JFileChooser.CANCEL_OPTION){
                //IJ.log("Save canceled !");
                return null;
            }
        
        
            

        } catch (Exception e) {
                e.printStackTrace();
        }
        return null;
                
    }
    
        
    }
    
    
    
    
    
    
    boolean loadSCMOScameraFiles(int width, int height,Roi r){
        
        //SCMOS MAP loading
        ImagePlus impvar=null;
        ImagePlus impoff=null;
        ImagePlus impgain=null;
        try{
            impvar=new ImagePlus(path_SCMOSvariance);
            impoff=new ImagePlus(path_SCMOSoffset);
            impgain=new ImagePlus(path_SCMOSgain);
        }catch(Exception e){}

        double [][] scmos=null;
        if ((impvar!=null)&&(impoff!=null)&&(impgain!=null)){
            
            int w1=impvar.getWidth();
            int w2=impoff.getWidth();
            int w3=impgain.getWidth();
            int h1=impvar.getHeight();
            int h2=impoff.getHeight();
            int h3=impgain.getHeight();
            
            if ((w1!=w2)||(w1!=w3)||(w2!=w3)||(h1!=h2)||(h1!=h3)||(h2!=h3)){
                IJ.log("ERROR : SCMOS images (variance/gain/offset) should have the same size");
                IJ.log("widths "+w1+"  "+w2+"  "+w3);
                IJ.log("heights "+h1+"  "+h2+"  "+h3);
                return false;
            }
            if ((w1!=width)||(h1!=height)){
                IJ.log("ERROR : SCMOS images parameters (variance/gain/offset) and opened image should have the same size");
                return false;
            }
            int startX=0;
            int startY=0;
            if (r!=null){
                Rectangle bound=r.getBounds();
                width=bound.width;
                height=bound.height;
                startX=bound.x;
                startY=bound.y;
            }
                
                
                
                
            
            scmosvargain=new double [width][height];
            scmosvariance=new double [width][height];
            scmosgain=new double [width][height];
            scmosoffset=new double [width][height];
            
            ImageProcessor ipvar=impvar.getProcessor();
            ImageProcessor ipoff=impoff.getProcessor();
            ImageProcessor ipgain=impgain.getProcessor();
            for (int p=0;p<width;p++){
                for (int pp=0;pp<height;pp++){
                    scmosvariance[p][pp]=ipvar.getPixelValue(startX+p, startY+pp);
                    scmosoffset[p][pp]=ipoff.getPixelValue(startX+p, startY+pp);
                    scmosgain[p][pp]=ipgain.getPixelValue(startX+p, startY+pp);
                    scmosvargain[p][pp]=scmosvariance[p][pp]/(scmosgain[p][pp]*scmosgain[p][pp]);
                }
            }
            return true;
        }
        else{
            IJ.log("ERROR, problem loading SCMOS images (variance/gain/offset)");
            return false;
        }
                
    }
    
    
       
    
    
    class MouseOptionLoad implements  MouseListener{
        
        String message;
        String actualPath;
        TextField textField;
        MouseOptionLoad(TextField textField,String actualPath,String message){
            this.message=message;
            this.actualPath=actualPath;
            this.textField=textField;
        }
        
        
        
    public void mouseClicked(MouseEvent e) 
    { 
     if (e.getClickCount() > 1) // Double-click 
     { 
      if (e.getSource() == textField) 
      { 
       String path=null;
       
        path = chooseLoadingPath(message); 
       if (path!=null){
            textField.setText(path); 
       }
      } 
     } 
    } 
    public void mousePressed(MouseEvent e) 
    { 

    } 

    public void mouseReleased(MouseEvent e) 
    { 

    } 

    public void mouseEntered(MouseEvent e) 
    { 

    } 

    public void mouseExited(MouseEvent e) 
    { 

    } 
    
    
    
    String chooseLoadingPath(String message){
        try {
            
            JFileChooser dialogue ;
            dialogue = new JFileChooser(new File(actualPath));
            File fichier;
            
            dialogue.setDialogTitle(message);
            int x=dialogue.showOpenDialog(null);
            if (x==JFileChooser.APPROVE_OPTION) {
                
                fichier = dialogue.getSelectedFile();
                
                OpenDialog.setDefaultDirectory(fichier.getParent());
                return fichier.getPath();
                
            }
            else if (x==JFileChooser.CANCEL_OPTION){
                //IJ.log("Load canceled !");
                return null;
            }
        
        
            

        } catch (Exception e) {
                e.printStackTrace();
        }
        return null;
                
    }
    
    
    }
    
    
    void printErrorMessage(){
        IJ.log("Sorry but it is not possible to launch many processes at the same time");
        IJ.log("If no process is running, please, restart ImageJ");
    }
    
    
    
    
    
    
    
    
    class ShowHistogramButton extends Button  implements ActionListener{
        Choice choicefield;
        TextField tfmin,tfmax;
        ShowHistogramButton(Choice choicefield,TextField tfmin,TextField tfmax){
            this.choicefield=choicefield;
            this.setLabel("show histogram");
            addActionListener(this);
            this.tfmax=tfmax;
            this.tfmin=tfmin;
        }
        
        public void actionPerformed(ActionEvent e) 
        { 
            //Execute when button is pressed 
            Statistic.computeHist(sl,choicefield.getSelectedIndex(),100);
            double [] mm=computeMinAndMax(choicefield.getSelectedIndex());
            
            tfmin.setText(""+mm[0]);
            tfmax.setText(""+mm[1]);
        } 
        
        
        private double [] computeMinAndMax(int idVariable){
            double mini=Double.MAX_VALUE;
            double maxi=Double.MIN_VALUE;
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    PLocalization pp = sl.fl.get(i).loc.get(j);
                    double val=pp.getValueVariable(idVariable);
                    if (val<mini){
                        mini=val;
                    }
                    if (val>maxi){
                        maxi=val;
                    }

                }
            }
            double [] d = new double[2];
            d[0]=mini;
            d[1]=maxi;
            return d;
        }
        
    }
    
    
    
    
    
    
    class PreviewButton extends Button  implements ActionListener{
        ImagePlus imppb;
        Vector<TextField> texts;
        Vector<Checkbox> bools;
        Vector<TextField> nums;
        PreviewButton(ImagePlus imp,Vector<TextField> texts,Vector<TextField> nums){
            this.texts=texts;
            this.bools=bools;
            this.nums=nums;
            this.imppb=imp;
            this.setLabel("preview");
            addActionListener(this);
        }
        
        public void actionPerformed(ActionEvent e) 
        { 
            //Execute when button is pressed 
            
            
            Action a = new Action();
            a.start();
            
                
            
        } 
        
        class Action extends Thread{
            
            public Action(){
                
            }
            
            public void run(){

                String path_calibration = texts.get(0).getText();



                double gain=0;

                if (!isSCMOS){
                    gain = (double)Double.parseDouble(nums.get(0).getText());
                }

                double nwat= (double)Double.parseDouble(nums.get(1).getText());

                double zfocus= (double)Double.parseDouble(nums.get(2).getText());

                int sizePatch = (int)Double.parseDouble(nums.get(3).getText());
                //nbStream = (int)gd.getNextNumber();
                //nbThread = (int)gd.getNextNumber();

                double axialRange= (double)Double.parseDouble(nums.get(4).getText());




                int photonThreshold = (int)Double.parseDouble(nums.get(5).getText());




            
            
            
            
            
                
                if (sizePatch%2==1){
                    sizePatch++;
                }
                
                
                
                
                boolean stacked=true;
                int nbImage=imppb.getNSlices();
                if (imppb.getNSlices()==1){
                    nbImage=imppb.getNFrames();
                    stacked=false;
                }
                ImageStack [] ims=null;
                if (stacked){
                    ims = new ImageStack[1];
                    ims[0]=imppb.getStack();

                }
                
                int select=imppb.getCurrentSlice();
                
                
                
                
                Roi r=imppb.getRoi();
                int width=imppb.getWidth();
                int height=imppb.getHeight();
                
                
                
                int startx=0;
                int starty=0;
                if (r!=null){
                    Rectangle bound=imppb.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                    startx=bound.x;
                    starty=bound.y;
                }
                if ((width<sizePatch)||(height<sizePatch)){
                    IJ.log("Process stopped: Localization impossible. The image, or selected ROI should be larger than patch size");
                    IJ.log("Please, select a large rectangle in the image before running the localization");
                    IJ.log("image width:"+width+" ; height:"+height+" ; sizePatch:"+sizePatch);
                    processing=new Boolean(false);
                    return;
                }
                
                
                ImageProcessor ip;
                if (stacked){
                    ip = ims[0].getProcessor(select);
                }
                else{
                    imppb.setSlice(select);
                    ip = imppb.getProcessor();
                }
                ImageProcessor ip2=new FloatProcessor(width,height);
                
                for (int i=0;i<width;i++){
                    for (int ii=0;ii<height;ii++){
                        ip2.putPixelValue(i, ii, ip.getPixelValue(i+startx, ii+starty));
                    }
                }
                
                ImagePlus impTest=new ImagePlus("test",ip2);
                
                if (gain==0||adu==0){
                    new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
                }
                else{
                    
                    
                    
                    
                    int sizeFFT=128;
                    sizeFFT=Math.max(sizeFFT,sizePatch*2);
                    sizeFFT=Math.min(sizeFFT,maxFFT);
                    
                    
                    
                    org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
                    if (!dp.loading){
                        IJ.log("impossible to load "+path_calibration);
                        processing=new Boolean(false);
                        return;
                    }
                    dp.setSizeoutput(sizePatch);
                    dp.setNwat(nwat);
                    dp.param.Zfocus=zfocus;
                    IJ.log("preview on current frame started");

                    org.pasteur.imagej.process.cpu.LocalizationPipeline lp;
                    if (isSCMOS){

                        lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,"",scmosvariance,scmosgain,scmosoffset,scmosvargain,true);
                    }
                    else{
                        lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,"",adu, gain,offset,true);
                    }
                    sl=new StackLocalization();
                    sl=lp.detectParticles(impTest,sl);
                    IJ.log("preview finished");
                }
            }  
            
        } 
        
        
    }
    
    
    
    
    
    
    
    class AutoFocusButton extends Button  implements ActionListener{
        ImagePlus imp;
        Vector<TextField> texts;
        Vector<Checkbox> bools;
        Vector<TextField> nums;
        AutoFocusButton(ImagePlus imp,Vector<TextField> texts,Vector<TextField> nums){
            this.texts=texts;
            this.bools=bools;
            this.nums=nums;
            this.imp=imp;
            this.setLabel("auto-focus");
            addActionListener(this);
        }
        
        public void actionPerformed(ActionEvent e) 
        { 
            //Execute when button is pressed 
            
            Action a = new Action();
            a.start();
            
                
            
        } 
        
        class Action extends Thread{
            
            public Action(){
                
            }
            
            public void run(){
                String path_calibration = texts.get(0).getText();


            
                double gain=0;

                if (!isSCMOS){
                    gain = (double)Double.parseDouble(nums.get(0).getText());
                }

                double nwat= (double)Double.parseDouble(nums.get(1).getText());

                double zfocustmp= (double)Double.parseDouble(nums.get(2).getText());

                int sizePatch = (int)Double.parseDouble(nums.get(3).getText());
                //nbStream = (int)gd.getNextNumber();
                //nbThread = (int)gd.getNextNumber();

                double axialRange= (double)Double.parseDouble(nums.get(4).getText());




                int photonThreshold = (int)Double.parseDouble(nums.get(5).getText());



                
            
                
                if (sizePatch%2==1){
                    sizePatch++;
                }
                
                int select=imp.getCurrentSlice();
                ImageProcessor ip = imp.getProcessor();
                ImagePlus impTest=new ImagePlus("test",ip);
                
                
                Roi r=imp.getRoi();
                int width=imp.getWidth();
                int height=imp.getHeight();
                
                if (r!=null){
                    Rectangle bound=imp.getRoi().getBounds();
                    width=bound.width;
                    height=bound.height;
                }
                if ((width<sizePatch)||(height<sizePatch)){
                    IJ.log("Process stopped: Localization impossible. The image, or selected ROI should be larger than patch size");
                    IJ.log("Please, select a large rectangle in the image before running the localization");
                    IJ.log("image width:"+width+" ; height:"+height+" ; sizePatch:"+sizePatch);
                    processing=new Boolean(false);
                    return;
                }
                
                
                
                
                
                
                if (gain==0||adu==0){
                    new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
                }
                else{
                    
                    
                    
                    
                    int sizeFFT=128;
                    sizeFFT=Math.max(sizeFFT,sizePatch*2);
                    sizeFFT=Math.min(sizeFFT,maxFFT);
                    
                    
                    
                    if (GPU_computation){
                        MyCudaStream.init(nbStream+1);
                        org.pasteur.imagej.process.gpu.DataPhase_ dp = new org.pasteur.imagej.process.gpu.DataPhase_(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(nwat);
                        dp.param.Zfocus=zfocustmp;
                        
                        
                        
                        IJ.log("auto-focus started on gpu");


                        org.pasteur.imagej.process.gpu.LocalizationAutoFocusPipeline_ lp;
                        if (isSCMOS){
                            lp=new org.pasteur.imagej.process.gpu.LocalizationAutoFocusPipeline_(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,scmosvariance,scmosgain,scmosoffset,scmosvargain,false);
                        }
                        else{
                            lp=new org.pasteur.imagej.process.gpu.LocalizationAutoFocusPipeline_(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,adu, gain,offset,false);
                        }
                        zfocustmp=lp.detectParticles(imp);
                        
                        
                        /*classes.testClass lp =new classes.testClass(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,rescale_slope, rescale_intercept);
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);*/



                        dp.free();
                    
                    
                        MyCudaStream.destroy();
                    }
                    else{
                        org.pasteur.imagej.process.cpu.DataPhase dp = new org.pasteur.imagej.process.cpu.DataPhase(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(nwat);
                        dp.param.Zfocus=zfocustmp;
                        IJ.log("auto-focus started on cpu : NOT YET IMPLEMENTED");
                        
                        org.pasteur.imagej.process.cpu.LocalizationPipeline lp;
                        if (isSCMOS){
                            
                            //lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain,false);
                        }
                        else{
                            //lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,adu, gain,offset,false);
                        }
                        
                        
                        
                    }
                    zfocus=zfocustmp;
                    nums.get(2).setText(""+zfocus);
                    IJ.log("auto-focus finished "+zfocus);
                }
            }
            
            
        }
        
        
    }
    
    
    
    
    
    
    
}
    
    
    
    
    