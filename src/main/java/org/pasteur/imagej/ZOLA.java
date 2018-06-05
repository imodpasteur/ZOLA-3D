/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej;



import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.process.*;
import org.pasteur.imagej.postprocess.*;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.cuda.*;




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
import java.awt.Font;


/**
 *
 * @author benoit
 */
public class ZOLA implements PlugIn  {
    
    boolean isSCMOS=false;
    
    int concecutiveFrameThreshold=1000;
    int smoothingFrameNumber=100;
            
            
    int maxFFT=140;
    double dualCamMaxDistanceMergingPlan=200;
    
    double dualCamMaxDistanceMergingZ=500;
    
    double maxDistanceMergingPlan=100;
    
    double maxDistanceMergingZ=200;
    
    boolean isFlipped=true;
    
    int sizeTextString=45;
    
    Prefs prefs = new Prefs(); 
    
    String path_SCMOSvariance="";
    String path_SCMOSoffset="";
    String path_SCMOSgain="";
    
    
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
    
    double parameterFit=1;
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
    int zernikeCoef=4;
    int orderFit=1;
    int iterationNumber=15;
    
    int idVariable=0;
    int binNumberVariable=100;
    
    double drift_sub_image_sizeUm=6;
    double drift_pixelsizeNM=100;
    int drift_bin=3;
    
    boolean GPU_computation;
       
    
    double sizeRendering=20;
            
    double maximumDistance;
    
    static Model_generator mg;
    
    static StackLocalization sl;
    static Boolean processing;
    
    
    
    double [][] scmosvargain=null;
    double [][] scmosvariance=null;
    double [][] scmosgain=null;
    double [][] scmosoffset=null;
    
    
    public ZOLA(){
        
        
    }
    
    
    
    
    public void run(String command){
        
        path_result=prefs.get("Zola.path_result", IJ.getDirectory("image")+"path_result.csv");
        
        path_localization=prefs.get("Zola.pathlocalization", IJ.getDirectory("image")+"localization_table.csv");
        path_localization2=prefs.get("Zola.pathlocalization2", IJ.getDirectory("image")+"localization_table.csv");
        path_calibration=prefs.get("Zola.pathcalibration", IJ.getDirectory("image")+"calibration_table.csv");
        path_registration=prefs.get("Zola.path_registration", IJ.getDirectory("image")+"calibration_table.csv");
        
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
        
        isFlipped=(boolean)((prefs.get("Zola.isFlipped", true)));
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
        parameterFit=prefs.get("Zola.parameterFit", 1);
        parameterFit2=prefs.get("Zola.parameterFit2", 1);
        wavelength=prefs.get("Zola.wavelength", .64);
        axialside=(int)prefs.get("Zola.axialside", 0);
        zernikeCoef=(int)prefs.get("Zola.zernikeCoef", 3);
        
        idVariable=(int)prefs.get("Zola.idVariable", 0);
        binNumberVariable=(int)prefs.get("Zola.binNumberVariable", 100);
        
        
        
    
    
        iterationNumber=(int)prefs.get("Zola.iterationNumber", 30);
        
        drift_sub_image_sizeUm=prefs.get("Zola.driftsubimagesizeUm", 6);
        drift_pixelsizeNM=prefs.get("Zola.driftpixelsizeNM", 50);
        drift_bin=(int)prefs.get("Zola.driftbin", 3);
        
        
        GPU_computation=(boolean)prefs.get("Zola.GPU_computation", false);
        
        sizeRendering=prefs.get("Zola.sizeRendering", 20);
        
        dualCamMaxDistanceMergingPlan=prefs.get("Zola.dualCamMaxDistanceMergingPlan", 250);
        
        dualCamMaxDistanceMergingZ=prefs.get("Zola.dualCamMaxDistanceMergingZ", 500);
        
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
        
        
        if (command.startsWith("Calibration")){
            
            calibration();
            
        }
        
        
        else if (command.startsWith("Filtering")){
            
            filtering();
            
        }
        
        
        
        
        else if (command.startsWith("CRLB_FromFile")){
            
            crlbfromfile();
            
            
        }
        else if (command.startsWith("CRLB")){
            
            crlb();
            
            
        }
        
        else if (command.startsWith("stathisto")){
            
            stathisto();
            
            
        }
        
        else if (command.startsWith("setcameraEMCCD")){
            
            setCameraEMCCD();
            
            
        }
        else if (command.startsWith("setcameraSCMOS")){
            
            setCameraSCMOS();
            
            
        }
        else if (command.startsWith("statmean")){
            
            statmean();
            
            
        }
        
        
        else if (command.startsWith("Condensation")){
            
            condensation();
            
            
        }
        
        
        else if (command.startsWith("test")){
            
            test();
            
            
        }
        
        
        else if (command.startsWith("driftCorrection")){
            
            driftCorrection();
            
            
        }
        
        else if (command.startsWith("resetDrift")){
            
            resetDriftCorrection();
            
            
        }
        
        
        
        
        else if (command.startsWith("Localization")){
            
            localization();
            
            
            
        }
        
        
        
        else if (command.startsWith("colorHist")){
            
            colorHist();
            
        }
        else if (command.startsWith("smartColorHist")){
            
            smartColorHist();
            
        }
        else if (command.startsWith("hist")){
            
            hist();
            
        }
        
        else if (command.startsWith("colorize")){
            
            colorizeHist();
            
        }
        
        
        
        else if (command.startsWith("Import")){
            importData();
            
        }
        else if (command.startsWith("Append")){
            appendData();
            
        }
        
        else if (command.startsWith("Export")){
            exportData();
            
        }
        
        
        else if (command.startsWith("phaseoptim")){
            phaseoptim();
            
        }
        
        
        else if (command.startsWith("getModels_GPUalloc")){
            getModels_GPUalloc();
        }
        else if (command.startsWith("getModels_GPUcomput")){
            getModels_GPUcomput();
        }
        else if (command.startsWith("getModels_GPUfree")){
            getModels_GPUfree();
        }
        
        
        
        
        prefs.set("Zola.path_SCMOSvariance", path_SCMOSvariance);
        prefs.set("Zola.path_SCMOSgain", path_SCMOSgain);
        prefs.set("Zola.path_SCMOSoffset", path_SCMOSoffset);
        prefs.set("Zola.GPU_computation", GPU_computation);
        prefs.set("Zola.pathlocalization", path_localization);
        prefs.set("Zola.pathlocalization2", path_localization2);
        prefs.set("Zola.path_registration", path_registration);
        if (path_calibration.length()>2){
            prefs.set("Zola.pathcalibration", path_calibration);
        }
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
        
        
        prefs.set("Zola.isFlipped", isFlipped);
        prefs.set("Zola.is3Drendering", is3Drendering);
        prefs.set("Zola.showLUT", showLUT);
        prefs.set("Zola.isCOLORrendering", isCOLORrendering);
        prefs.set("Zola.na", na);
        prefs.set("Zola.noil", noil);
        prefs.set("Zola.nwat", nwat);
        prefs.set("Zola.sigma", sigma);
        prefs.set("Zola.parameterFit", parameterFit);
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
        
        prefs.set("Zola.dualCamMaxDistanceMergingZ", dualCamMaxDistanceMergingZ);
        
        prefs.set("Zola.maxDistanceMergingPlan", maxDistanceMergingPlan);
        
        prefs.set("Zola.maxDistanceMergingZ", maxDistanceMergingZ);      
                
        prefs.savePreferences();
    }
    
    
    
    
    void localization(){
        
        
        nbStream=10;
        nbThread=120;  
            
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
//                int nn=(path.lastIndexOf("."));
//                if ((nn>path.length()-5)&&(nn>=0)){
//                    path=path.substring(0,nn);
//                }
                path=path+"ZOLA_localization_table.csv";
                
                
                GenericDialog gd = new GenericDialog("ZOLA: Localization");
                
                
                Font font = gd.getFont();
                Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
                
                gd.addCheckbox("Run_on_GPU :", GPU_computation);
                
                
                
                
                
                
                if (!isSCMOS){
                    gd.addMessage("Camera is EMCCD",fontBold);
                    gd.addMessage("ADU = "+adu);
                    gd.addMessage("Gain = "+gain);
                    gd.addMessage("Offset = "+offset);
                    
                }
                else{
                    if (scmosLoad){
                        gd.addMessage("Camera is SCMOS",fontBold);
                    }
                    else{
                        gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                    }
                    String sep=File.separator;
                    String [] var=this.path_SCMOSvariance.split(sep);
                    String [] off=this.path_SCMOSoffset.split(sep);
                    String [] gain=this.path_SCMOSgain.split(sep);
                    gd.addMessage("Offset file = "+off[var.length-1]);
                    gd.addMessage("Variance file = "+var[var.length-1]);
                    gd.addMessage("Gain file = "+gain[var.length-1]);
                    
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
                if ((nnn>=0)&&(nnn>path_localization.length()-5)){
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
                        DataPhase dp = new DataPhase(sizeFFT,path_calibration);
                        if (!dp.loading){
                            IJ.log("impossible to load "+path_calibration);
                            processing=new Boolean(false);
                            return;
                        }
                        dp.setSizeoutput(sizePatch);
                        dp.setNwat(nwat);
                        dp.param.Zfocus=zfocus;
                        IJ.log("Localization started");


                        LocalizationPipelineMany lp;
                        if (isSCMOS){
                            lp=new LocalizationPipelineMany(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain);
                        }
                        else{
                            lp=new LocalizationPipelineMany(dp,axialRange,photonThreshold,nbStream,nbThread,sizePatch,path_localization,adu, gain,offset);
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
                            
                            lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,scmosvariance,scmosgain,scmosoffset,scmosvargain);
                        }
                        else{
                            lp=new org.pasteur.imagej.process.cpu.LocalizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,adu, gain,offset);
                        }
                        sl=new StackLocalization();
                        sl=lp.detectParticles(imp,sl);
                        
                        
                    }
                    
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
                
        gd.addCheckbox("Run_on_GPU :", GPU_computation);
        
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
            DataPhase dp = new DataPhase(sizeFFT,path_calibration);
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
                CRLB crlb;
    //            if (dp2!=null){
    //                dp2.setSizeoutput(sizePatchPhaseRet);
    //                crlb =new CRLB(dp,dp2,-axialRange/2.,axialRange/2.,stepZloc,photonNumber,background,maxValuePlot); 
    //            }
    //            else{
                    crlb =new CRLB(dp,axialRange,stepZloc,photonNumber,background,-1); 
    //            }


                crlb.run(null);

            }

            MyCudaStream.destroy();
        }
        else{
            DataPhase dp = new DataPhase(sizeFFT,path_calibration);
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
                CRLB crlb;
    //            if (dp2!=null){
    //                dp2.setSizeoutput(sizePatchPhaseRet);
    //                crlb =new CRLB(dp,dp2,-axialRange/2.,axialRange/2.,stepZloc,photonNumber,background,maxValuePlot); 
    //            }
    //            else{
                    crlb =new CRLB(dp,axialRange,stepZloc,photonNumber,background,-1); 
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
    
    
    
    
    
    
    
    void crlbfromfile(){
        
        if ((processing !=null)&&(processing.booleanValue()==true)){
            printErrorMessage();
            return;
        }
        
        processing=new Boolean(true);
            
            
        

        MyCudaStream.init(1);


        GenericDialog gd = new GenericDialog("Fundamental limit (CRLB)");
        
        
        
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}

        String pathdir=IJ.getDirectory("image");
        if (pathdir==null){
            pathdir=IJ.getDirectory("current");
        }
        
        String pathResult=pathdir+"crlbtable.csv";
        
        
        
        String path_calibration2="";
        
        gd.addStringField("Calibration file:",path_calibration,sizeTextString); 
        
        gd.addStringField("Calibration_file_cam2 [optional]:",path_calibration2,sizeTextString); 
        
        
        gd.addMessage("Computation parameter",fontBold);

        
        gd.addNumericField("Size of patch: ", sizePatchPhaseRet,0);
        
        gd.addNumericField("Axial_range (µm): ", axialRange,2);

        gd.addMessage("CRLB parameters",fontBold);
        
        
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
        
        MyCudaStream.init(1);

        DataPhase dp = new DataPhase(sizeFFT,path_calibration);
        if (!dp.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        DataPhase dp2=null;
        if (path_calibration2.length()>3){
            dp2 = new DataPhase(sizeFFT,path_calibration2);
            if (!dp2.loading){
                IJ.log("impossible to load "+path_calibration2);
                processing=new Boolean(false);
                return;
            }
        }
        
        dp.setSizeoutput(sizePatchPhaseRet);
        dp.param.Zfocus=0;
        CRLB crlb; 
        String path=FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B");
//        if (dp2!=null){
//            String path2=classes.FileVectorLoader.chooseLoadingPath("Select a file with x , y , z , A , B for cam 2");
//            dp2.setSizeoutput(sizePatchPhaseRet);
//            crlb =new CRLB(dp,dp2,path,path2); 
//        }
//        else{
            crlb =new CRLB(dp,path,axialRange); 
//        }
        
        
        
        
        

        crlb.run(pathResult);

        dp.free();

        MyCudaStream.destroy();

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
            
            int zernikeCoef=15;
            String path=IJ.getDirectory("image");
            if (path==null){
                path=IJ.getDirectory("current");
            }
//            int nn=(path.lastIndexOf("."));
//            if ((nn>path.length()-5)&&(nn>=0)){
//                path=path.substring(0,nn);
//            }
            path=path+"phase_optimized.csv";
            

            GenericDialog gd = new GenericDialog("Phase optimization!");


            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}                
            gd.addMessage("This method is works only with GPU",fontBold);
            
            gd.addMessage("Image parameters",fontBold);
            
            gd.addNumericField("Pixel size (µm): ", xystep, 3);
            
            
            gd.addNumericField("Numerical aperture: ", na, 3);
            gd.addNumericField("Immersion_refractive index: ", noil, 3,6,"~1.518 for oil, ~1.33 for water");
            gd.addNumericField("Mounting_medium_refractive index: ", nwat, 3,6,"1.33 for water");
            
            gd.addNumericField("Distance_focus_to_coverslip (µm): ", zfocus, 3);
            gd.addNumericField("Wavelength: ", wavelength, 3);

            gd.addMessage("Optimization parameters",fontBold);
            
            
            gd.addNumericField("Axial_range (µm): ", axialRange,2);
            gd.addNumericField("Z_step (µm): ", stepZloc,2);
            
            gd.addNumericField("Photon number: ", photonNumber,0);
            gd.addNumericField("Background intensity: ", background,0);

            
            gd.addNumericField("Zernike coefficient #: ", zernikeCoef, 0);
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

                
                
                
                zernikeCoef = (int)gd.getNextNumber();
                iterationNumber = (int)gd.getNextNumber();
                
                pathResult = gd.getNextString(); 
            }
            else{
                processing=new Boolean(false);
                return;
            }
            
            
            int sizeFFT=128;

            MyCudaStream.init(1);
            int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+0);}
            
            PhaseOptimization po = new PhaseOptimization(sizeFFT,sizeFFT,xystep,axialRange, stepZloc,photonNumber, background,wavelength,noil,nwat,zfocus,na,zernikeCoef,pathResult);
            po.run(iterationNumber);
            po.free();
            
            
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
//            if ((nn>path.length()-5)&&(nn>=0)){
//                path=path.substring(0,nn);
//            }
            path=path+"ZOLA_calibration_PSF.csv";
            
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
            String [] fieldZernike = {"6","10","15","21","28","36","45"};
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
            
            gd.addCheckbox("Run_on_GPU :", GPU_computation);
            
            
                if (!isSCMOS){
                    gd.addMessage("Camera is EMCCD",fontBold);
                    gd.addMessage("ADU = "+adu);
                    gd.addMessage("Gain = "+gain);
                    gd.addMessage("Offset = "+offset);
                    
                }
                else{
                    
                    if (scmosLoad){
                        gd.addMessage("Camera is SCMOS",fontBold);
                    }
                    else{
                        gd.addMessage("Camera is SCMOS / problem with image sizes / ERROR will occur",fontBold);
                    }
                    String sep=File.separator;
                    String [] var=this.path_SCMOSvariance.split(sep);
                    String [] off=this.path_SCMOSoffset.split(sep);
                    String [] gain=this.path_SCMOSgain.split(sep);
                    gd.addMessage("Offset file = "+off[var.length-1]);
                    gd.addMessage("Variance file = "+var[var.length-1]);
                    gd.addMessage("Gain file = "+gain[var.length-1]);
                    
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
                
                
                xystep = (double)gd.getNextNumber()/1000;
                zstep = (double)gd.getNextNumber()/1000;
                
                axialside=gd.getNextChoiceIndex();
                
                na = (double)gd.getNextNumber();
                noil = (double)gd.getNextNumber();
                wavelength = (double)gd.getNextNumber();
                
                sizePatchPhaseRet = (int)gd.getNextNumber();
                
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
            if ((nnn>=0)&&(nnn>path_calibration.length()-5)){
                path_calibration=path_calibration.substring(0,nnn);
            }
            path_calibration+=".csv";
            
            if (sizePatchPhaseRet%2==1){
                sizePatchPhaseRet++;
            }
            
            int sizeFFT=128;
            sizeFFT=Math.max(sizeFFT,sizePatchPhaseRet*2);
            sizeFFT=Math.min(sizeFFT,maxFFT);
            
            
            
            
            if (gain==0||adu==0){
                new WaitForUserDialog("Error message", "gain has to be positive (default value = 1)").show();
            }
            else{
                IJ.log("ZOLA PSF modeling started");
                
                InitBackgroundAndPhotonNumberMany paramImage ;
                if (isSCMOS){
                    paramImage = new InitBackgroundAndPhotonNumberMany(image,xp,yp,sizePatchPhaseRet,this.scmosvariance,this.scmosgain,this.scmosoffset,this.scmosvargain,zstep);

                }
                else{
                    paramImage = new InitBackgroundAndPhotonNumberMany(image,xp,yp,sizePatchPhaseRet,adu,gain,offset,zstep);

                }
                long t=System.currentTimeMillis();
                
                if (GPU_computation){
                    //maybe test if image is a square
                    //IJ.log("init...");
                    MyCudaStream.init(1);
                    //IJ.log("init cuda ok "+xp+"  "+yp+"  "+sizePatchPhaseRet+"  "+rescale_slope+"  "+rescale_intercept+"  "+zstep+"  "+image);
                    
                    //IJ.log("param ok");
                    PhaseRetrievalProcessMany prp = new PhaseRetrievalProcessMany(sizeFFT,xystep,zstep,wavelength,noil,na,Integer.parseInt(fieldZernike[zernikeCoef]),paramImage,path_calibration,sigma,axialside);
                    //IJ.log("process many ok");
                    prp.run(iterationNumber);
                    //IJ.log("run ok");
                    prp.free();
                    MyCudaStream.destroy();
                }
                else{
                    
                    org.pasteur.imagej.process.cpu.PhaseRetrieval prp = new org.pasteur.imagej.process.cpu.PhaseRetrieval(sizeFFT,xystep,zstep,wavelength,noil,na,Integer.parseInt(fieldZernike[zernikeCoef]),paramImage,path_calibration,sigma,axialside);
                    prp.run(iterationNumber);
                }
                long t2=System.currentTimeMillis();
                IJ.log("elapsed time "+Math.floor((double)(t2-t)/(double)(1000*60))+" min "+Math.floor(60*((((double)(t2-t)/(double)(1000*60)))-Math.floor((double)(t2-t)/(double)(1000*60))))+" sec");
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
        
        ZRendering.smartColorRendering(null,sl, 20,1,true);
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
            gd.addCheckbox("Run_on_GPU :", GPU_computation);
            
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
                new DriftCorrection(sl,drift_pixelsizeNM, Integer.parseInt(binNumber[drift_bin]),drift_sub_image_sizeUm,path);
            }
            else{
                new org.pasteur.imagej.process.cpu.DriftCorrection(sl,drift_pixelsizeNM, Integer.parseInt(binNumber[drift_bin]),drift_sub_image_sizeUm);
            }
            ZRendering.smartColorRendering(null,sl, 20,1,true);
            
            IJ.log("Please, export the localization table to save drift correction");
            
            processing=new Boolean(false);
        }
        else{
            IJ.log("Please, import first a localization table");
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
                
                ZRendering.smartColorRendering(null,sl, 20,1,true);

                
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
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            double minX=Double.POSITIVE_INFINITY,maxX=Double.NEGATIVE_INFINITY,minY=Double.POSITIVE_INFINITY,maxY=Double.NEGATIVE_INFINITY,minZ=Double.POSITIVE_INFINITY,maxZ=Double.NEGATIVE_INFINITY;
            double maxCRLBX=Double.NEGATIVE_INFINITY,maxCRLBY=Double.NEGATIVE_INFINITY,maxCRLBZ=Double.NEGATIVE_INFINITY;
            double maxScore=Double.NEGATIVE_INFINITY;
            
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    PLocalization p = sl.fl.get(i).loc.get(j);
                    
                    
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
            
            
            
            gd.addNumericField("Max_chi2: ", maxScore,2,6,">1");
            
            

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
                ZRendering.smartColorRendering(null,sl, 20,1,true);

                
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
    
    
    
    
    
    
    
    
    
    
    void scatter(){
        
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

            GenericDialog gd = new GenericDialog("ZOLA: Render scatter plot");
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}

            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
            
            gd.addNumericField("Min_X: ", Math.floor(minX),0,6,"(nm)");
            gd.addNumericField("Max_X: ", Math.ceil(maxX),0,6,"(nm)");
            gd.addNumericField("Min_Y: ", Math.floor(minY),0,6,"(nm)");
            gd.addNumericField("Max_Y: ", Math.ceil(maxY),0,6,"(nm)");
            gd.addNumericField("Min_Z: ", Math.floor(minZ),0,6,"(nm)");
            gd.addNumericField("Max_Z: ", Math.ceil(maxZ),0,6,"(nm)");
            
            
            gd.addCheckbox("3D_rendering :", is3Drendering);
            gd.addCheckbox("Color_rendering :", isCOLORrendering);
            gd.showDialog();

            if (!gd.wasCanceled()){

                sizeRendering = (double)gd.getNextNumber();
                minX = (double)gd.getNextNumber();
                maxX = (double)gd.getNextNumber();
                minY = (double)gd.getNextNumber();
                maxY = (double)gd.getNextNumber();
                minZ = (double)gd.getNextNumber();
                maxZ = (double)gd.getNextNumber();
                is3Drendering = (boolean)gd.getNextBoolean();
                isCOLORrendering = (boolean)gd.getNextBoolean();
            }
            else{
                return;
            }

            if (is3Drendering){
                ZRendering.scatter3D(sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,isCOLORrendering);
            }
            else{
                ZRendering.scatter2D(sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,isCOLORrendering);
            }
        }
        else{
            IJ.log("please import a localization table first");
        }
        
    }
    
    void hist(){
        
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
            
            GenericDialog gd = new GenericDialog("ZOLA: Render histogram");
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
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
                
                ZRendering.hist3D(sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering);
            }
            else{
                ZRendering.hist2D(sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering);
            }
            
        }
        else{
            IJ.log("please import a localization table first");
        }
    }
    
    
    void smartColorHist(){
        
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
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
            gd.addNumericField("Shift_histogram: ", shiftrendering,0,6,"(pixels)");
            gd.addCheckbox("Show_calibration_bar:", showLUT);
            gd.showDialog();

            if (!gd.wasCanceled()){

                sizeRendering = (double)gd.getNextNumber();
                shiftrendering = (int)gd.getNextNumber();
                showLUT = (boolean)gd.getNextBoolean();
            }
            else{
                return;
            }
            if (shiftrendering<0){
                shiftrendering=0;
                IJ.log("render shift should be at least 0");
            }
            
            
            
            ZRendering.smartColorRendering(null,sl,sizeRendering,shiftrendering,showLUT);
            
            
        }
        else{
            IJ.log("please import a localization table first");
        }
    }
    
    
    
    
    
    
    void colorHist(){
        
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
            Font font = gd.getFont();
            Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
            gd.addNumericField("Pixel size: ", sizeRendering,1,6,"(nm)");
            
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

            if (!gd.wasCanceled()){

                sizeRendering = (double)gd.getNextNumber();
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
            
            
            
            ZRendering.colorRendering(null,sl,sizeRendering,minX,maxX,minY,maxY,minZ,maxZ,shiftrendering,showLUT);
            
            
        }
        else{
            IJ.log("please import a localization table first");
        }
    }
    
    
    
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
        
        
        org.pasteur.imagej.process.DataPhase dp = new org.pasteur.imagej.process.DataPhase(sizeFFT,path_calibration);
        if (!dp.loading){
            IJ.log("impossible to load "+path_calibration);
            processing=new Boolean(false);
            return;
        }
        dp.setSizeoutput(sizePatch);
        dp.setNwat(nwat);
        dp.setMany(model_number);
        
        
        
        
        
        
        mg=new Model_generator(dp,sizeImage);
        
        
        
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
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
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
            
            sl=new StackLocalization(path_localization);
            
        }
        
        
        ZRendering.smartColorRendering(null,sl, 20,1,true);
        
        
        IJ.log("localization table loaded");
        
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
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
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
            if (sl==null){
                sl=new StackLocalization(path_localization);
            }
            else{
                sl.append(path_localization);
            }
        }
        
        ZRendering.smartColorRendering(null,sl, 20,1,true);
        
        IJ.log("localization table loaded");
        
        processing=new Boolean(false);
    }
    
    
    
    
    void exportData(){
        
        
        TextField textField; 
        
        
        GenericDialog gd = new GenericDialog("ZOLA: Export localization table");
        Font font = gd.getFont();
        Font fontBold= gd.getFont();
            try{
                fontBold=new Font(font.getName(),Font.BOLD,font.getSize());
            }catch(Exception e){}
            
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
            if (sl!=null){
                sl.save(path_localization);
            }
        }
        
        IJ.log("localization table saved");
        
        
    }
    
    
    
    
    void test(){
        IJ.log("try to init GPU...");
        
        MyCudaStream.init(1);
        
        IJ.log("GPU init ok");
        
    }
    
    
    
    
    
    
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
    
    
    
    
    
    
    
    
    
}
