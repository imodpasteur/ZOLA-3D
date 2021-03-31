/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import jcuda.Pointer;
import org.pasteur.imagej.data.StackLocalization;
import org.pasteur.imagej.data.FrameLocalization;
import org.pasteur.imagej.data.PLocalization;
import org.pasteur.imagej.postprocess.ZRendering;
import org.pasteur.imagej.cuda.*; 
import org.pasteur.imagej.utils.*;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.ImageProcessor;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.awt.Rectangle;
import java.util.ArrayList;
import jcuda.runtime.JCuda;

import java.util.concurrent.locks.ReentrantLock;
import jcuda.runtime.cudaError;

/**
 *
 * @author benoit
 */
public class LocalizationPipelineHD_ {
     
    
    
    //static ArrayList<String> frameVariable = new ArrayList<String>();
    
    static ArrayList<Integer> freeLocalization = new ArrayList<Integer>();
    
    boolean show;
    
    String path_localization;
    
    StackLocalization stackloc;
            
    double rescale_slope; 
    double rescale_intercept;
    
    PlotResult pr ;
    
    ReentrantLock lock = new ReentrantLock();
    boolean predetectionFinished=false;
    int totNumberImageProcessed=0;
    int totNumberImagePreDetected=0;
    ArrayList<Integer> imageInBuffer = new ArrayList<Integer>();
    //ArrayList<String> frameNameInBuffer = new ArrayList<String>();
    int nblocalizationProcess;
    
    int sizePatch=64;
    
    int xinit;
    int yinit;
    int width;
    int height;
    int sizeImage;
    int sliceNumber;
    double thresholdCrossCorrelation=.1;
    
    int nbGlobalMask=100;
    int computedGlobalMask=0;
    double maxSNR=1; 
    
    //double [] rangePredetection;
    
    double [][][] image;
    
    
    double photonThreshold=500;
    DataPhase_  dp;
    
    int iterMaxDetection=30;
    int iterMaxLocalization=30;
    
    
    
    
    
    double axialRange;
    double minZ;
    double maxZ;
    double stepZ;
    
    
    int idLoc=0;
    SCMOScamera scmoscam=null;
    boolean isSCMOS;
    double adu=1;
    double gain=1;
    double offset=0;
    
    
    LocalizationThread [] lt;
    DataPhase_  [] dploc;//liste because multiple localizations
    Object [] monitor;
    
    public LocalizationPipelineHD_(DataPhase_ dp,double axialRange,int photonThreshold,int nblocalizationProcess,int sizePatch,String path_localization,double adu, double gain,double offset,boolean show){
        
        isSCMOS=false;
        this.adu=adu;
        this.gain=gain;
        this.offset=offset;
        
        this.rescale_slope=adu/gain;
        this.rescale_intercept=-offset*adu/gain;
        
        
        localizationPipelineMany(dp,axialRange,photonThreshold,nblocalizationProcess,sizePatch,path_localization,show);
        
    }
    
    public LocalizationPipelineHD_(DataPhase_ dp,double axialRange,int photonThreshold,int nblocalizationProcess,int nbThread,int sizePatch,String path_localization,SCMOScamera scmoscam,boolean show){
        
        isSCMOS=true;
        this.scmoscam=scmoscam;
        
        
        localizationPipelineMany(dp,axialRange,photonThreshold,nblocalizationProcess,sizePatch,path_localization,show);
        
    }
    
    
    
    
    
    
    
    
    void localizationPipelineMany(DataPhase_  dp,double axialRange,int photonThreshold,int nblocalizationProcess,int sizePatch,String path_localization,boolean show){
        
        /*this.scmos=scmos;
        if (scmos!=null){
            isSCMOS=true;
        }
        else{
            isSCMOS=false;
        }
        
        this.rescale_slope=rescale_slope;
        this.rescale_intercept=rescale_intercept;*/
        
        //frameVariable.add("name");
        //show=true;
        this.show=show;
        this.photonThreshold=photonThreshold;
        this.sizePatch=sizePatch;
        this.axialRange=axialRange;
        this.nblocalizationProcess=nblocalizationProcess;
        
        this.stepZ=dp.param.xystep*2;//Same precision according to Z than for X and Y during detection
        this.dp=dp;
        
        this.path_localization=path_localization;
        
        
        
        
    }
    
    
    
    
    
    
    
    public StackLocalization detectParticles(ImagePlus imp,StackLocalization stacklocinput){
        
        stackloc=stacklocinput;
        
        long timeBegin=System.currentTimeMillis();
        
        idLoc=0;
        
        
        IJ.showProgress(0);
        
        
        double posit=dp.param.Zfocus;
        //if (dp[i].param.noil!=dp[i].param.nwat){removed because maybe the center is not well positioned when noil=nwat

        SearchPSFcenter_ spc = new SearchPSFcenter_(dp,axialRange);
        posit=spc.getPosition();

        dp.param.ZfocusCenter=posit;
        
        
        
        
        
        this.minZ=dp.param.ZfocusCenter-axialRange/2.;
        this.maxZ=dp.param.ZfocusCenter+axialRange/2.;
//        IJ.write("\"id\",\"frame\",\"x [nm]\",\"y [nm]\",\"z [nm]\",\"intensity [photon]\",\"background\",\"likelihood\",\"residual\",\"crlbX\",\"crlbY\",\"crlbZ\",\"predetectionX\",\"predetectionY\",\"predetectionZ\"");
        
        
        
        
        
        
        
        
        sliceNumber=imp.getNSlices();
        if (imp.getNSlices()==1){
            sliceNumber=imp.getNFrames();
        }
        
        this.image=new double [sliceNumber][][];
        for (int i=0;i<sliceNumber;i++){
            
            
            image[i]=null;
            
        }
        
        
        
        
        
        
            
        width=imp.getWidth();
        height=imp.getHeight();
        xinit=0;
        yinit=0;
        Roi r=imp.getRoi();
        if (r!=null){
            Rectangle bound=imp.getRoi().getBounds();
            xinit=bound.x;
            yinit=bound.y;
            width=bound.width;
            height=bound.height;
        }
        
        
            
        
        sizeImage=Math.max(width,height);
        if (sizeImage%2==1){
            sizeImage++;
        }
        
            
        
        

        

        dploc=new DataPhase_[nblocalizationProcess];
                
        
        for (int ip=0;ip<nblocalizationProcess;ip++){
            dploc[ip]=new DataPhase_(dp,ip+1);
        }
        
        
        //dp.setMany(this.nbThread,isSCMOS);
        
        
        PreDetectionThread pdt = new PreDetectionThread(imp,sliceNumber);
        
        
        
        //Killer killer = new Killer(150000);
        //killer.start();

        
        
        
         
            
            
        
        pr = new PlotResult(30,width,height);//update every 60 second
        
        
        boolean saveOnTheFly=false;
        if (path_localization!=null){
            if (path_localization.length()>2){
                stackloc.saveOnTheFly(path_localization);
                saveOnTheFly=true;
            }
        }
        
        
        
        
        //Killer killer = new Killer(150000);
        
        
        GarbageLauncher gcl = new GarbageLauncher(1800000);//system gc every 30 min
        gcl.start();
        
        
        

        

        monitor = new Object[nblocalizationProcess];
        lt= new LocalizationThread[nblocalizationProcess];
        
        for (int ip=0;ip<nblocalizationProcess;ip++){
            monitor[ip] = new Object();
            //IJ.log("WARNING: dp should have another stream");
            lt[ip]=new LocalizationThread(monitor[ip],ip,dploc[ip]);
            
        }
        
        
        
        
        
        pdt.start();//main process
        
        for (int k=0;k<nblocalizationProcess;k++){
            synchronized(monitor[k]){
                lt[k].start();
                try{
                    monitor[k].wait();//to be sure the threads are launched at this stage and that wait() in the thread is called
                }catch(Exception er){IJ.log("oops: problem wait function for sync monitor 1");}
            }
        }
        
        //pr.start();
        IJ.log("WARNING: plot not started");
        
        
        //the program is running here
        
        
        try{
            pdt.join();
        }catch(Exception eeee){System.out.println("join pdt impossible");}
        
        
        for (int ip=0;ip<lt.length;ip++){
            if (IJ.escapePressed()){//if escape press -> kill localisation
                synchronized(monitor[ip]){
                    lt[ip].kill();
                    lt[ip].getMonitor().notify();
                }
            }
            try{
                lt[ip].join();
            }catch(Exception eeee){System.out.println("join lt impossible");}
             lt[ip].free();
             dploc[ip].free();
        }
        
        
        
        
        //free();
        
        
        pr.stop_();
        
        if (saveOnTheFly){
            stackloc.stopSaveOnTheFly();
        }
        
        gcl.stopRun();
        
        
        long timeEnd=System.currentTimeMillis();
        
        IJ.showProgress(0);
        
        
        
        IJ.log("elapsed time = "+((double)(timeEnd-timeBegin))/60000.+" min");
        IJ.log("localization number: "+this.idLoc);
        if (path_localization.length()>3){
            String ppath=path_localization;
            int nn=(ppath.lastIndexOf("."));
            if (nn>ppath.length()-5){
                ppath=ppath.substring(0,nn);
            }
            ppath=ppath+"_ZolaProtocol.json";
            //Protocol_localization.saveProtocolJSON(ppath);
            if (!isSCMOS){      
                Protocol_localization.saveProtocolJSON(ppath,dp.param.pathcalib, path_localization,(((double)(timeEnd-timeBegin))/60000.), dp.param.Zfocus,dp.param.nwat,sliceNumber,totNumberImageProcessed,IJ.escapePressed(),this.idLoc,sizePatch,this.isSCMOS,this.adu, this.gain,this.offset,"","","");
            }
            else{
                Protocol_localization.saveProtocolJSON(ppath,dp.param.pathcalib, path_localization,(((double)(timeEnd-timeBegin))/60000.), dp.param.Zfocus,dp.param.nwat,sliceNumber,totNumberImageProcessed,IJ.escapePressed(),this.idLoc,sizePatch,this.isSCMOS,1, 1,0,this.scmoscam.path_SCMOSvariance,this.scmoscam.path_SCMOSgain,this.scmoscam.path_SCMOSoffset);

            }
            //Protocol_localization.saveProtocolJSON(ppath,dp[0].param.pathcalib, path_localization,(((double)(timeEnd-timeBegin))/60000.), dp[0].param.Zfocus,dp[0].param.nwat,sliceNumber,totNumberImageProcessed,IJ.escapePressed(),this.idLoc,sizePatch,this.isSCMOS,this.adu, this.gain,this.offset,this.scmoscam.path_SCMOSvariance,this.scmoscam.path_SCMOSgain,this.scmoscam.path_SCMOSoffset);
            
            //this.save(path_localization, "Reconstruction\nCalibration path:"+dp[0].param.pathcalib+"\nLocalization path:"+path_localization+"\nReconstruction time:"+(((double)(timeEnd-timeBegin))/60000.)+" min\ndistance focus to coverslip:"+dp[0].param.Zfocus+"\nRefractive index medium:"+dp[0].param.nwat+"\nframe number:"+sliceNumber+"\nframe reconstructed:"+totNumberImageProcessed+"\nReconstruction aborted:"+IJ.escapePressed()+"\nlocalization number: "+this.idLoc+"\npatch size: "+sizePatch+"\n");
        }
        
        
        
        IJ.showStatus("localization finished");
        
        
        
        
         if (path_localization!=null){
            if (path_localization.length()>2){
                IJ.log("Localization file saved");
            }
            else{
                IJ.log("Localization file not saved...you can use Zola->export plugin to save it.");
            }
        }
        else{
            IJ.log("Localization file not saved...you can use Zola->export plugin to save it.");
        }
        freeLocalization.clear();
        IJ.resetEscape();
        return stackloc;
        
    }
    
    
    
    
    
    class PlotResult extends Thread{
        int time_s;
        ReentrantLock lockpr = new ReentrantLock();
        boolean running=true;
        ImagePlus impp=null;
        int sizePix=20;
        PlotResult(int time_s,int w,int h){
            impp=new ImagePlus();
            impp.setTitle("Live rendering 20nm/pixel");
            this.time_s=time_s;
        }
        public void run(){
            
            try{
                while (running){
                    Thread.sleep(time_s*1000);
                    if (running){
                        plot();
                    }
                }
                
            }catch(Exception e){}
        }
        
        public void launch_(){
            running=true;
            plot();
            start();
        }
        
        
        public void stop_(){
            running=false;
            plot();
        }
        
        public void plot(){
            try{
                lockpr.lock();
                impp=ZRendering.colorRendering(impp,stackloc, 20, 1,true);
            }
            finally{
                lockpr.unlock();
            }
        }
    }
    
    
    
    
    
    
    
    
    class Killer extends Thread{
        int time_ms;
        Killer(int time_ms){
            this.time_ms=time_ms;
        }
        public void run(){
            try{
                int nb=time_ms/10000;
                while (nb>0){
                    Thread.sleep(10000);
                    IJ.log("killer ... "+nb);
                    nb--;
                }
                
                System.exit(10);
            }catch(Exception e){}
        }
    }
    
    
    
    
    
    class GarbageLauncher extends Thread{
        int time_ms;
        boolean launched=true;
        GarbageLauncher(int time_ms){
            this.time_ms=time_ms;
            launched=true;
        }
        
        
        public void stopRun(){
            launched=false;
        }
        
        
        public void run(){
            int tms=(time_ms/1000)+1;
            try{
                loop:while (launched){
                    for (int i=0;i<1000;i++){
                        Thread.sleep(tms);
                        if (!launched){
                            break loop;
                        }
                    }
                    
                    System.gc();
                }
                
            }catch(Exception e){}
        }
    }
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    class PreDetectionThread extends Thread{
        
        
        int nbFrameParallel=1;
        boolean killer=false;
        ImagePlus imp;
        
        
        int pushedImage=0;
        int sliceNumber;
        
        Predetection_deconvolutionFastSkellamBckg_ predGPU;
        
        PreDetectionThread(ImagePlus imp,int sliceNumber){
            int sizePsftmp=dp.param.sizeoutput*2;
            IJ.log("paramsize "+dp.param.size);
            dp.setSizeoutput(sizePsftmp);
            this.sliceNumber=sliceNumber;
            nbFrameParallel=Math.min(nbFrameParallel, sliceNumber);
            IJ.log("zstep detection="+stepZ);
            this.imp=imp;
            //double stepZ=((double)maxZ-(double)minZ)/(double)16;
            //init with param of the first camera because predetection uses only one cam
            predGPU=new Predetection_deconvolutionFastSkellamBckg_(nbFrameParallel,width,height,dp,minZ, maxZ, stepZ,photonThreshold);
            //rangePredetection=predGPU.getRange();
            
            
            /*
            org.pasteur.imagej.process.gpu.Predetection_deconvolutionFastSkellamBckg_ pds = new org.pasteur.imagej.process.gpu.Predetection_deconvolutionFastSkellamBckg_(width,height, dp,mini, maxi, zstep,300);
        pds.pushImage(image[0]);
        pds.pushImage(image[1]);
        
        
        pds.deconvolution_RL();
            */
            
            
            
            
            
        }
        
        
        public void kill(){
            this.killer=true;
        }
        
        public void run(){
            
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
            
            
            loop:for (int z=1;z<=nbImage;z++){
                //IJ.log("image detection: "+(z-1));
                
                
                
                
                ////////////////////////////////////////////////////// WAIT free localization here
                
                
                    
                if (IJ.escapePressed()){
                    IJ.log("ESCAPE BUTTON PRESSED. LOCALIZATION WILL STOP SOON.");
                    totNumberImagePreDetected=pushedImage;
                    break loop;
                }
                
                
                
                ImageProcessor ip;
                String theName="default";
                if (stacked){
                    ip = ims[0].getProcessor(z);
                    theName = ims[0].getSliceLabel(z);
                }
                else{
                    imp.setSlice(z);
                    ip = imp.getProcessor();
                    theName=imp.getStack().getSliceLabel(z);
                }

                image[z-1]=new double[width][height];
                for (int ii=0,jj=xinit;ii<width;ii++,jj++){
                    for (int iii=0,jjj=yinit;iii<height;iii++,jjj++){
                            image[z-1][ii][iii]=ip.getPixelValue(jj, jjj);

                            if (isSCMOS){
                                IJ.log("WARNING:SCMOS");
                                image[z-1][ii][iii]=((image[z-1][ii][iii]-scmoscam.scmosoffset[ii][iii])/scmoscam.scmosgain[ii][iii])+scmoscam.scmosvargain[0][0];//here I add constant value for predetection for scmosvargain
                            }
                            else{
                                image[z-1][ii][iii]=image[z-1][ii][iii]*rescale_slope+rescale_intercept;
                            }
                            if (image[z-1][ii][iii]<0){
                                image[z-1][ii][iii]=0;
                            }
                    }
                }


                


                //ImageShow.imshow(image[z-1]," im2push ");
                //predetection on channel 0

                predGPU.pushImage(image[z-1]);
                pushedImage++;
                IJ.log("image pushed "+pushedImage);
                if (pushedImage<nbFrameParallel){
                    //do nothing -> we will process when nbFrameParallel==pushedImage
                }
                else if (pushedImage==sliceNumber){
                    //all images are pushed -> optimization 50 iterations
                    predGPU.deconvolution_RL(iterMaxDetection);
                    for (int i=0;i<nbFrameParallel;i++){
                        IJ.log("process end "+(pushedImage-nbFrameParallel+i));
                        this.popImage(pushedImage-nbFrameParallel+i);
                        
                    }
                }
                else if (pushedImage==nbFrameParallel){
                    //beginning -> iterMaxDetection iterations
                    predGPU.deconvolution_RL(iterMaxDetection);
                    //for (int i=0;i<nbFrameParallel;i++){
                        IJ.log("process start "+(pushedImage-nbFrameParallel));
                        this.popImage(pushedImage-nbFrameParallel);
                        
                    //}
                }
                else {
                    IJ.log("process middle "+(pushedImage-nbFrameParallel));
                    //otherwise: push / pop of images -> images are processed many time -> iterMaxDetection/nbFrameParallel iterations
                    predGPU.deconvolution_RL((iterMaxDetection/nbFrameParallel)+1);
                    this.popImage(pushedImage-nbFrameParallel);
                    
                }

                
                
                image[z-1]=null;
                
                ip=null;

            }
            predetectionFinished=true;
            
            
            endLoc();
            
            
            predGPU.free();
            IJ.log("end predetection");
        }
        
        
        //kill localisation task when all finished
        private void endLoc(){
            while (freeLocalization.size()!=nblocalizationProcess){
                IJ.log("-->"+freeLocalization.size()+"/"+nblocalizationProcess);
                try{
                    Thread.sleep(1000);
                }catch(Exception e){}
                if (this.killer){
                    return;
                }
            }
            for (int ip=0;ip<nblocalizationProcess;ip++){
                synchronized(monitor[ip]){
                    lt[ip].kill();
                    lt[ip].getMonitor().notify();
                }
            }
        }
        
        
        
        private void popImage(int numframe){
            
            while (freeLocalization.size()==0){
                try{
                    Thread.sleep(10);
                }catch(Exception e){}
                if (this.killer){
                    return;
                }
            }
            int index=freeLocalization.get(0);
            
            synchronized(monitor[index]) {
                
                
                JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(lt[index].getStreamId()));
                //It is important to set background, image and coordinate first... 
                lt[index].setBackground(predGPU.getPointerBackground());
                lt[index].setImage(predGPU.getPointerImage(),numframe);
                lt[index].setCoordinates(predGPU.popResult());
                try{
                    lock.lock(); 
                    freeLocalization.remove(0);
                }
                catch(Exception ee){
                    IJ.log("problem lock in LocalizationThread "+ee);//do it
                }
                finally{ 
                    lock.unlock();
                }


                //then localisation can be launched in parallel
                
                monitor[index].notify();
            }
            
        }
        
        
        
        
        
    }
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class LocalizationThread extends Thread{
        
        //String frameName="";
        int threadNumber;
        
        FrameLocalization fl ;
        Object monitor;
        
        boolean killer=false;
        LocalizationHD_ lochd;
        
        LocalizationThread(Object monitor,int threadNumber,DataPhase_ dp){
            
            this.monitor=monitor;
            
            try{
                lock.lock(); 
                freeLocalization.add(threadNumber);
            }
            catch(Exception ee){
                IJ.log("problem lock in LocalizationThread "+ee);//do it
            }
            finally{ 
                lock.unlock();
            }
            
            
            this.threadNumber=threadNumber;
            lochd=new LocalizationHD_(dp,200,width,height,minZ,maxZ);
        
        }
        
        public void free(){
            lochd.free();
        }
        
        public Object getMonitor(){
            return monitor;
        }
        
        public int getStreamId(){
            return lochd.getStreamId();
        }
        
        public void setBackground(Pointer device_bckg){
            lochd.setBackground(device_bckg);
            
        }
        
        public void setImage(Pointer device_image,int numframe){
            lochd.setImage(device_image,numframe);
        }
        
        public void setCoordinates(double [][] coord){
            lochd.setCoordinates(coord);
        }
        
        
        public void kill(){
            killer=true;
        }
        
         public void run(){
            synchronized(monitor) {
                monitor.notify();
                
                loop:while (true){
                
                    if (killer){//normally, it can never be true here because of synchronization
                        break loop;
                    }
                    
                    try{ 
                        monitor.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                    
                    if (killer){
                        break loop;
                    }
                    
                    
                    
                    //PROCESS HERE
                    FrameLocalization fl=lochd.optimize(iterMaxLocalization);
                    if (fl!=null){
                        stackloc.fl.add(fl);
                    }
                    
                    
                    //add avaibility
                    try{
                        lock.lock(); 
                        freeLocalization.add(threadNumber);
                    }
                    catch(Exception ee){
                        IJ.log("problem lock in LocalizationThread "+ee);//do it
                    }
                    finally{ 
                        lock.unlock();
                    }
                    
                    
                
                    
                }


            }
        }
        
        
    }
        
    
    
    



}
