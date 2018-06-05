/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;


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
import ij.process.ColorProcessor;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaError;

import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author benoit
 */
public class LocalizationPipelineMany {
     
    
    
    
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
    int nbThread=100;
    int nbPart=1;
    
    PolynomialFit [] pf=null;
    int sizePatch=32;
    
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
    
    double [][][] localMaxPosition;
    
    double [][][] image;
    
    double [] rangePredetection;
    double [][] psfPredetection;
    
    double photonThreshold=500;
    DataPhase [] dp;//liste because multiple camera possibly
    
    DataPhase [][] dpPart;
    LocalizationMany [][] loc;
    
    int iterMaxLocalization=20;
    
    
    double [][][] globalmask ;//entire stack
    double [][][] globalmask2 ;//entire stack
    
    
    double axialRange;
    double minZ;
    double maxZ;
    double stepZ;
    
    int dpnumber;
    
    int idLoc=0;
    
    double [][] scmosvariance;
    double [][] scmosgain;
    double [][] scmosoffset;
    double [][] scmosvargain;
    boolean isSCMOS;
    
    public LocalizationPipelineMany(DataPhase dp,double axialRange,int photonThreshold,int nbPart,int nbThread,int sizePatch,String path_localization,double adu, double gain,double offset){
        
        isSCMOS=false;
        this.rescale_slope=adu/gain;
        this.rescale_intercept=-offset*adu/gain;
        
        
        DataPhase [] ddp=new DataPhase[1];
        ddp[0]=dp;
        localizationPipelineMany(ddp,axialRange,photonThreshold,nbPart,nbThread,sizePatch,path_localization);
        
    }
    
    public LocalizationPipelineMany(DataPhase dp,double axialRange,int photonThreshold,int nbPart,int nbThread,int sizePatch,String path_localization,double [][] scmosvariance,double [][] scmosgain,double [][] scmosoffset,double [][] scmosvargain){
        
        
        isSCMOS=true;
        
        this.scmosvariance=scmosvariance;
        this.scmosoffset=scmosoffset;
        this.scmosgain=scmosgain;
        this.scmosvargain=scmosvargain;
        
        DataPhase [] ddp=new DataPhase[1];
        ddp[0]=dp;
        localizationPipelineMany(ddp,axialRange,photonThreshold,nbPart,nbThread,sizePatch,path_localization);
        
    }
    
    
    
    
    
    
    
    
    void localizationPipelineMany(DataPhase [] dp,double axialRange,int photonThreshold,int nbPart,int nbThread,int sizePatch,String path_localization){
        /*this.scmos=scmos;
        if (scmos!=null){
            isSCMOS=true;
        }
        else{
            isSCMOS=false;
        }
        
        this.rescale_slope=rescale_slope;
        this.rescale_intercept=rescale_intercept;*/
        
        this.photonThreshold=photonThreshold;
        this.sizePatch=sizePatch;
        this.axialRange=axialRange;
        this.nbPart=nbPart;
        this.nbThread=nbThread;
        
        this.stepZ=dp[0].param.xystep;//Same precision according to Z than for X and Y during detection
        dpnumber=dp.length;
        this.dp=dp;
        
        
        this.pf=new PolynomialFit[dpnumber];
        
        this.path_localization=path_localization;
        
        
        
        
    }
    
    
    
    
    
    
    
    public StackLocalization detectParticles(ImagePlus imp,StackLocalization stacklocinput){
        
        stackloc=stacklocinput;
        
        long timeBegin=System.currentTimeMillis();
        
        idLoc=0;
        
        
        IJ.showProgress(0);
        
        if (dp.length>1){
            IJ.log("error... dp[] should have length=1");
        }
        
        for (int i=0;i<dp.length;i++){
            double posit=dp[i].param.Zfocus;
            //if (dp[i].param.noil!=dp[i].param.nwat){removed because maybe the center is not well positioned when noil=nwat
            
            SearchPSFcenter spc = new SearchPSFcenter(dp[i],axialRange);
            posit=spc.getPosition();
            //}
            dp[i].param.ZfocusCenter=posit;
        }
        
        
        
        
        this.minZ=dp[0].param.ZfocusCenter-axialRange/2.;
        this.maxZ=dp[0].param.ZfocusCenter+axialRange/2.;
//        IJ.write("\"id\",\"frame\",\"x [nm]\",\"y [nm]\",\"z [nm]\",\"intensity [photon]\",\"background\",\"likelihood\",\"residual\",\"crlbX\",\"crlbY\",\"crlbZ\",\"predetectionX\",\"predetectionY\",\"predetectionZ\"");
        
        
        if (1!=dpnumber){
            IJ.log("error : each image source has to be associated to a retrieved phase");
            return null;
        }
        
        
        
        
        
        
        sliceNumber=imp.getNSlices();
        if (imp.getNSlices()==1){
            sliceNumber=imp.getNFrames();
        }
        
        
        this.localMaxPosition=new double [sliceNumber][][];
        this.image=new double [sliceNumber][][];
        for (int i=0;i<sliceNumber;i++){
            localMaxPosition[i]=null;
            
            image[i]=null;
            
        }
        
        
        //Initialize threads parameter
        for (int u=0;u<dp.length;u++){
            //very important to manage the end of the program (number of frames computed in parallel decreases)
            nbThread=Math.min((sliceNumber), nbThread);
        }
        
        
        nbPart=Math.min(sliceNumber/nbThread,nbPart);//manage wrong initialization compared to slice number
        
        
        
        
        
        int nbImage=imp.getNSlices();
        if (imp.getNSlices()==1){
            nbImage=imp.getNFrames();
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
        
            
         
        
        
        
        //as many dataPhase as there is available memory
        long [] memtot=new long[1];
        long [] memfree=new long[1];
        JCuda.cudaMemGetInfo(memfree, memtot);
        long memInit=memfree[0];
        DataPhase [] dptmp=new DataPhase[dp.length];
        for (int i=0;i<dp.length;i++){
            dptmp[i]=new DataPhase(dp[i],1);
            dptmp[i].setMany(this.nbThread,isSCMOS);
        }
        JCuda.cudaMemGetInfo(memfree, memtot);
        long memNeeded=(long)(1.5*(memInit-memfree[0]));//times 1.5 to not consume all memory
        
        
        
        PreDetectionThread pdt = new PreDetectionThread(imp);
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        
        nbPart=Math.min(1+(int)(memfree[0]/memNeeded),nbPart);
        
        dpPart = new DataPhase[nbPart][dp.length];
        dpPart[0]=dptmp;
        for (int up=1;up<nbPart;up++){
            for (int i=0;i<dp.length;i++){
                dpPart[up][i]=new DataPhase(dp[i],up+1);
                dpPart[up][i].setMany(this.nbThread,isSCMOS);
                //IJ.log("dpart posit "+dpPart[up][i].param.ZfocusCenter);
            }
        }
        
        
        
        int limit=Math.max(100, nbPart*nbThread+100);
        pdt.setLimit(limit);
        
        loc=new LocalizationMany[nbPart][nbThread];
        
        for (int up=0;up<nbPart;up++){
            for (int u=0;u<nbThread;u++){
                loc[up][u]=new LocalizationMany(1,dpPart[up],iterMaxLocalization,minZ-.3, maxZ+.3);
            }
        }
        
        
        
        //Killer killer = new Killer(150000);
        //killer.start();

        
        
        
                
                
            
            
            
        
        
               
        //StackParticle sp = new StackParticle(imp[0].getWidth(),imp[0].getHeight(),imp[0].getNSlices(),dp[0].param.xystep,false,0);
        
        
        
        pdt.start();
        
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
        
        computedGlobalMask=0;
        
        nbGlobalMask=Math.min(nbGlobalMask, nbImage);
        
        globalmask = new double [nbGlobalMask][width][height];
        globalmask2 = new double [nbGlobalMask][width][height];
        
        
        

        
        



        LocalizationThread [][] lt= new LocalizationThread[nbPart][nbThread];
        for (int ip=0;ip<lt.length;ip++){
            for (int i=0;i<lt[0].length;i++){
                lt[ip][i]=new LocalizationThread(ip,i);
            }
        }
        
        
        //launch step by step
        for (int ip=0;ip<lt.length;ip++){
            while(this.totNumberImagePreDetected<this.nbThread*(ip+1)){
                try{Thread.sleep(1000);}catch(Exception e){}//wait 1 second
            }
            
            for (int i=0;i<lt[0].length;i++){
                lt[ip][i].start();
            }
        }
        
        
        //the program continue after nbThread*lt.length images (all thread have readen at least one image) (important for end of task)
        
        for (int ip=0;ip<lt.length;ip++){
            for (int i=0;i<lt[0].length;i++){
                lt[ip][i].launch();
            }
        }
        
        
        
        try{
                pdt.join();
            }catch(Exception eeee){System.out.println("join pdt impossible");}
        
        for (int ip=0;ip<lt.length;ip++){
            for (int i=0;i<lt[0].length;i++){
                try{
                    lt[ip][i].join();
                }catch(Exception eeee){System.out.println("join lt impossible");}
            }
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
        
        this.save(path_localization, "Reconstruction\nCalibration path:"+dp[0].param.pathcalib+"\nLocalization path:"+path_localization+"\nReconstruction time:"+(((double)(timeEnd-timeBegin))/60000.)+" min\ndistance focus to coverslip:"+dp[0].param.Zfocus+"\nRefractive index medium:"+dp[0].param.nwat+"\nframe number:"+sliceNumber+"\nframe reconstructed:"+totNumberImageProcessed+"\nReconstruction aborted:"+IJ.escapePressed()+"\nlocalization number: "+this.idLoc+"\n");
        
        for (int up=0;up<nbPart;up++){
            for (int i=0;i<dp.length;i++){
                for (int u=0;u<100;u++){
                    JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(0));
                    for (int upt=0;upt<nbPart;upt++){
                        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(upt+1));
                    }
                }
                dpPart[up][i].free();
            }
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
                impp=ZRendering.smartColorRendering(impp,stackloc, 20, 1,true);
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
        
        ImagePlus imp;
        
        int maxRead=1000;//variable to block the program if the number of opened images is too big and use too much memory
        
        PredetectionGPU predGPU;
        PreDetectionThread(ImagePlus imp){
            
            this.imp=imp;
            
            //init with param of the first camera because predetection uses only one cam
            predGPU=new PredetectionGPU(sizeImage,dp[0],minZ, maxZ, stepZ,thresholdCrossCorrelation);
            psfPredetection=predGPU.getPSFNonNormalized();
            rangePredetection=predGPU.getRange();
            
            
            
            
            
            
        }
        
        void setLimit(int maxImageNumberInMemory){
            this.maxRead=maxImageNumberInMemory;
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
                lolo:while (imageInBuffer.size()>maxRead){
                    try{
                        Thread.sleep(10);
                        
                    }catch(Exception yy){}
                    if (IJ.escapePressed()){
                        break lolo;
                    }
                    if (z>nbImage){//test not possible
                        break loop;
                    }
                }
                if (IJ.escapePressed()){
                    IJ.log("ESCAPE BUTTON PRESSED. LOCALIZATION WILL STOP SOON.");
                    try{
                        lock.lock();
                        totNumberImagePreDetected=sliceNumber;
                        imageInBuffer.clear();
                        
                        
                    }
                    finally{
                        lock.unlock();
                    }
                    break loop;
                }
                
                ImageProcessor ip;
                if (stacked){
                    ip = ims[0].getProcessor(z);
                }
                else{
                    imp.setSlice(z);
                    ip = imp.getProcessor();
                }

                image[z-1]=new double[width][height];
                for (int ii=0,jj=xinit;ii<width;ii++,jj++){
                    for (int iii=0,jjj=yinit;iii<height;iii++,jjj++){
                            image[z-1][ii][iii]=ip.getPixelValue(jj, jjj);

                            if (isSCMOS){
                                image[z-1][ii][iii]=((image[z-1][ii][iii]-scmosoffset[ii][iii])/scmosgain[ii][iii])+scmosvargain[0][0];//here I add constant value for predetection for scmosvargain
                            }
                            else{
                                image[z-1][ii][iii]=image[z-1][ii][iii]*rescale_slope+rescale_intercept;
                            }
                            if (image[z-1][ii][iii]<0){
                                image[z-1][ii][iii]=0;
                            }
                    }
                }


                    


                //ImageShow.imshow(image[0][z-1]," im ");
                //predetection on channel 0

                predGPU.setImage(image[z-1]);
                
                
                localMaxPosition[z-1]=predGPU.convolveNormalizedInFourierDomain();

                
                    
                


                try{
                    lock.lock();
                    
                    if (isSCMOS){
                        for (int ii=0,jj=xinit;ii<width;ii++,jj++){
                            for (int iii=0,jjj=yinit;iii<height;iii++,jjj++){
                                image[z-1][ii][iii]=ip.getPixelValue(jj, jjj);
                                //do it again but not with constant value now for scmosvargain
                                image[z-1][ii][iii]=((image[z-1][ii][iii]-scmosoffset[ii][iii])/scmosgain[ii][iii])+scmosvargain[ii][iii];
                                
                            }
                        }
                    }
                    
                    imageInBuffer.add(z-1);
                    totNumberImagePreDetected++;
                }
                finally{
                    lock.unlock();
                }
                ip=null;

            }
            predetectionFinished=true;
            predGPU.free();
        }
        
        
        
        
    }
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class LocalizationThread extends Thread{
        
        
        boolean canBeLaunch=false;
        boolean stopReading=false;
        int [] positionPointer = new int[1];
        
        double [] psf;
        
        int imageNumber=0;
        
        
        boolean [][] mask ;
        double [][] patch ;
        double [][] patchscmos ;
        
        int threadNumber;
        int partNumber;
        
        FrameLocalization fl ;
        
        LocalizationThread(int partNumber,int threadNumber){
            
            mask = new boolean [width][height];
            
            patch =new double [dp.length][sizePatch*sizePatch];
            if (isSCMOS){
                patchscmos=new double [dp.length][sizePatch*sizePatch];
            }
            this.psf=new double[sizePatch*sizePatch];
            
            this.threadNumber=threadNumber;
            this.partNumber=partNumber;
        }
        
        public void launch(){
            
            canBeLaunch=true;
            
            
        }
        
        
        public void run(){
            
            mainlooper:while (totNumberImageProcessed<sliceNumber){
                
                
                if (IJ.escapePressed()){
                    for (int u=0;u<dp.length;u++){

                        dpPart[partNumber][u].modelMany.decrementNumberPSFToCompute();
                    }
                    stopReading=true;
                    break mainlooper;
                }
                
                
                
                boolean imageReaden=false;
                
                //block until all threads have readen 1 image
                //read image
                //IJ.log("begin lock readimage "+partNumber);
                try{
                    lock.lock();
                    if((imageInBuffer.size()>0)&&(!stopReading)){
                        imageReaden=true;
                        imageNumber=(int)imageInBuffer.get(0);
                        imageInBuffer.remove(0);
                    }
                }
                catch(Exception ee){
                    IJ.log("problem lock 2 in class localizationPipelineMany");//do it
                }
                finally{
                    lock.unlock();
                }
                

                //IJ.log("end lock readimage "+partNumber);
                
                
                if (imageReaden){
                    
                    fl=new FrameLocalization(imageNumber);
                    
                    for (int j=0;j<width;j++){
                        for (int jj=0;jj<height;jj++){
                            mask[j][jj]=false;
                        }
                    }

                    double [][] vect = localMaxPosition[imageNumber];
                    
//                    for (int i=0;i<Math.min(vect.length,10);i++){
//                        IJ.log("vect "+vect[i][0]+"   "+vect[i][1]+"  "+imageNumber);
//                    }
//                    
//                    IJ.log("vect length "+vect.length+"   "+nbMaxParticle);

                    
                    
                    
                    //IJ.log("maximum found: "+vect[0][0]+"/"+vect[0][2]+"/"+vect[0][3]+"   "+vect[1][0]+"/"+vect[1][2]+"   "+vect[2][0]+"/"+vect[2][2]);
                    
                    
                    //IJ.log("vect "+vect.length);
                    
                    
                    positionPointer[0]=0;
                    int numb=0;
                    loop:for (int n=0;n<vect.length;){
                        
                        if (IJ.escapePressed()){
                            break loop;
                        }
//                        int [] maxsP = myGetMaxPositions(vect,mask,positionPointer,width,height);
                        
                        int maxsP = myGetMaxPositions(vect,mask,positionPointer,width,height);
                        

//                        if (maxs!=null){
                        if (maxsP>=0){
                            int [] maxs=new int[4];
                            maxs[0]=(int)Math.round(vect[maxsP][0]);
                            maxs[1]=(int)Math.round(vect[maxsP][1]);
                            maxs[2]=(int)Math.round(vect[maxsP][2]);
                            maxs[3]=(int)Math.round(vect[maxsP][3]);
                            //IJ.log("maxs "+maxs[0]+"   "+maxs[1]+"  "+stream);
                            //maxs[2]=maxs[2]-1;
                            //IJ.log("position init : "+maxs[0]+"   "+maxs[1]+"   "+maxs[2]+"  "+maxs[3]+"   next Pointer:"+positionPointer[0]);
                            
                            boolean alreadyTreated=false;
                            loop1:for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
                                for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
                                    //for (int c=0;c<dp.length;c++){//channel
                                        if (mask[maxs[2]+j][maxs[3]+jj]==true){
                                            alreadyTreated=true;
                                            break loop1;
                                        }
                                        
                                        patch[0][i*sizePatch+ii]=image[imageNumber][maxs[2]+j][maxs[3]+jj];
                                        if (isSCMOS){
                                            patchscmos[0][i*sizePatch+ii]=scmosvargain[maxs[2]+j][maxs[3]+jj];
                                        }
                                        
                                    //}
                                }
                            }
                            if (alreadyTreated){
                                continue loop;
                            }
                            
                            
                            
                            //for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
                            //    for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
                            //        mask[imageNumber][maxs[2]+j][maxs[3]+jj]=true;//has been treated
                            //    }
                            //}

                            for (int u=0;u<psf.length;u++){
                                psf[u]=psfPredetection[maxs[1]][u];
                            }
                            //dp[0].psf.imshow(sizePatch, psf, "psfInit");



                            double [] AandB;
                            
                            
                            if (isSCMOS){
                                AandB=getAandB(patch[0], psf,patchscmos[0], .1);//compute A and B using threshold 10 percent of max of psf for channel 0
                            }
                            else{
                                AandB=getAandB(patch[0], psf, .1);
                            }
                            

                            //double snr=getSNR(patch[stream], psf);

                            //double photonBackgroundRatio=AandB[0]/AandB[1];
                            //IJ.log("A:"+AandB[0]+"   B:"+AandB[1]+"   "+photonBackgroundRatio);

                            boolean partTreated=false;


                            
                            //IJ.log("a&B "+AandB[0]+"  "+AandB[1]+"  "+AandB[2]);

                            if ((AandB[2]/AandB[1]>maxSNR)&&(AandB[1]>0)){



                                //IJ.log(""+imageNumber+"  "+AandB[2]+"  "+n+"   "+AandB[0]+"   "+AandB[1]);
                                
                                //if (AandB[0]>photonThreshold/2){
                                if (true){

                                /////////////////////////////
//                                for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
//                                    for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
//                                        //HERE ADD MODELcomparison
//                                        if (imageNumber<nbGlobalMask){
//                                            globalmask[imageNumber][maxs[2]+j][maxs[3]+jj]=psf[i*sizePatch+ii];
//                                            globalmask2[imageNumber][maxs[2]+j][maxs[3]+jj]=patch[0][i*sizePatch+ii];
//                                        }
//                                    }
//                                }
                                //globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=rangePredetection[maxs[1]];
                                /////////////////////////////

                                    //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda locPipeline loc 0 "+cudaResult);}
                                    
                                    
                                    if (isSCMOS){
                                        loc[partNumber][threadNumber].setSubWindowScmos(patch,patchscmos);
                                    }
                                    else{
                                        loc[partNumber][threadNumber].setSubWindow(patch);
                                    }

                                    //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda locPipeline loc 1 "+cudaResult);}
                                    
                                    
                                    //loc[partNumber][threadNumber].init(AandB[0], AandB[1], -vect[maxsP][5], -vect[maxsP][6], rangePredetection[maxs[1]]+vect[maxsP][4]);
                                    //IJ.log("A/B "+AandB[0]+"  "+AandB[1]);
                                    loc[partNumber][threadNumber].init(AandB[0], AandB[1], -vect[maxsP][5], -vect[maxsP][6], rangePredetection[maxs[1]]+vect[maxsP][4]);
                                    
                                    
                                    //moyOffset+=(-vect[maxsP][5]);
                                    
                                    //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda locPipeline loc 2 "+cudaResult);}

                                    //try{
                                    //    Thread.sleep(1000);
                                    //}catch(Exception egfregk){}
                                    boolean locate=true;
                                    //if (true)
                                    
                                    
                                    
                                    //ImageShow.imshow(loc[partNumber][threadNumber].mytest(),"mod "+rangePredetection[maxs[1]]);
                                    //double [] mod;
                                    //for (double zz=-1;zz<=10;zz+=1){
                                        //loc[partNumber][threadNumber].init(AandB[0], AandB[1], 0, 0, zz);
                                        //ImageShow.imshow(loc[partNumber][threadNumber].mytest(),"mod "+zz+"   "+rangePredetection[maxs[1]]);
                                    //}
                                    
                                    
                                    locate=loc[partNumber][threadNumber].localize();
                                    
                                    //else
                                    //    loc[partNumber][threadNumber].finishtest();
                                    
                                    if ((true)&&(locate)){



                                        

                                        double [] model=loc[partNumber][threadNumber].getPSF(0);//get model for channel 0
                                        //dp[0].psf.imshow(sizePatch, model, "modelFound");

                                        double maxModelThreshold=0;
                                        double background=loc[partNumber][threadNumber].getB()[0];
                                        for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
                                            for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
                                                if (maxModelThreshold<(model[i*sizePatch+ii]-background)){
                                                    maxModelThreshold=(model[i*sizePatch+ii]-background);
                                                }
                                            }
                                        }
                                        maxModelThreshold*=.2;
//                                        double std=validLocalization(patch[0],model,maxModelThreshold);

                                        //IJ.log("photons "+loc[threadNumber].getA()[0][0]+" / Th:"+photonThreshold);

                                        //IJ.log("photon : "+loc[partNumber][threadNumber].getA()[0][0]);
                                        if (loc[partNumber][threadNumber].getA()[0][0]>photonThreshold){
                                            
                                                //only displacement < than 3 pixels from the init
                                            if ((Math.abs(loc[partNumber][threadNumber].getX())<10*dpPart[partNumber][0].param.xystep)&&(Math.abs(loc[partNumber][threadNumber].getY())<10*dp[0].param.xystep)){
                                                int px_pix=(int)(maxs[2]-(loc[partNumber][threadNumber].getX()/dpPart[partNumber][0].param.xystep));
                                                int py_pix=(int)(maxs[3]-(loc[partNumber][threadNumber].getY()/dpPart[partNumber][0].param.xystep));
                                                
                                                
                                                if ((px_pix>=0)&&(px_pix<mask.length)&&(py_pix>=0)&&(py_pix<mask[0].length)&&(mask[px_pix][py_pix]==false)){
                                                    
                                                    
                                                    if ((loc[partNumber][threadNumber].getZ()>minZ)&&(loc[partNumber][threadNumber].getZ()<maxZ)){
                                                        
                                                        numb++;
                                                        n++;
                                                        
                                                        //IJ.log("zz "+rangePredetection[maxs[1]]+"  "+loc[partNumber][threadNumber].getZ());
                                                        double my_x=1000*(maxs[2]*dpPart[partNumber][0].param.xystep-loc[partNumber][threadNumber].getX());
                                                        double my_y=1000*(maxs[3]*dpPart[partNumber][0].param.xystep-loc[partNumber][threadNumber].getY());
                                                        double my_z=(1000*loc[partNumber][threadNumber].getZ());
                                                        double my_A=loc[partNumber][threadNumber].getA()[0][0];
                                                        double my_B=loc[partNumber][threadNumber].getB()[0];                                                        //double my_Score=(loc[partNumber][threadNumber].getGlobalLikelihood());
                                                        //double my_Score=(loc[partNumber][threadNumber].getGlobalLikelihood());
                                                        double my_Score= scoreCompute(patch[0],model,maxModelThreshold,background);
                                                        double my_crlbx=(1000*loc[partNumber][threadNumber].getCRLBX());
                                                        double my_crlby=(1000*loc[partNumber][threadNumber].getCRLBY());
                                                        double my_crlbz=(1000*loc[partNumber][threadNumber].getCRLBZ());

                                                        /*double spx=Math.floor(100* my_x)/100;

                                                        double spy=Math.floor(100* my_y)/100;

                                                        double spz=Math.floor(100* my_z)/100;

                                                        double spa=Math.floor(100* my_A)/100;

                                                        double spb=Math.floor(100* my_B)/100;

                                                        double spl=Math.floor(100* my_Score)/100;*/

    //                                                    double crlbX=0;
    //
    //                                                    double crlbY=0;
    //
    //                                                    double crlbZ=0;
    //                                                    double crlbX=loc[partNumber][threadNumber].getCRLBX(0);
    //        
    //                                                    double crlbY=loc[partNumber][threadNumber].getCRLBY(0);
    //        
    //                                                    double crlbZ=loc[partNumber][threadNumber].getCRLBZ(0);

//                                                        double crlbX=Math.floor(100* (1000*loc[partNumber][threadNumber].getCRLBX(0)))/100;
//
//                                                        double crlbY=Math.floor(100* (1000*loc[partNumber][threadNumber].getCRLBY(0)))/100;
//
//                                                        double crlbZ=Math.floor(100* (1000*loc[partNumber][threadNumber].getCRLBZ(0)))/100;
//
//                                                        double stdshow=Math.floor(100* (std))/100;
//
//                                                        double predetection = Math.floor(100*rangePredetection[maxs[1]])/100;
//                                                        IJ.write((idLoc)+","+(imageNumber)+","+(spx)+","+(spy)+","+(spz)+","+(spa)+","+(spb)+","+(spl)+","+(stdshow)+","+crlbX+","+crlbY+","+crlbZ+","+maxs[2]+","+maxs[3]+","+predetection);

                                                        PLocalization p = new PLocalization(idLoc,imageNumber,my_x,my_y,my_z,my_A,my_B,my_Score,my_crlbx,my_crlby,my_crlbz);
                                                        fl.loc.add(p); 
                                                        
                                                        
                                                         
                                                        //replacement in mask of the region corresponting to 20% of max PSF
                                                        try{
                                                        lock.lock();
                                                            idLoc++;
                                                        }
                                                        finally{
                                                            lock.unlock();
                                                        }

                                                        for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
                                                            for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
                                                                //HERE ADD MODELcomparison
                                                                if ((model[i*sizePatch+ii]-background)>maxModelThreshold){
                                                                    mask[maxs[2]+j][maxs[3]+jj]=true;//has been treated

                                                                }
                                                                if (imageNumber<nbGlobalMask){
                                                                    globalmask[imageNumber][maxs[2]+j][maxs[3]+jj]=model[i*sizePatch+ii];
                                                                    globalmask2[imageNumber][maxs[2]+j][maxs[3]+jj]=patch[0][i*sizePatch+ii];
                                                                }
                                                            }
                                                        }
                                                        if (imageNumber<nbGlobalMask){
                                                            //globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=AandB[0];
                                                            globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=my_Score;
                                                            //globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=rangePredetection[maxs[1]];
                                                            //globalmask[imageNumber][maxs[2]+1][maxs[3]+0]=my_z;
                                                        }
                                                        partTreated=true;
                                                    }
                                                }



                                            }
                                            
                                        }
                                        else{
                                            //IJ.log("low photon");
                                        }

                                    }

                                }
                                else{
                                    //IJ.log("part reject "+photonThreshold+"  "+AandB[0]);
                                }

                                   
                            }
                            if (partTreated==false){
                                //if non treated, delete 50% of max in patch to do not process it again
                                for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
                                    for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
                                        //HERE ADD MODELcomparison
                                        if ((patch[0][i*sizePatch+ii]>AandB[2]*.5)){
                                            mask[maxs[2]+j][maxs[3]+jj]=true;//has been treated

                                        }
                                    }
                                }
                            }


                            //in patch : compute RSB to reject or not
                            maxs=null;
                        }
                        else{
                            //no more particle
                            //IJ.log("break...");
                            break loop;
                        }

                        //iptest.putPixelValue(maxs[1], maxs[2], (float)maxs[0]);
                        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda locPipeline last 0 "+cudaResult);}
                    }
                    
                    
                    
                    while (!canBeLaunch){//wait end of first pass to continue
                        try{
                            Thread.sleep(1);
                        }
                        catch(Exception ee){
                            IJ.log("problem lock 2 in class localizationPipelineMany");//do it
                        }
                    }
                    
                    
                    //IJ.log("begin lock finishimage "+partNumber);
                    
                    try{
                        lock.lock();
                            image[imageNumber]=null;
                            
                        
                        stackloc.fl.add(fl);
                        localMaxPosition[imageNumber]=null;
                        
                        
                        
                        totNumberImageProcessed++;
                        
                        //IJ.log("to im process... "+totNumberImageProcessed+"  "+sliceNumber+"  "+(sliceNumber-totNumberImageProcessed)+"  id"+partNumber);
                        IJ.showProgress((double)totNumberImageProcessed/(double)(sliceNumber+1));
                        
                        if((sliceNumber-totNumberImageProcessed)<nbThread*nbPart){
                            
                            for (int u=0;u<dp.length;u++){
                                
                                dpPart[partNumber][u].modelMany.decrementNumberPSFToCompute();
                            }
                            stopReading=true;
                        }
                        
                        
                        if (imageNumber<nbGlobalMask){
                            computedGlobalMask++;
                            //IJ.log("comp "+computedGlobalMask);
                            if ( computedGlobalMask == (nbGlobalMask) ){
                                //ImageShow.imshow(globalmask,"model_frame1-100");
                                //ImageShow.imshow(globalmask2,"image_frame1-100");
                                pr.launch_();
                                globalmask=null;
                                globalmask2=null;
                            }
                        }
                        
                        
                    }
                    catch(Exception ee){
                        IJ.log("problem lock 3 in class localizationPipelineMany "+ee);//do it
                    }
                    finally{ 
                        lock.unlock();
                    }
                    
                    //IJ.log("end lock finishimage "+partNumber);
                    IJ.showStatus("image:"+imageNumber+" / localization #:"+numb);
                    //IJ.log(""+numb+" loc treated for image "+imageNumber+"   "+totNumberImageProcessed+"/"+totNumberImagePreDetected+"   "+partNumber);
                    
//                    for (int u=0;u<dp.length;u++){
//                        //very important to manage the end of the program (number of frames computed in parallel decreases)
//                        dpPart[partNumber][u].modelMany.updateNumberPSFToCompute(Math.min((sliceNumber-totNumberImageProcessed), nbThread));
//                    }
                }
                else{
                    try{//IJ.log("j attends mon image "+partNumber+"    imInBuffer"+imageInBuffer.size());
                        Thread.sleep(0);
                        
                    }
                    catch(Exception ee){
                        IJ.log("problem lock 2 in class localizationPipelineMany");//do it
                    }
                }
            }
            
        }
        
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda locPipeline last "+cudaResult);}

    }
        

        double validLocalization(double [] patch, double [] model,double maxModelThreshold){
            
            
            
            double sum=0;
            double count=0;
            for (int i=0;i<patch.length;i++){
                if (patch[i]>maxModelThreshold){
                    sum+=(patch[i]-model[i])*(patch[i]-model[i])/model[i];
                    count++;
                }
                
            }
            sum/=count;
            return sum;
        }
        
        
        public double scoreCompute(double [] patch, double [] model,double maxModelThreshold, double background){
            
            
            
            double sum=0;
            double sumModel=0;
            for (int i=0;i<patch.length;i++){
                //if (model[i]-background>maxModelThreshold){
                    sum+=(patch[i]-model[i])*(patch[i]-model[i])/model[i];
                    sumModel++;
                //}
                
            }
            sum/=sumModel;
            return sum;
        }
        
        
        /*public double scoreComputeProba(double [] patch, double [] model,double maxModelThreshold, double background){
            
            
            
            double sum=0;
            for (int i=0;i<patch.length;i++){
                //if (model[i]-background>maxModelThreshold){
                    sum+=(patch[i]-model[i])*(patch[i]-model[i])/model[i];
                    
                //}
                
            }
            return KhiTest.chisqr(patch.length*patch.length-5, sum);
        }*/
        
        
        
        

        private double [] getAandB(double [] patch, double [] psf, double thresholdPercentage){
            double maxPatch=Double.NEGATIVE_INFINITY;
            double maxi=Double.NEGATIVE_INFINITY;
            for (int i=0;i<patch.length;i++){
                if (psf[i]>maxi){
                    maxi=psf[i];
                }
                if (patch[i]>maxPatch){
                    maxPatch=patch[i];
                }
            }
            //IJ.log("maxi "+maxi);
            double threshold=maxi*thresholdPercentage;
            double countB=0;
            double B=0;
            double A=0;
            double countA=0;
            for (int i=0;i<patch.length;i++){
                    if (psf[i]<threshold){
                        B+=patch[i];
                        countB++;
                    }
                    else{
                        A+=patch[i];
                        countA++;
                    }
            }

            if (countB!=0){
                B/=countB;
            }
            else{
                B=0;
            }
            //IJ.log("count A "+countA+"  "+A+"  "+B);
            A-=countA*B;
            double [] res = {A,B,maxPatch};//last element=maxSNR
            return res;

        }

        
        
        

        private double [] getAandB(double [] patch, double [] psf, double [] scmos, double thresholdPercentage){
            double maxPatch=Double.NEGATIVE_INFINITY;
            double maxi=Double.NEGATIVE_INFINITY;
            for (int i=0;i<patch.length;i++){
                if (psf[i]>maxi){
                    maxi=psf[i];
                }
                if (patch[i]-scmos[i]>maxPatch){
                    maxPatch=patch[i]-scmos[i];
                }
            }
            //IJ.log("maxi "+maxi);
            double threshold=maxi*thresholdPercentage;
            double countB=0;
            double B=0;
            double A=0;
            double countA=0;
            for (int i=0;i<patch.length;i++){
                    if (psf[i]<threshold){
                        B+=patch[i]-scmos[i];
                        countB++;
                    }
                    else{
                        A+=patch[i]-scmos[i];
                        countA++;
                    }
            }

            if (countB!=0){
                B/=countB;
            }
            else{
                B=0;
            }
            //IJ.log("count A "+countA+"  "+A+"  "+B);
            A-=countA*B;
            double [] res = {A,B,maxPatch};//last element=maxSNR
            return res;

        }


        

        
        
        public int  myGetMaxPositions(double [][] v, boolean [][] mask,int [] pointerPosition,int width, int height)
        {
//            int [] results = new int [4];
            boolean found=false;
            int getPosit=-1;
            searchListSorted:for (int i=pointerPosition[0];(i<v.length);i++)
            {
                
                
                if (!mask[(int)Math.round(v[i][2])][(int)Math.round(v[i][3])])
                {
                    if (((int)Math.round(v[i][2])>sizePatch/2+1)&&((int)Math.round(v[i][3])>sizePatch/2+1)&&((int)Math.round(v[i][2])<width-sizePatch/2-1)&&((int)Math.round(v[i][3])<height-sizePatch/2-1)){
//                        results[0]=(int)Math.round(v[i][0]);
//                        results[1]=(int)Math.round(v[i][1]);
//                        results[2]=(int)Math.round(v[i][2]);
//                        results[3]=(int)Math.round(v[i][3]);
                        getPosit=i;
                        pointerPosition[0]=i+1;
                        found=true;
                        break searchListSorted;
                    }
                }
                
            }
            
            return getPosit;
//            if (found){
//                return results;		
//            }
//            else{
//                return null;
//            }
        }



        
        
    
    
        void free(){
            
            
            for (int u=0;u<dp.length;u++){
                //dp[u].free();
            }
            //MyCudaStream.destroy();
        }
    
    
        
    public void save(String path,String content){
        
        if (path!=null){
            try {
                int nn=(path.lastIndexOf("."));
                if (nn>path.length()-5){
                    path=path.substring(0,nn);
                }
                path=path+"_ZolaProtocol.txt";

                PrintWriter sortie;


                    sortie = new PrintWriter(new FileWriter(path, false));

                    sortie.println(""+content);
                            
                    

                    sortie.close();




            } catch (Exception e) {
                    e.printStackTrace();
                    IJ.log("error saving process");
            }
        }
        
    }
    
    
}
