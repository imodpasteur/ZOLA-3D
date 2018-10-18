/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.postprocess.ZRendering;
/**
 *
 * @author benoit
 */
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author benoit
 */
public class LocalizationPipeline {
     
    boolean show=false;
    String path_localization;
    
    StackLocalization stackloc;
            
    double rescale_slope; 
    double rescale_intercept;
    
    PlotResult pr ;
    ReentrantLock lockPred = new ReentrantLock();
    ReentrantLock lock = new ReentrantLock();
    boolean predetectionFinished=false;
    int totNumberImageProcessed=0;
    int totNumberImagePreDetected=0;
    ArrayList<Integer> imageInBuffer = new ArrayList<Integer>();
    int nbThread=Runtime.getRuntime().availableProcessors();
    
    PolynomialFit [] pf=null;
    int sizePatch=32;
    
    int xinit;
    int yinit;
    int width;
    int height;
    
    int sliceNumber;
    double thresholdCrossCorrelation=.1;
    
    int nbGlobalMask=100;
    int computedGlobalMask=0;
    double maxSNR=1; 
    
    double [][][] localMaxPosition;
    
    double [][][] image;
    
    double [] rangePredetection;
    double [][][] psfPredetection;
    
    double photonThreshold=500;
    DataPhase dp;//liste because multiple camera possibly
    
    Localization [] loc;
    
    int iterMaxLocalization=20;
    
    
    double [][][] globalmask=null ;//entire stack
    double [][][] globalmask2=null ;//entire stack
    
    
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
    
    public LocalizationPipeline(DataPhase dp,double axialRange,int photonThreshold,int sizePatch,String path_localization,double adu, double gain,double offset,boolean show){
    
        isSCMOS=false;
        this.rescale_slope=adu/gain;
        this.rescale_intercept=-offset*adu/gain;
        
        localizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,show);
        
    }
    
    public LocalizationPipeline(DataPhase dp,double axialRange,int photonThreshold,int sizePatch,String path_localization,double [][] scmosvariance,double [][] scmosgain,double [][] scmosoffset,double [][] scmosvargain,boolean show){
    
        isSCMOS=true;
        
        this.scmosvariance=scmosvariance;
        this.scmosoffset=scmosoffset;
        this.scmosgain=scmosgain;
        this.scmosvargain=scmosvargain;
        
        
        localizationPipeline(dp,axialRange,photonThreshold,sizePatch,path_localization,show);
        
    }
    
    void localizationPipeline(DataPhase dp,double axialRange,int photonThreshold,int sizePatch,String path_localization,boolean show){
        
        
        this.show=show;
        this.sizePatch=sizePatch;
        this.photonThreshold=photonThreshold;
        this.axialRange=axialRange;
        this.stepZ=dp.param.xystep;//Same precision according to Z than for X and Y during detection
        this.dp=dp;
        
        
        dpnumber=1;
        
          
        this.path_localization=path_localization;
        
    }
    
    
    
    
    
    
    public StackLocalization detectParticles(ImagePlus imp,StackLocalization stacklocinput){
        
       
        stackloc=stacklocinput;
        
        long timeBegin=System.currentTimeMillis();
        
        idLoc=0;
        
        
        IJ.showProgress(0);
        
        
        
        double posit=dp.param.Zfocus;
        //if (dp[i].param.noil!=dp[i].param.nwat){removed because maybe the center is not well positioned when noil=nwat

        SearchPSFcenter spc = new SearchPSFcenter(dp,axialRange);
        posit=spc.getPosition();
        //}
        dp.param.ZfocusCenter=posit;
        
        
        
        
        
        this.minZ=dp.param.ZfocusCenter-axialRange/2.;
        this.maxZ=dp.param.ZfocusCenter+axialRange/2.;
//        IJ.write("\"id\",\"frame\",\"x [nm]\",\"y [nm]\",\"z [nm]\",\"intensity [photon]\",\"background\",\"likelihood\",\"residual\",\"crlbX\",\"crlbY\",\"crlbZ\",\"predetectionX\",\"predetectionY\",\"predetectionZ\"");
        
        
        
        
        
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
        int cores = Runtime.getRuntime().availableProcessors();
        nbThread=cores;
        
        nbThread=Math.min((sliceNumber), nbThread);
        
        
        dp.setMany(nbThread);
        
        loc=new Localization[nbThread];
        for (int u=0;u<nbThread;u++){
            loc[u]=new Localization(1,dp,iterMaxLocalization,minZ-.3, maxZ+.3,u);
        }
        
        
        
        
        
        
        
        //Killer killer = new Killer(150000);
        //killer.start();

        
        
        
                
                
        
        
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
        
        
        
        
            
            
            
        
        
               
        //StackParticle sp = new StackParticle(imp[0].getWidth(),imp[0].getHeight(),imp[0].getNSlices(),dp[0].param.xystep,false,0);
        
        
        int limit=nbThread+100;
        PreDetectionThread [] pdt = new PreDetectionThread[nbThread];
        for (int i=0;i<nbThread;i++){
            pdt[i] = new PreDetectionThread(imp,limit,i,nbThread);
            pdt[i].start();
            
        }
        
        if (!show){
            pr = new PlotResult(1,width,height);//update every 30 second
        }
        
        boolean saveOnTheFly=false;
        if (path_localization!=null){
            if (path_localization.length()>2){
                stackloc.saveOnTheFly(path_localization);
                saveOnTheFly=true;
            }
        }
        
        
        
        
        //Killer killer = new Killer(150000);
        
        
        //GarbageLauncher gcl = new GarbageLauncher(1800000);//system gc every 30 min
        //gcl.start();
        
        
        computedGlobalMask=0;
        
        nbGlobalMask=Math.min(nbGlobalMask, nbImage);
        
        if (show){
            globalmask = new double [nbGlobalMask][width][height];
            globalmask2 = new double [nbGlobalMask][width][height];
        }
        
        
        

        
        


        
        LocalizationThread [] lt= new LocalizationThread[nbThread];
        for (int ip=0;ip<nbThread;ip++){
            lt[ip]=new LocalizationThread(ip);
            lt[ip].start();
            lt[ip].launch();
        }
        
        
        
        for (int ip=0;ip<nbThread;ip++){
            try{
                pdt[ip].join();
            }catch(Exception eeee){System.out.println("join pdt impossible");}
        }
        
        
        for (int ip=0;ip<lt.length;ip++){
            try{
                lt[ip].join();
            }catch(Exception eeee){System.out.println("join lt impossible");}
            
        }
        
        
        
        if (!show){
            pr.stop_();
        }
        
        if (saveOnTheFly){
            stackloc.stopSaveOnTheFly();
        }
        
        //gcl.stopRun();
        
        
        
        
        long timeEnd=System.currentTimeMillis();
        
        IJ.showProgress(0);
        
        
        
        IJ.log("elapsed time = "+((double)(timeEnd-timeBegin))/60000.+" min");
        IJ.log("localization number: "+this.idLoc);
        if (path_localization.length()>3){
            this.save(path_localization, "Reconstruction\nCalibration path:"+dp.param.pathcalib+"\nLocalization path:"+path_localization+"\nReconstruction time:"+(((double)(timeEnd-timeBegin))/60000.)+" min\ndistance focus to coverslip:"+dp.param.Zfocus+"\nRefractive index medium:"+dp.param.nwat+"\nframe number:"+sliceNumber+"\nframe reconstructed:"+totNumberImageProcessed+"\nReconstruction aborted:"+IJ.escapePressed()+"\nlocalization number: "+this.idLoc+"\npatch size: "+sizePatch+"\n");
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
        
        ImagePlus imp;
        int nbThread;
        int idThread;
        int maxRead=1000;//variable to block the program if the number of opened images is too big and use too much memory
        
        Predetection pred;
        PreDetectionThread(ImagePlus imp, int maxImageNumberInMemory,int idThread,int nbThread){
            this.nbThread=nbThread;
            this.idThread=idThread;
            this.imp=imp;
            
            this.maxRead=maxImageNumberInMemory;
            
            //init with param of the first camera because predetection uses only one cam
            
            pred=new Predetection(width,height,dp,minZ, maxZ, stepZ,thresholdCrossCorrelation);
            if (idThread==0){
                psfPredetection=pred.getPSFNonNormalized();
                rangePredetection=pred.getRange();
            }
            
        }
        
        
        
        
        
        public void run(){
            
            
            boolean stacked=true;
            int nbImage=imp.getNSlices();
            if (imp.getNSlices()==1){
                nbImage=imp.getNFrames();
                stacked=false;
            }
            ImageStack ims=null;
            if (stacked){
                ims=imp.getStack();
                
            }
            
            
            loop:for (int z=idThread+1;z<=nbImage;z+=nbThread){
                
                //IJ.log("image detection: "+(z-1));
                
                
                try{
                    lockPred.lock();
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
                }
                finally{
                    lockPred.unlock();
                }
                
                if (IJ.escapePressed()){
                    IJ.log("ESCAPE BUTTON PRESSED. LOCALIZATION WILL STOP SOON.");
                    try{
                        lockPred.lock();
                        imageInBuffer.clear();
                    }
                    finally{
                        lockPred.unlock();
                    }
                    break loop;
                }
                
                //IJ.log("process z "+z);
                   
                ImageProcessor ip;
                if (stacked){
                    ip = ims.getProcessor(z);
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

                pred.setImage(image[z-1]);
                
                
                localMaxPosition[z-1]=pred.convolveNormalizedInFourierDomain();

                
                    
                


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
        }
        
        
        
        
    }
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class LocalizationThread extends Thread{
        
        
        boolean canBeLaunch=false;
        boolean stopReading=false;
        int [] positionPointer = new int[1];
        
        double [][] psf;
        
        int imageNumber=0;
        
        
        boolean [][] mask ;
        double [][] patch ;
        double [][] patchscmos ;
        int threadNumber;
        
        FrameLocalization fl ;
        
        LocalizationThread(int threadNumber){
            
            mask = new boolean [width][height];
            
            patch =new double [sizePatch][sizePatch];
            if (isSCMOS){
                patchscmos=new double [sizePatch][sizePatch];
            }
            this.psf=new double[sizePatch][sizePatch];
            
            this.threadNumber=threadNumber;
        }
        
        public void launch(){
            
            canBeLaunch=true;
            
            
        }
        
        
        public void run(){
            
            mainlooper:while (totNumberImageProcessed<sliceNumber){
                
                
                if (IJ.escapePressed()){

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
                

                
                if (imageReaden){
                    
                    fl=new FrameLocalization(imageNumber);
                    
                    for (int j=0;j<width;j++){
                        for (int jj=0;jj<height;jj++){
                            mask[j][jj]=false;
                        }
                    }

                    double [][] vect = localMaxPosition[imageNumber];
                    
                    
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

                                    if (mask[maxs[2]+j][maxs[3]+jj]==true){
                                        alreadyTreated=true;
                                        break loop1;
                                    }
                                    patch[i][ii]=image[imageNumber][maxs[2]+j][maxs[3]+jj];
                                    if (isSCMOS){
                                        patchscmos[i][ii]=scmosvargain[maxs[2]+j][maxs[3]+jj];
                                    }
                                }
                            }
                            if (alreadyTreated){
                                continue loop;
                            }
                            
                            
                            

                            for (int u=0;u<psf.length;u++){
                                psf[u]=psfPredetection[maxs[1]][u];
                            }
                            //dp[0].psf.imshow(sizePatch, psf, "psfInit");



                            double [] AandB;
                            if (isSCMOS){
                                AandB=getAandB(patch, psf,patchscmos, .1);//compute A and B using threshold 10 percent of max of psf for channel 0
                            }
                            else{
                                AandB=getAandB(patch, psf, .1);//compute A and B using threshold 10 percent of max of psf for channel 0
                            }
                            //double snr=getSNR(patch[stream], psf);

                            //double photonBackgroundRatio=AandB[0]/AandB[1];
                            //IJ.log("A:"+AandB[0]+"   B:"+AandB[1]+"   "+photonBackgroundRatio);

                            boolean partTreated=false;


                            

                            if ((AandB[2]/AandB[1]>maxSNR)&&(AandB[1]>0)){


                                if (true){

                                
                                    
                                    
                                    loc[threadNumber].setSubWindow(0,0,patch);
                                    
                                    if (isSCMOS){
                                        loc[threadNumber].setSubWindowScmos(0,0,patchscmos);
                                    }
                                    
                                    loc[threadNumber].init(AandB[0], AandB[1], -vect[maxsP][5], -vect[maxsP][6], rangePredetection[maxs[1]]+vect[maxsP][4]);
                                    
                                    
                                    
                                    
                                    
                                    boolean locate=locate=loc[threadNumber].localize();
                                    
                                    
                                    
                                    if ((true)&&(locate)){



                                        

                                        double [][] model=loc[threadNumber].getPSF(0,0);//get model for channel 0 and cam 0
                                        //dp[0].psf.imshow(sizePatch, model, "modelFound");

                                        double maxModelThreshold=0;
                                        double background=loc[threadNumber].getB()[0];
                                        for (int i=0,j=-sizePatch/2;i<sizePatch;i++,j++){
                                            for (int ii=0,jj=-sizePatch/2;ii<sizePatch;ii++,jj++){
                                                if (maxModelThreshold<(model[i][ii]-background)){
                                                    maxModelThreshold=(model[i][ii]-background);
                                                }
                                            }
                                        }
                                        maxModelThreshold*=.2;
//                                        double std=validLocalization(patch[0],model,maxModelThreshold);

                                        //IJ.log("photons "+loc[threadNumber].getA()[0][0]+" / Th:"+photonThreshold);

                                        //IJ.log("photon : "+loc[partNumber][threadNumber].getA()[0][0]);
                                        if (loc[threadNumber].getA()[0][0]>photonThreshold){
                                            
                                                //only displacement < than 3 pixels from the init
                                            if ((Math.abs(loc[threadNumber].getX())<10*dp.param.xystep)&&(Math.abs(loc[threadNumber].getY())<10*dp.param.xystep)){
                                                int px_pix=(int)(maxs[2]-(loc[threadNumber].getX()/dp.param.xystep));
                                                int py_pix=(int)(maxs[3]-(loc[threadNumber].getY()/dp.param.xystep));
                                                
                                                
                                                if ((px_pix>=0)&&(px_pix<mask.length)&&(py_pix>=0)&&(py_pix<mask[0].length)&&(mask[px_pix][py_pix]==false)){
                                                    
                                                    
                                                    if ((loc[threadNumber].getZ()>minZ)&&(loc[threadNumber].getZ()<maxZ)){
                                                        
                                                        numb++;
                                                        n++;
                                                        //IJ.log("zz "+rangePredetection[maxs[1]]+"  "+loc[partNumber][threadNumber].getZ());
                                                        double my_x=1000*(maxs[2]*dp.param.xystep-loc[threadNumber].getX());
                                                        double my_y=1000*(maxs[3]*dp.param.xystep-loc[threadNumber].getY());
                                                        double my_z=(1000*loc[threadNumber].getZ());
                                                        double my_A=loc[threadNumber].getA()[0][0];
                                                        double my_B=loc[threadNumber].getB()[0];                                                        //double my_Score=(loc[partNumber][threadNumber].getGlobalLikelihood());
                                                        //double my_Score=(loc[partNumber][threadNumber].getGlobalLikelihood());
                                                        double my_Score= scoreCompute(patch,model,maxModelThreshold,background);
                                                        double my_crlbx=(1000*loc[threadNumber].getCRLBX());
                                                        double my_crlby=(1000*loc[threadNumber].getCRLBY());
                                                        double my_crlbz=(1000*loc[threadNumber].getCRLBZ());

                                                        
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
                                                                if ((model[i][ii]-background)>maxModelThreshold){
                                                                    mask[maxs[2]+j][maxs[3]+jj]=true;//has been treated

                                                                }
                                                                if (show){
                                                                    if (imageNumber<nbGlobalMask){
                                                                        globalmask[imageNumber][maxs[2]+j][maxs[3]+jj]=model[i][ii];
                                                                        globalmask2[imageNumber][maxs[2]+j][maxs[3]+jj]=patch[i][ii];
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        //if (imageNumber<nbGlobalMask){
                                                            //globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=AandB[0];
//                                                            globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=my_Score;
                                                            //globalmask[imageNumber][maxs[2]+0][maxs[3]+0]=rangePredetection[maxs[1]];
                                                            //globalmask[imageNumber][maxs[2]+1][maxs[3]+0]=my_z;
                                                        //}
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
                                        if ((patch[i][ii]>AandB[2]*.5)){
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
                        
                        if((sliceNumber-totNumberImageProcessed)<nbThread){
                              
                            
                            
                            stopReading=true;
                        }
                        
                        if (imageNumber<nbGlobalMask){
                            computedGlobalMask++;
                            
                            if ( computedGlobalMask == (nbGlobalMask) ){
                                if (show){
                                    ImageShow.imshow(globalmask,"model_frame");
                                    ImageShow.imshow(globalmask2,"image_frame");
                                    globalmask=null;
                                    globalmask2=null;
                                }
                                else{
                                    pr.launch_();
                                }
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
                    try{
                        Thread.sleep(500);
                        
                    }
                    catch(Exception ee){
                        IJ.log("problem lock 2 in class localizationPipelineMany");//do it
                    }
                }
            }
        }
        
        
    }
    
        

        double validLocalization(double [][] patch, double [][] model,double maxModelThreshold){
            
            
            
            double sum=0;
            double count=0;
            for (int i=0;i<patch.length;i++){
                for (int ii=0;ii<patch[i].length;ii++){
                    if (patch[i][ii]>maxModelThreshold){
                        sum+=(patch[i][ii]-model[i][ii])*(patch[i][ii]-model[i][ii])/model[i][ii];
                        count++;
                    }
                }
            }
            sum/=count;
            return sum;
        }
        
        
        public double scoreCompute(double [][] patch, double [][] model,double maxModelThreshold, double background){
            
            
            
            double sum=0;
            double sumModel=0;
            for (int i=0;i<patch.length;i++){
                for (int ii=0;ii<patch[i].length;ii++){
                //if (model[i]-background>maxModelThreshold){
                    sum+=(patch[i][ii]-model[i][ii])*(patch[i][ii]-model[i][ii])/model[i][ii];
                    sumModel++;
                //}
                }
            }
            sum/=sumModel;
            return sum;
        }

        

        private double [] getAandB(double [][] patch, double [][] psf, double thresholdPercentage){
            double maxPatch=Double.NEGATIVE_INFINITY;
            double maxi=Double.NEGATIVE_INFINITY;
            for (int i=0;i<patch.length;i++){
                for (int ii=0;ii<patch[i].length;ii++){
                    if (psf[i][ii]>maxi){
                        maxi=psf[i][ii];
                    }
                    if (patch[i][ii]>maxPatch){
                        maxPatch=patch[i][ii];
                    }
                }
            }
            //IJ.log("maxi "+maxi);
            double threshold=maxi*thresholdPercentage;
            double countB=0;
            double B=0;
            double A=0;
            double countA=0;
            for (int i=0;i<patch.length;i++){
                for (int ii=0;ii<patch[i].length;ii++){
                    if (psf[i][ii]<threshold){
                        B+=patch[i][ii];
                        countB++;
                    }
                    else{
                        A+=patch[i][ii];
                        countA++;
                    }
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

        
        
        

        private double [] getAandB(double [][] patch, double [][] psf, double [][] scmos, double thresholdPercentage){
            double maxPatch=Double.NEGATIVE_INFINITY;
            double maxi=Double.NEGATIVE_INFINITY;
            for (int i=0;i<patch.length;i++){
                for (int ii=0;ii<patch[i].length;ii++){
                    if (psf[i][ii]>maxi){
                        maxi=psf[i][ii];
                    }
                    if (patch[i][ii]-scmos[i][ii]>maxPatch){
                        maxPatch=patch[i][ii]-scmos[i][ii];
                    }
                }
            }
            //IJ.log("maxi "+maxi);
            double threshold=maxi*thresholdPercentage;
            double countB=0;
            double B=0;
            double A=0;
            double countA=0;
            for (int i=0;i<patch.length;i++){
                for (int ii=0;ii<patch[i].length;ii++){
                    if (psf[i][ii]<threshold){
                        B+=patch[i][ii]-scmos[i][ii];
                        countB++;
                    }
                    else{
                        A+=patch[i][ii]-scmos[i][ii];
                        countA++;
                    }
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
