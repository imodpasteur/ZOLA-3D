/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.cuda.*;

import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.cublasHandle;
import jcuda.jcusparse.cusparseHandle;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import java.util.concurrent.locks.ReentrantLock;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_N;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.ArrayList;
import jcuda.runtime.cudaError;
/**
 *
 * @author benoit
 */
public class Modelmany_double_ {
    
    
    
    
            
    
    Pointer device_subwindow;
    Pointer host_subwindow;
    double [] subwindow;
    
    
    
    Pointer device_subwindowSCMOS;
    Pointer host_subwindowSCMOS;
    double [] subwindowSCMOS;
    
    
    
    Pointer device_model;
    Pointer device_tmp;
    Pointer device_ones;
    Pointer device_AandB;
    
    Pointer device_psf;
    
    Pointer device_tmpLong;
    
    Pointer device_alpha;
    Pointer device_beta;
    Pointer device_likelihoodResult;
    
    
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;
    
    CUstream custream;
    
    PhaseParameters param;
    //PSFphaseJCudaFastDoubleMany psfMany;
    PSFmany_float_ psfMany_f;
    
    
    
    int nbParam=5;
    
    
    double [] X;
    double [] Y;
    double [] Z;
    double [] Zoil;
    double [] AandB;
    
    double [] likelihoodResult;
    Pointer host_likelihoodResult; 
    Pointer host_AandB;
    
    int numberPSFperModel;
    int numberModel;
    public int numberModelToCompute;
    int numberModelToCompute2;
    
    int sizeImage;
    
    //boolean parameterUpdated=true;
    boolean newImageSet=false;
    
    private final ReentrantLock lock = new ReentrantLock();
    //private final Condition cond = lock.newCondition();
    
    
    
    
    int numberPSFLaunched=0;//it is important to initialize this value before to start localize for ALL threads 
    int allPSFcomputed=0;
    boolean hasBeenComputed=false;
    
    
    ArrayList<Integer> freePosit = new ArrayList<Integer>();
    
    private  Object monitor = new Object();
    private  Object monitor2 = new Object();
    //private final  AtomicInteger atomic_numberPSFLaunched =  new AtomicInteger(0) ;
    //private  AtomicInteger atomic_allPSFcomputed =  new AtomicInteger(0) ;
    
    PSFmany_float_ psfMany_f_crlb;
    Pointer device_AandB_crlb;
    Pointer host_AandB_crlb;
    Pointer device_model_crlb;
    Pointer device_tmp_crlb;
    Pointer device_psf_crlb;
    double [] X_crlb;
    double [] Y_crlb;
    double [] Z_crlb;
    double [] Zoil_crlb;
    double [] AandB_crlb;
    Pointer device_crlbResult;
    Pointer device_ones_crlb;
    
    
            
    /*PSFphaseJCudaFastDoubleModelMany(PSFphaseJCudaFastDoubleMany psfMany,PhaseRetrievalParametersDouble param,int numberPSF){
        
        
        numberPSFLaunched=0;
        this.psfMany=psfMany;
        this.numberPSF=numberPSF;
        this.numberPSFToCompute=numberPSF;
        this.param=param;
        sizeImage=param.sizeoutput;
        
        for (int i=0;i<numberPSF;i++){
            freePosit.add(i);
        }
        
        construct();
        
    }*/
    
    //PrintTMP pt = new PrintTMP(2000);
    
    boolean isSCMOS=false;
    
    //single emitter fitting constructor
    Modelmany_double_(PSFmany_float_ psfMany_f,PSFmany_float_ psfMany_f_crlb,PhaseParameters param,int numberModel,boolean isSCMOS){
        
        modelmany_double_(psfMany_f,psfMany_f_crlb,param,numberModel,1,isSCMOS);
        
      
    }
    
    //multi emitter fitting constructor
    Modelmany_double_(PSFmany_float_ psfMany_f,PSFmany_float_ psfMany_f_crlb,PhaseParameters param,int numberModel,int numberPSFperModel,boolean isSCMOS){
        
        modelmany_double_(psfMany_f,psfMany_f_crlb,param,numberModel,numberPSFperModel,isSCMOS);
        
      
    }
    
    
    void modelmany_double_(PSFmany_float_ psfMany_f,PSFmany_float_ psfMany_f_crlb,PhaseParameters param,int numberModel,int numberPSFperModel,boolean isSCMOS){
        
        
        this.isSCMOS=isSCMOS;
        
        //atomic_numberPSFLaunched.set(0);
        //atomic_allPSFcomputed.set(0);
        numberPSFLaunched=0;
        this.psfMany_f=psfMany_f;
        this.psfMany_f_crlb=psfMany_f_crlb;
        this.numberModel=numberModel;
        this.numberPSFperModel=numberPSFperModel;
        this.numberModelToCompute=numberModel;
        this.param=param;
        sizeImage=param.sizeoutput;
        
        for (int i=0;i<numberModel;i++){
            freePosit.add(i);
        }
        
        
        
        construct();
        
        
        //pt.start();
        
    }
    
    public void resetNumberPSFToCompute(){
        
        synchronized(monitor) {
            for (int i=numberModelToCompute;i<numberModel;i++){
                freePosit.add(i);
            }
            numberModelToCompute=numberModel;
            monitor.notify();
        }
        
    }
    
    
    //at the end of the localization process, the number of molecule to localize can be < to numberPSF
    //so the number of localization done in parallel should decrease according the the available frame number for example
    public void decrementNumberPSFToCompute(){
        
        synchronized(monitor) {
            numberModelToCompute--;
            freePosit.remove(0);
            monitor.notify();
        }
        
    }
    
    public void decrementNumberPSFToCompute(int number){
        
        synchronized(monitor) {
            for (int i=0;i<number;i++){
                numberModelToCompute--;
                freePosit.remove(0);
            }
            monitor.notify();
        }
        
    }
    
    
    
    private void construct(){
        
        
        
        likelihoodResult = new double [numberModel];
        X = new double [numberModel*numberPSFperModel];
        Y = new double [numberModel*numberPSFperModel];
        Z = new double [numberModel*numberPSFperModel];
        Zoil = new double [numberModel*numberPSFperModel];
        AandB = new double [numberModel*numberPSFperModel*2];
        host_likelihoodResult = Pointer.to(likelihoodResult);
        
        
        
        
        subwindow=new double [sizeImage*sizeImage*numberModel];
        host_subwindow=Pointer.to(subwindow);
        device_subwindow = new Pointer();
        cudaMalloc(device_subwindow, sizeImage*sizeImage*numberModel * Sizeof.DOUBLE);
        
        
        if (isSCMOS){
            subwindowSCMOS=new double [sizeImage*sizeImage*numberModel];
            host_subwindowSCMOS=Pointer.to(subwindowSCMOS);
            device_subwindowSCMOS = new Pointer();
            cudaMalloc(device_subwindowSCMOS, sizeImage*sizeImage*numberModel * Sizeof.DOUBLE);
        }
        
        
        device_model = new Pointer();
        cudaMalloc(device_model, sizeImage*sizeImage*numberModel* Sizeof.DOUBLE);
        
        device_tmp = new Pointer();
        cudaMalloc(device_tmp, sizeImage*sizeImage*numberModel* Sizeof.DOUBLE);
        
        
        device_tmpLong = new Pointer();
        cudaMalloc(device_tmpLong, sizeImage*sizeImage*numberModel*numberPSFperModel* Sizeof.DOUBLE);
        
        
        
        
        
        device_AandB = new Pointer();
        cudaMalloc(device_AandB, 2*numberModel*numberPSFperModel* Sizeof.DOUBLE);
        
        
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.DOUBLE);
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.DOUBLE);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        
        
        
        double [] ones=new double[sizeImage*sizeImage*numberModel];
        for (int i=0;i<sizeImage*sizeImage*numberModel;i++){
            ones[i]=1;
        }
        Pointer h_ones=Pointer.to(ones);
        device_ones = new Pointer();
        JCuda.cudaMalloc(device_ones, sizeImage*sizeImage*numberModel *Sizeof.DOUBLE);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_ones, h_ones, sizeImage*sizeImage*numberModel*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        device_likelihoodResult=new Pointer();
        JCuda.cudaMalloc(device_likelihoodResult, numberModel* Sizeof.DOUBLE);
        
        handlecublas=MyCudaStream.getHandleCublas(param.stream);
        handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
        custream=MyCudaStream.getCUstream(param.stream);
        host_AandB=Pointer.to(AandB);
        
        
        
        
        
        //CRLB data//$$$$$$$$$$$$$$$$$$$$$$$$$
        
        device_AandB_crlb = new Pointer();
        cudaMalloc(device_AandB_crlb, 2*psfMany_f_crlb.numberPSF* Sizeof.DOUBLE);
        device_model_crlb = new Pointer();
        cudaMalloc(device_model_crlb, sizeImage*sizeImage*psfMany_f_crlb.numberPSF* Sizeof.DOUBLE);
        device_tmp_crlb = new Pointer();
        cudaMalloc(device_tmp_crlb, sizeImage*sizeImage*nbParam*nbParam* Sizeof.DOUBLE);
        
        X_crlb = new double [psfMany_f_crlb.numberPSF];
        Y_crlb = new double [psfMany_f_crlb.numberPSF];
        Zoil_crlb = new double [psfMany_f_crlb.numberPSF];
        Z_crlb = new double [psfMany_f_crlb.numberPSF];
        AandB_crlb = new double [psfMany_f_crlb.numberPSF*2];
        host_AandB_crlb=Pointer.to(AandB_crlb);
        double [] ones_crlb=new double[sizeImage*sizeImage*psfMany_f_crlb.numberPSF];
        for (int i=0;i<sizeImage*sizeImage*psfMany_f_crlb.numberPSF;i++){
            ones_crlb[i]=1;
        }
        Pointer h_ones_crlb=Pointer.to(ones_crlb);
        
        
        device_ones_crlb = new Pointer();
        JCuda.cudaMalloc(device_ones_crlb, sizeImage*sizeImage*psfMany_f_crlb.numberPSF *Sizeof.DOUBLE);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_ones_crlb, h_ones_crlb, sizeImage*sizeImage*psfMany_f_crlb.numberPSF*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        device_crlbResult=new Pointer();
        JCuda.cudaMalloc(device_crlbResult, nbParam*nbParam* Sizeof.DOUBLE);
        
    }
    
    
    public void setX(int id,double value){
        this.X[id]=value;
    }
    
    public void setY(int id,double value){
        this.Y[id]=value;
    }
    
    public void setZ(int id,double value){
        this.Z[id]=value;
    }
    public void setZoil(int id,double value){
        this.Zoil[id]=value;
    }
    
    public void setA(int id,double value){
        this.AandB[id*2]=value;
        //parameterUpdated=false;
    }
    
    public void setB(int id,double value){
        this.AandB[id*2+1]=value;
        //parameterUpdated=false;
    }
    
    
    void updateAandB(){
        int cudaResult = cudaMemcpyAsync(device_AandB, host_AandB, 2*numberModel*numberPSFperModel*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem cpy A and B");}
        //parameterUpdated=true;
    }
    
    
    
    
    public int setSubWindow(double [] subwindow){
        
        if (sizeImage*sizeImage!=subwindow.length){
            IJ.log("error in image dimension (localizationMany class)   "+subwindow.length+"   "+(sizeImage*sizeImage));
            return -1;
        }
        int id=-1;
        
        synchronized(monitor) {
            while (freePosit.size()==0){//no available position
                //no wait if (frame number == thread number) // this is only to manage the case: (frame number < thread number)
                try{Thread.sleep(0, 5000);}catch(Exception ee){}
                IJ.log("WARNING: no available position setSubWinBlocked "+(param.stream-1)+"   freePosit:"+freePosit.size());
            }
            
            newImageSet=true;
            id=freePosit.get(0);
            freePosit.remove(0);
            for (int i=0;i<subwindow.length;i++){
                this.subwindow[id*subwindow.length+i]=subwindow[i];
            }
        }
        
        
        return id;
    }
    
    
    
    
    
    
    
    public int setSubWindowScmos(double [] subwindow,double [] subwindowSCMOS){
        
        if (sizeImage*sizeImage!=subwindow.length){
            IJ.log("error in image dimension (localizationMany class)   "+subwindow.length+"   "+(sizeImage*sizeImage));
            return -1;
        }
        int id=-1;
        
        synchronized(monitor) {
            while (freePosit.size()==0){//no available position
                //no wait if (frame number == thread number) // this is only to manage the case: (frame number < thread number)
                try{Thread.sleep(0, 500000);}catch(Exception ee){}
                IJ.log("setSubWinBlocked "+(param.stream-1)+"   freePosit:"+freePosit.size());
            }
            
            newImageSet=true;
            id=freePosit.get(0);
            freePosit.remove(0);
            for (int i=0;i<subwindow.length;i++){
                this.subwindow[id*subwindow.length+i]=subwindow[i];
                if (isSCMOS){
                    this.subwindowSCMOS[id*subwindowSCMOS.length+i]=subwindowSCMOS[i];
                }
                else{
                    IJ.log("ERROR... you should not call setSubWindowScmos function if isSCMOS==false");
                }
            }
        }
        
        
        return id;
    }
    
    
    
    
    
    
    
    
    
    public void freePosit(int id){
        
        
        synchronized(monitor) {
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));
            freePosit.add(id);
        }
        /*try{
            lock.lock();
            
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));
            freePosit.add(id);
        }
        catch(Exception urf){}
        finally{lock.unlock();}*/
    }
    
    void computeLikelihood(){
        int   cudaResult;
        if (newImageSet){
            
            cudaResult=cudaMemcpyAsync(device_subwindow, host_subwindow, sizeImage*sizeImage*numberModel*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem cpy subwindow");}
            
            if (isSCMOS){
                cudaResult=cudaMemcpyAsync(device_subwindowSCMOS, host_subwindowSCMOS, sizeImage*sizeImage*numberModel*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem cpy subwindow");}
            
            }
            
            newImageSet=false;
        }
        
        //if (!parameterUpdated){
            updateAandB();
        //}
        
        
        psfMany_f.computePSF(X, Y, Zoil,Z);
        
        
        this.device_psf=psfMany_f.getPointerPSF();
        
        
        
        
        //compute model
        if (numberPSFperModel==1){
            if (isSCMOS){
                MyVecDouble.addPhotonsAndBackgroundMany_scmos(custream, sizeImage*sizeImage*numberModel, sizeImage*sizeImage, device_model, device_psf, device_AandB,device_subwindowSCMOS);
            }
            else{
                MyVecDouble.addPhotonsAndBackgroundMany(custream, sizeImage*sizeImage*numberModel, sizeImage*sizeImage, device_model, device_psf, device_AandB);
            }
        }
        else{
            if (isSCMOS){
                MyVecDouble.addPhotonsAndBackgroundManyReshuffle_scmos(custream, sizeImage*sizeImage*numberModel*numberPSFperModel, sizeImage*sizeImage,numberPSFperModel, device_tmpLong, device_psf, device_AandB,device_subwindowSCMOS);
            }
            else{
                MyVecDouble.addPhotonsAndBackgroundManyReshuffle(custream, sizeImage*sizeImage*numberModel*numberPSFperModel, sizeImage*sizeImage,numberPSFperModel, device_tmpLong, device_psf, device_AandB);
            }
            
            
            //psfMany_f.imshow(sizeImage*sizeImage*numberModel*numberPSFperModel, sizeImage, device_psf, "psf", "DOUBLE");
            //psfMany_f.imshow(sizeImage*sizeImage*numberModel*numberPSFperModel, sizeImage, device_tmpLong, "psfLong", "DOUBLE");
            
            jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,sizeImage*sizeImage*numberModel,numberPSFperModel,device_alpha,device_tmpLong,sizeImage*sizeImage*numberModel,device_ones,1,device_beta,device_model,1);
        
            //psfMany_f.imshow(sizeImage*sizeImage*numberModel, sizeImage, device_model, "model", "DOUBLE");
            //IJ.log("here convert device_tmpLong to device_model");
        }
        
        
        
        //compute Likelihood
        MyVecDouble.computePoissonLikelihood(custream,sizeImage*sizeImage*numberModel, device_tmp, device_subwindow, device_model);

        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizeImage*sizeImage,numberModel,device_alpha,device_tmp,sizeImage*sizeImage,device_ones,1,device_beta,this.device_likelihoodResult,1);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR Dgemv likelihood");}
        
        
        cudaResult=cudaMemcpyAsync(host_likelihoodResult, device_likelihoodResult, numberModel*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        
        
        
        
        
        
    }
    
    
    
    
    
    public void getModel(int id,double [] model){
        
        int cudaResult=cudaMemcpyAsync(Pointer.to(model), device_model.withByteOffset(sizeImage*sizeImage*id*Sizeof.DOUBLE), sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem cpy getModel");}
        
    }
    
    
    
    
    public void computeFisherMatrix(int idimage,double x, double y, double zoil, double z, double a, double b,double h,double [][] fisherMatrix){
        
        //try{
        //    lock.lock();
            
        synchronized(monitor) {
            
            
            int id=0;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=1;
            this.X_crlb[id]=x-h;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=2;
            this.X_crlb[id]=x+h;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=3;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y-h;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=4;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y+h;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=5;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z-h;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=6;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z+h;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b;
            
            id=7;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a-h*1000;
            this.AandB_crlb[id*2+1]=b;
            
            id=8;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a+h*1000;
            this.AandB_crlb[id*2+1]=b;
            
            id=9;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b-h*10;
            
            id=10;
            this.X_crlb[id]=x;
            this.Y_crlb[id]=y;
            this.Z_crlb[id]=z;
            this.Zoil_crlb[id]=zoil;
            this.AandB_crlb[id*2]=a;
            this.AandB_crlb[id*2+1]=b+h*10;
            
            
            
            int cudaResult = cudaMemcpyAsync(device_AandB_crlb, host_AandB_crlb, 2*psfMany_f_crlb.numberPSF*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
            if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem cpy A and B");}
            
            
            
            
            
            
            psfMany_f_crlb.computePSF(X_crlb, Y_crlb, Zoil_crlb, Z_crlb);
            this.device_psf_crlb=psfMany_f_crlb.getPointerPSF();
            
            
            
            
            MyVecDouble.addPhotonsAndBackgroundMany(custream, sizeImage*sizeImage*psfMany_f_crlb.numberPSF, sizeImage*sizeImage, device_model_crlb, device_psf_crlb, device_AandB_crlb);
            //psfMany_f_crlb.imshow(sizeImage*sizeImage*11, sizeImage, device_model_crlb, "model", "DOUBLE");
             
            double [] resultCRLBCuda=new double [nbParam*nbParam];
            
            MyVecDouble.computeCRLB(custream,sizeImage*sizeImage*nbParam*nbParam,nbParam,device_tmp_crlb,device_model_crlb,h);
            //MyVecDouble.computeCRLBwrong(custream,sizeImage*sizeImage*nbParam*nbParam,nbParam,device_tmp_crlb,device_model_crlb,this.device_subwindow.withByteOffset(sizeImage*sizeImage*idimage*Sizeof.DOUBLE),h);
            
            //psfMany_f_crlb.imshow(sizeImage*sizeImage*nbParam*nbParam, sizeImage, device_tmp_crlb, "tmp", "DOUBLE");
            
            
            cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizeImage*sizeImage,nbParam*nbParam,device_alpha,device_tmp_crlb,sizeImage*sizeImage,device_ones_crlb,1,device_beta,this.device_crlbResult,1);
            if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR Dgemv likelihood");}
            
            Pointer crlbRes=Pointer.to(resultCRLBCuda);
            cudaResult=cudaMemcpyAsync(crlbRes, device_crlbResult, nbParam*nbParam*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
            
            //double [][] mat=new double [nbParam][nbParam];
            for (int u=0;u<nbParam;u++){
                for (int uu=0;uu<nbParam;uu++){
                    fisherMatrix[u][uu]=resultCRLBCuda[u+uu*nbParam];
                    //mat[u][uu]=resultCRLBCuda[u+uu*nbParam];
                }
            }
            
            resultCRLBCuda=null;
            JCuda.cudaFree(crlbRes);
            
            
            
            
            
        }
        //catch(Exception urf){}
        //finally{lock.unlock();}
        
    }
    
    
    
    
    
    
    public void setParameters(int id, double x, double y, double zoil, double z, double a, double b){
        
        
        //IJ.log("get lik "+id+" / "+numberPSFToCompute);
        //IJ.log("get lik in PSFphase... "+a+"  "+b+"  "+x+"  "+y+"  "+z);
        this.setA(id, a);
        this.setB(id, b);
        this.setX(id, x);
        this.setY(id, y);
        this.setZoil(id, zoil);
        this.setZ(id, z);
        
        
        
        
    }
    
    
    
    //set parameter for each model (id) (the dim of vectors should be: numberPSFperModel
    public void setParameters(int id, double [] x, double [] y, double zoil, double [] z, double [] a, double [] b){
        
        
        //IJ.log("get lik "+id+" / "+numberPSFToCompute);
        //IJ.log("get lik in PSFphase... "+a+"  "+b+"  "+x+"  "+y+"  "+z);
        for (int i=0;i<x.length;i++){
            this.setA(id*this.numberPSFperModel+i, a[i]);
            this.setB(id*this.numberPSFperModel+i, b[i]);
            this.setX(id*this.numberPSFperModel+i, x[i]);
            this.setY(id*this.numberPSFperModel+i, y[i]);
            this.setZoil(id*this.numberPSFperModel+i, zoil);
            this.setZ(id*this.numberPSFperModel+i, z[i]);
            
            
        }
        
        
        
        
        
        
        
    }
    
    
    
    public double getLikelihood(int id, double x, double y, double zoil, double z, double a, double b){
        
        //IJ.log("get lik "+id+" / "+numberPSFToCompute);
        //IJ.log("get lik in PSFphase... "+a+"  "+b+"  "+x+"  "+y+"  "+z);
        this.setA(id, a);
        this.setB(id, b);
        this.setX(id, x);
        this.setY(id, y);
        this.setZoil(id, zoil);
        this.setZ(id, z);
        
        
        
        return getLikelihood(id);
        
    }
    
    
    
    
    
    public double getLikelihood(int id){
        
        
        
        synchronized(monitor) {
            hasBeenComputed=false; 
            while (!hasBeenComputed){
                
                numberPSFLaunched++;
                
                int atom=numberPSFLaunched;
                if (((atom)==this.numberModelToCompute)&&(!hasBeenComputed)){//no need to use numberPSFperModel here
                    
                    hasBeenComputed=true;
                    computeLikelihood();
                    numberModelToCompute2=this.numberModelToCompute;
                    allPSFcomputed=0;
                    monitor.notifyAll();
                }
                
                
                
                
                while (atom!=this.numberModelToCompute){
                    
                    try{ 
                        monitor.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    if (!hasBeenComputed){
                        
                        //IJ.log("notify "+(param.stream-1)+"   "+numberPSFLaunched+"/"+numberPSFToCompute);
                        
                    }
                    atom=numberPSFLaunched;
                }

                
                    
                if (!hasBeenComputed){//possible if numberPSFToCompute decreases at the end of the frames (notify called)
                    //IJ.log("free thread");
                    numberPSFLaunched--;
                    //we decrement and continue loop
                }
                  
            }
        }
        
        
        
        
        //wait finish to put back atomic_numberPSFLaunched=0
          
        synchronized(monitor2) {
            allPSFcomputed++;
            //IJ.log("All "+allPSFcomputed);
            if ((allPSFcomputed)==numberModelToCompute2){
                numberPSFLaunched=0;
                //IJ.log("notif ok ");
                monitor2.notifyAll();
            }
            else{
                while (allPSFcomputed!=numberModelToCompute2){
                    
                    try{
                        monitor2.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                }
            }
        }
        
        //IJ.log("get lik ok... "+id);
        
        return likelihoodResult[id];
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public void free(){
        //IJ.log("free PSF_model many called");
        JCuda.cudaFree(device_subwindow);
        JCuda.cudaFree(device_model);
        JCuda.cudaFree(device_tmp);
        JCuda.cudaFree(device_tmpLong);
        JCuda.cudaFree(device_ones);
        JCuda.cudaFree(device_AandB);
        
        JCuda.cudaFree(device_alpha);
        
        JCuda.cudaFree(device_beta);
        JCuda.cudaFree(device_likelihoodResult);
        
        
        JCuda.cudaFree( device_AandB_crlb);
        JCuda.cudaFree( device_model_crlb);
        JCuda.cudaFree( device_tmp_crlb);
        JCuda.cudaFree( device_crlbResult);
        JCuda.cudaFree( device_ones_crlb);
        
        JCuda.cudaFree( host_subwindow);
        JCuda.cudaFree( host_likelihoodResult);
        JCuda.cudaFree( host_AandB);
        JCuda.cudaFree( host_AandB_crlb);
        
        
        if (isSCMOS){
            JCuda.cudaFree(device_subwindowSCMOS);
        }
    
            
        //pt.stopRun();
    
    }
    
    
    
    
    
    
    
    class PrintTMP extends Thread{
        int time_ms;
        boolean launched=true;
        PrintTMP(int time_ms){
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
                    
                    IJ.log("while: "+(param.stream-1)+"   "+numberPSFLaunched+"/"+numberModelToCompute);
                }
                
            }catch(Exception e){}
        }
    }
    
    
}
