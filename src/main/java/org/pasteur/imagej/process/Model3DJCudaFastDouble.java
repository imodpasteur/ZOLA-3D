/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;
import org.pasteur.imagej.cuda.*;

import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas2;
import static  jcuda.jcublas.JCublas2.cublasCreate;
import jcuda.jcublas.cublasHandle;
import jcuda.jcublas.cublasOperation;
import static jcuda.jcusparse.JCusparse.cusparseCreate;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaError;




 
 
import static jcuda.jcublas.JCublas2.cublasCreate;
import static jcuda.jcublas.JCublas2.cublasDestroy;
import static jcuda.jcublas.JCublas2.cublasDgemm;
import static jcuda.jcublas.JCublas2.cublasGetVector;
import static jcuda.jcublas.JCublas2.cublasSetPointerMode;
import static jcuda.jcublas.JCublas2.cublasSetVector;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_N;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.jcublas.cublasPointerMode.CUBLAS_POINTER_MODE_DEVICE;
import static jcuda.runtime.JCuda.cudaFree;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.cudaMemcpyKind.*;
 
import java.util.Random;
 
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.JCublas2;
import jcuda.jcublas.cublasHandle;
import jcuda.jcusparse.cusparseHandle;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import jcuda.runtime.cudaMemcpyKind;

/**
/**
 *
 * @author benoit
 */
public class Model3DJCudaFastDouble {
    
    
    
    boolean isAmplitude=false ;
    boolean isBackground = false;
    boolean isBackground2D = false;
    boolean isSCMOS = false;
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;
    
    CUstream custream;
    
    Pointer devresone;
    double [] resone;
    Pointer hostresone;
    
    Pointer devresstack;
    double [] resstack;
    Pointer hostresstack;
    
    
    Pointer device_Amplitude;
    Pointer device_Background;
    Pointer device_Background_2D;
    Pointer device_SCMOS;
    double [] background2D;
    double [] scmos;
    
    Pointer device_input;
    //Pointer device_logI;
    Pointer device_model;
    Pointer device_tmp;
    Pointer device_psf;
    Pointer device_ones;
    
    Pointer device_alpha;
    Pointer device_beta;
    
    
    Pointer host_input;
    double [] res;
    double [] input;
    int sizepow;
    int size;
    int nbImage;
    PhaseRetrievalParametersDouble param;
    public Model3DJCudaFastDouble(PhaseRetrievalParametersDouble param, double [][][] inputImage, double [][] scmos){
        this.param=param;
        
        if ((inputImage[0].length!=inputImage[0][0].length)){
            IJ.log("ERROR, image should be square (PhaseJCudaFast class)");return;
        }
        if ((inputImage[0].length!=param.sizeoutput)){
            IJ.log("ERROR images size should be the same :   sizOutput:"+param.sizeoutput+"   sizImage:"+inputImage[0].length);
        }
        
        this.size=inputImage[0].length;
        sizepow=inputImage[0].length*inputImage[0][0].length;
        nbImage=inputImage.length;
        input=new double [nbImage*sizepow];
        background2D=new double [sizepow];
        if (scmos!=null){
            isSCMOS=true;
            this.scmos=new double [sizepow];
        }
        else{
            isSCMOS=false;
            this.scmos=null;
        }
        //double [] logI=new double [nbImage*sizepow];
        //double ll=0;
        for (int rz=0;rz<nbImage;rz++){
            for (int i=0;i<size;i++){
                for (int ii=0;ii<size;ii++){
                    input[rz*sizepow+i*size+ii]=inputImage[rz][i][ii];
                    /*logI[rz*sizepow+i*size+ii]=0;
                    loop:for (double u=2;u<inputImage[rz][i][ii];u++){
                        logI[rz*sizepow+i*size+ii]+=Math.log(u);
                        if (inputImage[rz][i][ii]>2000){
                            break loop;//can be very long if higher than 5000 and it does not change the result. This is added to deal with eventually bad double precision
                        }
                    }
                    ll+=logI[rz*sizepow+i*size+ii];*/
                }
            }
        }
        if (isSCMOS){
            for (int i=0;i<size;i++){
                for (int ii=0;ii<size;ii++){
                    this.scmos[i*size+ii]=scmos[i][ii];
                }
            }
        }
        //IJ.log("sum of log : "+ll);
        Pointer hh=Pointer.to(input);
        device_input=new Pointer();
        JCuda.cudaMalloc(device_input, sizepow*nbImage *Sizeof.DOUBLE);
        this.setHost2Device(device_input, hh, sizepow*nbImage, Sizeof.DOUBLE);
        
        
        device_Amplitude=new Pointer();
        JCuda.cudaMalloc(device_Amplitude, nbImage *Sizeof.DOUBLE);
        
        device_Background=new Pointer();
        JCuda.cudaMalloc(device_Background, nbImage *Sizeof.DOUBLE);
        
        
        device_Background_2D=new Pointer();
        JCuda.cudaMalloc(device_Background_2D, sizepow *Sizeof.DOUBLE);
        if (isSCMOS){
            device_SCMOS=new Pointer();
            JCuda.cudaMalloc(device_SCMOS, sizepow *Sizeof.DOUBLE);
            this.setHost2Device(device_SCMOS, Pointer.to(this.scmos), sizepow, Sizeof.DOUBLE);
            int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setbckg "+cudaResult);}
        }
        
//        Pointer hhlogI=Pointer.to(logI);
//        device_logI=new Pointer();
//        JCuda.cudaMalloc(device_logI, sizepow*nbImage *Sizeof.DOUBLE);
//        this.setHost2Device(device_logI, hhlogI, sizepow*nbImage, Sizeof.DOUBLE);
        
        double [] ones=new double[sizepow*nbImage];
        for (int i=0;i<sizepow*nbImage;i++){
            ones[i]=1;
        }
        Pointer h_ones=Pointer.to(ones);
        
        device_ones = new Pointer();
        JCuda.cudaMalloc(device_ones, sizepow*nbImage *Sizeof.DOUBLE);
        setHost2Device(device_ones,h_ones,sizepow*nbImage,Sizeof.DOUBLE);
        
        device_model=new Pointer();
        device_psf=new Pointer();
        device_tmp=new Pointer();
        JCuda.cudaMalloc(device_model, sizepow*nbImage *Sizeof.DOUBLE);
        JCuda.cudaMalloc(device_tmp, sizepow*nbImage *Sizeof.DOUBLE);
        JCuda.cudaMalloc(device_psf, sizepow*nbImage *Sizeof.DOUBLE);
        
        
        res=new double[sizepow*nbImage];
        
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.DOUBLE);
        cudaMemcpyAsync(device_alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.DOUBLE);
        cudaMemcpyAsync(device_beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        
        
        handlecublas=MyCudaStream.getHandleCublas(param.stream);
        handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
        custream=MyCudaStream.getCUstream(param.stream);
        
        devresone=new Pointer();
        JCuda.cudaMalloc(devresone, 1 * Sizeof.DOUBLE);
        resone = new double [1];
        hostresone = Pointer.to(resone);
        
        devresstack=new Pointer();
        JCuda.cudaMalloc(devresstack, nbImage * Sizeof.DOUBLE);
        resstack = new double [nbImage];
        hostresstack = Pointer.to(resstack);
        
    }
    
    
    
    
    public void updateImage(double [][][] newimage){
        double [] input=new double [nbImage*sizepow];
        
        for (int rz=0;rz<nbImage;rz++){
            for (int i=0;i<size;i++){
                for (int ii=0;ii<size;ii++){
                    input[rz*sizepow+i*size+ii]=newimage[rz][i][ii];
                    
                }
            }
        }
        Pointer hh=Pointer.to(input);
        this.setHost2Device(device_input, hh, sizepow*nbImage, Sizeof.DOUBLE);
    }
    
    
    
    
    public void free(){
        
        
        JCuda.cudaFree(devresone);
        JCuda.cudaFree(devresstack);
        JCuda.cudaFree(device_Amplitude);
        JCuda.cudaFree(device_Background);
        JCuda.cudaFree(device_Background_2D);
        if (this.scmos!=null){
            JCuda.cudaFree(device_SCMOS);
        }
        JCuda.cudaFree(device_input);
        //JCuda.cudaFree(device_logI);
        JCuda.cudaFree(device_model);
        JCuda.cudaFree(device_tmp);
        JCuda.cudaFree(device_psf);
        JCuda.cudaFree(device_ones);
        JCuda.cudaFree(device_alpha);
        JCuda.cudaFree(device_beta);
        
    }
    
    
    public void setPSF(Pointer device_psf,int zposit){
        
        
        int cudaResult=JCuda.cudaMemcpy(this.device_psf.withByteOffset(param.sizeoutput*param.sizeoutput * Sizeof.DOUBLE*zposit), device_psf, param.sizeoutput*param.sizeoutput * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR copy PSF cuda device 2 device setpsf model");return ;
        }
    }
    
    
    
    
    
    
    public void setPSFMany(Pointer device_psf){
        
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult);}
        
         cudaResult=JCuda.cudaMemcpy(this.device_psf, device_psf, param.sizeoutput*param.sizeoutput * Sizeof.DOUBLE*this.nbImage, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR copy PSF cuda device 2 device setpsf model");return ;
        }
    }
    
    
    
    
    
    
    
    
    public void computeModel(double A_photons, int s, double B_background){
        MyVecDouble.mulScalar(custream,sizepow, device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_psf.withByteOffset(sizepow*s*Sizeof.DOUBLE),A_photons);//square
        MyVecDouble.addScalar(custream,sizepow, device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE), B_background);
        if (this.isSCMOS){
            MyVecDouble.add(custream,sizepow, device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE),device_SCMOS);
        }
    }
    
    
    public void setAmplitude(double [] A_photons){
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setbckg "+cudaResult);}
        this.setHost2Device(device_Amplitude, Pointer.to(A_photons), nbImage, Sizeof.DOUBLE);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setphoton "+cudaResult);}
        
        isAmplitude=true;
    }
    
    public void computeModel(double B_background){
        if (isAmplitude){
            if (isSCMOS){
                MyVecDouble.computeModelMany1_scmos(custream,sizepow*this.nbImage,sizepow,device_model,device_psf,device_Amplitude,B_background,device_SCMOS);
            }
            else{
                MyVecDouble.computeModelMany1(custream,sizepow*this.nbImage,sizepow,device_model,device_psf,device_Amplitude,B_background);
            }
        }
        else{
            IJ.log("oops... problem model computation... amplitude not set");
        }
//        for (int i=0;i<A_photons.length;i++){
//            MyVecDouble.mulScalar(custream,sizepow, device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE), device_psf.withByteOffset(sizepow*i*Sizeof.DOUBLE),A_photons[i]);//square
//        }
//        MyVecDouble.addScalar(custream,sizepow*this.nbImage, device_model, device_model, B_background);
    }
    
    public void setBackground(double [] B_background){
        this.setHost2Device(device_Background, Pointer.to(B_background), nbImage, Sizeof.DOUBLE);
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setbckg "+cudaResult);}
        isBackground=true;
    }
    
    public void setBackground2D(double [][] B_background_2D){
        
            for (int i=0;i<size;i++){
                for (int ii=0;ii<size;ii++){
                    background2D[i*size+ii]=B_background_2D[i][ii];
                }
            }
        this.setHost2Device(device_Background_2D, Pointer.to(background2D), sizepow, Sizeof.DOUBLE);
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setbckg "+cudaResult);}
        isBackground2D=true;
    }
    public void computeModel(){
        
        if (isBackground&&isAmplitude){
            if (isSCMOS){
                MyVecDouble.computeModelMany2_scmos(custream,sizepow*this.nbImage,sizepow,device_model,device_psf,device_Amplitude,device_Background,device_SCMOS);
            }
            else{
                MyVecDouble.computeModelMany2(custream,sizepow*this.nbImage,sizepow,device_model,device_psf,device_Amplitude,device_Background);
            }
        }
        else{
            IJ.log("oops... problem model computation... amplitude or background not set");
        }
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setbckg "+cudaResult);}
        
//        for (int i=0;i<A_photons.length;i++){
//            MyVecDouble.mulScalar(custream,sizepow, device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE), device_psf.withByteOffset(sizepow*i*Sizeof.DOUBLE),A_photons[i]);//square
//            MyVecDouble.addScalar(custream,sizepow, device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE),B_background[i]);//square
//        }
        
    }
    
    
    public void computeModel_background2D(){
        
        if (isBackground2D&&isAmplitude){
            if (isSCMOS){
                MyVecDouble.computeModelMany3_scmos(custream,sizepow*this.nbImage,sizepow,device_model,device_psf,device_Amplitude,device_Background_2D,this.device_SCMOS);
            }
            else{
                MyVecDouble.computeModelMany3(custream,sizepow*this.nbImage,sizepow,device_model,device_psf,device_Amplitude,device_Background_2D);
            }
        }
        else{
            IJ.log("oops... problem model computation... amplitude or background not set");
        }
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda model3d setbckg "+cudaResult);}
        
//        for (int i=0;i<A_photons.length;i++){
//            MyVecDouble.mulScalar(custream,sizepow, device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE), device_psf.withByteOffset(sizepow*i*Sizeof.DOUBLE),A_photons[i]);//square
//            MyVecDouble.addScalar(custream,sizepow, device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*i*Sizeof.DOUBLE),B_background[i]);//square
//        }
        
    }
    
    public double [][][] getModel(){
        Pointer hostoutput=Pointer.to(res);//output point to res
        
        this.setDevice2Host(hostoutput, device_model,sizepow*nbImage,Sizeof.DOUBLE);
        
        //2D conversion:
        
        double [][][] res3D = new double [this.nbImage][this.size][this.size];
        for (int z=0;z<nbImage;z++){
            for (int i=0;i<size;i++){
                for (int ii=0;ii<size;ii++){
                    res3D[z][i][ii]=res[z*sizepow+i*size+ii];
                }
            }
        }
        //long tttttt=System.currentTimeMillis();
        
        return res3D;
    }
    
    
    
    public double [][][] getTmp(){
        Pointer hostoutput=Pointer.to(res);//output point to res
        
        this.setDevice2Host(hostoutput, device_tmp,sizepow*nbImage,Sizeof.DOUBLE);
        
        //2D conversion:
        
        double [][][] res3D = new double [this.nbImage][this.size][this.size];
        for (int z=0;z<nbImage;z++){
            for (int i=0;i<size;i++){
                for (int ii=0;ii<size;ii++){
                    res3D[z][i][ii]=res[z*sizepow+i*size+ii];
                }
            }
        }
        //long tttttt=System.currentTimeMillis();
        
        return res3D;
    }
    
    
    
    public double getLikelihood(int method){
        if (method==0){
            return getLikelihoodPoisson();
        }
        else{
            return getLikelihoodMeanSquare();
        }
    }
    public double getLikelihood(int method, int s){
        if (method==0){
            return getLikelihoodPoisson(s);
        }
        else{
            return getLikelihoodMeanSquare(s);
        }
    }
    public double [] getLikelihoodStack(int method){
        if (method==0){
            return getLikelihoodStackPoisson();
        }
        else{
            return getLikelihoodStackMeanSquare();
        }
    }
    
    
    public double getLikelihoodPoisson(){
//        MyVecDouble.log(custream,sizepow*nbImage, device_tmp, device_model);
//        MyVecDouble.mul(custream,sizepow*nbImage, device_tmp, device_tmp, device_input);
//        //MyVecDouble.sub(custream,sizepow*nbImage, device_tmp, device_tmp, device_logI);
//        MyVecDouble.sub(custream,sizepow*nbImage, device_tmp, device_model, device_tmp);
        
        MyVecDouble.computePoissonLikelihood(custream, sizepow*nbImage, device_tmp, device_input, device_model);
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizepow*nbImage,1,device_alpha,device_tmp,sizepow*nbImage,device_ones,1,device_beta,devresone,1);
        
        
        cudaMemcpyAsync(hostresone, devresone, 1*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        return resone[0];
    }
    
    
    
    
    
    public double getLikelihoodPoisson(int s){ 
//        MyVecDouble.log(custream,sizepow, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE));
//        MyVecDouble.mul(custream,sizepow, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_input.withByteOffset(sizepow*s*Sizeof.DOUBLE));
//        //MyVecDouble.sub(custream,sizepow, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_logI.withByteOffset(sizepow*s*Sizeof.DOUBLE));
//        MyVecDouble.sub(custream,sizepow, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE));
        
        MyVecDouble.computePoissonLikelihood(custream, sizepow*nbImage, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_input.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE));
        
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizepow,1,device_alpha,device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE),sizepow,device_ones.withByteOffset(sizepow*s*Sizeof.DOUBLE),1,device_beta,devresone,1);
        
        
        cudaMemcpyAsync(hostresone, devresone, 1*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        return resone[0];
        
        
    }
    
    
    public double [] getLikelihoodStackPoisson(){ 

        MyVecDouble.computePoissonLikelihood(custream, sizepow*nbImage, device_tmp, device_input, device_model);
        
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizepow,nbImage,device_alpha,device_tmp,sizepow,device_ones,1,device_beta,devresstack,1);
        
        
        cudaMemcpyAsync(hostresstack, devresstack, nbImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        double [] res = new double [resstack.length];
        for (int i=0;i<resstack.length;i++){
            res[i]=resstack[i];
        }
        return res;
        
        
    }
    
    
    
    public double getLikelihoodMeanSquare(){
        
        MyVecDouble.sub(custream,sizepow*nbImage, device_tmp, device_model,device_input);
        MyVecDouble.mul(custream,sizepow*nbImage, device_tmp, device_tmp, device_tmp);
        
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizepow*nbImage,1,device_alpha,device_tmp,sizepow*nbImage,device_ones,1,device_beta,devresone,1);
        
        
        cudaMemcpyAsync(hostresone, devresone, 1*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        return resone[0];
    }
    
    public double [] getLikelihoodStackMeanSquare(){ 

        MyVecDouble.sub(custream,sizepow*nbImage, device_tmp, device_model,device_input);
        MyVecDouble.mul(custream,sizepow*nbImage, device_tmp, device_tmp, device_tmp);
        
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizepow,nbImage,device_alpha,device_tmp,sizepow,device_ones,1,device_beta,devresstack,1);
        
        
        cudaMemcpyAsync(hostresstack, devresstack, nbImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        double [] res = new double [resstack.length];
        for (int i=0;i<resstack.length;i++){
            res[i]=resstack[i];
        }
        return res;
        
        
    }
    
    
    
    public double getLikelihoodMeanSquare(int s){
        MyVecDouble.sub(custream,sizepow, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_model.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_input.withByteOffset(sizepow*s*Sizeof.DOUBLE));
        MyVecDouble.mul(custream,sizepow, device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE), device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE));
        
        
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizepow,1,device_alpha,device_tmp.withByteOffset(sizepow*s*Sizeof.DOUBLE),sizepow,device_ones.withByteOffset(sizepow*s*Sizeof.DOUBLE),1,device_beta,devresone,1);
        
        
        cudaMemcpyAsync(hostresone, devresone, 1*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        
        return resone[0];
        
    }
    
    
    
    
    void setHost2Device(Pointer device,Pointer host,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpyAsync(device, host, size * sizeElement, cudaMemcpyKind.cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR host2device cuda");return ;
        }
    }
    
    
    void setDevice2Host(Pointer host,Pointer device,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpyAsync(host, device, size * sizeElement, cudaMemcpyKind.cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR device2host cuda-");return ;
        }
    }
    
    
    
}
