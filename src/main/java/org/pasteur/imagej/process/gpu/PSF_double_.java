/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.utils.ImageShow;
import jcuda.*;
import jcuda.runtime.*;
import jcuda.jcufft.*;
import jcuda.jcufft.JCufft;
import jcuda.jcusparse.*;
import jcuda.jcublas.*;
import static jcuda.jcusparse.JCusparse.*;
import static jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import static jcuda.jcusparse.cusparseMatrixType.CUSPARSE_MATRIX_TYPE_GENERAL;
import static jcuda.jcusparse.cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE;
import static jcuda.runtime.JCuda.*;
import static jcuda.runtime.cudaMemcpyKind.*;
import ij.ImageStack;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import jcuda.driver.CUstream;
import static jcuda.jcublas.JCublas2.cublasCreate;
import static jcuda.jcublas.JCublas2.cublasSetPointerMode;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.jcublas.cublasPointerMode.CUBLAS_POINTER_MODE_DEVICE;
/**
 *
 * @author benoit
 */
public class PSF_double_ {
    
    
    cufftHandle plan;
    cufftHandle plan2;
    
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;
    CUstream custream;
            
            
    Pointer  device_gaussian;
    Pointer  device_gaussian2bleSize;
    
    Pointer device_phase;
    Pointer device_pupil;
    public Pointer device_kx;
    Pointer device_ky;
    Pointer device_kz;
    Pointer device_kz_oil;
    Pointer device_kz_is_imaginary;
    Pointer device_kz_oil_is_imaginary;
    
    
    Pointer device_psf;
    Pointer device_psfFFT;
    //Pointer device_psfsingle;
    
    
    
    
    
    Pointer device_phtmpFull;
    Pointer device_realtmpFull;
    Pointer device_tmpFull;
    Pointer device_imagtmpFull;
    Pointer device_phtmp;
    Pointer device_realtmp;
    Pointer device_tmp;
    Pointer device_imagtmp;
    Pointer device_fftdata;
    Pointer device_sparseIndexOdd;
    Pointer device_sparseIndexEven;
    Pointer device_sparseIndexShift2D;
    
    Pointer device_sparseIndexOddDisk;
    Pointer device_sparseIndexEvenDisk;
    Pointer device_sparseIndexOddShift2D;
    Pointer device_sparseIndexEvenShift2D;
    Pointer device_sparseIndexOutput;
    
    
    Pointer device_sparseIndexOddShift2DOutput;
    Pointer device_sparseIndexEvenShift2DOutput;
    Pointer device_sparseIndexEvenShift2DOutputNext;
    Pointer device_sparseIndexOddShift2DOutputNext;
    
    Pointer host_phase;
    Pointer host_pupil;
    Pointer host_kx;
    Pointer host_ky;
    Pointer host_kz;
    Pointer host_kz_oil;
    Pointer host_kz_is_imaginary;
    Pointer host_kz_oil_is_imaginary;
    Pointer host_sparseIndexOdd;
    Pointer host_sparseIndexEven;
    Pointer host_sparseIndexShift2D;
    Pointer host_sparseIndexOutput;
    
    Pointer host_sparseIndexOddDisk;
    Pointer host_sparseIndexEvenDisk;
    
    Pointer host_sparseIndexOddShift2D;
    Pointer host_sparseIndexEvenShift2D;
    
    Pointer host_sparseIndexOddShift2DOutput;
    Pointer host_sparseIndexEvenShift2DOutput;
    Pointer host_sparseIndexEvenShift2DOutputNext;
    Pointer host_sparseIndexOddShift2DOutputNext;
    
    int [] sparseIndexOdd;
    int [] sparseIndexEven;
    int [] sparseIndexShift2D;
    int [] sparseIndexEvenDisk;
    int [] sparseIndexOddDisk;
    int [] sparseIndexEvenShift2D;
    int [] sparseIndexOddShift2D;
    int [] sparseIndexOutput;
    int [] sparseIndexOddShift2DOutput;
    int [] sparseIndexEvenShift2DOutput;
    int [] sparseIndexEvenShift2DOutputNext;
    int [] sparseIndexOddShift2DOutputNext;
    double [] phase ;
    double [] pupil ;
    double [] kx ;
    double [] ky ;
    double [] kz ;
    double [] kz_oil ;
    double [] kz_is_imaginary;
    double [] kz_oil_is_imaginary ;
    boolean success=true;
    int sizepow;
    
    double [] matrix;
    
    double [] somResOne;
    Pointer hostSomOne;
    
    double [][] res2D ;
    double [] res1D ;
    double [] res;
    GaussianKernel_ gk;  
    PhaseParameters param;
    
    int cudaResult;
    int sizeoutput1;
    int sizeoutput2;
    public PSF_double_(PhaseParameters param){
        sizeoutput1=param.sizeoutput;
        sizeoutput2=param.sizeoutput*2;
        if (sizeoutput2>param.size){
            sizeoutput2=param.size;
        }
        gk = new GaussianKernel_(sizeoutput2,param.sizeoutput,param.sigmaGaussianKernel,param.stream);
        
        this.param=param;
    
        
        this.sizepow=param.size*param.size;
        
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 1 "+cudaResult+"   "+param.stream);}
        
        
        
        phase = new double[param.sizeDisk];
        pupil = new double[param.sizeDisk];
        kx = new double[param.sizeDisk];
        ky = new double[param.sizeDisk];
        kz = new double[param.sizeDisk];
        kz_oil = new double[param.sizeDisk];
        kz_is_imaginary = new double[param.sizeDisk];
        kz_oil_is_imaginary = new double[param.sizeDisk];
        
        
        sparseIndexOddDisk = new int[param.sizeDisk];
        sparseIndexEvenDisk = new int[param.sizeDisk];
        
        sparseIndexOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutput = new int[sizeoutput2*sizeoutput2];
        
        for (int i=0,j=((param.size/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            for (int ii=0,jj=((param.size/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                sparseIndexOutput[i*sizeoutput2+ii]=j*param.size+jj;
            }
        }
        
        double left=(param.nwat/(param.wavelength))*(param.nwat/(param.wavelength));
        double left_oil=(param.noil/param.wavelength)*(param.noil/param.wavelength);
        double center=param.centerFourierImage;
        int id=0;
        for (int i=0;i<param.size;i++){
            for (int ii=0;ii<param.size;ii++){
                double disk=((((i-center)/(param.size*param.xystep))*((i-center)/(param.size*param.xystep))+((ii-center)/(param.size*param.xystep))*((ii-center)/(param.size*param.xystep))));
                
                if (disk<=param.ringsize){
                    
                    
                    //double iPositShifted=((i+size/2)%size);
                    //double iiPositShifted=((ii+size/2)%size);
                    
                    kx[id]=(i-center)*(2*Math.PI)/(param.xystep*param.size);
                    ky[id]=(ii-center)*(2*Math.PI)/(param.xystep*param.size);
                    
                    
                    double kxx=((i-center)/(param.xystep*(param.size)))*((i-center)/(param.xystep*(param.size)));
                    double kyy=((ii-center)/(param.xystep*(param.size)))*((ii-center)/(param.xystep*(param.size)));


                    //double left=4*Math.PI*Math.PI*k*k;

                    double right=(kxx+kyy);


                    if (left>right){
                        kz[id]=2*Math.PI*Math.sqrt(left-right);
                        kz_is_imaginary[id]=1;
                    }
                    else{
                        kz[id]=0;//-2*Math.PI*Math.sqrt(Math.abs(left-right));
                        kz_is_imaginary[id]=0;
                    }
                    
                    if (left_oil>right){
                        kz_oil[id]=2*Math.PI*Math.sqrt(left_oil-right);
                        kz_oil_is_imaginary[id]=1;
                    }
                    else{
                        kz_oil[id]=0;//2*Math.PI*Math.sqrt(Math.abs(left_oil-right));
                        kz_oil_is_imaginary[id]=0;
                    }
                    
                    if (param.withApoFactor){
                        if (left_oil>right){
                            pupil[id]=(double)(1/Math.pow(1-(right/left_oil),.25));//with apodization factor
                        }
                        else{
                            pupil[id]=0;
                        }
                    }
                    else{
                        pupil[id]=1/(double)Math.sqrt(param.sizeDisk);//like that -> final sum=1 (unuseful actually because we normalize at the end)
                    }
                    
                    
                    phase[id]=0;
                    
                    
                    
                    sparseIndexOddDisk[id]=((((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size))*2)+1;
                    sparseIndexEvenDisk[id]=((((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size))*2);
                    id++;
                }
                
                
                
            }
        }
        
        
        
        
        
        
        
        
        
        
        
        sparseIndexOdd = new int[sizepow];
        sparseIndexEven = new int[sizepow];
        sparseIndexOddShift2D = new int[sizepow];
        sparseIndexEvenShift2D = new int[sizepow];
        sparseIndexShift2D = new int[sizepow];
        for (int i=0;i<sizepow;i++){
            sparseIndexOdd[i]=(i*2)+1;
            sparseIndexEven[i]=(i*2);
        }
        for (int i=0;i<param.size;i++){
            for (int ii=0;ii<param.size;ii++){
                sparseIndexShift2D[i*param.size+ii]=(((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size));
            }
        }
        for (int i=0;i<param.size;i++){
            //String s="";
            //String s2="";
            for (int ii=0;ii<param.size;ii++){
                
                //sparseIndexOdd[i*param.size+ii]=  i*param.size*2  +  (ii*2  +  1) ;
                //sparseIndexEven[i*param.size+ii]= i*param.size*2  +   ii*2 ;
                //s+=sparseIndexEven[i*param.size+ii]+"  ";
                sparseIndexOddShift2D[i*param.size+ii]=2*(((i+param.size/2)%param.size)*param.size  + (((ii)+param.size/2)%param.size))+1;
                sparseIndexEvenShift2D[i*param.size+ii]=2*(((i+param.size/2)%param.size)*param.size + (((ii)+param.size/2)%param.size));
                //s2+=sparseIndexEvenShift2D[i*param.size+ii]+"  ";
            }
            //IJ.log("s "+s);
            //IJ.log("                                                                "+s2);
        }
        
        
        for (int i=0,j=((param.size/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            //String s="";
            for (int ii=0,jj=((param.size/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                sparseIndexEvenShift2DOutput[i*sizeoutput2+ii]=sparseIndexEvenShift2D[j*param.size+jj];
                sparseIndexOddShift2DOutput[i*sizeoutput2+ii]=sparseIndexOddShift2D[j*param.size+jj];
                //s+=sparseIndexEvenShift2DOutput[i*param.size+ii]+"  ";
            }
            //IJ.log("ss "+s);
        }
        
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                
                sparseIndexEvenShift2DOutputNext[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2 + (((ii)+sizeoutput2/2)%sizeoutput2));
                sparseIndexOddShift2DOutputNext[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2 + (((ii)+sizeoutput2/2)%sizeoutput2))+1;
                
            }
        }
        //this.imshow(sizeoutput2, sparseIndexEvenShift2DOutputNext, "toto even");
        //this.imshow(sizeoutput2, sparseIndexOddShift2DOutputNext, "toto odd");
        
        //res[i][ii]=matrix[(i+w/2)%w][(ii+h/2)%h];
        
        
        
        
        
        
        host_phase= Pointer.to(phase);
        host_pupil= Pointer.to(pupil);
        host_kx= Pointer.to(kx);
        host_ky= Pointer.to(ky);
        host_kz= Pointer.to(kz);
        
        host_kz_oil= Pointer.to(kz_oil);
        host_kz_is_imaginary= Pointer.to(kz_is_imaginary);
        host_kz_oil_is_imaginary= Pointer.to(kz_oil_is_imaginary);
        
        host_sparseIndexOdd= Pointer.to(sparseIndexOdd);
        host_sparseIndexEven= Pointer.to(sparseIndexEven);
        host_sparseIndexOddDisk= Pointer.to(sparseIndexOddDisk);
        host_sparseIndexEvenDisk= Pointer.to(sparseIndexEvenDisk);
        host_sparseIndexOddShift2D= Pointer.to(sparseIndexOddShift2D);
        host_sparseIndexEvenShift2D= Pointer.to(sparseIndexEvenShift2D);
        host_sparseIndexShift2D= Pointer.to(sparseIndexShift2D);
        host_sparseIndexOutput= Pointer.to(sparseIndexOutput);
        host_sparseIndexOddShift2DOutput= Pointer.to(sparseIndexOddShift2DOutput);
        host_sparseIndexEvenShift2DOutput= Pointer.to(sparseIndexEvenShift2DOutput);
        host_sparseIndexEvenShift2DOutputNext= Pointer.to(sparseIndexEvenShift2DOutputNext);
        
        host_sparseIndexOddShift2DOutputNext= Pointer.to(sparseIndexOddShift2DOutputNext);
        
        
        
        
        
        deviceConstructor();
        
        setHost2Device(device_phase,host_phase,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_pupil,host_pupil,param.sizeDisk,Sizeof.DOUBLE);
        
        
        
        setHost2Device(device_kx,host_kx,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_ky,host_ky,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_kz_oil,host_kz_oil,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_kz,host_kz,param.sizeDisk,Sizeof.DOUBLE);
        
        setHost2Device(device_kz_is_imaginary,host_kz_is_imaginary,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_kz_oil_is_imaginary,host_kz_oil_is_imaginary,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_sparseIndexOutput,host_sparseIndexOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2DOutput,host_sparseIndexOddShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2DOutput,host_sparseIndexEvenShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2DOutputNext,host_sparseIndexEvenShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2DOutputNext,host_sparseIndexOddShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddDisk,host_sparseIndexOddDisk,param.sizeDisk,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenDisk,host_sparseIndexEvenDisk,param.sizeDisk,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2D,host_sparseIndexOddShift2D,sizepow,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2D,host_sparseIndexEvenShift2D,sizepow,Sizeof.INT);
        setHost2Device(device_sparseIndexOdd,host_sparseIndexOdd,sizepow,Sizeof.INT);
        setHost2Device(device_sparseIndexEven,host_sparseIndexEven,sizepow,Sizeof.INT);
        setHost2Device(device_sparseIndexShift2D,host_sparseIndexShift2D,sizepow,Sizeof.INT);
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 2 "+cudaResult+"   "+param.stream);}
        
        
        res2D = new double [sizeoutput2][sizeoutput2];
        res1D = new double [sizeoutput2*sizeoutput2];
        res = new double [sizeoutput2*sizeoutput2];
        
        
        
        
        plan = new cufftHandle();
        JCufft.cufftPlan2d(plan, param.size,param.size, cufftType.CUFFT_Z2Z);
        JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(param.stream));
        
        plan2 = new cufftHandle();
        JCufft.cufftPlan2d(plan2, sizeoutput2,sizeoutput2, cufftType.CUFFT_Z2Z);
        JCufft.cufftSetStream(plan2, MyCudaStream.getCudaStream_t(param.stream));
        
        handlecublas=MyCudaStream.getHandleCublas(param.stream);
        handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
        custream=MyCudaStream.getCUstream(param.stream);
        
        matrix=new double[sizeoutput2*sizeoutput2];
        
         cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 3 "+cudaResult+"   "+param.stream);}
        
        
        
        this.computeGaussianKernel(param.sigmaGaussianKernel);
        
         cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 4 "+cudaResult+"   "+param.stream);}
        
        
    }

    
    
    void deviceConstructor(){
        
        somResOne = new double [1];
        hostSomOne=Pointer.to(somResOne);
    
        device_phase = new Pointer();
        device_pupil = new Pointer();
        
        int cudaResult = JCuda.cudaMalloc(device_phase, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
            
        cudaResult = JCuda.cudaMalloc(device_pupil, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_kx = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kx, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_ky = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_ky, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_kz = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_kz_oil = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz_oil, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_kz_is_imaginary = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz_is_imaginary, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_kz_oil_is_imaginary = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz_oil_is_imaginary, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        device_phtmp = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_phtmp, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_tmp = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_realtmp = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_realtmp, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        device_imagtmp = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_imagtmp, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        device_psf = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psf, sizeoutput1 * sizeoutput1 *Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_psfFFT = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfFFT, 2*sizeoutput2 * sizeoutput2 *Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        /*device_psfsingle = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfsingle, sizeoutput2 * sizeoutput2 *Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }*/
        
        
        
        device_phtmpFull = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_phtmpFull, sizepow * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_tmpFull = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmpFull, sizepow * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_realtmpFull = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_realtmpFull, sizepow * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        device_imagtmpFull = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_imagtmpFull, sizepow * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        device_fftdata = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_fftdata, sizepow*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        
        device_sparseIndexOddDisk = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddDisk, param.sizeDisk * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexEvenShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_sparseIndexEvenShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_sparseIndexOddShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_sparseIndexOddShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexEvenDisk = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenDisk, param.sizeDisk * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        
        device_sparseIndexOdd = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOdd, sizepow * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexEven = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEven, sizepow * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        
        device_sparseIndexOddShift2D = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2D, sizepow * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexEvenShift2D = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2D, sizepow * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR PSFPhaseJCudaFastDouble malloc cuda");return ;
        }
        
        
        device_sparseIndexShift2D = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexShift2D, sizepow * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_gaussian=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian, sizeoutput2* sizeoutput2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR PSFPhaseJCudaFastDouble malloc cuda gaussian kernel");return ;
        }
        device_gaussian2bleSize=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian2bleSize, sizeoutput2* sizeoutput2*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR PSFPhaseJCudaFastDouble malloc cuda gaussian kernel");return ;
        }
    }
    
    
    
    void setHost2Device(Pointer device,Pointer host,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpyAsync(device, host, size * sizeElement, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR PSFPhaseJCudaFastDouble host2device cuda");return ;
        }
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF computePSF() 0 "+cudaResult);}
        
    }
    
    
    void setDevice2Host(Pointer host,Pointer device,int size,int sizeElement){
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF computePSF() 0 "+cudaResult);}
        
        int cudaResult = JCuda.cudaMemcpyAsync(host, device, size * sizeElement, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR PSFPhaseJCudaFastDouble device2host cuda "+cudaResult);return ;
        }
    }
    
    
    
    
    
    
    
    public void setSizeoutput(int sizeoutput){
        
        
        
        sizeoutput2=sizeoutput*2;
        sizeoutput1=sizeoutput;
        
    
        if (sizeoutput2>param.size){
            sizeoutput2=param.size;
            
        }
        
        if (sizeoutput2<=0){
            sizeoutput2=param.size;
        }
        
        if (sizeoutput1>param.size){
            sizeoutput1=param.size;
            
        }
        
        if (sizeoutput1<=0){
            sizeoutput1=param.size;
        }
        
        param.sizeoutput=sizeoutput1;
        
        
        JCuda.cudaFree(device_psf);
        JCuda.cudaFree(device_psfFFT);
        JCuda.cudaFree(device_sparseIndexOutput);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput);
        JCuda.cudaFree(device_gaussian);
        JCuda.cudaFree(device_gaussian2bleSize);
        
        
        
        
        
        sparseIndexOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutput = new int[sizeoutput2*sizeoutput2];
        
        for (int i=0,j=((param.size/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            for (int ii=0,jj=((param.size/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                sparseIndexOutput[i*sizeoutput2+ii]=j*param.size+jj;
            }
        }
        
        for (int i=0,j=((param.size/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            for (int ii=0,jj=((param.size/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                sparseIndexEvenShift2DOutput[i*sizeoutput2+ii]=sparseIndexEvenShift2D[j*param.size+jj];
                sparseIndexOddShift2DOutput[i*sizeoutput2+ii]=sparseIndexOddShift2D[j*param.size+jj];
            }
        }
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                
                sparseIndexEvenShift2DOutputNext[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2 + (((ii)+sizeoutput2/2)%sizeoutput2));
                sparseIndexOddShift2DOutputNext[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2 + (((ii)+sizeoutput2/2)%sizeoutput2))+1;
                
            }
        }
        
        host_sparseIndexOutput= Pointer.to(sparseIndexOutput);
        host_sparseIndexOddShift2DOutput= Pointer.to(sparseIndexOddShift2DOutput);
        host_sparseIndexEvenShift2DOutput= Pointer.to(sparseIndexEvenShift2DOutput);
        host_sparseIndexEvenShift2DOutputNext= Pointer.to(sparseIndexEvenShift2DOutputNext);
        host_sparseIndexOddShift2DOutputNext= Pointer.to(sparseIndexOddShift2DOutputNext);
        device_sparseIndexOutput = new Pointer();
        int cudaResult = JCuda.cudaMalloc(device_sparseIndexOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexEvenShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexEvenShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_sparseIndexOddShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        device_sparseIndexOddShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        setHost2Device(device_sparseIndexOutput,host_sparseIndexOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2DOutput,host_sparseIndexOddShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2DOutput,host_sparseIndexEvenShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2DOutputNext,host_sparseIndexEvenShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2DOutputNext,host_sparseIndexOddShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
        
        
        
        device_gaussian=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian, sizeoutput2* sizeoutput2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        device_gaussian2bleSize=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian2bleSize, sizeoutput2* sizeoutput2*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        device_psf = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psf, sizeoutput1 *sizeoutput1 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        device_psfFFT = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfFFT, 2*sizeoutput2 *sizeoutput2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }
        
        /*device_psfsingle = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfsingle, sizeoutput *sizeoutput * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }*/
        
        
        
        
        res2D = new double [sizeoutput2][sizeoutput2];
        res1D = new double [sizeoutput2*sizeoutput2];
        res = new double [sizeoutput2*sizeoutput2];
        gk.free();
        
        gk = new GaussianKernel_(sizeoutput2,sizeoutput,param.sigmaGaussianKernel,param.stream);
        //gk = new GaussianKernel(sizeoutput2,1,param.stream);
        
        matrix=new double[sizeoutput2*sizeoutput2];
        this.computeGaussianKernel(param.sigmaGaussianKernel);
    }
    
    
    
    
    
    
    
    //needed when nwat change
    public void resetKz(){
        double left=(param.nwat/(param.wavelength))*(param.nwat/(param.wavelength));
        double center=param.centerFourierImage;
        int id=0;
        for (int i=0;i<param.size;i++){
            for (int ii=0;ii<param.size;ii++){
                double disk=((((i-center)/(param.size*param.xystep))*((i-center)/(param.size*param.xystep))+((ii-center)/(param.size*param.xystep))*((ii-center)/(param.size*param.xystep))));
                
                if (disk<=param.ringsize){
                    
                    
                    //double iPositShifted=((i+size/2)%size);
                    //double iiPositShifted=((ii+size/2)%size);
                    
                    
                    
                    
                    double kxx=((i-center)/(param.xystep*(param.size)))*((i-center)/(param.xystep*(param.size)));
                    double kyy=((ii-center)/(param.xystep*(param.size)))*((ii-center)/(param.xystep*(param.size)));


                    //double left=4*Math.PI*Math.PI*k*k;

                    double right=(kxx+kyy);


                    if (left>right){
                        kz[id]=2*Math.PI*Math.sqrt(left-right);
                        kz_is_imaginary[id]=1;
                    }
                    else{
                        kz[id]=2*Math.PI*Math.sqrt(Math.abs(left-right));
                        kz_is_imaginary[id]=0;
                    }
                    
                    id++;
                }
                
                
                
            }
        }
        
        host_kz_is_imaginary= Pointer.to(kz_is_imaginary);
        host_kz= Pointer.to(kz);
        
        setHost2Device(device_kz,host_kz,param.sizeDisk,Sizeof.DOUBLE);
        setHost2Device(device_kz_is_imaginary,host_kz_is_imaginary,param.sizeDisk,Sizeof.DOUBLE);
        
    }
    
    
    public void setSigma(double sigmaInPixels){
        this.computeGaussianKernel(sigmaInPixels);
    }
            
    public void computeGaussianKernel(double sigmaInPixels){
        
        double sigmaFourier=(1./(sigmaInPixels))*sizeoutput2/(Math.PI*2);
        double sigpow=sigmaFourier*sigmaFourier;
        double c=sizeoutput2/2;
        double j,jj;
        for (int i=0;i<sizeoutput2;i++){
            j=i;
            for (int ii=0;ii<sizeoutput2;ii++){
                jj=ii;
                matrix[i*sizeoutput2+ii]=Math.exp(-.5*((j-c)*(j-c)+(jj-c)*(jj-c))/(sigpow));
            }
        }
        cudaMemcpyAsync(device_tmp, Pointer.to(matrix), sizeoutput2*sizeoutput2*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF comp gauss ker 1 "+cudaResult+"   "+param.stream);}
        //MyCudaStream.init(5);
        cusparseDgthr(handlecusparse, sizeoutput2*sizeoutput2, device_tmp, device_gaussian,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF comp gauss ker 2 "+cudaResult+"   "+param.stream);}
        
        cusparseDsctr(handlecusparse, sizeoutput2*sizeoutput2, device_tmp, device_sparseIndexEvenShift2DOutputNext,device_gaussian2bleSize, CUSPARSE_INDEX_BASE_ZERO);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF comp gauss ker 3 "+cudaResult+"   "+param.stream);}
        
        cusparseDsctr(handlecusparse, sizeoutput2*sizeoutput2, device_tmp, device_sparseIndexOddShift2DOutputNext,device_gaussian2bleSize, CUSPARSE_INDEX_BASE_ZERO);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF comp gauss ker 4 "+cudaResult+"   "+param.stream);}
        
    }
    
    public double [][] getPhase(){
        double [][] mat = new double [param.size][param.size];
        
        this.setDevice2Host(host_phase,device_phase, param.sizeDisk, Sizeof.DOUBLE);

        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=phase[i];
        }
            
        return mat;
    }
    
    public double [][] getPupil(){
        double [][] mat = new double [param.size][param.size];
        
        this.setDevice2Host( host_pupil,device_pupil, param.sizeDisk, Sizeof.DOUBLE);
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=pupil[i];
        }
        
        return mat;
    }
    
    public double [][] getKx(){
        this.setDevice2Host( host_kx,device_kx, param.sizeDisk, Sizeof.DOUBLE);
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kx[i];
        }
        return mat;
    }
    
    public double [][] getKy(){
        this.setDevice2Host( host_ky,device_ky, param.sizeDisk, Sizeof.DOUBLE);
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=ky[i];
        }
        return mat;
    }
    
    public double [][] getKz(){
        this.setDevice2Host( host_kz,device_kz, param.sizeDisk, Sizeof.DOUBLE);
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kz[i];
        }
        return mat;
    }
    
    
    
    
    
    public void updatePhase(double [][] mat){
            for (int i=0;i<param.sizeDisk;i++){
                phase[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
            }
            host_phase= Pointer.to(phase);
            setHost2Device(device_phase,host_phase,param.sizeDisk,Sizeof.DOUBLE);
        
        
    }
    
    
//    public void updatePupil(double [][][] mat){
//        for (int or=0;or<param;or++){
//            for (int i=0;i<param.sizeDisk;i++){
//                pupil[i]=mat[or][param.disk2D[i][0]][param.disk2D[i][1]];
//            }
//            host_pupil= Pointer.to(pupil);
//            setHost2Device(device_pupil[or],host_pupil,param.sizeDisk,Sizeof.DOUBLE);
//        }
//        
//    }
    
    
    
    
    public void updatePupil(double [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            pupil[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_pupil= Pointer.to(pupil);
        setHost2Device(device_pupil,host_pupil,param.sizeDisk,Sizeof.DOUBLE);
        
        
    }
    
    
    
    public void updatePhase(Pointer device_p){
        //this.device_phase=device_p;
        
        int cudaResult = JCuda.cudaMemcpyAsync(device_phase, device_p, param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR host2device cuda");return ;
        }
        
    }
    
//    public void updatePupil(Pointer [] device_p){
//        this.device_pupil=device_p;
//        
//    }
    
//    public void updatePupil(Pointer device_p){
//        this.device_pupil=device_p;
//        
//    }
    
    
    
//    public void updateKz(Pointer device_p){
//        this.device_kz=device_p;
//        //this.imshow(device_kz, "kz");
//        
//    }
    
    
    
    
    
    
    
    public void updateKx(double [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            kx[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_kx= Pointer.to(kx);
        setHost2Device(device_kx,host_kx,param.sizeDisk,Sizeof.DOUBLE);
        
    }
    
    
    
    
    
    public void updateKy(double [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            ky[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_ky= Pointer.to(ky);
        setHost2Device(device_ky,host_ky,param.sizeDisk,Sizeof.DOUBLE);
        
    }
    
    
    public void updateKz(double [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            kz[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_kz= Pointer.to(kz);
        setHost2Device(device_kz,host_kz,param.sizeDisk,Sizeof.DOUBLE);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    public void free(){
        
        JCuda.cudaFree(device_phase);
        
        JCuda.cudaFree(device_pupil);
        JCuda.cudaFree(device_kx);
        JCuda.cudaFree(device_ky);
        JCuda.cudaFree(device_kz);
        JCuda.cudaFree(device_kz_oil);
        
        JCuda.cudaFree(device_kz_is_imaginary);
        JCuda.cudaFree(device_kz_oil_is_imaginary);
        
        JCuda.cudaFree(device_phtmp);
        JCuda.cudaFree(device_realtmp);
        JCuda.cudaFree(device_tmp);
        JCuda.cudaFree(device_imagtmp);
        JCuda.cudaFree(device_phtmpFull);
        JCuda.cudaFree(device_realtmpFull);
        JCuda.cudaFree(device_tmpFull);
        JCuda.cudaFree(device_imagtmpFull);
        JCuda.cudaFree(device_fftdata);
        
        JCuda.cudaFree(device_sparseIndexOdd);
        JCuda.cudaFree(device_sparseIndexEven);
        JCuda.cudaFree(device_sparseIndexOddDisk);
        JCuda.cudaFree(device_sparseIndexEvenDisk);
        JCuda.cudaFree(device_sparseIndexOddShift2D);
        JCuda.cudaFree(device_sparseIndexEvenShift2D);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput);
        
        JCuda.cudaFree(device_sparseIndexShift2D);
        
        JCuda.cudaFree(device_gaussian);
        JCuda.cudaFree(device_gaussian2bleSize);
        
        JCuda.cudaFree(device_psfFFT);
        JCuda.cudaFree(device_psf);
        //JCuda.cudaFree(device_psfsingle);
        JCuda.cudaFree(device_sparseIndexOutput);
        
        JCufft.cufftDestroy(plan);
        plan=null;
        JCufft.cufftDestroy(plan2);
        plan2=null;
        gk.free();
        
    }
    
    
    
    
    
    
    
    private double [] shift2D(double[] matrix){
        
        double [] res = new double [param.size*param.size];
        for (int i=0;i<param.size;i++){
            for (int ii=0;ii<param.size;ii++){
                res[i*param.size+ii]=matrix[((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size)];
            }
        }
        return res;
        
    }
    
    
    
    
    
    public void imshow(int largersize,int [] im,String title){
        ImageProcessor ip = new FloatProcessor(im.length/largersize,largersize);
        for (int i=0;i<im.length/largersize;i++){
            for (int ii=0;ii<largersize;ii++){
                //ip.putPixelValue(i, ii, i*largersize+ii);
                ip.putPixelValue(i, ii, (float)im[i*largersize+ii]);
            }
        }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    public void imshow(int largersize,double [] im,String title){
        ImageProcessor ip = new FloatProcessor(im.length/largersize,largersize);
        for (int i=0;i<im.length/largersize;i++){
            for (int ii=0;ii<largersize;ii++){
                
                
                //ip.putPixelValue(i, ii, i*largersize+ii);
                ip.putPixelValue(i, ii, (float)im[i*largersize+ii]);
            }
        }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    public void imshow(int largersize,float [] im,String title){
        ImageProcessor ip = new FloatProcessor(im.length/largersize,largersize);
        for (int i=0;i<im.length/largersize;i++){
            for (int ii=0;ii<largersize;ii++){
                
                
                //ip.putPixelValue(i, ii, i*largersize+ii);
                ip.putPixelValue(i, ii, (float)im[i*largersize+ii]);
            }
        }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    public void imshow(int sizeTotal, int largersize,Pointer device,String title,String type){
        IJ.log("imshow is time consuming");
        if (type.startsWith("INT")){
            int [] vect = new int[sizeTotal];
            Pointer p =Pointer.to(vect);
            this.setDevice2Host(p,device,sizeTotal,Sizeof.INT);
            this.imshow(largersize,vect, title);
        }
        else if (type.startsWith("DOUBLE")){
            double [] vect = new double[sizeTotal];
            Pointer p =Pointer.to(vect);
            this.setDevice2Host(p,device,sizeTotal,Sizeof.DOUBLE);
            this.imshow(largersize,vect, title);
        }
        else if (type.startsWith("FLOAT")){
            float [] vect = new float[sizeTotal];
            Pointer p =Pointer.to(vect);
            this.setDevice2Host(p,device,sizeTotal,Sizeof.FLOAT);
            this.imshow(largersize,vect, title);
        }
    }
    
    
    
    public void imshow(Pointer deviceDiskShape,String title){
        double [] kk=new double [this.kz.length];
        Pointer hostkk = Pointer.to(kk);
        this.setDevice2Host( hostkk,deviceDiskShape, param.sizeDisk, Sizeof.DOUBLE);
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kk[i];
        }
        ImageShow.imshow(mat,title);
    }
    
    
    public void computeFourierMask(double z){
        
        JCuda.cudaMemsetAsync(device_imagtmp, 0, param.sizeDisk * Sizeof.DOUBLE, MyCudaStream.getCudaStream_t(param.stream));
        
        
        
        
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_phtmp, device_kz, z);
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_phtmp, device_phtmp, param.getweightZ());
        MyVecDouble.add(custream,param.sizeDisk, device_phtmp, device_phtmp, device_phase);
        
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_realtmp, device_kz_oil, 0);
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_realtmp, device_realtmp, param.getweightZ());
        MyVecDouble.sub(custream,param.sizeDisk, device_phase, device_phase, device_realtmp);
        
        MyVecDouble.cos(custream,param.sizeDisk, device_realtmp, device_phtmp);
        MyVecDouble.mul(custream,param.sizeDisk, device_realtmp, device_realtmp,device_pupil);
        
        
        MyVecDouble.sin(custream,param.sizeDisk, device_imagtmp, device_phtmp);
        MyVecDouble.mul(custream,param.sizeDisk, device_imagtmp, device_imagtmp,device_pupil);
        
        
        this.imshow(device_realtmp, "realtmp");
        this.imshow(device_imagtmp, "imagtmp");
        
        
        
    }
    
    
    
    
    
    
    public void updateSigmaGaussianKernel(double sigma){
        gk.setSigma(sigma);
        //IJ.log("sig set ok");
        this.setSigma(sigma);
    }
    
    
    
    
    
     
     
    public void computePSF(double x,double y,double zoil,double zwat){
        
        
        /*zwat+=param.Zfocus;
        zwat*=(param.nwat/(param.noil))*(param.nwat/(param.noil));*/
        
        
        //zoil*=param.noil/param.nwat;
        /*MyVecDouble.computePSF_phase(custream, param.sizeDisk, device_realtmp, device_imagtmp, device_kx, device_ky, device_kz, device_pupil, device_phase, x, y, z*param.weightZ);
        
        
        cudaResult = JCuda.cudaMemsetAsync(device_fftdata, 0, sizepow*2 * Sizeof.DOUBLE, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cudaMemset cuda 1");}
        
        
        cudaResult = cusparseDsctr(handlecusparse, param.sizeDisk, device_realtmp, device_sparseIndexEvenDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR dsctr cuda 3");}
        cudaResult = cusparseDsctr(handlecusparse, param.sizeDisk, device_imagtmp, device_sparseIndexOddDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR dsctr cuda 4");}
        */
        
        
        //IJ.log("zwat "+zwat+"  zoil "+zoil+"  "+param.nwat/(param.noil));
        
        JCuda.cudaMemsetAsync(device_fftdata, 0, sizepow*2 * Sizeof.DOUBLE, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cudaMemset cuda 1");}
        
        MyVecDouble.computePSF_phaseNwithOil(custream, param.sizeDisk, device_kx, device_ky, device_kz,device_kz_is_imaginary, device_kz_oil,device_kz_oil_is_imaginary, device_pupil, device_phase, x, y,zwat, zoil,device_sparseIndexEvenDisk,device_sparseIndexOddDisk,device_fftdata);
        
        //MyVecDouble.computePSF_phaseN(custream, param.sizeDisk, device_kx, device_ky, device_kz, device_pupil, device_phase, x,y,zoil       ,device_sparseIndexEvenDisk,device_sparseIndexOddDisk,device_fftdata);
        
        //this.imshow(param.size*param.size*2,param.size*2, device_fftdata, "device fft","DOUBLE");
        
        JCufft.cufftExecZ2Z(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 5");}
         
        //this.imshow(param.size*param.size*2,param.size*2, device_fftdata, "device fft","DOUBLE");
        
        //cudaResult = cusparseDgthr(handlecusparse, sizeoutput2*sizeoutput2, device_fftdata, device_realtmpFull,device_sparseIndexOddShift2DOutput, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 6");}
        //cudaResult = cusparseDgthr(handlecusparse, sizeoutput2*sizeoutput2, device_fftdata, device_imagtmpFull,device_sparseIndexEvenShift2DOutput, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 6");}
        
        
        
        //this.imshow(sizeoutput2*sizeoutput2,sizeoutput2, device_sparseIndexEvenShift2DOutput, "device_sparseIndexEvenShift2DOutput","INT");
        
        //this.imshow(sizeoutput2*sizeoutput2,sizeoutput2, device_sparseIndexEvenShift2DOutputNext, "device_sparseIndexEvenShift2DOutputNext","INT");
        //this.imshow(sizeoutput2*sizeoutput2,sizeoutput2, device_sparseIndexOddShift2DOutputNext, "device_sparseIndexOddShift2DOutputNext","INT");
        
        
        
        //MyVecDouble.computePSF_signalN(custream, sizeoutput2*sizeoutput2, device_psf, Math.sqrt(sizepow),device_sparseIndexEvenShift2DOutput,device_sparseIndexOddShift2DOutput,device_fftdata );
        
        
        
        MyVecDouble.computePSF_signalN2(custream, sizeoutput2*sizeoutput2, Math.sqrt(sizepow),device_sparseIndexEvenShift2DOutput,device_sparseIndexOddShift2DOutput,device_fftdata,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext,device_psfFFT);
        
        
        //this.imshow(sizeoutput2*sizeoutput2*2,sizeoutput2*2, device_psfFFT, "device fft 2","DOUBLE");
        //this.imshow(sizeoutput2*sizeoutput2,sizeoutput2, device_sparseIndexOddShift2DOutputNext, "AFT device_sparseIndexOddShift2DOutputNext","INT");
        
        
        //MyVecDouble.computePSF_signal(custream, sizeoutput2*sizeoutput2, device_psf, device_realtmpFull, device_imagtmpFull, Math.sqrt(sizepow));
       
        
        
        //FILTERING
        
       
        gk.filterPSF(device_psfFFT,device_psf);
        
        
        
        
        JCublas2.cublasDasum(handlecublas,sizeoutput1*sizeoutput1,device_psf,1,device_tmp);//norm sum=1
        
        
        JCuda.cudaMemcpyAsync(hostSomOne, device_tmp, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));;if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda PSF computePSF()2 "+cudaResult);}
        
        
        
        
        MyVecDouble.divScalar(custream,sizeoutput1*sizeoutput1, device_psf, device_psf,somResOne[0]);//square
        
        
        
        
    }
    
    
    
    
    
    
    
    
    public void computePSFNonOptim(double x,double y,double z){
        
        
        
        cudaResult = JCuda.cudaMemsetAsync(device_imagtmp, 0, param.sizeDisk * Sizeof.DOUBLE, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
        
        
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_phtmp, device_kz, z);
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_phtmp, device_phtmp, param.weightZ);
        MyVecDouble.add(custream,param.sizeDisk, device_phtmp, device_phtmp, device_phase);
        
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_tmp, device_kx, x);
        MyVecDouble.add(custream,param.sizeDisk, device_phtmp, device_phtmp, device_tmp);
        
        MyVecDouble.mulScalar(custream,param.sizeDisk, device_tmp, device_ky, y);
        MyVecDouble.add(custream,param.sizeDisk, device_phtmp, device_phtmp, device_tmp);
        
        
        
        MyVecDouble.cos(custream,param.sizeDisk, device_realtmp, device_phtmp);
        MyVecDouble.mul(custream,param.sizeDisk, device_realtmp, device_realtmp,device_pupil);
        
        
        MyVecDouble.sin(custream,param.sizeDisk, device_imagtmp, device_phtmp);
        MyVecDouble.mul(custream,param.sizeDisk, device_imagtmp, device_imagtmp,device_pupil);
        
        
        
        
        long x1=System.currentTimeMillis();
        
        
        cudaResult = JCuda.cudaMemsetAsync(device_fftdata, 0, sizepow*2 * Sizeof.DOUBLE, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cudaMemset cuda 1");}
        
        cudaResult = cusparseDsctr(handlecusparse, param.sizeDisk, device_realtmp, device_sparseIndexEvenDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR dsctr cuda 3");}
        cudaResult = cusparseDsctr(handlecusparse, param.sizeDisk, device_imagtmp, device_sparseIndexOddDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR dsctr cuda 4");}
        
        
        
        long x2=System.currentTimeMillis();
        
        
        cudaResult = JCufft.cufftExecZ2Z(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 5");}
        
        
        long x3=System.currentTimeMillis();
        
        
        
        MyVecDouble.divScalar(custream,sizepow*2, device_fftdata, device_fftdata, Math.sqrt(sizepow));//not necessary...TO REMOVE ?????
        
        
        //gather complex
        
        cudaResult = cusparseDgthr(handlecusparse, sizepow, device_fftdata, device_realtmpFull,device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 6");}
        cudaResult = cusparseDgthr(handlecusparse, sizepow, device_fftdata, device_imagtmpFull,device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 7");}
        
                
        
        cudaResult = cusparseDgthr(handlecusparse, sizepow, device_realtmpFull, device_tmpFull,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 8");}
        cudaResult = cusparseDgthr(handlecusparse, sizepow, device_imagtmpFull, device_phtmpFull,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 9");}
        
        
        //this.imshow(sizepow,(int)Math.sqrt(sizepow), device_tmpFull, "device_tmpFull","DOUBLE");
        
        MyVecDouble.mul(custream,sizepow, device_realtmpFull, device_tmpFull,device_tmpFull);//square
        MyVecDouble.mul(custream,sizepow, device_imagtmpFull, device_phtmpFull,device_phtmpFull);//square
        MyVecDouble.add(custream,sizepow, device_tmpFull, device_realtmpFull, device_imagtmpFull);
        
        
        
       cudaResult =  cusparseDgthr(handlecusparse, sizeoutput2*sizeoutput2, device_tmpFull, device_psf,device_sparseIndexOutput, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 10");}
        
        
        gk.filter(device_psf);
        
        
        
        JCublas2.cublasDasum(handlecublas,sizeoutput2*sizeoutput2,device_psf,1,device_tmp);//norm sum=1
        
        cudaResult=JCuda.cudaMemcpyAsync(hostSomOne, device_tmp, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));;if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda PSF computePSF()2 "+cudaResult);}
        
        MyVecDouble.divScalar(custream,sizeoutput2*sizeoutput2, device_psf, device_psf,somResOne[0]);//square
        
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF computePSF() 3 "+cudaResult);}
        
        
        //this.imshow(sizepow,(int)Math.sqrt(sizepow), device_psf, "device_psfAfter","DOUBLE");
        
    }
    
    
    
    
    
    
    public double [][] getPSF(){
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF getPSF([][]) "+cudaResult);}
        Pointer hostoutput=Pointer.to(res);//output point to res
        //IJ.log("res "+res.length+"   ");
        this.setDevice2Host(hostoutput, device_psf,sizeoutput1*sizeoutput1,Sizeof.DOUBLE);
        
        //2D conversion:
        
        res2D = new double [sizeoutput1][sizeoutput1];//already done in setsizeoutput
        for (int i=0;i<sizeoutput1;i++){
            for (int ii=0;ii<sizeoutput1;ii++){
                res2D[i][ii]=res[i*sizeoutput1+ii];
            }
        }
        //long tttttt=System.currentTimeMillis();
        
        return res2D;
    }
    
    
    
    
    
    public double [] getPSF1D(){
        Pointer hostoutput=Pointer.to(res);//output point to res
        //IJ.log("res "+res.length+"   ");
        this.setDevice2Host(hostoutput, device_psf,sizeoutput1*sizeoutput1,Sizeof.DOUBLE);
        
        //2D conversion:
        
        res1D = new double [sizeoutput1*sizeoutput1];
        for (int i=0;i<sizeoutput1;i++){
            for (int ii=0;ii<sizeoutput1;ii++){
                res1D[i*sizeoutput1+ii]=res[i*sizeoutput1+ii];
            }
        }
        //long tttttt=System.currentTimeMillis();
        
        return res1D;
    }
    
    public Pointer getPointerPSF(){
        return device_psf;
    }
    
    
    public Pointer getPointerkx(){
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 1 "+cudaResult+"   "+param.stream);}
        
        return device_kx;
    }
    
    public Pointer getPointerky(){
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 1 "+cudaResult+"   "+param.stream);}
        
        return device_ky;
    }
    
    public Pointer getPointerkz(){
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 1 "+cudaResult+"   "+param.stream);}
        return device_kz;
    }
    
    public Pointer getPointerkzOil(){
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF 1 "+cudaResult+"   "+param.stream);}
        return device_kz_oil;
    }
    
    
    public Pointer getPointerPupil(){
        return device_pupil;
    }
    
//    //weight the pupil term by term with the kernel pointer (use for filtering with gaussian kernel)
//    public void weightPupil(Pointer device_weight){
//        MyVecDouble.mul(custream,param.sizeDisk, device_pupil, device_pupil, device_weight);
//
//        MyVecDouble.mul(custream,param.sizeDisk, device_pupil, device_pupil, device_pupil);
//        
//        Pointer p = new Pointer();
//        JCuda.cudaMalloc(p, 1 * Sizeof.DOUBLE);
//        double [] som = new double [1];
//        JCublas2.cublasDasum(handlecublas,param.sizeDisk,device_pupil,1,p);//norm sum=1
//        JCuda.cudaMemcpyAsync(Pointer.to(som), p, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
//        
//        MyVecDouble.divScalar(custream,param.sizeDisk, device_pupil, device_pupil,som[0]);//square
//        MyVecDouble.sqrt(custream,param.sizeDisk, device_pupil, device_pupil);
//        
//    }
    
    
    
    
    //weight the pupil term by term with the kernel pointer (use for filtering with gaussian kernel)
//    public void normalizePupilToSum1(){
//        
//        MyVecDouble.mul(param.sizeDisk, device_pupil, device_pupil, device_pupil);
//        
//        Pointer p = new Pointer();
//        JCuda.cudaMalloc(p, 1 * Sizeof.DOUBLE);
//        double [] som = new double [1];
//        JCublas2.cublasDasum(handlecublas,param.sizeDisk,device_pupil,1,p);//norm sum=1
//        JCuda.cudaMemcpyAsync(Pointer.to(som), p, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);
//        
//        MyVecDouble.divScalar(param.sizeDisk, device_pupil, device_pupil,som);//square
//        MyVecDouble.sqrt(param.sizeDisk, device_pupil, device_pupil);
//    }
    
    
    
    
    public Pointer getPointerSparseIndexShift2D(){
        return this.device_sparseIndexShift2D;
    }
    
    
    
    
    
    public Pointer getPointerSparseIndexEven(){
        return this.device_sparseIndexEven;
    }
    
    public Pointer getPointerSparseIndexOdd(){
        return this.device_sparseIndexOdd;
    }
    
    
    public Pointer getPointerFourierMaskReal(){
        return device_realtmp;
    }
    public Pointer getPointerFourierMaskImag(){
        return device_imagtmp;
    }
    
    public void setPointerFourierMaskReal(Pointer data){
        JCuda.cudaMemcpyAsync(device_realtmp, data, param.sizeDisk * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));
    }
    
    public void setPointerFourierMaskImag(Pointer data){
        JCuda.cudaMemcpyAsync(device_imagtmp, data, param.sizeDisk * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));
    }
    
    
    public double [][][] shift2D(double[][][] matrix){
        int w=matrix[0].length;
        int h=matrix[0][0].length;
        if ((w%2==0)&&(h%2==0)){
            double [][][] res = new double [matrix.length][w][h];
            for (int u=0;u<matrix.length;u++){

                
                for (int i=0;i<w;i++){
                    for (int ii=0;ii<h;ii++){
                        res[u][i][ii]=matrix[u][(i+w/2)%w][(ii+h/2)%h];
                        //res[i][ii]=matrix[i][ii];
                    }
                }
                
                
            }
            return res;
        }
        else{
            IJ.log("problem image size has to be even");
            return null;
        }
    }
    
    
    
    
    public void save(String path){
        
        double [][][] im = new double [5][param.size][param.size];
        
        im[0]=this.getPhase();
        im[1]=this.getPupil();
        im[2]=this.getKx();
        im[3]=this.getKy();
        im[4]=this.getKz();
        //im=this.shift2D(im);
        
        for (int i=0;i<im[3].length;i++){
            for (int ii=0;ii<im[3][0].length;ii++){
                im[4][i][ii]*=param.weightZ;
            }
        }
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[k][i][ii]);
                }
            }
            if (k==0)
                ims.addSlice("phase",ip);
            else if (k==1)
                ims.addSlice("magnitude",ip);
            else if (k==2)
                ims.addSlice("kx",ip);
            else if (k==3)
                ims.addSlice("ky",ip);
            else if (k==4)
                ims.addSlice("kz",ip);
            else{
                IJ.log("curious add");
            }
        }
        ImagePlus imp=new ImagePlus("my_phase",ims);
        imp.show();
        if (path!=null){
            IJ.save(imp, path+".tif");
        }
        else{
            IJ.log("null path...save not possible");
        }
    
    }
    
    
    
}