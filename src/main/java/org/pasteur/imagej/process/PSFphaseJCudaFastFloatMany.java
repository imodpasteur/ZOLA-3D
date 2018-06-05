/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;

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
public class PSFphaseJCudaFastFloatMany {
    
    
    float [] position;
    
    Pointer device_alpha;
    Pointer device_beta;
    Pointer device_ones;
    Pointer device_sumRes;
    
    Pointer device_position;
    Pointer host_position;
    int currentPosition=0;
    
    cufftHandle plan;
    cufftHandle plan2;
    
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;
    CUstream custream;
            
            
    public Pointer device_phase;
    Pointer device_pupil;
    public Pointer device_kx;
    Pointer device_ky;
    Pointer device_kz;
    Pointer device_kz_oil;
    Pointer device_kz_is_imaginary;
    Pointer device_kz_oil_is_imaginary;
    //Pointer device_kz_right;
    //Pointer device_kz_left;
    
    
    Pointer device_psf;
    Pointer device_psf_double;
    Pointer device_psfFFT;
    //Pointer device_psfsingle;
    
    Pointer  device_gaussian;
    Pointer  device_gaussian2bleSize;
    
    
    
    //Pointer device_phtmpFull;
    //Pointer device_realtmpFull;
    Pointer device_tmpFull;
    //Pointer device_imagtmpFull;
    //Pointer device_phtmp;
    //Pointer device_realtmp;
    //Pointer device_tmp;
    //Pointer device_imagtmp;
    Pointer device_fftdata;
    //Pointer device_sparseIndexOdd;
    //Pointer device_sparseIndexEven;
    //Pointer device_sparseIndexShift2D;
    
    Pointer device_sparseIndexOddDisk;
    Pointer device_sparseIndexEvenDisk;
    //Pointer device_sparseIndexOddShift2D;
    //Pointer device_sparseIndexEvenShift2D;
    //Pointer device_sparseIndexOutput;
    
    Pointer device_sparseIndexOddShift2DOutput;
    Pointer device_sparseIndexEvenShift2DOutput;
    Pointer device_sparseIndexOddShift2DOutput1;
    Pointer device_sparseIndexEvenShift2DOutput1;
    Pointer device_sparseIndexEvenShift2DOutputNext;
    Pointer device_sparseIndexOddShift2DOutputNext;
    Pointer device_sparseIndexShift2DOutput;
    
    Pointer host_phase;
    Pointer host_pupil;
    Pointer host_kx;
    Pointer host_ky;
    Pointer host_kz;
    Pointer host_kz_oil;
    Pointer host_kz_is_imaginary;
    Pointer host_kz_oil_is_imaginary;
    
    
    Pointer host_sparseIndexOddDisk;
    Pointer host_sparseIndexEvenDisk;
    
    
    
    Pointer host_sparseIndexOddShift2DOutput;
    Pointer host_sparseIndexEvenShift2DOutput;
    Pointer host_sparseIndexEvenShift2DOutputNext;
    Pointer host_sparseIndexOddShift2DOutputNext;
    Pointer host_sparseIndexEvenShift2DOutput1;
    Pointer host_sparseIndexOddShift2DOutput1;
    
    Pointer host_sparseIndexShift2DOutput;
    
    //int [] sparseIndexOdd;
    //int [] sparseIndexEven;
    //int [] sparseIndexShift2D;
    int [] sparseIndexEvenDisk;
    int [] sparseIndexOddDisk;
    int [] sparseIndexEvenShift2D;
    int [] sparseIndexOddShift2D;
    //int [] sparseIndexOutput;
    int [] sparseIndexOddShift2DOutput;
    int [] sparseIndexEvenShift2DOutput;
    int [] sparseIndexEvenShift2DOutputNext;
    int [] sparseIndexOddShift2DOutputNext;
    int [] sparseIndexEvenShift2DOutput1;
    int [] sparseIndexOddShift2DOutput1;
    int [] sparseIndexShift2DOutput;
    float [] phase ;
    float [] pupil ;
    float [] kx ;
    float [] ky ;
    float [] kz ;
    float [] kz_oil ;
    float [] kz_is_imaginary ;
    float [] kz_oil_is_imaginary ;
    //double [] kz_right ;
    //double [] kz_left ;
    boolean success=true;
    int sizepow;
    
    float [] matrix;
    
    float [] somResOne;
    Pointer hostSomOne;
    
    float [][] res2D ;
    float [] res1D ;
    float [] res;
    PhaseRetrievalParametersDouble param;
    
    int sizeoutput1;
    int sizeoutput2;
    int cudaResult;
    int numberPSF;
    
    
    
    public PSFphaseJCudaFastFloatMany(PhaseRetrievalParametersDouble param,int numberPSF){
        sizeoutput1=param.sizeoutput;
        sizeoutput2=param.sizeoutput*2;
        if (sizeoutput2>param.size){
            sizeoutput2=param.size;
        }
        this.numberPSF=numberPSF;
        position=new float [4*numberPSF];
        host_position=Pointer.to(position);
        currentPosition=0;
        
        this.param=param;
    
        
        this.sizepow=param.size*param.size;
        
        
        
        
        phase = new float[param.sizeDisk];
        pupil = new float[param.sizeDisk];
        kx = new float[param.sizeDisk];
        ky = new float[param.sizeDisk];
        kz = new float[param.sizeDisk];
        kz_oil = new float[param.sizeDisk];
        kz_is_imaginary = new float[param.sizeDisk];
        kz_oil_is_imaginary = new float[param.sizeDisk];
        
        sparseIndexOddDisk = new int[param.sizeDisk];
        sparseIndexEvenDisk = new int[param.sizeDisk];
        
        
                
        //sparseIndexOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutput1 = new int[sizeoutput1*sizeoutput1];
        sparseIndexOddShift2DOutput1 = new int[sizeoutput1*sizeoutput1];
        for (int i=0,j=((param.size/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            for (int ii=0,jj=((param.size/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                //sparseIndexOutput[i*sizeoutput2+ii]=j*param.size+jj;
            }
        }
        
        double left=((param.nwat/(param.wavelength))*(param.nwat/(param.wavelength)));
        double left_oil=((param.noil/param.wavelength)*(param.noil/param.wavelength));
        
        float center=(float)param.centerFourierImage;
        int id=0;
        for (int i=0;i<param.size;i++){
            for (int ii=0;ii<param.size;ii++){
                float disk=(float)((((i-center)/(param.size*param.xystep))*((i-center)/(param.size*param.xystep))+((ii-center)/(param.size*param.xystep))*((ii-center)/(param.size*param.xystep))));
                
                if (disk<=param.ringsize){
                    
                    
                    //double iPositShifted=((i+size/2)%size);
                    //double iiPositShifted=((ii+size/2)%size);
                    
                    kx[id]=(float)((i-center)*(2*Math.PI)/(param.xystep*param.size));
                    ky[id]=(float)((ii-center)*(2*Math.PI)/(param.xystep*param.size));
                    
                    
                    double kxx=((i-center)/(param.xystep*(param.size)))*((i-center)/(param.xystep*(param.size)));
                    double kyy=((ii-center)/(param.xystep*(param.size)))*((ii-center)/(param.xystep*(param.size)));


                    //double left=4*Math.PI*Math.PI*k*k;

                    float right=(float)(kxx+kyy);

                    
                    if (left>right){
                        kz[id]=(float)(2*Math.PI*Math.sqrt(left-right));
                        kz_is_imaginary[id]=(float)1;
                    }
                    else{
                        kz[id]=0;//(float)(2*Math.PI*Math.sqrt(Math.abs(left-right)));
                        kz_is_imaginary[id]=(float)0;
                    }
                    
                    if (left_oil>right){
                        kz_oil[id]=(float)(2*Math.PI*Math.sqrt(left_oil-right));
                        kz_oil_is_imaginary[id]=(float)1;
                    }
                    else{
                        kz_oil[id]=0;//(float)(2*Math.PI*Math.sqrt(Math.abs(left_oil-right)));
                        kz_oil_is_imaginary[id]=(float)0;
                        
                    }
                    pupil[id]=(float)(1/(double)Math.sqrt(param.sizeDisk));//like that -> final sum=1
                    phase[id]=0;
                    
                    
                    
                    sparseIndexOddDisk[id]=((((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size))*2)+1;
                    sparseIndexEvenDisk[id]=((((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size))*2);
                    id++;
                }
                
                
                
            }
        }
        
        
        
        
        
        
        
        
        
        
        
        //sparseIndexOdd = new int[sizepow];
        //sparseIndexEven = new int[sizepow];
        sparseIndexOddShift2D = new int[sizepow];
        sparseIndexEvenShift2D = new int[sizepow];
        int [] sparseIndexOddShift2D1 = new int[sizeoutput2*sizeoutput2];
        int [] sparseIndexEvenShift2D1 = new int[sizeoutput2*sizeoutput2];
        //sparseIndexShift2D = new int[sizepow];
        sparseIndexShift2DOutput = new int[sizeoutput2*sizeoutput2];
        for (int i=0;i<sizepow;i++){
            //sparseIndexOdd[i]=(i*2)+1;
            //sparseIndexEven[i]=(i*2);
        }
        for (int i=0;i<param.size;i++){
            for (int ii=0;ii<param.size;ii++){
                //sparseIndexShift2D[i*param.size+ii]=(((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size));
            }
        }
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                sparseIndexShift2DOutput[i*sizeoutput2+ii]=(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2+((ii+sizeoutput2/2)%sizeoutput2));
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
        
        
        
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                sparseIndexOddShift2D1[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2  + (((ii)+sizeoutput2/2)%sizeoutput2))+1;
                sparseIndexEvenShift2D1[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2 + (((ii)+sizeoutput2/2)%sizeoutput2));
            }
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
        
        for (int i=0,j=((sizeoutput2/2)-(sizeoutput1/2));i<sizeoutput1;i++,j++){
            //String s="";
            for (int ii=0,jj=((sizeoutput2/2)-(sizeoutput1/2));ii<sizeoutput1;ii++,jj++){
                sparseIndexEvenShift2DOutput1[i*sizeoutput1+ii]=sparseIndexEvenShift2D1[j*sizeoutput2+jj];
                sparseIndexOddShift2DOutput1[i*sizeoutput1+ii]=sparseIndexOddShift2D1[j*sizeoutput2+jj];
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        host_phase= Pointer.to(phase);
        host_pupil= Pointer.to(pupil);
        host_kx= Pointer.to(kx);
        host_ky= Pointer.to(ky);
        host_kz= Pointer.to(kz);
        host_kz_oil= Pointer.to(kz_oil);
        host_kz_is_imaginary= Pointer.to(kz_is_imaginary);
        host_kz_oil_is_imaginary= Pointer.to(kz_oil_is_imaginary);
        //host_kz_left= Pointer.to(kz_left);
        //host_kz_right= Pointer.to(kz_right);
        
        //host_sparseIndexOdd= Pointer.to(sparseIndexOdd);
        //host_sparseIndexEven= Pointer.to(sparseIndexEven);
        host_sparseIndexOddDisk= Pointer.to(sparseIndexOddDisk);
        host_sparseIndexEvenDisk= Pointer.to(sparseIndexEvenDisk);
        //host_sparseIndexOddShift2D= Pointer.to(sparseIndexOddShift2D);
        //host_sparseIndexEvenShift2D= Pointer.to(sparseIndexEvenShift2D);
        //host_sparseIndexShift2D= Pointer.to(sparseIndexShift2D);
        host_sparseIndexShift2DOutput= Pointer.to(sparseIndexShift2DOutput);
        //host_sparseIndexOutput= Pointer.to(sparseIndexOutput);
        host_sparseIndexOddShift2DOutput= Pointer.to(sparseIndexOddShift2DOutput);
        host_sparseIndexEvenShift2DOutput= Pointer.to(sparseIndexEvenShift2DOutput);
        host_sparseIndexEvenShift2DOutputNext= Pointer.to(sparseIndexEvenShift2DOutputNext);
        
        host_sparseIndexOddShift2DOutputNext= Pointer.to(sparseIndexOddShift2DOutputNext);
        
        host_sparseIndexEvenShift2DOutput1= Pointer.to(sparseIndexEvenShift2DOutput1);
        
        host_sparseIndexOddShift2DOutput1= Pointer.to(sparseIndexOddShift2DOutput1);
        
        
        
        
        
        boolean construct=deviceConstructor();
        
        if (construct){
            setHost2Device(device_phase,host_phase,param.sizeDisk,Sizeof.FLOAT);
            setHost2Device(device_pupil,host_pupil,param.sizeDisk,Sizeof.FLOAT);



            setHost2Device(device_kx,host_kx,param.sizeDisk,Sizeof.FLOAT);
            setHost2Device(device_ky,host_ky,param.sizeDisk,Sizeof.FLOAT);
            setHost2Device(device_kz,host_kz,param.sizeDisk,Sizeof.FLOAT);
            setHost2Device(device_kz_oil,host_kz_oil,param.sizeDisk,Sizeof.FLOAT);

            setHost2Device(device_kz_is_imaginary,host_kz_is_imaginary,param.sizeDisk,Sizeof.FLOAT);
            setHost2Device(device_kz_oil_is_imaginary,host_kz_oil_is_imaginary,param.sizeDisk,Sizeof.FLOAT);

            //setHost2Device(device_kz_left,host_kz_left,param.sizeDisk,Sizeof.FLOAT);
            //setHost2Device(device_kz_right,host_kz_right,param.sizeDisk,Sizeof.FLOAT);
            //setHost2Device(device_sparseIndexOutput,host_sparseIndexOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
            setHost2Device(device_sparseIndexOddShift2DOutput,host_sparseIndexOddShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
            setHost2Device(device_sparseIndexEvenShift2DOutput,host_sparseIndexEvenShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
            setHost2Device(device_sparseIndexEvenShift2DOutputNext,host_sparseIndexEvenShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
            setHost2Device(device_sparseIndexOddShift2DOutputNext,host_sparseIndexOddShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
            setHost2Device(device_sparseIndexEvenShift2DOutput1,host_sparseIndexEvenShift2DOutput1,sizeoutput1*sizeoutput1,Sizeof.INT);
            setHost2Device(device_sparseIndexOddShift2DOutput1,host_sparseIndexOddShift2DOutput1,sizeoutput1*sizeoutput1,Sizeof.INT);
            setHost2Device(device_sparseIndexOddDisk,host_sparseIndexOddDisk,param.sizeDisk,Sizeof.INT);
            setHost2Device(device_sparseIndexEvenDisk,host_sparseIndexEvenDisk,param.sizeDisk,Sizeof.INT);
            //setHost2Device(device_sparseIndexOddShift2D,host_sparseIndexOddShift2D,sizepow,Sizeof.INT);
            //setHost2Device(device_sparseIndexEvenShift2D,host_sparseIndexEvenShift2D,sizepow,Sizeof.INT);
            //setHost2Device(device_sparseIndexOdd,host_sparseIndexOdd,sizepow,Sizeof.INT);
            //setHost2Device(device_sparseIndexEven,host_sparseIndexEven,sizepow,Sizeof.INT);
            //setHost2Device(device_sparseIndexShift2D,host_sparseIndexShift2D,sizepow,Sizeof.INT);
            
            
            
            setHost2Device(device_sparseIndexShift2DOutput,host_sparseIndexShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);

            res2D = new float [sizeoutput2][sizeoutput2];
            res1D = new float [sizeoutput2*sizeoutput2];
            res = new float [sizeoutput2*sizeoutput2];


            plan = new cufftHandle();
            plan2 = new cufftHandle();
            int [] size= new int[2];
            size[0]=param.size;
            size[1]=param.size;
            int [] sizeoutput= new int[2];
            sizeoutput[0]=sizeoutput2;
            sizeoutput[1]=sizeoutput2;
            //JCufft.cufftPlan2d(plan, param.size,param.size, cufftType.CUFFT_Z2Z);
            JCufft.cufftPlanMany(plan, 2,size, size,1,param.size*param.size, size,1,param.size*param.size,cufftType.CUFFT_C2C,numberPSF);
            JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(param.stream));

            JCufft.cufftPlanMany(plan2, 2,sizeoutput, sizeoutput,1,sizeoutput2*sizeoutput2, sizeoutput,1,sizeoutput2*sizeoutput2,cufftType.CUFFT_C2C,numberPSF);
            JCufft.cufftSetStream(plan2, MyCudaStream.getCudaStream_t(param.stream));

            handlecublas=MyCudaStream.getHandleCublas(param.stream);
            handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
            custream=MyCudaStream.getCUstream(param.stream);

            device_sumRes = new Pointer();
            cudaMalloc(device_sumRes, numberPSF * Sizeof.FLOAT);
            device_alpha = new Pointer();
            cudaMalloc(device_alpha, 1 * Sizeof.FLOAT);
            cudaMemcpyAsync(device_alpha, Pointer.to(new float[]{1}), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
            device_beta = new Pointer();
            cudaMalloc(device_beta, 1 * Sizeof.FLOAT);
            cudaMemcpyAsync(device_beta, Pointer.to(new float[]{0}), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
            float [] ones=new float[sizeoutput1*sizeoutput1*numberPSF];
            for (int i=0;i<sizeoutput1*sizeoutput1*numberPSF;i++){
                ones[i]=1;
            }
            Pointer h_ones=Pointer.to(ones);
            device_ones = new Pointer();
            JCuda.cudaMalloc(device_ones, sizeoutput1*sizeoutput1*numberPSF *Sizeof.FLOAT);
            //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
            cudaMemcpyAsync(device_ones, h_ones, sizeoutput1*sizeoutput1*numberPSF*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));

            matrix=new float[sizeoutput2*sizeoutput2];
            this.computeGaussianKernel(param.sigmaGaussianKernel);

        }
        
        
        
        
        
    }

    
    
    boolean deviceConstructor(){
        
        
        long sizeMemoryNeeded=(long)(param.size*param.size*5)+(long)(numberPSF * Sizeof.FLOAT+sizeoutput2*sizeoutput2*numberPSF *Sizeof.FLOAT)+(long)(2*Sizeof.FLOAT+sizeoutput2* sizeoutput2*2 * Sizeof.FLOAT)+(long)(sizeoutput2* sizeoutput2 * Sizeof.FLOAT)+(long)(sizeoutput2* sizeoutput2 * Sizeof.INT+param.sizeDisk * Sizeof.INT)+(long)(4*sizeoutput2* sizeoutput2* Sizeof.INT)+(long)(param.sizeDisk * Sizeof.INT+numberPSF*sizepow*2 * Sizeof.FLOAT+sizepow * Sizeof.FLOAT)+(long)(3*sizeoutput2 * sizeoutput2*numberPSF *Sizeof.FLOAT)+(long)(sizeoutput2 * sizeoutput2*numberPSF *Sizeof.DOUBLE)+(long)(8*param.sizeDisk * Sizeof.FLOAT+4*this.numberPSF * Sizeof.FLOAT);
        
        long [] memtot=new long[1];
        long [] memfree=new long[1];
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        
        if (memfree[0]-100<(sizeMemoryNeeded)){
            IJ.log("Problem: not enough GPU memory to allocate memory");
            return false;
        }
        
        
        somResOne = new float [numberPSF];
        hostSomOne=Pointer.to(somResOne);
    
        device_phase = new Pointer();
        device_pupil = new Pointer();
        
        this.device_position=new Pointer();
        
        int cudaResult = JCuda.cudaMalloc(device_position, 4*this.numberPSF * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){success=false;IJ.log("ERROR malloc cuda");return false;}
        
        cudaResult = JCuda.cudaMalloc(device_phase, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return false;
        }
            
        cudaResult = JCuda.cudaMalloc(device_pupil, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_kx = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kx, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_ky = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_ky, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_kz = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_kz_oil = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz_oil, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        device_kz_is_imaginary = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz_is_imaginary, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_kz_oil_is_imaginary = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_kz_oil_is_imaginary, param.sizeDisk * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
//        device_kz_right = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_kz_right, param.sizeDisk * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        device_kz_left = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_kz_left, param.sizeDisk * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        
//        device_phtmp = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_phtmp, param.sizeDisk * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        device_tmp = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_tmp, param.sizeDisk * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        device_realtmp = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_realtmp, param.sizeDisk * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        
//        device_imagtmp = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_imagtmp, param.sizeDisk * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
        
        
        device_psf = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psf, sizeoutput1 * sizeoutput1*numberPSF *Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_psf_double = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psf_double, sizeoutput1 * sizeoutput1*numberPSF *Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_psfFFT = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfFFT, 2*sizeoutput2 * sizeoutput2*numberPSF *Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        /*device_psfsingle = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfsingle, sizeoutput2 * sizeoutput2 *Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return ;
        }*/
        
        
        
//        device_phtmpFull = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_phtmpFull, sizepow * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
        
        
        device_tmpFull = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmpFull, sizepow * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
//        
//        device_realtmpFull = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_realtmpFull, sizepow * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        
//        device_imagtmpFull = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_imagtmpFull, sizepow * Sizeof.FLOAT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
        
        
        
        device_fftdata = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_fftdata, numberPSF*sizepow*2 * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return false ;
        }
        
        
        
        device_sparseIndexOddDisk = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddDisk, param.sizeDisk * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
//        device_sparseIndexOutput = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_sparseIndexOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
        
        
        
        device_sparseIndexEvenShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_sparseIndexEvenShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_sparseIndexOddShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        device_sparseIndexEvenShift2DOutput1 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutput1, sizeoutput1* sizeoutput1* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_sparseIndexOddShift2DOutput1 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutput1, sizeoutput1* sizeoutput1* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        
        device_sparseIndexOddShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        device_sparseIndexEvenDisk = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenDisk, param.sizeDisk * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        
        
//        device_sparseIndexOdd = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_sparseIndexOdd, sizepow * Sizeof.INT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        device_sparseIndexEven = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_sparseIndexEven, sizepow * Sizeof.INT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        
//        device_sparseIndexOddShift2D = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2D, sizepow * Sizeof.INT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        device_sparseIndexEvenShift2D = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2D, sizepow * Sizeof.INT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
//        
//        device_sparseIndexShift2D = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_sparseIndexShift2D, sizepow * Sizeof.INT);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            success=false;IJ.log("ERROR malloc cuda");return ;
//        }
//        
        
        
        device_sparseIndexShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexShift2DOutput, sizeoutput2* sizeoutput2 * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        device_gaussian=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian,sizeoutput2* sizeoutput2 * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return  false;
        }
        
        device_gaussian2bleSize=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian2bleSize, sizeoutput2* sizeoutput2*2 * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return  false;
        }
        
        return true;
    }
    
    
    
    void setHost2Device(Pointer device,Pointer host,int size,int sizeElement){
        
        int cudaResult = JCuda.cudaMemcpyAsync(device, host, size * sizeElement, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR host2device cuda");return ;
        }
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF computePSF() 0 "+cudaResult);}
        
    }
    
    
    void setDevice2Host(Pointer host,Pointer device,int size,int sizeElement){
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF computePSF() 0 "+cudaResult);}
        
        int cudaResult = JCuda.cudaMemcpyAsync(host, device, size * sizeElement, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR device2host cuda "+cudaResult);return ;
        }
    }
    
    
    
    
    
    
    
    public boolean setSizeoutput(int sizeoutput){
        
        
        long sizeMemoryRemoved=(long)(this.sizeoutput2 *this.sizeoutput2*numberPSF * Sizeof.DOUBLE)+(long)(3*this.sizeoutput2 *this.sizeoutput2*numberPSF * Sizeof.FLOAT)+(long)(this.sizeoutput2* this.sizeoutput2*2 * Sizeof.FLOAT)+(long)(this.sizeoutput2* this.sizeoutput2 * Sizeof.FLOAT)+(long)(this.sizeoutput2*this.sizeoutput2*numberPSF *Sizeof.FLOAT)+(long)(5*this.sizeoutput2* this.sizeoutput2 * Sizeof.INT);
        long sizeMemoryNeeded=(long)(sizeoutput*2 *sizeoutput*2*numberPSF * Sizeof.DOUBLE)+(long)(3*sizeoutput*2 *sizeoutput*2*numberPSF * Sizeof.FLOAT)+(long)(sizeoutput*2* sizeoutput*2*2 * Sizeof.FLOAT)+(long)(sizeoutput*2* sizeoutput*2 * Sizeof.FLOAT)+(long)(sizeoutput*2*sizeoutput*2*numberPSF *Sizeof.FLOAT)+(long)(5*sizeoutput*2* sizeoutput*2 * Sizeof.INT);
        
        long [] memtot=new long[1];
        long [] memfree=new long[1];
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        
        if ((memfree[0]-100+sizeMemoryRemoved)<(sizeMemoryNeeded)){
            IJ.log("Problem: not enough GPU memory to allocate memory");
            return false;
        }
        
        
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
        
        
        JCuda.cudaFree(device_psf_double);
        JCuda.cudaFree(device_psf);
        JCuda.cudaFree(device_psfFFT);
        
        //JCuda.cudaFree(device_sparseIndexOutput);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput1);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput1);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput);
        JCuda.cudaFree(device_sparseIndexShift2DOutput);
        JCuda.cudaFree(device_gaussian);
        JCuda.cudaFree(device_gaussian2bleSize);
        JCuda.cudaFree(device_ones);
        
        
        device_sparseIndexShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexShift2DOutput, sizeoutput2* sizeoutput2 * Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return false;
        }
        
        
        sparseIndexEvenShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutputNext = new int[sizeoutput2*sizeoutput2];
        sparseIndexEvenShift2DOutput1 = new int[sizeoutput1*sizeoutput1];
        sparseIndexOddShift2DOutput1 = new int[sizeoutput1*sizeoutput1];
        sparseIndexEvenShift2DOutput = new int[sizeoutput2*sizeoutput2];
        sparseIndexOddShift2DOutput = new int[sizeoutput2*sizeoutput2];
        
        
        
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                sparseIndexShift2DOutput[i*sizeoutput2+ii]=(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2+((ii+sizeoutput2/2)%sizeoutput2));
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
        
        int [] sparseIndexOddShift2D1 = new int[sizeoutput2*sizeoutput2];
        int [] sparseIndexEvenShift2D1 = new int[sizeoutput2*sizeoutput2];
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                sparseIndexOddShift2D1[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2  + (((ii)+sizeoutput2/2)%sizeoutput2))+1;
                sparseIndexEvenShift2D1[i*sizeoutput2+ii]=2*(((i+sizeoutput2/2)%sizeoutput2)*sizeoutput2 + (((ii)+sizeoutput2/2)%sizeoutput2));
            }
        }
        
        for (int i=0,j=((sizeoutput2/2)-(sizeoutput1/2));i<sizeoutput1;i++,j++){
            for (int ii=0,jj=((sizeoutput2/2)-(sizeoutput1/2));ii<sizeoutput1;ii++,jj++){
                sparseIndexEvenShift2DOutput1[i*sizeoutput1+ii]=sparseIndexEvenShift2D1[j*sizeoutput2+jj];
                sparseIndexOddShift2DOutput1[i*sizeoutput1+ii]=sparseIndexOddShift2D1[j*sizeoutput2+jj];
            }
        }
        
        host_sparseIndexShift2DOutput= Pointer.to(sparseIndexShift2DOutput);
        
        host_sparseIndexOddShift2DOutput= Pointer.to(sparseIndexOddShift2DOutput);
        host_sparseIndexEvenShift2DOutput= Pointer.to(sparseIndexEvenShift2DOutput);
        host_sparseIndexEvenShift2DOutputNext= Pointer.to(sparseIndexEvenShift2DOutputNext);
        host_sparseIndexOddShift2DOutputNext= Pointer.to(sparseIndexOddShift2DOutputNext);
        host_sparseIndexEvenShift2DOutput1= Pointer.to(sparseIndexEvenShift2DOutput1);
        host_sparseIndexOddShift2DOutput1= Pointer.to(sparseIndexOddShift2DOutput1);

        
        
        device_sparseIndexEvenShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return false ;
        }
        
        device_sparseIndexEvenShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_sparseIndexOddShift2DOutputNext = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutputNext, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        device_sparseIndexEvenShift2DOutput1 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexEvenShift2DOutput1, sizeoutput1* sizeoutput1* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_sparseIndexOddShift2DOutput1 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutput1, sizeoutput1* sizeoutput1* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        device_sparseIndexOddShift2DOutput = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_sparseIndexOddShift2DOutput, sizeoutput2* sizeoutput2* Sizeof.INT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return false ;
        }
        
        
        
        
        setHost2Device(device_sparseIndexOddShift2DOutput,host_sparseIndexOddShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2DOutput,host_sparseIndexEvenShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexEvenShift2DOutputNext,host_sparseIndexEvenShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2DOutputNext,host_sparseIndexOddShift2DOutputNext,sizeoutput2*sizeoutput2,Sizeof.INT);
        
        setHost2Device(device_sparseIndexEvenShift2DOutput1,host_sparseIndexEvenShift2DOutput1,sizeoutput2*sizeoutput2,Sizeof.INT);
        setHost2Device(device_sparseIndexOddShift2DOutput1,host_sparseIndexOddShift2DOutput1,sizeoutput2*sizeoutput2,Sizeof.INT);
        
        setHost2Device(device_sparseIndexShift2DOutput,host_sparseIndexShift2DOutput,sizeoutput2*sizeoutput2,Sizeof.INT);
        
        
        
        float [] ones=new float[sizeoutput1*sizeoutput1*numberPSF];
        for (int i=0;i<sizeoutput1*sizeoutput1*numberPSF;i++){
            ones[i]=1;
        }
        Pointer h_ones=Pointer.to(ones);
        device_ones = new Pointer();
        JCuda.cudaMalloc(device_ones, sizeoutput1*sizeoutput1*numberPSF *Sizeof.FLOAT);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_ones, h_ones, sizeoutput1*sizeoutput1*numberPSF*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        
        
        
        device_gaussian=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian, sizeoutput2* sizeoutput2 * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return  false;
        }
        device_gaussian2bleSize=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian2bleSize, sizeoutput2* sizeoutput2*2 * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return  false;
        }
        
        
        
        device_psf = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psf, sizeoutput1 *sizeoutput1*numberPSF * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_psf_double = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psf_double, sizeoutput1 *sizeoutput1*numberPSF * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        device_psfFFT = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_psfFFT, 2*sizeoutput2 *sizeoutput2*numberPSF * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            success=false;IJ.log("ERROR malloc cuda");return  false;
        }
        
        
        
        
        
        
        res2D = new float [sizeoutput1][sizeoutput1];
        res1D = new float [sizeoutput1*sizeoutput1];
        res = new float [sizeoutput1*sizeoutput1];
        
        matrix=new float[sizeoutput2*sizeoutput2];
        this.computeGaussianKernel(param.sigmaGaussianKernel);
        
        JCufft.cufftDestroy(plan2);
        plan2=null;
        
        plan2 = new cufftHandle();
        int [] size= new int[2];
        size[0]=param.size;
        size[1]=param.size;
        int [] sizeoutputi= new int[2];
        sizeoutputi[0]=sizeoutput2;
        sizeoutputi[1]=sizeoutput2;
        
        JCufft.cufftPlanMany(plan2, 2,sizeoutputi, sizeoutputi,1,sizeoutput2*sizeoutput2, sizeoutputi,1,sizeoutput2*sizeoutput2,cufftType.CUFFT_C2C,numberPSF);
        JCufft.cufftSetStream(plan2, MyCudaStream.getCudaStream_t(param.stream));
        return true;
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
                        kz[id]=(float)(2*Math.PI*Math.sqrt(left-right));
                        kz_is_imaginary[id]=1;
                    }
                    else{
                        kz[id]=(float)(2*Math.PI*Math.sqrt(Math.abs(left-right)));
                        kz_is_imaginary[id]=0;
                    }
                    
                    id++;
                }
                
                
                
            }
        }
        
        host_kz_is_imaginary= Pointer.to(kz_is_imaginary);
        host_kz= Pointer.to(kz);
        
        setHost2Device(device_kz,host_kz,param.sizeDisk,Sizeof.FLOAT);
        setHost2Device(device_kz_is_imaginary,host_kz_is_imaginary,param.sizeDisk,Sizeof.FLOAT);
        
    }
    
    public void setSigma(double sigmaInPixels){
        this.computeGaussianKernel(sigmaInPixels);
    }
            
    public void computeGaussianKernel(double sigmaInPixels){
        //IJ.log("update sigma "+sigmaInPixels);
        double sigmaFourier=(1./(sigmaInPixels))*sizeoutput2/(Math.PI*2);
        double sigpow=sigmaFourier*sigmaFourier;
        double c=sizeoutput2/2;
        double j,jj;
        for (int i=0;i<sizeoutput2;i++){
            j=i;
            for (int ii=0;ii<sizeoutput2;ii++){
                jj=ii;
                matrix[i*sizeoutput2+ii]=(float)(Math.exp(-.5*((j-c)*(j-c)+(jj-c)*(jj-c))/(sigpow)));
            }
        }
        cudaMemcpyAsync(device_tmpFull, Pointer.to(matrix), sizeoutput2*sizeoutput2*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(param.stream));
        
        cusparseSgthr(handlecusparse, sizeoutput2*sizeoutput2, device_tmpFull, device_gaussian,device_sparseIndexShift2DOutput, CUSPARSE_INDEX_BASE_ZERO);
        
        cusparseSsctr(handlecusparse, sizeoutput2*sizeoutput2, device_tmpFull, device_sparseIndexEvenShift2DOutputNext,device_gaussian2bleSize, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSsctr(handlecusparse, sizeoutput2*sizeoutput2, device_tmpFull, device_sparseIndexOddShift2DOutputNext,device_gaussian2bleSize, CUSPARSE_INDEX_BASE_ZERO);
        
        
        
    }
    
    
    
    
    public float [][] getPhase(){
        float [][] mat = new float [param.size][param.size];
        
        this.setDevice2Host(host_phase,device_phase, param.sizeDisk, Sizeof.FLOAT);

        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=phase[i];
        }
            
        return mat;
    }
    
    
    public double [][] getPhaseDouble(){
        double [][] mat = new double [param.size][param.size];
        
        this.setDevice2Host(host_phase,device_phase, param.sizeDisk, Sizeof.FLOAT);

        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=phase[i];
        }
            
        return mat;
    }
    
    public float [][] getPupil(){
        float [][] mat = new float [param.size][param.size];
        
        this.setDevice2Host( host_pupil,device_pupil, param.sizeDisk, Sizeof.FLOAT);
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=pupil[i];
        }
        
        return mat;
    }
    
    public float [][] getKx(){
        this.setDevice2Host( host_kx,device_kx, param.sizeDisk, Sizeof.FLOAT);
        float [][] mat = new float [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kx[i];
        }
        return mat;
    }
    
    public float [][] getKy(){
        this.setDevice2Host( host_ky,device_ky, param.sizeDisk, Sizeof.FLOAT);
        float [][] mat = new float [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=ky[i];
        }
        return mat;
    }
    
    public float [][] getKz(){
        this.setDevice2Host( host_kz,device_kz, param.sizeDisk, Sizeof.FLOAT);
        float [][] mat = new float [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kz[i];
        }
        return mat;
    }
    
    
    
    
    
    public void updatePhase(float [][] mat){
            for (int i=0;i<param.sizeDisk;i++){
                phase[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
            }
            host_phase= Pointer.to(phase);
            setHost2Device(device_phase,host_phase,param.sizeDisk,Sizeof.FLOAT);
        
        
    }
    
    
//    public void updatePupil(double [][][] mat){
//        for (int or=0;or<param;or++){
//            for (int i=0;i<param.sizeDisk;i++){
//                pupil[i]=mat[or][param.disk2D[i][0]][param.disk2D[i][1]];
//            }
//            host_pupil= Pointer.to(pupil);
//            setHost2Device(device_pupil[or],host_pupil,param.sizeDisk,Sizeof.FLOAT);
//        }
//        
//    }
    
    
    
    
    public void updatePupil(float [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            pupil[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_pupil= Pointer.to(pupil);
        setHost2Device(device_pupil,host_pupil,param.sizeDisk,Sizeof.FLOAT);
        
        
    }
    
    
    
    public void updatePhase(Pointer device_p){
        
        //this.imshowDouble(device_p, "phaseDouble");
        
        //double -> float
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image_0 "+cudaResult);}
        MyVecDouble.double2float(custream,param.sizeDisk, device_phase, device_p);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image_1 "+cudaResult);}
        
        //this.imshowFloat(device_phase, "phaseFloat");
        
        //this.imshowFloat(device_phase, "phaseInit");
    }
    
    
//    public void updatePupil(Pointer device_p){
//        this.device_pupil=device_p;
//        
//    }
//    
//    
//    
//    public void updateKz(Pointer device_p){
//        this.device_kz=device_p;
//        //this.imshow(device_kz, "kz");
//        
//    }
    
    
    
    
    
    
    
    public void updateKx(float [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            kx[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_kx= Pointer.to(kx);
        setHost2Device(device_kx,host_kx,param.sizeDisk,Sizeof.FLOAT);
        
    }
    
    
    
    
    
    public void updateKy(float [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            ky[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_ky= Pointer.to(ky);
        setHost2Device(device_ky,host_ky,param.sizeDisk,Sizeof.FLOAT);
        
    }
    
    
    public void updateKz(float [][] mat){
        for (int i=0;i<param.sizeDisk;i++){
            kz[i]=mat[param.disk2D[i][0]][param.disk2D[i][1]];
        }
        host_kz= Pointer.to(kz);
        setHost2Device(device_kz,host_kz,param.sizeDisk,Sizeof.FLOAT);
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    public void free(){
        
        //IJ.log("free PSF called");
        
        JCuda.cudaFree(device_phase);
        
        JCuda.cudaFree(device_pupil);
        
        JCuda.cudaFree(device_position);
        
        JCuda.cudaFree(device_kx);
        JCuda.cudaFree(device_ky);
        JCuda.cudaFree(device_kz);
        JCuda.cudaFree(device_kz_oil);
        JCuda.cudaFree(device_kz_is_imaginary);
        JCuda.cudaFree(device_kz_oil_is_imaginary);
        JCuda.cudaFree(device_tmpFull);
//        JCuda.cudaFree(device_imagtmpFull);
        JCuda.cudaFree(device_fftdata);
        
        
        JCuda.cudaFree(device_sparseIndexOddDisk);
        JCuda.cudaFree(device_sparseIndexEvenDisk);
        
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutputNext);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput1);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput1);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput);
        
        
        JCuda.cudaFree(device_sparseIndexShift2DOutput);
        JCuda.cudaFree(device_gaussian);
        JCuda.cudaFree(device_gaussian2bleSize);
        
        JCuda.cudaFree(device_psfFFT);
        JCuda.cudaFree(device_psf);
        JCuda.cudaFree(device_psf_double);
        
        JCuda.cudaFree(device_ones);
        JCuda.cudaFree(device_alpha);
        JCuda.cudaFree(device_beta);
        
        JCuda.cudaFree(device_sumRes);
        
        
    JCuda.cudaFreeHost(host_sparseIndexShift2DOutput);
    
    JCuda.cudaFreeHost(host_sparseIndexOddShift2DOutput);
    JCuda.cudaFreeHost(host_sparseIndexEvenShift2DOutput);
    JCuda.cudaFreeHost(host_sparseIndexEvenShift2DOutputNext);
    JCuda.cudaFreeHost(host_sparseIndexOddShift2DOutputNext);
    
    JCuda.cudaFreeHost(host_sparseIndexEvenShift2DOutput1);
    JCuda.cudaFreeHost(host_sparseIndexOddShift2DOutput1);
    
    JCuda.cudaFreeHost(host_sparseIndexOddDisk);
    JCuda.cudaFreeHost(host_sparseIndexEvenDisk);
    
    JCuda.cudaFreeHost(host_kz_oil_is_imaginary);
    JCuda.cudaFreeHost(host_kz_is_imaginary);
    JCuda.cudaFreeHost(host_kz_oil);
    JCuda.cudaFreeHost(host_kz);
    JCuda.cudaFreeHost(host_ky);
    JCuda.cudaFreeHost(host_kx);
    JCuda.cudaFreeHost(host_pupil);
    JCuda.cudaFreeHost(host_phase);
    JCuda.cudaFreeHost(host_position);
        
        
        
        JCufft.cufftDestroy(plan2);
        plan2=null;
        
        JCufft.cufftDestroy(plan);
        plan=null;
        
        
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
    
    
    
    public void imshowDouble(Pointer deviceDiskShape,String title){
        double [] kk=new double [this.kz.length];
        Pointer hostkk = Pointer.to(kk);
        this.setDevice2Host( hostkk,deviceDiskShape, param.sizeDisk, Sizeof.DOUBLE);
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kk[i];
        }
        ImageShow.imshow(mat,title);
    }
    
    
    public void imshowFloat(Pointer deviceDiskShape,String title){
        float [] kk=new float [this.kz.length];
        Pointer hostkk = Pointer.to(kk);
        this.setDevice2Host( hostkk,deviceDiskShape, param.sizeDisk, Sizeof.FLOAT);
        float [][] mat = new float [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kk[i];
        }
        ImageShow.imshow(mat,title);
    }
    
    
    
    
    public void updateSigmaGaussianKernel(double sigma){
        this.setSigma(sigma);
    }
    
    
    
    
     
     
    public void computePSF(double [] x, double [] y, double [] zoil, double [] z){
        
        
        
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult+"  "+param.stream);}
        
        
        cudaResult=JCuda.cudaMemsetAsync(device_fftdata, 0, this.numberPSF*sizepow*2 * Sizeof.FLOAT, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem set 0");}
        for (int i=0;i<x.length;i++){
            this.position[i]=(float)x[i];
            this.position[i+numberPSF]=(float)y[i];
            this.position[i+2*numberPSF]=(float)(z[i]*param.weightZ);
            this.position[i+3*numberPSF]=(float)(zoil[i]);
            
        }
        
        
        
        setHost2Device(device_position,host_position,4*numberPSF,Sizeof.FLOAT);
        
        MyVecDouble.computePSF_phaseNManywithOil_f(custream, param.sizeDisk*numberPSF,param.sizeDisk,this.sizepow*2, device_kx, device_ky, device_kz,device_kz_is_imaginary, device_kz_oil, device_kz_oil_is_imaginary, device_pupil, device_phase, device_position,device_sparseIndexEvenDisk,device_sparseIndexOddDisk,device_fftdata,numberPSF);
        
        
        
        
        
        //imshow(this.sizepow*2*numberPSF, (int)Math.sqrt(this.sizepow)*2,device_fftdata,"fft"+this.sizepow,"FLOAT");
        
        cudaResult=JCufft.cufftExecC2C(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 5");}
         
        
        
        
        
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR FFTinvFirst PSFmany");}
        
        MyVecDouble.computePSF_signalN2Many_f(custream, numberPSF*sizeoutput2*sizeoutput2,sizeoutput2*sizeoutput2,sizepow*2, (float)Math.sqrt(sizepow),device_sparseIndexEvenShift2DOutput,device_sparseIndexOddShift2DOutput,device_fftdata,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext,device_psfFFT);
        
        
        ////FILTERING
        
        cudaResult=JCufft.cufftExecC2C(plan2, device_psfFFT,device_psfFFT,JCufft.CUFFT_FORWARD);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR FFTfor PSFmany");}
        
        
        
        MyVecDouble.mulMany_f(custream, sizeoutput2*sizeoutput2*2*this.numberPSF, sizeoutput2*sizeoutput2*2, device_psfFFT, device_psfFFT, device_gaussian2bleSize);
        
        
        
        cudaResult=JCufft.cufftExecC2C(plan2, device_psfFFT, device_psfFFT, JCufft.CUFFT_INVERSE);

        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR FFTinv PSFmany");}
        
        
        //MyVecDouble.computePSF_signalNsqrtMany_f(custream, sizeoutput1*sizeoutput1*this.numberPSF,sizeoutput1*sizeoutput1, device_psf, device_psfFFT, sizeoutput2*sizeoutput2,device_sparseIndexEvenShift2DOutput1,device_sparseIndexOddShift2DOutput1);
       
        MyVecDouble.computePSF_signalNsqrtMany_fcrop(custream, sizeoutput1*sizeoutput1*this.numberPSF,sizeoutput1*sizeoutput1,sizeoutput2*sizeoutput2, device_psf, device_psfFFT, sizeoutput2*sizeoutput2,device_sparseIndexEvenShift2DOutput1,device_sparseIndexOddShift2DOutput1);
       
        
        cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_T,sizeoutput1*sizeoutput1,numberPSF,device_alpha,device_psf,sizeoutput1*sizeoutput1,device_ones,1,device_beta,device_sumRes,1);
        
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cublasDgemv PSFmany");}
        //JCuda.cudaMemcpyAsync(hostSomOne, device_sumRes, numberPSF * Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));;if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda PSF computePSF()2 "+cudaResult);}
        
        //IJ.log("somResOne "+somResOne[0]+"  "+somResOne[1]);
        
        MyVecDouble.divScalarMany_f(custream,sizeoutput1*sizeoutput1*numberPSF,sizeoutput1*sizeoutput1, device_psf ,device_psf_double , device_psf,device_sumRes);
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult);}
        
        
        
    }
    
    
    
    
    
    
    
     /*
     
    public void computePSFwithoutOilRefIndexMismatch(double [] x, double [] y, double [] z){
        
        
        
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult+"  "+param.stream);}
        
        
        cudaResult=JCuda.cudaMemsetAsync(device_fftdata, 0, this.numberPSF*sizepow*2 * Sizeof.FLOAT, MyCudaStream.getCudaStream_t(param.stream));
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR mem set 0");}
        for (int i=0;i<x.length;i++){
            this.position[i]=(float)x[i];
            this.position[i+numberPSF]=(float)y[i];
            this.position[i+2*numberPSF]=(float)(z[i]*param.weightZ);
        }
        
        
        setHost2Device(device_position,host_position,3*numberPSF,Sizeof.FLOAT);
        
        MyVecDouble.computePSF_phaseNMany_f(custream, param.sizeDisk*numberPSF,param.sizeDisk,this.sizepow*2, device_kx, device_ky, device_kz, device_pupil, device_phase, device_position,device_sparseIndexEvenDisk,device_sparseIndexOddDisk,device_fftdata,numberPSF);
        
        
        //imshow(this.sizepow*2*numberPSF, (int)Math.sqrt(this.sizepow)*2,device_fftdata,"fft"+this.sizepow,"FLOAT");
        
        cudaResult=JCufft.cufftExecC2C(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 5");}
         
        
        
        
        
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR FFTinvFirst PSFmany");}
        
        MyVecDouble.computePSF_signalN2Many_f(custream, numberPSF*sizeoutput2*sizeoutput2,sizeoutput2*sizeoutput2,sizepow*2, (float)Math.sqrt(sizepow),device_sparseIndexEvenShift2DOutput,device_sparseIndexOddShift2DOutput,device_fftdata,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext,device_psfFFT);
        
        
        ////FILTERING
        
        cudaResult=JCufft.cufftExecC2C(plan2, device_psfFFT,device_psfFFT,JCufft.CUFFT_FORWARD);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR FFTfor PSFmany");}
        
        
        
        MyVecDouble.mulMany_f(custream, sizeoutput2*sizeoutput2*2*this.numberPSF, sizeoutput2*sizeoutput2*2, device_psfFFT, device_psfFFT, device_gaussian2bleSize);
        
        
        
        cudaResult=JCufft.cufftExecC2C(plan2, device_psfFFT, device_psfFFT, JCufft.CUFFT_INVERSE);

        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR FFTinv PSFmany");}
        
        
        MyVecDouble.computePSF_signalNsqrtMany_f(custream, sizeoutput2*sizeoutput2*this.numberPSF,sizeoutput2*sizeoutput2, device_psf, device_psfFFT, sizeoutput2*sizeoutput2,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext);
       
        
        cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_T,sizeoutput2*sizeoutput2,numberPSF,device_alpha,device_psf,sizeoutput2*sizeoutput2,device_ones,1,device_beta,device_sumRes,1);
        
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cublasDgemv PSFmany");}
        //JCuda.cudaMemcpyAsync(hostSomOne, device_sumRes, numberPSF * Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));;if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda PSF computePSF()2 "+cudaResult);}
        
        //IJ.log("somResOne "+somResOne[0]+"  "+somResOne[1]);
        
        MyVecDouble.divScalarMany_f(custream,sizeoutput2*sizeoutput2*numberPSF,sizeoutput2*sizeoutput2, device_psf ,device_psf_double , device_psf,device_sumRes);
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult);}
        
        
        
    }
    
    */
    
    
    
    
    
    
    
    
    public float [][] getPSF(int number){
        
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF getPSF([][]) "+cudaResult);}
        Pointer hostoutput=Pointer.to(res);//output point to res
        //IJ.log("res "+res.length+"   ");
        this.setDevice2Host(hostoutput, device_psf.withByteOffset(number*sizeoutput1*sizeoutput1*Sizeof.FLOAT),sizeoutput1*sizeoutput1,Sizeof.FLOAT);
        
        //2D conversion:
        
        res2D = new float [sizeoutput1][sizeoutput1];//already done in setsizeoutput
        for (int i=0;i<sizeoutput1;i++){
            for (int ii=0;ii<sizeoutput1;ii++){
                res2D[i][ii]=res[i*sizeoutput1+ii];
            }
        }
        //long tttttt=System.currentTimeMillis();
        
        return res2D;
    }
    
    
    
    
    public float [] getAllPSF(){
        
        float [] res = new float[numberPSF*sizeoutput1*sizeoutput1];
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF getPSF([][]) "+cudaResult);}
        Pointer hostoutput=Pointer.to(res);//output point to res
        //IJ.log("res "+res.length+"   ");
        this.setDevice2Host(hostoutput, device_psf,numberPSF*sizeoutput1*sizeoutput1,Sizeof.FLOAT);
        
        //2D conversion:
        
        
        return res;
    }
    
    /*public double [] getPSF1D(int number){
        Pointer hostoutput=Pointer.to(res);//output point to res
        //IJ.log("res "+res.length+"   ");
        this.setDevice2Host(hostoutput, device_psf.withByteOffset(number*sizeoutput2*sizeoutput2*Sizeof.FLOAT),sizeoutput2*sizeoutput2,Sizeof.FLOAT);
        
        //2D conversion:
        
        res1D = new double [sizeoutput2*sizeoutput2];
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                res1D[i*sizeoutput2+ii]=res[i*sizeoutput2+ii];
            }
        }
        //long tttttt=System.currentTimeMillis();
        
        return res1D;
    }*/
    
    public Pointer getPointerPSF(int number){
        return device_psf_double.withByteOffset(number*sizeoutput2*sizeoutput2*Sizeof.DOUBLE);
    }
    
    public Pointer getPointerPSF(){
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult);}
        
        return device_psf_double;
    }
    
    
    public Pointer getPointerPSFfloat(){
        
        
        return device_psf;
    }
    
    
    
    
    public void save(String path){
        
        float [][][] im = new float [5][param.size][param.size];
        
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