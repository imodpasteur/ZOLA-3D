/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.utils.ImageShow;
import ij.IJ;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcufft.JCufft;
import jcuda.jcufft.cufftHandle;
import jcuda.jcufft.cufftType;
import jcuda.jcublas.JCublas2;
import static  jcuda.jcublas.JCublas2.cublasCreate;
import jcuda.jcublas.cublasHandle;
import jcuda.jcublas.cublasOperation;
import static jcuda.jcusparse.JCusparse.cusparseSgthr;
import static jcuda.jcusparse.JCusparse.cusparseSsctr;
import static jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import jcuda.runtime.cudaError;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import jcuda.jcublas.JCublas;
import jcuda.jcublas.JCublas2;

import static jcuda.jcublas.JCublas2.cublasCreate;
import static jcuda.jcublas.JCublas2.cublasDestroy;
import static jcuda.jcublas.JCublas2.cublasDgemm;
import static jcuda.jcublas.JCublas2.cublasGetVector;
import static jcuda.jcublas.JCublas2.cublasSetPointerMode;
import static jcuda.jcublas.JCublas2.cublasSetVector;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_N;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_C;
import static jcuda.jcublas.cublasPointerMode.CUBLAS_POINTER_MODE_DEVICE;
import jcuda.runtime.cudaMemcpyKind;
/**
 *
 * @author benoit
 */
public class Predetection_deconvolutionTest_ {
    
    double RSB=3;
    double photonThreshold=0;
    float weightBackground=(float)1.05;
    int cudaResult;
    cufftHandle plan;
    
    cufftHandle planMany;
    
    cublasHandle  handlecublas;
    
    
    
    
    
    Pointer device_fullImage;
    Pointer device_background;
    
    Pointer device_correlImage;
    
    
    Pointer device_sparseIndexPaddingShift2DEven;
    
    Pointer device_o;
    Pointer device_convFFT;
    
    Pointer device_fullImageFFT;
    
    Pointer device_dicFFT;
    Pointer device_dicFFT_mirrored;
    
    
    Pointer device_tmpFFT;
    
    Pointer  device_beta;
    Pointer  device_alpha;
    Pointer  device_a;
    
    
    
    
    //new pointers
    
    Pointer device_M;
    Pointer device_tmp;
    Pointer device_tmp2;
    Pointer device_PSF;
    Pointer device_PSF_mirrored;
    
    
    
    
    Pointer host_resultCorrelation;
    float [] resultCorrelation_f;
    float [][][] resultCorrelation;
    
    
    
    float [] fullImage;
    float [] fullO;
    float [] fullB;
    
    int sizePSF;
    
    Pointer host_fullImage;
    
    //double background_average;
    
    int width;
    int height;
    int sizeFullImage;
    int sizeFullImageX;
    int sizeFullImageY;
    int totalSize;
    int totalSizeUnpadded;
    DataPhase_ dparam;
    
    public double [] range;
    
    double [][] psf;
    
    double [][] resvect ;
    
    double [][] matrix2_detection;
    double [][][] matrix3_detection ;
    
    int streamId;
    
    //double thresholdDeconvolution=.2;
    
    double [] tmp_b0;
    double [][] tmp_b1;
    
    double zstep;
    double xystep;
    
    
    public Predetection_deconvolutionTest_(int sizeFullImage, DataPhase_ dparam,double mini, double maxi, double step,double photonThreshold){
        this.photonThreshold=photonThreshold;
        double initStep=step;
        
        if (sizeFullImage%2==1){
            sizeFullImage++;
        }
        
        
        this.sizeFullImage=sizeFullImage;
        
        this.sizeFullImageX=sizeFullImage;
        this.sizeFullImageY=sizeFullImage;
        
        this.streamId=dparam.getStream();
        
        int sizeFullImageFFT=sizeFullImage*2;
        totalSize=sizeFullImageFFT*sizeFullImageFFT;
        totalSizeUnpadded=sizeFullImage*sizeFullImage;
        
        
        
        long [] memtot=new long[1];
        long [] memfree=new long[1];
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        
        
        int number=(int)Math.ceil((maxi-mini)/step);
        range = new double [number];
        double z=mini;
        
        for (int i=0;i<range.length;i++){
            range[i]=z;
            z+=step;
        }
        
        
        long nb=totalSizeUnpadded+sizeFullImage*sizeFullImage+totalSize*2 ;
        
        long nb2=(totalSize*2*range.length)*2+sizeFullImage*sizeFullImage*range.length;
        
        
        long neededMemory=((long)3*(nb+nb2) * (long)Sizeof.FLOAT);//*2 memory mo free
        
        while (neededMemory>memfree[0]){
            step*=1.1;
            number=(int)Math.ceil((maxi-mini)/step);
            range = new double [number];
            z=mini;
            for (int i=0;i<range.length;i++){
                range[i]=z;
                z+=step;
            }
            nb2=(totalSize*2*range.length)*2+sizeFullImage*sizeFullImage*range.length;
            neededMemory=((long)2*(nb+nb2) * (long)Sizeof.FLOAT);//*2 memory mo free
            
        }
        
        
        
        if (step>initStep){
            IJ.log("Due to available GPU memory, z detection is performed with a precision of "+step+" µm");
        }
        
        zstep=step;
        xystep=dparam.param.xystep;
        
        
        
        
        fullImage=new float [totalSizeUnpadded];
        
        fullO=new float [totalSizeUnpadded*range.length];
        fullB=new float [totalSizeUnpadded];
        host_fullImage=Pointer.to(fullImage);
        
        
        matrix3_detection = new double [range.length][sizeFullImage][sizeFullImage];
        matrix2_detection = new double [sizeFullImage][sizeFullImage];
        
        int [] sparseIndexPaddingShift2DEven = new int[sizeFullImage*sizeFullImage];
        int [] sparseIndexPadding = new int[sizeFullImage*sizeFullImage];
        int [] sparseIndexOdd = new int[sizeFullImageFFT*sizeFullImageFFT];
        int [] sparseIndexEven = new int[sizeFullImageFFT*sizeFullImageFFT];
        int [] sparseIndexShift2D = new int[sizeFullImageFFT*sizeFullImageFFT];
        for (int i=0;i<sizeFullImage;i++){
            for (int ii=0;ii<sizeFullImage;ii++){
                sparseIndexPadding[i*sizeFullImage+ii]=(i+sizeFullImage/2)*sizeFullImageFFT+(ii+sizeFullImage/2);
                
                int x=((i+sizeFullImage/2));
                int y=((ii+sizeFullImage/2));
                int p=(((x+sizeFullImageFFT/2)%sizeFullImageFFT)*sizeFullImageFFT+((y+sizeFullImageFFT/2)%sizeFullImageFFT));
                sparseIndexPaddingShift2DEven[i*sizeFullImage+ii]=p*2;
            }
        }
        
        
        for (int i=0;i<sizeFullImageFFT;i++){
            for (int ii=0;ii<sizeFullImageFFT;ii++){
                sparseIndexShift2D[i*sizeFullImageFFT+ii]=(((i+sizeFullImageFFT/2)%sizeFullImageFFT)*sizeFullImageFFT+((ii+sizeFullImageFFT/2)%sizeFullImageFFT));
                
                
           
            }
        }
        
        
        for (int i=0;i<sizeFullImageFFT*sizeFullImageFFT;i++){
            sparseIndexOdd[i]=(i*2)+1;
            sparseIndexEven[i]=(i*2);
        }
        
        
        
        sizePSF=dparam.param.sizeoutput;
        
       
        
        device_sparseIndexPaddingShift2DEven=new Pointer();
        cudaResult =cudaMalloc(device_sparseIndexPaddingShift2DEven, totalSizeUnpadded * Sizeof.INT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 6 "+cudaResult);}
        cudaResult =cudaMemcpyAsync(device_sparseIndexPaddingShift2DEven, Pointer.to(sparseIndexPaddingShift2DEven), totalSizeUnpadded*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
        
        
        
        
        
        
        
        
        
        
        
        
        
        device_background = new Pointer();
        cudaResult =cudaMalloc(device_background, sizeFullImage*sizeFullImage * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        resultCorrelation_f=new float[sizeFullImage*sizeFullImage*range.length];
        resultCorrelation=null;
        host_resultCorrelation = Pointer.to(resultCorrelation_f);
        
        device_correlImage = new Pointer();
        cudaResult =cudaMalloc(device_correlImage, sizeFullImage*sizeFullImage*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        device_fullImage = new Pointer();
        cudaResult =cudaMalloc(device_fullImage, sizeFullImage*sizeFullImage * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_o = new Pointer();
        cudaResult =cudaMalloc(device_o, sizeFullImage*sizeFullImage*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_M = new Pointer();
        cudaResult =cudaMalloc(device_M, sizeFullImage*sizeFullImage * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_tmp = new Pointer();
        cudaResult =cudaMalloc(device_tmp, sizeFullImage*sizeFullImage * sizePSF*sizePSF * range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        device_tmp2 = new Pointer();
        cudaResult =cudaMalloc(device_tmp2, sizeFullImage*sizeFullImage * range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_PSF = new Pointer();
        cudaResult =cudaMalloc(device_PSF, sizePSF*sizePSF*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        device_PSF_mirrored = new Pointer();
        cudaResult =cudaMalloc(device_PSF_mirrored, sizePSF*sizePSF*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        
        
        
        plan = new cufftHandle();
        
        cudaResult =JCufft.cufftPlan2d(plan, sizeFullImageFFT,sizeFullImageFFT, cufftType.CUFFT_C2C);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 11 "+cudaResult);}
        
        cudaResult =JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}
        
        
        handlecublas=MyCudaStream.getHandleCublas(streamId);
        
        planMany = new cufftHandle();
        
        int [] size= new int[2];
        size[0]=sizeFullImageFFT;
        size[1]=sizeFullImageFFT;
        
        //JCufft.cufftPlan2d(plan, param.size,param.size, cufftType.CUFFT_Z2Z);
        JCufft.cufftPlanMany(planMany, 2,size, size,1,sizeFullImageFFT*sizeFullImageFFT, size,1,sizeFullImageFFT*sizeFullImageFFT,cufftType.CUFFT_C2C,range.length);
        
        cudaResult =JCufft.cufftSetStream(planMany, MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}
        
        
        this.dparam=dparam;
        
        
        
        
        //test();
       
        
        
        
        device_fullImageFFT = new Pointer();
        cudaResult =cudaMalloc(device_fullImageFFT, totalSize*2 * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 181 "+cudaResult);}
        
        
        
        
        
        device_dicFFT = new Pointer();
        cudaResult =cudaMalloc(device_dicFFT, totalSize*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        device_dicFFT_mirrored = new Pointer();
        cudaResult =cudaMalloc(device_dicFFT_mirrored, totalSize*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        
        device_convFFT = new Pointer();
        cudaResult =cudaMalloc(device_convFFT, totalSize*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        
        device_tmpFFT = new Pointer();
        cudaResult =cudaMalloc(device_tmpFFT, totalSize*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 184 "+cudaResult);}
        
        
        
        double [] one = new double [1];
        one[0]=1;
        Pointer hostOne=Pointer.to(one);
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.FLOAT);
        float [] alpha = new float[1];
        alpha[0]=1;
        cudaMemcpyAsync(device_alpha, Pointer.to(alpha), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(streamId));
        
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.FLOAT);
        float [] beta = new float[1];
        beta[0]=0;
        cudaMemcpyAsync(device_beta, Pointer.to(beta), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(streamId));
        
        int sizeDeviceA=(int)Math.max(range.length,sizePSF*sizePSF)*totalSizeUnpadded;
        device_a = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_a, sizeDeviceA * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda build 7");return ;
        }
        float [] a=new float [sizeDeviceA];
        for (int i=0;i<sizeDeviceA;i++){
            a[i]=1;
        }
        Pointer host_a=Pointer.to(a);
        JCuda.cudaMemcpyAsync(device_a, host_a, sizeDeviceA* Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(streamId));
        
        
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        //IJ.log("warning no dictionnary");
        
         makeDictionnary();
        
         
         
         
    }
    
    
    public int getSizeImageCrossCorrel(){
        return this.sizeFullImage;
    }
    
    
        
        
    
    
    
    private void makeDictionnary(){
        
        int sizeFullImageFFT=sizeFullImage*2;
        int [] sparseIndexPaddingShift2DEvenDic = new int[sizeFullImage*sizeFullImage];
        for (int i=0;i<dparam.param.sizeoutput;i++){
            for (int ii=0;ii<dparam.param.sizeoutput;ii++){
                
                int x=((i+sizeFullImageFFT/2-dparam.param.sizeoutput/2));
                int y=((ii+sizeFullImageFFT/2-dparam.param.sizeoutput/2));
                int p=(((x+sizeFullImageFFT/2)%sizeFullImageFFT)*sizeFullImageFFT+((y+sizeFullImageFFT/2)%sizeFullImageFFT));
                sparseIndexPaddingShift2DEvenDic[i*dparam.param.sizeoutput+ii]=p*2;
                
            }
        }
        
        Pointer device_sparseIndexPaddingShift2DEvenDic;
        device_sparseIndexPaddingShift2DEvenDic=new Pointer();
        cudaResult =cudaMalloc(device_sparseIndexPaddingShift2DEvenDic, dparam.param.sizeoutput*dparam.param.sizeoutput * Sizeof.INT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 6 "+cudaResult);}
        cudaResult =cudaMemcpyAsync(device_sparseIndexPaddingShift2DEvenDic, Pointer.to(sparseIndexPaddingShift2DEvenDic), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
        
        
        Pointer device_dic;
        device_dic = new Pointer();
        cudaResult =cudaMalloc(device_dic, dparam.param.sizeoutput*dparam.param.sizeoutput*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 182 "+cudaResult);}
        
        
        cudaResult =JCuda.cudaMemsetAsync(device_dicFFT, 0, this.totalSize*2*range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));
        cudaResult =JCuda.cudaMemsetAsync(device_dicFFT_mirrored, 0, this.totalSize*2*range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));
        //double [][][] psftmp = new double [range.length][dparam.param.sizeoutput][dparam.param.sizeoutput];
        double [] psf_d = new double [dparam.param.sizeoutput*dparam.param.sizeoutput];
        float [] psf_f = new float [dparam.param.sizeoutput*dparam.param.sizeoutput];
        float [] psf_f_mirrored = new float [dparam.param.sizeoutput*dparam.param.sizeoutput];
        psf= new double [range.length][dparam.param.sizeoutput*dparam.param.sizeoutput];
        for (int i=0;i<range.length;i++){
        //for (int i=0;i<1;i++){
            //IJ.log("range PSFmaker "+range[i]+"  "+dparam.param.Zfocus);
            
            
            dparam.psf.computePSF(0,0,dparam.param.Zfocus,range[i]);
            //if (i==0){ImageShow.imshow(dparam.psf.getPSF(),"psf");}
                    
            //conversion float
            psf_d=dparam.psf.getPSF1D();
            double mean=0;
            for (int u=0;u<psf_d.length;u++){
                psf[i][u]=psf_d[u];//*(int)Math.sqrt(this.totalSize);
                psf_f[u]=(float)psf_d[u];//*(int)Math.sqrt(this.totalSize);//multiply by padded width to avoid to do it after FFT
            }
            
            for (int u=0;u<dparam.param.sizeoutput;u++){
                for (int uu=0;uu<dparam.param.sizeoutput;uu++){
                    psf_f_mirrored[u*dparam.param.sizeoutput+uu]=psf_f[(dparam.param.sizeoutput-u-1)*dparam.param.sizeoutput+(dparam.param.sizeoutput-uu-1)];
                    //psf_f_mirrored[u*dparam.param.sizeoutput+uu]=psf_f[u*dparam.param.sizeoutput+(uu)];
                }
            }
            
            
            //dparam.psf.imshow(dparam.param.sizeoutput, psf_d, "DOUBLE");
            //if (i==0){dparam.psf.imshow(dparam.param.sizeoutput, psf_f, "FLOAT");}
            
            //if (i==0){dparam.psf.imshow(dparam.param.sizeoutput, psf_d, "DOUBLE");}
        
            
            //dparamFullImage.psf.imshow(this.totalSizeUnpadded,(int)Math.sqrt(this.totalSizeUnpadded), dparamFullImage.psf.getPointerPSF(), "ImagePSFgen","DOUBLE");
            
            cudaResult =cudaMemcpyAsync(device_dic, Pointer.to(psf_f), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            
            cudaResult =cudaMemcpyAsync(device_PSF.withByteOffset(i*dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT), Pointer.to(psf_f), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            
            
            
            //dparam.psf.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput,dparam.param.sizeoutput, device_dic, "dic","FLOAT");
            
            //dparam.psf.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput,dparam.param.sizeoutput, device_sparseIndexPaddingShift2DEvenDic, "sparse","INT");
        
            cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), dparam.param.sizeoutput*dparam.param.sizeoutput, device_dic, device_sparseIndexPaddingShift2DEvenDic,device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
            
            
            
            cudaResult =JCufft.cufftExecC2C(plan, device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 28 "+cudaResult);}

            //if (i==0){dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "dic1_"+range[i],"FLOAT");}
        
            
            //if (i==0){dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "fftImage2","FLOAT");}
        
            
            //IJ.log("make dic ok");
            
            
            
            
            
            cudaResult =cudaMemcpyAsync(device_dic, Pointer.to(psf_f_mirrored), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            cudaResult =cudaMemcpyAsync(device_PSF_mirrored.withByteOffset(i*dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT), Pointer.to(psf_f_mirrored), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), dparam.param.sizeoutput*dparam.param.sizeoutput, device_dic, device_sparseIndexPaddingShift2DEvenDic,device_dicFFT_mirrored.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
            cudaResult =JCufft.cufftExecC2C(plan, device_dicFFT_mirrored.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),device_dicFFT_mirrored.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 28 "+cudaResult);}
            
            //if (i==0){dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT_mirrored.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "fftImage2","FLOAT");}
        
        
        }
        
        JCuda.cudaFree(device_dic) ;
        JCuda.cudaFree(device_sparseIndexPaddingShift2DEvenDic) ;
        MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(streamId), this.totalSize*2*this.range.length, device_dicFFT, device_dicFFT, (float)(1./Math.sqrt(totalSize)));//here we multiply with sqrt(total size) to avoid to do it after
        //MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(streamId), this.totalSize*2*this.range.length, device_dicFFT_mirrored, device_dicFFT_mirrored, (float)(Math.sqrt(totalSize)));
        
    }
        
    
    
    public double [][] getPSFNonNormalized(){
        double [][] psf = new double [this.psf.length][this.psf[0].length];
        for (int i=0;i<psf.length;i++){
            for (int ii=0;ii<this.psf[0].length;ii++){
                psf[i][ii]=this.psf[i][ii];
            }
        }
        return psf;
    }
    
    
    
    public double [] getRange(){
        double [] range = new double [this.range.length];
        for (int i=0;i<this.range.length;i++){
            
                range[i]=this.range[i];
            
        }
        return range;
    }
    
    
    
    public double getStep(){
        return zstep;
    }
    
    
    public boolean setOinit(){
        
        double [][][] o = new double [range.length][this.sizeFullImage][this.sizeFullImage];
        IJ.log("o length "+o.length);
        o[0][3][3]=1;
        o[1][3][3]=0;
        o[0][3][2]=0;
        try{
            o[0][12][12]=1;
        o[1][12][12]=0;
        o[0][12][11]=0;
        }catch(Exception rthurir){}
        for (int i=0;i<o.length;i++){
            
            for (int ii=0;ii<o[i].length;ii++){
                for (int iii=0;iii<o[i][ii].length;iii++){
                    float v=(float)o[i][ii][iii];
                    fullO[i*o[i].length*o[i][ii].length+ii*o[i][ii].length+iii]=v+(float)0;
                }
            }
            
        }
        Pointer host_o=Pointer.to(fullO);
        cudaResult =cudaMemcpyAsync(device_o,host_o, totalSizeUnpadded*range.length*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        //dparam.psf.imshow(totalSizeUnpadded*range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "input","FLOAT");//show first image
        
        
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
//    public boolean setOinit(double [][][] o){
//        
//        for (int i=0;i<o.length;i++){
//            for (int ii=0;ii<o[i].length;ii++){
//                for (int iii=0;iii<o[i][ii].length;iii++){
//                    fullO[i*o[i].length*o[ii].length+ii*o[ii].length+iii]=(float)o[i][ii][iii];
//                }
//            }
//        }
//        Pointer host_o=Pointer.to(fullO);
//        cudaResult =cudaMemcpyAsync(device_o,host_o, totalSizeUnpadded*range.length*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
//        
//        
//        
//        return (cudaResult == cudaError.cudaSuccess);
//    }
    
    
    
    
    
    public boolean setBackgroundinit(float meanValue){
        
        double [][] b = new double [this.sizeFullImage][this.sizeFullImage];
        
        
        for (int i=0;i<sizeFullImage;i++){
            for (int ii=0;ii<sizeFullImage;ii++){
                fullB[i*sizeFullImage+ii]=meanValue;
            }
        }
        Pointer host_b=Pointer.to(fullB);
        cudaResult =cudaMemcpyAsync(device_background,host_b, totalSizeUnpadded*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        //dparam.psf.imshow(totalSizeUnpadded,sizeFullImage, device_o, "dic","FLOAT");//show first image
        
        
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    public boolean setBackgroundinit(double [][] b){
        
        for (int i=0;i<b.length;i++){
            for (int ii=0;ii<b[i].length;ii++){
                float v=(float)b[i][ii];
                fullB[i*b[i].length+ii]=v;
            }
        }
        Pointer host_b=Pointer.to(fullB);
        cudaResult =cudaMemcpyAsync(device_background,host_b, totalSizeUnpadded*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        //dparam.psf.imshow(totalSizeUnpadded,sizeFullImage, device_o, "dic","FLOAT");//show first image
        
        
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    public boolean setBackgroundMorphologyCPU(double [][] image){
        
        long t0 = System.currentTimeMillis();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        if (tmp_b0==null){
            tmp_b0=new double [sizePSF*sizePSF];
        }
        if (tmp_b1==null){
            tmp_b1=new double [sizeFullImage][sizeFullImage];
        }
        
        int x,y;
        for (int i=0;i<image.length;i++){
            for (int ii=0;ii<image[i].length;ii++){
                int nb=0;
                for (int j=-this.sizePSF/2;j<this.sizePSF/2;j++){
                    x=i+j;
                    if ((x>=0)&&(x<image.length)){
                        for (int jj=-this.sizePSF/2;jj<this.sizePSF/2;jj++){
                            y=ii+jj;
                            if ((y>=0)&&(y<image.length)){
                                tmp_b0[nb]=image[x][y];
                                nb++;
                            }
                        }
                        
                    }
                }
                Arrays.sort(tmp_b0,0,nb);
                tmp_b1[i][ii]=tmp_b0[nb/10];//first decile
                
            }
            
        }
        //ImageShow.imshow(tmp_b1,"tmp_b1");
        
        for (int i=0;i<image.length;i++){
            for (int ii=0;ii<image[i].length;ii++){
                int nb=0;
                for (int j=-this.sizePSF/2;j<this.sizePSF/2;j++){
                    x=i+j;
                    if ((x>=0)&&(x<image.length)){
                        for (int jj=-this.sizePSF/2;jj<this.sizePSF/2;jj++){
                            y=ii+jj;
                            if ((y>=0)&&(y<image.length)){
                                tmp_b0[nb]=tmp_b1[x][y];
                                nb++;
                            }
                        }
                        
                    }
                }
                Arrays.sort(tmp_b0,0,nb);
                fullB[i*sizeFullImage+ii]=(float)tmp_b0[9*nb/10];//last decile
                
            }
            
        }
        
        Pointer host_b=Pointer.to(fullB);
        cudaResult =cudaMemcpyAsync(device_background,host_b, totalSizeUnpadded*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        long t1 = System.currentTimeMillis();
        IJ.log("bckg CPU "+(t1-t0));
        dparam.psf.imshow(totalSizeUnpadded,sizeFullImage, device_background, "bckgCPU","FLOAT");//show first image
        
        
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    
    public boolean setBackgroundMorphologyGPU(double [][] image){
        
        long t0 = System.currentTimeMillis();
JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
                
        int smoothsize=(int)Math.sqrt(sizePSF);
        if (smoothsize%2==0){
            smoothsize+=1;
        }
        if (smoothsize<3){
            smoothsize=3;
        }
        //remove SLIGHTLY noise
        MyVecDouble.localMean(MyCudaStream.getCUstream(streamId), this.totalSizeUnpadded, this.sizeFullImageX, this.sizeFullImageY, smoothsize, device_tmp, device_fullImage);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMinimum(MyCudaStream.getCUstream(streamId), this.totalSizeUnpadded, this.sizeFullImageX, this.sizeFullImageY, sizePSF, device_tmp2, device_tmp);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMean(MyCudaStream.getCUstream(streamId), this.totalSizeUnpadded, this.sizeFullImageX, this.sizeFullImageY, smoothsize, device_tmp, device_tmp2);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMaximum(MyCudaStream.getCUstream(streamId), this.totalSizeUnpadded, this.sizeFullImageX, this.sizeFullImageY, sizePSF, device_tmp2, device_tmp);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMean(MyCudaStream.getCUstream(streamId), this.totalSizeUnpadded, this.sizeFullImageX, this.sizeFullImageY, smoothsize, device_background, device_tmp2);//here we multiply with sqrt(total size) to avoid to do it after
        

        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}

        
        long t1 = System.currentTimeMillis();
        IJ.log("bckg GPU "+(t1-t0));
        dparam.psf.imshow(totalSizeUnpadded,sizeFullImage, device_fullImage, "device_image","FLOAT");//show first image
        dparam.psf.imshow(totalSizeUnpadded,sizeFullImage, device_tmp2, "device_tmp2","FLOAT");//show first image
        dparam.psf.imshow(totalSizeUnpadded,sizeFullImage, device_background, "device_bckgGPU","FLOAT");//show first image
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    
    public boolean setImage(double [][] image){
        
        
        this.width=image.length;
        this.height=image[0].length;
        
//        double mean=0;
//        for (int i=0,j=0;i<width;i++){
//            for (int ii=0;ii<height;ii++){
//                mean+=(float)image[i][ii];
//            }
//        }
//        mean/=width*height;
//        double std=Math.sqrt(mean);
//        float count=0;
//        float meanBckg=0;
//        for (int i=0,j=0;i<width;i++){
//            for (int ii=0;ii<height;ii++){
//                if ((image[i][ii]>mean-3*std)&&(image[i][ii]<mean+3*std)){
//                     meanBckg+=(float)image[i][ii];
//                     count++;
//                }
//                
//               
//            }
//        }
//        meanBckg/=count;
        //no need to set it at 0
        cudaResult=JCuda.cudaMemsetAsync(device_fullImageFFT, 0, this.totalSize*2 * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
       
        if (Math.max(width, height)+(Math.max(width, height)%2)!=sizeFullImage){
            IJ.log("OOPS: problems are coming: image size should be set for FFT first "+width+"  "+height+"  "+sizeFullImage);
        }
        for (int i=0,j=0;i<sizeFullImage;i++){
            for (int ii=0;ii<sizeFullImage;ii++){
                if ((i<width)&&(ii<height)){
                    
                    fullImage[j++]=(float)image[i][ii];
                }
                else{
                    fullImage[j++]=0;//meanBckg;
                }
            }
        }
        
        
        
        cudaResult=cudaMemcpyAsync(device_fullImage, host_fullImage, totalSizeUnpadded*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult+"  "+totalSizeUnpadded);}
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        
        setBackgroundMorphologyCPU(image);
        
        setBackgroundMorphologyGPU(image);
        
        //setBackgroundinit(meanBckg);
        
       
       cudaResult=cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), totalSizeUnpadded, device_fullImage, this.device_sparseIndexPaddingShift2DEven,device_fullImageFFT, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 3");}
        
       cudaResult= JCufft.cufftExecC2C(plan, device_fullImageFFT,device_fullImageFFT,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}
       
       
       //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(this.totalSizeUnpadded), device_fullImage, "Image","FLOAT");
        
       
       
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    
    
    
    public double [][] deconvolution_RL(int iter_max){
        if (this.resultCorrelation==null){
            resultCorrelation=new float[range.length][sizeFullImage][sizeFullImage];
        }
        return deconvolution_RL(iter_max,resultCorrelation);
    }
    
    
    
    
    public double [][] deconvolution_RL(int iter_max,float [][][] resultCorrelation){
        //iter_max=1;//30;
        setOinit();IJ.log("WARNING: o init manually");
        //IJ.log("WARNING: iter deconv fixed value");
        //long t0=System.currentTimeMillis();
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        float power=(float).125;
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
            
        //dparam.psf.imshow(totalSizeUnpadded*this.range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "oooooINIT ","FLOAT");
        iter_max=1;
        if (true)
                return null;
        MyVecDouble.mulScalarFloat(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, device_background, device_background, (float)weightBackground);
        //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_background, "device_background mul","FLOAT");
        toto: for (int t=0;t<iter_max;t++){
            
            //IJ.log("iter "+(iter_max-t)+" ");
            for (int y=0;y<1;y++){//range.length;y++){
            //for (int y=range.length-1;y>=0;y--){//range.length;y++){
                
                cudaResult=JCuda.cudaMemsetAsync(device_tmp2, 0, sizeFullImageX*sizeFullImageY*range.length* Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
                //cudaResult=JCuda.cudaMemsetAsync(device_tmp, 0, sizeFullImageX*sizeFullImageY* sizePSF* sizePSF* Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
                dparam.psf.imshow(sizeFullImage*sizeFullImage,sizeFullImage, device_o, "device_o","FLOAT");
                dparam.psf.imshow(sizePSF*sizePSF,sizePSF, device_PSF, "device_PSF","FLOAT");
                //System.out.println("tutu");
                cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}

                long t0=System.currentTimeMillis();
                for (int yuyu=0;yuyu<100;yuyu++){  
                cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
                
                MyVecDouble.manualFilteringStackedFast(MyCudaStream.getCUstream(streamId), sizeFullImageX*sizeFullImageY*sizePSF*sizePSF*range.length,sizeFullImageX,sizeFullImageY,sizePSF,range.length,device_tmp, device_o,device_PSF);
                
                
                cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_T,sizePSF*sizePSF,sizeFullImageX*sizeFullImageY*range.length,device_alpha,device_tmp,sizePSF*sizePSF,device_a,1,device_beta,device_tmp2,1);
                
                cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_N,totalSizeUnpadded,range.length,device_alpha,device_tmp2,sizeFullImageX*sizeFullImageY,device_a,1,device_beta,device_tmp,1);
                
                //divide I by (o*psf+b):
                MyVecDouble.addanddivide(MyCudaStream.getCUstream(streamId), sizeFullImageX*sizeFullImageY, device_M, device_fullImage , device_tmp,device_background);

                }
                
                long t1=System.currentTimeMillis();
                
               IJ.log("TIME A "+(t1-t0));
                
                
                //dparam.psf.imshow(sizePSF*sizePSF,sizePSF, device_PSF, "device_PSF","FLOAT");
                //dparam.psf.imshow(sizeFullImageX*sizeFullImageY*sizePSF*sizePSF,sizePSF*sizePSF, device_tmp, "device_tmp","FLOAT");
                //dparam.psf.imshow(sizeFullImageX*sizeFullImageY,sizeFullImageY, device_M, "device_M","FLOAT");
                
//                cudaResult=JCuda.cudaMemsetAsync(device_convFFT, 0, this.totalSize*2*this.range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
//
//                //dparam.psf.imshow(sizeFullImage*sizeFullImage,sizeFullImage, device_sparseIndexPaddingShift2DEven, "device_sparseIndexPaddingShift2DEven","INT");
//
//
//                //dparam.psf.imshow(sizeFullImage*sizeFullImage,sizeFullImage, device_o, "device_oInit","FLOAT");
//                
//                IJ.log("SIZE "+totalSize*2*3+"  "+(int)Math.sqrt(this.totalSize)*2);
//
//                IJ.log("loop can be condensed by implementing cuda operation");
//                for (int i=0;i<range.length;i++){
//                    cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), sizeFullImage*sizeFullImage, device_o.withByteOffset(i*this.totalSizeUnpadded*Sizeof.FLOAT), device_sparseIndexPaddingShift2DEven,device_convFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
//                }
//                
//                if (y==0){dparam.psf.imshow(totalSize*2*5,(int)Math.sqrt(this.totalSize)*2, device_convFFT, "CUSPARSE device_convFFT","FLOAT");}
//                
                
                
                
                cudaResult=JCuda.cudaMemsetAsync(device_convFFT, 0, this.totalSize*2*this.range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
                
                //dparam.psf.imshow(totalSize*2*5,(int)Math.sqrt(this.totalSize)*2, device_convFFT, "RESET device_convFFT","FLOAT");
                
                long t2=System.currentTimeMillis();
                
            for (int yuyu=0;yuyu<100;yuyu++){    
                cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
                
                
                //dparam.psf.imshow(sizeFullImage*sizeFullImage*range.length,(int)this.sizeFullImage, device_o, "device_o","FLOAT");
                
                MyVecDouble.mycusparsemoduloSsctr(MyCudaStream.getCUstream(streamId), sizeFullImage*sizeFullImage*range.length,totalSize*2,sizeFullImage*sizeFullImage,device_convFFT, device_o, device_sparseIndexPaddingShift2DEven);
                
                
                //
                
                
                

                cudaResult= JCufft.cufftExecC2C(planMany, device_convFFT,device_convFFT,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}

                //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_convFFT, "device_convFFT","FLOAT");

                //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT, "device_dicFFT","FLOAT");


                MyVecDouble.complexeMulKernel(MyCudaStream.getCUstream(streamId), totalSize*range.length, totalSize, device_tmpFFT, device_convFFT, device_dicFFT);

                
                JCufft.cufftExecC2C(planMany, device_tmpFFT,device_tmpFFT,JCufft.CUFFT_INVERSE);

                
                MyVecDouble.makeResultCorrelation(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, totalSizeUnpadded,totalSize*2, device_correlImage, device_tmpFFT, this.device_sparseIndexPaddingShift2DEven);


                MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, device_correlImage, device_correlImage, (float)Math.sqrt(totalSize));

                
                
                

                //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "device_correlImage","FLOAT");

                //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage.withByteOffset(5*this.totalSizeUnpadded*Sizeof.FLOAT), "device_correlImage5","FLOAT");



                //sum psf*o
                //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
                cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_N,totalSizeUnpadded,range.length,device_alpha,device_correlImage,totalSizeUnpadded,device_a,1,device_beta,device_tmpFFT,1);
                //if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda deconv "+cudaResult+"   "+streamId);}
                //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}

                //NOW, IMAGES ARE NOT STACKED

                //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_tmpFFT, "p*o "+y,"FLOAT");
                //if (y==0){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_tmpFFT, "p*o "+y,"FLOAT");}//#################### GOOD TO SHOW

                //divide I by (o*psf+b):
                MyVecDouble.addanddivide(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, device_convFFT, device_fullImage , device_tmpFFT,device_background);

                cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
            }
                dparam.psf.imshow(sizeFullImageX*sizeFullImageY,sizeFullImageY, device_convFFT, "device_convFFT","FLOAT");
                
                long t3=System.currentTimeMillis();
                
                IJ.log("TIME B "+(t3-t2));

                //FFT computation
                cudaResult=JCuda.cudaMemsetAsync(device_tmpFFT, 0, this.totalSize*2 * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
                
                
                //if (y==9){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), "o_mid "+y,"FLOAT");}
                //if (y==0){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), "o ("+y,"FLOAT");}

                //here, device_convFFT corresponds to i/(p*o+b)
               //this was written to test FFT : cudaResult=cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), totalSizeUnpadded, device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), this.device_sparseIndexPaddingShift2DEven,device_tmpFFT, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 3");}
               cudaResult=cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), totalSizeUnpadded, device_convFFT, this.device_sparseIndexPaddingShift2DEven,device_tmpFFT, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 3");}

               cudaResult= JCufft.cufftExecC2C(plan, device_tmpFFT,device_tmpFFT,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}

               //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_tmpFFT, "device_tmpFFT 0 "+y,"FLOAT");
               
               //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_dicFFT.withByteOffset(y*this.totalSize*2*Sizeof.FLOAT), "mulPSF 0 "+y,"FLOAT");
               
               MyVecDouble.complexeMulKernel(MyCudaStream.getCUstream(streamId), totalSize, totalSize,  device_convFFT,device_tmpFFT, device_dicFFT_mirrored.withByteOffset(y*this.totalSize*2*Sizeof.FLOAT));

               //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_convFFT, "device_convFFT 0 "+y,"FLOAT");
               
                JCufft.cufftExecC2C(plan, device_convFFT,device_convFFT,JCufft.CUFFT_INVERSE);

                //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_convFFT, "device_convFFT 1 "+y,"FLOAT");
                
                MyVecDouble.makeResultCorrelation(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, totalSizeUnpadded,totalSize*2, device_correlImage, device_convFFT, this.device_sparseIndexPaddingShift2DEven);
                
                //here, device_correlImage contains (i/(p*o+b)*p*)
                
                //if (y==9){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "(i/(p*o+b)*p*)  "+t,"FLOAT");}
                
                //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "device_correlImage  "+y,"FLOAT");
                
                //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "device_correlImage "+y,"FLOAT");

                MyVecDouble.mul_fl_pow(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), device_correlImage,power);

                
                
            }
            
            //dparam.psf.imshow(totalSizeUnpadded*this.range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "ooooo "+t,"FLOAT");
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
            power*=2;
            /*if (t>2){
                power*=2;
            }
            if (t>4){
                power*=2;
            }
            if (t>8){
                power*=2;
            }*/
            power=Math.min(power, 4);
            
        }
        
        
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        //long t1=System.currentTimeMillis();
        //IJ.log("time = "+(t1-t0));
        
        //dparam.psf.imshow(totalSizeUnpadded*range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "input","FLOAT");//show first image
        
        
        
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
          
        cudaResult =cudaMemcpyAsync(this.host_resultCorrelation, device_o, range.length*totalSizeUnpadded*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
        
        for (int i=0;i<this.sizeFullImage;i++){
            for (int ii=0;ii<this.sizeFullImage;ii++){
                for (int z=0;z<range.length;z++){
                    resultCorrelation[z][i][ii]=(float)((double)this.resultCorrelation_f[z*sizeFullImage*sizeFullImage+i*sizeFullImage+ii]);
                }
                
            }   
        } 
        //ImageShow.imshow(resultCorrelation,"resultCorrelation");
         
        
        
        
        
        double [][] vect= getVectMaxPosition(resultCorrelation);
        
        Arrays.sort(vect, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o2[0]).compareTo(o1[0]);
            }
        });
        
        
        
        
        
        /*ImageShow.imshow(resultCorrelation, "o");
        
        double [][] tmp = new double[width][height];
        for (int i=0;i<this.sizeFullImage;i++){
            for (int ii=0;ii<this.sizeFullImage;ii++){
                for (int z=0;z<range.length;z++){
                    resultCorrelation[z][i][ii]=0;
                }
                tmp[i][ii]=0;//this.thresholdDeconvolution;//background value
            }   
        } 
        
        double [][][] psfShow = new double [range.length][this.sizePSF][this.sizePSF];
        for (int z=0;z<range.length;z++){
                for (int i=0;i<sizePSF;i++){
                    for (int ii=0;ii<sizePSF;ii++){
                        psfShow[z][i][ii]=this.psf[z][i*sizePSF+ii];
                    }
                }
            }
        ImageShow.imshow(psfShow,"psf");
        for (int i=0;i<vect.length;i++){
            resultCorrelation[(int)vect[i][1]][(int)vect[i][2]][(int)vect[i][3]]=(float)(vect[i][0]);
            for (int j=-this.sizePSF/2,k=0;j<sizePSF/2;j++){
                for (int jj=-this.sizePSF/2;jj<sizePSF/2;jj++,k++){
                    if (((int)vect[i][3]+jj>=0)&&((int)vect[i][3]+jj<height)&&((int)vect[i][2]+j>=0)&&((int)vect[i][2]+j<width)){
                        double xxx=this.psf[(int)vect[i][1]][k];
                        
                        tmp[(int)vect[i][2]+j][(int)vect[i][3]+jj]+=(float)(vect[i][0]*this.psf[(int)vect[i][1]][k]);
                    }
                }
            }
        }
        
        
        ImageShow.imshow(tmp,"model");
        
        ImageShow.imshow(resultCorrelation,"resultCorrelation END");
        
        IJ.log("len "+vect.length);
        */
        
        
        return vect;
        
    }
    
    
    
    
    public void runCrossCorrelationOnGPU(){
        
        cudaResult=JCuda.cudaMemsetAsync(device_convFFT, 0, totalSizeUnpadded*range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
                
        
        MyVecDouble.complexeConjugateKernel(MyCudaStream.getCUstream(streamId), totalSize*range.length, totalSize, device_tmpFFT, device_fullImageFFT, device_dicFFT);
        
        
        JCufft.cufftExecC2C(planMany, device_tmpFFT,device_tmpFFT,JCufft.CUFFT_INVERSE);
        
        float [] tmp = new float[1];
        
        //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_background, "device_background ","FLOAT");
        //MyVecDouble.mulScalarFloat(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, device_background, device_background, (float).001);

        //dparam.psf.imshow(totalSize*2*range.length,(int)Math.sqrt(this.totalSize), device_dicFFT, "dic2","FLOAT");
        
        //MyVecDouble.makeResultCorrelation(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, totalSizeUnpadded,totalSize*2, device_correlImage, device_tmpFFT, this.device_sparseIndexPaddingShift2DEven);
        //use the following to use cross correl as input of deconvolution
        MyVecDouble.makeResultCorrelationNormalized(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, totalSizeUnpadded,totalSize*2, device_o, device_tmpFFT, this.device_sparseIndexPaddingShift2DEven,(float)0.999*(float)((float)Math.sqrt(totalSizeUnpadded)/256.)*this.sizePSF*this.sizePSF,device_background,(float)0.1);
        
        //dparam.psf.imshow(totalSizeUnpadded*range.length,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "cross correl init","FLOAT");//show first image
        //dparam.psf.imshow(totalSizeUnpadded*range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "cross_correl","FLOAT");//show first image
        
        //compute average 
        //cudaResult=JCublas2.cublasDasum(handlecublas,totalSizeUnpadded*range.length,device_tmpFFT,1,device_o);
//        cudaResult=JCublas2.cublasSasum(handlecublas,totalSizeUnpadded*range.length,device_o,1,device_convFFT);
//        
//        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
//        cudaResult =cudaMemcpyAsync( Pointer.to(tmp),device_convFFT, 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
//        IJ.log("sum="+tmp[0]+"   "+(tmp[0]/(totalSizeUnpadded*range.length))+"   "+totalSizeUnpadded+"  "+range.length);
//        //here, device_tmpFFT[0] contains the sum
//        tmp[0]=0;
//        cudaResult =cudaMemcpyAsync( device_convFFT,Pointer.to(tmp), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
//        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        //MyVecDouble.subtractMeanWithSumAsInputWithPositiveConstraint(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, device_o, device_o, device_convFFT,(float)0.001);
        //MyVecDouble.sub(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, device_o, device_o, device_background);
        
    }
    
    
    
    
    
    
        
        double [][] getVectMaxPosition(float [][][] resultConvolution){
            
            
            
            int decalFilter=2;
            ArrayList<double []> al = new ArrayList<double []>();
            
            int sizeEdge2remove=decalFilter+1;//sizePSF/2
            
            
            for (int r=0;r<resultConvolution.length;r++){
                for (int i=Math.max(decalFilter,sizeEdge2remove);i<width-Math.max(decalFilter,sizeEdge2remove);i++){
                    for (int j=Math.max(decalFilter,sizeEdge2remove);j<height-Math.max(decalFilter,sizeEdge2remove);j++){
                        matrix3_detection[r][i][j]=Double.NEGATIVE_INFINITY;
                        if ((resultConvolution[r][i][j]>0)){
                            
                            for (int a=-decalFilter;a<=decalFilter;a++){
                                if ((r+a>0)&&(r+a<resultConvolution.length)){
                                    matrix3_detection[r][i][j]=Math.max(resultConvolution[r+a][i][j],matrix3_detection[r][i][j]);
                                }
                            }
                        }
                    }
                }
            }
            
            
            for (int r=0;r<resultConvolution.length;r++){
                for (int i=Math.max(decalFilter,sizeEdge2remove);i<width-Math.max(decalFilter,sizeEdge2remove);i++){
                    for (int j=Math.max(decalFilter,sizeEdge2remove);j<height-Math.max(decalFilter,sizeEdge2remove);j++){
                        matrix2_detection[i][j]=Double.NEGATIVE_INFINITY;
                        if ((resultConvolution[r][i][j]>0)){
                            for (int a=-decalFilter;a<=decalFilter;a++){
                                matrix2_detection[i][j]=Math.max(matrix3_detection[r][i][j+a],matrix2_detection[i][j]);
                            }
                        }
                    }
                }
                double max=0;
                //double min=0;
                for (int j=Math.max(decalFilter,sizeEdge2remove);j<height-Math.max(decalFilter,sizeEdge2remove);j++){
                    for (int i=Math.max(decalFilter,sizeEdge2remove);i<width-Math.max(decalFilter,sizeEdge2remove);i++){
                        if ((resultConvolution[r][i][j]>this.fullB[i*this.sizeFullImage+j]/100)){//check if local maximum higher than avg background/2 (to be faster at the end)
                            
                            if (resultConvolution[r][i][j]==matrix2_detection[i][j]){
                                max=Double.NEGATIVE_INFINITY;
                                for (int a=-decalFilter;a<=decalFilter;a++){
                                    max=Math.max(matrix2_detection[i+a][j],max);
                                }
                                if (max==resultConvolution[r][i][j]){
                                    
//                                    min=Double.POSITIVE_INFINITY;
//                                    for (int a=-decalFilter;a<=decalFilter;a++){
//                                        for (int aa=-decalFilter;aa<=decalFilter;aa++){
//                                            min=Math.min(resultConvolution[r][i+a][j+aa],min);
//                                        }
//                                    }
//                                    if((min>0)){
                                        double [] machin = new double[8];
                                        
                                        
                                        double topZ=0;
                                        if ((r>1)&&(r<resultConvolution.length-1)){
                                            topZ=topOfParabola(-1.*zstep,resultConvolution[r-1][i][j],0,resultConvolution[r][i][j],1.*zstep,resultConvolution[r+1][i][j]);
                                        }
                                        
                                        double topX=0;
                                        if ((i>1)&&(i<resultConvolution[0].length-1)){
                                            topX=topOfParabola(-1.*xystep,resultConvolution[r][i-1][j],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i+1][j]);
                                        }
                                        
                                        double topY=0;
                                        if ((j>1)&&(j<resultConvolution[0][0].length-1)){
                                            topY=topOfParabola(-1.*xystep,resultConvolution[r][i][j-1],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i][j+1]);
                                        }
                                        
                                        
                                        machin[0]=resultConvolution[r][i][j];
                                        machin[1]=r;
                                        machin[2]=i;
                                        machin[3]=j;
                                        machin[4]=topZ;
                                        machin[5]=topX;
                                        machin[6]=topY;
                                        machin[7]=this.fullB[i*this.sizeFullImage+j];
                                        al.add(machin);
//                                    }
                                }
                                
                                
                            }
                        }
                    }
                }
            }
            
            //compute photon sum instead of max value so that photon number is well estimated
            for (int i=al.size()-1;i>=0;i--){
                double sum=0;
                for (int j=-decalFilter;j<decalFilter;j++){
                    for (int jj=-decalFilter;jj<decalFilter;jj++){
                        for (int jjj=-decalFilter;jjj<decalFilter;jjj++){
                            if ((al.get(i)[1]+j>=0)&&(al.get(i)[1]+j<resultConvolution.length)&&(al.get(i)[2]+jj>=0)&&(al.get(i)[2]+jj<width)&&(al.get(i)[3]+jjj>=0)&&(al.get(i)[3]+jjj<height)){
                                sum+=resultConvolution[(int)al.get(i)[1]+j][(int)al.get(i)[2]+jj][(int)al.get(i)[3]+jjj];
                            }
                        }
                    }
                }
                IJ.log("al.get(i)[7] "+al.get(i)[7]+"  "+sum+"  "+i);
                if ((sum<al.get(i)[7]*RSB)||(sum<photonThreshold)){
                    al.remove(i);
                }
                else{
                    al.get(i)[0]=sum;
                }
            }
            
            
            
            double [][] vect = new double [al.size()][8];
            for (int i=0;i<al.size();i++){
                vect[i]=al.get(i);
            }
            
            al.clear();
            
            
            
            
            return vect;
        }
        
        
        
        
        double [][] getVectMaxPositionNewOld(float [][][] resultConvolution){
            double thresholdDeconvolution=500;
            int decalFilter=2;
            ArrayList<double []> al = new ArrayList<double []>();
            if (Math.random()<.001){
                IJ.log("speed could be slightly improved in Predetection_deconvolution(getVectMaxPosition)");
                //remove tmp
            }
            
            double [][] tmp = new double [width][height];
            double [][][] tmp0 = new double [resultConvolution.length][width][height];
            
            for (int r=0;r<resultConvolution.length;r++){
                for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                    for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                        if ((resultConvolution[r][i][j]>0)){
                            tmp0[r][i][j]=Double.NEGATIVE_INFINITY;
                            for (int a=-decalFilter;a<=decalFilter;a++){
                                if ((r+a>0)&&(r+a<resultConvolution.length)){
                                    tmp0[r][i][j]=Math.max(resultConvolution[r+a][i][j],tmp0[r][i][j]);
                                }
                            }
                        }
                    }
                }
            }
            
            
            for (int r=0;r<resultConvolution.length;r++){
                for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                    for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                        if ((resultConvolution[r][i][j]>0)){
                            tmp[i][j]=Double.NEGATIVE_INFINITY;
                            for (int a=-decalFilter;a<=decalFilter;a++){
                                tmp[i][j]=Math.max(tmp0[r][i][j+a],tmp[i][j]);
                            }
                        }
                    }
                }
                double max=0;
                //double min=0;
                for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                    for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                        if ((resultConvolution[r][i][j]>thresholdDeconvolution)){
                            
                            if (resultConvolution[r][i][j]==tmp[i][j]){
                                max=Double.NEGATIVE_INFINITY;
                                for (int a=-decalFilter;a<=decalFilter;a++){
                                    max=Math.max(tmp[i+a][j],max);
                                }
                                if (max==resultConvolution[r][i][j]){
                                    
//                                    min=Double.POSITIVE_INFINITY;
//                                    for (int a=-decalFilter;a<=decalFilter;a++){
//                                        for (int aa=-decalFilter;aa<=decalFilter;aa++){
//                                            min=Math.min(resultConvolution[r][i+a][j+aa],min);
//                                        }
//                                    }
//                                    if((min>0)){
                                        double [] machin = new double[7];
                                        
                                        
                                        double topZ=0;
                                        if ((r>1)&&(r<resultConvolution.length-1)){
                                            topZ=topOfParabola(-1.*zstep,resultConvolution[r-1][i][j],0,resultConvolution[r][i][j],1.*zstep,resultConvolution[r+1][i][j]);
                                        }
                                        
                                        double topX=0;
                                        if ((i>1)&&(i<resultConvolution[0].length-1)){
                                            topX=topOfParabola(-1.*xystep,resultConvolution[r][i-1][j],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i+1][j]);
                                        }
                                        
                                        double topY=0;
                                        if ((j>1)&&(j<resultConvolution[0][0].length-1)){
                                            topY=topOfParabola(-1.*xystep,resultConvolution[r][i][j-1],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i][j+1]);
                                        }
                                        
                                        
                                        machin[0]=resultConvolution[r][i][j];
                                        machin[1]=r;
                                        machin[2]=i;
                                        machin[3]=j;
                                        machin[4]=topZ;
                                        machin[5]=topX;
                                        machin[6]=topY;
                                        al.add(machin);
//                                    }
                                }
                                
                                
                            }
                        }
                    }
                }
            }
            
            double [][] vect = new double [al.size()][7];
            for (int i=0;i<al.size();i++){
                vect[i]=al.get(i);
            }
            
            al.clear();
            
            //compute photon sum instead of max value so that photon number is well estimated
            for (int i=0;i<vect.length;i++){
                double sum=0;
                for (int j=-decalFilter;j<decalFilter;j++){
                    for (int jj=-decalFilter;jj<decalFilter;jj++){
                        for (int jjj=-decalFilter;jjj<decalFilter;jjj++){
                            if ((vect[i][1]+j>=0)&&(vect[i][1]+j<resultConvolution.length)&&(vect[i][2]+jj>=0)&&(vect[i][2]+jj<width)&&(vect[i][3]+jjj>=0)&&(vect[i][3]+jjj<height)){
                                sum+=resultConvolution[(int)vect[i][1]+j][(int)vect[i][2]+jj][(int)vect[i][3]+jjj];
                            }
                        }
                    }
                }
                vect[i][0]=sum;
            }
            
            
            return vect;
        }

        
        double [][] getVectMaxPositionOld(float [][][] resultConvolution){
        double thresholdDeconvolution=500;

            int decalFilter=2;
            ArrayList<double []> al = new ArrayList<double []>();
            double [][] tmp = new double [width][height];
            for (int r=0;r<resultConvolution.length;r++){
                for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                    for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                        if ((resultConvolution[r][i][j]>0)){
                            tmp[i][j]=Double.NEGATIVE_INFINITY;
                            for (int a=-decalFilter;a<=decalFilter;a++){
                                tmp[i][j]=Math.max(resultConvolution[r][i][j+a],tmp[i][j]);
                            }
                        }
                    }
                }
                double max=0;
                //double min=0;
                for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                    for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                        if ((resultConvolution[r][i][j]>thresholdDeconvolution)){
                            
                            if (resultConvolution[r][i][j]==tmp[i][j]){
                                max=Double.NEGATIVE_INFINITY;
                                for (int a=-decalFilter;a<=decalFilter;a++){
                                    max=Math.max(tmp[i+a][j],max);
                                }
                                if (max==resultConvolution[r][i][j]){
                                    
//                                    min=Double.POSITIVE_INFINITY;
//                                    for (int a=-decalFilter;a<=decalFilter;a++){
//                                        for (int aa=-decalFilter;aa<=decalFilter;aa++){
//                                            min=Math.min(resultConvolution[r][i+a][j+aa],min);
//                                        }
//                                    }
//                                    if((min>0)){
                                        double [] machin = new double[8];
                                        
                                        
                                        double topZ=0;
                                        if ((r>1)&&(r<resultConvolution.length-1)){
                                            topZ=topOfParabola(-1.*zstep,resultConvolution[r-1][i][j],0,resultConvolution[r][i][j],1.*zstep,resultConvolution[r+1][i][j]);
                                        }
                                        
                                        double topX=0;
                                        if ((i>1)&&(i<resultConvolution[0].length-1)){
                                            topX=topOfParabola(-1.*xystep,resultConvolution[r][i-1][j],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i+1][j]);
                                        }
                                        
                                        double topY=0;
                                        if ((j>1)&&(j<resultConvolution[0][0].length-1)){
                                            topY=topOfParabola(-1.*xystep,resultConvolution[r][i][j-1],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i][j+1]);
                                        }
                                        
                                        
                                        machin[0]=resultConvolution[r][i][j];//photon number
                                        machin[1]=r;//z
                                        machin[2]=i;//x
                                        machin[3]=j;//y
                                        machin[4]=topZ;//z subpixelic
                                        machin[5]=topX;//x subpixelic
                                        machin[6]=topY;//y subpixelic
                                        
                                        al.add(machin);
//                                    }
                                }
                                
                                
                            }
                        }
                    }
                }
            }
            
            double [][] vect = new double [al.size()][7];
            for (int i=0;i<al.size();i++){
                vect[i]=al.get(i);
            }
            
            al.clear();
            return vect;
        }


    
    
        private double topOfParabola(double xa,double ya,double xb,double yb,double xc,double yc){
            double a=(yc-ya)/((xc-xa)*(xc-xb))-(yb-ya)/((xb-xa)*(xc-xb));
            double b=((yb-ya)/(xb-xa))-a*(xb+xa);
            return(-b/(2*a));
        }
        
        
        
        double [][] getVectMaxPositionOldOld(float [][][] resultConvolution){
            double thresholdDeconvolution=500;


            int decalFilter=3;
            ArrayList<double []> al = new ArrayList<double []>();

            for (int r=0;r<resultConvolution.length;r++){
                toto:for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                    for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                        if ((resultConvolution[r][i][j]>thresholdDeconvolution)){
                            boolean ok=true;
                            search:for (int a=-decalFilter;a<=decalFilter;a++){
                            for (int aa=-decalFilter;aa<=decalFilter;aa++){
                                    if ((a!=0)||(aa!=0)){
                                        if (resultConvolution[r][i][j]<=resultConvolution[r][(a+i)][(aa+j)]){
                                            ok=false;
                                            break search;
                                        }
                                    }
                                }
                            }

                            if (ok){
                                double [] machin = new double[4];
                                machin[0]=resultConvolution[r][i][j];
                                machin[1]=r;
                                machin[2]=i;
                                machin[3]=j;
                                al.add(machin);
                            }
                        }
                    }
                }
            }
            
            double [][] vect = new double [al.size()][4];
            for (int i=0;i<al.size();i++){
                vect[i]=al.get(i);
            }
            al.clear();
            return vect;
        }


    
    
        
        
        
        
    
    
    
    
    
    
    public void free(){
        
        
        JCuda.cudaFree(device_tmp) ;
        JCuda.cudaFree(device_tmp2) ;
        JCuda.cudaFree(device_M) ;
        JCuda.cudaFree(device_PSF) ;
        JCuda.cudaFree(device_PSF_mirrored) ;
        
        JCuda.cudaFree(device_background) ;
        
        JCuda.cudaFree(device_o) ;

        JCuda.cudaFree(device_fullImage) ;


        JCuda.cudaFree(device_a) ;
        JCuda.cudaFree(device_alpha) ;
        JCuda.cudaFree(device_beta) ;

        JCuda.cudaFree(device_correlImage) ;

        JCuda.cudaFree(device_sparseIndexPaddingShift2DEven) ;


        JCuda.cudaFree(device_fullImageFFT) ;

        JCuda.cudaFree(device_dicFFT) ;


        JCuda.cudaFree(device_tmpFFT) ;

        JCufft.cufftDestroy(plan);
        plan=null; 
        JCufft.cufftDestroy(planMany);
        planMany=null; 

        JCuda.cudaFree(this.host_fullImage) ;
        JCuda.cudaFree(this.host_resultCorrelation) ;
    }
    
}
