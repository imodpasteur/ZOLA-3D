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
public class Predetection_deconvolutionFast_ {
    
    ArrayList<Integer> freeFramePosition = new ArrayList<Integer>();
    
    
    int cudaResult;
    cufftHandle plan;
    
    cufftHandle planMany;
    cufftHandle planDepthMany;
    cublasHandle  handlecublas;
    
    
    int nbFrame=1;
    int nb_iter=50;
    
    Pointer device_image;
    Pointer device_background;
    
    Pointer device_tmp1;
    Pointer device_tmp2;
    
    
    Pointer device_sparseIndexPaddingShift2DEven;
    
    
    
    
    Pointer device_o;
    
    Pointer device_imageFFT;
    
    Pointer device_dicFFT;
    Pointer device_dicFFT_mirrored;
    
    
    Pointer device_tmpFFT1;
    Pointer device_tmpFFT2;
    
    
    Pointer  device_beta;
    Pointer  device_alpha;
    Pointer  device_a;
    
    
    
    
    
    
    
    
    float [][] image;
    float [][] o;
    float [][] bckg;
    
    
    Pointer [] host_image;
    
    Pointer [] host_bckg;
    
    int sizePSF;
    
    
    
    int width;
    int height;
    int depth;
    int wh;
    
    
    int widthFFT;
    int heightFFT;
    int whFFT;
    
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
    
    
    public Predetection_deconvolutionFast_(int width,int height, DataPhase_ dparam,double mini, double maxi, double step ){
        
        
        //step*=5;IJ.log("WARNING WARNING: deconv : step mult by 5");
        
        
        for (int i=0;i<nbFrame;i++){
            freeFramePosition.add(i);
        }
        
        double initStep=step;
        
        
        
        this.streamId=dparam.getStream();
        
        this.width=width;
        this.height=height;
        wh=width*height;
        
        widthFFT=width*2;
        heightFFT=height*2;
        whFFT=widthFFT*heightFFT;
        
        
        
        long [] memtot=new long[1];
        long [] memfree=new long[1];
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        
        step=.5;IJ.log("WARNING, step fixed to .5");
        
        int number=(int)Math.ceil((maxi-mini)/step);
        
        
        
        range = new double [number];
        double z=mini;
        
        for (int i=0;i<range.length;i++){
            range[i]=z;
            z+=step;
        }
        
        
        
        depth=range.length;
        
        
        
        if (step>initStep){
            IJ.log("Due to available GPU memory, z detection is performed with a precision of "+step+" µm");
        }
        
        zstep=step;
        xystep=dparam.param.xystep;
        
        
        
        
        image=new float [nbFrame][wh];
        host_image = new Pointer[nbFrame];
        for (int i=0;i<nbFrame;i++){
            host_image[i]=Pointer.to(image[i]);
        }
        
        
        o=new float [nbFrame][wh*depth];
        bckg=new float [nbFrame][wh];
        
        host_bckg = new Pointer[nbFrame];
        for (int i=0;i<nbFrame;i++){
            host_bckg[i]=Pointer.to(image[i]);
        }
        
        
        int [] sparseIndexPaddingShift2DEven = new int[wh];
        int [] sparseIndexPadding = new int[wh];
        int [] sparseIndexOdd = new int[whFFT];
        int [] sparseIndexEven = new int[whFFT];
        int [] sparseIndexShift2D = new int[whFFT];
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                sparseIndexPadding[i*height+ii]=(i+width/2)*heightFFT+(ii+height/2);
                
                int x=((i+width/2));
                int y=((ii+height/2));
                int p=(((x+widthFFT/2)%widthFFT)*heightFFT+((y+heightFFT/2)%heightFFT));
                sparseIndexPaddingShift2DEven[i*height+ii]=p*2;
            }
        }
        
        
        for (int i=0;i<widthFFT;i++){
            for (int ii=0;ii<heightFFT;ii++){
                sparseIndexShift2D[i*heightFFT+ii]=(((i+widthFFT/2)%widthFFT)*heightFFT+((ii+heightFFT/2)%heightFFT));
                
                
           
            }
        }
        
        
        for (int i=0;i<whFFT;i++){
            sparseIndexOdd[i]=(i*2)+1;
            sparseIndexEven[i]=(i*2);
        }
        
        
        
        sizePSF=dparam.param.sizeoutput;
        
       
        
        device_sparseIndexPaddingShift2DEven=new Pointer();
        cudaResult =cudaMalloc(device_sparseIndexPaddingShift2DEven, wh * Sizeof.INT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 6 "+cudaResult);}
        cudaResult =cudaMemcpyAsync(device_sparseIndexPaddingShift2DEven, Pointer.to(sparseIndexPaddingShift2DEven), wh*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
        
        
        
        
        
        
        
        
        
        
        
        
        
        device_background = new Pointer();
        cudaResult =cudaMalloc(device_background, wh * nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        
        
        device_tmp1 = new Pointer();
        cudaResult =cudaMalloc(device_tmp1, wh*depth*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        device_tmp2 = new Pointer();
        cudaResult =cudaMalloc(device_tmp2, wh*depth*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_tmpFFT1 = new Pointer();
        cudaResult =cudaMalloc(device_tmpFFT1, whFFT*2*depth*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        device_tmpFFT2 = new Pointer();
        cudaResult =cudaMalloc(device_tmpFFT2, whFFT*2*depth*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_image = new Pointer();
        cudaResult =cudaMalloc(device_image, wh*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        device_o = new Pointer();
        cudaResult =cudaMalloc(device_o, wh*range.length*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        
        
        
        
        plan = new cufftHandle();
        
        cudaResult =JCufft.cufftPlan2d(plan, widthFFT,heightFFT, cufftType.CUFFT_C2C);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 11 "+cudaResult);}
        
        cudaResult =JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}
        
        
        handlecublas=MyCudaStream.getHandleCublas(streamId);
        
        planDepthMany = new cufftHandle();
        
        int [] size= new int[2];
        size[0]=widthFFT;
        size[1]=heightFFT;
        
        //JCufft.cufftPlan2d(plan, param.size,param.size, cufftType.CUFFT_Z2Z);
        JCufft.cufftPlanMany(planDepthMany, 2,size, size,1,widthFFT*heightFFT, size,1,whFFT,cufftType.CUFFT_C2C,depth*nbFrame);
        
        cudaResult =JCufft.cufftSetStream(planDepthMany, MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}
        
        
        planMany = new cufftHandle();
        
        
        //JCufft.cufftPlan2d(plan, param.size,param.size, cufftType.CUFFT_Z2Z);
        JCufft.cufftPlanMany(planMany, 2,size, size,1,widthFFT*heightFFT, size,1,whFFT,cufftType.CUFFT_C2C,nbFrame);
        
        cudaResult =JCufft.cufftSetStream(planMany, MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}
        
        
        this.dparam=dparam;
        
        
        
        
        //test();
       
        
        
        
        device_imageFFT = new Pointer();
        cudaResult =cudaMalloc(device_imageFFT, whFFT*2*nbFrame * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 181 "+cudaResult);}
        
        
        
        
        device_dicFFT = new Pointer();
        cudaResult =cudaMalloc(device_dicFFT, whFFT*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        device_dicFFT_mirrored = new Pointer();
        cudaResult =cudaMalloc(device_dicFFT_mirrored, whFFT*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 183 "+cudaResult);}
        
        
        
        
        
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
        
        
        device_a = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_a, depth*wh*nbFrame * Sizeof.FLOAT);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda build 7");return ;
        }
        float [] a=new float [depth*wh*nbFrame];
        for (int i=0;i<depth*wh*nbFrame;i++){
            a[i]=1;
        }
        Pointer host_a=Pointer.to(a);
        JCuda.cudaMemcpyAsync(device_a, host_a, depth*wh*nbFrame* Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(streamId));
        
        
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        //IJ.log("warning no dictionnary");
        
         makeDictionnary();
        
         
         
         
    }
    
    
    
        
        
    
    
    
    private void makeDictionnary(){
        
        int [] sparseIndexPaddingShift2DEvenDic = new int[wh];
        for (int i=0;i<dparam.param.sizeoutput;i++){
            for (int ii=0;ii<dparam.param.sizeoutput;ii++){
                
                int x=((i+widthFFT/2-dparam.param.sizeoutput/2));
                int y=((ii+heightFFT/2-dparam.param.sizeoutput/2));
                int p=(((x+widthFFT/2)%widthFFT)*heightFFT+((y+heightFFT/2)%heightFFT));
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
        
        
        cudaResult =JCuda.cudaMemsetAsync(device_dicFFT, 0, whFFT*2*depth * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));
        cudaResult =JCuda.cudaMemsetAsync(device_dicFFT_mirrored, 0, whFFT*2*depth * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));
        //double [][][] psftmp = new double [range.length][dparam.param.sizeoutput][dparam.param.sizeoutput];
        double [] psf_d = new double [dparam.param.sizeoutput*dparam.param.sizeoutput];
        float [] psf_f = new float [dparam.param.sizeoutput*dparam.param.sizeoutput];
        float [] psf_f_mirrored = new float [dparam.param.sizeoutput*dparam.param.sizeoutput];
        psf= new double [range.length][dparam.param.sizeoutput*dparam.param.sizeoutput];
        
        double [] psftmp= new double [range.length*dparam.param.sizeoutput*dparam.param.sizeoutput];
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
                    psftmp[i*dparam.param.sizeoutput*dparam.param.sizeoutput+u*dparam.param.sizeoutput+uu]=psf_f[u*dparam.param.sizeoutput+uu];
                    psf_f_mirrored[u*dparam.param.sizeoutput+uu]=psf_f[(dparam.param.sizeoutput-u-1)*dparam.param.sizeoutput+(dparam.param.sizeoutput-uu-1)];
                    //psf_f_mirrored[u*dparam.param.sizeoutput+uu]=psf_f[u*dparam.param.sizeoutput+(uu)];
                }
            }
            
            
            //if (i==0){dparam.psf.imshow(dparam.param.sizeoutput, psf_f, "FLOAT");}
            
            //if (i==0){dparam.psf.imshow(dparam.param.sizeoutput, psf_d, "DOUBLE");}
        
            
            //dparamFullImage.psf.imshow(this.totalSizeUnpadded,(int)Math.sqrt(this.totalSizeUnpadded), dparamFullImage.psf.getPointerPSF(), "ImagePSFgen","DOUBLE");
            
            cudaResult =cudaMemcpyAsync(device_dic, Pointer.to(psf_f), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            
            
            
            //dparam.psf.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput,dparam.param.sizeoutput, device_dic, "dic","FLOAT");
            
            //dparam.psf.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput,dparam.param.sizeoutput, device_sparseIndexPaddingShift2DEvenDic, "sparse","INT");
        
            cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), dparam.param.sizeoutput*dparam.param.sizeoutput, device_dic, device_sparseIndexPaddingShift2DEvenDic,device_dicFFT.withByteOffset(i*whFFT*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
            
            
            
            cudaResult =JCufft.cufftExecC2C(plan, device_dicFFT.withByteOffset(i*whFFT*2*Sizeof.FLOAT),device_dicFFT.withByteOffset(i*whFFT*2*Sizeof.FLOAT),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 28 "+cudaResult);}

            //if (i==0){dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "dic1_"+range[i],"FLOAT");}
        
            
            //if (i==0){dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "fftImage2","FLOAT");}
        
            
            //IJ.log("make dic ok");
            
            
            
            
            
            cudaResult =cudaMemcpyAsync(device_dic, Pointer.to(psf_f_mirrored), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            
            cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), dparam.param.sizeoutput*dparam.param.sizeoutput, device_dic, device_sparseIndexPaddingShift2DEvenDic,device_dicFFT_mirrored.withByteOffset(i*whFFT*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
            cudaResult =JCufft.cufftExecC2C(plan, device_dicFFT_mirrored.withByteOffset(i*whFFT*2*Sizeof.FLOAT),device_dicFFT_mirrored.withByteOffset(i*this.whFFT*2*Sizeof.FLOAT),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 28 "+cudaResult);}
            
            //if (i==0){dparam.psf.imshow(whFFT*2,widthFFT, device_dicFFT_mirrored.withByteOffset(i*whFFT*2*Sizeof.FLOAT), "fftImage2","FLOAT");}
        
        
        }
        
        
        dparam.psf.imshow(dparam.param.sizeoutput, psftmp, "psf");
        
        
        JCuda.cudaFree(device_dic) ;
        JCuda.cudaFree(device_sparseIndexPaddingShift2DEvenDic) ;
        MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(streamId), whFFT*2*this.range.length, device_dicFFT, device_dicFFT, (float)(1./Math.sqrt(1)));//here we multiply with sqrt(total size) to avoid to do it after
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
    
    
    
    
    
    
    
    public boolean setBackgroundMorphology(int indexFrame){
        
        
        int smoothsize=(int)Math.sqrt(sizePSF);
        if (smoothsize%2==0){
            smoothsize+=1;
        }
        if (smoothsize<3){
            smoothsize=3;
        }
        //remove SLIGHTLY noise
        MyVecDouble.localMean(MyCudaStream.getCUstream(streamId), wh, width, height, smoothsize, this.device_tmp1, device_image.withByteOffset(wh*indexFrame*Sizeof.FLOAT));//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMinimum(MyCudaStream.getCUstream(streamId), this.wh, this.width, this.height, sizePSF, this.device_tmp2, device_tmp1);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMean(MyCudaStream.getCUstream(streamId), this.wh, this.width, this.height, smoothsize, device_tmp1, device_tmp2);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMaximum(MyCudaStream.getCUstream(streamId), this.wh, this.width, this.height, sizePSF, device_tmp2, device_tmp1);//here we multiply with sqrt(total size) to avoid to do it after
        
        MyVecDouble.localMean(MyCudaStream.getCUstream(streamId), this.wh, this.width, this.height, smoothsize, device_background.withByteOffset(wh*indexFrame*Sizeof.FLOAT), device_tmp2);//here we multiply with sqrt(total size) to avoid to do it after
        
        host_bckg[indexFrame]=Pointer.to(bckg[indexFrame]);
        cudaResult =cudaMemcpyAsync(host_bckg[indexFrame],device_background.withByteOffset(wh*indexFrame*Sizeof.FLOAT), wh*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
        
        //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_background, "device_bckg","FLOAT");
        
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    public boolean setImage(double [][] image){
        if (freeFramePosition.size()==0){
            IJ.log("ERROR: no free position to process image for deconvolution");
            return false;
        }
        int index=freeFramePosition.get(0);
        freeFramePosition.remove(0);
        
        this.width=image.length;
        this.height=image[0].length;
        
        //cudaResult=JCuda.cudaMemsetAsync(device_imageFFT.withByteOffset(whFFT*2*index*Sizeof.FLOAT), 0, whFFT*2 * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
       
        for (int i=0,j=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                this.image[index][j++]=(float)image[i][ii];
                if ((i<width)&&(ii<height)){
                    
                    
                }
                else{
                    this.image[index][j++]=0;//meanBckg;
                }
            }
        }
        
        
        cudaResult=cudaMemcpyAsync(device_image.withByteOffset(wh*index*Sizeof.FLOAT), this.host_image[index], wh*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        
        setBackgroundMorphology(index);
        
       
       cudaResult=cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), wh, device_image.withByteOffset(wh*index*Sizeof.FLOAT), this.device_sparseIndexPaddingShift2DEven,device_imageFFT.withByteOffset(whFFT*2*index*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 3");}
        
       cudaResult= JCufft.cufftExecC2C(plan, device_imageFFT.withByteOffset((whFFT*2*index*Sizeof.FLOAT)),device_imageFFT.withByteOffset((whFFT*2*index*Sizeof.FLOAT)),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}
       
       
       
       //dparam.psf.imshow(wh*nbFrame,height, device_image, "Image","FLOAT");
       //dparam.psf.imshow(wh*nbFrame,height, device_background, "bckg","FLOAT");
        
       JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        
       
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    
    public boolean setOinit(){
        
        double [][][][] o = new double [nbFrame][range.length][this.width][this.height];
        float [] fullO = new float[depth*nbFrame*width*height];
        
        try{
        o[0][0][60][60]=2000;
        o[0][1][20][20]=1500;
        o[1][0][20][61]=1000;
        }catch(Exception e){IJ.log("point not set");}
        
        for (int i=0;i<o.length;i++){
            
            for (int ii=0;ii<o[i].length;ii++){
                for (int iii=0;iii<o[i][ii].length;iii++){
                    for (int iiii=0;iiii<o[i][ii][iii].length;iiii++){
                
                        float v=(float)o[i][ii][iii][iiii];
                        fullO[i*o[i].length*o[i][ii].length*o[i][ii][iii].length+ii*o[i][ii].length*o[i][ii][iii].length+iii*o[i][ii][iii].length+iiii]=v+(float).1;
                    }
                }
            }
            
        }
        Pointer host_o=Pointer.to(fullO);
        cudaResult =cudaMemcpyAsync(device_o,host_o, wh*nbFrame*range.length*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        //dparam.psf.imshow(wh*range.length,(int)Math.sqrt(wh), device_o, "input","FLOAT");//show first image
        
        
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    
    
    //oIsNull: set it to false if o is already initialized
    //oIsNull: set it to true to initialise it to fixed positive value (FAST)
    public double [][] deconvolution_RL(){
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        
        //setOinit();
        
        MyVecDouble.initializeVectorToValue(MyCudaStream.getCUstream(streamId), wh*depth*nbFrame, device_o, (float)0.1);
            
        
        

        //sum psf*o
        
                
                
        float power=(float).125;
        //dparam.psf.imshow(totalSizeUnpadded*this.range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "oooooINIT ","FLOAT");
        
        
        long t1=System.currentTimeMillis();
        
        //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_background, "device_background mul","FLOAT");
        toto:for (int t=0;t<this.nb_iter;t++){
            
            
            //dparam.psf.imshow(wh*2,width, device_o.withByteOffset(wh*2*0*Sizeof.FLOAT), "init "+t,"FLOAT");
            
            
            
            
            
            cudaResult=JCuda.cudaMemsetAsync(device_tmpFFT1, 0, (long)this.whFFT*(long)2*(long)this.range.length*(long)nbFrame * (long)Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}



            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}

            MyVecDouble.mycusparsemoduloSsctrMany(MyCudaStream.getCUstream(streamId), (long)wh*(long)range.length*(long)nbFrame,nbFrame,depth,whFFT*2,width,height,device_tmpFFT1, device_o, device_sparseIndexPaddingShift2DEven);
            
            
            
            
            cudaResult= JCufft.cufftExecC2C(planDepthMany, device_tmpFFT1,device_tmpFFT1,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}

            
            
            
            MyVecDouble.complexeMulKernelMany(MyCudaStream.getCUstream(streamId), (long)whFFT*(long)depth*(long)nbFrame, depth,whFFT, device_tmpFFT2, device_tmpFFT1, device_dicFFT);
            
            
            
            
            JCufft.cufftExecC2C(planDepthMany, device_tmpFFT2,device_tmpFFT2,JCufft.CUFFT_INVERSE);

            
            
            
            MyVecDouble.makeResultCorrelationMany(MyCudaStream.getCUstream(streamId), (long)wh*(long)range.length*(long)nbFrame, nbFrame,depth,wh,whFFT*2, device_tmpFFT1, device_tmpFFT2, this.device_sparseIndexPaddingShift2DEven);

            
            //dparam.psf.imshow(wh*range.length*2,(int)Math.sqrt(wh), device_tmpFFT1, "device_tmpFFT1","FLOAT");//show first image
            
            
            MyVecDouble.turnMatrixMany(MyCudaStream.getCUstream(streamId), (long)wh*(long)depth*(long)nbFrame, nbFrame,depth,wh, device_tmpFFT2, device_tmpFFT1);

            if (nb_iter-1==t)
                dparam.psf.imshow(wh*range.length,(int)Math.sqrt(wh), device_tmpFFT1, "device_tmpFFT1","FLOAT");//show first image
            
            


            cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_N,wh*nbFrame,range.length,device_alpha,device_tmpFFT2,wh*nbFrame,device_a,1,device_beta,device_tmpFFT1,1);
            //cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_T,range.length,wh*nbFrame,device_alpha,device_tmpFFT2,wh*nbFrame,device_a,1,device_beta,device_tmpFFT1,1);
            
            
            if (nb_iter-1==t)
            dparam.psf.imshow(wh,height, device_tmpFFT1, "sum "+t,"FLOAT");
            
            
            
            //dparam.psf.imshow(wh*nbFrame,width, device_tmpFFT2, "device_tmpFFT1","FLOAT");//show first image
            
            
            
            //dparam.psf.imshow(wh*nbFrame,height, device_tmpFFT1, "data","FLOAT");//show first image

            
            //the following performs I/(p*o+b)
            MyVecDouble.addanddivide(MyCudaStream.getCUstream(streamId), wh*nbFrame, device_tmpFFT2, device_image , device_tmpFFT1,device_background);
            
            //dparam.psf.imshow(wh*nbFrame,height, device_tmpFFT2, "div","FLOAT");//show first image

            //dparam.psf.imshow(wh*nbFrame,height, device_image, "device_image","FLOAT");//show first image
            
            //dparam.psf.imshow(wh*nbFrame,height, device_tmpFFT2, "divFirst","FLOAT");//show first image

            
            
            MyVecDouble.copyMany(MyCudaStream.getCUstream(streamId), wh*nbFrame*depth,wh,depth,nbFrame, device_tmpFFT1, device_tmpFFT2 );

                
            //dparam.psf.imshow(wh*nbFrame*depth,height, device_tmpFFT1, "copy","FLOAT");//show first image

            
            
            
            
            //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}

            //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}

            //dparam.psf.imshow(wh*nbFrame,height, device_tmpFFT2, "devided","FLOAT");//show first image
            

            //FFT computation
            cudaResult=JCuda.cudaMemsetAsync(device_tmpFFT2, 0, (long)this.whFFT*(long)2*(long)this.range.length*(long)nbFrame * (long)Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}

            
            //if (y==9){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), "o_mid "+y,"FLOAT");}
            //if (y==0){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), "o ("+y,"FLOAT");}

            //here, device_convFFT corresponds to i/(p*o+b)
           //this was written to test FFT : cudaResult=cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), totalSizeUnpadded, device_o.withByteOffset(y*this.totalSizeUnpadded*Sizeof.FLOAT), this.device_sparseIndexPaddingShift2DEven,device_tmpFFT, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 3");}
           
           
           
           
           MyVecDouble.mycusparsemoduloSsctrMany(MyCudaStream.getCUstream(streamId), (long)wh*(long)range.length*(long)nbFrame,nbFrame,depth,whFFT*2,width,height,device_tmpFFT2, device_tmpFFT1, device_sparseIndexPaddingShift2DEven);
            
           
           //dparam.psf.imshow(whFFT*2*nbFrame*depth,heightFFT*2, device_tmpFFT2, "sparse","FLOAT");//show first image
           
           
           
           cudaResult= JCufft.cufftExecC2C(planDepthMany, device_tmpFFT2,device_tmpFFT2,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}

           //dparam.psf.imshow(whFFT*2*nbFrame,heightFFT*2, device_tmpFFT1, "fft","FLOAT");//show first image
           

           //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_dicFFT.withByteOffset(y*this.totalSize*2*Sizeof.FLOAT), "mulPSF 0 "+y,"FLOAT");

           MyVecDouble.complexeMulKernelMany(MyCudaStream.getCUstream(streamId), (long)whFFT*(long)depth*(long)nbFrame, depth,whFFT, device_tmpFFT1,device_tmpFFT2, device_dicFFT_mirrored);
            
            
           
           
           //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_convFFT, "device_convFFT 0 "+y,"FLOAT");

            JCufft.cufftExecC2C(planDepthMany, device_tmpFFT1,device_tmpFFT1,JCufft.CUFFT_INVERSE);

            //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(totalSize), device_convFFT, "device_convFFT 1 "+y,"FLOAT");


            MyVecDouble.makeResultCorrelationMany(MyCudaStream.getCUstream(streamId), (long)wh*(long)range.length*(long)nbFrame, nbFrame,depth,wh,whFFT*2, device_tmpFFT2, device_tmpFFT1, this.device_sparseIndexPaddingShift2DEven);
            
            //dparam.psf.imshow(wh*nbFrame*depth,height, device_tmpFFT2, "weight","FLOAT");//show first image
           
            
            
            //here, device_correlImage contains (i/(p*o+b)*p*)

            //if (y==9){dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "(i/(p*o+b)*p*)  "+t,"FLOAT");}

            //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "device_correlImage  "+y,"FLOAT");

            //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_correlImage, "device_correlImage "+y,"FLOAT");

            

            MyVecDouble.mul_fl_pow(MyCudaStream.getCUstream(streamId), (long)wh*(long)range.length*(long)nbFrame, device_o, device_o, device_tmpFFT2,power);

            
            
            
            
            
            
            
            //dparam.psf.imshow(totalSizeUnpadded*this.range.length,(int)Math.sqrt(totalSizeUnpadded), device_o, "ooooo "+t,"FLOAT");
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
            power*=2;
            
            power=Math.min(power, 3);
            
        }
        
        
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        long t2=System.currentTimeMillis();
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}


        IJ.log("time  "+(t2-t1));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        //long t1=System.currentTimeMillis();
        //IJ.log("time = "+(t1-t0));
        
        dparam.psf.imshow(wh*range.length*nbFrame,height, device_o, "o","FLOAT");//show first image
        
        dparam.psf.imshow(wh,height, device_background, "bckg","FLOAT");
        
        
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+streamId);}
        
        
        
        
        return null;
        
    }
    
    
    /*
    
    public void runCrossCorrelationOnGPU(){
        
        //cudaResult=JCuda.cudaMemsetAsync(device_convFFT, 0, totalSizeUnpadded*range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
               
        
        MyVecDouble.complexeConjugateKernel(MyCudaStream.getCUstream(streamId), totalSize*range.length, totalSize, device_tmpFFT, device_fullImageFFT, device_dicFFT);
        
        
        JCufft.cufftExecC2C(planMany, device_tmpFFT,device_tmpFFT,JCufft.CUFFT_INVERSE);
        
        float [] tmp = new float[1];
        
        //dparam.psf.imshow(totalSizeUnpadded,(int)Math.sqrt(totalSizeUnpadded), device_background, "device_background ","FLOAT");
        //MyVecDouble.mulScalarFloat(MyCudaStream.getCUstream(streamId), totalSizeUnpadded, device_background, device_background, (float).001);

        //dparam.psf.imshow(totalSize*2*range.length,(int)Math.sqrt(this.totalSize), device_dicFFT, "dic2","FLOAT");
        
        //MyVecDouble.makeResultCorrelation(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, totalSizeUnpadded,totalSize*2, device_correlImage, device_tmpFFT, this.device_sparseIndexPaddingShift2DEven);
        //use the following to use cross correl as input of deconvolution
        MyVecDouble.makeResultCorrelationNormalized(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, totalSizeUnpadded,totalSize*2, device_o, device_tmpFFT, this.device_sparseIndexPaddingShift2DEven,(float)0.999*(float)this.sizePSF,device_background,(float)0.1);
        
        IJ.log("This function should not be used because it is not normalized cross correlation (normalization is expensive on GPU)");
        
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
                
                if ((sum<al.get(i)[7]*RSB)||(sum<photonThreshold)){
                    al.remove(i);
                }
                else{
                    //IJ.log("al.get(i)[7] "+al.get(i)[7]+"  "+sum+"  "+i+"  al: "+(al.get(i)[7]*RSB));
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


    
    
        
        
        
        
    
    
    
    */
    
    
    public void free(){
        JCuda.cudaFree(device_background) ;
        
        JCuda.cudaFree(device_o) ;
        JCuda.cudaFree(device_tmp1) ;
        JCuda.cudaFree(device_tmp2) ;
        
        JCuda.cudaFree(device_tmpFFT1) ;
        JCuda.cudaFree(device_tmpFFT2) ;



        JCuda.cudaFree(device_a) ;
        JCuda.cudaFree(device_alpha) ;
        JCuda.cudaFree(device_beta) ;


        JCuda.cudaFree(device_sparseIndexPaddingShift2DEven) ;


        JCuda.cudaFree(device_imageFFT) ;
        
        JCuda.cudaFree(device_image) ;

        JCuda.cudaFree(device_dicFFT) ;
        
        JCuda.cudaFree(device_dicFFT_mirrored) ;



        JCufft.cufftDestroy(plan);
        plan=null; 
        JCufft.cufftDestroy(planDepthMany);
        planDepthMany=null; 
        
        JCufft.cufftDestroy(planMany);
        planMany=null; 

    }
    
}
