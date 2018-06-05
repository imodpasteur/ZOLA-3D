/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;

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
/**
 *
 * @author benoit
 */
public class PredetectionGPU {
    
    int cudaResult;
    cufftHandle plan;
    
    cufftHandle planMany;
    
    Pointer device_fullImage;
    
    Pointer device_correlImage;
    
    
    Pointer device_sparseIndexPaddingShift2DEven;
    
    
    
    
    Pointer device_fullImageFFT;
    
    Pointer device_dicFFT;
    
    
    Pointer device_tmpFFT;
    
    
    Pointer host_resultCorrelation;
    float [] resultCorrelation_f;
    float [][][] resultCorrelation;
    
    
    static double [] varPSF;
    
    
    
    float [] fullImage;
    double [] meanImage;
    double [] sumImage;
    double [] sumSquareImage;
    double [] sumImageSecond;
    double [] normImage;
    double [] sumSquareImageSecond;
    double [] varImage;
    
    int sizePSF;
    
    Pointer host_fullImage;
    
    
    int width;
    int height;
    int sizeFullImage;
    int totalSize;
    int totalSizeUnpadded;
    DataPhase dparam;
    
    public double [] range;
    
    double [][] psf;
    
    double [][] resvect ;
    
    int streamId;
    
    double thresholdCrossCorrelation=.2;
    
    
    double zstep;
    double xystep;
    public PredetectionGPU(int sizeFullImage, DataPhase dparam,double mini, double maxi, double step,double thresholdCrossCorrelation){

        double initStep=step;
        
        if (sizeFullImage%2==1){
            sizeFullImage++;
        }
        
        
        this.sizeFullImage=sizeFullImage;
        this.thresholdCrossCorrelation=thresholdCrossCorrelation;
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
        
        
        
        meanImage=new double [totalSizeUnpadded];
        sumImage=new double [totalSizeUnpadded];
        sumSquareImage=new double [totalSizeUnpadded];
        sumImageSecond=new double [totalSizeUnpadded];
        normImage=new double [totalSizeUnpadded];
        sumSquareImageSecond=new double [totalSizeUnpadded];
        varImage=new double [totalSizeUnpadded];
        fullImage=new float [totalSizeUnpadded];
        
        host_fullImage=Pointer.to(fullImage);
        
        
        
        
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
        
        
        for (int ii=0;ii<sizeFullImage;ii++){
            for (int i=0;i<sizeFullImage;i++){
                int valx1=sizePSF/2;
                int valx2=sizePSF/2;
                if (i<sizePSF/2-1){
                    valx1=i+1;
                }
                if (i>sizeFullImage-sizePSF/2-2){
                    valx2=sizeFullImage-i-1;
                }
                int valy1=sizePSF/2;
                int valy2=sizePSF/2;
                if (ii<sizePSF/2-1){
                    valy1=ii+1;
                }
                if (ii>sizeFullImage-sizePSF/2-2){
                    valy2=sizeFullImage-ii-2;
                }
                normImage[i*sizeFullImage+ii]=(valx1+valx2)*(valy1+valy2);
            }
        }
       
        
        device_sparseIndexPaddingShift2DEven=new Pointer();
        cudaResult =cudaMalloc(device_sparseIndexPaddingShift2DEven, totalSizeUnpadded * Sizeof.INT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 6 "+cudaResult);}
        cudaResult =cudaMemcpyAsync(device_sparseIndexPaddingShift2DEven, Pointer.to(sparseIndexPaddingShift2DEven), totalSizeUnpadded*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
        
        
        
        
        
        
        
        device_fullImage = new Pointer();
        cudaResult =cudaMalloc(device_fullImage, sizeFullImage*sizeFullImage * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        resultCorrelation_f=new float[sizeFullImage*sizeFullImage*range.length];
        resultCorrelation=new float[range.length][sizeFullImage][sizeFullImage];
        host_resultCorrelation = Pointer.to(resultCorrelation_f);
        
        device_correlImage = new Pointer();
        cudaResult =cudaMalloc(device_correlImage, sizeFullImage*sizeFullImage*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 8 "+cudaResult);}
        
        
        
        
        
        plan = new cufftHandle();
        
        cudaResult =JCufft.cufftPlan2d(plan, sizeFullImageFFT,sizeFullImageFFT, cufftType.CUFFT_C2C);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 11 "+cudaResult);}
        
        cudaResult =JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}
        
        
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
        
        
        device_tmpFFT = new Pointer();
        cudaResult =cudaMalloc(device_tmpFFT, totalSize*2*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 184 "+cudaResult);}
        
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        //IJ.log("warning no dictionnary");
        
         makeDictionnary();
        
         
         
         
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
        
        
        varPSF=new double [range.length];
        cudaResult =JCuda.cudaMemsetAsync(device_dicFFT, 0, this.totalSize*2*range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));
        //double [][][] psftmp = new double [range.length][dparam.param.sizeoutput][dparam.param.sizeoutput];
        double [] psf_d = new double [dparam.param.sizeoutput*dparam.param.sizeoutput];
        float [] psf_f = new float [dparam.param.sizeoutput*dparam.param.sizeoutput];
        psf= new double [range.length][dparam.param.sizeoutput*dparam.param.sizeoutput];
        for (int i=0;i<range.length;i++){
        //for (int i=0;i<1;i++){
            //IJ.log("range PSFmaker "+range[i]+"  "+dparam.param.Zfocus);
            
            
            dparam.psf.computePSF(0,0,dparam.param.Zfocus,range[i]);
            
            //conversion float
            psf_d=dparam.psf.getPSF1D();
            double mean=0;
            for (int u=0;u<psf_d.length;u++){
                psf[i][u]=psf_d[u];
                mean+=psf_d[u];
            }
            mean/=(double)psf_d.length;
            varPSF[i]=0;
            for (int u=0;u<psf_d.length;u++){
                psf_d[u]-=mean;
                psf_f[u]=(float)psf_d[u];
                varPSF[i]+=psf_d[u]*psf_d[u];
            }
            varPSF[i]/=(double)(psf_d.length);
            
            double 
            //dparam.psf.imshow(dparam.param.sizeoutput, psf_d, "DOUBLE");
            //dparam.psf.imshow(dparam.param.sizeoutput, psf_f, "FLOAT");
        
            
            //dparamFullImage.psf.imshow(this.totalSizeUnpadded,(int)Math.sqrt(this.totalSizeUnpadded), dparamFullImage.psf.getPointerPSF(), "ImagePSFgen","DOUBLE");
            
            cudaResult =cudaMemcpyAsync(device_dic, Pointer.to(psf_f), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
            
            //dparam.psf.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput,dparam.param.sizeoutput, device_dic, "dic","FLOAT");
            
            //dparam.psf.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput,dparam.param.sizeoutput, device_sparseIndexPaddingShift2DEvenDic, "sparse","INT");
        
            cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), dparam.param.sizeoutput*dparam.param.sizeoutput, device_dic, device_sparseIndexPaddingShift2DEvenDic,device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
            
            
            
            cudaResult =JCufft.cufftExecC2C(plan, device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 28 "+cudaResult);}

            //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "dic1_"+range[i],"FLOAT");
        
            
            //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), "fftImage2","FLOAT");
        
            
            //IJ.log("make dic ok");
        
        }
        
        JCuda.cudaFree(device_dic) ;
        JCuda.cudaFree(device_sparseIndexPaddingShift2DEvenDic) ;
        MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(streamId), this.totalSize*2*this.range.length, device_dicFFT, device_dicFFT, (float)Math.sqrt(totalSize));
        
        
        
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
    
    public boolean setImage(double [][] image){
        
        
        this.width=image.length;
        this.height=image[0].length;
        
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
                    fullImage[j++]=0;
                }
            }
        }
        
        
        
        
        
        cudaResult=cudaMemcpyAsync(device_fullImage, host_fullImage, totalSizeUnpadded*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult+"  "+totalSizeUnpadded);}
        
        
        
        
        
        
        
        cudaResult=cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), totalSizeUnpadded, device_fullImage, this.device_sparseIndexPaddingShift2DEven,device_fullImageFFT, CUSPARSE_INDEX_BASE_ZERO);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cusparseDgthr cuda 3");}
        
       cudaResult= JCufft.cufftExecC2C(plan, device_fullImageFFT,device_fullImageFFT,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 6");}
       
       
       //dparam.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_fullImageFFT, "fftImage","FLOAT");
        
       
       //double [][] tmp0 = new double [sizeFullImage][sizeFullImage];
       //double [][] tmp1 = new double [sizeFullImage][sizeFullImage];
       //double [][] tmp2 = new double [sizeFullImage][sizeFullImage];
       
        for (int i=0;i<sizeFullImage;i++){
            sumImage[i*sizeFullImage+0]=0;
            sumSquareImage[i*sizeFullImage+0]=0;
            for (int t=0;t<sizePSF/2+1;t++){//before 0 it is outside
                if (t<sizeFullImage){
                    sumImage[i*sizeFullImage+0]+=fullImage[i*sizeFullImage+t];
                    sumSquareImage[i*sizeFullImage+0]+=fullImage[i*sizeFullImage+t]*fullImage[i*sizeFullImage+t];
                }
            }
        
            for (int ii=1;ii<sizeFullImage;ii++){
                sumImage[i*sizeFullImage+ii]=sumImage[i*sizeFullImage+ii-1];
                if (ii-sizePSF/2>=0){
                    sumImage[i*sizeFullImage+ii]-=fullImage[i*sizeFullImage+ii-sizePSF/2];
                }
                else{
                    //outside
                }
                if (ii+sizePSF/2<sizeFullImage){
                    sumImage[i*sizeFullImage+ii]+=fullImage[i*sizeFullImage+ii+sizePSF/2];
                }
                else{
                    //outside
                }
                sumSquareImage[i*sizeFullImage+ii]=sumSquareImage[i*sizeFullImage+ii-1];
                if (ii-sizePSF/2>=0){
                    sumSquareImage[i*sizeFullImage+ii]-=fullImage[i*sizeFullImage+ii-sizePSF/2]*fullImage[i*sizeFullImage+ii-sizePSF/2];
                }
                if (ii+sizePSF/2<sizeFullImage){
                    sumSquareImage[i*sizeFullImage+ii]+=fullImage[i*sizeFullImage+ii+sizePSF/2]*fullImage[i*sizeFullImage+ii+sizePSF/2];
                }
                
            }
        }
        
        for (int ii=0;ii<sizeFullImage;ii++){
            sumImageSecond[0*sizeFullImage+ii]=0;
            sumSquareImageSecond[0*sizeFullImage+ii]=0;
            for (int t=0;t<sizePSF/2+1;t++){
                if (t<sizeFullImage){
                    sumImageSecond[0*sizeFullImage+ii]+=sumImage[t*sizeFullImage+ii];
                    sumSquareImageSecond[0*sizeFullImage+ii]+=sumSquareImage[t*sizeFullImage+ii];
                }
            }
            for (int i=1;i<sizeFullImage;i++){
                sumImageSecond[i*sizeFullImage+ii]=sumImageSecond[(i-1)*sizeFullImage+ii];
                if (i-sizePSF/2>=0){
                    sumImageSecond[i*sizeFullImage+ii]-=sumImage[(i-sizePSF/2)*sizeFullImage+ii];
                }
                if (i+sizePSF/2<sizeFullImage){
                    sumImageSecond[i*sizeFullImage+ii]+=sumImage[(i+sizePSF/2)*sizeFullImage+ii];
                }
                
                sumSquareImageSecond[i*sizeFullImage+ii]=sumSquareImageSecond[(i-1)*sizeFullImage+ii];
                if (i-sizePSF/2>=0){
                    sumSquareImageSecond[i*sizeFullImage+ii]-=sumSquareImage[(i-sizePSF/2)*sizeFullImage+ii];
                }
                if (i+sizePSF/2<sizeFullImage){
                    sumSquareImageSecond[i*sizeFullImage+ii]+=sumSquareImage[(i+sizePSF/2)*sizeFullImage+ii];
                }
                
            }
            for (int i=0;i<sizeFullImage;i++){
                meanImage[i*sizeFullImage+ii]=sumImageSecond[i*sizeFullImage+ii]/normImage[i*sizeFullImage+ii];
            }
        }
        
        
        
        for (int i=0;i<sizeFullImage;i++){
            for (int ii=0;ii<sizeFullImage;ii++){
                varImage[i*sizeFullImage+ii]=(sumSquareImageSecond[i*sizeFullImage+ii]/normImage[i*sizeFullImage+ii])-meanImage[i*sizeFullImage+ii]*meanImage[i*sizeFullImage+ii];
            }
        }
        
        //for (int ii=0;ii<sizeFullImage;ii++){
        //    for (int i=0;i<sizeFullImage;i++){
        //        tmp0[i][ii]=meanImage[i*sizeFullImage+ii];
        //        tmp1[i][ii]=varImage[i*sizeFullImage+ii];
        //        tmp2[i][ii]=normImage[i*sizeFullImage+ii];
        //    }
        //}
        //classes.ImageShow.imshow(tmp0,"mean_image");
        //classes.ImageShow.imshow(tmp1,"var_image");
        //classes.ImageShow.imshow(tmp2,"norm_image");
        //dparam.psf.imshow(sizeFullImage, sumSquareImageSecond, "sumSquareImageSecond");
        
        //dparam.psf.imshow(sizeFullImage, varImage, "myVar");
       
       
        return (cudaResult == cudaError.cudaSuccess);
    }
    
    
    
    
    
    
    
    
    public double [][] convolveNormalizedInFourierDomain(){
        
        
        
        MyVecDouble.complexeConjugateKernel(MyCudaStream.getCUstream(streamId), totalSize*range.length, totalSize, device_tmpFFT, device_fullImageFFT, device_dicFFT);
        
        
        JCufft.cufftExecC2C(planMany, device_tmpFFT,device_tmpFFT,JCufft.CUFFT_INVERSE);
        
        
        
        //dparam.psf.imshow(totalSize*2*range.length,(int)Math.sqrt(this.totalSize), device_dicFFT, "dic2","FLOAT");
        
        MyVecDouble.makeResultCorrelation(MyCudaStream.getCUstream(streamId), totalSizeUnpadded*range.length, totalSizeUnpadded,totalSize*2, device_correlImage, device_tmpFFT, this.device_sparseIndexPaddingShift2DEven);
        
        //dparam.psf.imshow(totalSizeUnpadded*range.length,(int)Math.sqrt(this.totalSizeUnpadded), device_correlImage, "device_correlImage","FLOAT");
        
        
        cudaResult =cudaMemcpyAsync(this.host_resultCorrelation, device_correlImage, totalSizeUnpadded*range.length*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
        
        
        for (int z=0;z<range.length;z++){
            for (int i=0;i<this.width;i++){
                for (int ii=0;ii<this.height;ii++){
                    if (Math.sqrt(varPSF[z]*varImage[i*sizeFullImage+ii])>0){
                        this.resultCorrelation[z][i][ii]=(float)(((double)this.resultCorrelation_f[z*sizeFullImage*sizeFullImage+i*sizeFullImage+ii])/Math.sqrt(varPSF[z]*varImage[i*sizeFullImage+ii]));
                    }
                    else if (varPSF[z]<0){
                        this.resultCorrelation[z][i][ii]=-1;
                    }
                    else if (varImage[i*sizeFullImage+ii]<0){
                        this.resultCorrelation[z][i][ii]=-2;
                    }
                }   
            } 
        }
        
        
        //ImageShow.imshow(resultCorrelation, "cross");
        
        
        
        
        double [][] vect= getVectMaxPosition(resultCorrelation);
        
        Arrays.sort(vect, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o2[0]).compareTo(o1[0]);
            }
        });
        return vect;
    }
    
    
    
    
    
        
        double [][] getVectMaxPosition(float [][][] resultConvolution){


            
            int decalFilter=3;
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
                        if ((resultConvolution[r][i][j]>thresholdCrossCorrelation)){
                            
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
                                        float offsetX=0;
                                        float offsetY=0;
                                        float offsetZ=0;
                                        float count=0;
                                        int shift=0;
                                        for (int a=0;a<=0;a++){
                                            if ((r+a>=0)&&(r+a<resultConvolution.length)){
                                                for (int aa=-shift;aa<=shift;aa++){
                                                    if ((i+aa>=0)&&(i+aa<resultConvolution[r].length)){
                                                        for (int aaa=-shift;aaa<=shift;aaa++){
                                                            if ((j+aaa>=0)&&(j+aaa<resultConvolution[r][i].length)){
                                                                offsetZ+=((float)(a*zstep))*resultConvolution[r+a][i+aa][j+aaa];
                                                                offsetX+=((float)(aa*xystep))*resultConvolution[r+a][i+aa][j+aaa];
                                                                offsetY+=((float)(aaa*xystep))*resultConvolution[r+a][i+aa][j+aaa];
                                                                count+=resultConvolution[r+a][i+aa][j+aaa];
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        offsetX/=count;
                                        offsetY/=count;
                                        offsetZ/=count;
                                        machin[0]=resultConvolution[r][i][j];
                                        machin[1]=r;
                                        machin[2]=i;
                                        machin[3]=j;
                                        machin[4]=offsetZ;
                                        machin[5]=offsetX;
                                        machin[6]=offsetY;
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


    
    
        
        
        
        
        double [][] getVectMaxPositionOld(float [][][] resultConvolution){



            int decalFilter=3;
            ArrayList<double []> al = new ArrayList<double []>();

            for (int r=0;r<resultConvolution.length;r++){
                toto:for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                    for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                        if ((resultConvolution[r][i][j]>thresholdCrossCorrelation)){
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
         JCuda.cudaFree(device_fullImage) ;
         
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
