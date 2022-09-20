/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.cublasHandle;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import jcuda.runtime.cudaError;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice;
import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.data.FrameLocalization;
import org.pasteur.imagej.data.PLocalization;
import org.pasteur.imagej.utils.Matrixe;



/**
 *
 * @author benoit
 */
public class LocalizationHD_ {
    
    int numframe;
    cublasHandle  handlecublas;
    
    
    int index_photon=4;
    int index_X=5;//int x position in image
    int index_Y=6;//int y position in image
    int index_x=0;//float x shift in (µm)
    int index_y=1;//float y shift in (µm)
    int index_z=2;//float z shift in (µm)
    int index_f=3;//focus in (µm)
    int index_crlb_x=7;
    int index_crlb_y=8;
    int index_crlb_z=9;
    int index_chi2=10;
    
    
    float [] position;
    
    int sizePatch;
    
    int cudaResult;
    
    Pointer device_image;
    
    Pointer device_model;
    
    Pointer device_likelihood;
    
    Pointer device_crop;
    Pointer device_crop_partial_x;
    Pointer device_crop_partial_y;
    Pointer device_crop_partial_z;
    Pointer device_crop_partial_photon;
    Pointer device_fisher;
    Pointer device_fisher_sum;
    Pointer device_crop_tmp_bckg;
    Pointer device_crop_tmp1;
    Pointer device_crop_tmp2;
    
    Pointer device_cropFloat;
    
    Pointer device_background;
    
    Pointer device_PSF;
    
    Pointer device_parameters;
    
    Pointer device_parameterSave;
    
    Pointer device_gamma;
    Pointer device_sum_gamma;
    float [] sum_gamma = new float[1];
    Pointer device_likelihoodResult;
    Pointer device_likelihoodResult1;
    Pointer device_likelihoodResult2;
    Pointer device_likelihoodResult3;
    
    Pointer device_invInfo;
    
    Pointer device_ones;
    Pointer device_alpha;
    Pointer device_beta;
    
    Pointer device_ones_f;
    Pointer device_alpha_f;
    Pointer device_beta_f;
    
    Pointer host_position;
    int wh;
    int width,height;
    
    DataPhase_ dp;
    
    double [] A;
    double [] B;
    double [] X;
    double [] Y;
    int [] Xint;//X position = (Xint-X)
    int [] Yint;
    double [] Z;
    double [] Zfocus;
    int streamId;
    
    int maxMoleculePerImage;
    int numberToCompute;//number of psf to compute for current image (should be <=maxMoleculePerImage)
    
    double [] rangePredetection;
    
    float [] fisher ;
    double [][] fisher2=new double[4][4];
    double [] std =new double[4];
    double [] score;//chi2
    
    double minZ;
    double maxZ;
            
    public LocalizationHD_(DataPhase_ dp,int maxMoleculePerImage,int width,int height,double minZ,double maxZ){
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.streamId=dp.getStream();
        this.width=width;
        this.height=height;
        this.wh=width*height;
        this.dp=dp;
        
        sizePatch=dp.param.sizeoutput;
        
        this.maxMoleculePerImage=maxMoleculePerImage;
        
        dp.setManyPSF(maxMoleculePerImage);
        
        Zfocus=new double[maxMoleculePerImage];
        for (int i=0;i<maxMoleculePerImage;i++){
            Zfocus[i]=(float)dp.param.Zfocus;
        }
        
        
        X = new double [maxMoleculePerImage];
        Y = new double [maxMoleculePerImage];
        Z = new double [maxMoleculePerImage];
        Xint = new int [maxMoleculePerImage];
        Yint = new int [maxMoleculePerImage];
        
        fisher = new float[maxMoleculePerImage*16];
        score= new double[maxMoleculePerImage];
        position=new float [11*maxMoleculePerImage];
        
        host_position=Pointer.to(position);
        
        
        device_image = new Pointer();
        cudaResult =cudaMalloc(device_image, wh * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 0 "+cudaResult);}
        
        device_likelihood = new Pointer();
        cudaResult =cudaMalloc(device_likelihood, wh * Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 0 "+cudaResult);}
        
        device_model = new Pointer();
        cudaResult =cudaMalloc(device_model, wh * Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 0 "+cudaResult);}
        
        device_background = new Pointer();
        cudaResult =cudaMalloc(device_background, wh * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 1 "+cudaResult);}
        
        device_parameters=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_parameters, 8*maxMoleculePerImage * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}

        device_gamma=new Pointer();//used to weight x,y,z,photon during update
        cudaResult = JCuda.cudaMalloc(device_gamma, 6*maxMoleculePerImage * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}

        device_sum_gamma=new Pointer();//used to weight x,y,z,photon during update
        cudaResult = JCuda.cudaMalloc(device_sum_gamma, 1 * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}

        
        
        
        
        
        device_crop=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        device_crop_tmp_bckg=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_tmp_bckg, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        device_crop_tmp1=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_tmp1, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        device_crop_tmp2=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_tmp2, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        
        device_crop_partial_x=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_partial_x, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        device_crop_partial_y=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_partial_y, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        
        device_crop_partial_z=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_partial_z, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        device_crop_partial_photon=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_crop_partial_photon, wh *maxMoleculePerImage* Sizeof.DOUBLE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        
        //4*4 fisher matrix
        device_fisher=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_fisher, wh *maxMoleculePerImage*16* Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        device_fisher_sum=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_fisher_sum, maxMoleculePerImage*16* Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
        
        
        
        device_cropFloat=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_cropFloat, wh *maxMoleculePerImage* Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}

        
        
        double [] ones=new double[sizePatch*sizePatch];
        for (int i=0;i<sizePatch*sizePatch;i++){
            ones[i]=1;
        }
        Pointer h_ones=Pointer.to(ones);
        device_ones = new Pointer();
        JCuda.cudaMalloc(device_ones, sizePatch*sizePatch *Sizeof.DOUBLE);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_ones, h_ones, sizePatch*sizePatch*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(dp.param.stream));
        
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.DOUBLE);
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(dp.param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.DOUBLE);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(dp.param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        
        
        
        
        device_invInfo =  new Pointer();
        cudaMalloc(device_invInfo, maxMoleculePerImage * Sizeof.INT);
        
        
        float [] ones_f=new float[sizePatch*sizePatch];
        for (int i=0;i<sizePatch*sizePatch;i++){
            ones_f[i]=1;
        }
        Pointer h_ones_f=Pointer.to(ones_f);
        device_ones_f = new Pointer();
        JCuda.cudaMalloc(device_ones_f, sizePatch*sizePatch *Sizeof.FLOAT);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_ones_f, h_ones_f, sizePatch*sizePatch*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(dp.param.stream));
        
        device_alpha_f = new Pointer();
        cudaMalloc(device_alpha_f, 1 * Sizeof.FLOAT);
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_alpha_f, Pointer.to(new float[]{(float)1.0}), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(dp.param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        device_beta_f = new Pointer();
        cudaMalloc(device_beta_f, 1 * Sizeof.FLOAT);
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_beta_f, Pointer.to(new float[]{(float)0.0}), 1*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(dp.param.stream));
        //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        
        
        
        
        
        device_likelihoodResult=new Pointer();
        JCuda.cudaMalloc(device_likelihoodResult, maxMoleculePerImage* Sizeof.DOUBLE);
        
        device_likelihoodResult1=new Pointer();
        JCuda.cudaMalloc(device_likelihoodResult1, maxMoleculePerImage* Sizeof.DOUBLE);
        
        device_likelihoodResult2=new Pointer();
        JCuda.cudaMalloc(device_likelihoodResult2, maxMoleculePerImage* Sizeof.DOUBLE);
        
        device_likelihoodResult3=new Pointer();
        JCuda.cudaMalloc(device_likelihoodResult3, maxMoleculePerImage* Sizeof.DOUBLE);
        
        device_parameterSave=new Pointer();
        JCuda.cudaMalloc(device_parameterSave, maxMoleculePerImage* Sizeof.FLOAT);
        
        
        
        handlecublas=MyCudaStream.getHandleCublas(dp.param.stream);
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        
        
    }
    
    public void setRangePredetection(double [] range){
        this.rangePredetection=new double [range.length];
        for (int i=0;i<range.length;i++){
            rangePredetection[i]=range[i];
        }
        
        
    }
    
    public void setImage(Pointer device_image,int numframe){
        this.numframe=numframe;
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        cudaResult=cudaMemcpyAsync(this.device_image, device_image, wh*Sizeof.FLOAT, cudaMemcpyDeviceToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
    }
    
    public void setCoordinates(double [][] vect){
        this.numberToCompute=Math.min(vect.length,maxMoleculePerImage);
        
        IJ.log("number"+numberToCompute+"  "+vect.length+"  "+maxMoleculePerImage);
        for (int i=0;i<numberToCompute;i++){
            /*this.A[i]=vect[i][0];
            this.B[i]=0;*/
            this.X[i]=-vect[i][5];
            this.Y[i]=-vect[i][6];
            this.Xint[i]=(int)vect[i][2];
            this.Yint[i]=(int)vect[i][3];
            this.Z[i]=vect[i][1]+vect[i][4];
            
        }
        
        
        for (int i=0;i<numberToCompute;i++){
            this.position[i+index_x*numberToCompute]=(float)-vect[i][5];//x relative
            this.position[i+index_y*numberToCompute]=(float)-vect[i][6];//y relative
            this.position[i+index_z*numberToCompute]=(float)(vect[i][1]+vect[i][4]);//z
            this.position[i+index_f*numberToCompute]=(float)(dp.param.Zfocus);//z focus
            this.position[i+index_photon*numberToCompute]=(float)(vect[i][0]);//photon
            this.position[i+index_X*numberToCompute]=(float)(vect[i][2]-sizePatch/2);//X position (int)
            this.position[i+index_Y*numberToCompute]=(float)(vect[i][3]-sizePatch/2);//Y position (int)
            
            IJ.log("zz:"+(vect[i][1]+vect[i][4])+"    photon:"+(vect[i][0]));
            //this.position[i+7*numberToCompute]=(float)(0);//bckg
            //IJ.log("i"+i+"    X"+vect[i][2]+"  Y"+vect[i][3]+"    Z"+(vect[i][1]+vect[i][4])+"   "+vect[i][0]);
            
        }
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        cudaMemcpyAsync(this.device_parameters, host_position, 8*numberToCompute*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
    }
    
    
    
    public void setBackground(Pointer device_bckg){
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        cudaResult=cudaMemcpyAsync(this.device_background, device_bckg, wh*Sizeof.FLOAT, cudaMemcpyDeviceToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
        //dp.psf.imshow(wh,height, device_bckg, "bckg","FLOAT");//show first image
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
    }
    
    
    public int getStreamId(){
        return dp.param.stream;
    }
    
    public FrameLocalization optimize(int nb_iter){
        
            
        //IJ.log("optimize");
        if (numberToCompute>0){
            //needed to compute likelihood (has to be launched after setImage and setCoordinates)
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            MyVecDouble.cropFromImageFloat(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch,sizePatch,numberToCompute,device_cropFloat,width,height,device_image,device_parameters,index_X,index_Y);//6 and 7 are index of X/Y patch starting
            //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_cropFloat, "image","FLOAT");//show first image


            FrameLocalization fl = new FrameLocalization(numframe);


            float h=(float)0.0001;
            //nb_iter=0;IJ.log("ERROR: nbiter=0");
            //start with z because z predetection may not be very precise
            
            
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            cudaMemcpyAsync(host_position,this.device_parameters, 5*numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            
            
            
            for (int i=0;i<numberToCompute;i++){
                if (true){//if ((Math.abs(position[i+index_x*numberToCompute])<dp.param.xystep*5)&&(Math.abs(position[i+index_y*numberToCompute])<dp.param.xystep*5)&&(position[i+index_z*numberToCompute]>minZ)&&(position[i+index_z*numberToCompute]<maxZ)){
                    
//                    if (i==11)
//                        IJ.log("XYZXYZ "+i+"   "+position[i+index_X*numberToCompute]+"    "+position[i+index_Y*numberToCompute]+"   "+position[i+index_z*numberToCompute]+"    "+position[i+index_photon*numberToCompute]);
                }
                                                        
            }
            
            
            updatePhotonParameter((float).001);
            updatePhotonParameter((float).001);
            updatePhotonParameter((float).001);
            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ   "+zzz[0]);}
            updateParameter(index_z,h);
            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ   "+zzz[0]);}
            updatePhotonParameter((float).001);
            updateParameter(index_z,h);
            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ   "+zzz[0]);}
            updatePhotonParameter((float).001);
            updateParameter(index_z,h);
            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ   "+zzz[0]);}
            updateParameter(index_z,h);
            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ   "+zzz[0]);}
            updateParameter(index_z,h);
            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ   "+zzz[0]);}
            
            for (int y=0;y<nb_iter;y++){
                
                            
                updateParameter(index_z,h);
                updatePhotonParameter((float).1);
                            if (false){float [] zzz=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(zzz), device_parameters.withByteOffset(index_z*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            float [] ppp=new float[numberToCompute];
                            cudaMemcpyAsync(Pointer.to(ppp), device_parameters.withByteOffset(index_photon*numberToCompute*Sizeof.FLOAT), numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
                            IJ.log("ZZZ ("+y+")  "+zzz[0]+"   "+ppp[0]);}
                
                updateParameter(index_x,h);
                updateParameter(index_y,h);
                
                            
                
            }
            
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            cudaMemcpyAsync(host_position,this.device_parameters, 5*numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            
            
            computeFisher();
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            int count=0;
            
            dp.psf.imshow(wh,height, device_model, "model ML","DOUBLE");//show first image
            
            for (int i=0;i<numberToCompute;i++){
                if (true){//((Math.abs(position[i+index_x*numberToCompute])<dp.param.xystep*5)&&(Math.abs(position[i+index_y*numberToCompute])<dp.param.xystep*5)&&(position[i+index_z*numberToCompute]>minZ)&&(position[i+index_z*numberToCompute]<maxZ)){
                    count++;
                    double my_x=1000*(position[i+index_X*numberToCompute]*dp.param.xystep+(dp.param.sizeoutput/2)*dp.param.xystep-position[i+index_x*numberToCompute]);
                    double my_y=1000*(position[i+index_Y*numberToCompute]*dp.param.xystep+(dp.param.sizeoutput/2)*dp.param.xystep-position[i+index_y*numberToCompute]);
                    double my_z=(1000*position[i+index_z*numberToCompute]);
                    //IJ.log("final Z "+my_z+"   "+position[i+index_photon*numberToCompute]);
                    double my_A=position[i+index_photon*numberToCompute];
                    double my_B=0;                                                        
                    double my_Score= position[i+index_chi2*numberToCompute];
                    double my_crlbx=(1000*position[i+index_crlb_x*numberToCompute]);
                    double my_crlby=(1000*position[i+index_crlb_y*numberToCompute]);
                    double my_crlbz=(1000*position[i+index_crlb_z*numberToCompute]);
                    PLocalization p = new PLocalization(i,numframe,my_x,my_y,my_z,my_A,my_B,my_Score,my_crlbx,my_crlby,my_crlbz);
                    //p.addListOfVariables_String(frameVariable);
                    //p.setValueOtherVariable_String(0, frameName);
                    fl.loc.add(p); 
                }
                                                        
            }
            IJ.log("frame:"+fl.numFrame+"    loc:"+count);
            return fl;
        }
        else{
            return null;
        }
    }
    
    
    //before calling this function, device_likelihoodResult2 has to be computed (current likelihood)
    //this function update device_likelihoodResult2 for the next call
    private void updateParameter(int indexParameter,float h){
        
        
        
        
        
        float maxJump=(float)(dp.param.xystep/2);//if x,y or z index -> max grad jump = pixel size / 2
        float minJump=(float)0.0000000001;
        
        if (indexParameter==4){//if photon index -> max grad jump = infinity
            maxJump=1000000000;
            
        }
        
        MyVecDouble.initializeVectorToValue(MyCudaStream.getCUstream(streamId), 6*maxMoleculePerImage, device_gamma, (float)10);
        
        
        
        
        
        
        this.computeLikelihoodAndBackground(device_likelihoodResult2);
        
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,indexParameter,-h,device_parameters);//6 and 7 are index of X/Y patch starting
        this.computeLikelihoodWithPredefinedBackground(device_likelihoodResult1);

        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,indexParameter,h*2,device_parameters);//6 and 7 are index of X/Y patch starting
        this.computeLikelihoodWithPredefinedBackground(device_likelihoodResult3);

        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,indexParameter,-h,device_parameters);//6 and 7 are index of X/Y patch starting

//        if (indexParameter==4){
//                                            double [] lik1=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik1), device_likelihoodResult1, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            double [] lik2=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik2), device_likelihoodResult2, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            double [] lik3=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik3), device_likelihoodResult3, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            IJ.log("lik ("+indexParameter+")    "+0+"  "+lik1[0]+"  "+lik2[0]+"  "+lik3[0]);
//
//                                            {float [] p=new float[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(p), device_parameters.withByteOffset(indexParameter*numberToCompute*Sizeof.FLOAT), 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//                                            IJ.log("param="+p[0]);}
//
//        }
        //gamma is divided by 10 at each iteration ; gamma = 0 if update is ok --> no more update
        for (int precision_loop=0;precision_loop<5;precision_loop++){
            
            MyVecDouble.updateParameter(MyCudaStream.getCUstream(streamId), numberToCompute,indexParameter,h,device_likelihoodResult1,device_likelihoodResult2,device_likelihoodResult3,device_parameters,device_parameterSave,device_gamma,minJump,maxJump);//6 and 7 are index of X/Y patch starting
            this.computeLikelihoodWithPredefinedBackground(device_likelihoodResult);
            
//            if (indexParameter==4){
//                                            double [] lik=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik), device_likelihoodResult, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            IJ.log("new lik ("+precision_loop+")    "+0+"  "+lik[0]);
//                                            
//                                            {float [] p=new float[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(p), device_parameters.withByteOffset(indexParameter*numberToCompute*Sizeof.FLOAT), 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//                                            IJ.log("param="+p[0]);}
//            }
            MyVecDouble.checkLikelihood(MyCudaStream.getCUstream(streamId), numberToCompute,indexParameter,device_likelihoodResult2,device_likelihoodResult,device_parameters,device_parameterSave,device_gamma);//6 and 7 are index of X/Y patch starting
            
            
            
//            float [] lik = new float[this.maxMoleculePerImage];
//            //print likelihood
//            cudaMemcpyAsync(Pointer.to(lik), device_likelihoodResult2, this.numberToCompute*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//            if (numberToCompute>=3){
//                IJ.log("likA "+lik[0]+"   "+lik[1]+"  "+lik[2]);
//            }
            
            
            
            //check gamma: may improve speed
            cudaResult=jcuda.jcublas.JCublas2.cublasSasum(handlecublas,numberToCompute,device_gamma.withByteOffset(indexParameter*numberToCompute*Sizeof.FLOAT),1,device_sum_gamma);
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            cudaMemcpyAsync(Pointer.to(sum_gamma), device_sum_gamma, 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            
            
            if (sum_gamma[0]<.0000000000001){//may run faster
                break;
            }
            
            
            
            
            //if device_likelihoodResult<device_likelihoodResult2 -> device_likelihoodResult2=device_likelihoodResult
            
            
        }
        
    }
    
    
    
    //before calling this function, device_likelihoodResult2 has to be computed (current likelihood)
    //this function update device_likelihoodResult2 for the next call
    private void updatePhotonParameter(float h){
        float maxJump=(float)1000000000;//if x,y or z index -> max grad jump = pixel size / 2
        float minJump=(float)0.0000000001;
        
        
        
        MyVecDouble.initializeVectorToValue(MyCudaStream.getCUstream(streamId), 6*maxMoleculePerImage, device_gamma, (float)1);
        
        
        
        
        
        
        this.computeLikelihoodAndBackground(device_likelihoodResult2);
        
        
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,-h,device_parameters);//6 and 7 are index of X/Y patch starting
        this.computeLikelihoodWithPredefinedBackgroundAndPredefinedPsf(device_likelihoodResult1);

        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,h*2,device_parameters);//6 and 7 are index of X/Y patch starting
        this.computeLikelihoodWithPredefinedBackgroundAndPredefinedPsf(device_likelihoodResult3);

        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,-h,device_parameters);//6 and 7 are index of X/Y patch starting
//        {
//                                            double [] lik1=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik1), device_likelihoodResult1, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            double [] lik2=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik2), device_likelihoodResult2, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            double [] lik3=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik3), device_likelihoodResult3, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            IJ.log("lik   "+11+"  "+lik1[11]+"  "+lik2[11]+"  "+lik3[11]);
//
//                                            {float [] p=new float[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(p), device_parameters.withByteOffset(index_photon*numberToCompute*Sizeof.FLOAT), 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//                                            IJ.log("param="+p[11]);}
//
//        }
        //gamma is divided by 10 at each iteration ; gamma = 0 if update is ok --> no more update
        for (int precision_loop=0;precision_loop<5;precision_loop++){
            
            MyVecDouble.updateParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,h,device_likelihoodResult1,device_likelihoodResult2,device_likelihoodResult3,device_parameters,device_parameterSave,device_gamma,minJump,maxJump);//6 and 7 are index of X/Y patch starting
            this.computeLikelihoodWithPredefinedBackground(device_likelihoodResult);
            
//            {
//                                            double [] lik=new double[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(lik), device_likelihoodResult, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//
//                                            IJ.log("new lik ("+precision_loop+")    "+11+"  "+lik[11]);
//                                            
//                                            {float [] p=new float[numberToCompute];
//                                            cudaMemcpyAsync(Pointer.to(p), device_parameters.withByteOffset(index_photon*numberToCompute*Sizeof.FLOAT), 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
//                                            IJ.log("param="+p[11]);}
//            }
            MyVecDouble.checkLikelihood(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,device_likelihoodResult2,device_likelihoodResult,device_parameters,device_parameterSave,device_gamma);//6 and 7 are index of X/Y patch starting
            
            
            
            
            
            //check gamma: may improve speed
            cudaResult=jcuda.jcublas.JCublas2.cublasSasum(handlecublas,numberToCompute,device_gamma.withByteOffset(index_photon*numberToCompute*Sizeof.FLOAT),1,device_sum_gamma);
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            cudaMemcpyAsync(Pointer.to(sum_gamma), device_sum_gamma, 1*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
            if (sum_gamma[0]<.0000000000001){//may run faster
                break;
            }
            
            
            
            //if device_likelihoodResult<device_likelihoodResult2 -> device_likelihoodResult2=device_likelihoodResult
            
            
        }
        
        MyVecDouble.max(MyCudaStream.getCUstream(streamId), numberToCompute, device_parameters.withByteOffset(index_photon*numberToCompute*Sizeof.FLOAT), (float)0);
        
        
    }
    
    
    private void computeLikelihoodWithPredefinedBackground(Pointer device_lik){
        computeLikelihoodWithPredefinedBackground(device_lik,false);
    }
    
    private void computeLikelihoodWithPredefinedBackground(Pointer device_lik,boolean show){
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        
        //here, device_cropFloat has to  contain image crops
        
        MyVecDouble.computeLikelihoodAndModelwithPhotonNumberAndBackground(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop,device_cropFloat,device_crop_tmp1,device_PSF,device_parameters,index_photon, device_crop_tmp_bckg);//6 and 7 are index of X/Y patch starting
        
        //here, device_crop_tmp1 is the model
         
        //here, device_crop is the likelihood
        
        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizePatch*sizePatch,numberToCompute,device_alpha,device_crop,sizePatch*sizePatch,device_ones,1,device_beta,device_lik,1);
        

        if (show){
            //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop, "likelihood 2","DOUBLE");//show first image
            
            dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_tmp1, "model 2","DOUBLE");//show first image
            
            
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        }
    }
    
    //same function but we do not recompute PSF that is unchanged (faster)
    private void computeLikelihoodWithPredefinedBackgroundAndPredefinedPsf(Pointer device_lik){
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        
        //here, device_cropFloat has to  contain image crops
        
        MyVecDouble.computeLikelihoodAndModelwithPhotonNumberAndBackground(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop,device_cropFloat,device_crop_tmp1,device_PSF,device_parameters,index_photon, device_crop_tmp_bckg);//6 and 7 are index of X/Y patch starting
        
        //here, device_crop_tmp1 is the model
         
        //here, device_crop is the likelihood
        
        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizePatch*sizePatch,numberToCompute,device_alpha,device_crop,sizePatch*sizePatch,device_ones,1,device_beta,device_lik,1);
        
    }
    
    
    private void computeLikelihoodAndBackground(Pointer device_lik){
        computeLikelihoodAndBackground(device_lik,false);
    }
    
    private void computeLikelihoodAndBackground(Pointer device_lik,boolean show){
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_PSF, "psf","FLOAT");//show first image
        //dp.psf.imshow(8*numberToCompute,numberToCompute, device_parameters, "device_parameters","FLOAT");//show first image

        


        //ici, je devrais faire l'inverse: mettre les psf dans une image plutot que d'utiliser des crops de l'image
        //Je pense qu'il faut pas recalculer le background: utiliser celui de deconvolution

        //dp.psf.imshow(wh,height, device_image, "imageinput","FLOAT");


        MyVecDouble.computeModelAndLikelihood(MyCudaStream.getCUstream(streamId), wh,width,height,sizePatch,sizePatch,numberToCompute,device_likelihood,device_model,device_image,device_PSF,device_parameters,index_X,index_Y,index_photon,device_background);//6 and 7 are index of X/Y patch starting
        //dp.psf.imshow(wh,height, device_model, "model","DOUBLE");

        //dp.psf.imshow(wh,height, device_likelihood, "lik","FLOAT");

        MyVecDouble.cropFromImage(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch,sizePatch,numberToCompute,device_crop,width,height,device_likelihood,device_parameters,index_X,index_Y);//6 and 7 are index of X/Y patch starting
        
        
        
        
        
        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizePatch*sizePatch,numberToCompute,device_alpha,device_crop,sizePatch*sizePatch,device_ones,1,device_beta,device_lik,1);
        
        MyVecDouble.cropFromImage(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch,sizePatch,numberToCompute,device_crop_tmp1,width,height,device_model,device_parameters,index_X,index_Y);//6 and 7 are index of X/Y patch starting
        MyVecDouble.subtractModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp_bckg,device_crop_tmp1,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        //here, device_crop_tmp_bckg contains background + overlapping emitters only
        
        if (show){
            //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop, "likelihood 1","DOUBLE");//show first image
            
            //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_PSF, "psf","FLOAT");
            
            MyVecDouble.cropFromImage(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch,sizePatch,numberToCompute,device_crop,width,height,device_model,device_parameters,index_X,index_Y);//6 and 7 are index of X/Y patch starting
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            
            dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop, "model 1","DOUBLE");//show first image
            JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            
            
            //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_tmp_bckg, "bckg 1","DOUBLE");//show first image

            //JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
            
            
        }
        
        

    }
    
    
    private void computeFisher(){
        
        
        float h=(float).0005;
        float h_photon=(float).0005;
        
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        
        
        
        
        
        //model with all psf
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelAndLikelihood(MyCudaStream.getCUstream(streamId), wh,width,height,sizePatch,sizePatch,numberToCompute,device_likelihood,device_model,device_image,device_PSF,device_parameters,index_X,index_Y,index_photon,device_background);//6 and 7 are index of X/Y patch starting
        MyVecDouble.cropFromImage(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch,sizePatch,numberToCompute,device_crop,width,height,device_model,device_parameters,index_X,index_Y);//6 and 7 are index of X/Y patch starting

        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop, "model","DOUBLE");//show first image
        
        

        //here, device_crop contains crops with all psfs
        //JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));IJ.log("fisher computed");JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        //MyVecDouble.subtractModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp_bckg,device_crop,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        
        
        //here, device_crop_tmp_bckg contains background + overlapping emitters only
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_tmp_bckg, "bckg","DOUBLE");//show first image
        
        
        
        //first: compute chi2:
        MyVecDouble.cropFromImageFloat(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch,sizePatch,numberToCompute,device_cropFloat,width,height,device_image,device_parameters,index_X,index_Y);//6 and 7 are index of X/Y patch starting
        MyVecDouble.chi2(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute,sizePatch*sizePatch,device_crop_tmp1,device_crop,device_cropFloat);//
        
        
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop, "device_crop","DOUBLE");//show first image
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_cropFloat, "device_cropFloat","FLOAT");//show first image
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_tmp1, "chi2","DOUBLE");//show first image
        
        
        
        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,sizePatch*sizePatch,numberToCompute,device_alpha,device_crop_tmp1,sizePatch*sizePatch,device_ones,1,device_beta,device_crop_tmp2,1);
        cudaMemcpyAsync(Pointer.to(score), device_crop_tmp2, numberToCompute*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
        //here: score1contains chi2 sum
        
        
        
        
        
        
        
        //partial photon
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,-h_photon,device_parameters);//6 and 7 are index of X/Y patch starting
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp1,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,+2*h_photon,device_parameters);//6 and 7 are index of X/Y patch starting
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp2,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        MyVecDouble.partialModel(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute,device_crop_partial_photon,device_crop_tmp2,device_crop_tmp1,h_photon);
        //put back
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_photon,-h_photon,device_parameters);//6 and 7 are index of X/Y patch starting
        
        
        
        
        
        
        //partial X
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_x,-h,device_parameters);//6 and 7 are index of X/Y patch starting
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp1,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_tmp1, "X0","DOUBLE");//show first image
        
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_x,2*h,device_parameters);//6 and 7 are index of X/Y patch starting
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp2,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_tmp2, "X1","DOUBLE");//show first image
                
        MyVecDouble.partialModel(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute,device_crop_partial_x,device_crop_tmp2,device_crop_tmp1,h);
        //put back
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_x,-h,device_parameters);
        
        
        
        
        
        //partial Y
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_y,-h,device_parameters);//6 and 7 are index of X/Y patch starting
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp1,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_y,2*h,device_parameters);//6 and 7 are index of X/Y patch starting
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp2,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        MyVecDouble.partialModel(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute,device_crop_partial_y,device_crop_tmp2,device_crop_tmp1,h);
        //put back
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_y,-h,device_parameters);
        
        
        
        
        
        
        
        
        //partial Z
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_z,-h,device_parameters);//6 and 7 are index of X/Y patch starting
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp1,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_z,2*h,device_parameters);//6 and 7 are index of X/Y patch starting
        dp.psf_fMany.computePSF(device_parameters,numberToCompute);
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        device_PSF=dp.psf_fMany.getPointerPSFfloat();
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        MyVecDouble.computeModelwithPhotonNumber(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute, numberToCompute,device_crop_tmp2,device_PSF,device_parameters,index_photon);//6 and 7 are index of X/Y patch starting
        
        MyVecDouble.partialModel(MyCudaStream.getCUstream(streamId),sizePatch*sizePatch*numberToCompute,device_crop_partial_z,device_crop_tmp2,device_crop_tmp1,h);
        //put back
        MyVecDouble.shiftParameter(MyCudaStream.getCUstream(streamId), numberToCompute,index_z,-h,device_parameters);
        
        
        
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_partial_photon, "ph","DOUBLE");//show first image
        
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_partial_x, "X","DOUBLE");//show first image
        
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_partial_y, "Y","DOUBLE");//show first image
        
        //dp.psf.imshow(sizePatch*sizePatch*numberToCompute,sizePatch, device_crop_partial_z, "Z","DOUBLE");//show first image
        
        MyVecDouble.computeFisherMatrix(MyCudaStream.getCUstream(streamId), sizePatch*sizePatch*numberToCompute*16,sizePatch*sizePatch,numberToCompute,device_fisher,device_crop,device_crop_partial_x,device_crop_partial_y,device_crop_partial_z,device_crop_partial_photon);
        
        cudaResult=jcuda.jcublas.JCublas2.cublasSgemv(handlecublas,CUBLAS_OP_T,sizePatch*sizePatch,16*numberToCompute,device_alpha_f,device_fisher,sizePatch*sizePatch,device_ones_f,1,device_beta_f,device_fisher_sum,1);
        
        //dp.psf.imshow(16*sizePatch*sizePatch*2,sizePatch, device_fisher, "device_fisher","FLOAT");
        //dp.psf.imshow(16*numberToCompute,4, device_fisher_sum, "device_fisher_sum","FLOAT");
        
        //this does not work
        //jcuda.jcublas.JCublas2.cublasSmatinvBatched(handlecublas,4,device_fisher_sum,4,device_fisher,  4,device_invInfo, 1);
        //dp.psf.imshow(16*numberToCompute,4, device_fisher, "INV","FLOAT");
        //dp.psf.imshow(numberToCompute,1, device_invInfo, "info","INT");
        
        
        cudaMemcpyAsync(Pointer.to(fisher), device_fisher_sum, numberToCompute*16*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cumemcpy cuda 0 "+cudaResult);}
        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));
        
        for (int i=0;i<numberToCompute;i++){
            
            for (int u=0;u<4;u++){
                for (int uu=0;uu<4;uu++){
                    fisher2[u][uu]=fisher[i*16+u*4+uu];
                }
            }
        
            Matrixe mat = new Matrixe(fisher2);
            try{
                mat=Matrixe.inverse(mat);
                fisher2=mat.getMatrixe();
                for (int u=0;u<4;u++){
                    std[u]=Math.sqrt(fisher2[u][u]);
                }
            }catch(Exception ee){//if not inversible:
                //dont take into account covar if non inversible
                for (int u=0;u<4;u++){
                    if (fisher2[u][u]!=0){
                        std[u]=Math.sqrt(1/fisher2[u][u]);
                    }
                    else{
                        std[u]=Double.MAX_VALUE;
                    }
                }
                
            }
            
            mat=null;
            
            
            this.position[i+index_crlb_x*numberToCompute]=(float)std[0];
            this.position[i+index_crlb_y*numberToCompute]=(float)std[1];
            this.position[i+index_crlb_z*numberToCompute]=(float)std[2];
            
            this.position[i+index_chi2*numberToCompute]=(float)score[i];
                

            
        }
        
        
        
        
        
        
        
        
        
        
    }
    
    
    
    
    
    
    
    public void free(){
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        JCuda.cudaFree(device_parameterSave) ; 
        JCuda.cudaFree(device_likelihoodResult3) ; 
        JCuda.cudaFree(device_likelihoodResult2) ; 
        JCuda.cudaFree(device_likelihoodResult1) ; 
        JCuda.cudaFree(device_likelihoodResult) ; 
        JCuda.cudaFree(device_beta_f) ; 
        JCuda.cudaFree(device_ones_f) ; 
        JCuda.cudaFree(device_alpha_f) ; 
        JCuda.cudaFree(device_invInfo) ; 
        JCuda.cudaFree(device_alpha) ; 
        JCuda.cudaFree(device_beta) ; 
        JCuda.cudaFree(device_ones) ; 
        JCuda.cudaFree(device_cropFloat) ; 
        JCuda.cudaFree(device_fisher_sum) ; 
        JCuda.cudaFree(device_fisher) ; 
        JCuda.cudaFree(device_crop_partial_photon) ; 
        JCuda.cudaFree(device_crop_partial_z) ; 
        JCuda.cudaFree(device_crop_partial_y) ; 
        JCuda.cudaFree(device_crop_partial_x) ; 
        JCuda.cudaFree(device_crop_tmp2) ; 
        JCuda.cudaFree(device_crop_tmp1) ; 
        JCuda.cudaFree(device_crop) ; 
        JCuda.cudaFree(device_crop_tmp_bckg) ; 
        JCuda.cudaFree(device_sum_gamma) ; 
        JCuda.cudaFree(device_gamma) ; 
        JCuda.cudaFree(device_parameters) ; 
        JCuda.cudaFree(device_background) ; 
        JCuda.cudaFree(device_model) ; 
        JCuda.cudaFree(device_likelihood) ; 
        JCuda.cudaFree(device_image) ; 
        
    }
    
    
}
