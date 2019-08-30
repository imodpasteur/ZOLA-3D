/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

/**
 *
 * @author benoit
 */

import ij.ImageStack;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.process.*;
import org.pasteur.imagej.utils.*;
import java.util.Arrays;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.JCublas;
import jcuda.jcublas.JCublas2;
import jcuda.jcublas.cublasHandle;
import jcuda.jcufft.JCufft;
import jcuda.jcufft.cufftHandle;
import jcuda.jcufft.cufftType;
import static jcuda.jcusparse.JCusparse.cusparseDgthr;
import static jcuda.jcusparse.JCusparse.cusparseDsctr;
import jcuda.jcusparse.cusparseHandle;
import static jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpy;
import jcuda.runtime.cudaError;
import jcuda.runtime.cudaMemcpyKind;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;


/**
 *
 * @author benoit
 */
public class GS_ {
    
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;

    CUstream custream;
    
    cufftHandle plan;
    
    //sizeImage*sizeFFT
    Pointer [] device_image;
    Pointer [] device_psf;
    Pointer [] device_real_stack;
    Pointer [] device_imag_stack;
    Pointer device_realSpacial;
    Pointer device_imagSpacial;
    Pointer device_tmpSpacial;
    
    Pointer device_sparseIndexEven;
    Pointer device_sparseIndexOdd;
    Pointer device_sparseIndexShift2D;
    
    //sizeFFT*sizeFFT*2
    Pointer device_fftdata;
    
    //disk size
    Pointer device_phase;
    Pointer device_pupil;
    Pointer [] device_phase_stack;
    Pointer [] device_pupil_stack;
    Pointer device_kx;
    Pointer device_ky;
    Pointer device_kz;
    Pointer device_kz_oil;
    Pointer device_kz_is_imaginary;
    Pointer device_kz_oil_is_imaginary;
    Pointer device_realDisk;
    Pointer device_imagDisk;
    Pointer device_tmpDisk;
    Pointer device_sparseIndexOddDisk;
    Pointer device_sparseIndexEvenDisk;
    
    
    
    
    
    
    
    
    
    PhaseParameters param;
    
    int sizeImage;
    int sizeDisk;
    
    double [][] pile;
    double [] outputDisk;
    public double [] deltaZ;
    
    
    
    int orderPupil=0;
    int axialside;
    public GS_(double [][][] stack,int sizeFFT,int order,double xystep,double zstep,double wavelength,double noil,double nwat,double zfocus,double na,double sigmaGaussianKernel,int zernikeCoefNumber,int axialside,boolean withApoFactor){
        param=new PhaseParameters(sizeFFT,sizeFFT,1,order,xystep,zstep,wavelength,noil,na,sigmaGaussianKernel,withApoFactor);
        this.sizeImage=sizeFFT;
        this.sizeDisk=param.sizeDisk;
        gs(stack);
    }
    
    
    public GS_(double [][][] stack,PhaseParameters param,int axialside){
        this.param=param;
        this.sizeImage=param.size;
        this.sizeDisk=param.sizeDisk;
        this.axialside=axialside;
        gs(stack);
    }
        
        
    private void gs(double [][][] stack){
        
        
        
        //param.weightZ=1.3;
        
        
        //IJ.log("stack length "+stack.length);
        double centerstack=stack.length/2;
        deltaZ=new double [stack.length];
        //IJ.log("center is : "+center+"  !!!!!!!!!!§§§§§§§§§§§§§§§!!!!!!!");
        for (int s=0;s<stack.length;s++){
            if (axialside==0){
                deltaZ[s]=(((double)s)-centerstack)*param.zstep;
            }
            else{
                deltaZ[s]=-(((double)s)-centerstack)*param.zstep;
            }
            
        }
        
        
        //stack=this.normalize(stack);
        
        //resCrossCorel=new double [stack.length-1][sizeImage][sizeImage];
        
        pile=new double [stack.length][sizeImage*sizeImage];
        
        outputDisk=new double[sizeDisk];
        
        
        //device in spacial domain
        device_image=new Pointer[stack.length];
        device_psf=new Pointer[stack.length];
        device_real_stack=new Pointer[stack.length];
        device_imag_stack=new Pointer[stack.length];
        
        device_realSpacial=new Pointer();
        device_imagSpacial=new Pointer();
        device_tmpSpacial=new Pointer();
        
        
        
        if (JCuda.cudaMalloc(device_realSpacial, sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        if (JCuda.cudaMalloc(device_imagSpacial, sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        if (JCuda.cudaMalloc(device_tmpSpacial, sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        
        
        for (int i=0;i<pile.length;i++){
            for (int ii=0;ii<stack[i].length;ii++){
                for (int iii=0;iii<stack[i][ii].length;iii++){
                    pile[i][(ii+sizeImage/2-stack[i].length/2)*sizeImage+(iii+sizeImage/2-stack[i][ii].length/2)]=stack[i][ii][iii];
                }
            }
            device_image[i]=new Pointer();
            device_psf[i]=new Pointer();
            device_real_stack[i]=new Pointer();
            device_imag_stack[i]=new Pointer();
            if (JCuda.cudaMalloc(device_image[i], sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
            if (JCuda.cudaMalloc(device_psf[i], sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
            if (JCuda.cudaMalloc(device_real_stack[i], sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
            if (JCuda.cudaMalloc(device_imag_stack[i], sizeImage*sizeImage * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
            
            
            cudaMemcpy(device_image[i], Pointer.to(pile[i]), sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
            JCuda.cudaMemset(device_real_stack[i], 0, sizeImage*sizeImage * Sizeof.DOUBLE);
            JCuda.cudaMemset(device_imag_stack[i], 0, sizeImage*sizeImage * Sizeof.DOUBLE);
        }
        
        //device in FFT domain
        
        
        
        
        
        double [] phase = new double[param.sizeDisk];
        double [] pupil = new double[param.sizeDisk];
        double [] kx = new double[param.sizeDisk];
        double [] ky = new double[param.sizeDisk];
        double [] kz = new double[param.sizeDisk];
        double [] kz_oil = new double[param.sizeDisk];
        double [] kz_is_imaginary = new double[param.sizeDisk];
        double [] kz_oil_is_imaginary = new double[param.sizeDisk];
        int [] sparseIndexOddDisk = new int[param.sizeDisk];
        int [] sparseIndexEvenDisk = new int[param.sizeDisk];
        
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
                        kz[id]=0;//2*Math.PI*Math.sqrt(Math.abs(left-right));
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
                        pupil[id]=1/(double)Math.sqrt(param.sizeDisk_cpu);//like that -> final sum=1 (unuseful actually because we normalize at the end)
                    }
                    
                    phase[id]=0;
                    
                    
                    
                    sparseIndexOddDisk[id]=((((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size))*2)+1;
                    sparseIndexEvenDisk[id]=((((i+param.size/2)%param.size)*param.size+((ii+param.size/2)%param.size))*2);
                    id++;
                }
                
                
                
            }
        }
        
        device_sparseIndexEvenDisk=new Pointer();
        device_sparseIndexOddDisk=new Pointer();
        
        
        
        if (JCuda.cudaMalloc(device_sparseIndexOddDisk, param.sizeDisk * Sizeof.INT) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_sparseIndexOddDisk, Pointer.to(sparseIndexOddDisk), param.sizeDisk*Sizeof.INT, cudaMemcpyHostToDevice);
        
        if (JCuda.cudaMalloc(device_sparseIndexEvenDisk, param.sizeDisk * Sizeof.INT) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_sparseIndexEvenDisk, Pointer.to(sparseIndexEvenDisk), param.sizeDisk*Sizeof.INT, cudaMemcpyHostToDevice);
        
        
        device_tmpDisk=new Pointer();
        if (JCuda.cudaMalloc(device_tmpDisk, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        
        device_realDisk=new Pointer();
        if (JCuda.cudaMalloc(device_realDisk, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        
        device_imagDisk=new Pointer();
        if (JCuda.cudaMalloc(device_imagDisk, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        
        device_kz=new Pointer();
        if (JCuda.cudaMalloc(device_kz, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_kz, Pointer.to(kz), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        device_kz_oil=new Pointer();
        if (JCuda.cudaMalloc(device_kz_oil, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_kz_oil, Pointer.to(kz_oil), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        device_kz_is_imaginary=new Pointer();
        if (JCuda.cudaMalloc(device_kz_is_imaginary, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_kz_is_imaginary, Pointer.to(kz_is_imaginary), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        device_kz_oil_is_imaginary=new Pointer();
        if (JCuda.cudaMalloc(device_kz_oil_is_imaginary, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_kz_oil_is_imaginary, Pointer.to(kz_oil_is_imaginary), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        
        
        device_kx=new Pointer();
        if (JCuda.cudaMalloc(device_kx, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_kx, Pointer.to(kx), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        
        device_ky=new Pointer();
        if (JCuda.cudaMalloc(device_ky, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_ky, Pointer.to(ky), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        
        device_phase=new Pointer();
        if (JCuda.cudaMalloc(device_phase, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_phase, Pointer.to(phase), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        
        
        device_pupil=new Pointer();
        if (JCuda.cudaMalloc(device_pupil, param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_pupil, Pointer.to(pupil), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        
        device_fftdata=new Pointer();
        if (JCuda.cudaMalloc(device_fftdata, sizeImage*sizeImage*2 * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        
        
        
        device_phase_stack=new Pointer[stack.length];
        device_pupil_stack=new Pointer[stack.length];
        for (int i=0;i<stack.length;i++){
            device_phase_stack[i]=new Pointer();
            if (JCuda.cudaMalloc(device_phase_stack[i], param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        
            device_pupil_stack[i]=new Pointer();
            if (JCuda.cudaMalloc(device_pupil_stack[i], param.sizeDisk * Sizeof.DOUBLE) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
            cudaMemcpy(device_pupil_stack[i], Pointer.to(pupil), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        }
        
        
        int [] sparseIndexOdd = new int[sizeImage*sizeImage];
        int [] sparseIndexEven = new int[sizeImage*sizeImage];
        int [] sparseIndexShift2D = new int[sizeImage*sizeImage];
        
        for (int i=0;i<sizeImage*sizeImage;i++){
            sparseIndexOdd[i]=(i*2)+1;
            sparseIndexEven[i]=(i*2);
        }
        for (int i=0;i<sizeImage;i++){
            for (int ii=0;ii<sizeImage;ii++){
                sparseIndexShift2D[i*sizeImage+ii]=(((i+sizeImage/2)%sizeImage)*sizeImage+((ii+sizeImage/2)%sizeImage));
            }
        }
        
        device_sparseIndexEven=new Pointer();
        if (JCuda.cudaMalloc(device_sparseIndexEven, sizeImage*sizeImage * Sizeof.INT) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_sparseIndexEven, Pointer.to(sparseIndexEven), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice);
        
        
        device_sparseIndexOdd=new Pointer();
        if (JCuda.cudaMalloc(device_sparseIndexOdd, sizeImage*sizeImage * Sizeof.INT) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_sparseIndexOdd, Pointer.to(sparseIndexOdd), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice);
        
        device_sparseIndexShift2D=new Pointer();
        if (JCuda.cudaMalloc(device_sparseIndexShift2D, sizeImage*sizeImage * Sizeof.INT) != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda gaussian kernel");return ;}
        cudaMemcpy(device_sparseIndexShift2D, Pointer.to(sparseIndexShift2D), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice);
        
        
        plan = new cufftHandle();
        JCufft.cufftPlan2d(plan, sizeImage,sizeImage, cufftType.CUFFT_Z2Z);
        JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(param.stream));
        
        handlecublas=MyCudaStream.getHandleCublas(param.stream);
        handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
        custream=MyCudaStream.getCUstream(param.stream);
        
        
        
        
        
        
        init();
        
        
    }
    
    
    
    
    
    
    
    public void run(int iter){
        for (int i=0;i<iter;i++){
            
            this.computeFFT();//
            this.computeMeanPhase();
            //this.computeMeanPupil();
            //this.computeMeanPupilTmp();
            this.resetPupil();
            /*if (i<iter/4.){
                this.setPupilToRing();
            }
            else {
                this.computeMeanPupil();
            }*/
            
            
            /*else if(i<3.*iter/4.){
                this.computeMeanPupil();
                this.computeMeanPupilFitGaussian();
            }
            else{
                this.computeMeanPupil();
                this.computeMeanPupilFitGaussian();
            }*/
            //setPupilToRing();
            //
            //ImageShow.imshow(this.getPupilStack(),"pupil");
            
            this.computeInvFFT();
            
            
        }
        //ImageShow.imshow(this.getPhaseStack(),"phaseFFT");
        //ImageShow.imshow(this.getPupilStack(),"pupil");
        
        
        //ImageShow.imshow(this.getPupil(),"pupil");
        //ImageShow.imshow(this.getPSF(),"psf");   
        //getPhaseZernikeFitCoefs(15);
    }
    
    
    private void init(){
        //JCuda.cudaMemset(device_pupil, 0, sizeDisk * Sizeof.DOUBLE);
        JCuda.cudaMemset(device_phase, 0, sizeDisk * Sizeof.DOUBLE);
        //MyVecDouble.addScalar(custream,sizeDisk, device_pupil, device_pupil, 1/(double)Math.sqrt(param.sizeDisk));
        
        for (int i=0;i<device_image.length;i++){
            JCuda.cudaMemset(device_fftdata, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE);
            JCuda.cudaMemset(device_phase_stack[i], 0, sizeDisk * Sizeof.DOUBLE);
            
            
            MyVecDouble.mulScalar(custream,sizeDisk, device_tmpDisk, device_kz_oil, +this.deltaZ[i]);
            MyVecDouble.add(custream,sizeDisk, this.device_phase_stack[i], device_tmpDisk, this.device_phase_stack[i]);
            
            
            cudaMemcpy(device_pupil_stack[i],device_pupil,  sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
            //compute abs of magnitude
            MyVecDouble.fabs(custream,sizeDisk, device_tmpDisk, device_pupil_stack[i]);
            //this.imshow(device_tmpDisk,"device_tmpDisk");
            //compute cos and sin
            MyVecDouble.cos(custream,sizeDisk, device_realDisk, device_phase_stack[i]);
            MyVecDouble.sin(custream,sizeDisk, device_imagDisk, device_phase_stack[i]);
            
            MyVecDouble.mul(custream,sizeDisk, device_realDisk,device_realDisk, device_tmpDisk);
            MyVecDouble.mul(custream,sizeDisk, device_imagDisk,device_imagDisk, device_tmpDisk);
            
            
            
            cusparseDsctr(handlecusparse, param.sizeDisk, device_realDisk, device_sparseIndexEvenDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDsctr(handlecusparse, param.sizeDisk, device_imagDisk, device_sparseIndexOddDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);

            JCufft.cufftExecZ2Z(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);
            
            
            MyVecDouble.divScalar(custream,sizeImage*sizeImage*2, device_fftdata, device_fftdata, ((double)sizeImage));
        
        
            //gather complex
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_fftdata, device_realSpacial,device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_fftdata, device_imagSpacial,device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);

            //this.imshow(sizeImage*sizeImage,sizeImage,device_realSpacial,"device_realSpacial after","DOUBLE");

            
            //this.imshow(sizeImage*sizeImage,sizeImage,device_imagSpacial,"device_imagSpacial after","DOUBLE");
            
            
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_realSpacial, device_real_stack[i],device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_imagSpacial, device_imag_stack[i],device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);

            
            
            

            MyVecDouble.mul(custream,sizeImage*sizeImage, device_realSpacial, device_real_stack[i],device_real_stack[i]);//square
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imagSpacial, device_imag_stack[i],device_imag_stack[i]);//square
            MyVecDouble.add(custream,sizeImage*sizeImage, device_psf[i], device_realSpacial, device_imagSpacial);
            
            double [] imag = new double [sizeImage*sizeImage];
            cudaMemcpy(Pointer.to(imag), device_imag_stack[i], sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            double [] real = new double [sizeImage*sizeImage];
            cudaMemcpy(Pointer.to(real), device_real_stack[i], sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeImage*sizeImage;ii++){
                imag[ii]=Math.atan2(imag[ii], real[ii]);
            }
            cudaMemcpy(device_imagSpacial,Pointer.to(imag),  sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
            
            //do not manage NaN values
            //VecDouble.div(sizeImage*sizeImage, this.device_imagSpacial, device_imag_stack[i], device_real_stack[i]);
            //VecDouble.atan(sizeImage*sizeImage, this.device_imagSpacial, this.device_imagSpacial);
            
            
            //magnitude replacement
            
            MyVecDouble.cos(custream,sizeImage*sizeImage, device_tmpSpacial, device_imagSpacial);//cosin
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_real_stack[i], device_image[i], device_tmpSpacial);
            
            MyVecDouble.sin(custream,sizeImage*sizeImage, device_tmpSpacial, device_imagSpacial);//sin
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag_stack[i], device_image[i], device_tmpSpacial);
            
            
        }
        
    }
    
    
    
    private void computeFFT(){
        double [] imag = new double [sizeDisk];
        double [] real = new double [sizeDisk];
        for (int i=0;i<device_image.length;i++){
            //IJ.log("set magn and phase here, not real and imag part");
            
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_real_stack[i], device_realSpacial,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_imag_stack[i], device_imagSpacial,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
        
            
            cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_realSpacial, device_sparseIndexEven,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);
            
            cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_imagSpacial, device_sparseIndexOdd,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);
            
            
            JCufft.cufftExecZ2Z(plan, device_fftdata,device_fftdata,JCufft.CUFFT_FORWARD);
            
            MyVecDouble.divScalar(custream,2*sizeImage*sizeImage, device_fftdata, device_fftdata, sizeImage);
            
            
            /*
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_fftdata,this.device_realSpacial, device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_fftdata,this.device_imagSpacial, device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);
            
            //this.imshow(sizeImage * sizeImage, sizeImage, device_realSpacial, "real after", "DOUBLE");

            //this.imshow(sizeImage * sizeImage, sizeImage, device_imagSpacial, "imag after", "DOUBLE");
            
            MyVecDouble.div(custream,sizeImage*sizeImage, device_tmpSpacial, device_imagSpacial, device_realSpacial);
            MyVecDouble.atan(custream,sizeImage*sizeImage, device_tmpSpacial, device_tmpSpacial);
            */
            
            
            //this.imshow(sizeImage*sizeImage*2,sizeImage*2,device_fftdata, "fft","DOUBLE");
            
            cusparseDgthr(handlecusparse, param.sizeDisk, device_fftdata,device_realDisk, device_sparseIndexEvenDisk, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, param.sizeDisk, device_fftdata,device_imagDisk, device_sparseIndexOddDisk, CUSPARSE_INDEX_BASE_ZERO);
            
            //compute pupil
            MyVecDouble.mul(custream,sizeDisk,this.device_pupil_stack[i],device_realDisk,device_realDisk);
            MyVecDouble.mul(custream,sizeDisk,device_tmpDisk,device_imagDisk,device_imagDisk);
            MyVecDouble.add(custream,sizeDisk,device_tmpDisk,this.device_pupil_stack[i],device_tmpDisk);
            MyVecDouble.sqrt(custream,sizeDisk,this.device_pupil_stack[i],device_tmpDisk);
            
            
            //this.imshow(sizeImage*sizeImage,sizeImage,device_imagSpacial, "dev input","DOUBLE");
            //this.imshow(device_realDisk, "real");
            
            //compute atan:
            
            
            cudaMemcpy(Pointer.to(imag), device_imagDisk, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            
            cudaMemcpy(Pointer.to(real), device_realDisk, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeDisk;ii++){
                imag[ii]=Math.atan2(imag[ii], real[ii]);
            }
            cudaMemcpy(device_phase_stack[i],Pointer.to(imag),  sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
            
            //this.imshow(sizeImage * sizeImage, sizeImage, device_phase_stack[i], "device_phase_stack[i]", "DOUBLE");
            
            //VecDouble.div(sizeDisk, this.device_phase_stack[i], device_imagDisk, device_realDisk);
            //VecDouble.atan(sizeDisk, this.device_phase_stack[i], this.device_phase_stack[i]);
            
            //this.imshow(this.device_phase_stack[i], "phase"+i);
        
            
        }
        
        
        
    }
    
    
    private void resetPupil(){
        
        for (int i=0;i<device_image.length;i++){
            cudaMemcpy(device_pupil_stack[i],device_pupil,  sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
        }
    }
    
    
    private void setPupilToRing(){
        JCuda.cudaMemset(device_pupil, 0, sizeDisk * Sizeof.DOUBLE);
        MyVecDouble.addScalar(custream,sizeDisk, device_pupil, device_pupil, 1/(double)Math.sqrt(param.sizeDisk));
        for (int i=0;i<device_image.length;i++){
            cudaMemcpy(device_pupil_stack[i],device_pupil,  sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
        }
    }
    
    
    
    private void computeMeanPupil(){
        JCuda.cudaMemset(device_pupil, 0, sizeDisk * Sizeof.DOUBLE);
        for (int i=0;i<device_image.length;i++){
            MyVecDouble.add(custream,sizeDisk, device_pupil, device_pupil, device_pupil_stack[i]);
        }
        //normalization sum=1
        MyVecDouble.mul(custream,sizeDisk, device_pupil, device_pupil, device_pupil);
        
        //double som=JCublas.cublasDasum(sizeDisk,device_pupil,1);//norm sum=1
        
        Pointer p = new Pointer();
        JCuda.cudaMalloc(p, 1 * Sizeof.DOUBLE);
        double [] som = new double [1];
        JCublas2.cublasDasum(handlecublas,sizeDisk,device_pupil,1,p);//norm sum=1
        JCuda.cudaMemcpy(Pointer.to(som), p, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);
        
        MyVecDouble.divScalar(custream,sizeDisk, device_pupil, device_pupil,som[0]);//square
        MyVecDouble.sqrt(custream,sizeDisk, device_pupil, device_pupil);
        for (int i=0;i<device_image.length;i++){
            cudaMemcpy(device_pupil_stack[i], device_pupil, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
        
        }
        JCuda.cudaFree(p);
    }
    
    private void computeMeanPupilTmp(){
        int num=10;
        for (int i=0;i<device_image.length;i+=num){
            JCuda.cudaMemset(device_pupil, 0, sizeDisk * Sizeof.DOUBLE);
            for (int ii=i;ii<i+num;ii++){
                if (ii<device_image.length){
                    MyVecDouble.add(custream,sizeDisk, device_pupil, device_pupil, device_pupil_stack[ii]);
                }
            }
            MyVecDouble.mul(custream,sizeDisk, device_pupil, device_pupil, device_pupil);
            double som=JCublas.cublasDasum(sizeDisk,device_pupil,1);//norm sum=1
            MyVecDouble.divScalar(custream,sizeDisk, device_pupil, device_pupil,som);//square
            MyVecDouble.sqrt(custream,sizeDisk, device_pupil, device_pupil);
            
            for (int ii=i;ii<i+num;ii++){
                if (ii<device_image.length){
                    cudaMemcpy(device_pupil_stack[ii], device_pupil, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
                }
            }
        }
        
        
    }
    
    
    
    
    
    /*
    private void computeMeanPhase(){
        
        double [] kz = new double [sizeDisk];
        double [] ph = new double [sizeDisk];
        double [] mc = new double [sizeDisk];
        double [] ms = new double [sizeDisk];
        for (int ii=0;ii<sizeDisk;ii++){
            mc[ii]=0;
                    ms[ii]=0;
        }
        cudaMemcpy(Pointer.to(kz), device_kz, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
        for (int i=0;i<device_image.length;i++){
        //{int i=device_image.length/2;
            //unpropagation
            
            cudaMemcpy(Pointer.to(ph), device_phase_stack[i], sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeDisk;ii++){
                
                    ph[ii]-=this.deltaZ[i]*param.weightZ*kz[ii];
                    mc[ii]+=Math.cos(ph[ii]);
                    ms[ii]+=Math.sin(ph[ii]);
                
            }
            //this.imshow(sizeImage, mc, "cos"+i);
                //this.imshow(sizeImage, ms, "sin"+i);
        }
        for (int ii=0;ii<sizeDisk;ii++){
                
                    
                    mc[ii]/=(double)device_image.length;
                    ms[ii]/=(double)device_image.length;
                ph[ii]=Math.signum(ms[ii])*Math.acos(mc[ii]);
            }
        
        //here, the mean phase is in device_cosphase
        
        
        //propagation
        
        for (int i=0;i<device_image.length;i++){
        //{int i=3;
            //IJ.log("this.deltaZ[i]   "+this.deltaZ[i]);
            for (int ii=0;ii<sizeDisk;ii++){
                mc[ii]=ph[ii]+this.deltaZ[i]*param.weightZ*kz[ii];
                //mc[ii]=ph[ii];
            }
            cudaMemcpy(this.device_phase_stack[i], Pointer.to(mc), sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
            
            
        }
        
    }*/
    
    
    
    private void computeMeanPhase(){
        JCuda.cudaMemset(device_realDisk, 0, sizeDisk * Sizeof.DOUBLE);//cosin
        JCuda.cudaMemset(device_imagDisk, 0, sizeDisk * Sizeof.DOUBLE);//sin
        for (int i=0;i<device_image.length;i++){
            //unpropagation
            MyVecDouble.mulScalar(custream,sizeDisk, device_phase, device_kz_oil, +this.deltaZ[i]*1);
            MyVecDouble.add(custream,sizeDisk, device_phase, device_phase, this.device_phase_stack[i]);
            
            MyVecDouble.cos(custream,sizeDisk, device_tmpDisk, device_phase);//cosin
            //this.imshow(device_tmpDisk, "cos"+i);
            MyVecDouble.add(custream,sizeDisk, device_realDisk, device_realDisk, device_tmpDisk);
            
            MyVecDouble.sin(custream,sizeDisk, device_tmpDisk, device_phase);//sin
            MyVecDouble.add(custream,sizeDisk, device_imagDisk, device_imagDisk, device_tmpDisk);
        }
        //normalization
        MyVecDouble.divScalar(custream,sizeDisk, device_realDisk, device_realDisk, device_image.length);
        MyVecDouble.divScalar(custream,sizeDisk, device_imagDisk, device_imagDisk, device_image.length);
        
        
        //compute sign of sin function
        MyVecDouble.gteScalar(custream,sizeDisk, device_tmpDisk, device_imagDisk, 0);
        MyVecDouble.subScalar(custream,sizeDisk, device_tmpDisk, device_tmpDisk, .5);
        MyVecDouble.mulScalar(custream,sizeDisk, device_tmpDisk, device_tmpDisk, 2.);
        //compute arc cos
        MyVecDouble.acos(custream,sizeDisk, device_realDisk, device_realDisk);
        MyVecDouble.mul(custream,sizeDisk, device_phase, device_tmpDisk, device_realDisk);
        //here, the mean phase is in device_phase
        //this.imshow(device_phase, "phase");
        //propagation
        for (int i=0;i<device_image.length;i++){
            MyVecDouble.mulScalar(custream,sizeDisk, device_tmpDisk, device_kz_oil, -this.deltaZ[i]*1);
            MyVecDouble.add(custream,sizeDisk, this.device_phase_stack[i], device_tmpDisk, device_phase);
            
            
            
            
            
        }
    }
    
    
    
    
    private void computeInvFFT(){
        double [] imag = new double [sizeImage*sizeImage];
        double [] real = new double [sizeImage*sizeImage];
        for (int i=0;i<device_image.length;i++){
            JCuda.cudaMemset(device_fftdata, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE);
            //compute abs of magnitude
            MyVecDouble.fabs(custream,sizeDisk, device_tmpDisk, device_pupil_stack[i]);
            
            //compute cos and sin
            MyVecDouble.cos(custream,sizeDisk, device_realDisk, device_phase_stack[i]);
            MyVecDouble.sin(custream,sizeDisk, device_imagDisk, device_phase_stack[i]);
            
            MyVecDouble.mul(custream,sizeDisk, device_realDisk,device_realDisk, device_tmpDisk);
            MyVecDouble.mul(custream,sizeDisk, device_imagDisk,device_imagDisk, device_tmpDisk);
            
            cusparseDsctr(handlecusparse, param.sizeDisk, device_realDisk, device_sparseIndexEvenDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDsctr(handlecusparse, param.sizeDisk, device_imagDisk, device_sparseIndexOddDisk,device_fftdata, CUSPARSE_INDEX_BASE_ZERO);

            //this.imshow(sizeImage*sizeImage*2,sizeImage*2,device_fftdata, "fft","DOUBLE");
            JCufft.cufftExecZ2Z(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);
            
            
            
            
            MyVecDouble.divScalar(custream,sizeImage*sizeImage*2, device_fftdata, device_fftdata, sizeImage);
        
        
            //gather complex
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_fftdata, device_realSpacial,device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_fftdata, device_imagSpacial,device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);



            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_realSpacial, device_real_stack[i],device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_imagSpacial, device_imag_stack[i],device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);

            
            

            MyVecDouble.mul(custream,sizeImage*sizeImage, device_realSpacial, device_real_stack[i],device_real_stack[i]);//square
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imagSpacial, device_imag_stack[i],device_imag_stack[i]);//square
            MyVecDouble.add(custream,sizeImage*sizeImage, device_psf[i], device_realSpacial, device_imagSpacial);
            
            
            
            cudaMemcpy(Pointer.to(imag), device_imag_stack[i], sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            
            cudaMemcpy(Pointer.to(real), device_real_stack[i], sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeImage*sizeImage;ii++){
                imag[ii]=Math.atan2(imag[ii], real[ii]);
            }
            cudaMemcpy(device_imagSpacial,Pointer.to(imag),  sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
            
            //do not manage NaN
            //VecDouble.div(sizeImage*sizeImage, this.device_imagSpacial, device_imag_stack[i], device_real_stack[i]);
            //VecDouble.atan(sizeImage*sizeImage, this.device_imagSpacial, this.device_imagSpacial);
            
            //replacement with magnitude of image input
            
            MyVecDouble.cos(custream,sizeImage*sizeImage, device_tmpSpacial, device_imagSpacial);//cosin
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_real_stack[i], device_image[i], device_tmpSpacial);
            
            MyVecDouble.sin(custream,sizeImage*sizeImage, device_tmpSpacial, device_imagSpacial);//sin
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag_stack[i], device_image[i], device_tmpSpacial);
            
        }
        
        
    }
    
    
    double [][][] getPSF(){
        double [][][] stackres =new double [this.device_image.length][sizeImage][sizeImage];
        for (int i=0;i<this.device_image.length;i++){
            cudaMemcpy(Pointer.to(pile[i]), device_psf[i], sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeImage;ii++){
                for (int iii=0;iii<sizeImage;iii++){
                    stackres[i][ii][iii]=pile[i][ii*sizeImage+iii];
                }
            }
        }
        return stackres;
    }
    
    
    double [][][] getImage(){
        double [][][] stackres =new double [this.device_image.length][sizeImage][sizeImage];
        for (int i=0;i<this.device_image.length;i++){
            cudaMemcpy(Pointer.to(pile[i]), device_image[i], sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeImage;ii++){
                for (int iii=0;iii<sizeImage;iii++){
                    stackres[i][ii][iii]=pile[i][ii*sizeImage+iii];
                }
            }
        }
        return stackres;
    }
    
    double [][][] getPhaseStack(){
        double [][][] stackres =new double [this.device_phase_stack.length][sizeImage][sizeImage];
        for (int i=0;i<this.device_phase_stack.length;i++){
            cudaMemcpy(Pointer.to(outputDisk), device_phase_stack[i], sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeDisk;ii++){
                stackres[i][param.disk2D[ii][0]][param.disk2D[ii][1]]=outputDisk[ii];
                
            }
        }
        return stackres;
    }
    
    
    double [][][] getPupilStack(){
        double [][][] stackres =new double [this.device_pupil_stack.length][sizeImage][sizeImage];
        for (int i=0;i<this.device_pupil_stack.length;i++){
            cudaMemcpy(Pointer.to(outputDisk), device_pupil_stack[i], sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
            for (int ii=0;ii<sizeDisk;ii++){
                stackres[i][param.disk2D[ii][0]][param.disk2D[ii][1]]=outputDisk[ii];
                
            }
        }
        return stackres;
    }
    
    
    public double [][] getPhase(){
        double [][] res =new double [sizeImage][sizeImage];
        cudaMemcpy(Pointer.to(outputDisk), device_phase, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
        for (int ii=0;ii<sizeDisk;ii++){
            res[param.disk2D[ii][0]][param.disk2D[ii][1]]=outputDisk[ii];
        }
        return res;
    }
    
    
    
    
    
    public double [][] getPupil(){
        double [][] res =new double [sizeImage][sizeImage];
        cudaMemcpy(Pointer.to(outputDisk), device_pupil, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
        for (int ii=0;ii<sizeDisk;ii++){
            res[param.disk2D[ii][0]][param.disk2D[ii][1]]=outputDisk[ii];
        }
        return res;
    }
    
    
    
    
    
    public double [] getPhaseZernikeFitCoefs(int coefNumber){
        
        Zernike z = new Zernike(this.sizeImage, param.sizeRadiusRingPixel,coefNumber);
        double [][] myPhase=new double [sizeImage][sizeImage];
        cudaMemcpy(Pointer.to(outputDisk), device_phase, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
        for (int ii=0;ii<sizeDisk;ii++){
            myPhase[param.disk2D[ii][0]][param.disk2D[ii][1]]=outputDisk[ii];
        }
        
        ImageShow.imshow(myPhase,"phaseGS");
        
        double [][] myPupil=new double [sizeImage][sizeImage];
        cudaMemcpy(Pointer.to(outputDisk), device_pupil, sizeDisk*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
        for (int ii=0;ii<sizeDisk;ii++){
            myPupil[param.disk2D[ii][0]][param.disk2D[ii][1]]=outputDisk[ii];
        }
        
        double [] w=z.fitCosSin(myPhase, myPupil);
        
        return w;
    }
    
    
    
    
    
    
    
    
    
    
    
    public double [][][] normalize(double [][][] image){
        double tailleFiltre=1;
        int decX=(sizeImage-image[0].length)/2;
        int decY=(sizeImage-image[0].length)/2;
        GaussianKernel_ gk = new GaussianKernel_(this.sizeImage,this.sizeImage,tailleFiltre,0);
        double [][][] res=new double [image.length][sizeImage][sizeImage];
        double [][][] tmp=new double [image.length][sizeImage][sizeImage];
        boolean [][][] isPSF=new boolean [image.length][sizeImage][sizeImage];
        double [][][] tmpFilt=new double [image.length][sizeImage][sizeImage];
        double [] photonNumber=new double [image.length];
        double [] photonNumberAfterFilt=new double [image.length];
        for (int z=0;z<image.length;z++){
            photonNumber[z]=0;
            //anscomb transform
            for (int i=0;i<image[z].length;i++){
                for (int ii=0;ii<image[z][0].length;ii++){
                    photonNumber[z]+=image[z][i][ii];
                    tmp[z][i+decX][ii+decY]=2.*Math.sqrt(image[z][i][ii]+3./8.);
                }
            }
            
            gk.setImage(tmp[z]);
            //gk.filter();
            gk.getImage(tmpFilt[z]);
            

            
            
        }
        
        
        double meanback=0;
        double count=0;
        for (int i=0;i<image.length;i++){

            //first edge
            for (int ii=0;ii<image[0].length;ii+=image[0].length-1){
                for (int iii=0;iii<image[0][0].length;iii++){
                    meanback+=tmpFilt[i][ii+decX][iii+decY];
                    count++;
                }
            }
            for (int ii=1;ii<image[0].length-1;ii++){
                for (int iii=0;iii<image[0][0].length;iii+=image[0][0].length-1){
                    meanback+=tmpFilt[i][ii+decX][iii+decY];
                    count++;
                }
            }
            //second edge
            for (int ii=1;ii<image[0].length-1;ii++){
                for (int iii=1;iii<image[0][0].length-1;iii+=image[0][0].length-3){
                    meanback+=tmpFilt[i][ii+decX][iii+decY];
                    count++;
                }
            }
            for (int ii=1;ii<image[0].length;ii+=image[0].length-3){
                for (int iii=2;iii<image[0][0].length;iii++){
                    meanback+=tmpFilt[i][ii+decX][iii+decY];
                    count++;
                }
            }
        }
        meanback/=count;
        
        double std=0;
        
        for (int i=0;i<image.length;i++){
            for (int ii=0;ii<image[0].length;ii+=image[0].length-1){
                for (int iii=0;iii<image[0][0].length;iii++){
                    std+=(tmpFilt[i][ii+decX][iii+decY]-meanback)*(tmpFilt[i][ii+decX][iii+decY]-meanback);
                }
            }
            for (int ii=1;ii<image[0].length-1;ii++){
                for (int iii=0;iii<image[0][0].length;iii+=image[0][0].length-1){
                    std+=(tmpFilt[i][ii+decX][iii+decY]-meanback)*(tmpFilt[i][ii+decX][iii+decY]-meanback);
                }
            }
            for (int ii=1;ii<image[0].length-1;ii++){
                for (int iii=1;iii<image[0][0].length-1;iii+=image[0][0].length-3){
                    std+=(tmpFilt[i][ii+decX][iii+decY]-meanback)*(tmpFilt[i][ii+decX][iii+decY]-meanback);
                }
            }
            for (int ii=1;ii<image[0].length;ii+=image[0].length-3){
                for (int iii=2;iii<image[0][0].length;iii++){
                    std+=(tmpFilt[i][ii+decX][iii+decY]-meanback)*(tmpFilt[i][ii+decX][iii+decY]-meanback);
                }
            }
        }
        std/=count;
        std=Math.sqrt(std);
        
        meanback=meanback+std*2.355/2;
        
        double [] vect4errode=new double [(((int)Math.ceil(tailleFiltre))*2+1)*(((int)Math.ceil(tailleFiltre))*2+1)];
        //IJ.log("meanBack + stdBack  "+meanback);
        
        
        //ImageShow.imshow(tmpFilt,"im filtered");
        
        for (int i=0;i<tmpFilt.length;i++){
            //IJ.log("mean back B:"+meanback+"   "+std+"   mean+5std "+(meanback+5*std)+"   "+count);
            for (int ii=0;ii<tmpFilt[0].length;ii++){
                for (int iii=0;iii<tmpFilt[0][0].length;iii++){
                    tmp[i][ii][iii]=-1;
                    isPSF[i][ii][iii]=false;
                    tmpFilt[i][ii][iii]-=meanback;
                    if (tmpFilt[i][ii][iii]<=0){
                        tmpFilt[i][ii][iii]=-1;
                    }
                }
            }
            //set borders to 0
            for (int ii=0;ii<sizeImage;ii+=sizeImage-1){
                for (int iii=0;iii<image[0][0].length;iii++){
                    tmpFilt[i][ii][iii]=-1;
                }
            }
            for (int ii=1;ii<sizeImage-1;ii++){
                for (int iii=0;iii<sizeImage;iii+=sizeImage-1){
                    tmpFilt[i][ii][iii]=-1;
                }
            }
            for (int ii=1;ii<sizeImage-1;ii++){
                for (int iii=1;iii<sizeImage-1;iii+=sizeImage-3){
                    tmpFilt[i][ii][iii]=-1;
                }
            }
            for (int ii=1;ii<sizeImage;ii+=sizeImage-3){
                for (int iii=2;iii<sizeImage;iii++){
                    tmpFilt[i][ii][iii]=-1;
                }
            }
        
            //IJ.log("mean back B:"+meanback+"   "+std+"   mean+5std "+(meanback+5*std)+"   "+count);
            for (int ii=2;ii<sizeImage-2;ii++){
                for (int iii=2;iii<sizeImage-2;iii++){
                    if (tmpFilt[i][ii][iii]>0){
                        for (int a=-((int)Math.ceil(tailleFiltre)),k=0;a<=((int)Math.ceil(tailleFiltre));a++){
                            for (int aa=-((int)Math.ceil(tailleFiltre));aa<=((int)Math.ceil(tailleFiltre));aa++,k++){
                                vect4errode[k]=tmpFilt[i][ii+a][iii+aa];
                            }
                        }
                        Arrays.sort(vect4errode);
                        tmp[i][ii][iii]=vect4errode[0];
                        
                    }
                    
                    if(tmp[i][ii][iii]>0){
                        isPSF[i][ii][iii]=true;
                    }
                    
                }
            }
        }
        
        
        
        
        
        
        gk.free();
        
        
        double newMeanBcg=0;
        double countnew=0;
        //compute new mean background
        for (int z=0;z<image.length;z++){
            //anscomb transform
            for (int i=0;i<image[z].length;i++){
                for (int ii=0;ii<image[z][0].length;ii++){
                    
                    res[z][i+decX][ii+decY]=Math.sqrt(image[z][i][ii]);
                    
                    if (!isPSF[z][i+decX][ii+decY]){
                        newMeanBcg+=res[z][i+decX][ii+decY];
                        countnew++;
                    }
                    
                }
            }
        }
        double newstd=0;
        newMeanBcg/=countnew;
        for (int z=0;z<image.length;z++){
            //anscomb transform
            for (int i=0;i<image[z].length;i++){
                for (int ii=0;ii<image[z][0].length;ii++){
                    if (!isPSF[z][i+decX][ii+decY]){
                        newstd+=(res[z][i+decX][ii+decY]-newMeanBcg)*(res[z][i+decX][ii+decY]-newMeanBcg);
                        
                    }
                    
                }
            }
        }
        newstd/=countnew;
        newstd=Math.sqrt(newstd)*2.355/2;
        //IJ.log("new mean : "+newMeanBcg+"  new std "+newstd);
        for (int z=0;z<image.length;z++){
            photonNumberAfterFilt[z]=0;
            //anscomb transform
            for (int i=0;i<sizeImage;i++){
                for (int ii=0;ii<sizeImage;ii++){
                    
                    if (isPSF[z][i][ii]){
                        res[z][i][ii]-=(newMeanBcg+newstd);
                        photonNumberAfterFilt[z]+=res[z][i][ii];
                    }
                    else{
                        res[z][i][ii]=0;
                    }
                    if (res[z][i][ii]<0){
                        res[z][i][ii]=0;
                    }
                    
                }
            }
            //IJ.write(""+photonNumber[z]+"  "+photonNumberAfterFilt[z]);
            
        }
        
        ImageShow.imshow(res,"res");
        
        
        return res;
    }
    
    
    public double [][][] normalizeOld(double [][][] image){
        
        double [][][] im=new double[image.length][image[0].length][image[0][0].length];
        
        
        im = new double [image.length][this.sizeImage][this.sizeImage];
        for (int k=0;k<image.length;k++){
            
            for (int i=0;i<image[0].length;i++){
                for (int ii=0;ii<image[0][0].length;ii++){
                    if ((sizeImage/2-image[0].length/2+i>=0)&&(sizeImage/2-image[0].length/2+i<sizeImage)&&(sizeImage/2-image[0][0].length/2+i>=0)&&(sizeImage/2-image[0][0].length/2+i<sizeImage)){
                        im[k][sizeImage/2-image[0].length/2+i][sizeImage/2-image[0][0].length/2+ii]=image[k][i][ii];
                    }
                }
            }
        }
        
        double meanback=0;
        double count=0;
        for (int i=0;i<image.length;i++){
            
            
            for (int ii=0;ii<image[0].length;ii+=image[0].length-1){
                for (int iii=0;iii<image[0][0].length;iii++){
                    meanback+=image[i][ii][iii];
                    count++;
                }
            }
            for (int ii=0;ii<image[0].length;ii++){
                for (int iii=0;iii<image[0][0].length;iii+=image[0][0].length-1){
                    meanback+=image[i][ii][iii];
                    count++;
                }
            }
            for (int ii=1;ii<image[0].length;ii+=image[0].length-2){
                for (int iii=1;iii<image[0][0].length-1;iii+=image[0][0].length-1){
                    meanback+=image[i][ii][iii];
                    count++;
                }
            }
            for (int ii=1;ii<image[0].length-1;ii++){
                for (int iii=1;iii<image[0][0].length;iii+=image[0][0].length-2){
                    meanback+=image[i][ii][iii];
                    count++;
                }
            }
        }
        meanback/=count;
        double std=0;
        
        for (int i=0;i<image.length;i++){
            for (int ii=0;ii<image[0].length;ii+=image[0].length-1){
                for (int iii=0;iii<image[0][0].length;iii+=image[0][0].length-1){
                    std+=(image[i][ii][iii]-meanback)*(image[i][ii][iii]-meanback);
                }
            }
            for (int ii=0;ii<image[0].length;ii++){
                for (int iii=0;iii<image[0][0].length;iii+=image[0][0].length-1){
                    std+=(image[i][ii][iii]-meanback)*(image[i][ii][iii]-meanback);
                }
            }
            for (int ii=1;ii<image[0].length;ii+=image[0].length-2){
                for (int iii=1;iii<image[0][0].length-1;iii+=image[0][0].length-1){
                    std+=(image[i][ii][iii]-meanback)*(image[i][ii][iii]-meanback);
                }
            }
            for (int ii=1;ii<image[0].length-1;ii++){
                for (int iii=1;iii<image[0][0].length;iii+=image[0][0].length-2){
                    std+=(image[i][ii][iii]-meanback)*(image[i][ii][iii]-meanback);
                }
            }
        }
        std/=count;
        std=Math.sqrt(std);
        meanback=meanback+5*std;
        //IJ.log("meanBack + stdBack  "+meanback);
        for (int i=0;i<im.length;i++){
            //IJ.log("mean back B:"+meanback+"   "+std+"   mean+5std "+(meanback+5*std)+"   "+count);
            for (int ii=0;ii<im[0].length;ii++){
                for (int iii=0;iii<im[0][0].length;iii++){
                    im[i][ii][iii]-=meanback;
                    if (im[i][ii][iii]<0){
                        im[i][ii][iii]=0;
                    }
                }
            }
        }
        
        //normalize
        for (int i=0;i<im.length;i++){
            double som=0;
            for (int ii=0;ii<im[0].length;ii++){
                for (int iii=0;iii<im[0][0].length;iii++){
                    som+=im[i][ii][iii];
                }
            }
            //IJ.log("sum A:"+som);
            for (int ii=0;ii<im[0].length;ii++){
                for (int iii=0;iii<im[0][0].length;iii++){
                    im[i][ii][iii]/=som;
                    im[i][ii][iii]=Math.sqrt(im[i][ii][iii]);
                }
            }
        }
        
        return im;
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
    }
    
    
    
    public void imshow(Pointer deviceDiskShape,String title){
        double [] kk=new double [sizeDisk];
        Pointer hostkk = Pointer.to(kk);
        this.setDevice2Host( hostkk,deviceDiskShape, param.sizeDisk, Sizeof.DOUBLE);
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=kk[i];
        }
        ImageShow.imshow(mat,title);
    }
    
    
    void setHost2Device(Pointer device,Pointer host,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpy(device, host, size * sizeElement, cudaMemcpyKind.cudaMemcpyHostToDevice);
        
    }
    
    
    void setDevice2Host(Pointer host,Pointer device,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpy(host, device, size * sizeElement, cudaMemcpyKind.cudaMemcpyDeviceToHost);
        
    }
    
    
    
    
    
    
    public void free(){
        for (int i=0;i<device_image.length;i++){
            JCuda.cudaFree(device_image[i]);
            JCuda.cudaFree(device_real_stack[i]);
            JCuda.cudaFree(device_imag_stack[i]);
            JCuda.cudaFree(device_phase_stack[i]);
            JCuda.cudaFree(device_pupil_stack[i]);
            
        }        
        
        //sizeFFT*sizeFFT
        JCuda.cudaFree( device_realSpacial);
        JCuda.cudaFree( device_imagSpacial);
        JCuda.cudaFree( device_tmpSpacial);

        JCuda.cudaFree( device_sparseIndexEven);
        JCuda.cudaFree( device_sparseIndexOdd);
        JCuda.cudaFree( device_sparseIndexShift2D);

        //sizeFFT*sizeFFT*2
        JCuda.cudaFree( device_fftdata);

        //disk size
        JCuda.cudaFree( device_phase);
        JCuda.cudaFree( device_pupil);
        
        JCuda.cudaFree( device_kx);
        JCuda.cudaFree( device_ky);
        JCuda.cudaFree( device_kz);
        JCuda.cudaFree( device_realDisk);
        JCuda.cudaFree( device_imagDisk);
        JCuda.cudaFree( device_tmpDisk);
        JCuda.cudaFree( device_sparseIndexOddDisk);
        JCuda.cudaFree( device_sparseIndexEvenDisk);
        
        
        JCufft.cufftDestroy(plan);
        plan=null;
        
        
    }
    
}