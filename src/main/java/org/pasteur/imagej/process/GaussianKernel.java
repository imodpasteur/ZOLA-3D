/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;


import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.utils.ImageShow;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.JCublas;
import jcuda.jcublas.JCublas2;
import jcuda.jcublas.cublasHandle;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaError;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import jcuda.jcufft.*;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.jcusparse.JCusparse.cusparseDgthr;
import static jcuda.jcusparse.JCusparse.cusparseDsctr;
import jcuda.jcusparse.cusparseHandle;
import static jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import jcuda.runtime.cudaMemcpyKind;
/**
 *
 * @author benoit
 */
public class GaussianKernel {
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;

    CUstream custream;
    cufftHandle plan;
    cufftHandle planD;
    int streamId=0;
    Pointer  device_gaussian;
    Pointer  device_gaussian2bleSize;
    
    Pointer  device_image;
    
    Pointer  device_tmp;
    Pointer  device_tmp2;
    
    Pointer  device_real;
    Pointer  device_imag;
    
    
    Pointer device_sparseIndexZero;
    Pointer device_sparseIndexEven;
    Pointer device_sparseIndexOdd;
    Pointer device_sparseIndexShift2D;
    Pointer device_sparseIndexOddShift2D;
    Pointer device_sparseIndexEvenShift2D;
    
    
    Pointer device_sparseIndexOddShift2DOutput;
    Pointer device_sparseIndexEvenShift2DOutput;
    
    Pointer host_sparseIndexOddShift2DOutput;
    Pointer host_sparseIndexEvenShift2DOutput;
    
    int [] sparseIndexOddShift2DOutput;
    int [] sparseIndexEvenShift2DOutput; 
    
    double [] matrix;
    int sizeImage;
    int sizeImageOutput;
    double sigmaFourier;
    int streamNumber;
    double sigmaInPixels;
    /*PhaseRetrievalParametersDouble param;
    
    public GaussianKernel(PhaseRetrievalParametersDouble param){
        this.param=param;
        
        
        matrix=new double[param.sizeDisk];
        device_gaussian=new Pointer();
        int cudaResult = JCuda.cudaMalloc(device_gaussian, param.sizeDisk * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
    }*/
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public GaussianKernel(int sizeImage,int sizeImageOutput, double sigmaInPixels,int streamNumber){
        this.sizeImage=sizeImage;
        this.sizeImageOutput=sizeImageOutput;
        matrix=new double[sizeImage*sizeImage];
        device_gaussian=new Pointer();
        int cudaResult = JCuda.cudaMalloc(device_gaussian, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        device_gaussian2bleSize=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_gaussian2bleSize, sizeImage*sizeImage*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        device_image=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_image, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        device_tmp=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        
        device_real=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_real, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        device_imag=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_imag, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        device_tmp2=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp2, sizeImage*sizeImage*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        
        int [] sparseIndexZero = new int[sizeImage*sizeImage];
        int [] sparseIndexOdd = new int[sizeImage*sizeImage];
        int [] sparseIndexEven = new int[sizeImage*sizeImage];
        int [] sparseIndexEvenShift2D = new int[sizeImage*sizeImage];
        int [] sparseIndexOddShift2D = new int[sizeImage*sizeImage];
        int [] sparseIndexShift2D = new int[sizeImage*sizeImage];
        
        for (int i=0;i<sizeImage*sizeImage;i++){
            sparseIndexZero[i]=0;
            sparseIndexOdd[i]=(i*2)+1;
            sparseIndexEven[i]=(i*2);
        }
        for (int i=0;i<sizeImage;i++){
            for (int ii=0;ii<sizeImage;ii++){
                sparseIndexShift2D[i*sizeImage+ii]=(((i+sizeImage/2)%sizeImage)*sizeImage+((ii+sizeImage/2)%sizeImage));
            }
        }
        for (int i=0;i<sizeImage;i++){
            for (int ii=0;ii<sizeImage;ii++){
                
                sparseIndexOddShift2D[i*sizeImage+ii]=2*(((i+sizeImage/2)%sizeImage)*sizeImage  + (((ii)+sizeImage/2)%sizeImage))+1;
                sparseIndexEvenShift2D[i*sizeImage+ii]=2*(((i+sizeImage/2)%sizeImage)*sizeImage + (((ii)+sizeImage/2)%sizeImage));
                
            }
        }
        
        sparseIndexEvenShift2DOutput = new int[sizeImageOutput*sizeImageOutput];
        sparseIndexOddShift2DOutput = new int[sizeImageOutput*sizeImageOutput];
        
        for (int i=0,j=((sizeImage/2)-(sizeImageOutput/2));i<sizeImageOutput;i++,j++){
            //String s="";
            for (int ii=0,jj=((sizeImage/2)-(sizeImageOutput/2));ii<sizeImageOutput;ii++,jj++){
                sparseIndexEvenShift2DOutput[i*sizeImageOutput+ii]=sparseIndexEvenShift2D[j*sizeImage+jj];
                sparseIndexOddShift2DOutput[i*sizeImageOutput+ii]=sparseIndexOddShift2D[j*sizeImage+jj];
                //s+=sparseIndexEvenShift2DOutput[i*param.size+ii]+"  ";
            }
            //IJ.log("ss "+s);
        }
        
        
        device_sparseIndexZero=new Pointer();
        cudaMalloc(device_sparseIndexZero, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexZero, Pointer.to(sparseIndexZero), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        
        device_sparseIndexEven=new Pointer();
        cudaMalloc(device_sparseIndexEven, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexEven, Pointer.to(sparseIndexEven), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        device_sparseIndexEvenShift2D=new Pointer();
        cudaMalloc(device_sparseIndexEvenShift2D, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexEvenShift2D, Pointer.to(sparseIndexEvenShift2D), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        
        device_sparseIndexOddShift2D=new Pointer();
        cudaMalloc(device_sparseIndexOddShift2D, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexOddShift2D, Pointer.to(sparseIndexOddShift2D), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        device_sparseIndexEvenShift2DOutput=new Pointer();
        cudaMalloc(device_sparseIndexEvenShift2DOutput, sizeImageOutput*sizeImageOutput * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexEvenShift2DOutput, Pointer.to(sparseIndexEvenShift2DOutput), sizeImageOutput*sizeImageOutput*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        
        device_sparseIndexOddShift2DOutput=new Pointer();
        cudaMalloc(device_sparseIndexOddShift2DOutput, sizeImageOutput*sizeImageOutput * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexOddShift2DOutput, Pointer.to(sparseIndexOddShift2DOutput), sizeImageOutput*sizeImageOutput*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        
        device_sparseIndexOdd=new Pointer();
        cudaMalloc(device_sparseIndexOdd, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexOdd, Pointer.to(sparseIndexOdd), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        device_sparseIndexShift2D=new Pointer();
        cudaMalloc(device_sparseIndexShift2D, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexShift2D, Pointer.to(sparseIndexShift2D), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        this.streamNumber=streamNumber;
        this.streamId=streamNumber;
        //IJ.log("GaussianKernel : launch stream # "+streamNumber);
        
        handlecublas=MyCudaStream.getHandleCublas(streamNumber);
        handlecusparse=MyCudaStream.getHandleCuSparse(streamNumber);
        custream=MyCudaStream.getCUstream(streamNumber);
        
        this.computeGaussianKernel(sigmaInPixels);
        
        plan = new cufftHandle();
        JCufft.cufftPlan2d(plan, sizeImage,sizeImage, cufftType.CUFFT_Z2Z);
        
        JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(streamNumber));
        
        planD = new cufftHandle();
        JCufft.cufftPlan2d(planD, sizeImage,sizeImage, cufftType.CUFFT_D2Z);
        
        JCufft.cufftSetStream(planD, MyCudaStream.getCudaStream_t(streamNumber));
        
        
    }
    
    
    /*
    
    public GaussianKernel(int sizeImage){
        this.sizeImage=sizeImage;
        matrix=new double[sizeImage*sizeImage];
        device_gaussian=new Pointer();
        int cudaResult = JCuda.cudaMalloc(device_gaussian, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        device_image=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_image, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        device_tmp=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        
        device_real=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_real, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        device_imag=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_imag, sizeImage*sizeImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
        }
        
        device_tmp2=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp2, sizeImage*sizeImage*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda gaussian kernel");return ;
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
        cudaMalloc(device_sparseIndexEven, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexEven, Pointer.to(sparseIndexEven), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        
        device_sparseIndexOdd=new Pointer();
        cudaMalloc(device_sparseIndexOdd, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexOdd, Pointer.to(sparseIndexOdd), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        device_sparseIndexShift2D=new Pointer();
        cudaMalloc(device_sparseIndexShift2D, sizeImage*sizeImage * Sizeof.INT);
        cudaMemcpyAsync(device_sparseIndexShift2D, Pointer.to(sparseIndexShift2D), sizeImage*sizeImage*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        
        
        plan = new cufftHandle();
        JCufft.cufftPlan2d(plan, sizeImage,sizeImage, cufftType.CUFFT_Z2Z);
        
    }*/
    
    
    
    
    public void setSigma(double sigmaInPixels){
        this.computeGaussianKernel(sigmaInPixels);
    }
    
    
    public double getSigma(){
        return (sigmaInPixels);
    }
            
    
    
    /*public void computeGaussianKernel(double sigma){
        
        sigmaFourier=(1./(sigma/param.xystep))*param.size/(Math.PI*2);
        double sigpow=sigmaFourier*sigmaFourier;
        double c=param.size/2;
        for (int i=0;i<param.sizeDisk;i++){
            matrix[i]=Math.exp(-.5*((param.disk2D[i][0]-c)*(param.disk2D[i][0]-c)+(param.disk2D[i][1]-c)*(param.disk2D[i][1]-c))/(sigpow));
        }
        cudaMemcpy(device_gaussian, Pointer.to(matrix), param.sizeDisk*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
    }*/
    
    public double [][] getKernel(){
        double [][] mat = new double [sizeImage][sizeImage];
        double [] vect = new double [sizeImage*sizeImage];
        
        this.setDevice2Host( Pointer.to(vect),device_gaussian, sizeImage*sizeImage, Sizeof.DOUBLE);
        for (int i=0;i<sizeImage;i++){
            for (int ii=0;ii<sizeImage;ii++){
                mat[i][ii]=vect[i*sizeImage+ii];
            }
        }
        
        
        return mat;
    }
    
    public void computeGaussianKernel(double sigmaInPixels){
        this.sigmaInPixels=sigmaInPixels;
        sigmaFourier=(1./(sigmaInPixels))*sizeImage/(Math.PI*2);
        double sigpow=sigmaFourier*sigmaFourier;
        double c=sizeImage/2;
        double j,jj;
        for (int i=0;i<sizeImage;i++){
            j=i;
            for (int ii=0;ii<sizeImage;ii++){
                jj=ii;
            matrix[i*sizeImage+ii]=Math.exp(-.5*((j-c)*(j-c)+(jj-c)*(jj-c))/(sigpow));
            }
        }
        cudaMemcpyAsync(device_tmp, Pointer.to(matrix), sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp, device_gaussian,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_tmp, device_sparseIndexEvenShift2D,device_gaussian2bleSize, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_tmp, device_sparseIndexOddShift2D,device_gaussian2bleSize, CUSPARSE_INDEX_BASE_ZERO);
    
        
        
    }
    
    
//    public void compute2DGaussianKernel(double sigmaInPixelsX,double sigmaInPixelsY,double angle){
//        if ((sigmaInPixelsX==0)||(sigmaInPixelsY==0)){
//            IJ.log("problem gaussian kernel: sigma=0");
//            return;
//        }
//        double sigmaFourierX=(1./(sigmaInPixelsX))*sizeImage/(Math.PI*2);
//        double sigpowX=sigmaFourierX*sigmaFourierX;
//        double sigmaFourierY=(1./(sigmaInPixelsY))*sizeImage/(Math.PI*2);
//        double sigpowY=sigmaFourierY*sigmaFourierY;
//        
//        
//        double c=(double)(sizeImage)/2.-.5;
//        double j,jj;
//        /*double det=sigpowX*sigpowY-sigpowXY*sigpowXY;
//        if (det==0){
//            det=sigpowX*sigpowY;
//            if (det==0){
//                IJ.log("problem gaussian kernel: matrix singular: determinant=0");
//            }
//        }
//        IJ.log("det "+det+"  "+sigpowX+"  "+sigpowY+"  "+sigpowXY);
//        double [][] matvarcovar=new double [2][2];
//        matvarcovar[0][0]=sigpowY/det;
//        matvarcovar[0][1]=-sigpowXY/det;
//        matvarcovar[1][0]=-sigpowXY/det;
//        matvarcovar[1][1]=sigpowX/det;*/
//        
//        double a,b,aa,bb,x;
//        for (int i=0;i<sizeImage;i++){
//            j=i;
//            for (int ii=0;ii<sizeImage;ii++){
//                jj=ii;
//                aa=j-c;
//                bb=jj-c;
//                a=Math.cos(angle)*aa-Math.sin(angle)*bb;
//                b=Math.sin(angle)*aa+Math.cos(angle)*bb;
//                //x=a*(matvarcovar[0][0]*a+matvarcovar[1][0]*b)+b*(matvarcovar[0][1]*a+matvarcovar[1][1]*b);
//                x=(a*a/sigpowX+b*b/sigpowY);
//                matrix[i*sizeImage+ii]=Math.exp(-.5*x);
//            }
//        }
//        cudaMemcpyAsync(device_tmp, Pointer.to(matrix), sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
//        cudaMemcpyAsync(device_gaussian, Pointer.to(matrix), sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
//        //cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp, device_gaussian,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
//    }
    
    
    public Pointer getGaussianKernel(){
        return this.device_gaussian;
    }
    
    
    public void setImage(Pointer image){
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel Set Image "+cudaResult);}
        cudaMemcpyAsync(device_image, image, sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice,MyCudaStream.getCudaStream_t(streamId));
    }
    
    
    public void setImage(double [][] image){
        if ((image.length!=this.sizeImage)&&(image[0].length!=this.sizeImage)){
            IJ.log("error different sizes. gaussian kernel filtering impossible !");
            return;
        }
        double [] vect = new double [sizeImage*sizeImage];
        for (int i=0;i<sizeImage;i++){
            for (int ii=0;ii<sizeImage;ii++){
                vect[i*sizeImage+ii]=image[i][ii];
            }
        }
        
        cudaMemcpyAsync(device_image, Pointer.to(vect), sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
    }
    
    
    public Pointer getImage(){
        return device_image;
    }
    
    public void getImage(double [][] image){
        if ((image.length!=this.sizeImage)&&(image[0].length!=this.sizeImage)){
            IJ.log("error different sizes. gaussian kernel filtering impossible !");
            return;
        }
        double [] vect = new double [sizeImage*sizeImage];
        cudaMemcpyAsync(Pointer.to(vect),device_image, sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
        //this.imshow(sizeImage*sizeImage,sizeImage,device_image, "im", "DOUBLE");
        for (int i=0;i<sizeImage;i++){
            for (int ii=0;ii<sizeImage;ii++){
                image[i][ii]=vect[i*sizeImage+ii];
            }
        }
        
    }
    
    
    public void free(){
        //IJ.log("free gaussian kernel");
        
        
        
        
        JCuda.cudaFree(device_sparseIndexZero);
        JCuda.cudaFree(device_gaussian);
        JCuda.cudaFree(device_gaussian2bleSize);
        JCuda.cudaFree(device_image);
        JCuda.cudaFree(device_tmp);
        JCuda.cudaFree(device_tmp2);
        JCuda.cudaFree(device_real);
        JCuda.cudaFree(device_imag);
        JCuda.cudaFree(device_sparseIndexOdd);
        JCuda.cudaFree(device_sparseIndexEven);
        JCuda.cudaFree(device_sparseIndexShift2D);
        JCuda.cudaFree(device_sparseIndexEvenShift2D);
        JCuda.cudaFree(device_sparseIndexOddShift2D);
        JCuda.cudaFree(device_sparseIndexEvenShift2DOutput);
        JCuda.cudaFree(device_sparseIndexOddShift2DOutput);
        
        JCufft.cufftDestroy(plan);
        plan=null;
        JCufft.cufftDestroy(planD);
        planD=null;
    }
    
    
    
    public void filter(){
        JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
        //this.imshow(sizeImage*sizeImage,sizeImage,device_image, "device_image","DOUBLE");
        
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_image, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            
        
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_tmp, device_sparseIndexEven,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);

        
        JCufft.cufftExecZ2Z(plan, device_tmp2,device_tmp2,JCufft.CUFFT_FORWARD);
        
        MyVecDouble.divScalar(custream,2*sizeImage*sizeImage, device_tmp2, device_tmp2, sizeImage);
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_real, device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_imag, device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);
        
        
        
        //////////////////////////////////////////////////////////
        
        
        MyVecDouble.mul(custream,sizeImage*sizeImage, device_real,device_real, device_gaussian);
        MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag,device_imag, device_gaussian);
        
        
        
        //////////////////////////////////////////////////

        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_real, device_sparseIndexEven,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_imag, device_sparseIndexOdd,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);

        
        
        
        JCufft.cufftExecZ2Z(plan, device_tmp2, device_tmp2, JCufft.CUFFT_INVERSE);




        MyVecDouble.divScalar(custream,sizeImage*sizeImage*2, device_tmp2, device_tmp2, sizeImage);
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_real,device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_imag,device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);



        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_real, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_imag, device_image,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);




        MyVecDouble.mul(custream,sizeImage*sizeImage, device_real, device_tmp,device_tmp);//square
        MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag, device_image,device_image);//square
        MyVecDouble.add(custream,sizeImage*sizeImage, device_image, device_real, device_imag);
        MyVecDouble.sqrt(custream,sizeImage*sizeImage, device_image, device_image);
        
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel filter 3 "+cudaResult);}
    }
    
    
    
    
    
    
    
    
    
    public void filterPSF(Pointer device_imageFullinput, Pointer device_PSFoutput){
        //IJ.log("GaussianKernel streamId "+streamId);
        
        /*
        JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
        
        
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_image, device_sparseIndexEvenShift2D,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
        */
        //this.imshow(sizeImage*sizeImage*2,sizeImage*2,device_imageFullinput, "fft","DOUBLE");
        
        JCufft.cufftExecZ2Z(plan, device_imageFullinput,device_tmp2,JCufft.CUFFT_FORWARD);
        
        
        
        
        
        
        MyVecDouble.mul(custream,sizeImage*sizeImage*2, device_tmp2,device_tmp2, device_gaussian2bleSize);
        
        
        
        JCufft.cufftExecZ2Z(plan, device_tmp2, device_tmp2, JCufft.CUFFT_INVERSE);

        
        /*
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_real,device_sparseIndexEvenShift2D, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_imag,device_sparseIndexOddShift2D, CUSPARSE_INDEX_BASE_ZERO);


        
        
        MyVecDouble.computePSF_signalsqrt(custream, sizeImage*sizeImage, device_image, device_real, device_imag, sizeImage*sizeImage);
        */
        MyVecDouble.computePSF_signalNsqrt(custream, sizeImageOutput*sizeImageOutput, device_PSFoutput, device_tmp2, sizeImage*sizeImage,device_sparseIndexEvenShift2DOutput,device_sparseIndexOddShift2DOutput);
        
    }
    
    
    
    
    
    
    
//    public void filterPSFNormalized(Pointer device_imageFullinput, Pointer device_PSFoutput){
//        //IJ.log("GaussianKernel streamId "+streamId);
//        
//        /*
//        JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
//        
//        
//        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_image, device_sparseIndexEvenShift2D,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
//        */
//        //this.imshow(sizeImage*sizeImage*2,sizeImage*2,device_imageFullinput, "fft","DOUBLE");
//        
//        JCufft.cufftExecZ2Z(plan, device_imageFullinput,device_tmp2,JCufft.CUFFT_FORWARD);
//        
//        
//        
//        
//        
//        
//        MyVecDouble.mul(custream,sizeImage*sizeImage*2, device_tmp2,device_tmp2, device_gaussian2bleSize);
//        
//        
//        JCufft.cufftExecZ2Z(plan, device_tmp2, device_tmp2, JCufft.CUFFT_INVERSE);
//
//
//        /*
//        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_real,device_sparseIndexEvenShift2D, CUSPARSE_INDEX_BASE_ZERO);
//        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_imag,device_sparseIndexOddShift2D, CUSPARSE_INDEX_BASE_ZERO);
//
//
//
//        
//        MyVecDouble.computePSF_signalsqrt(custream, sizeImage*sizeImage, device_image, device_real, device_imag, sizeImage*sizeImage);
//        */
//        
//        
//        MyVecDouble.computePSF_signalNsqrtNormalized(custream, sizeImage*sizeImage, device_PSFoutput, device_tmp2, sizeImage*sizeImage,device_sparseIndexEvenShift2D,device_sparseIndexOddShift2D,device_tmp);
//        double [] ttt = new double [10];
//        JCuda.cudaMemcpyAsync(Pointer.to(ttt), device_tmp, 10 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(streamId));
//        for (int i=0;i<10;i++){
//            IJ.log("i "+ttt[i]);
//        }
//        
//        
//    }
    
    
    
    
    
    
    
    
    public void filter(Pointer device_image){
        //IJ.log("GaussianKernel streamId "+streamId);
        
        
        
        JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
        
        
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_image, device_sparseIndexEvenShift2D,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);

        
        JCufft.cufftExecZ2Z(plan, device_tmp2,device_tmp2,JCufft.CUFFT_FORWARD);
        
        
        
        
        //this.imshow(sizeImage*sizeImage*2,sizeImage*2,device_tmp2, "fft","DOUBLE");
        
        MyVecDouble.mul(custream,sizeImage*sizeImage*2, device_tmp2,device_tmp2, device_gaussian2bleSize);
        
        
        JCufft.cufftExecZ2Z(plan, device_tmp2, device_tmp2, JCufft.CUFFT_INVERSE);


        /*
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_real,device_sparseIndexEvenShift2D, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_imag,device_sparseIndexOddShift2D, CUSPARSE_INDEX_BASE_ZERO);



        
        MyVecDouble.computePSF_signalsqrt(custream, sizeImage*sizeImage, device_image, device_real, device_imag, sizeImage*sizeImage);
        */
        MyVecDouble.computePSF_signalNsqrt(custream, sizeImage*sizeImage, device_image, device_tmp2, sizeImage*sizeImage,device_sparseIndexEvenShift2D,device_sparseIndexOddShift2D);
        
    }
    
    
    
    
    
    
    public void filterNonOptim(Pointer device_image){
        //IJ.log("GaussianKernel streamId "+streamId);
        JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_image, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            
        
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_tmp, device_sparseIndexEven,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);

        
        JCufft.cufftExecZ2Z(plan, device_tmp2,device_tmp2,JCufft.CUFFT_FORWARD);
        
        MyVecDouble.divScalar(custream,2*sizeImage*sizeImage, device_tmp2, device_tmp2, sizeImage);
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_real, device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_imag, device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);
        
        
        
        //////////////////////////////////////////////////////////
        
        
        if (false){
            //compute pupil
            MyVecDouble.mul(custream,sizeImage*sizeImage,device_tmp,device_real,device_real);
            MyVecDouble.mul(custream,sizeImage*sizeImage,device_image,device_imag,device_imag);
            MyVecDouble.add(custream,sizeImage*sizeImage,device_tmp,device_tmp,device_image);
            MyVecDouble.sqrt(custream,sizeImage*sizeImage,device_tmp,device_tmp);



            MyVecDouble.mul(custream,sizeImage*sizeImage, device_tmp, device_tmp, device_gaussian);//filtering


            double [] imag = new double [sizeImage*sizeImage];
            double [] real = new double [sizeImage*sizeImage];
            //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel filter 0 "+cudaResult);}
            cudaMemcpyAsync(Pointer.to(imag), device_imag, sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
            //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel filter 1 "+cudaResult);}
            cudaMemcpyAsync(Pointer.to(real), device_real, sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
            for (int ii=0;ii<sizeImage*sizeImage;ii++){
                imag[ii]=Math.atan2(imag[ii], real[ii]);
            }
            cudaMemcpyAsync(device_image,Pointer.to(imag),  sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));

            //cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel filter 2 "+cudaResult);}


            JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
            //compute abs of magnitude
            MyVecDouble.fabs(custream,sizeImage*sizeImage, device_tmp, device_tmp);



            MyVecDouble.cos(custream,sizeImage*sizeImage, device_real, device_image);
            MyVecDouble.sin(custream,sizeImage*sizeImage, device_imag, device_image);

            MyVecDouble.mul(custream,sizeImage*sizeImage, device_real,device_real, device_tmp);
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag,device_imag, device_tmp);
        
        }
        else{
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_real,device_real, device_gaussian);
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag,device_imag, device_gaussian);
        }
        
        //////////////////////////////////////////////////

        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_real, device_sparseIndexEven,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_imag, device_sparseIndexOdd,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);

        //this.imshow(sizeImage*sizeImage*2,sizeImage*2,device_fftdata, "fft","DOUBLE");
        JCufft.cufftExecZ2Z(plan, device_tmp2, device_tmp2, JCufft.CUFFT_INVERSE);




        MyVecDouble.divScalar(custream,sizeImage*sizeImage*2, device_tmp2, device_tmp2, sizeImage);
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_real,device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2, device_imag,device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);



        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_real, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_imag, device_image,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);




        MyVecDouble.mul(custream,sizeImage*sizeImage, device_real, device_tmp,device_tmp);//square
        MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag, device_image,device_image);//square
        MyVecDouble.add(custream,sizeImage*sizeImage, device_image, device_real, device_imag);
        MyVecDouble.sqrt(custream,sizeImage*sizeImage, device_image, device_image);
        
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda gaussianKernel filter 3 "+cudaResult);}
    }
    
    
    
    
    public void filterOldDontWork(){
        
        JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
        //this.imshow(sizeImage*sizeImage, sizeImage, device_image, "imageinit", "DOUBLE");
        
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_image, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
        
        //this.imshow(sizeImage*sizeImage, sizeImage, device_tmp, "imageshifted", "DOUBLE");
        
        //dparamFullImage.psf.imshow(totalSize,(int)Math.sqrt(this.totalSize), device_fullImagePadded2, "ImagePaddedShifted","DOUBLE");
        
        cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_tmp, device_sparseIndexEven,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
        
        //this.imshow(sizeImage*sizeImage*2, sizeImage, device_tmp2, "imageForFFTeven", "DOUBLE");
        
        //dparamFullImage.psf.imshow(totalSize*2,(int)Math.sqrt(this.totalSize), device_fftdata, "fftImage1","DOUBLE");
        
        JCufft.cufftExecZ2Z(plan, device_tmp2,device_tmp2,JCufft.CUFFT_FORWARD);
        
        MyVecDouble.divScalar(custream,sizeImage*sizeImage*2, device_tmp2, device_tmp2, sizeImage);
        
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_real, device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
        cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_imag, device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);
        
        
        
        //this.imshow(sizeImage*sizeImage, sizeImage, device_real, "imagereal", "DOUBLE");
        //this.imshow(sizeImage*sizeImage, sizeImage, device_imag, "imageimag", "DOUBLE");
        
        if (false){
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_tmp, device_real, device_real);
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_image, device_imag, device_imag);
            MyVecDouble.add(custream,sizeImage*sizeImage, device_tmp, device_tmp, device_image);     
            MyVecDouble.sqrt(custream,sizeImage*sizeImage, device_tmp, device_tmp);  //magnitude 
            this.imshow(sizeImage*sizeImage, sizeImage, device_tmp, "device_tmp", "DOUBLE");
            
            
            
            /*MyVecDouble.div(custream,sizeImage*sizeImage, device_image, device_imag, device_real);
            MyVecDouble.atan(custream,sizeImage, device_image, device_image);//phase*/
            
            
            
            double [] imag = new double [sizeImage*sizeImage];
            cudaMemcpyAsync(Pointer.to(imag), device_imag, sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
            double [] real = new double [sizeImage*sizeImage];
            cudaMemcpyAsync(Pointer.to(real), device_real, sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(streamId));
            for (int ii=0;ii<sizeImage*sizeImage;ii++){
                imag[ii]=Math.atan2(imag[ii], real[ii]);
                //IJ.log("imag "+imag[ii]);
            }
            cudaMemcpyAsync(device_image,Pointer.to(imag),  sizeImage*sizeImage*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));


            //this.imshow(sizeImage*sizeImage, sizeImage, device_tmp, "magn", "DOUBLE");
            //this.imshow(sizeImage*sizeImage, sizeImage, device_image, "phase", "DOUBLE");



            //////////////////MyVecDouble.mul(custream,sizeImage*sizeImage, device_tmp, device_tmp, device_gaussian);//filtering

            //this.imshow(sizeImage*sizeImage, sizeImage, device_gaussian, "gaussian", "DOUBLE");
            //this.imshow(sizeImage*sizeImage, sizeImage, device_tmp, "filtered", "DOUBLE");


            MyVecDouble.cos(custream,sizeImage*sizeImage, device_tmp2, device_image);
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_real, device_tmp, device_tmp2);

            MyVecDouble.sin(custream,sizeImage*sizeImage, device_tmp2, device_image);;
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag, device_tmp, device_tmp2);
            
            
            
        }
        else{
            //MyVecDouble.mul(custream,sizeImage*sizeImage, device_real, device_real, device_gaussian);//filtering
           // MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag, device_imag, device_gaussian);//filtering
            
            //this.imshow(sizeImage*sizeImage, sizeImage, device_gaussian,"matrix","DOUBLE");
        }
        
            //this.imshow(sizeImage*sizeImage, sizeImage, device_gaussian, "gauss", "DOUBLE");
            //this.imshow(sizeImage*sizeImage, sizeImage, device_real, "real", "DOUBLE");
            //this.imshow(sizeImage*sizeImage, sizeImage, device_imag, "imag", "DOUBLE");
        
            JCuda.cudaMemsetAsync(device_tmp2, 0, sizeImage*sizeImage*2 * Sizeof.DOUBLE,MyCudaStream.getCudaStream_t(streamId));
            
            cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_real, device_sparseIndexEven,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDsctr(handlecusparse, sizeImage*sizeImage, device_imag, device_sparseIndexOdd,device_tmp2, CUSPARSE_INDEX_BASE_ZERO);
            
            JCufft.cufftExecZ2Z(plan, device_tmp2,device_tmp2,JCufft.CUFFT_INVERSE);
            
            MyVecDouble.divScalar(custream,sizeImage*sizeImage*2, device_tmp2, device_tmp2, sizeImage);

            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_real, device_sparseIndexEven, CUSPARSE_INDEX_BASE_ZERO);
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_tmp2,device_imag, device_sparseIndexOdd, CUSPARSE_INDEX_BASE_ZERO);
        
            
            
            //this.imshow(sizeImage*sizeImage, sizeImage, device_real, "realaftfft", "DOUBLE");
            //this.imshow(sizeImage*sizeImage, sizeImage, device_imag, "imagaftfft", "DOUBLE");
        
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_real, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_real, device_tmp,device_tmp);//square
            
            cusparseDgthr(handlecusparse, sizeImage*sizeImage, device_imag, device_tmp,device_sparseIndexShift2D, CUSPARSE_INDEX_BASE_ZERO);
            MyVecDouble.mul(custream,sizeImage*sizeImage, device_imag, device_tmp,device_tmp);//square
            
            MyVecDouble.add(custream,sizeImage*sizeImage, device_tmp, device_real, device_image);//magnitude
            
            
            
            MyVecDouble.sqrt(custream,sizeImage*sizeImage, device_image, device_tmp);
                
            
            
            if (false){//normalize sum=1
            
                //double som=JCublas.cublasDasum(sizeImage*sizeImage,device_tmp,1);//norm sum=1
                
                Pointer p = new Pointer();
                JCuda.cudaMalloc(p, 1 * Sizeof.DOUBLE);
                double [] som = new double [1];
                JCublas2.cublasDasum(handlecublas,sizeImage*sizeImage,device_tmp,1,p);//norm sum=1
                JCuda.cudaMemcpy(Pointer.to(som), p, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);
                
                MyVecDouble.divScalar(custream,sizeImage*sizeImage, device_image, device_tmp,som[0]);//square
                MyVecDouble.sqrt(custream,sizeImage*sizeImage, device_image, device_image);
                JCuda.cudaFree(p);
            }
        
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
    
    
    
    
    
    void setHost2Device(Pointer device,Pointer host,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpy(device, host, size * sizeElement, cudaMemcpyKind.cudaMemcpyHostToDevice);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR GaussianKernel host2device cuda");return ;
        }
    }
    
    
    void setDevice2Host(Pointer host,Pointer device,int size,int sizeElement){
        int cudaResult = JCuda.cudaMemcpy(host, device, size * sizeElement, cudaMemcpyKind.cudaMemcpyDeviceToHost);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR GaussianKernel device2host cuda "+cudaResult);return ;
        }
    }
    
    
}
