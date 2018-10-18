/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;

/**
 *
 * @author benoit
 */


import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.utils.FastFourierTransform;
import ij.IJ;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.pasteur.imagej.utils.ImageShow;





/**
 *
 * @author benoit
 */
public class PSFPhase {
    
    FastFourierTransform fft;
    FastFourierTransform fft2;
    
    int [] sparseIndexDiskX;
    int [] sparseIndexDiskY;
    
    int [][] sparseIndexShift2DOutputX;
    int [][] sparseIndexShift2DOutputY;
    
    int [][] sparseIndexShift2DOutputNextX;
    int [][] sparseIndexShift2DOutputNextY;
    
    int [][] sparseIndexShift2DOutputLastX;
    int [][] sparseIndexShift2DOutputLastY;
    
    /*int [] sparseIndexOdd;
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
    int [] sparseIndexOddShift2DOutputNext;*/
    double [] phase ;
    double [] pupil ;
    double [] kx ;
    double [] ky ;
    double [] kz ;
    double [] kz_oil ;
    double [] kz_is_imaginary;
    double [] kz_oil_is_imaginary ;
    boolean success=true;
    
    
    double [] matrix;
    
    double [][] gaussianMatrix;
    
    double [][] real;
    double [][] imag;
    
    double [][] realtmp;//just pointer
    double [][] imagtmp;//just pointer
    
    double [][] realtmp2;
    double [][] imagtmp2;
    
    double [][] magn;
    
    double [] somResOne;
    
    double [][] res2D ;
    double [] res1D ;
    double [] res;
    PhaseParameters param;
    
    int cudaResult;
    int sizeoutput1;
    int sizeoutput2;
    public PSFPhase(PhaseParameters param){
        sizeoutput1=param.sizeoutput;
        sizeoutput2=param.sizeoutput*2;
        if (sizeoutput2>param.size_cpu){
            sizeoutput2=param.size_cpu;
        }
        
        
        sizeoutput2=(int)Math.pow(2, (int)Math.ceil(Math.log(sizeoutput2)/Math.log(2)));
        
        this.param=param;
    
        
        
        
        fft = new FastFourierTransform(param.size_cpu,param.size_cpu);
        fft2 = new FastFourierTransform(sizeoutput2,sizeoutput2);
        
        phase = new double[param.sizeDisk_cpu];
        pupil = new double[param.sizeDisk_cpu];
        kx = new double[param.sizeDisk_cpu];
        ky = new double[param.sizeDisk_cpu];
        kz = new double[param.sizeDisk_cpu];
        kz_oil = new double[param.sizeDisk_cpu];
        kz_is_imaginary = new double[param.sizeDisk_cpu];
        kz_oil_is_imaginary = new double[param.sizeDisk_cpu];
        
        
        sparseIndexDiskX = new int[param.sizeDisk_cpu];
        sparseIndexDiskY = new int[param.sizeDisk_cpu];
        sparseIndexShift2DOutputX = new int[sizeoutput2][sizeoutput2];
        sparseIndexShift2DOutputY = new int[sizeoutput2][sizeoutput2];
        sparseIndexShift2DOutputNextX = new int[sizeoutput2][sizeoutput2];
        sparseIndexShift2DOutputNextY = new int[sizeoutput2][sizeoutput2];
        
        sparseIndexShift2DOutputLastX = new int[sizeoutput1][sizeoutput1];
        sparseIndexShift2DOutputLastY = new int[sizeoutput1][sizeoutput1];
        
        
        
        double left=(param.nwat/(param.wavelength))*(param.nwat/(param.wavelength));
        double left_oil=(param.noil/param.wavelength)*(param.noil/param.wavelength);
        double center=param.centerFourierImage_cpu;
        int id=0;
        for (int i=0;i<param.size_cpu;i++){
            for (int ii=0;ii<param.size_cpu;ii++){
                double disk=((((i-center)/(param.size_cpu*param.xystep))*((i-center)/(param.size_cpu*param.xystep))+((ii-center)/(param.size_cpu*param.xystep))*((ii-center)/(param.size_cpu*param.xystep))));
                
                if (disk<=param.ringsize){
                    
                    
                    //double iPositShifted=((i+size/2)%size);
                    //double iiPositShifted=((ii+size/2)%size);
                    
                    kx[id]=(i-center)*(2*Math.PI)/(param.xystep*param.size_cpu);
                    ky[id]=(ii-center)*(2*Math.PI)/(param.xystep*param.size_cpu);
                    
                    
                    double kxx=((i-center)/(param.xystep*(param.size_cpu)))*((i-center)/(param.xystep*(param.size_cpu)));
                    double kyy=((ii-center)/(param.xystep*(param.size_cpu)))*((ii-center)/(param.xystep*(param.size_cpu)));


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
                    
                    if (left_oil>right){
                        kz_oil[id]=2*Math.PI*Math.sqrt(left_oil-right);
                        kz_oil_is_imaginary[id]=1;
                    }
                    else{
                        kz_oil[id]=2*Math.PI*Math.sqrt(Math.abs(left_oil-right));
                        kz_oil_is_imaginary[id]=0;
                    }
                    
                    if (param.withApoFactor){
                        if (left_oil>right){
                            pupil[id]=(float)(1/Math.pow(1-(right/left_oil),.25));//with apodization factor
                        }
                        else{
                            pupil[id]=0;
                        }
                    }
                    else{
                        pupil[id]=1/(double)Math.sqrt(param.sizeDisk_cpu);//like that -> final sum=1 (unuseful actually because we normalize at the end)
                    }
                    
                    phase[id]=0;
                    
                    sparseIndexDiskX[id]=((i+param.size_cpu/2)%param.size_cpu);
                    sparseIndexDiskY[id]=((ii+param.size_cpu/2)%param.size_cpu);
                    
                    
                    id++;
                }
                
                
                
            }
        }
        
        
        real = new double[param.size_cpu][param.size_cpu];
        imag = new double[param.size_cpu][param.size_cpu];
        realtmp2 = new double[sizeoutput2][sizeoutput2];
        imagtmp2 = new double[sizeoutput2][sizeoutput2];
        magn = new double[sizeoutput1][sizeoutput1];
        for (int i=0;i<param.size_cpu;i++){
            for (int ii=0;ii<param.size_cpu;ii++){
                real[i][ii]=0;
                imag[i][ii]=0;
            }
        }
        
        
        
        
        
        
        
        
        for (int i=0,j=((param.size_cpu/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            
            for (int ii=0,jj=((param.size_cpu/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                sparseIndexShift2DOutputX[i][ii]=((j+param.size_cpu/2)%param.size_cpu);
                sparseIndexShift2DOutputY[i][ii]=((jj+param.size_cpu/2)%param.size_cpu);
                sparseIndexShift2DOutputNextX[i][ii]=((i+sizeoutput2/2)%sizeoutput2);
                sparseIndexShift2DOutputNextY[i][ii]=(((ii)+sizeoutput2/2)%sizeoutput2);
                
                
            }
            
        }
        
        for (int i=0,j=((sizeoutput2/2)-(sizeoutput1/2));i<sizeoutput1;i++,j++){
            
            for (int ii=0,jj=((sizeoutput2/2)-(sizeoutput1/2));ii<sizeoutput1;ii++,jj++){
                sparseIndexShift2DOutputLastX[i][ii]=((j+sizeoutput2/2)%sizeoutput2);
                sparseIndexShift2DOutputLastY[i][ii]=(((jj)+sizeoutput2/2)%sizeoutput2);
                
                
            }
        }
        
        
        
        
        
        res2D = new double [sizeoutput2][sizeoutput2];
        res1D = new double [sizeoutput2*sizeoutput2];
        res = new double [sizeoutput2*sizeoutput2];
        
        
        
        this.gaussianMatrix=new double[sizeoutput2][sizeoutput2];
        
        this.computeGaussianKernel(param.sigmaGaussianKernel);
         
    }

    
    
    
    
    
    
    
    public void setSizeoutput(int sizeoutput){
        
        
        
        sizeoutput2=sizeoutput*2;
        sizeoutput1=sizeoutput;
        
        
    
        if (sizeoutput2>param.size_cpu){
            sizeoutput2=param.size_cpu;
            
        }
        
        if (sizeoutput2<=0){
            sizeoutput2=param.size_cpu;
        }
        
        if (sizeoutput1>param.size_cpu){
            sizeoutput1=param.size_cpu;
            
        }
        
        if (sizeoutput1<=0){
            sizeoutput1=param.size_cpu;
        }
        
        param.sizeoutput=sizeoutput1;
        sizeoutput2=(int)Math.pow(2, (int)Math.ceil(Math.log(sizeoutput2)/Math.log(2)));
        
        
        
        sparseIndexShift2DOutputX = new int[sizeoutput2][sizeoutput2];
        sparseIndexShift2DOutputY = new int[sizeoutput2][sizeoutput2];
        
        sparseIndexShift2DOutputNextX = new int[sizeoutput2][sizeoutput2];
        sparseIndexShift2DOutputNextY = new int[sizeoutput2][sizeoutput2];
        
        sparseIndexShift2DOutputLastX = new int[sizeoutput1][sizeoutput1];
        sparseIndexShift2DOutputLastY = new int[sizeoutput1][sizeoutput1];
        
        for (int i=0,j=((param.size_cpu/2)-(sizeoutput2/2));i<sizeoutput2;i++,j++){
            for (int ii=0,jj=((param.size_cpu/2)-(sizeoutput2/2));ii<sizeoutput2;ii++,jj++){
                sparseIndexShift2DOutputX[i][ii]=((j+param.size_cpu/2)%param.size_cpu);
                sparseIndexShift2DOutputY[i][ii]=((jj+param.size_cpu/2)%param.size_cpu);
                sparseIndexShift2DOutputNextX[i][ii]=((i+sizeoutput2/2)%sizeoutput2);
                sparseIndexShift2DOutputNextY[i][ii]=(((ii)+sizeoutput2/2)%sizeoutput2);
                
            }
        }
        
        for (int i=0,j=((sizeoutput2/2)-(sizeoutput1/2));i<sizeoutput1;i++,j++){
            for (int ii=0,jj=((sizeoutput2/2)-(sizeoutput1/2));ii<sizeoutput1;ii++,jj++){
                sparseIndexShift2DOutputLastX[i][ii]=((j+sizeoutput2/2)%sizeoutput2);
                sparseIndexShift2DOutputLastY[i][ii]=(((jj)+sizeoutput2/2)%sizeoutput2);
                
            }
        }
        
        res2D = new double [sizeoutput2][sizeoutput2];
        res1D = new double [sizeoutput2*sizeoutput2];
        res = new double [sizeoutput2*sizeoutput2];
        
        
        realtmp2 = new double[sizeoutput2][sizeoutput2];
        imagtmp2 = new double[sizeoutput2][sizeoutput2];
        magn = new double[sizeoutput1][sizeoutput1];
        
        this.gaussianMatrix=new double[sizeoutput2][sizeoutput2];
        this.computeGaussianKernel(param.sigmaGaussianKernel);
        fft2=null;
        fft2 = new FastFourierTransform(sizeoutput2,sizeoutput2);
        
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
                gaussianMatrix[i][ii]=(float)(Math.exp(-.5*((j-c)*(j-c)+(jj-c)*(jj-c))/(sigpow)));
                
            }
        }
        
        
        
    }
    
    
    
    public void updateSigmaGaussianKernel(double sigma){
        this.setSigma(sigma);
    }
    
    
    
    
    
    
    
    public void updatePhase(double [] mat){
            for (int i=0;i<param.sizeDisk_cpu;i++){
                phase[i]=mat[i];
            }
            
        
        
    }
    
    
    
    
    
    
    public void updatePhasePointer(double [] phase){
        this.phase=phase;
        
    }
    
    
    
    
    
    
    //needed when nwat change
    public void resetKz(){
        double left=(param.nwat/(param.wavelength))*(param.nwat/(param.wavelength));
        double center_cpu=param.centerFourierImage_cpu;
        int id=0;
        for (int i=0;i<param.size_cpu;i++){
            for (int ii=0;ii<param.size_cpu;ii++){
                double disk=((((i-center_cpu)/(param.size_cpu*param.xystep))*((i-center_cpu)/(param.size_cpu*param.xystep))+((ii-center_cpu)/(param.size_cpu*param.xystep))*((ii-center_cpu)/(param.size_cpu*param.xystep))));
                
                if (disk<=param.ringsize){
                    
                    
                    //double iPositShifted=((i+size/2)%size);
                    //double iiPositShifted=((ii+size/2)%size);
                    
                    
                    
                    
                    double kxx=((i-center_cpu)/(param.xystep*(param.size_cpu)))*((i-center_cpu)/(param.xystep*(param.size_cpu)));
                    double kyy=((ii-center_cpu)/(param.xystep*(param.size_cpu)))*((ii-center_cpu)/(param.xystep*(param.size_cpu)));


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
    
    
    
    public void imshow(double [] disk,String title){
        double [][] mat = new double [param.size_cpu][param.size_cpu];
        for (int i=0;i<param.sizeDisk_cpu;i++){
            mat[param.disk2D_cpu[i][0]][param.disk2D_cpu[i][1]]=disk[i];
        }
        ImageShow.imshow(mat,title);
    }
    
    
    public double [][] getPSFPointer(){
        return magn;
    }
    
    public double [][] getPSF(){
        double [][] psf = new double [this.magn.length][this.magn[0].length];
        for (int i=0;i<magn.length;i++){
            for (int ii=0;ii<magn[0].length;ii++){
                psf[i][ii]=magn[i][ii];
            }
        }
        return psf;
    }
    
    public double [] getKxPointer(){
        return kx;
    }
    
    public double [] getKyPointer(){
        return ky;
    }
    
    public double [] getKzPointer(){
        return kz;
    }
    
    public void imshowPhase(){
        this.imshow(phase, "Phase");
    }
    
    public void computePSF(double x,double y,double zoil,double zwat){
        
        
        double v=0;
        
        double p,q,r;
        for (int i=0;i<param.disk2D_cpu.length;i++){
            p=kx[i]*x + ky[i]*y + phase[i];
            q= p + kz[i]*zwat - kz_oil[i]*zoil;
            v+=q;
            r= p + kz[i]*zwat*kz_is_imaginary[i] - kz_oil[i]*zoil*kz_oil_is_imaginary[i];
            real[sparseIndexDiskX[i]][sparseIndexDiskY[i]]=pupil[i]*Math.cos(q);
            imag[sparseIndexDiskX[i]][sparseIndexDiskY[i]]=pupil[i]*Math.sin(r);
        }
        
        
        fft.setReal(real);
        fft.setImag(imag);
        fft.ifft2D();
        realtmp=fft.getPointerRealOut2D();
        imagtmp=fft.getPointerImagOut2D();
        
        
        //JCufft.cufftExecZ2Z(plan, device_fftdata, device_fftdata, JCufft.CUFFT_INVERSE);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR fft cuda 5");}
        
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                p=realtmp[sparseIndexShift2DOutputX[i][ii]][sparseIndexShift2DOutputY[i][ii]];
                q=imagtmp[sparseIndexShift2DOutputX[i][ii]][sparseIndexShift2DOutputY[i][ii]];
                realtmp2[this.sparseIndexShift2DOutputNextX[i][ii]][this.sparseIndexShift2DOutputNextY[i][ii]]=p*p+q*q;
                imagtmp2[this.sparseIndexShift2DOutputNextX[i][ii]][this.sparseIndexShift2DOutputNextY[i][ii]]=0;
            }
        }
        
        
        
        fft2.setReal(realtmp2);
        fft2.setImag(imagtmp2);
        fft2.fft2D();
        realtmp=fft2.getPointerRealOut2D();
        imagtmp=fft2.getPointerImagOut2D();
        for (int i=0;i<sizeoutput2;i++){
            for (int ii=0;ii<sizeoutput2;ii++){
                realtmp[i][ii]*=this.gaussianMatrix[sparseIndexShift2DOutputNextX[i][ii]][sparseIndexShift2DOutputNextY[i][ii]];
                imagtmp[i][ii]*=this.gaussianMatrix[sparseIndexShift2DOutputNextX[i][ii]][sparseIndexShift2DOutputNextY[i][ii]];
            }
        }
        
        fft2.setReal(realtmp);
        fft2.setImag(imagtmp);
        fft2.ifft2D();
        
        realtmp=fft2.getPointerRealOut2D();
        imagtmp=fft2.getPointerImagOut2D();
        r=0;
        for (int i=0;i<sizeoutput1;i++){
            for (int ii=0;ii<sizeoutput1;ii++){
                p=realtmp[sparseIndexShift2DOutputLastX[i][ii]][sparseIndexShift2DOutputLastY[i][ii]];
                q=imagtmp[sparseIndexShift2DOutputLastX[i][ii]][sparseIndexShift2DOutputLastY[i][ii]];
                magn[i][ii]=Math.sqrt(p*p+q*q);
                r+=magn[i][ii];
            }
        }
        
        for (int i=0;i<sizeoutput1;i++){
            for (int ii=0;ii<sizeoutput1;ii++){
                magn[i][ii]/=r;
            }
        }
        
        
        
    }
    
    
    
    
    
    
}
