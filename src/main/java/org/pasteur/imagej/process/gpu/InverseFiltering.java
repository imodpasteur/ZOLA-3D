/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import ij.IJ;
import org.pasteur.imagej.utils.*;
import java.util.Random;

/**
 *
 * @author benoit
 */
public class InverseFiltering {
    
    
    
    
    double axialRange; 
    double stepZ;
    
    DataPhase_ dp;
    
    
    public double [] range;
    
    
    public InverseFiltering(DataPhase_ dp,double axialRange, double stepZ){
        this.dp=dp;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        
        
    }
    
    
    double [][] pad(double [][] im){
        Random rand = new Random();
        double mean=mean(im);
        double var=var(im);
        double epsilon=0;
        
        double [][] r = new double [im.length*2][im[0].length*2];
        if((mean>0)&&(var>0)){
            for (int i=0;i<im.length*2;i++){
                for (int ii=0;ii<im[0].length*2;ii++){
                    do{
                        //r[i][ii]=rand.nextGaussian()*Math.sqrt(var)+mean;
                        r[i][ii]=mean+rand.nextGaussian()*.00001;
                    }while(r[i][ii]<0);
                }
            }
        }
        else{
            IJ.log("error in padding");
        }
        for (int i=0;i<im.length;i++){
            for (int ii=0;ii<im[0].length;ii++){
                r[i+im.length/2][ii+im[i].length/2]=im[i][ii];
            }
        }
        
        return r;
    }
    
    
    
    
    public double mean(double [][] x){
        double mean=0;
        
        double num=(x.length*x[0].length);
        for (int i=0;i<x.length;i++){
            for (int ii=0;ii<x[0].length;ii++){
                mean+=x[i][ii];
            }
        }
        mean/=num;
        
        return mean;
    }
    
    
    public double var(double [][] x){
        double mean=mean(x);
        double var=0;
        double num=(x.length*x[0].length);
        for (int i=0;i<x.length;i++){
            for (int ii=0;ii<x[0].length;ii++){
                var+=(x[i][ii]-mean)*(x[i][ii]-mean);
            }
        }
        var/=num;
        return var;
    }
    
    
    public double [][] addnoise(double [][] x,double photon,double bckg,double sigma){
        double [][] y = new double[x.length][x[0].length];
        Random r = new Random();
        for (int i=0;i<x.length;i++){
            for (int ii=0;ii<x[0].length;ii++){
                y[i][ii]=x[i][ii]*photon+r.nextGaussian()*sigma+bckg;
            }
        }
        
        return y;
    }
    
    
    //see https://cnx.org/contents/BQy9SIix@1.2:BCzesWfo@2/Deconvolution-with-Inverse-and-Weiner-Filters
    public void run(double [][][] image){
        
        dp.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        double position=spsfc.getPosition();
        
        double minZ=position-axialRange/2;
        double maxZ=position+axialRange/2;
        
        int w,h;
        
        Cross predGPU=new Cross(dp.param.sizeoutput,dp,minZ, maxZ, stepZ,.3);
        double [][][] psf=predGPU.getPSFNonNormalized3D();
        double []range=predGPU.getRange();
        double [][] padded=pad(psf[0]);
        w=padded.length;
        h=padded[0].length;
        double [][][] ResultInvFilt = new double[psf.length][][];
        
        for (int pp=0;pp<psf.length;pp++){
        //for (int pp=0;pp<1;pp++){
            

            double [][] a = new double[w][h];
            double [][] b = new double[w][h];
            double [][] c = new double[w][h];
            double [][] d = new double[w][h];


            double [][] imageA = new double[w][h];
            double [][] imageB = new double[w][h];
            double [][] imageC = new double[w][h];
            double [][] imageD = new double[w][h];

            FastFourierTransform fftFilter = new FastFourierTransform(w,h);
            FastFourierTransform fftImage = new FastFourierTransform(w,h);
            
            fftFilter.setReal(pad(psf[pp]));
            fftFilter.setImag(a);
            fftFilter.fft2D();
            a=fftFilter.getPointerRealOut2D();
            b=fftFilter.getPointerImagOut2D();
            //ImageShow.imshow(a,"Filt_FFT_Re");
            //ImageShow.imshow(b,"Filt_FFT_Im");
            
            
            
            
            
            

            //fftImage.setReal(addnoise(pad(image[0]),2000,10,2));
            fftImage.setReal(pad(image[0]));
            fftImage.setImag(imageB);
            fftImage.fft2D();
            imageA=fftImage.getPointerRealOut2D();
            imageB=fftImage.getPointerImagOut2D();
            
            
            
            double gamma=100000;
            double sig =2;// var(image[0]);
            IJ.log("sig "+sig);
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    double tmp=(a[i][ii]*a[i][ii]+b[i][ii]*b[i][ii]);
                    int method=1;
                    if (method==1){
                        if (tmp!=0){
                            if (Math.sqrt((a[i][ii]*a[i][ii]+b[i][ii]*b[i][ii])/(tmp*tmp))<gamma){;
                                c[i][ii]=a[i][ii]/tmp;
                                d[i][ii]=-b[i][ii]/tmp;
                            }
                            else{
                                c[i][ii]=gamma*(a[i][ii]/tmp)/Math.sqrt((a[i][ii]*a[i][ii]+b[i][ii]*b[i][ii])/(tmp*tmp));
                                d[i][ii]=gamma*(-b[i][ii]/tmp)/Math.sqrt((a[i][ii]*a[i][ii]+b[i][ii]*b[i][ii])/(tmp*tmp));
                            }


                        }
                        else{
                            IJ.log("WARNING zero division");
                        }
                    }
                    if (method==2){//wiener
                        double S=(imageA[i][ii]*imageA[i][ii]+imageB[i][ii]*imageB[i][ii]);
                        if (Math.random()<.01){
                            IJ.log("weight "+(w*w*h*h)*sig+"   "+S);
                        }
                        
                        //c[i][ii]=(a[i][ii]*S/(tmp*S+(w*w*h*h)*sig));
                        //d[i][ii]=(-b[i][ii]*S/(tmp*S+(w*w*h*h)*sig));
                        c[i][ii]=1/a[i][ii];//(a[i][ii]/tmp);
                        d[i][ii]=1/b[i][ii];//(b[i][ii]/tmp);
                    }
                    
                }
            }
            //ImageShow.imshow(c,"Filt_FFT_Re_inv");
            //ImageShow.imshow(d,"Filt_FFT_Im_inv");


            //ImageShow.imshow(imageA,"Image_FFT_Re");
            //ImageShow.imshow(imageB,"Image_FFT_Im");


            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    //imageC[i][ii]=imageA[i][ii]*c[i][ii]-imageB[i][ii]*d[i][ii];
                    //imageD[i][ii]=imageB[i][ii]*c[i][ii]+imageA[i][ii]*d[i][ii];
                    imageC[i][ii]=c[i][ii];
                    imageD[i][ii]=d[i][ii];

                }
            }



            fftFilter.setReal(imageC);
            fftFilter.setImag(imageD);

            fftFilter.ifft2D();
            fftFilter.shiftOutput();
            a=fftFilter.getPointerRealOut2D();
            b=fftFilter.getPointerImagOut2D();
            //ImageShow.imshow(a,"realFinal");
            //ImageShow.imshow(b,"imagFinal");

            double [][] res = new double [w][h];
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    res[i][ii]=Math.sqrt(a[i][ii]*a[i][ii]+b[i][ii]*b[i][ii]);
                }
            }
            
            
            ResultInvFilt[pp]=res;
            

        }
        ImageShow.imshow(ResultInvFilt,"result inverse filtering");
        
        
        predGPU.free();
        
        
        
    }
    
    
    
    
    
    
    
}
