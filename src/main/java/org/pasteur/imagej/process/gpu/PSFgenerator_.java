/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import org.pasteur.imagej.process.InitBackgroundAndPhotonNumber;
import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.cuda.*;
import jcuda.runtime.JCuda;

import jcuda.runtime.cudaError;
import ij.IJ;
import java.awt.Color;
import ij.gui.Plot;
import jcuda.Pointer;
import org.pasteur.imagej.utils.Zernike;

/**
 *
 * @author benoit
 */

public class PSFgenerator_ {
    //int stream=0;//not parallel process: iterative
    double [] deltaZ;
    //double [][][][] imageBlur;
    int methodLikelihood=0;
    int center;
    DataPhase_ dparam;
    
    double [][] x;
    double [][] y;
    double [][] z;
    int nbSlice;
    int nbstack;
    
    double imageLength;
    
    //double [] fit_b_0;
    //double [] fit_b_1;
    double [] fit_a;
    double [][] fit_a_each;
    double [][][] background;
    
    double epsilon=.00001;//stop criterion
    
    double [][] registrationStack;//z,x,y registration for each stack
    
    String path_calibration;
    
    int fitOrderBackground=2;
    double sigma;
    double [] zwat ;
    
    double dx=0;//drift X supposed linear
    double dy=0;//drift Y supposed linear
    double dz=0;//drift Z supposed linear
    
    double maxDrift=10;//max drift between each frame (nm)
    double axialRange,stepZ;
    
    public PSFgenerator_(int sizeFFT,double xystep,double zstep,double wavelength,double noil,double na ,String path_calibration,double sigma,int axialside,boolean withApoFactor,double axialRange, double stepZ){
        this.sigma=sigma;
        
        
        dparam = new DataPhase_(sizeFFT,sizeFFT,0,xystep,zstep,wavelength,noil,na,1.0,withApoFactor);
        dparam.setNwat(dparam.param.noil);//ca ne change rien normalement car bille collée à lamelle
        
        dparam.param.Zfocus=0;
        //astuce for registration
        
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        
        registrationStack=new double [3][nbstack];
        
        center=sizeFFT/2;
        
        
        
        deltaZ=new double [nbSlice];
        //IJ.log("center is : "+center+"  !!!!!!!!!!§§§§§§§§§§§§§§§!!!!!!!");
        for (int s=0;s<nbSlice;s++){
            if (axialside==0){
                deltaZ[s]=(((double)s)-center)*dparam.param.zstep;
            }
            else{
                deltaZ[s]=-(((double)s)-center)*dparam.param.zstep;
            }
            
        }
        
        zwat = new double [nbSlice];
        for (int i=0;i<nbSlice;i++){
            zwat[i]=0;
        }
        
        x=new double[nbstack][nbSlice];//registration of each stack
        y=new double[nbstack][nbSlice];//registration of each stack
        z=new double[nbstack][nbSlice];//registration of each stack
        
        
        
        this.path_calibration=path_calibration;
    }
    
    
    
    private void tophat(){
        
        //compute phase:
        int nbDataPerImage=dparam.param.sizeDisk;

        double maxDist=dparam.param.sizeRadiusRingPixel;

        double center=dparam.param.centerFourierImage;
        double dist;
        
        double [] phase = new double [nbDataPerImage];
        
        for (int i=0;i<nbDataPerImage;i++){
            int x=dparam.param.disk2D[i][0];
            int y=dparam.param.disk2D[i][1];


            dist=Math.sqrt(((double)x-center)*((double)x-center)+((double)y-center)*((double)y-center));
            

            if (dist<0.8*maxDist){
                phase[i]=0;
            }
            else {
                phase[i]=3.141592;
            }
            
        }
        
        dparam.phaseNonZer.setValuesPhase(phase);

        Pointer devicePhase=dparam.phaseNonZer.getPointerPhase();
        
        dparam.psf.updatePhase(devicePhase);
    }
    
    
    
    
    
    private void fermat(){
        
        
        Zernike z = new Zernike(dparam.param.size,dparam.param.sizeRadiusRingPixel,28);
        double [][][] zer=z.Z;
        ImageShow.imshow(zer,"zer");
        
        ImageShow.imshow(zer[5],"zer");
        //compute phase:
        int nbDataPerImage=dparam.param.sizeDisk;

        double a=dparam.param.sizeRadiusRingPixel;

        double center=dparam.param.centerFourierImage;
        double dist;
        
        double [] phase = new double [nbDataPerImage];
        double order=1;
        for (int i=0;i<nbDataPerImage;i++){
            int x=dparam.param.disk2D[i][0];
            int y=dparam.param.disk2D[i][1];
            double xx=((double)x-center);
            double yy=((double)y-center);
            double angle=Math.atan2(yy, xx);
            double r=xx*xx+yy*yy;
            
            
            
            // fermat
            double dn=Double.POSITIVE_INFINITY;
            double dp=Double.POSITIVE_INFINITY;
            
            for (double o=order%(int)order;o<=order;o++){
                
                double fp=a*a*(angle);
                double fn=a*a*(angle+Math.PI);
                if (Double.isNaN(fp)){
                    IJ.log(""+angle);
                }
                if (r<fp){
                    dp=fp-r;
                }

                if (r<fn){
                    dn=fn-r;

                }
                //IJ.log("fn fp "+fn+"  "+fp);
                if (dn<dp){
                    phase[i]=0;//-Math.PI*Math.sqrt(r)/a;
                }
                else{
                    phase[i]=1.5*Math.PI/2;
                }
                
                
            
            }
            //phase[i]=angle;
            
        }
        
        dparam.phaseNonZer.setValuesPhase(phase);

        Pointer devicePhase=dparam.phaseNonZer.getPointerPhase();
        
        dparam.psf.updatePhase(devicePhase);
    }
    
    
    
    
    
    
    
    private void split(){
        
        //compute phase:
        int nbDataPerImage=dparam.param.sizeDisk;

        double a=dparam.param.sizeRadiusRingPixel;

        double center=dparam.param.centerFourierImage;
        double dist;
        
        double [] phase = new double [nbDataPerImage];
        double order=1;
        for (int i=0;i<nbDataPerImage;i++){
            int x=dparam.param.disk2D[i][0];
            int y=dparam.param.disk2D[i][1];
            double xx=((double)x-center);
            double yy=((double)y-center);
            double angle=Math.atan2(yy, xx);
            double r=xx*xx+yy*yy;
            
            
            
            
            if ((xx>0)&&(yy>0)){
                phase[i]=angle;//-Math.PI*Math.sqrt(r)/a;
            }
            else if ((xx<0)&&(yy<0)){
                phase[i]=angle;//-Math.PI*Math.sqrt(r)/a;
            }
            else{
                phase[i]=Math.PI/2+angle;
            }
                
                
            
            
            //phase[i]=angle;
            
        }
        
        dparam.phaseNonZer.setValuesPhase(phase);

        Pointer devicePhase=dparam.phaseNonZer.getPointerPhase();
        
        dparam.psf.updatePhase(devicePhase);
    }
    
    
    
    
    
    
    public void run(){
        
        
        dparam.param.sigmaGaussianKernel=0;
        
        
        
        
        //this.tophat();
        //this.fermat();
        
        this.split();
        // compute psf
        
        
        
        dparam.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dparam,axialRange);
        double position=spsfc.getPosition();
        IJ.log("position"+position);
        int nb=(int)((axialRange)/stepZ)+1;
        
        
        double [][][] psf = new double [nb][][];
        
        
        int k=0;
        
        for (double u=-axialRange/2+position;k<nb;u+=stepZ,k++){
            dparam.psf.computePSF(0, 0,this.dparam.param.Zfocus,u);
            psf[k]=dparam.psf.getPSF();
        }
        
        ImageShow.imshow(psf,"psf");
        
        
        
        double [][] phaseres = dparam.psf.getPhase();
        
        
        ImageShow.imshow(phaseres,"phaseres");
        
        
    }


    public void free(){
        dparam.free();
    }



}
    
    
    

 