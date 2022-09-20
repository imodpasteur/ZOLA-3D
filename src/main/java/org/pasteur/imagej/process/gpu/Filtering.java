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
public class Filtering {
    
    
    
    
    double axialRange; 
    double stepZ;
    double minPosizionZ;
    DataPhase_ dp;
    
    
    public double [] range;
    
    
    public Filtering(DataPhase_ dp,double minPosizionZ,double axialRange, double stepZ){
        this.dp=dp;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        this.minPosizionZ=minPosizionZ;
        
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
    
    
    
    public void run(double [][][] image){
        
        int sizeFilt=dp.param.sizeoutput;
        
        dp.psf.resetKz();
        
        //SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        //double position=spsfc.getPosition();
        
        double minZ=minPosizionZ;
        double maxZ=minPosizionZ+axialRange;
        
        int w,h;
        IJ.log(""+minZ+"  "+maxZ+"  "+stepZ);
        Cross predGPU=new Cross(dp.param.sizeoutput,dp,minZ, maxZ, stepZ,.3);
        
        double []range=predGPU.getRange();
        //double [][] padded=pad(psf[0]);
        
        
        
        w=image[0].length;
        h=image[0][0].length;
        float [][][] ResultFilt = new float[range.length][sizeFilt][sizeFilt];
        
        float [][][] ResultFinal = new float[range.length][w][h];
        
        
        double [][] tmp = new double[sizeFilt][sizeFilt];
        
        for (int i=0;i<Math.ceil((double)w/(double)sizeFilt);i+=1){
            for (int ii=0;ii<Math.ceil((double)h/(double)sizeFilt);ii+=1){
                for (int u=0;u<sizeFilt;u++){
                    for (int uu=0;uu<sizeFilt;uu++){
                        if ((i*sizeFilt+u)<w && (ii*sizeFilt+uu)<h){
                            tmp[u][uu]=(double)image[0][(int)(i*sizeFilt)+u][(int)(ii*sizeFilt)+uu];
                        }
                    }
                }
                predGPU.setImage(tmp);
                predGPU.convolveNormalizedInFourierDomain(ResultFilt);
                for (int y=0;y<range.length;y++){
                    for (int u=0;u<sizeFilt;u++){
                        for (int uu=0;uu<sizeFilt;uu++){
                            if ((i*sizeFilt+u)<w && (ii*sizeFilt+uu)<h){
                                ResultFinal[y][(int)(i*sizeFilt)+u][(int)(ii*sizeFilt)+uu]=ResultFilt[y][u][uu];
                            }
                        }
                    }
                }
                ImageShow.imshow(ResultFilt,"rf "+i+"  "+ii);
                
            }
            
        }
        
        
        double [][][] psf=predGPU.getPSFNonNormalized3D();
        ImageShow.imshow(psf,"psf");
        
        
        ImageShow.imshow(ResultFinal,"result filtering");
        
        
        predGPU.free();
        
        
        
    }
    
    
    
    
    
    
    
}
