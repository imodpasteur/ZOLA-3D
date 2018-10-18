/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;





import org.pasteur.imagej.utils.Matrixe;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.cuda.MyCudaStream;
import org.pasteur.imagej.process.gpu.Model3D_;
import jcuda.runtime.JCuda;

import jcuda.runtime.cudaError;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import jcuda.Pointer;


/**
 *
 * @author benoit
 */
public class ZernikePhaseOptimization_ {
    
    //int stream=0;//not parallel process: iterative
    double [] deltaZ;
    
    int methodLikelihood=0;
    
    DataPhase_ dparam;
    
    
    
    double axialRange; 
    double stepZ;
    
    double [] zCRLB;
    
    
    String path_result;
    
    float photonNumber; 
    float background;
    double nwat;
    double Zfocus;
    public ZernikePhaseOptimization_(int sizeFFT,int sizePatch,double xystep,double axialRange, double stepZ,double photonNumber, double background,double wavelength,double noil,double nwat,double zfocus,double na,int zernikeCoefNumber,String path_result,boolean withApoFactor){
        
        this.photonNumber=(float)photonNumber;
        this.background=(float)background;
        
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        
        dparam = new DataPhase_(sizeFFT,sizePatch,0,xystep,stepZ,wavelength,noil,na,1,zernikeCoefNumber,withApoFactor);
        this.nwat=nwat;
        this.Zfocus=zfocus;
        
        int nb=(int)((axialRange)/stepZ)+1;
        
        
        
        
        deltaZ=new double [nb];
        //IJ.log("center is : "+center+"  !!!!!!!!!!§§§§§§§§§§§§§§§!!!!!!!");
        for (int s=0;s<nb;s++){
            deltaZ[s]=0-(axialRange/2.)+dparam.param.zstep*(double)s;
            //IJ.log("deltaZ "+deltaZ[s]);
        }
        
        
        
        dparam.setMany(nb*7);//images of 1 stack are computed in parallel * 7 for fisher matrix
        
        
        
        
        
        this.path_result=path_result;
        
        
    }
    
    
    
    
    
    public void run(int nbIter){
        
        
        
        dparam.setNwat(dparam.param.noil);
        dparam.param.Zfocus=0;
        //IJ.log("zstep "+dparam.param.zstep);
        
        this.phase_optim_zernike(nbIter,true);
        
        
        dparam.setNwat(nwat);
        dparam.param.Zfocus=this.Zfocus;
                
        
        if (dparam.param.nwat!=dparam.param.noil){
            IJ.log("process ran a second time to deal with refractive index mismatch");
            SearchPSFcenter_ spsfc= new SearchPSFcenter_(dparam,axialRange);
            double position=spsfc.getPosition();
            for (int s=0;s<deltaZ.length;s++){
                deltaZ[s]=position-(axialRange/2.)+dparam.param.zstep*(double)s;
                //IJ.log("deltaZ "+deltaZ[s]);
            }

            this.phase_optim_zernike(nbIter/2,false);
        }
        
        
        float [][][] im = new float [deltaZ.length][][];
        
        for (int i=0;i<deltaZ.length;i++){
            im[i]=dparam.psf_fMany.getPSF(i*7);
        }
        
        ImageShow.imshow(im,"PSF");
        
        
        /*double [][] rescrlb=compute_CRLB((float)0.05);
        for (int i=0;i<rescrlb.length;i++){
            String s="";
            for (int ii=0;ii<3;ii++){
                s+=rescrlb[i][ii]+",";
            }
            IJ.write(""+s);
        }*/
        
        if (path_result.length()>2){
            dparam.save(path_result);
        }
        else{
            IJ.log("calibration file not saved...you should select a path to save it");
        }
        
    }
    
    private void phase_optim_zernike(int iterations,boolean init){
        
        
        
        double hdec=.05;
        
        
        
        
        double lik=Double.MAX_VALUE;
        int [] compl=this.dparam.phaseZer.complexity;
        
        int maxCompl=0;
        for (int i=0;i<compl.length;i++){
            if (maxCompl<compl[i]){
                maxCompl=compl[i];
            }
        }
        
        
        IJ.showProgress(.02);
        if (init){
            IJ.log("initialization started");
            initPhaseZernike(hdec,15);
        }
        
        


        
        IJ.log("optimization started");
        
        
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
        //ImageShow.imshow(dparam.psf_fMany.getPhase(),"phase "+0);
        
        for (int t=0;t<iterations;t++){
            
            
            IJ.log("remaining iterations: "+(iterations-t));
            IJ.showProgress(((double)t/(double)iterations));
            
            //if (t%2==0){
            //    IJ.log("remaining iterations: "+(int)((double)iterations-(double)t/2));
            //}
            //IJ.showProgress(.02+0.5*((double)t/(double)iterations));
            
            if (t<iterations/3){
                lik=this.updatePhaseZernike(hdec,10);
                
            }
            else if (t<iterations/2){
                lik=this.updatePhaseZernike(hdec,15);
                //lik=this.updatePhaseZernike(hdec,15);
                
            }
            else{
                lik=this.updatePhaseZernike(hdec,dparam.phaseZer.numCoef);
                //lik=this.updatePhaseZernike(hdec,dparam.phaseZer.numCoef);
                
            }
            
            
            
        
            
        }
        
        
        
        
        /*for (int t=0;t<iterations;t++){
            
            if (t%2==0){
                IJ.log("remaining iterations: "+(int)((double)iterations/2-(double)t/2));
            }
            IJ.showProgress(0.5+0.5*((double)t/(double)iterations));
            if (t<iterations/3){
                lik=this.updatePhaseZernike_std(hdec,10);
                
            }
            else if (t<iterations/2){
                lik=this.updatePhaseZernike_std(hdec,15);
                //lik=this.updatePhaseZernike(hdec,15);
                
            }
            else{
                lik=this.updatePhaseZernike_std(hdec,dparam.phaseZer.numCoef);
                //lik=this.updatePhaseZernike(hdec,dparam.phaseZer.numCoef);
                
            }
            
            
        
            
        }*/
        IJ.showProgress(1);
        
        ImageShow.imshow(dparam.psf_fMany.getPhase(),"phase");
        
        //initPhasePistonOneValue(hdec);
        /*for (int i=0;i<iterations;i++){
            IJ.log("iter..."+(iterations-i));
            this.updatePhasePistonOneValue(hdec);
        }*/
        
        //ImageShow.imshow(dparam.psf_fMany.getPhase(),"phase after modif");
        
        
        
        double [][] crlb=this.compute_CRLB((float)hdec);
        double [] xcrlb= new double [crlb.length];
        double [] ycrlb= new double [crlb.length];
        double [] zcrlb= new double [crlb.length];
        
        double meanX=0;
        double meanY=0;
        double meanZ=0;
        
        for (int i=0;i<crlb.length;i++){
            xcrlb[i]=crlb[i][0];
            ycrlb[i]=crlb[i][1];
            zcrlb[i]=crlb[i][2];
            meanX+=xcrlb[i];
            meanY+=ycrlb[i];
            meanZ+=zcrlb[i];
        }
        meanX/=(double)crlb.length;
        meanY/=(double)crlb.length;
        meanZ/=(double)crlb.length;
        
        //IJ.write(""+(deltaZ[deltaZ.length-1]-deltaZ[0])+" "+(deltaZ[1]-deltaZ[0])+" "+dparam.phaseZer.numCoef+" "+this.photonNumber+" "+this.background+" "+meanX+" "+meanY+" "+meanZ);
        
        //this.plot(deltaZ, xcrlb,"CRLB (X)","Z (µm)","CRLB (µm)");
        //this.plot(deltaZ, ycrlb,"CRLB (Y)","Z (µm)","CRLB (µm)");
        //this.plot(deltaZ, zcrlb,"CRLB (Z)","Z (µm)","CRLB (µm)");
        
        
        
    }
    
    
    
    public void plot(double [] x, double [] y,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel);
        p.setFont(0, 18);
        double xmin=Double.MAX_VALUE;
        double xmax=Double.NEGATIVE_INFINITY;
        double ymin=Double.MAX_VALUE;
        double ymax=Double.NEGATIVE_INFINITY;
        for (int i=0;i<x.length;i++){
            if (xmin>x[i]){
                xmin=x[i];
            }
            if (xmax<x[i]){
                xmax=x[i];
            }
            if (ymin>y[i]){
                ymin=y[i];
            }
            if (ymax<y[i]){
                ymax=y[i];
            }
        }
        if (ymin<.025){
            ymin=0;
            ymax=.025;
        }
        p.setLimits(xmin, xmax, ymin, ymax);
        for (int ii=0;ii<x.length;ii++){
            
            p.setColor(Color.red);
            p.add("CIRCLE",  x,y);
            //p.show();
            p.setLineWidth(2);
        }
        p.show();
        
    }
    
    
    
    private double updatePhasePistonOneValue(double h){
        
        
        
        dparam.phaseZer.setA(0, 1);
        
        
        double lastLikelihood=-1;
        
        for (int i=0;i<dparam.phaseZer.nbDataPerImage;i++){
        //for (int i=0;i<10;i++){
            
            double value=dparam.phaseZer.getPistonValue(i);
            
            
            
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

            double lik2=this.compute_sumstd_CRLB((float)h);
            
            
            
            

            dparam.phaseZer.modifyOnePixelOfPiston(i, value-h);
            
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            double lik1=this.compute_sumstd_CRLB((float)h);






            

            dparam.phaseZer.modifyOnePixelOfPiston(i, value+h);
            
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            double lik3=this.compute_sumstd_CRLB((float)h);







            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)/(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }


            //IJ.log("lik("+i+") "+lik1+"  "+lik2+"  "+lik3+"  "+grad);


            boolean found=false;
            loop:for (double gamma=1;gamma>.02;gamma/=10){
                
                
                dparam.phaseZer.modifyOnePixelOfPiston(i, value-gamma*grad);
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                double lik=this.compute_sumstd_CRLB((float)h);
                
                
                
                
                if (lik<lik2){
                    lastLikelihood=lik;
                    found=true;
                    break loop;
                }
                else{
                    dparam.phaseZer.modifyOnePixelOfPiston(i, value);
                }
            }
            
            if (!found){
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            }
            

        }
        
        return lastLikelihood;
        
    }
    
    
    
    
    
    
    private double initPhasePistonOneValue(double h){
        
        
        
        dparam.phaseZer.setA(0, 1);
        
        
        double lastLikelihood=-1;
        
        for (int i=0;i<dparam.phaseZer.nbDataPerImage;i+=1){
        //for (int i=0;i<10;i++){
            IJ.log(""+(dparam.phaseZer.nbDataPerImage-i));
            double lik=Double.MAX_VALUE;
            
            double v=dparam.phaseZer.getA(i);
            
            double pi=Math.PI;
            
            double finVal=-pi;
            for (double value=v-pi;value<v+pi;value+=.1){
                dparam.phaseZer.modifyOnePixelOfPiston(i, value);
            
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                double liktmp=this.compute_sumstd_CRLB((float)h);
                if (liktmp<lik){
                    lik=liktmp;
                    finVal=value;
                }
            }
            
            dparam.phaseZer.modifyOnePixelOfPiston(i, finVal);
            

        }
        
        return lastLikelihood;
        
    }
    
    
    
    
    
    private double updatePhaseZernike_std(double h,int numbCoef){
        
        
        if (numbCoef>dparam.phaseZer.numCoef){
            numbCoef=dparam.phaseZer.numCoef;
        }
        
        double lastLikelihood=-1;
        
       
        for (int p=3;p<numbCoef;p++){
            
            if (true){
            
                
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, -h));
                
                double lik1=this.compute_sumstd_CRLB((float)h);
                
                
                

                
                
                
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                
                double lik2=this.compute_sumstd_CRLB((float)h);
                
                
                
                
                
                
                
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, h));
                
                double lik3=this.compute_sumstd_CRLB((float)h);
                
                
                
                
                
                
                
                double save=dparam.phaseZer.getA(p);
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                
                //IJ.log("lik("+p+") "+lik1+"  "+lik2+"  "+lik3+"  "+grad);
                
                
                boolean found=false;
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    
                    
                    
                    dparam.phaseZer.setA(p, save-gamma*grad);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    double lik=this.compute_sumstd_CRLB((float)h);
                
                    
                    
                    
                    if (lik<lik2){
                        lastLikelihood=lik;
                        found=true;
                        break loop;
                    }
                    else{
                        dparam.phaseZer.setA(p, save);
                    }
                }
                if (!found){
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    
                    
                }
            }
        }
        
        return lastLikelihood;
        
    }
    
    
    
    
    
    
    
    private double updatePhaseZernike(double h,int numbCoef){
        
        if (numbCoef>dparam.phaseZer.numCoef){
            numbCoef=dparam.phaseZer.numCoef;
        }
        
        
        double lastLikelihood=-1;
        
       
        for (int p=3;p<numbCoef;p++){
            
            if ((p!=4)&&(p!=12)){
            
                
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, -h));
                
                double lik1=this.compute_sum_CRLB((float)h);
                
                
                

                
                
                
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                
                double lik2=this.compute_sum_CRLB((float)h);
                
                
                
                
                
                
                
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, h));
                
                double lik3=this.compute_sum_CRLB((float)h);
                
                
                
                
                
                
                
                double save=dparam.phaseZer.getA(p);
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                
                //IJ.log("lik("+p+") "+lik1+"  "+lik2+"  "+lik3+"  "+grad);
                
                
                boolean found=false;
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    
                    
                    
                    dparam.phaseZer.setA(p, save-gamma*grad);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    double lik=this.compute_sum_CRLB((float)h);
                
                    
                    
                    
                    if (lik<lik2){
                        lastLikelihood=lik;
                        found=true;
                        break loop;
                    }
                    else{
                        dparam.phaseZer.setA(p, save);
                    }
                }
                if (!found){
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    
                    
                }
            }
        }
        
        return lastLikelihood;
        
    }
    
    
    
    
    
    
    
    private void initPhaseZernike(double h,int numbCoef){
        
        if (numbCoef>dparam.phaseZer.numCoef){
            numbCoef=dparam.phaseZer.numCoef;
        }
        
        
        for (int p=3;p<numbCoef;p++){
            if (p%5==0){
                IJ.log("step "+(int)(p/5)+"/"+(int)(numbCoef/5));
            }
            //IJ.log("init coef "+p);
            if ((p!=4)&&(p!=12)){
                dparam.phaseZer.setA(p, 0);
                
                double likelihood=Double.MAX_VALUE;
                double val=0;
                for (double up=-5*Math.PI;up<=5*Math.PI;up+=Math.PI){
                    //IJ.log("up_init "+up);
                    if (up!=0){
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, up));
                        double lik1=this.compute_sum_CRLB((float)h);
                        if (likelihood>lik1){
                            likelihood=lik1;
                            val=up;
                        }
                    }
                    //IJ.log("up "+up+"  "+lik1);
                }
                
                
                for (double up=-Math.PI;up<=Math.PI;up+=Math.PI/4){
                    //IJ.log("up_init "+up);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, up));
                    double lik1=this.compute_sum_CRLB((float)h);
                    if (likelihood>lik1){
                        likelihood=lik1;
                        val=up;
                    }
                    //IJ.log("up "+up+"  "+lik1);
                }
                
                dparam.phaseZer.setA(p, val);
                
                
                
                
                //dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                
            }
        }
        
        
    }
    
    
    
    double compute_sumstd_CRLB(float hdec){
        double [][] crlb=compute_CRLB( hdec);
        double som=0;
        double somX=0;
        double somY=0;
        for (int i=0;i<crlb.length;i++){
            som+=crlb[i][2];
            somX+=crlb[i][0];
            somY+=crlb[i][1];
        }
        
        /////////////////begin modif
        double som1=som/(double)crlb.length;
        double som1X=somX/(double)crlb.length;
        double som1Y=somY/(double)crlb.length;
        double mean=som1;
        double meanX=som1X;
        double meanY=som1Y;
        double som2=0;
        double som2X=0;
        double som2Y=0;
        for (int i=0;i<crlb.length;i++){
            som2+=(crlb[i][2]-mean)*(crlb[i][2]-mean);
            som2X+=(crlb[i][0]-meanX)*(crlb[i][0]-meanX);
            som2Y+=(crlb[i][1]-meanY)*(crlb[i][1]-meanY);
        }
        som2/=(double)crlb.length;
        som2X/=(double)crlb.length;
        som2Y/=(double)crlb.length;
        
        
        som=som1+Math.sqrt(som2);
        
        ////////////////end modif
        
        
        //add X & Y
        //som+=som1X+Math.sqrt(som2X);
        //som+=som1Y+Math.sqrt(som2Y);
        
        
        //IJ.log("X Y Z");
        
        return som;
    }
    
    
    
    
    
    double compute_sum_CRLB(float hdec){
        double [][] crlb=compute_CRLB( hdec);
        double som=0;
        double somX=0;
        double somY=0;
        for (int i=0;i<crlb.length;i++){
            som+=crlb[i][2];
            somX+=crlb[i][0];
            somY+=crlb[i][1];
        }
        
        som/=(double)crlb.length;
        somX/=(double)crlb.length;
        somY/=(double)crlb.length;
        
        
        ////////////////begin modif to delete... max optim
        /*double max=Double.NEGATIVE_INFINITY;
        for (int i=0;i<crlbZ.length;i++){
            if (crlbZ[i][2]>max){
                max=crlbZ[i][2];
            }
        }
        som=max;*/
        ////////////////end modif to delete
        
        
        
        //add X & Y
        //som+=somX+somY;
        
        return som;
    }
    
    
    
    
    double [][] compute_CRLB(float hdec){
        
        
        int nbParam=5;
        
        //dparam.psf.computePSF(dparam.param.centerX,dparam.param.centerY,z);.
        
        double [][] crlb=new double[deltaZ.length][nbParam];
        
        double [] x = new double [deltaZ.length*7];
        double [] y = new double [deltaZ.length*7];
        double [] z = new double [deltaZ.length*7];
        double [] zoil = new double [deltaZ.length*7];
        
        for (int i=0;i<deltaZ.length;i++){
            x[i*7+0]=0;
            y[i*7+0]=0;
            z[i*7+0]=deltaZ[i];
            
            
            
            x[i*7+1]=+hdec;
            y[i*7+1]=0;
            z[i*7+1]=deltaZ[i];
            
            x[i*7+2]=-hdec;
            y[i*7+2]=0;
            z[i*7+2]=deltaZ[i];
            
            
            
            x[i*7+3]=0;
            y[i*7+3]=+hdec;
            z[i*7+3]=deltaZ[i];
            
            x[i*7+4]=0;
            y[i*7+4]=-hdec;
            z[i*7+4]=deltaZ[i];
            
            
            x[i*7+5]=0;
            y[i*7+5]=0;
            z[i*7+5]=deltaZ[i]+hdec;
            
            x[i*7+6]=0;
            y[i*7+6]=0;
            z[i*7+6]=deltaZ[i]-hdec;
            
            
        }
        
        for (int i=0;i<deltaZ.length;i++){
            for (int ii=0;ii<7;ii++){
                zoil[i*7+ii]=dparam.param.Zfocus;
            }
        }
        
        dparam.psf_fMany.computePSF(x, y, zoil,z);
        
        
        
        
        double [] u=new double [nbParam];
        double [] v=new double [nbParam];
        
        for (int t=0;t<deltaZ.length;t++){
            u[0]=0;
            u[1]=0;
            u[2]=deltaZ[t];
            u[3]=photonNumber;
            u[4]=background;
            v[0]=0;
            v[1]=0;
            v[2]=deltaZ[t];
            v[3]=photonNumber;
            v[4]=background;
            
            
            
            float [][] f1 =dparam.psf_fMany.getPSF(t*7+0);
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f1[i][ii]=f1[i][ii]*this.photonNumber+this.background;
                }
            }
            
            
            float [][][] f0=new float [nbParam][f1.length][f1[0].length];
            float [][][] f2=new float [nbParam][f1.length][f1[0].length];

                        
            float [][] x2 =dparam.psf_fMany.getPSF(t*7+1);
            
            float [][] x0 =dparam.psf_fMany.getPSF(t*7+2);
            for (int i=0;i<x0.length;i++){
                for (int ii=0;ii<x0[0].length;ii++){
                    f0[0][i][ii]=x0[i][ii]*this.photonNumber+this.background;
                    f2[0][i][ii]=x2[i][ii]*this.photonNumber+this.background;
                }
            }

            
            float [][] y2 =dparam.psf_fMany.getPSF(t*7+3);
            
            float [][] y0 =dparam.psf_fMany.getPSF(t*7+4);
            for (int i=0;i<y0.length;i++){
                for (int ii=0;ii<y0[0].length;ii++){
                    f0[1][i][ii]=y0[i][ii]*this.photonNumber+this.background;
                    f2[1][i][ii]=y2[i][ii]*this.photonNumber+this.background;
                }
            }

            
            float [][] z2 =dparam.psf_fMany.getPSF(t*7+5);
            
            float [][] z0 =dparam.psf_fMany.getPSF(t*7+6);
            for (int i=0;i<z0.length;i++){
                for (int ii=0;ii<z0[0].length;ii++){
                    f0[2][i][ii]=z0[i][ii]*this.photonNumber+this.background;
                    f2[2][i][ii]=z2[i][ii]*this.photonNumber+this.background;
                }
            }

            //double  [][] a2 =new double[f1.length][f1[0].length];
            //double  [][] a0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[3][i][ii]=f1[i][ii]*(this.photonNumber-hdec)+this.background;
                    f2[3][i][ii]=f1[i][ii]*(this.photonNumber+hdec)+this.background;
                }
            }

            //double  [][] b2 =new double[f1.length][f1[0].length];
            //double  [][] b0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[4][i][ii]=f1[i][ii]*this.photonNumber+(this.background-hdec);
                    f2[4][i][ii]=f1[i][ii]*this.photonNumber+(this.background+hdec);
                }
            }
            
            
            
            //compute fisher
            double d1,d2;
            double [][] I=new double [nbParam][nbParam];
            for (int p=0;p<nbParam;p++){
                for (int pp=0;pp<nbParam;pp++){
                    I[p][pp]=0;
                    for (int i=0;i<f1.length;i++){
                        for (int ii=0;ii<f1[0].length;ii++){
                            d1=(f2[p][i][ii]-f0[p][i][ii])/(2*Math.abs(hdec));
                            d2=(f2[pp][i][ii]-f0[pp][i][ii])/(2*Math.abs(hdec));
                            I[p][pp]+=(1/f1[i][ii])*d1*d2;
                        }
                    }
                }
            }

            double [] std =new double[nbParam];
            boolean inversible=true;
            Matrixe mat = new Matrixe(I);
            try{
                mat=Matrixe.inverse(mat);
                I=mat.getMatrixe();
            }catch(Exception ee){
                
                //dont take into account covar if non inversible
                IJ.log("fisher matrix non inversible at z="+z);
                for (int i=0;i<nbParam;i++){
                    if (I[i][i]!=0){
                        std[i]=Math.sqrt(1/I[i][i]);
                    }
                    else{
                        std[i]=Double.MAX_VALUE;
                    }
                }
                inversible=false;

            }

            
            if (inversible){
                for (int i=0;i<nbParam;i++){
                    std[i]=Math.sqrt(I[i][i]);
                }
            }
            
            crlb[t]=std;
            
        }
        
        return crlb;
        
    }
    
    
    public void free(){ 
        dparam.free();
    }
    
    
    
}
