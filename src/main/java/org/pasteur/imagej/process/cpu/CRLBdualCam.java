/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;


import org.pasteur.imagej.utils.*;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;

/**
 *
 * @author benoit
 */
public class CRLBdualCam {
    
    
    double axialRange; 
    double stepZ;
    
    DataPhase dp1;
    DataPhase dp2;
    
    
    double [] xCRLB;
    double [] yCRLB;
    double [] zCRLB;
    
    double photonNumber1; 
    double background1;
    double photonNumber2; 
    double background2;
    
    double maxValuePlotUm=-1;
    double [][] matrixParameter1;
    double [][] matrixParameter2;
    
    
    public CRLBdualCam(DataPhase dp1,DataPhase dp2,double axialRange, double stepZ,double photonNumber1, double background1,double photonNumber2, double background2,double maxValuePlotUm){
        this.dp1=dp1;
        this.dp2=dp2;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        this.photonNumber1=photonNumber1;
        this.background1=background1;
        this.photonNumber2=photonNumber2;
        this.background2=background2;
        this.maxValuePlotUm=maxValuePlotUm;
        matrixParameter1=null;
        matrixParameter2=null;
    }
    
    public CRLBdualCam(DataPhase dp1,DataPhase dp2,String path1,String path2,double axialRange){
        
        this.axialRange=axialRange;
        this.dp1=dp1;
        this.dp2=dp2;
        matrixParameter1=FileVectorLoader.getTableFromFile(path1, ",");
        matrixParameter2=FileVectorLoader.getTableFromFile(path2, ",");
        
        
        
        
    }
    
    
    
    public void run(){
        run(null);
    }
    
    public void run(String path){
        
        dp1.psf.resetKz();
        
        SearchPSFcenter spsfc1= new SearchPSFcenter(dp1,axialRange);
        double position1=spsfc1.getPosition();
        
        dp2.psf.resetKz();
        
        SearchPSFcenter spsfc2= new SearchPSFcenter(dp2,axialRange);
        double position2=spsfc2.getPosition();
        
        
        if (matrixParameter1==null){
            
            int nb=(int)((axialRange)/stepZ)+1;
            if (nb>0){
                double [] x_abs=new double[nb];
                xCRLB=new double[nb];
                yCRLB=new double[nb];
                zCRLB=new double[nb];
                double [][] crlb=new double[nb][5];
                int k=0;
                double [][][] im1 = new double [nb][][];
                double [][][] im2 = new double [nb][][];


                double meanCRLBX=0;
                double meanCRLBY=0;
                double meanCRLBZ=0;
                double count=0;
                for (double u=-axialRange/2+position1;k<nb;u+=stepZ,k++){
                    dp1.psf.computePSF(0, 0,this.dp1.param.Zfocus,u);
                    double [][] ppsf1=dp1.psf.getPSF();
                    im1[k]=new double[ppsf1.length][];
                    for (int j=0;j<ppsf1.length;j++){
                        im1[k][j]=new double[ppsf1[j].length];
                        for (int jj=0;jj<ppsf1[j].length;jj++){
                            im1[k][j][jj]=ppsf1[j][jj];
                        }
                    }

                    dp2.psf.computePSF(0, 0,this.dp2.param.Zfocus,u);
                    double [][] ppsf2=dp2.psf.getPSF();
                    im2[k]=new double[ppsf2.length][];
                    for (int j=0;j<ppsf2.length;j++){
                        im2[k][j]=new double[ppsf2[j].length];
                        for (int jj=0;jj<ppsf2[j].length;jj++){
                            im2[k][j][jj]=ppsf2[j][jj];
                        }
                    }

                    IJ.log("ERROR: second position might be not good :");
                    double [] res=this.computeFisher(0, 0, u,0, 0, u, .0005);

                    xCRLB[k]=res[0];
                    yCRLB[k]=res[1];
                    zCRLB[k]=res[2];
                    crlb[k][0]=u;
                    crlb[k][1]=res[0];
                    crlb[k][2]=res[1];
                    crlb[k][3]=res[2];
                    crlb[k][4]=Math.max(Math.max(res[0],res[1]),res[2]);
                    x_abs[k]=u;
                    if ((nb-k)%10==0){
                        //IJ.log("remaining : "+((nb-k)/10));
                        IJ.showProgress(1-(((double)(nb-k))/100.));
                    }
                    //IJ.log("crlb "+zCRLB[k]);
                    meanCRLBX+=xCRLB[k];
                    meanCRLBY+=yCRLB[k];
                    meanCRLBZ+=zCRLB[k];
                    count++;
                }
                meanCRLBX/=count;
                meanCRLBY/=count;
                meanCRLBZ/=count;

                //IJ.write(""+(axialRange)+" "+(stepZ)+" "+dp.phaseZer.numCoef+" "+this.photonNumber+" "+this.background+" "+meanCRLBX+" "+meanCRLBY+" "+meanCRLBZ);
                this.dp1.psf.imshowPhase();
                ImageShow.imshow(im1,"PSF 1 model");
                this.dp2.psf.imshowPhase();
                ImageShow.imshow(im1,"PSF 2 model");

                if (path!=null){
                    if (path.length()>2){
                        FileVectorLoader.saveTableInFile(path, crlb, ",");
                    }
                    else{
                        IJ.log("CRLB file not saved...you should select a path to save it");
                    }

                }
                if (nb>1){
                    for (int i=0;i<xCRLB.length;i++){
                        xCRLB[i]*=1000;
                        yCRLB[i]*=1000;
                        zCRLB[i]*=1000;

                    }
                    this.plot(x_abs, xCRLB,"CRLB(X)","Z (µm)","sigma X (nm)");
                    this.plot(x_abs, yCRLB,"CRLB(Y)","Z (µm)","sigma Y (nm)");
                    this.plot(x_abs, zCRLB,"CRLB(Z)","Z (µm)","sigma Z (nm)");
                }

            }
            else{
                IJ.log("oops: wrong range");
            }
        }
        else{
            int nb=matrixParameter1.length;
            if (nb>0){
                double [] x_abs=new double[nb];
                xCRLB=new double[nb];
                yCRLB=new double[nb];
                zCRLB=new double[nb];
                double [][] crlb=new double[nb][5];
                int k=0;
                double [][][] im1 = new double [nb][][];
                
                double [][][] im2 = new double [nb][][];
                
                
                
                for (int i=0;i<matrixParameter1.length;i++){
                //for (double u=minZ;k<nb;u+=stepZ,k++){
                    //IJ.log("i "+i+"  "+matrixParameter[i][0]+"  "+matrixParameter[i][1]+"  "+matrixParameter[i][2]+"  "+matrixParameter[i][3]+"  "+matrixParameter[i][4]);
                    this.background1=matrixParameter1[i][4];
                    this.photonNumber1=matrixParameter1[i][3];
                    this.background2=matrixParameter2[i][4];
                    this.photonNumber2=matrixParameter2[i][3];
                    
                    dp1.psf.computePSF(matrixParameter1[i][0], matrixParameter1[i][1], dp1.param.Zfocus,matrixParameter1[i][2]);
                    im1[k]=dp1.psf.getPSF();
                    
                    dp2.psf.computePSF(matrixParameter2[i][0], matrixParameter2[i][1], dp2.param.Zfocus,matrixParameter2[i][2]);
                    im2[k]=dp2.psf.getPSF();
                    
                    
                    
                    double [] res;
                    
                    res=this.computeFisher(matrixParameter1[i][0], matrixParameter1[i][1], matrixParameter1[i][2],matrixParameter2[i][0], matrixParameter2[i][1], matrixParameter2[i][2], .0005);
                    
                    xCRLB[k]=res[0];
                    yCRLB[k]=res[1];
                    zCRLB[k]=res[2];
                    crlb[k][0]=matrixParameter1[i][2];
                    crlb[k][1]=res[0];
                    crlb[k][2]=res[1];
                    crlb[k][3]=res[2];
                    crlb[k][4]=Math.max(Math.max(res[0],res[1]),res[2]);
                    x_abs[k]=matrixParameter1[i][2];
                    if ((nb-k)%10==0){
                        IJ.showProgress(1-(((double)(nb-k))/100.));
                    }
                    //IJ.log("crlb "+zCRLB[k]);
                    k++;

                }
                //ImageShow.imshow(this.dp.psf.getPhase(),"Phase");
                ImageShow.imshow(im1,"PSF model");
                ImageShow.imshow(im2,"PSF model");
                
                if (path!=null){
                    if (path.length()>2){
                        FileVectorLoader.saveTableInFile(path, crlb, ",");
                    }
                    else{
                        IJ.log("CRLB file not saved...you should select a path to save it");
                    }
                }

                if (nb>1){
                    for (int i=0;i<xCRLB.length;i++){
                        xCRLB[i]*=1000;
                        yCRLB[i]*=1000;
                        zCRLB[i]*=1000;
                        
                    }
                    this.plot(x_abs, xCRLB,"CRLB(X)","Z (µm)","sigma X (nm)");
                    this.plot(x_abs, yCRLB,"CRLB(Y)","Z (µm)","sigma Y (nm)");
                    this.plot(x_abs, zCRLB,"CRLB(Z)","Z (µm)","sigma Z (nm)");
                }

            }
            else{
                IJ.log("oops: wrong range");
            }
        }
        
    }
            
    
    
    
    
    
    
          
    double [] computeFisher(double x, double y, double z,double xc, double yc, double zc,double hdec){
        
        
        int nbParam=5;
        
        double [][] I=new double [nbParam][nbParam];
        
        double [][] I1=new double [nbParam][nbParam];
        
        double [][] I2=new double [nbParam][nbParam];
        
        
        {
            dp1.psf.computePSF(x, y, dp1.param.Zfocus,z);
            double [][] f1 =dp1.psf.getPSF();
            
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f1[i][ii]=f1[i][ii]*this.photonNumber1+this.background1;
                }
            }


            double [][][] f0=new double [nbParam][f1.length][f1[0].length];
            double [][][] f2=new double [nbParam][f1.length][f1[0].length];

            dp1.psf.computePSF(x+hdec, y, dp1.param.Zfocus,z);
            double [][] x2 =dp1.psf.getPSF();
            
            dp1.psf.computePSF(x-hdec, y, dp1.param.Zfocus,z);
            double [][] x0 =dp1.psf.getPSF();
            for (int i=0;i<x0.length;i++){
                for (int ii=0;ii<x0[0].length;ii++){
                    f0[0][i][ii]=x0[i][ii]*this.photonNumber1+this.background1;
                    f2[0][i][ii]=x2[i][ii]*this.photonNumber1+this.background1;
                }
            }
            dp1.psf.computePSF(x, y+hdec, dp1.param.Zfocus,z);
            double [][] y2 =dp1.psf.getPSF();
            dp1.psf.computePSF(x, y-hdec, dp1.param.Zfocus,z);
            double [][] y0 =dp1.psf.getPSF();
            for (int i=0;i<y0.length;i++){
                for (int ii=0;ii<y0[0].length;ii++){
                    f0[1][i][ii]=y0[i][ii]*this.photonNumber1+this.background1;
                    f2[1][i][ii]=y2[i][ii]*this.photonNumber1+this.background1;
                }
            }

            dp1.psf.computePSF(x, y, dp1.param.Zfocus,z+hdec);
            double [][] z2 =dp1.psf.getPSF();
            dp1.psf.computePSF(x, y, dp1.param.Zfocus,z-hdec);
            double [][] z0 =dp1.psf.getPSF();
            for (int i=0;i<z0.length;i++){
                for (int ii=0;ii<z0[0].length;ii++){
                    f0[2][i][ii]=z0[i][ii]*this.photonNumber1+this.background1;
                    f2[2][i][ii]=z2[i][ii]*this.photonNumber1+this.background1;
                }
            }

            double  [][] a2 =new double[f1.length][f1[0].length];
            double  [][] a0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[3][i][ii]=f1[i][ii]*(this.photonNumber1-hdec)+this.background1;
                    f2[3][i][ii]=f1[i][ii]*(this.photonNumber1+hdec)+this.background1;
                }
            }

            double  [][] b2 =new double[f1.length][f1[0].length];
            double  [][] b0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[4][i][ii]=f1[i][ii]*this.photonNumber1+(this.background1-hdec);
                    f2[4][i][ii]=f1[i][ii]*this.photonNumber1+(this.background1+hdec);
                }
            }
            
            
            
            
            

            //compute fisher
            double d1,d2;
            for (int p=0;p<nbParam;p++){
                for (int pp=0;pp<nbParam;pp++){
                    I1[p][pp]=0;
                    for (int i=0;i<f1.length;i++){
                        for (int ii=0;ii<f1[0].length;ii++){
                            d1=(f2[p][i][ii]-f0[p][i][ii])/(2*Math.abs(hdec));
                            d2=(f2[pp][i][ii]-f0[pp][i][ii])/(2*Math.abs(hdec));
                            I1[p][pp]+=(1/f1[i][ii])*d1*d2;
                        }
                    }
                }
            }
        }
        
        
        {
            dp2.psf.computePSF(xc, yc, dp2.param.Zfocus,zc);
            double [][] f1 =dp2.psf.getPSF();
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f1[i][ii]=f1[i][ii]*this.photonNumber2+this.background2;
                }
            }


            double [][][] f0=new double [nbParam][f1.length][f1[0].length];
            double [][][] f2=new double [nbParam][f1.length][f1[0].length];

            dp2.psf.computePSF(xc+hdec, yc, dp2.param.Zfocus,zc);
            double [][] x2 =dp2.psf.getPSF();
            
            dp2.psf.computePSF(xc-hdec, yc, dp2.param.Zfocus,zc);
            double [][] x0 =dp2.psf.getPSF();
            for (int i=0;i<x0.length;i++){
                for (int ii=0;ii<x0[0].length;ii++){
                    f0[0][i][ii]=x0[i][ii]*this.photonNumber2+this.background2;
                    f2[0][i][ii]=x2[i][ii]*this.photonNumber2+this.background2;
                }
            }
            dp2.psf.computePSF(xc, yc+hdec, dp2.param.Zfocus,zc);
            double [][] y2 =dp2.psf.getPSF();
            dp2.psf.computePSF(xc, yc-hdec, dp2.param.Zfocus,zc);
            double [][] y0 =dp2.psf.getPSF();
            for (int i=0;i<y0.length;i++){
                for (int ii=0;ii<y0[0].length;ii++){
                    f0[1][i][ii]=y0[i][ii]*this.photonNumber2+this.background2;
                    f2[1][i][ii]=y2[i][ii]*this.photonNumber2+this.background2;
                }
            }

            dp2.psf.computePSF(xc, yc, dp2.param.Zfocus,zc+hdec);
            double [][] z2 =dp2.psf.getPSF();
            dp2.psf.computePSF(xc, yc, dp2.param.Zfocus,zc-hdec);
            double [][] z0 =dp2.psf.getPSF();
            for (int i=0;i<z0.length;i++){
                for (int ii=0;ii<z0[0].length;ii++){
                    f0[2][i][ii]=z0[i][ii]*this.photonNumber2+this.background2;
                    f2[2][i][ii]=z2[i][ii]*this.photonNumber2+this.background2;
                }
            }

            double  [][] a2 =new double[f1.length][f1[0].length];
            double  [][] a0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[3][i][ii]=f1[i][ii]*(this.photonNumber2-hdec)+this.background2;
                    f2[3][i][ii]=f1[i][ii]*(this.photonNumber2+hdec)+this.background2;
                }
            }

            double  [][] b2 =new double[f1.length][f1[0].length];
            double  [][] b0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[4][i][ii]=f1[i][ii]*this.photonNumber2+(this.background2-hdec);
                    f2[4][i][ii]=f1[i][ii]*this.photonNumber2+(this.background2+hdec);
                }
            }
            
            
            
            
            

            //compute fisher
            double d1,d2;
            for (int p=0;p<nbParam;p++){
                for (int pp=0;pp<nbParam;pp++){
                    I2[p][pp]=0;
                    for (int i=0;i<f1.length;i++){
                        for (int ii=0;ii<f1[0].length;ii++){
                            d1=(f2[p][i][ii]-f0[p][i][ii])/(2*Math.abs(hdec));
                            d2=(f2[pp][i][ii]-f0[pp][i][ii])/(2*Math.abs(hdec));
                            I2[p][pp]+=(1/f1[i][ii])*d1*d2;
                        }
                    }
                }
            }
        }
        
        
        for (int p=0;p<nbParam;p++){
            //String s="";
            for (int pp=0;pp<nbParam;pp++){
                I[p][pp]=I2[p][pp]+I1[p][pp];
                //s+=I[p][pp]+"  ";
            }
            //IJ.log("I:   "+s);
        }
        
        
        
        //////////////////////////////////////////////////////
//        for (int p=0;p<nbParam;p++){
//            for (int pp=0;pp<nbParam;pp++){
//                if (p!=pp){
//                    I[p][pp]=0;
//                }
//            }
//        }
//        IJ.log("WARNING: Covariance modified");
        //////////////////////////////////////////////////////
        
        
        Matrixe mat = new Matrixe(I);
        try{
            mat=Matrixe.inverse(mat);
            I=mat.getMatrixe();
        }catch(Exception ee){
            double [] std =new double[nbParam];
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

            return std;
        }
        
        
        
        
        
        double [] std =new double[nbParam];
        
        for (int i=0;i<nbParam;i++){
            std[i]=Math.sqrt(I[i][i]);
        }
        
        
        
        
        
        
        return std;
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
        
        ymin=0;
        if (maxValuePlotUm>0){
            ymax=maxValuePlotUm;
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
            
}
