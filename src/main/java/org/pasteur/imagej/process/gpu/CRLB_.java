/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;


import org.pasteur.imagej.utils.FileVectorLoader;
import org.pasteur.imagej.utils.Matrixe;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.process.gpu.SearchPSFcenter_;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;

/**
 *
 * @author benoit
 */
public class CRLB_ {
    
    double [][] matrixParameter;
    
    double axialRange; 
    double stepZ;
    
    DataPhase_ dp;
    
    
    double [] xCRLB;
    double [] yCRLB;
    double [] zCRLB;
    
    double photonNumber; 
    double background;
    
    double maxValuePlotUm=-1;
    
    
    public CRLB_(DataPhase_ dp,double axialRange, double stepZ,double photonNumber, double background,double maxValuePlotUm){
        this.dp=dp;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        this.photonNumber=photonNumber;
        this.background=background;
        this.maxValuePlotUm=maxValuePlotUm;
        matrixParameter=null;
        
    }
    
    
    
    public CRLB_(DataPhase_ dp,double axialRange, double stepZ,double photonNumber, double background){
        this.dp=dp;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        this.photonNumber=photonNumber;
        this.background=background;
        matrixParameter=null;
    }
    
    
    
    
    public CRLB_(DataPhase_ dp,String path,double axialRange){
        
        this.axialRange=axialRange;
        this.dp=dp;
        matrixParameter=FileVectorLoader.getTableFromFile(path, ",");
        
        
        
    }
    
    public void run(){
        run(null);
    }
    
    public void run(String path){
        
        dp.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        double position=spsfc.getPosition();
        
        IJ.log("position "+position);
        
        boolean computeSTD=false;
        
        if (matrixParameter==null){
            
            int nb=(int)((axialRange)/stepZ)+1;
            if (nb>0){
                double [] x_abs=new double[nb];
                xCRLB=new double[nb];
                yCRLB=new double[nb];
                zCRLB=new double[nb];
                double [][] crlb=new double[nb][5];
                int k=0;
                double [][][] im = new double [nb][][];
                
                double [] std_PSFX = null;
                double [] std_PSFY = null;
                double [] std_PSFAVG = null;
                if (computeSTD==true){
                    std_PSFX = new double [nb];
                    std_PSFY = new double [nb];
                    std_PSFAVG = new double [nb];
                }
                
                double [][][] im2=null;
                
                double meanCRLBX=0;
                double meanCRLBY=0;
                double meanCRLBZ=0;
                double count=0;
                for (double u=-axialRange/2+position;k<nb;u+=stepZ,k++){
                    dp.psf.computePSF(0, 0,this.dp.param.Zfocus,u);
                    im[k]=dp.psf.getPSF();
                    
                    if (computeSTD){
                        double meanX=0;
                        double meanY=0;
                        double sum=0;
                        for (int up1=0;up1<im[k].length;up1++){
                            for (int up2=0;up2<im[k][up1].length;up2++){
                                meanX+=up1*im[k][up1][up2];
                                meanY+=up2*im[k][up1][up2];
                                sum+=im[k][up1][up2];
                            }
                        }
                        meanX/=sum;
                        meanY/=sum;

                        std_PSFX[k]=0;
                        std_PSFY[k]=0;
                        for (int up1=0;up1<im[k].length;up1++){
                            for (int up2=0;up2<im[k][up1].length;up2++){
                                std_PSFX[k]+=im[k][up1][up2]*(up1-meanX)*(up1-meanX);
                                std_PSFY[k]+=im[k][up1][up2]*(up2-meanY)*(up2-meanY);

                            }
                        }
                        std_PSFAVG[k]=std_PSFX[k]+std_PSFY[k];
                        std_PSFX[k]=Math.sqrt(std_PSFX[k]/sum);
                        std_PSFY[k]=Math.sqrt(std_PSFY[k]/sum);
                        std_PSFAVG[k]=Math.sqrt(std_PSFAVG[k]/sum);
                    }
                    
                    double [] res=this.computeFisher(0, 0, u, .0005);
                    
                    xCRLB[k]=res[0];
                    yCRLB[k]=res[1];
                    zCRLB[k]=res[2];
                    //IJ.write(""+u+","+res[0]+","+res[1]+","+res[2]);
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
        
                ImageShow.imshow(this.dp.psf.getPhase(),"Phase");
                ImageShow.imshow(this.dp.psf.getPupil(),"Pupil");
                ImageShow.imshow(im,"PSF model");
                
                if (path!=null){
                    if (path.length()>2){
                        FileVectorLoader.saveTableInFile(path, crlb, ",");
                    }
                    else{
                        IJ.log("CRLB file not saved...you should select a path to save it");
                    }
                    
                }
                if (nb>1){
                    
                    double [][] tmp= new double [xCRLB.length][4];
                    
                    for (int i=0;i<xCRLB.length;i++){
                        
                        xCRLB[i]*=1000;
                        yCRLB[i]*=1000;
                        zCRLB[i]*=1000;
                        tmp[i][0]=x_abs[i];
                        tmp[i][1]=xCRLB[i];
                        tmp[i][2]=yCRLB[i];
                        tmp[i][3]=zCRLB[i];
                        
                    }
                    //FileVectorLoader.saveTableInFile("/home/benoit/tmp/TET_photon="+this.photonNumber+"_bckg="+this.background+".csv", tmp, ",");
                    this.plot(x_abs, xCRLB,"CRLB(X)","Z (µm)","sigma X (nm)");
                    this.plot(x_abs, yCRLB,"CRLB(Y)","Z (µm)","sigma Y (nm)");
                    this.plot(x_abs, zCRLB,"CRLB(Z)","Z (µm)","sigma Z (nm)");
                    if (computeSTD){
                        this.plot(x_abs, std_PSFX,"std(Z)","Z (µm)","std PSF X (pix)");
                        this.plot(x_abs, std_PSFY,"std(Z)","Z (µm)","std PSF Y (pix)");
                        this.plot(x_abs, std_PSFAVG,"std(Z)","Z (µm)","std PSF Average (pix)");
                    }
                }

            }
            else{
                IJ.log("oops: wrong range");
            }
        }
        else{
            
            int nb=matrixParameter.length;
            if (nb>0){
                double [] x_abs=new double[nb];
                xCRLB=new double[nb];
                yCRLB=new double[nb];
                zCRLB=new double[nb];
                double [][] crlb=new double[nb][5];
                int k=0;
                double [][][] im = new double [nb][][];
                
                double [][][] im2=null;
                
                
                for (int i=0;i<matrixParameter.length;i++){
                //for (double u=minZ;k<nb;u+=stepZ,k++){
                    //IJ.log("i "+i+"  "+matrixParameter[i][0]+"  "+matrixParameter[i][1]+"  "+matrixParameter[i][2]+"  "+matrixParameter[i][3]+"  "+matrixParameter[i][4]);
                    this.background=matrixParameter[i][4];
                    this.photonNumber=matrixParameter[i][3];
                    
                    dp.psf.computePSF(matrixParameter[i][0], matrixParameter[i][1], dp.param.Zfocus,matrixParameter[i][2]);
                    im[k]=dp.psf.getPSF();
                    
                    
                    
                    double [] res;
                    
                    res=this.computeFisher(matrixParameter[i][0], matrixParameter[i][1], matrixParameter[i][2], .0005);
                    
                    xCRLB[k]=res[0];
                    yCRLB[k]=res[1];
                    zCRLB[k]=res[2];
                    crlb[k][0]=matrixParameter[i][2];
                    crlb[k][1]=res[0];
                    crlb[k][2]=res[1];
                    crlb[k][3]=res[2];
                    crlb[k][4]=Math.max(Math.max(res[0],res[1]),res[2]);
                    x_abs[k]=matrixParameter[i][2];
                    if ((nb-k)%10==0){
                        IJ.showProgress(1-(((double)(nb-k))/100.));
                    }
                    //IJ.log("crlb "+zCRLB[k]);
                    k++;

                }
                //ImageShow.imshow(this.dp.psf.getPhase(),"Phase");
                ImageShow.imshow(im,"PSF model");
                
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
        IJ.showProgress(0);
    }
            
    
    
          
    double [] computeFisher(double x, double y, double z,double hdec){
        
        
        int nbParam=5;
        
        
        
        double [][] I=new double [nbParam][nbParam];
        
        
        {
            dp.psf.computePSF(x, y, dp.param.Zfocus,z);
            double [][] f1 =dp.psf.getPSF();
            
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f1[i][ii]=f1[i][ii]*this.photonNumber+this.background;
                }
            }


            double [][][] f0=new double [nbParam][f1.length][f1[0].length];
            double [][][] f2=new double [nbParam][f1.length][f1[0].length];

            dp.psf.computePSF(x+hdec, y, dp.param.Zfocus,z);
            double [][] x2 =dp.psf.getPSF();
            
            dp.psf.computePSF(x-hdec, y, dp.param.Zfocus,z);
            double [][] x0 =dp.psf.getPSF();
            for (int i=0;i<x0.length;i++){
                for (int ii=0;ii<x0[0].length;ii++){
                    f0[0][i][ii]=x0[i][ii]*this.photonNumber+this.background;
                    f2[0][i][ii]=x2[i][ii]*this.photonNumber+this.background;
                }
            }
            dp.psf.computePSF(x, y+hdec, dp.param.Zfocus,z);
            double [][] y2 =dp.psf.getPSF();
            dp.psf.computePSF(x, y-hdec, dp.param.Zfocus,z);
            double [][] y0 =dp.psf.getPSF();
            for (int i=0;i<y0.length;i++){
                for (int ii=0;ii<y0[0].length;ii++){
                    f0[1][i][ii]=y0[i][ii]*this.photonNumber+this.background;
                    f2[1][i][ii]=y2[i][ii]*this.photonNumber+this.background;
                }
            }

            dp.psf.computePSF(x, y, dp.param.Zfocus,z+hdec);
            double [][] z2 =dp.psf.getPSF();
            dp.psf.computePSF(x, y, dp.param.Zfocus,z-hdec);
            double [][] z0 =dp.psf.getPSF();
            for (int i=0;i<z0.length;i++){
                for (int ii=0;ii<z0[0].length;ii++){
                    f0[2][i][ii]=z0[i][ii]*this.photonNumber+this.background;
                    f2[2][i][ii]=z2[i][ii]*this.photonNumber+this.background;
                }
            }

            double  [][] a2 =new double[f1.length][f1[0].length];
            double  [][] a0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[3][i][ii]=f1[i][ii]*(this.photonNumber-hdec)+this.background;
                    f2[3][i][ii]=f1[i][ii]*(this.photonNumber+hdec)+this.background;
                }
            }

            double  [][] b2 =new double[f1.length][f1[0].length];
            double  [][] b0 =new double[f1.length][f1[0].length];
            for (int i=0;i<f1.length;i++){
                for (int ii=0;ii<f1[0].length;ii++){
                    f0[4][i][ii]=f1[i][ii]*this.photonNumber+(this.background-hdec);
                    f2[4][i][ii]=f1[i][ii]*this.photonNumber+(this.background+hdec);
                }
            }

            //compute fisher
            double d1,d2;
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
        for (int p=0;p<nbParam;p++){
            String s="";
            for (int pp=0;pp<nbParam;pp++){
                s+=I[p][pp]+" ";
                    
                
            }
            IJ.log("I  "+s);
        }
        //////////////////////////////////////////////////////
        IJ.log("");
        
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //perform 3D localization precision (in case of 3D localization: stack of images)
    public void run3D(String path){
        
        
        
        dp.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        double position=spsfc.getPosition();
        
        //IJ.log("position "+position);
        
        boolean computeSTD=false;


        int nb=(int)((axialRange)/stepZ)+1;
        if (nb>0){
            double [] x_abs=new double[nb];
            xCRLB=new double[nb];
            yCRLB=new double[nb];
            zCRLB=new double[nb];
            double [][] crlb=new double[nb][5];
            int k=0;
            double [][][] im = new double [nb][][];

            double [] std_PSFX = null;
            double [] std_PSFY = null;
            double [] std_PSFAVG = null;
            if (computeSTD==true){
                std_PSFX = new double [nb];
                std_PSFY = new double [nb];
                std_PSFAVG = new double [nb];
            }

            double [][][] im2=null;

            double meanCRLBX=0;
            double meanCRLBY=0;
            double meanCRLBZ=0;
            double count=0;
            
            double start_obj=-axialRange/2+position;
            double end_obj=axialRange/2+position;
            
            double [] res=this.computeFisher3D(0., 0., position,start_obj,end_obj,stepZ, .0005);
            
            for (double u=start_obj;u<=end_obj;u+=stepZ,k++){
                dp.psf.computePSF(0, 0,u,this.dp.param.Zfocus);
                im[k]=dp.psf.getPSF();

            }
            

            //IJ.write(""+(axialRange)+" "+(stepZ)+" "+dp.phaseZer.numCoef+" "+this.photonNumber+" "+this.background+" "+meanCRLBX+" "+meanCRLBY+" "+meanCRLBZ);

            ImageShow.imshow(this.dp.psf.getPhase(),"Phase");
            ImageShow.imshow(this.dp.psf.getPupil(),"Pupil");
            ImageShow.imshow(im,"PSF model");
            
            IJ.write("CRLB X: "+res[0]+"   CRLB Y: "+res[1]+"   CRLB Z: "+res[2]);
            

        

            
        }
        
    }
    
    
    
    
    
    
    
    
    
    
    //simulate 3D stack of images and estimate CRLB for 1 molecule using all images
    double [] computeFisher3D(double x, double y, double fix_z,double start_z_obj,double end_z_obj,double z_step_obj,double hdec){
        
        
        int nbParam=5;
        
        double [][] I=new double [nbParam][nbParam];
        for (double z_obj=start_z_obj;z_obj<=end_z_obj;z_obj+=z_step_obj){
        

            {
                dp.psf.computePSF(x, y, z_obj,fix_z);
                double [][] f1 =dp.psf.getPSF();

                for (int i=0;i<f1.length;i++){
                    for (int ii=0;ii<f1[0].length;ii++){
                        f1[i][ii]=f1[i][ii]*this.photonNumber+this.background;
                    }
                }


                double [][][] f0=new double [nbParam][f1.length][f1[0].length];
                double [][][] f2=new double [nbParam][f1.length][f1[0].length];

                dp.psf.computePSF(x+hdec, y, z_obj,fix_z);
                double [][] x2 =dp.psf.getPSF();

                dp.psf.computePSF(x-hdec, y, z_obj,fix_z);
                double [][] x0 =dp.psf.getPSF();
                for (int i=0;i<x0.length;i++){
                    for (int ii=0;ii<x0[0].length;ii++){
                        f0[0][i][ii]=x0[i][ii]*this.photonNumber+this.background;
                        f2[0][i][ii]=x2[i][ii]*this.photonNumber+this.background;
                    }
                }
                dp.psf.computePSF(x, y+hdec, z_obj,fix_z);
                double [][] y2 =dp.psf.getPSF();
                dp.psf.computePSF(x, y-hdec, z_obj,fix_z);
                double [][] y0 =dp.psf.getPSF();
                for (int i=0;i<y0.length;i++){
                    for (int ii=0;ii<y0[0].length;ii++){
                        f0[1][i][ii]=y0[i][ii]*this.photonNumber+this.background;
                        f2[1][i][ii]=y2[i][ii]*this.photonNumber+this.background;
                    }
                }

                dp.psf.computePSF(x, y, z_obj,fix_z+hdec);
                double [][] z2 =dp.psf.getPSF();
                dp.psf.computePSF(x, y, z_obj,fix_z-hdec);
                double [][] z0 =dp.psf.getPSF();
                for (int i=0;i<z0.length;i++){
                    for (int ii=0;ii<z0[0].length;ii++){
                        f0[2][i][ii]=z0[i][ii]*this.photonNumber+this.background;
                        f2[2][i][ii]=z2[i][ii]*this.photonNumber+this.background;
                    }
                }

                double  [][] a2 =new double[f1.length][f1[0].length];
                double  [][] a0 =new double[f1.length][f1[0].length];
                for (int i=0;i<f1.length;i++){
                    for (int ii=0;ii<f1[0].length;ii++){
                        f0[3][i][ii]=f1[i][ii]*(this.photonNumber-hdec)+this.background;
                        f2[3][i][ii]=f1[i][ii]*(this.photonNumber+hdec)+this.background;
                    }
                }

                double  [][] b2 =new double[f1.length][f1[0].length];
                double  [][] b0 =new double[f1.length][f1[0].length];
                for (int i=0;i<f1.length;i++){
                    for (int ii=0;ii<f1[0].length;ii++){
                        f0[4][i][ii]=f1[i][ii]*this.photonNumber+(this.background-hdec);
                        f2[4][i][ii]=f1[i][ii]*this.photonNumber+(this.background+hdec);
                    }
                }

                //compute fisher
                double d1,d2;
                for (int p=0;p<nbParam;p++){
                    for (int pp=0;pp<nbParam;pp++){
                        
                        for (int i=0;i<f1.length;i++){
                            for (int ii=0;ii<f1[0].length;ii++){
                                d1=(f2[p][i][ii]-f0[p][i][ii])/(2*Math.abs(hdec));
                                d2=(f2[pp][i][ii]-f0[pp][i][ii])/(2*Math.abs(hdec));
                                I[p][pp]+=(1/f1[i][ii])*d1*d2;
                            }
                        }
                    }
                }
            }

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
        for (int p=0;p<nbParam;p++){
            String s="";
            for (int pp=0;pp<nbParam;pp++){
                s+=I[p][pp]+" ";
                    
                
            }
            IJ.log("I  "+s);
        }
        //////////////////////////////////////////////////////
        IJ.log("");
        
        Matrixe mat = new Matrixe(I);
        try{
            mat=Matrixe.inverse(mat);
            I=mat.getMatrixe();
        }catch(Exception ee){
            double [] std =new double[nbParam];
            //dont take into account covar if non inversible
            IJ.log("fisher matrix non inversible at z="+fix_z);
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
