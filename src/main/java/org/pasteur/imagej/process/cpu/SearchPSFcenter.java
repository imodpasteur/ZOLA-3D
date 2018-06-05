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
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.utils.Matrixe;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;

/**
 *
 * @author benoit
 */
public class SearchPSFcenter {
    DataPhase dp;
    double axialRange;
    double zfocus;
    double photonNumber=3000;
    double background=50;
    double epsilon=0.00001;//0.01nm
    public SearchPSFcenter(DataPhase dp,double axialRange){
        this.axialRange=axialRange;
        this.dp=new DataPhase(dp.param.size,dp.param.size,dp);//need same size as size FFT (large size) to compute min CRLB
        
       
        
    }
    
    public double getPosition(){
        if (axialRange<=0){
            IJ.log("oops, axial range should be positive");
            return dp.param.Zfocus;
        }
        
        
        this.zfocus=dp.param.Zfocus;
        
        double stepZ=dp.param.xystep*5;//first, range = .5µm
        double minZ=zfocus-2*axialRange-Math.abs(zfocus);
        double maxZ=zfocus+Math.abs(zfocus)+2*axialRange;
        double position;
        
        
        position=getMinPositionCRLB(minZ, maxZ, stepZ);
        
        double oldPosit=Double.POSITIVE_INFINITY;
        int iter=0;
        loop:for (int i=0;i<10;i++){
            
            minZ=position-stepZ*5;
            maxZ=position+stepZ*5;
            stepZ=stepZ/2;
            
            position=getMinPositionCRLB(minZ, maxZ, stepZ);
            
            
            if (Math.abs(oldPosit-position)<epsilon){
                break loop;
            }
            
            oldPosit=position;
            iter++;
        }
        
        
        position=getMinRangeCRLB(position-axialRange*2, position+axialRange*2, dp.param.xystep);//then, get best range around best position integrated with Zstep = pixelsize 
        
        return position;
    }
    
    
    
    double getMinRangeCRLB(double minZ, double maxZ, double stepZ){
        
        
        axialRange/=2;

        int nb=0;
        
        
/////////////////////
//        for (double u=minZ;u<maxZ;u+=stepZ){
//            nb++;
//        }
//        double [] crlbtmp=new double [nb];
//        double [] x_abs=new double[nb];
//        double [][][] im = new double [nb][][];
//        nb=0;
//        int k=0;
//////////////////////
        
        
        
        for (double u=minZ-axialRange/2;u<=maxZ+axialRange/2;u+=stepZ){
            nb++;
        }
        double [] crlb = new double [nb];
        
        ///////////////////////////////////////////////////////////
//        double [] crlbX = new double [nb];
        //////////////////////////////////////////////////////////
        nb=0;
        for (double u=minZ-axialRange/2;u<=maxZ+axialRange/2;u+=stepZ){
            double [] res=this.computeFisher(0, 0, u, .005);
            crlb[nb]=Math.sqrt(res[0]*res[0]+res[1]*res[1]);
//            crlbX[nb]=u;
//            IJ.log("crlb[nb] "+nb+"  "+crlb[nb]+"  "+crlb.length+"  "+crlbX.length);
            
            nb++;
        }
        
        double minivalue=Double.POSITIVE_INFINITY;
        double miniposit=0;
        nb=(int)(axialRange/(2*stepZ));
//        IJ.log("NB INIT "+nb);
        for (double u=minZ;u<maxZ;u+=stepZ){
            double CRLB=0;
            int nb2=0;
            for (double i=-axialRange/2;i<=axialRange/2;i++){
                CRLB+=crlb[nb+nb2];
                nb2++;
            }
            
///////////////////////
//            dp.psf.computePSF(0, 0, dp.param.Zfocus,u);
//            im[k]=dp.psf.getPSF();
//            crlbtmp[k]=CRLB;
//            x_abs[k]=u;
//            k++;
////////////////////
            
            
            if (minivalue>CRLB){
                minivalue=CRLB;
                miniposit=u;
            }
            nb++;
            
        }
        
////////////////////
//        ImageShow.imshow(im,"psf");
//        this.plot(x_abs, crlbtmp,"CRLB(X)","Z (µm)","sigma (µm), min="+miniposit);
////////////////////
        
        axialRange*=2;
        
        return miniposit;
        
    }
    
    
    
    
    
    double getMinPositionCRLB(double minZ, double maxZ, double stepZ){
        
        
//        int nb=0;
//        for (double u=minZ;u<maxZ;u+=stepZ){
//            nb++;
//        }
//        double [] crlbtmp=new double [nb];
//        double [] x_abs=new double[nb];
//        double [][][] im = new double [nb][][];
//        
        
        
        
        double minivalue=Double.POSITIVE_INFINITY;
        double miniposit=0;
        for (double u=minZ;u<maxZ;u+=stepZ){
            
            double [] res=this.computeFisher(0, 0, u, .005);
            double CRLB=Math.sqrt(res[0]*res[0]+res[1]*res[1]);
            
//            dp.psf.computePSF(0, 0, 0,u);
//            im[k]=dp.psf.getPSF();
//            crlbtmp[k]=CRLB;
//            x_abs[k]=u;
//            k++;
            
            if (minivalue>CRLB){
                minivalue=CRLB;
                miniposit=u;
            }
            
        }
        
//        ImageShow.imshow(im,"psf");
//        this.plot(x_abs, crlbtmp,"CRLB(X)","Z (µm)","sigma (µm)");
        
        
        return miniposit;
        
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
                    f0[3][i][ii]=f1[i][ii]*(this.photonNumber-hdec*100)+this.background;
                    f2[3][i][ii]=f1[i][ii]*(this.photonNumber+hdec*100)+this.background;
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
                    //IJ.log("I "+p+"  "+pp+"  "+I[p][pp]);
                }
            }
            /*ImageShow.imshow(f2,"f2");
            ImageShow.imshow(f0,"f0");
            ImageShow.imshow(f1,"f1");
            try{
                Thread.sleep(10000);
            }catch(Exception rrr){}*/
        }
        
        
        
        
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
            IJ.log(""+ylabel+"  "+i+"  "+y[i]);
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
