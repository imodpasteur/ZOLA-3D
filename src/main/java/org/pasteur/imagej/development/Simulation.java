/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.development;

import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.process.gpu.SearchPSFcenter_;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;

/**
 *
 * @author benoit
 */
public class Simulation {
    
    
    double minX;
    double maxX;
    double minY;
    double maxY;
    //double minZ;
    //double maxZ;
    
    double axialRange;
    
    double stepZ;
    
    double stepX;
    
    double stepY;
    
    int copies;
    
    DataPhase_ dp;
    
    
    double photonNumber; 
    double background;
    double backSecOrder;
    int patchSize;
    
    public Simulation(int copies,DataPhase_ dp,double axialRange, double stepZ, double stepX, double stepY,double photonNumber, double background,double backSecOrder,int patchSize){
        this.backSecOrder=backSecOrder;
        this.dp=dp;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        
        this.copies=copies;
        
        this.stepX=stepX;
        this.stepY=stepY;
        
        this.photonNumber=photonNumber;
        this.background=background;
        this.patchSize=patchSize;
    }
    
    
    /*public Simulation(int copies,DataPhase dp1,DataPhase dp2,double minZ, double maxZ, double stepZ, double stepX, double stepY,double photonNumber, double background,double backSecOrder,int patchSize){
        this.backSecOrder=backSecOrder;
        this.dp=new DataPhase[2];
        this.dp[0]=dp1;
        this.dp[1]=dp2;
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.stepZ=stepZ;
        
        this.copies=copies;
        
        this.stepX=stepX;
        this.stepY=stepY;
        
        this.photonNumber=photonNumber;
        this.background=background;
        
        this.patchSize=patchSize;
    }*/
    
    
    public Simulation(int copies,DataPhase_ dp,double minX, double maxX,double minY, double maxY,double axialRange, double stepZ, double stepX, double stepY,double photonNumber, double background,double backSecOrder,int patchSize){
        this.backSecOrder=backSecOrder;
        this.dp=dp;
        this.minX=minX;
        this.maxX=maxX;
        this.minY=minY;
        this.maxY=maxY;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        
        this.copies=copies;
        
        this.stepX=stepX;
        this.stepY=stepY;
        
        this.photonNumber=photonNumber;
        this.background=background;
        this.patchSize=patchSize;
    }
    
    /*public Simulation(int copies,DataPhase dp1,DataPhase dp2,double minX, double maxX,double minY, double maxY,double minZ, double maxZ, double stepZ, double stepX, double stepY,double photonNumber, double background,double backSecOrder,int patchSize){
        this.backSecOrder=backSecOrder;
        this.dp=new DataPhase[2];
        this.dp[0]=dp1;
        this.dp[1]=dp2;
        this.minX=minX;
        this.maxX=maxX;
        this.minY=minY;
        this.maxY=maxY;
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.stepZ=stepZ;
        
        this.copies=copies;
        
        this.stepX=stepX;
        this.stepY=stepY;
        
        this.photonNumber=photonNumber;
        this.background=background;
        this.patchSize=patchSize;
    }*/
    
    
    public void run(){
        
        dp.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        double position=spsfc.getPosition();
        
        
        
        
        PolynomialFit pf = new PolynomialFit(2,1,2);
        pf.log();
        pf.a[0][0]=background;
        pf.a[0][2]=backSecOrder;
        pf.a[0][5]=backSecOrder;
        pf.log();
        double [][] background=new double[patchSize][patchSize];
        for (int i=0;i<patchSize;i++){
            for (int ii=0;ii<patchSize;ii++){
                double [] v = new double [2];
                v[0]=i-patchSize/2;
                v[1]=ii-patchSize/2;
                v=pf.transform(v);
                background[i][ii]=v[0];
            }
        }
        ImageShow.imshow(background,"background");
        
        int nb=(int)((axialRange)/stepZ)+1;
        
        if (nb>0){
            
            int k=0;
            double [][] psf;
            float [][] psfmany;
            double [][][] im = new double [nb*copies][patchSize][patchSize];
            double [][][] immany = new double [nb*copies][patchSize][patchSize];
            double [][][] imnoise = new double [nb*copies][patchSize][patchSize];
            dp.setMany(1);
            
            double [] xposit=new double [1];
            double [] yposit=new double [1];
            double [] zoilposit=new double [1];
            double [] zposit=new double [1];
            for (double u=-axialRange/2+position,uu=-.0,uuu=-.0;k<nb;u+=stepZ,uu+=stepX,uuu+=stepY,k++){
            //for (double u=minZ,uu=-.0,uuu=-.0;k<nb;u+=stepZ,uu+=stepX,uuu+=stepY,k++){
                IJ.showProgress(k, nb);
                xposit[0]=uu;
                yposit[0]=uuu;
                zoilposit[0]=this.dp.param.Zfocus;
                zposit[0]=u;
                dp.psf_fMany.computePSF(xposit, yposit, zoilposit, zposit);
                dp.psf.computePSF(uu, uuu, zoilposit[0],u);
                psf=dp.psf.getPSF();
                psfmany=dp.psf_fMany.getPSF(0);
                if ((patchSize>=psf.length)&&(patchSize>=psf[0].length)){
                    for (int y=0;y<copies;y++){
                        int decX=0;
                        int decY=0;
                        if ((patchSize>psf.length)&&(patchSize>psf[0].length)){
                            decX=(patchSize-psf.length)/2;
                            decY=(patchSize-psf[0].length)/2;
                            for (int z=0;z<patchSize;z++){
                                for (int zz=0;zz<patchSize;zz++){
                                    im[k*copies+y][z][zz]=background[z][zz];
                                    immany[k*copies+y][z][zz]=background[z][zz];
                                    imnoise[k*copies+y][z][zz]=background[z][zz];
                                }
                            }
                        }

                        for (int z=0;z<psf.length;z++){
                            for (int zz=0;zz<psf[z].length;zz++){
                                im[k*copies+y][z+decX][zz+decY]=psf[z][zz]*photonNumber+background[z][zz];
                                immany[k*copies+y][z+decX][zz+decY]=(double)psfmany[z][zz]*photonNumber+background[z][zz];
                                imnoise[k*copies+y][z+decX][zz+decY]=this.poissonValueLog(im[k*copies+y][z+decX][zz+decY]);
                            }
                        }
                        for (int z=0;z<patchSize;z++){
                            for (int zz=0;zz<patchSize;zz++){
                                imnoise[k*copies+y][z][zz]=this.poissonValueLog(im[k*copies+y][z][zz]);
                            }
                        }
                    }
                }
                else{
                    IJ.log("problem patch size image simulation");
                }
                    
                

                
            }
            
            

            ImageShow.imshow(im,"PSF/woNoise/"+copies+"copies/cam");
            ImageShow.imshow(immany,"PSFmany/woNoise/"+copies+"copies/cam");
            ImageShow.imshow(imnoise,"PSF/wNoise/"+copies+"copies/cam");
            
            

        }
        else{
            IJ.log("oops: wrong range");
        }
        
        
    }
            
    
    
    
    /*
    
    public void run2(){
        
        
        classes.PolynomialFit pf = new classes.PolynomialFit(2,1,2);
        pf.log();
        pf.a[0][0]=background;
        pf.a[0][2]=backSecOrder;
        pf.a[0][5]=backSecOrder;
        pf.log();
        double [][] background=new double[patchSize][patchSize];
        for (int i=0;i<patchSize;i++){
            for (int ii=0;ii<patchSize;ii++){
                double [] v = new double [2];
                v[0]=i-patchSize/2;
                v[1]=ii-patchSize/2;
                v=pf.transform(v);
                background[i][ii]=v[0];
            }
        }
        classes.ImageShow.imshow(background,"background");
        
        
        int nbX=(int)((maxX-minX)/stepX)+1;
        int nbY=(int)((maxY-minY)/stepY)+1;
        int nbZ=(int)((maxZ-minZ)/stepZ)+1;
        if ((nbX>0)&&(nbY>0)&&(nbY>0)){
            
            for (int p=0;p<dp.length;p++){
                dp[p].setMany(nbX);
                
            }
            double [] xposit=new double [nbX];
            double [] yposit=new double [nbX];
            double [] zoilposit=new double [nbX];
            double [] zposit=new double [nbX];
                        
            int k=0,kk=0,kkk=0;
            double [][] psf;
            float [][] psfmany;
            double [][][][] im = new double [dp.length][nbZ*nbX*nbY*copies][patchSize][patchSize];
            double [][][][] immany = new double [dp.length][nbZ*nbX*nbY*copies][patchSize][patchSize];
            double [][][][] imnoise = new double [dp.length][nbZ*nbX*nbY*copies][patchSize][patchSize];
            int kid=0;
            for (double u=minZ;k<nbZ;u+=stepZ,k++){
                IJ.showProgress(k, nbZ);
                kk=0;
                for (double uu=minY;kk<nbY;uu+=stepY,kk++){
                    String s="";
                    kkk=0;
                        
                    for (double uuu=minX;kkk<nbX;uuu+=stepX,uuu+=stepX,kkk++){
                        xposit[kkk]=uuu;
                        yposit[kkk]=uu;
                        zposit[kkk]=u;
                        zoilposit[kkk]=0;
                    }
                    for (int p=0;p<dp.length;p++){
                        dp[p].psf_fMany.computePSF(xposit, yposit, zoilposit, zposit);
                    }
                    for (double uuu=minX;kkk<nbX;uuu+=stepX,uuu+=stepX,kkk++){
                        
                        s+="x:"+((int)(uuu*1000))+" y:"+((int)(uu*1000))+" z:"+((int)(u*1000))+"\n";
                        
                        
                        for (int p=0;p<dp.length;p++){
                            
                            
                            
                            psfmany=dp[p].psf_fMany.getPSF(kkk);
                            
                            dp[p].psf.computePSF(uu, uuu,0, u);
                            psf=dp[p].psf.getPSF();
                            if ((patchSize>=psf.length)&&(patchSize>=psf[0].length)){
                                for (int y=0;y<copies;y++){
                                    int decX=0;
                                    int decY=0;
                                    if ((patchSize>psf.length)&&(patchSize>psf[0].length)){
                                        decX=(patchSize-psf.length)/2;
                                        decY=(patchSize-psf[0].length)/2;
                                        for (int z=0;z<patchSize;z++){
                                            for (int zz=0;zz<patchSize;zz++){
                                                im[p][kid*copies+y][z][zz]=background[z][zz];
                                                immany[p][kid*copies+y][z][zz]=background[z][zz];
                                                imnoise[p][kid*copies+y][z][zz]=background[z][zz];
                                            }
                                        }
                                    }

                                    for (int z=0;z<psf.length;z++){
                                        for (int zz=0;zz<psf[z].length;zz++){
                                            im[p][kid*copies+y][z+decX][zz+decY]=psf[z][zz]*photonNumber+background[z][zz];
                                            immany[p][kid*copies+y][z+decX][zz+decY]=psfmany[z][zz]*photonNumber+background[z][zz];
                                            imnoise[p][kid*copies+y][z+decX][zz+decY]=this.poissonValueLog(im[p][kid*copies+y][z+decX][zz+decY]);
                                        }
                                    }
                                    for (int z=0;z<patchSize;z++){
                                        for (int zz=0;zz<patchSize;zz++){
                                            imnoise[p][kid*copies+y][z][zz]=this.poissonValueLog(im[p][kid*copies+y][z][zz]);
                                        }
                                    }
                                }
                            }
                            else{
                                IJ.log("problem patch size image simulation");
                            }

                        }
                        kid++;
                    }
                    IJ.write(""+s);
                }

                
            }
            
            
            
            for (int p=0;p<dp.length;p++){
                ImageShow.imshow(im[p],"PSF/woNoise/"+copies+"copies/cam"+p);
                ImageShow.imshow(immany[p],"PSFmany/woNoise/"+copies+"copies/cam"+p);
                ImageShow.imshow(imnoise[p],"PSF/wNoise/"+copies+"copies/cam"+p);
            }
            

        }
        else{
            IJ.log("oops: wrong range");
        }
        
        
    }
            
    
    */
    
    
    int poissonValueLog(double pixVal) {
        
      java.util.Random r = new java.util.Random();
      double L = (-(pixVal));
      int k = 0;
      double p = 0;
      do {
         k++;
         // Generate uniform random number u in [0,1] and let p ← p × u.
         p += Math.log(r.nextDouble());
      } while (p >= L);
      return (k - 1);
   }
    
    
    
    
            
}
