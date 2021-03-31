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
import java.util.Random;

/**
 *
 * @author benoit
 */
public class SimulationOneBead {
    
    
    
    //double minZ;
    //double maxZ;
    
    
    
    
    DataPhase_ dp;
    
    
    double photonNumber; 
    double background;
    double backSecOrder;
    int patchSize;
    double axialRange;
    
    Random random;
    
    public SimulationOneBead(DataPhase_ dp,double axialRange,double photonNumber, double background,int patchSize){
        this.backSecOrder=backSecOrder;
        this.dp=dp;
        
        random = new Random(0);
        this.axialRange=axialRange;
        this.photonNumber=photonNumber;
        this.background=background;
        this.patchSize=patchSize;
    }
    
    
    public void run(){
        
        dp.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        double position=spsfc.getPosition();
        
        
        
        
        double [][] background=new double[patchSize][patchSize];
        
        int cameraFlip=+1;//  1 for cam 1   //   -1 for cam 2

        int nb=18;
        
        double theposition=.9*(axialRange/2.);
        
        
        double [] zposition = new double [nb];
        zposition[0]=position+(theposition)*cameraFlip;
        zposition[1]=position+(theposition)*cameraFlip;
        zposition[2]=position+(theposition)*cameraFlip;
        zposition[3]=position+(theposition)*cameraFlip;
        zposition[4]=position+(-theposition)*cameraFlip;
        zposition[5]=position+(-theposition)*cameraFlip;
        zposition[6]=position+(-theposition)*cameraFlip;
        zposition[7]=position+(-theposition)*cameraFlip;
        
        zposition[8]=position;
        zposition[9]=position;
        
        zposition[10]=position+(-theposition)*cameraFlip;
        zposition[11]=position+(-theposition)*cameraFlip;
        zposition[12]=position+(-theposition)*cameraFlip;
        zposition[13]=position+(-theposition)*cameraFlip;
        zposition[14]=position+(theposition)*cameraFlip;
        zposition[15]=position+(theposition)*cameraFlip;
        zposition[16]=position+(theposition)*cameraFlip;
        zposition[17]=position+(theposition)*cameraFlip;
        
        
        
        int [] xshift = new int [nb];
        xshift[0]=patchSize;
        xshift[1]=patchSize;
        xshift[2]=7*patchSize;
        xshift[3]=7*patchSize;
        xshift[4]=patchSize;
        xshift[5]=4*patchSize;
        xshift[6]=7*patchSize;
        xshift[7]=4*patchSize;
        xshift[8]=4*patchSize;//middle
        xshift[9]=4*patchSize;//close to middle
        
        xshift[10]=2*patchSize;
        xshift[11]=2*patchSize;
        xshift[12]=6*patchSize;
        xshift[13]=6*patchSize;
        xshift[14]=2*patchSize;
        xshift[15]=4*patchSize;
        xshift[16]=6*patchSize;
        xshift[17]=4*patchSize;
        
        int [] yshift = new int [nb];
        yshift[0]=patchSize;
        yshift[1]=7*patchSize;
        yshift[2]=patchSize;
        yshift[3]=7*patchSize;
        yshift[4]=4*patchSize;
        yshift[5]=patchSize;
        yshift[6]=4*patchSize;
        yshift[7]=7*patchSize;
        yshift[8]=4*patchSize;//middle
        yshift[9]=4*patchSize;//close to middle
        
        yshift[10]=2*patchSize;
        yshift[11]=6*patchSize;
        yshift[12]=2*patchSize;
        yshift[13]=6*patchSize;
        yshift[14]=4*patchSize;
        yshift[15]=2*patchSize;
        yshift[16]=4*patchSize;
        yshift[17]=6*patchSize;
        
        int nbDistance=100;
        int [] dist = new int [nbDistance];
        for (int i=0;i<nbDistance;i++){
            dist[i]=i+1;
        }
        
        
        
        if (nb>0){
            
            
            double [][] psf;
            float [][] psfmany;
            
            double [][][] imFull = new double [nbDistance][patchSize*9][patchSize*9];
            double [][][] imFullNoise = new double [nbDistance][patchSize*9][patchSize*9];
            for (int ic=0;ic<imFull.length;ic++){
                for (int i=0;i<imFull[ic].length;i++){
                    for (int ii=0;ii<imFull[ic][i].length;ii++){
                        imFull[ic][i][ii]=this.background;
                    }
                }
            }
            
            double [][][] im = new double [nb*nbDistance][patchSize][patchSize];
            double [][][] immany = new double [nb*nbDistance][patchSize][patchSize];
            //double [][][] imnoise = new double [nb*nbDistance][patchSize][patchSize];
            dp.setMany(1);
            
            double [] xposit=new double [1];
            double [] yposit=new double [1];
            double [] zoilposit=new double [1];
            double [] zposit=new double [1];
            
            
            
            for (int k=0;k<nb;k++){
                
                
                
                
            //for (double u=minZ,uu=-.0,uuu=-.0;k<nb;u+=stepZ,uu+=stepX,uuu+=stepY,k++){
                IJ.showProgress(k, nb);
                
                xposit[0]=0;
                yposit[0]=0;
                zoilposit[0]=this.dp.param.Zfocus;
                zposit[0]=zposition[k];
                dp.psf_fMany.computePSF(xposit, yposit, zoilposit, zposit);
                dp.psf.computePSF(xposit[0], yposit[0], zoilposit[0],zposit[0]);
                psf=dp.psf.getPSF();
                psfmany=dp.psf_fMany.getPSF(0);
                
                if ((patchSize>=psf.length)&&(patchSize>=psf[0].length)){
                    for (int y=0;y<nbDistance;y++){
                        
                        
                        double zum=(axialRange)*(((double)y/(double)nbDistance)-.5);
                        IJ.log("zum "+zum);
                
                
                
                        int decX=0;
                        int decY=0;
                        if ((patchSize>psf.length)&&(patchSize>psf[0].length)){
                            decX=(patchSize-psf.length)/2;
                            decY=(patchSize-psf[0].length)/2;
                            for (int z=0;z<patchSize;z++){
                                for (int zz=0;zz<patchSize;zz++){
                                    im[k*nbDistance+y][z][zz]=background[z][zz];
                                    immany[k*nbDistance+y][z][zz]=background[z][zz];
                                    //imnoise[k*nbDistance+y][z][zz]=background[z][zz];
                                }
                            }
                        }
                        
                        
                        
                        
                        if (k==8){
                            //double randZ=radius*random.nextDouble();
                            dp.psf.computePSF(xposit[0], yposit[0], zoilposit[0],position+zum*cameraFlip);
                            psf=dp.psf.getPSF();
                        }
                        
                        for (int z=0;z<psf.length;z++){
                            for (int zz=0;zz<psf[z].length;zz++){
                                if (k!=9){
                                    
                                    imFull[y][xshift[k]+z+decX][yshift[k]+zz+decY]+=(double)psf[z][zz]*photonNumber;
                                    
                                }
                            }
                        }
                    }
                }
                else{
                    IJ.log("problem patch size image simulation");
                }
                    
                

                
            }
            IJ.log("adding noise... please wait");
            for (int ic=0;ic<imFull.length;ic++){
                IJ.log(""+(imFull.length-ic)+" ...");
                for (int i=0;i<imFull[ic].length;i++){
                    for (int ii=0;ii<imFull[ic][i].length;ii++){
                        imFullNoise[ic][i][ii]=this.poissonValueLog(imFull[ic][i][ii]);
                    }
                }
            }
            IJ.log("top 4");
            ImageShow.imshow(im,"PSF/woNoise/cam");
            ImageShow.imshow(imFull,"image/woNoise/cam");
            ImageShow.imshow(imFullNoise,"image/wNoise/cam");
            
            

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
