/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.pasteur.imagej.data.StackLocalization;
import org.pasteur.imagej.utils.FastFourierTransform;
import org.pasteur.imagej.utils.ImageShow;

/**
 *
 * @author benoit
 */
public class RegistrationCrossCorrel {
    
    StackLocalization sl1;
    StackLocalization sl2;
    int pixelSizeNM;
    public RegistrationCrossCorrel(StackLocalization sl1,StackLocalization sl2,int pixelSizeNM){
        this.sl1=sl1;
        this.sl2=sl2;
        this.pixelSizeNM=pixelSizeNM;
        
    }
    
    //return shift X and Y
    public double [] getShift(){
        
        
        double maxX=Double.NEGATIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        double minZ=Double.NEGATIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        for (int i=0;i<sl1.fl.size();i++){
            for (int j=0;j<sl1.fl.get(i).loc.size();j++){
                if (sl1.fl.get(i).loc.get(j).X>maxX){
                    maxX=sl1.fl.get(i).loc.get(j).X;
                }
                if (sl1.fl.get(i).loc.get(j).Y>maxY){
                    maxY=sl1.fl.get(i).loc.get(j).Y;
                }
            }
        }
        for (int i=0;i<sl2.fl.size();i++){
            for (int j=0;j<sl2.fl.get(i).loc.size();j++){
                if (sl2.fl.get(i).loc.get(j).X>maxX){
                    maxX=sl2.fl.get(i).loc.get(j).X;
                }
                if (sl2.fl.get(i).loc.get(j).Y>maxY){
                    maxY=sl2.fl.get(i).loc.get(j).Y;
                }
            }
        }
        
        double [][] im1=makePaddedImage(sl1,pixelSizeNM,maxX,maxY,1);
        double [][] im2=makePaddedImage(sl2,pixelSizeNM,maxX,maxY,1);
        int w=im1.length;
        int h=im1[0].length;
        double [][] tmp1a=new double [w][h];
        double [][] tmp2a=new double [w][h];
        double [][] tmp1b=new double [w][h];
        double [][] tmp2b=new double [w][h];
        double [][] tmp1=new double [w][h];
        double [][] tmp2=new double [w][h];
        
        FastFourierTransform fft = new FastFourierTransform(im1.length,im1[0].length);
        fft.setReal(im1);
        fft.setImag(tmp1a);
        fft.fft2D();
        double [][] pt1=fft.getPointerRealOut2D();
        double [][] pt2=fft.getPointerImagOut2D();
        for (int i=0;i<w;i++){
            for (int j=0;j<h;j++){
                tmp1a[i][j]=pt1[i][j];
                tmp2a[i][j]=pt2[i][j];
            }
        }
        
        
        fft.setReal(im2);
        fft.setImag(tmp1b);
        fft.fft2D();
        pt1=fft.getPointerRealOut2D();
        pt2=fft.getPointerImagOut2D();
        for (int i=0;i<w;i++){
            for (int j=0;j<h;j++){
                tmp1b[i][j]=pt1[i][j];
                tmp2b[i][j]=pt2[i][j];
            }
        }
        
        
        for (int i=0;i<w;i++){
            for (int j=0;j<h;j++){
                tmp1[i][j] = tmp1a[i][j]*tmp1b[i][j]+tmp2a[i][j]*tmp2b[i][j];
                tmp2[i][j] = -tmp1a[i][j]*tmp2b[i][j]+tmp2a[i][j]*tmp1b[i][j];
            }
        }
        
        
        
        
        fft.setReal(tmp1);
        fft.setImag(tmp2);
        
        fft.ifft2D();
        
        
        pt1=fft.getPointerRealOut2D();
        for (int i=0;i<w;i++){
            for (int j=0;j<h;j++){
                tmp1[i][j]=pt1[i][j];
            }
        }
        for (int i=0;i<w;i++){
            for (int j=0;j<h;j++){
                tmp2[i][j] = tmp1[(i+w/2)%w][(j+h/2)%h];
            }
        }
        //ImageShow.imshow(tmp2);
        
        int maxPositX=0;
        int maxPositY=0;
        double maxi=tmp2[0][0];
        for (int i=0;i<w;i++){
            for (int j=0;j<h;j++){
                if (tmp2[i][j]>maxi){
                    maxi=tmp2[i][j];
                    maxPositX=i;
                    maxPositY=j;
                }
            }
        }
        
        
        double topX=topOfParabola(-1.*pixelSizeNM,tmp2[maxPositX-1][maxPositY],0,tmp2[maxPositX][maxPositY],1.*pixelSizeNM,tmp2[maxPositX+1][maxPositY]);
        double topY=topOfParabola(-1.*pixelSizeNM,tmp2[maxPositX][maxPositY-1],0,tmp2[maxPositX][maxPositY],1.*pixelSizeNM,tmp2[maxPositX][maxPositY+1]);
        
        double shiftX=(w/2)*pixelSizeNM-(maxPositX*pixelSizeNM+topX);
        double shiftY=(h/2)*pixelSizeNM-(maxPositY*pixelSizeNM+topY);
        IJ.log("Lateral of cam 2: "+shiftX+"  "+shiftY);
        double [] result = new double[2];
        result[0]=shiftX;
        result[1]=shiftY;
        return result;
        
    }
    
    
    
    
    
    
    


    
    
    private double topOfParabola(double xa,double ya,double xb,double yb,double xc,double yc){
        double a=(yc-ya)/((xc-xa)*(xc-xb))-(yb-ya)/((xb-xa)*(xc-xb));
        double b=((yb-ya)/(xb-xa))-a*(xb+xa);
        return(-b/(2*a));
    }
    
    
    
    
    double [][] makePaddedImage(StackLocalization sl,double pixelsizeNM,double maxX,double maxY,int shift){
        
        
        
        
        
        double binColor=255;
        int width=(int)Math.ceil((maxX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY)/pixelsizeNM);
        int depth=1;
        width+=shift*2;
        height+=shift*2;
        int nbShift=1+shift*2;
        int [][] shifT = new int[1+shift*2][1+shift*2];
        
        //create first element
        shifT[0][0]=1;
        //create first line
        for (int i=1;i<=shift;i++){
            shifT[0][i]=shifT[0][i-1]+1;
        }
        for (int i=shift+1;i<1+shift*2;i++){
            shifT[0][i]=shifT[0][i-1]-1;
        }
        //create each column
        for (int ii=0;ii<nbShift;ii++){
            for (int i=1;i<=shift;i++){
                shifT[i][ii]=shifT[i-1][ii]+shifT[0][ii];
            }
        }
        for (int ii=0;ii<nbShift;ii++){
            for (int i=shift+1;i<1+shift*2;i++){
                shifT[i][ii]=shifT[i-1][ii]-shifT[0][ii];
            }
        }
        
        
        double x;
        double y;
        double z;
        
        width=2*(int)Math.pow(2,(Math.ceil(Math.log(Math.max(width, 1))/Math.log(2))));
        height=2*(int)Math.pow(2,(Math.ceil(Math.log(Math.max(height, 1))/Math.log(2))));
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        
        
        double [][] image = new double [width][height];
        
        
        int xx,yy;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    xx=(int)((x)/pixelsizeNM)+shift;
                    yy=(int)((y)/pixelsizeNM)+shift;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if ((xx-shift>=0)&&(yy-shift>=0)&&(xx+shift<width)&&(yy+shift<height)){
                        for (int a=-shift,aa=0;a<=shift;a++,aa++){
                            for (int b=-shift,bb=0;b<=shift;b++,bb++){
                                if ((xx+a+width/4>=0)&&(xx+a+width/4<width)&&(yy+b+height/4>=0)&&(yy+b+height/4<height)){
                                    image[xx+a+width/4][yy+b+height/4]+=shifT[aa][bb];
                                }
                            }
                        }
                    }
                }
            }
        }
        return image;
    }
    
}
