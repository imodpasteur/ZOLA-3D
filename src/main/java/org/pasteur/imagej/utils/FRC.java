/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

/**
 *
 * @author benoit
 */
public class FRC {
    
    double [][] A,B;
    int w,h;
    double pixSizeXY;
    public FRC(double [][]A, double [][]B,double pixSizeXY){
        
        if ((A.length!=B.length)||(A[0].length!=B[0].length)){
            IJ.log("Sorry but the 2 images should have the same size");
            return;
        }
        
        this.A=A;
        this.B=B;
        this.w=A.length;
        this.h=A[0].length;
        this.pixSizeXY=pixSizeXY;
    }
    
    
    
    
    private void addMatriceTo(double [][] mat, double [][] toAdd){
        for (int i=0;i<mat.length;i++){
            for (int ii=0;ii<mat[0].length;ii++){
                    mat[i][ii]+=toAdd[i][ii];
            }
        }
    }
    
    private void setMatriceTo(double [][] mat, double [][] toSet){
        for (int i=0;i<mat.length;i++){
            for (int ii=0;ii<mat[0].length;ii++){
                    mat[i][ii]=toSet[i][ii];
            }
        }
    }
    
    
    
    private void shift2D(double[][] output,double [][] input){
        int w=input.length;
        int h=input[0].length;
        if ((w%2==0)&&(h%2==0)){
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                        output[i][ii]=input[(i+w/2)%w][(ii+h/2)%h];
                }
            }
        }
        else{
            return;
        }
    }
    
    
    
    private double [][] makeDistanceVector(int size){
        int w=size;
        int h=size;
        double cx=w/2;
        double cy=h/2;
        double [][] v = new double [h*w][3];//vector of dist/x/y
        int k=0;
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                v[k][0]=Math.sqrt((double)((i-cx)*(i-cx)+(ii-cy)*(ii-cy)));
                v[k][1]=i;
                v[k][2]=ii;
                k++;
            }
        }
        
        
        Arrays.sort(v, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o1[0]).compareTo(o2[0]);
            }
        });
        
        return v;
        
    }
    
    
    public void run(){
        
        int size=1024;
        IJ.log("SIZE "+size);
        double [][] weight = new double [w][h];
        double xtmp,ytmp;
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                if ((i<w*1./8.)||(i>=w*7./8.)){
                    xtmp=Math.sin(4*Math.PI*(double)i/(double)w);
                    xtmp*=xtmp;
                }
                else{
                    xtmp=1;
                }
                if ((ii<h*1./8.)||(ii>=h*7./8.)){
                    ytmp=Math.sin(4*Math.PI*(double)ii/(double)h);
                    ytmp*=ytmp;
                }
                else{
                    ytmp=1;
                }
                
                weight[i][ii]=1;//xtmp*ytmp*ztmp;
                
            }
        }
        
        
        
        int splitX=(int)Math.ceil((double)w/(double)(size/2));
        int splitY=(int)Math.ceil((double)h/(double)(size/2));
        double [][] A = new double [size][size];
        double [][] B = new double [size][size];
        double [][] C = new double [size][size];
        double [][] D = new double [size][size];
        double [][] E = new double [size][size];
        double [][] F = new double [size][size];
        double [][] G = new double [size][size];
        double [][] H = new double [size][size];
        double [][] tmpointer;
        FastFourierTransform fft1 = new FastFourierTransform(size,size);
        FastFourierTransform fft2 = new FastFourierTransform(size,size);
        
        loop:for (int x=0;x<splitX;x++){
            
            for (int y=0;y<splitY;y++){
                IJ.showProgress((double)(x*splitY+y)/(double)(splitX*splitY));
                //Image copy including padding
                for (int i=0;i<size/2;i++){
                    for (int ii=0;ii<size/2;ii++){
                        if (((size/2)*x+i<w)&&((size/2)*y+ii<h)){
                            A[i+size/4][ii+size/4]=this.A[(size/2)*x+i][(size/2)*y+ii]*weight[(size/2)*x+i][(size/2)*y+ii];
                            B[i+size/4][ii+size/4]=this.B[(size/2)*x+i][(size/2)*y+ii]*weight[(size/2)*x+i][(size/2)*y+ii];
                        }
                        
                    }
                }

                fft1.setReal(A);
                fft1.setImag(C);

                fft2.setReal(B);
                fft2.setImag(D);

                fft1.fft2D();
                fft2.fft2D();


                tmpointer=fft1.getPointerRealOut2D();
                addMatriceTo(E,tmpointer);
                tmpointer=fft1.getPointerRealOut2D();
                addMatriceTo(F,tmpointer);

                tmpointer=fft2.getPointerRealOut2D();
                addMatriceTo(G,tmpointer);
                tmpointer=fft2.getPointerRealOut2D();
                addMatriceTo(H,tmpointer);
                    
                
            }
        }
        
        
        
        
        IJ.log("FFT OK");
        
        
        shift2D(A,E);
        shift2D(B,F);
        
        IJ.log("SHIFT image 1 OK");
        
        shift2D(C,G);
        shift2D(D,H);
        
        IJ.log("SHIFT image 2 OK");
        
        
        
        

        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                for (int iii=0;iii<size;iii++){
                    E[i][ii]=(A[i][ii]*A[i][ii]+B[i][ii]*B[i][ii]);
                    F[i][ii]=(C[i][ii]*C[i][ii]+D[i][ii]*D[i][ii]);
                    G[i][ii]=(A[i][ii]*C[i][ii]+B[i][ii]*D[i][ii]);
                    H[i][ii]=(B[i][ii]*C[i][ii]-A[i][ii]*D[i][ii]);
                }
            }
        }
        
        IJ.showProgress(.99);
        
        IJ.log("NORMALIZED CORRELATION OK");
//        ImageShow.imshow(A,"A");
//        ImageShow.imshow(B,"B");
//        ImageShow.imshow(C,"C");
//        ImageShow.imshow(D,"D");
//        
//        ImageShow.imshow(E,"E");
//        ImageShow.imshow(F,"F");
//        ImageShow.imshow(G,"G");
//        ImageShow.imshow(H,"H");
        
        
        double [][] v=makeDistanceVector(size);
        //A contains distance matrix
        IJ.log("DIST MATRIX OK");
        
        
        
        double step=1;//pixel step to convert in vector
        int sizeVect=(int)Math.ceil(size*Math.sqrt(2)/(2*step));//sqrt(2) because 2D
        double [] vectE = new double[sizeVect];
        double [] vectF = new double[sizeVect];
        double [] vectG = new double[sizeVect];
        double [] vectH = new double[sizeVect];
        double [] vectN = new double[sizeVect];//number of pixels
        
        
        //sum over distances to convert 2D images to vectors
        for (int i=0;i<v.length;i++){
            int bin=(int)(v[i][0]/step);
            
            
            vectE[bin]+=E[(int)v[i][1]][(int)v[i][2]];
            vectF[bin]+=F[(int)v[i][1]][(int)v[i][2]];
            vectG[bin]+=G[(int)v[i][1]][(int)v[i][2]];
            vectH[bin]+=H[(int)v[i][1]][(int)v[i][2]];
            vectN[bin]++;
            
        }
        
        
        
        
        
        double [] vectFinal = new double[sizeVect];
        double [] xaxisXY = new double[sizeVect];
        for (int i=0;i<sizeVect;i++){
            vectFinal[i]=(vectG[i]/Math.sqrt(vectE[i]*vectF[i]));// - 2/Math.sqrt(vectN[i]/2);
            xaxisXY[i]=i*step/(size*this.pixSizeXY);
        }
        
        double [] vectFinalSmooth = new double[sizeVect];
        
        double m,c;
        for (int i=0;i<sizeVect;i++){
            m=0;
            c=0;
            for (int u=-5;u<=5;u++){
                if ((u+i>=0)&&(u+i<sizeVect)){
                    m+=vectFinal[i+u];
                    c++;
                }
            }
            vectFinalSmooth[i]=m/c;
        }
        
        double threshold=1./7.;
        
        int xresol=0;
        while (vectFinalSmooth[xresol]>threshold){
            xresol++;
        }
        IJ.log("r "+xaxisXY[xresol]+"  "+xresol);
        double resolutionXY=xaxisXY[xresol];
        
        
        plot(xaxisXY,vectFinal,vectFinalSmooth,xaxisXY[xresol],threshold,"FRC","1/nm","FRC");
        
        
        IJ.showProgress(1);
        
        IJ.log("Resolution XY: "+(1/resolutionXY));
        
        
        //well done
        
        
    }
    
    
    
    
    
    
    public void plot(double [] x, double [] y, double [] y2,double lineX,double lineY,String title,String xlabel,String ylabel){
        
        double xMin=x[0]; double xMax=x[0]; double yMin=y[0]; double yMax=y[0];
        
        for (int i=0;i<y.length;i++){
            if (x[i]<xMin){
                xMin=x[i];
            }
            if (x[i]>xMax){
                xMax=x[i];
            }
            if (y[i]<yMin){
                yMin=y[i];
            }
            if (y[i]>yMax){
                yMax=y[i];
            }
            
        }
        
        
        Plot p = new Plot(""+title,xlabel,ylabel);
        p.setLimits(xMin-.01*(xMax-xMin), xMax+.01*(xMax-xMin), yMin-.01*(yMax-yMin), yMax+.01*(yMax-yMin));
        p.setLineWidth(2);
        p.setFont(0, 18);
        for (int ii=0;ii<x.length;ii++){
            
            p.setColor(Color.blue);
            p.add("CIRCLE",  x,y);
            //p.show();


            
        }
        p.setColor(Color.red);
        p.drawLine(0, lineY, lineX, lineY);
        p.drawLine(lineX, 0, lineX, lineY);
        p.addPoints(x, y2,Plot.LINE);
        p.show();
        
    }
    
    
    
    
    
}
