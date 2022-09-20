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
public class FSC {
    
    double [][][] A,B;
    int w,h,d;
    double pixSizeXY, pixSizeZ;
    public FSC(double [][][] A, double [][][] B,double pixSizeXY,double pixSizeZ){
        
        if ((A.length!=B.length)||(A[0].length!=B[0].length)||(A[0][0].length!=B[0][0].length)){
            IJ.log("Sorry but the 2 images should have the same size");
            return;
        }
        
        this.A=A;
        this.B=B;
        this.d=A.length;
        this.w=A[0].length;
        this.h=A[0][0].length;
        this.pixSizeXY=pixSizeXY;
        this.pixSizeZ=pixSizeZ;
    }
    
    
    
    
    private void addMatriceTo(double [][][] mat, double [][][] toAdd){
        for (int i=0;i<mat.length;i++){
            for (int ii=0;ii<mat[0].length;ii++){
                for (int iii=0;iii<mat[0][0].length;iii++){
                    mat[i][ii][iii]+=toAdd[i][ii][iii];
                }
            }
        }
    }
    
    private void setMatriceTo(double [][][] mat, double [][][] toSet){
        for (int i=0;i<mat.length;i++){
            for (int ii=0;ii<mat[0].length;ii++){
                for (int iii=0;iii<mat[0][0].length;iii++){
                    mat[i][ii][iii]=toSet[i][ii][iii];
                }
            }
        }
    }
    
    
    
    private void shift3D(double[][][] output,double [][][] input){
        int d=input.length;
        int w=input[0].length;
        int h=input[0][0].length;
        if ((w%2==0)&&(h%2==0)&&(d%2==0)){
            for (int i=0;i<d;i++){
                for (int ii=0;ii<w;ii++){
                    for (int iii=0;iii<h;iii++){
                        output[i][ii][iii]=input[(i+d/2)%d][(ii+w/2)%w][(iii+h/2)%h];
                    }
                }
            }
        }
        else{
            return;
        }
    }
    
    
    
    private double [][] makeDistanceVector(int size){
        int d=size;
        int w=size;
        int h=size;
        double cx=w/2;
        double cy=h/2;
        double cz=d/2;
        double [][] v = new double [d*h*w][4];//vector of dist/x/y/z
        int k=0;
        for (int i=0;i<d;i++){
            for (int ii=0;ii<w;ii++){
                for (int iii=0;iii<h;iii++){
                    v[k][0]=Math.sqrt((double)((i-cz)*(i-cz)+(ii-cx)*(ii-cx)+(iii-cy)*(iii-cy)));
                    v[k][1]=i;
                    v[k][2]=ii;
                    v[k][3]=iii;
                    k++;
                }
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
        
        int size=256;
        IJ.log("SIZE "+size);
        double [][][] weight = new double [d][w][h];
        double xtmp,ytmp,ztmp;
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                for (int iii=0;iii<d;iii++){
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
                    if ((iii<d*1./8.)||(iii>=d*7./8.)){
                        ztmp=Math.sin(4*Math.PI*(double)iii/(double)d);
                        ztmp*=ztmp;
                    }
                    else{
                        ztmp=1;
                    }
                    weight[iii][i][ii]=xtmp*ytmp*ztmp;
                }
            }
        }
        
        
        
        int splitX=(int)Math.ceil((double)w/(double)(size/2));
        int splitY=(int)Math.ceil((double)h/(double)(size/2));
        int splitZ=(int)Math.ceil((double)d/(double)(size/2));
        double [][][] A = new double [size][size][size];
        double [][][] B = new double [size][size][size];
        double [][][] C = new double [size][size][size];
        double [][][] D = new double [size][size][size];
        double [][][] E = new double [size][size][size];
        double [][][] F = new double [size][size][size];
        double [][][] G = new double [size][size][size];
        double [][][] H = new double [size][size][size];
        double [][][] numeratorI = new double [size][size][size];
        double [][][] tmpointer;
        FastFourierTransform fft1 = new FastFourierTransform(size,size,size);
        FastFourierTransform fft2 = new FastFourierTransform(size,size,size);
        
        loop:for (int x=0;x<splitX;x++){
            
            for (int y=0;y<splitY;y++){
                IJ.showProgress((double)(x*splitY+y)/(double)(splitX*splitY));
                for (int z=0;z<splitZ;z++){
                    //Image copy including padding
                    for (int i=0;i<size/2;i++){
                        for (int ii=0;ii<size/2;ii++){
                            for (int iii=0;iii<size/2;iii++){
                                if (((size/2)*x+i<w)&&((size/2)*y+ii<h)&&((size/2)*z+iii<d)){
                                    A[iii+size/4][i+size/4][ii+size/4]=this.A[(size/2)*z+iii][(size/2)*x+i][(size/2)*y+ii]*weight[(size/2)*z+iii][(size/2)*x+i][(size/2)*y+ii];
                                    B[iii+size/4][i+size/4][ii+size/4]=this.B[(size/2)*z+iii][(size/2)*x+i][(size/2)*y+ii]*weight[(size/2)*z+iii][(size/2)*x+i][(size/2)*y+ii];
                                }
                            }
                        }
                    }
                    
                    fft1.setReal(A);
                    fft1.setImag(C);
                    
                    fft2.setReal(B);
                    fft2.setImag(D);
                    
                    fft1.fft3D();
                    fft2.fft3D();
                    
                    
                    tmpointer=fft1.getPointerRealOut3D();
                    addMatriceTo(E,tmpointer);
                    tmpointer=fft1.getPointerImagOut3D();
                    addMatriceTo(F,tmpointer);
                    
                    tmpointer=fft2.getPointerRealOut3D();
                    addMatriceTo(G,tmpointer);
                    tmpointer=fft2.getPointerImagOut3D();
                    addMatriceTo(H,tmpointer);
                    
                }
            }
        }
        
        
        
        
        IJ.log("FFT OK");
        
        
        shift3D(A,E);
        shift3D(B,F);
        
        IJ.log("SHIFT image 1 OK");
        
        shift3D(C,G);
        shift3D(D,H);
        
        IJ.log("SHIFT image 2 OK");
        
        
        
        

        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                for (int iii=0;iii<size;iii++){
                    E[iii][i][ii]=(A[iii][i][ii]*A[iii][i][ii]+B[iii][i][ii]*B[iii][i][ii]);
                    F[iii][i][ii]=(C[iii][i][ii]*C[iii][i][ii]+D[iii][i][ii]*D[iii][i][ii]);
                    G[iii][i][ii]=(A[iii][i][ii]*C[iii][i][ii]+B[iii][i][ii]*D[iii][i][ii]);
                    H[iii][i][ii]=(B[iii][i][ii]*C[iii][i][ii]-A[iii][i][ii]*D[iii][i][ii]);
                    numeratorI[iii][i][ii]=(G[iii][i][ii]*G[iii][i][ii]+H[iii][i][ii]*H[iii][i][ii]);
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
        
        
        
        double step=2;//pixel step to convert in vector
        int sizeVect=(int)Math.ceil(size*Math.sqrt(3)/(2*step));//sqrt(3) because 3D
        double [] vectE = new double[sizeVect];
        double [] vectF = new double[sizeVect];
        double [] vectG = new double[sizeVect];
        double [] vectH = new double[sizeVect];
        double [] numerator = new double[sizeVect];
        double [] vectN = new double[sizeVect];//number of pixels
        
        
        //sum over distances to convert 3D images to vectors
        for (int i=0;i<v.length;i++){
            int bin=(int)(v[i][0]/step);
            
            
            vectE[bin]+=E[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]];
            vectF[bin]+=F[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]];
            vectG[bin]+=G[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]];
            vectH[bin]+=H[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]];
            vectN[bin]++;
            numerator[bin]+=G[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]]*G[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]]+H[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]]*H[(int)v[i][1]][(int)v[i][2]][(int)v[i][3]];
        }
        
        
        
        
        
        double [] vectFinal = new double[sizeVect];
        double [] xaxisXY = new double[sizeVect];
        double [] xaxisZ = new double[sizeVect];
        for (int i=0;i<sizeVect;i++){
            vectFinal[i]=(Math.sqrt(numerator[i])/Math.sqrt(vectE[i]*vectF[i]));// - 2/Math.sqrt(vectN[i]/2);
            xaxisXY[i]=i*step/((size)*this.pixSizeXY);//(size/2) because we use half of image ?
            xaxisZ[i]=i*step/((size)*this.pixSizeZ);
            
        }
        
        
        
        int xresol=0;
        while (vectFinal[xresol]>(1./7.)){
            xresol++;
        }
        IJ.log("r "+xaxisXY[xresol]+"  "+xresol);
        double resolutionXY=xaxisXY[xresol];
        double resolutionZ=xaxisZ[xresol];
        
        
        plot(xaxisXY,vectFinal,xaxisXY[xresol],vectFinal[xresol],"FSC","1/nm","Lateral FSC");
        
        
        IJ.showProgress(1);
        
        IJ.log("Resolution XY: "+(1/resolutionXY));
        IJ.log("Resolution Z: "+(1/resolutionZ));
        
        
        //well done
        
        
    }
    
    
    
    
    
    
    public void plot(double [] x, double [] y,double lineX,double lineY,String title,String xlabel,String ylabel){
        
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
        
        p.show();
        
    }
    
    
    
    
    
}
