/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

/**
 *
 * @author benoit
 */
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 *
 * @author benoit
 */
public class FourierTransform {
    int nbThread=1;//no thread....
    public FourierTransform(int nothing){
        
    }
    
    public FourierTransform(){
        
    }
    public void dft(double [][] realIn,double [][] imagIn,double [][] realOut,double [][] imagOut){
        
        int w=realIn.length;
        int h=realIn[0].length;
        DFT [] dft = new DFT[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            dft[i]=new DFT(realIn,imagIn,realOut,imagOut,start,end,true);
            dft[i].run();
        }
        
        
        
    }
    
    
    
    
    

    
    
    
    public void fft(float [][][] rIn,float [][][] iIn,float [][][] realOut,float [][][] imagOut){
        float [][][] realIn = new float [rIn.length][rIn[0].length][rIn[0][0].length];
        float [][][] imagIn = new float [iIn.length][iIn[0].length][iIn[0][0].length];
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                for (int iii=0;iii<rIn[0][0].length;iii++){
                    realIn[i][ii][iii]=rIn[i][ii][iii];
                    imagIn[i][ii][iii]=iIn[i][ii][iii];
                }
            }
        }
        int w=realIn.length;
        int h=realIn[0].length;
        int d=realIn[0][0].length;
        FastFT3D [] fft = new FastFT3D[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            fft[i]=new FastFT3D(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].run();
        }
        
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        
        step=h/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=h;
            }
            fft[i]=new FastFT3D(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].changeSecondStep();
            fft[i].run();
        }
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            fft[i]=new FastFT3D(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].changeThirdStep();
            fft[i].run();
        }
        
        
    }
    
    
    
    
    
    
    public void ifft(float [][][] rIn,float [][][] iIn,float [][][] realOut,float [][][] imagOut){
        float [][][] realIn = new float [rIn.length][rIn[0].length][rIn[0][0].length];
        float [][][] imagIn = new float [iIn.length][iIn[0].length][iIn[0][0].length];
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                for (int iii=0;iii<rIn[0][0].length;iii++){
                    realIn[i][ii][iii]=rIn[i][ii][iii];
                    imagIn[i][ii][iii]=iIn[i][ii][iii];
                }
            }
        }
        int w=realIn.length;
        int h=realIn[0].length;
        int d=realIn[0][0].length;
        FastFT3D [] fft = new FastFT3D[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            
            //IJ.log("i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            
            fft[i]=new FastFT3D(realIn,imagIn,realOut,imagOut,start,end,false);
            fft[i].run();
        }
        
        
        
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        
        step=h/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=h;
            }
            fft[i]=new FastFT3D(realIn,imagIn,realOut,imagOut,start,end,false);
            fft[i].changeSecondStep();
            fft[i].run();
        }
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            fft[i]=new FastFT3D(realIn,imagIn,realOut,imagOut,start,end,false);
            fft[i].changeThirdStep();
            fft[i].run();
        }
        
        
    }
    
    
    
    
    
    
    
    public void fft(double [][] rIn,double [][] iIn,double [][] realOut,double [][] imagOut){
        double [][] realIn = new double [rIn.length][rIn[0].length];
        double [][] imagIn = new double [iIn.length][iIn[0].length];
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                realIn[i][ii]=rIn[i][ii];
                imagIn[i][ii]=iIn[i][ii];
            }
        }
        int w=realIn.length;
        int h=realIn[0].length;
        FastFT [] fft = new FastFT[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            //IJ.log("i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            fft[i]=new FastFT(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].run();
        }
        
        
        
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        step=h/nbThread;
        for (int i=0;i<nbThread;i++){
            //IJ.log("next i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=h;
            }
            fft[i]=new FastFT(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].changeSecondStep();
            fft[i].run();
        }
        //IJ.log("end ");
        
        
    }
    
    
    
    public void ffttest(double [][] rIn,double [][] iIn,double [][] realOut,double [][] imagOut){
        double [][] realIn = new double [rIn.length][rIn[0].length];
        double [][] imagIn = new double [iIn.length][iIn[0].length];
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                realIn[i][ii]=rIn[i][ii];
                imagIn[i][ii]=iIn[i][ii];
            }
        }
        int w=realIn.length;
        int h=realIn[0].length;
        FastFT [] fft = new FastFT[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            //IJ.log("i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            fft[i]=new FastFT(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].run();
        }
        
        
        
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        step=h/nbThread;
        for (int i=0;i<nbThread;i++){
            //IJ.log("next i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=h;
            }
            fft[i]=new FastFT(realIn,imagIn,realOut,imagOut,start,end,true);
            fft[i].changeSecondStep();
            fft[i].run();
        }
        //IJ.log("end ");
        
        
    }
    
    
    public void ifft(double [][] rIn,double [][] iIn,double [][] realOut,double [][] imagOut){
        double [][] realIn = new double [rIn.length][rIn[0].length];
        double [][] imagIn = new double [iIn.length][iIn[0].length];
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                realIn[i][ii]=rIn[i][ii];
                imagIn[i][ii]=iIn[i][ii];
            }
        }
        int w=realIn.length;
        int h=realIn[0].length;
        FastFT [] fft = new FastFT[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            //IJ.log("i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            fft[i]=new FastFT(realIn,imagIn,realOut,imagOut,start,end,false);
            fft[i].run();
        }
        
        
        realIn=realOut;
        imagIn=imagOut;
        
        step=h/nbThread;
        for (int i=0;i<nbThread;i++){
            //IJ.log("next i "+i);
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=h;
            }
            fft[i]=new FastFT(realIn,imagIn,realOut,imagOut,start,end,false);
            fft[i].changeSecondStep();
            fft[i].run();
        }
        //IJ.log("end ");
        
        
    }
    
    
    
    
    
    
    public void idft(double [][] realIn,double [][] imagIn,double [][] realOut,double [][] imagOut){
        
        int w=realIn.length;
        int h=realIn[0].length;
        DFT [] dft = new DFT[nbThread];
        int step=w/nbThread;
        for (int i=0;i<nbThread;i++){
            int start=i*step;
            int end=(i+1)*step;
            if (i==nbThread-1){
                end=w;
            }
            dft[i]=new DFT(realIn,imagIn,realOut,imagOut,start,end,false);
            dft[i].run();
        }
        
        
        
    }
    
    public double [][] shift2D(double[][] matrix){
        int w=matrix.length;
        int h=matrix[0].length;
        if ((w%2==0)&&(h%2==0)){
            double [][] res = new double [w][h];
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    res[i][ii]=matrix[(i+w/2)%w][(ii+h/2)%h];
                    //res[i][ii]=matrix[i][ii];
                }
            }
            return res;
        }
        else{
            IJ.log("problem image size has to be even");
            return null;
        }
    }
    
    
    
    
    
    public void shift2D(double[][] matrix,double [][] result){
        int w=matrix.length;
        int h=matrix[0].length;
        if ((w%2==0)&&(h%2==0)){
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    result[i][ii]=matrix[(i+w/2)%w][(ii+h/2)%h];
                    //res[i][ii]=matrix[i][ii];
                }
            }
        }
        else{
            IJ.log("problem image size has to be even");
            return;
        }
    }
    
    
    
    public double [][] copy(double[][] matrix){
        int w=matrix.length;
        int h=matrix[0].length;
        
        double [][] res = new double [w][h];
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                res[i][ii]=matrix[i][ii];
            }
        }
        return res;
        
    }
    
    public double [][] minus(double[][] matrix){
        int w=matrix.length;
        int h=matrix[0].length;
        
        double [][] res = new double [w][h];
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                res[i][ii]=-matrix[i][ii];
            }
        }
        return res;
        
    }
    
    
    public double [][][] shift2D(double[][][] matrix){
        int w=matrix[0].length;
        int h=matrix[0][0].length;
        if ((w%2==0)&&(h%2==0)){
            double [][][] res = new double [matrix.length][w][h];
            for (int u=0;u<matrix.length;u++){

                
                for (int i=0;i<w;i++){
                    for (int ii=0;ii<h;ii++){
                        res[u][i][ii]=matrix[u][(i+w/2)%w][(ii+h/2)%h];
                        //res[i][ii]=matrix[i][ii];
                    }
                }
                
                
            }
            return res;
        }
        else{
            IJ.log("problem image size has to be even");
            return null;
        }
    }
    
    
    
    
    
    
    public double [][][] shift3D(double[][][] matrix){
        int d=matrix.length;
        int w=matrix[0].length;
        int h=matrix[0][0].length;
        if ((w%2==0)&&(h%2==0)&&(d%2==0)){
            double [][][] res = new double [d][w][h];
            for (int u=0;u<matrix.length;u++){

                for (int i=0;i<w;i++){
                    for (int ii=0;ii<w;ii++){
                        for (int iii=0;iii<h;iii++){
                            res[i][ii][iii]=matrix[(i+d/2)%d][(ii+w/2)%w][(iii+h/2)%h];
                        }
                    }
                }
                
                
            }
            return res;
        }
        else{
            IJ.log("problem image size has to be even");
            return null;
        }
    }
    
    
    
    public void shift3D(double[][][] matrix,double [][][] result){
        int d=matrix.length;
        int w=matrix[0].length;
        int h=matrix[0][0].length;
        if ((w%2==0)&&(h%2==0)&&(d%2==0)){
            for (int u=0;u<matrix.length;u++){
                for (int i=0;i<w;i++){
                    for (int ii=0;ii<w;ii++){
                        for (int iii=0;iii<h;iii++){
                            result[i][ii][iii]=matrix[(i+d/2)%d][(ii+w/2)%w][(iii+h/2)%h];
                        }
                    }
                }
                
                
            }
        }
        else{
            IJ.log("problem image size has to be even");
        }
    }
    
    
    
    
    
    
    public float [][][] shift3D(float[][][] matrix){
        int d=matrix.length;
        int w=matrix[0].length;
        int h=matrix[0][0].length;
        if ((w%2==0)&&(h%2==0)&&(d%2==0)){
            float [][][] res = new float [d][w][h];
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    for (int iii=0;iii<d;iii++){
                        res[i][ii][iii]=matrix[(i+w/2)%w][(ii+h/2)%h][(iii+d/2)%d];
                    }
                }
            }
            return res;
        }
        else{
            IJ.log("problem image size has to be even");
            return null;
        }
    }
    
    
    
    public void shift3D(float[][][] matrix,float [][][] result){
        int w=matrix.length;
        int h=matrix[0].length;
        int d=matrix[0][0].length;
        if ((w%2==0)&&(h%2==0)&&(d%2==0)){
            for (int i=0;i<w;i++){
                for (int ii=0;ii<h;ii++){
                    for (int iii=0;iii<d;iii++){
                        result[i][ii][iii]=matrix[(i+w/2)%w][(ii+h/2)%h][(iii+d/2)%d];
                    }
                }
            }
        }
        else{
            IJ.log("problem image size has to be even");
        }
    }
    
    
    
    
    public double [][] getMagnitude(double [][] real,double [][] imag){
        
        int w=real.length;
        int h=real[0].length;
        double [][] res=new double [w][h];
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                res[i][ii]=Math.sqrt(real[i][ii]*real[i][ii]+imag[i][ii]*imag[i][ii]);
            }
        }
        return res;
    }
    
    public void getMagnitude(double [][] real,double [][] imag,double [][] result){
        
        int w=real.length;
        int h=real[0].length;
        
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                result[i][ii]=Math.sqrt(real[i][ii]*real[i][ii]+imag[i][ii]*imag[i][ii]);
            }
        }
    }
    
    public double [][] getPhase(double [][] real,double [][] imag){
        
        int w=real.length;
        int h=real[0].length;
        double [][] res=new double [w][h];
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                res[i][ii]=Math.atan2(imag[i][ii],real[i][ii]);
                /*if (real[i][ii]>0)
                    res[i][ii]=Math.atan(imag[i][ii]/real[i][ii]);
                else if (real[i][ii]<0)
                    res[i][ii]=Math.atan(imag[i][ii]/real[i][ii])+Math.PI;
                else if ((real[i][ii]<0)&&(imag[i][ii]<0))
                    res[i][ii]=Math.atan(imag[i][ii]/real[i][ii])-Math.PI;
                else if ((real[i][ii]==0)&&(imag[i][ii]>0))
                    res[i][ii]=Math.PI/2;
                else if ((real[i][ii]==0)&&(imag[i][ii]<0))
                    res[i][ii]=-Math.PI/2;
                else if ((real[i][ii]==0)&&(imag[i][ii]==0))
                    res[i][ii]=Double.NaN;*/
            }
        }
        return res;
    }
    
    
    public void getPhase(double [][] real,double [][] imag,double [][] result){
        
        int w=real.length;
        int h=real[0].length;
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                result[i][ii]=Math.atan2(imag[i][ii],real[i][ii]);
                /*if (real[i][ii]>0)
                    res[i][ii]=Math.atan(imag[i][ii]/real[i][ii]);
                else if (real[i][ii]<0)
                    res[i][ii]=Math.atan(imag[i][ii]/real[i][ii])+Math.PI;
                else if ((real[i][ii]<0)&&(imag[i][ii]<0))
                    res[i][ii]=Math.atan(imag[i][ii]/real[i][ii])-Math.PI;
                else if ((real[i][ii]==0)&&(imag[i][ii]>0))
                    res[i][ii]=Math.PI/2;
                else if ((real[i][ii]==0)&&(imag[i][ii]<0))
                    res[i][ii]=-Math.PI/2;
                else if ((real[i][ii]==0)&&(imag[i][ii]==0))
                    res[i][ii]=Double.NaN;*/
            }
        }
    }
    
    
    public double [][] getReal(double [][] magn,double [][] phase){
        
        int w=magn.length;
        int h=magn[0].length;
        double [][] res=new double [w][h];
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                res[i][ii]=Math.abs(magn[i][ii])*(Math.cos(phase[i][ii]));
            }
        }
        return res;
    }
    
    
    public void getReal(double [][] magn,double [][] phase,double [][] result){
        
        int w=magn.length;
        int h=magn[0].length;
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                result[i][ii]=Math.abs(magn[i][ii])*(Math.cos(phase[i][ii]));
            }
        }
    }
    
    
    public double [][] getImag(double [][] magn,double [][] phase){
        
        int w=magn.length;
        int h=magn[0].length;
        double [][] res=new double [w][h];
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                res[i][ii]=Math.abs(magn[i][ii])*(Math.sin(phase[i][ii]));
            }
        }
        return res;
    }
    
    public void getImag(double [][] magn,double [][] phase,double [][] result){
        
        int w=magn.length;
        int h=magn[0].length;
        for (int i=0;i<w;i++){
            for (int ii=0;ii<h;ii++){
                result[i][ii]=Math.abs(magn[i][ii])*(Math.sin(phase[i][ii]));
            }
        }
    }

    
    
    class DFT{
        double [][] realIn;double [][] imagIn;double [][] realOut;double [][] imagOut;
        int startW;int endW;
        boolean direct;
        DFT(double [][] realIn,double [][] imagIn,double [][] realOut,double [][] imagOut,int startW,int endW,boolean direct){
            this.realIn=realIn;
            this.imagIn=imagIn;
            this.realOut=realOut;
            this.imagOut=imagOut;
            this.direct=direct;
            this.startW=startW;
            this.endW=endW;
                    
            
        }
        
        public void run(){
        
            int w=realIn.length;
            int h=realIn[0].length;
            double dw=(double)w;
            double dh=(double)h;
            double di=(double)startW;
            if (direct){
                for (int i=startW;i<endW;i++,di++){
                    double dii=0;
                    for (int ii=0;ii<h;ii++,dii++){
                        realOut[i][ii]=0;
                        imagOut[i][ii]=0;
                        double dj=0;
                        for (int j=0;j<w;j++,dj++){
                            double djj=0;
                            for (int jj=0;jj<h;jj++,djj++){
                                realOut[i][ii]+=realIn[j][jj]*Math.cos(-2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                                imagOut[i][ii]+=realIn[j][jj]*Math.sin(-2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                                realOut[i][ii]-=imagIn[j][jj]*Math.sin(-2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                                imagOut[i][ii]+=imagIn[j][jj]*Math.cos(-2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                            }
                        }
                    }
                }
            }
            else{
                for (int i=startW;i<endW;i++,di++){
                    double dii=0;
                    for (int ii=0;ii<h;ii++,dii++){
                        realOut[i][ii]=0;
                        imagOut[i][ii]=0;
                        double dj=0;
                        for (int j=0;j<w;j++,dj++){
                            double djj=0;
                            for (int jj=0;jj<h;jj++,djj++){
                                realOut[i][ii]+=realIn[j][jj]*Math.cos(2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                                imagOut[i][ii]+=realIn[j][jj]*Math.sin(2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                                realOut[i][ii]-=imagIn[j][jj]*Math.sin(2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                                imagOut[i][ii]+=imagIn[j][jj]*Math.cos(2*Math.PI*((di*dj/dw)+(dii*djj/dh)))/ Math.sqrt(dw*dh);
                            }
                        }
                    }
                }
            }
        }
    }
    
    void imshow(double [][] im,String title){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[i][ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    class FastFT{
        double [][] realIn;double [][] imagIn;double [][] realOut;double [][] imagOut;
        int start;int end;
        boolean direct;
        boolean firstPart=true;
        boolean test;
        FastFT(double [][] realIn,double [][] imagIn,double [][] realOut,double [][] imagOut,int start,int end,boolean direct){
            this.realIn=realIn;
            this.imagIn=imagIn;
            this.realOut=realOut;
            this.imagOut=imagOut;
            this.direct=direct;
            this.start=start;
            this.end=end;
            this.test=false;
            
        }
        
        FastFT(double [][] realIn,double [][] imagIn,double [][] realOut,double [][] imagOut,int start,int end,boolean direct,boolean test){
            this.realIn=realIn;
            this.imagIn=imagIn;
            this.realOut=realOut;
            this.imagOut=imagOut;
            this.direct=direct;
            this.start=start;
            this.end=end;
            this.test=test;
            
        }
        
        public void changeSecondStep(){
            firstPart=false;
        }
        
        public void run(){
            if (firstPart){
                run1();
            }
            else{
                run2();
            }
        }
        //first part
        public void run1(){
            
            int w=realIn.length;
            int h=realIn[0].length;
            double[] vr = new double[h]; 
            double[] vi = new double[h]; 
            double[] vfr = new double[h]; 
            double[] vfi = new double[h]; 
            FastFT1D fft1d = new FastFT1D(test);
            if (direct){
                for (int i=start;i<end;i++){
                    for (int ii = 0; ii < h; ii++){
                        vr[ii] = (realIn[i][ii]);
                        //IJ.log("realIn "+realIn[i][ii]);
                        vi[ii] = (imagIn[i][ii]);
                    }
                    fft1d.fft(vr,vi,vfr,vfi); 
                    for (int n = 0; n < h; n++){
                            realOut[i][n] = (vfr[n]);
                            imagOut[i][n] = (vfi[n]);
                    }
                            
                }
                
                
            }
            else{
                for (int i=start;i<end;i++){
                    for (int ii = 0; ii < h; ii++){
                        vr[ii] = (realIn[i][ii]);
                        vi[ii] = (imagIn[i][ii]);
                    }
                    fft1d.ifft(vr,vi,vfr,vfi); 
                    for (int n = 0; n < h; n++){
                            realOut[i][n] = (vfr[n]);
                            imagOut[i][n] = (vfi[n]);
                    }
                            
                }
                
                
            }
        }
        
        
        //second part
        public void run2(){
            
            int w=realIn.length;
            int h=realIn[0].length;
            double[] vr = new double[w]; 
            double[] vi = new double[w]; 
            double[] vfr = new double[w]; 
            double[] vfi = new double[w]; 
            FastFT1D fft1d = new FastFT1D(test);
            if (direct){
                for (int i=start;i<end;i++){
                    for (int ii = 0; ii < w; ii++){
                        vr[ii] = (realIn[ii][i]);
                        vi[ii] = (imagIn[ii][i]);
                    }
                    fft1d.fft(vr,vi,vfr,vfi); 
                    for (int n = 0; n < w; n++){
                            realOut[n][i] = (vfr[n]);
                            imagOut[n][i] = (vfi[n]);
                    }
                }
                
            }
            else{
                for (int i=start;i<end;i++){
                    for (int ii = 0; ii < w; ii++){
                        vr[ii] = (realIn[ii][i]);
                        vi[ii] = (imagIn[ii][i]);
                    }
                    fft1d.ifft(vr,vi,vfr,vfi); 
                    for (int n = 0; n < w; n++){
                            realOut[n][i] = (vfr[n]);
                            imagOut[n][i] = (vfi[n]);
                    }
                }
            }
        }
    }
    
    
    
    
    
    
    
    
    
    class FastFT3D{
        float [][][] realIn;float [][][] imagIn;float [][][] realOut;float [][][] imagOut;
        boolean direct;
        int start;int end;
        boolean firstPart=true;
        boolean secondPart=false;
        boolean test;
        FastFT3D(float [][][] realIn,float [][][] imagIn,float [][][] realOut,float [][][] imagOut,int start,int end,boolean direct){
            this.realIn=realIn;
            this.imagIn=imagIn;
            this.realOut=realOut;
            this.imagOut=imagOut;
            this.direct=direct;
            this.start=start;
            this.end=end;
            this.test=false;
            
        }
        
        FastFT3D(float [][][] realIn,float [][][] imagIn,float [][][] realOut,float [][][] imagOut,int start,int end,boolean direct,boolean test){
            this.realIn=realIn;
            this.imagIn=imagIn;
            this.realOut=realOut;
            this.imagOut=imagOut;
            this.direct=direct;
            this.start=start;
            this.end=end;
            this.test=test;
            
        }
        
        
        public void changeFirstStep(){
            firstPart=true;
            secondPart=false;
        }
        public void changeSecondStep(){
            firstPart=false;
            secondPart=true;
        }
        public void changeThirdStep(){
            firstPart=false;
            secondPart=false;
        }
        
        public void run(){
            if (firstPart){
                run1();
                
            }
            else if (secondPart){
                run2();
            }
            else{
                run3();
            }
        }
        //first part
        public void run1(){
            
            int w=realIn.length;
            int h=realIn[0].length;
            int d=realIn[0][0].length;
            double[] vr = new double[d]; 
            double[] vi = new double[d]; 
            double[] vfr = new double[d]; 
            double[] vfi = new double[d]; 
            
            FastFT1D fft1d = new FastFT1D(test);
            if (direct){
                for (int i=start;i<end;i++){
                    for (int ii = 0; ii < h; ii++){
                        for (int iii = 0; iii < d; iii++){
                            vr[iii] = (double)(realIn[i][ii][iii]);
                            //IJ.log("realIn "+realIn[i][ii]);
                            vi[iii] = (double)(imagIn[i][ii][iii]);
                        }

                        fft1d.fft(vr,vi,vfr,vfi); 
                        for (int n = 0; n < d; n++){
                                realOut[i][ii][n] = (float)(vfr[n]);
                                imagOut[i][ii][n] = (float)(vfi[n]);
                        }
                    }  
                }
                
                
            }
            else{
                for (int i=start;i<end;i++){
                    for (int ii = 0; ii < h; ii++){
                        for (int iii = 0; iii < d; iii++){
                            vr[iii] = (double)(realIn[i][ii][iii]);
                            //IJ.log("realIn "+realIn[i][ii]);
                            vi[iii] = (double)(imagIn[i][ii][iii]);
                        }

                        fft1d.ifft(vr,vi,vfr,vfi); 
                        for (int n = 0; n < d; n++){
                                realOut[i][ii][n] = (float)(vfr[n]);
                                imagOut[i][ii][n] = (float)(vfi[n]);
                        }
                    }   
                }
                
                
            }
        }
        
        
        //second part
        public void run2(){
            
            int w=realIn.length;
            int h=realIn[0].length;
            int d=realIn[0][0].length;
            double[] vr = new double[w]; 
            double[] vi = new double[w]; 
            double[] vfr = new double[w]; 
            double[] vfi = new double[w]; 
            FastFT1D fft1d = new FastFT1D(test);
            if (direct){
                for (int i=start;i<end;i++){
                    for (int iii = 0; iii < d; iii++){
                        for (int ii = 0; ii < w; ii++){
                            vr[ii] = (double)(realIn[ii][i][iii]);
                            vi[ii] = (double)(imagIn[ii][i][iii]);
                        }
                        fft1d.fft(vr,vi,vfr,vfi); 
                        for (int n = 0; n < w; n++){
                                realOut[n][i][iii] = (float)(vfr[n]);
                                imagOut[n][i][iii] = (float)(vfi[n]);
                        }
                    }
                }
                
            }
            else{
                for (int i=start;i<end;i++){
                    for (int iii = 0; iii < d; iii++){
                        for (int ii = 0; ii < w; ii++){
                            vr[ii] = (double)(realIn[ii][i][iii]);
                            vi[ii] = (double)(imagIn[ii][i][iii]);
                        }
                        fft1d.ifft(vr,vi,vfr,vfi); 
                        for (int n = 0; n < w; n++){
                                realOut[n][i][iii] = (float)(vfr[n]);
                                imagOut[n][i][iii] = (float)(vfi[n]);
                        }
                    }
                }
            }
        }
        
        
        //second part
        public void run3(){
            
            int w=realIn.length;
            int h=realIn[0].length;
            int d=realIn[0][0].length;
            double[] vr = new double[h]; 
            double[] vi = new double[h]; 
            double[] vfr = new double[h]; 
            double[] vfi = new double[h]; 
            
            FastFT1D fft1d = new FastFT1D(test);
            if (direct){
                for (int i=start;i<end;i++){
                    for (int iii = 0; iii < d; iii++){
                        for (int ii = 0; ii < h; ii++){
                            vr[ii] = (double)(realIn[i][ii][iii]);
                            vi[ii] = (double)(imagIn[i][ii][iii]);
                        }
                        fft1d.fft(vr,vi,vfr,vfi); 
                        for (int n = 0; n < h; n++){
                                realOut[i][n][iii] = (float)(vfr[n]);
                                imagOut[i][n][iii] = (float)(vfi[n]);
                        }
                    }
                }
                
            }
            else{
                for (int i=start;i<end;i++){
                    for (int iii = 0; iii < d; iii++){
                        for (int ii = 0; ii < h; ii++){
                            vr[ii] = (double)(realIn[i][ii][iii]);
                            vi[ii] = (double)(imagIn[i][ii][iii]);
                        }
                        fft1d.ifft(vr,vi,vfr,vfi); 
                        for (int n = 0; n < h; n++){
                                realOut[i][n][iii] = (float)(vfr[n]);
                                imagOut[i][n][iii] = (float)(vfi[n]);
                        }
                    }
                }
            }
        }
    }
    
    
    
    
    
    public class FastFT1D {
        boolean test;
        double TWOPI=2*Math.PI;
	/** constructor for the use of the 1D transformation routines */ 
	public FastFT1D() {
		test=false;
	}
        
        
        public FastFT1D(boolean test) {
		this.test=test;
	}
        
        /** constructor for the use of the 1D transformation routines */ 
	public void fft(double [] realIn,double [] imagIn,double [] realOut,double [] imagOut) {
            
            boolean allZero=true;
            loop:for (int i=0;i<realIn.length;i++){
                if ((realIn[i]!=0)||(imagIn[i]!=0)){
                    allZero=false;
                    break loop;
                }
            }
            if (!allZero){
                this.transform(realIn,imagIn,realOut,imagOut,1);
            }
            else{
                for (int i=0;i<realIn.length;i++){
                    realOut[i]=0;
                    imagOut[i]=0;
                }
            }
	}
        
        /** constructor for the use of the 1D transformation routines */ 
	public void ifft(double [] realIn,double [] imagIn,double [] realOut,double [] imagOut) {
            
            boolean allZero=true;
            loop:for (int i=0;i<realIn.length;i++){
                if ((realIn[i]!=0)||(imagIn[i]!=0)){
                    allZero=false;
                    break loop;
                }
            }
            if (!allZero){
                this.transform(realIn,imagIn,realOut,imagOut,-1);
            }
            else{
                for (int i=0;i<realIn.length;i++){
                    realOut[i]=0;
                    imagOut[i]=0;
                }
            }
	}
        
        
        /** iterative version of fft */  
	private void transform(double [] realIn,double [] imagIn,double [] realOut,double [] imagOut,int inverseFFT) {
            
                int N = realIn.length*2+1;
                double [] x = new double[N];
                for (int t=0;t<realIn.length;t++){
                    x[2*t+1]=realIn[t];
                    x[2*t+2]=imagIn[t];
                }
                
                
                fft(x,realIn.length,inverseFFT);
                
                
                for (int t=0;t<realIn.length;t++){
                    realOut[t]=x[2*t+1]/Math.sqrt(realIn.length);
                    imagOut[t]=x[2*t+2]/Math.sqrt(realIn.length);
                }
	}
        
        
        // compute the FFT of x[], assuming its length is a power of 2
    public void fft(double [] data, int nn, int isign){
        int n, mmax, m, j, istep, i;
        double wtemp, wr, wpr, wpi, wi, theta;
        double tempr, tempi;
   
        n = nn << 1;
        j = 1;
        for (i = 1; i < n; i += 2) {
        if (j > i) {
            tempr = data[j];     
            data[j] = data[i];     
            data[i] = tempr;
            
            tempr = data[j+1]; 
            data[j+1] = data[i+1]; 
            data[i+1] = tempr;
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = 2*mmax;
        theta = -TWOPI/(isign*mmax);
        
        wtemp = Math.sin(0.5*theta);
        wpi = Math.sin(theta);
        
        wpr = -2.0*wtemp*wtemp;
        
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j =i + mmax;
                tempr = wr*data[j]   - wi*data[j+1];
                tempi = wr*data[j+1] + wi*data[j];
                data[j]   = data[i]   - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }
        mmax = istep;
    }


    
    

    }
    
    
    

    

}
}

