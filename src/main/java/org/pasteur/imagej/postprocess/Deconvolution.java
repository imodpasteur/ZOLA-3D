/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;
//import org.pasteur.imagej.utils.*;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.pasteur.imagej.data.*;
import static org.pasteur.imagej.postprocess.ZRendering.nameScatterPlot;
/**
 *
 * @author benoit
 */
public class Deconvolution {
    
    
    
    int [] indexStack;
    int [] indexFrame;
    double [] X;
    double [] Y;
    double [] Z;
    double [] stdX;
    double [] stdY;
    double [] stdZ;
    double [] varX;
    double [] varY;
    double [] varZ;
    
    double currentX;
    double currentY;
    double currentZ;
    int currentIndex;
    
    double [] newX;
    double [] newY;
    double [] newZ;
    
    int threadNumber=4;
    Gaussian [] gaussian;
    Object [] monitor1;
    Object monitor0;
    int cpt;
    double sumGaussian;
    int computation=0;
    
    StackLocalization sl;
    
    public Deconvolution(StackLocalization sl){
        this.sl=sl;
        int number=0;
        for (int i=0;i<sl.fl.size();i++){
            number+=sl.fl.get(i).loc.size();
        }
        indexStack = new int[number];
        indexFrame = new int[number];
        X = new double[number];
        Y = new double[number];
        Z = new double[number];
        stdX = new double[number];
        stdY = new double[number];
        stdZ = new double[number];
        varX = new double[number];
        varY = new double[number];
        varZ = new double[number];
        newX = new double[number];
        newY = new double[number];
        newZ = new double[number];
        for (int i=0,k=0;i<sl.fl.size();i++){
            for (int ii=0;ii<sl.fl.get(i).loc.size();ii++,k++){
                indexStack[k]=i;
                indexFrame[k]=ii;
                X[k]=sl.fl.get(i).loc.get(ii).X;
                Y[k]=sl.fl.get(i).loc.get(ii).Y;
                Z[k]=sl.fl.get(i).loc.get(ii).Z;
                stdX[k]=sl.fl.get(i).loc.get(ii).crlb_X;
                stdY[k]=sl.fl.get(i).loc.get(ii).crlb_Y;
                stdZ[k]=sl.fl.get(i).loc.get(ii).crlb_Z;
                varX[k]=sl.fl.get(i).loc.get(ii).crlb_X*sl.fl.get(i).loc.get(ii).crlb_X;
                varY[k]=sl.fl.get(i).loc.get(ii).crlb_Y*sl.fl.get(i).loc.get(ii).crlb_Y;
                varZ[k]=sl.fl.get(i).loc.get(ii).crlb_Z*sl.fl.get(i).loc.get(ii).crlb_Z;
                newX[k]=sl.fl.get(i).loc.get(ii).X;
                newY[k]=sl.fl.get(i).loc.get(ii).Y;
                newZ[k]=sl.fl.get(i).loc.get(ii).Z;
            }
        }
        monitor1 = new Object[threadNumber];
        monitor0 = new Object();
        gaussian = new Gaussian[threadNumber];
        for (int i=0;i<threadNumber;i++){
            monitor1[i] = new Object();
            gaussian[i]= new Gaussian(monitor1[i]);
        }
        for (int k=0;k<threadNumber;k++){
            synchronized(monitor1[k]){
                gaussian[k].start();
                try{
                    monitor1[k].wait();//to be sure the threads are launched at this stage and that wait() in the thread is called
                }catch(Exception er){IJ.log("oops: problem wait function for sync monitor 1");}
            }
        }
        
        int numberInThread=(number/threadNumber);
        for (int i=0;i<threadNumber;i++){
            int [] index=null;
            if (numberInThread==0){//if number < threadNumber --> 1 thread only
                if (i==0){
                    index = new int[number];
                    for (int j=0;j<number;j++){
                        index[j]=j;
                    }
                }
            }
            else{
                int end=((i+1)*numberInThread);
                if (i==threadNumber-1){
                    end=number;
                }
                index = new int[end-(i*numberInThread)];
                for (int j=(i*numberInThread),jj=0;j<end;j++,jj++){
                    index[jj]=j;
                }
            }
            gaussian[i].setIndex(index);
        }
        
        //showMap2D(2);
        
        
        
    }
    
    
    void showMap2D(double pixelsizeNM){
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        for (int i=0;i<X.length;i++){
            if (Z[i]>maxZ){
                maxZ=Z[i];
            }
            if (Z[i]<minZ){
                minZ=Z[i];
            }
            if (X[i]>maxX){
                maxX=X[i];
            }
            if (X[i]<minX){
                minX=X[i];
            }
            if (Y[i]>maxY){
                maxY=Y[i];
            }
            if (Y[i]<minY){
                minY=Y[i];
            }
        }
        
        
        
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        
        
        
        
        double x;
        double y;
        double z;
        
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        int depth=(int)Math.ceil((maxZ-minZ)/pixelsizeNM);
        
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        depth=Math.max(depth, 1);
        
        
        ImageStack ims = new ImageStack(width,height);
        
        ImageProcessor [] ip = new FloatProcessor[depth];
        for (int t=0;t<depth;t++){
            ip[t]= new FloatProcessor(width,height);
//            
        }
        for (int i=0;i<width;i++){
            IJ.log("remaining: "+(width-i));
            for (int ii=0;ii<height;ii++){
                for (int iii=0;iii<depth;iii++){
                    
                    currentX=minX+pixelsizeNM*(double)i;
                    currentY=minY+pixelsizeNM*(double)ii;
                    currentZ=minZ+pixelsizeNM*(double)iii;
                    
                    computation=0;
                    cpt=0;
                    sumGaussian=0;
                    synchronized(monitor0){//parrallel gaussian derivative computation
                        for (int k=0;k<threadNumber;k++){
                            synchronized(monitor1[k]){
                                monitor1[k].notify();
                            }
                        }
                        try{
                            monitor0.wait();
                        }catch(Exception ee){IJ.log("error wait function "+ee);}
                    }
                    
                    ip[iii].putPixelValue(i, ii, sumGaussian);
                    
                }
            }
        }
        
        
        for (int t=0;t<depth;t++){
            ims.addSlice("z="+(((double)t)*pixelsizeNM), ip[t]);
        }
        ImagePlus imp = new ImagePlus("MAP "+pixelsizeNM+"nm per px",ims);
        imp.show();
        
        
    }
    
    
    StackLocalization finish(){
        
        ZRendering.colorRendering(null, sl, 10, 0, false);
        for (int k=0;k<threadNumber;k++){
            gaussian[k].kill();
            synchronized(monitor1[k]){
                monitor1[k].notify();
            }
            try{
                gaussian[k].join();
            }catch(Exception eee){IJ.log("Thread impossible to join() "+eee);}
        }
        
        
        for (int j=0;j<X.length;j++){
            sl.fl.get(this.indexStack[j]).loc.get(indexFrame[j]).X=newX[j];
            sl.fl.get(this.indexStack[j]).loc.get(indexFrame[j]).Y=newY[j];
            sl.fl.get(this.indexStack[j]).loc.get(indexFrame[j]).Z=newZ[j];
        }
        //ZRendering.colorRendering(null, sl, 10, 0, false);
        return sl;
    }
    
    
    
    public StackLocalization run(int iter){
        
        int emitterNumber=X.length;
        
        
        
        
        
            
        for (int j=0;j<emitterNumber;j++){
            IJ.log("remaining iteration "+(emitterNumber-j));
            currentX=newX[j];
            currentY=newY[j];
            currentZ=newZ[j];
            currentIndex=j;
            
            for (int i=0;i<iter;i++){
                
                
                computation=1;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                newX[j]+=this.sumGaussian;
                currentX=newX[j];
                
                computation=2;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                newY[j]+=this.sumGaussian;
                currentY=newY[j];
                
                
                computation=3;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                newZ[j]+=this.sumGaussian;
                currentZ=newZ[j];
                
            }
            
            
            
            
            
        }
        
        
        return finish();
        
        
    }
    
    
    
    
    
    
    
    
    
    public StackLocalization runMulti(int iter){
        
        int emitterNumber=X.length;
        double alpha=1;
        
        double step=1;//(nm)
        
        
            
        for (int j=0;j<emitterNumber;j++){
            IJ.log("remaining iteration "+(emitterNumber-j));
            currentX=newX[j];
            currentY=newY[j];
            currentZ=newZ[j];
            currentIndex=-1;
            double Gelse;
            double Gelseprime;
            double Gelsesecond;
            double Gcurrent;
            double Gcurrentprime;
            double Gcurrentsecond;
            double grad;
            for (int i=0;i<10;i++){
                
                
                computation=4;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                Gelse=sumGaussian;
                
                computation=5;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                Gelseprime=sumGaussian;
                
                
                computation=8;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                Gelsesecond=sumGaussian;
                
                Gcurrent=Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) ));
                Gcurrentprime=(-(currentX-X[j])/varX[j])*Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) ));
                Gcurrentsecond=(-1/varX[j])*Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) )) + (-(currentX-X[j])/varX[j])*Gcurrentprime;
                
                
                grad=(Gcurrent*Gelseprime+Gcurrentprime*Gelse)/(Gcurrentprime*Gelseprime+Gcurrentsecond*Gelse+Gelsesecond*Gcurrent+Gelseprime*Gcurrentprime);
                
                while (grad>step){
                    grad/=10;
                }
                while (grad<-step){
                    grad/=10;
                }
                
                
                newX[j]-=alpha*grad;
                currentX=newX[j];
                
                
                
                
                
                
                
                computation=4;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                Gelse=sumGaussian;
                
                computation=6;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                Gelseprime=sumGaussian;
                
                
                computation=9;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                Gelsesecond=sumGaussian;
                
                Gcurrent=Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) ));
                Gcurrentprime=(-(currentY-Y[j])/varY[j])*Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) ));
                Gcurrentsecond=(-1/varY[j])*Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) )) + (-(currentY-Y[j])/varY[j])*Gcurrentprime;
                
                
                grad=(Gcurrent*Gelseprime+Gcurrentprime*Gelse)/(Gcurrentprime*Gelseprime+Gcurrentsecond*Gelse+Gelsesecond*Gcurrent+Gelseprime*Gcurrentprime);
                
                while (grad>step){
                    grad/=10;
                }
                while (grad<-step){
                    grad/=10;
                }
                
                IJ.log("gr "+grad);
                
                newY[j]-=alpha*(grad);
                currentY=newY[j];
                
                
                
                
                
                
                
                
                computation=4;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                Gelse=sumGaussian;
                
                computation=7;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                Gelseprime=sumGaussian;
                
                
                computation=10;
                cpt=0;
                sumGaussian=0;
                synchronized(monitor0){//parrallel gaussian derivative computation
                    for (int k=0;k<threadNumber;k++){
                        synchronized(monitor1[k]){
                            monitor1[k].notify();
                        }
                    }
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                }
                //IJ.log("partial X "+sumGaussian+"  "+newX[j]);
                Gelsesecond=sumGaussian;
                
                Gcurrent=Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) ));
                Gcurrentprime=(-(currentZ-Z[j])/varZ[j])*Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) ));
                Gcurrentsecond=(-1/varZ[j])*Math.exp(-0.5*(  ((currentX-X[j])*(currentX-X[j])/varX[j]) + ((currentY-Y[j])*(currentY-Y[j])/varY[j])   +  ((currentZ-Z[j])*(currentZ-Z[j])/varZ[j]) )) + (-(currentZ-Z[j])/varZ[j])*Gcurrentprime;
                
                grad=(Gcurrent*Gelseprime+Gcurrentprime*Gelse)/(Gcurrentprime*Gelseprime+Gcurrentsecond*Gelse+Gelsesecond*Gcurrent+Gelseprime*Gcurrentprime);
                
                while (grad>step){
                    grad/=10;
                }
                while (grad<-step){
                    grad/=10;
                }
                
                newZ[j]-=alpha*(grad);
                currentZ=newZ[j];
                
                
                
            }
            
            
            
            
            
        }
        
        
        return finish();
        
        
    }
    
    
    
    class Gaussian extends Thread{
        
        int [] index=null;
        boolean killer=false;
        Object monitor;
        
        Gaussian(Object monitor){
            this.monitor=monitor;
        }
        
        public void kill(){
            killer=true;
        }
        
        
        
        public void setIndex(int [] index){
            this.index=index;
        }
        
        
        
        public void run(){
            synchronized(monitor) {
                monitor.notify();
                
                loop:while (true){
                
                    
                    try{ 
                        monitor.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                    
                    if (killer){
                        break loop;
                    }
                    
                    //HERE: COMPUTATION
                    //computation= { 0:likelihood ; 1:partialX ; 2:partialY ; 3:partialZ}
                    double sum=0;
                    if (index!=null){
                        if (computation==0){
                            for (int i=0;i<index.length;i++){
                                sum+=Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                            }
                        }
                        else if (computation==1){
                            for (int i=0;i<index.length;i++){
                                sum+=(-(currentX-X[index[i]])/varX[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                            }
                        }
                        else if (computation==2){
                            for (int i=0;i<index.length;i++){
                                sum+=(-(currentY-Y[index[i]])/varY[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                            }
                        }
                        else if (computation==3){
                            for (int i=0;i<index.length;i++){
                                sum+=(-(currentZ-Z[index[i]])/varZ[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                            }
                        }
                        //WITHOUT CURRENT INDEX
                        else if (computation==4){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                        else if (computation==5){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=(-(currentX-X[index[i]])/varX[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                        else if (computation==6){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=(-(currentY-Y[index[i]])/varY[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                        else if (computation==7){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=(-(currentZ-Z[index[i]])/varZ[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                        
                        else if (computation==8){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=(-(currentX-X[index[i]])/varX[index[i]])*(-(currentX-X[index[i]])/varX[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                    sum+=(-1/varX[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                        else if (computation==9){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=(-(currentY-Y[index[i]])/varY[index[i]])*(-(currentY-Y[index[i]])/varY[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                    sum+=(-1/varY[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                        else if (computation==10){
                            for (int i=0;i<index.length;i++){
                                if (index[i]!=currentIndex){
                                    sum+=(-(currentZ-Z[index[i]])/varZ[index[i]])*(-(currentZ-Z[index[i]])/varZ[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                    sum+=(-1/varZ[index[i]])*Math.exp(-0.5*(  ((currentX-X[index[i]])*(currentX-X[index[i]])/varX[index[i]]) + ((currentY-Y[index[i]])*(currentY-Y[index[i]])/varY[index[i]])   +  ((currentZ-Z[index[i]])*(currentZ-Z[index[i]])/varZ[index[i]]) ));
                                }
                            }
                        }
                    }
                    
                    
                    
                    synchronized(monitor0) {
                        sumGaussian+=sum;
                        cpt++;
                        if (cpt==threadNumber){
                            monitor0.notify();
                        }
                    }
                    
                
                    if (killer){
                        break loop;
                    }
                }


            }
        }
        
        
        
    }
    
    
    
}
