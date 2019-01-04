/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;


import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.utils.Matrixe;
import jcuda.Pointer;
import ij.IJ;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.cublasHandle;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import jcuda.jcusparse.cusparseHandle;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import jcuda.runtime.cudaError;
import jcuda.runtime.cudaMemcpyKind;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice;
import java.util.Random;

/** 
 *
 * @author benoit
 * One LocalizationMany per thread
 * LocalizationMany could potentially deal with multi cam and multi frame
 */
public class LocalizationMultiEmitter_ {
    
    Random random=new Random();
    
    
    int cpt=0;
    Object [][] monitor1;
    Object [][] monitor2;
    Object monitor0;
    Object monitor02;
    int nbProcess;//number of phase (1 per camera)
    //boolean [][] toBeblocked;
    ComputeLikelihood [][] cl;
    ComputeLikelihoodWithoutPSFupdating [][] clwo;
    
    
    double [][][] fisherMatrix;
    
    double xsave;
    
    boolean withregistration=false;
    double positionRefX;
    double positionRefY;
    double decX;
    double decY;
    double dx1;
    double dy1;
    double dz1;
    double dx2;
    double dy2;
    double dz2;
    PolynomialFit pf;
    
    
    double amplitudeMutation=1;//during mutation, position change at maximum 1000nm
    double probaPhotonChange=.1;//during mutation, photon number can change at maximum 10 %
    
    double maxShift=50;//At one step, the maximum gradient shift is 50nm only
    
    double [][] thetaXhusband;
    double [][] thetaXwife;
    double [][] thetaYhusband;
    double [][] thetaYwife;
    double [][] thetaZhusband;
    double [][] thetaZwife;
    double [][][] thetaAhusband;
    double [][][] thetaAwife;
    double [][][] thetaBhusband;
    double [][][] thetaBwife;
    
    double [] meanX;
    double [] meanY;
    double [] meanZ;
    double [][] meanA;
    double [][] meanB;
    
    double [][][] save;
    
    double [] stdX;
    double [] stdY;
    double [] stdZ;
    double [][] stdA;
    double [][] stdB;
    
    double bestLikelihood;
    double [] bestX;
    double [] bestY;
    double [] bestZ;
    double [][] bestA;
    double [][] bestB;
    
    
    double [][] likelihood;//the likelihood for each camera and each individual
    double [][] sortedLikelihood;//the likelihood sorted in increasing order according to each individual ; second dim=2: likelihood&index&cumul
    double [][] thetaX;
    double [][] thetaY;
    double [][] thetaZ;
    double [] thetaZfocus;
    double [][][] thetaA;
    double [][][] thetaB;
    
    double [] x2, y2, z2;
    double [] v = new double [3];
    
    
    
    double h=.0005;//.5 nm ; do not reduce it because we use float instead of double for PSF generator
    
    double epsilon=.0001; // distance in micrometer (stop criterion)
    
    DataPhase_ [] dparam;
    int width;//width of image
    int height;//height of image
    double [] subwindow;// 4 dimensions in the vect: []nbProcess ; size;width;height
    int totalsize;
    int iterMax;
    //for one cam
    double photonThreshold=0;
    
    double minZ;
    double maxZ;
    
    int [][] id;//vector for each camera...
    
    
    double [][] modelPSF;
    
    double [] resCRLB;
            
    //in case of multiframe fitting, the likelihood is just the sum of likelihoods according to each image...
    
    int numberPSFperModel;
    int popsize;
    int popsizeSelection;//corresponds to number of peaple used for selection (popsizeSelection<=popsize)
    
    public LocalizationMultiEmitter_(DataPhase_ dparam,int iterMax,double minZ, double maxZ,int popsize,int popsizeSelection,int numberPSFperModel){
        
        this.dparam = new DataPhase_[1];
        this.dparam[0] = dparam;
        this.localizationMultiEmitter_(iterMax, minZ, maxZ, popsize,popsizeSelection, numberPSFperModel);
            
        
    }
    
    
    
    // for multiple cam
    public LocalizationMultiEmitter_(DataPhase_ [] dparam, int iterMax,double minZ, double maxZ,int popsize,int popsizeSelection,int numberPSFperModel){
        
        this.dparam = dparam;
        this.localizationMultiEmitter_(iterMax, minZ, maxZ, popsize,popsizeSelection, numberPSFperModel);
    }




    private void localizationMultiEmitter_(int iterMax,double minZ, double maxZ,int popsize,int popsizeSelection,int numberPSFperModel){
        
        
        
        this.popsize=popsize;
        this.numberPSFperModel=numberPSFperModel;
        this.popsizeSelection=popsizeSelection;
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.iterMax=iterMax;
        nbProcess=dparam.length;
        this.width=dparam[0].param.sizeoutput;//width & height are supposed to be the same for each source
        this.height=dparam[0].param.sizeoutput;
        maxShift=dparam[0].param.xystep/2.;
        
        id = new int [nbProcess][popsize];
        resCRLB=new double [3];
        for (int i=0;i<nbProcess;i++){
            for (int pop=0;pop<popsize;pop++){
                id[i][pop]=0;
            }
        }
        fisherMatrix=new double [nbProcess][5][5];
        modelPSF=new double [nbProcess][width*height];
        
        
        thetaAhusband=new double[this.nbProcess][popsize/2][numberPSFperModel];
        thetaBhusband=new double[this.nbProcess][popsize/2][numberPSFperModel];
        thetaXhusband=new double[popsize/2][numberPSFperModel];
        thetaYhusband=new double[popsize/2][numberPSFperModel];
        thetaZhusband=new double[popsize/2][numberPSFperModel];
        thetaAwife=new double[this.nbProcess][popsize/2][numberPSFperModel];
        thetaBwife=new double[this.nbProcess][popsize/2][numberPSFperModel];
        thetaXwife=new double[popsize/2][numberPSFperModel];
        thetaYwife=new double[popsize/2][numberPSFperModel];
        thetaZwife=new double[popsize/2][numberPSFperModel];
        
        
        bestA=new double[this.nbProcess][numberPSFperModel];
        bestB=new double[this.nbProcess][numberPSFperModel];
        bestX=new double[numberPSFperModel];
        bestY=new double[numberPSFperModel];
        bestZ=new double[numberPSFperModel];
        bestLikelihood=Double.POSITIVE_INFINITY;
        
        save=new double[this.nbProcess][popsize][numberPSFperModel];
        
        meanA=new double[this.nbProcess][numberPSFperModel];
        meanB=new double[this.nbProcess][numberPSFperModel];
        meanX=new double[numberPSFperModel];
        meanY=new double[numberPSFperModel];
        meanZ=new double[numberPSFperModel];
        stdA=new double[this.nbProcess][numberPSFperModel];
        stdB=new double[this.nbProcess][numberPSFperModel];
        stdX=new double[numberPSFperModel];
        stdY=new double[numberPSFperModel];
        stdZ=new double[numberPSFperModel];
        
        thetaA=new double[this.nbProcess][popsize][numberPSFperModel];
        thetaB=new double[this.nbProcess][popsize][numberPSFperModel];
        
        thetaX=new double[popsize][numberPSFperModel];
        thetaY=new double[popsize][numberPSFperModel];
        thetaZ=new double[popsize][numberPSFperModel];
        
        thetaZfocus=new double[nbProcess];
        for (int i=0;i<nbProcess;i++){
            thetaZfocus[i]=dparam[i].param.Zfocus;
        }
        
        
            
        monitor1 = new Object[nbProcess][popsize];
        monitor2 = new Object[nbProcess][popsize];
        monitor0 = new Object();
        monitor02 = new Object();

        //toBeblocked=new boolean [nbProcess][popsize];
        likelihood=new double [nbProcess][popsize];
        sortedLikelihood=new double [popsize][3];
        cl = new ComputeLikelihood[nbProcess][popsize];
        clwo = new ComputeLikelihoodWithoutPSFupdating[nbProcess][popsize];
        for (int i=0;i<nbProcess;i++){
            for (int pop=0;pop<popsize;pop++){
                monitor1[i][pop] = new Object();
                monitor2[i][pop] = new Object();
                //toBeblocked[i][pop]=true;
                cl[i][pop]=new ComputeLikelihood(monitor1[i][pop],i,pop);
                clwo[i][pop]=new ComputeLikelihoodWithoutPSFupdating(monitor2[i][pop],i,pop);
            }
        }



        for (int k=0;k<nbProcess;k++){
            for (int pop=0;pop<popsize;pop++){
                synchronized(monitor1[k][pop]){
                    cl[k][pop].start();
                    try{
                        monitor1[k][pop].wait();//to be sure the threads are launched at this stage and that wait() in the thread is called
                    }catch(Exception er){IJ.log("oops: problem wait function for sync monitor 1");}
                }
            }
        }
        
        for (int k=0;k<nbProcess;k++){
            for (int pop=0;pop<popsize;pop++){
                synchronized(monitor2[k][pop]){
                    clwo[k][pop].start();
                    try{
                        monitor2[k][pop].wait();//to be sure the threads are launched at this stage and that wait() in the thread is called
                    }catch(Exception er){IJ.log("oops: problem wait function for sync monitor 1");}
                }
            }
        }
            
        
        
    }
    
    
    
    
    
    
    
    
    
    public void setMinimumPhotonNumber(double value){
        this.photonThreshold=value;
        
        
        
    }
    
    
    public double [] getA(int idParam){
        return thetaA[idParam][(int)sortedLikelihood[0][1]];
    }
    public double [] getB(int idParam){
        return thetaB[idParam][(int)sortedLikelihood[0][1]];
    }
    public double [] getX(){
        return thetaX[(int)sortedLikelihood[0][1]];
    }
    public double [] getY(){
        return thetaY[(int)sortedLikelihood[0][1]];
    }
    public double [] getZ(){
        return thetaZ[(int)sortedLikelihood[0][1]];   }
    
    
    
    public double getLikelihood(){
        return (sortedLikelihood[0][0]);
    }
    
    
    
    
    
    public void init(double a, double b, double x, double y, double z){
        for (int pop=0;pop<popsize;pop++){
            for (int j=0;j<this.numberPSFperModel;j++){
                this.thetaX[pop][j]=x;
                this.thetaY[pop][j]=y;
                this.thetaZ[pop][j]=z;
                for (int i=0;i<this.nbProcess;i++){
                    this.thetaA[i][pop][j]=a;
                    this.thetaB[i][pop][j]=b;
                }

            }
        }
        
        
        
    }
    
    public void init(double [] a, double [] b, double x, double y, double z){
        for (int pop=0;pop<popsize;pop++){
            for (int j=0;j<this.numberPSFperModel;j++){
                this.thetaX[pop][j]=x;
                this.thetaY[pop][j]=y;
                this.thetaZ[pop][j]=z;
                for (int i=0;i<this.nbProcess;i++){
                    this.thetaA[i][pop][j]=a[i];
                    this.thetaB[i][pop][j]=b[i];
                }

            }
        }
        
    }
    
    
    
    //input corresponds to a list of possible parameters
    public void init(double [] aList, double [] bList, double [] xList, double [] yList, double [] zList,double randomPercentageAB, double sigmaRandomXYZ){
        
        
        
        int nb=aList.length;
        
        if (this.nbProcess==1){
            
            for (int pop=0;pop<popsize;pop++){
                int posit = (int)Math.random()*nb;
                for (int j=0;j<this.numberPSFperModel;j++){
                    
                    this.thetaX[pop][j]=xList[posit];
                    this.thetaY[pop][j]=yList[posit];
                    this.thetaZ[pop][j]=zList[posit];
                    this.thetaA[0][pop][j]=aList[posit];
                    if (j==0){
                        this.thetaB[0][pop][j]=bList[posit];
                    }
                    else{
                        this.thetaB[0][pop][j]=0;
                    }
                }
            }
        }
        else{
            IJ.log("initialization function is not developped for dealing with multi camera yet");
        }
        
    }
    
    
    
    
    
    //input corresponds to a list of possible parameters
    public void initSorted(double [] aList, double [] bList, double [] xList, double [] yList, double [] zList,double randomPercentageAB, double sigmaRandomXYZ){
        
        
        int [] order=getSortedOrder(xList,yList,zList);
        //IJ.log("ord "+order.length+"  "+xList.length);
        //first order the initial list to create a path
        int nb=aList.length;
        int shortNB=Math.min(numberPSFperModel, nb);
        ArrayList<Integer> list = new ArrayList<Integer>();//used to generate random list
        for (int i=0;i<nb;i++){
            list.add(i);
        }
        int [] positions = new int[shortNB];
        if (this.nbProcess==1){
            
            for (int pop=0;pop<popsize;pop++){
                Collections.shuffle(list);//generate random serie
                for (int i=0;i<shortNB;i++){
                    positions[i]=list.get(i);
                }
                Arrays.sort(positions);//sort it to have roughtly the same ordering for each individual -> better after to compute mean according to population
                double k=0;
                for (int j=0;j<this.numberPSFperModel;j++){
                    int posit=order[positions[(int)k]];
                    
                    this.thetaX[pop][j]=xList[posit]+random.nextGaussian()*sigmaRandomXYZ;
                    this.thetaY[pop][j]=yList[posit]+random.nextGaussian()*sigmaRandomXYZ;
                    this.thetaZ[pop][j]=zList[posit]+random.nextGaussian()*sigmaRandomXYZ;
                    this.thetaA[0][pop][j]=aList[posit]+random.nextGaussian()*aList[posit]*randomPercentageAB;
                    if (j==0){
                        this.thetaB[0][pop][j]=bList[posit]+random.nextGaussian()*bList[posit]*randomPercentageAB;
                    }
                    else{
                        this.thetaB[0][pop][j]=0;
                    }
                    k+=((double)shortNB/(double)numberPSFperModel);
                }
            }
        }
        else{
            IJ.log("initialization function is not developped for dealing with multi camera yet");
        }
        
    }
    
    
    
    
    
    
    
    //input corresponds to a list of possible parameters
    public void initSorted(double [] aList, double [] bList, double [] xList, double [] yList, double [] zList){
        
        
        int [] order=getSortedOrder(xList,yList,zList);
        //IJ.log("ord "+order.length+"  "+xList.length);
        //first order the initial list to create a path
        int nb=aList.length;
        int shortNB=Math.min(numberPSFperModel, nb);
        ArrayList<Integer> list = new ArrayList<Integer>();//used to generate random list
        for (int i=0;i<nb;i++){
            list.add(i);
        }
        int [] positions = new int[shortNB];
        if (this.nbProcess==1){
            
            for (int pop=0;pop<popsize;pop++){
                Collections.shuffle(list);//generate random serie
                for (int i=0;i<shortNB;i++){
                    positions[i]=list.get(i);
                }
                Arrays.sort(positions);//sort it to have roughtly the same ordering for each individual -> better after to compute mean according to population
                double k=0;
                for (int j=0;j<this.numberPSFperModel;j++){
                    int posit=order[positions[(int)k]];
                    
                    this.thetaX[pop][j]=xList[posit];
                    this.thetaY[pop][j]=yList[posit];
                    this.thetaZ[pop][j]=zList[posit];
                    this.thetaA[0][pop][j]=aList[posit];
                    if (j==0){
                        this.thetaB[0][pop][j]=bList[posit];
                    }
                    else{
                        this.thetaB[0][pop][j]=0;
                    }
                    k+=((double)shortNB/(double)numberPSFperModel);
                }
            }
        }
        else{
            IJ.log("initialization function is not developped for dealing with multi camera yet");
        }
        
    }
    
    
    
    //this function return a list of position corresponding to the path of sorted 3D positions
    private int [] getSortedOrder(double [] x, double [] y, double [] z){
        int nb=x.length;
        if (nb>1){
            double [][] dist = new double [(nb*(nb+1))/2-nb][3];
            
            for (int i=0,k=0;i<nb-1;i++){
                for (int ii=i+1;ii<nb;ii++,k++){
                    dist[k][0] = Math.sqrt((x[i]-x[ii])*(x[i]-x[ii])+(y[i]-y[ii])*(y[i]-y[ii])+(z[i]-z[ii])*(z[i]-z[ii]));
                    dist[k][1] = i;
                    dist[k][2] = ii;
                }
            }
            Arrays.sort(dist, new Comparator<double[]>() {
                @Override
                public int compare(double[] o1, double[] o2) {
                    return ((Double) o1[0]).compareTo(o2[0]);
                }
            });
            
            
            double [] connexionNumber = new double [nb];
            //ArrayList order = new ArrayList();
            int [][] segment = new int [nb-1][2];
            for (int i=0,k=0;k<nb-1;i++){
                int p1=(int)dist[i][1];
                int p2=(int)dist[i][2];
                if ((connexionNumber[p1]<2)&&(connexionNumber[p2]<2)){
                    connexionNumber[p1]++;
                    connexionNumber[p2]++;
                    segment[k][0]=p1;
                    segment[k][1]=p2;
                    k++;
                }
            }
            int start =0;//search start of path
            for (int i=0;i<connexionNumber.length;i++){
                if (connexionNumber[i]<2){
                    start=i;
                    break;
                }
            }
            int [] ordered = new int[nb];
            boolean [] treated= new boolean[nb-1];
            ordered[0]=start;
            for (int i=1;i<nb;i++){
                search:for (int j=0;j<nb-1;j++){
                    if (!treated[j]){
                        if (ordered[i-1]==segment[j][0]){
                            ordered[i]=segment[j][1];
                            treated[j]=true;
                            break search;
                        }
                        if (ordered[i-1]==segment[j][1]){
                            ordered[i]=segment[j][0];
                            treated[j]=true;
                            break search;
                        }
                        
                    }
                }
            }
            return ordered;
                    
        }    
        else{
            return new int[]{0};
        }
        
    }
    
    
    
    
    
    //positionRef: position of center of subwindow
    //dec : decalage entre les 2 patchs
    //dx, dy,dz : drifts
    //pf : PolynomialFit function for registration
    public void setRegistrationParameters(double positionRefX,double positionRefY,double decX,double decY,double dx1,double dy1,double dz1,double dx2,double dy2,double dz2,PolynomialFit pf,double xsave){
        withregistration=true;
        this.decX=decX;
        this.decY=decY;
        this.positionRefX=positionRefX;
        this.positionRefY=positionRefY;
        this.dx1=dx1;
        this.dy1=dy1;
        this.dz1=dz1;
        this.dx2=dx2;
        this.dy2=dy2;
        this.dz2=dz2;
        this.pf=pf;
        this.xsave=xsave;
    }
    
    
    
    public void setSubWindow(double [] subwindow){
        //could be improved
        for (int pop=0;pop<popsize;pop++){
            id[0][pop]=dparam[0].modelMany.setSubWindow(subwindow);
        }
        
    }
    
    
    
    
    
    public void setSubWindow(double [][] subwindow){
        
        for (int u=0;u<subwindow.length;u++){
            //IJ.log("subWin "+id.length+"  "+subwindow.length+"  "+dparam.length+"  "+dparam[u].modelMany);
            for (int pop=0;pop<popsize;pop++){
                id[u][pop]=dparam[u].modelMany.setSubWindow(subwindow[u]);
            }
        }
        
        
    }
    
    public void setSubWindowScmos(double [] subwindow,double [] subwindowSCMOS){
        for (int pop=0;pop<popsize;pop++){
            id[0][pop]=dparam[0].modelMany.setSubWindowScmos(subwindow,subwindowSCMOS);
        }
        
        
        
    }
    
    
    public void setSubWindowScmos(double [][] subwindow,double [][] subwindowSCMOS){
        
        for (int u=0;u<subwindow.length;u++){
            for (int pop=0;pop<popsize;pop++){
                id[u][pop]=dparam[u].modelMany.setSubWindowScmos(subwindow[u],subwindowSCMOS[u]);
            }
        }
        
        
    }
    
    
    
    
    
    
    public void finish(){
        
        
        this.bestLikelihood=Double.POSITIVE_INFINITY;
            
        for (int p=0;p<nbProcess;p++){
            dparam[p].modelMany.getModel(id[p][(int)this.sortedLikelihood[0][1]],modelPSF[p]);
        }
        
        for (int p=0;p<nbProcess;p++){
            for (int pop=0;pop<popsize;pop++){
                dparam[p].modelMany.freePosit(id[p][pop]);
            }
            
        }
        
        
        
        
        
    }
    
    
    
    
    public void kill(){
        for (int k=0;k<nbProcess;k++){
            for (int pop=0;pop<popsize;pop++){
                cl[k][pop].kill();
                synchronized(monitor1[k][pop]){
                    monitor1[k][pop].notify();
                }
                try{
                    cl[k][pop].join();
                }catch(Exception eee){IJ.log("Thread impossible to join() "+eee);}
                
                clwo[k][pop].kill();
                synchronized(monitor2[k][pop]){
                    monitor2[k][pop].notify();
                }
                try{
                    clwo[k][pop].join();
                }catch(Exception eee){IJ.log("Thread impossible to join() "+eee);}
            
            }
        }
    }
    
    
    
    
    
    
    
    
    public double getCRLBX(){
        
        return resCRLB[0];
            
    }
    
    
    
    
    public double getCRLBY(){
        
        return resCRLB[1];
            
    }
    
    
    
    
    public double getCRLBZ(){
        
        return resCRLB[2];
            
    }
    
    
    public double [] getPSF(int idCamera){
        
        return modelPSF[idCamera];
            
    }
    
    void testPhotonNumber(){
        for (int i=0;i<this.nbProcess;i++){
            for (int pop=0;pop<popsize;pop++){
                for (int j=0;j<this.numberPSFperModel;j++){
                    this.thetaA[i][pop][j]=Math.max(1, this.thetaA[i][pop][j]);
                    this.thetaB[i][pop][j]=Math.max(.01, thetaB[i][pop][j]);
                }
            }
        }
    }
    
    
    
    
    public boolean localizeCrossEntropy(){
        
        
        IJ.log("loc...");
                
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeB();
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeB();
        testPhotonNumber();
        
        this.computeLikelihood();
        this.sortLikelihood();
        
        for (int i=0;i<this.popsize;i++){
            //IJ.log("value pop("+i+")  =  "+thetaA[0][(int)sortedLikelihood[i][1]][0]+"  "+thetaB[0][(int)sortedLikelihood[i][1]][0]+"  "+thetaX[(int)sortedLikelihood[i][1]][0]+"  "+thetaY[(int)sortedLikelihood[i][1]][0]+"  "+thetaZ[(int)sortedLikelihood[i][1]][0]+"  "+likelihood[0][(int)sortedLikelihood[i][1]]);
        }
        
        //for (int t=0;t<iterMax;t++){
        for (int t=0;t<10;t++){
            this.computeMeanStdOverPopulation();
            for (int i=0;i<this.numberPSFperModel;i++){
                //IJ.log("mean ("+i+")  =  "+meanA[0][i]+"  "+meanB[0][i]+"  "+meanX[i]+"  "+meanY[i]+"  "+meanZ[i]+"  ");
                //IJ.log("std ("+i+")  =  "+stdA[0][i]+"  "+stdB[0][i]+"  "+stdX[i]+"  "+stdY[i]+"  "+stdZ[i]+"  ");
            }
            this.makeNewPopulationSelection();
            
            this.computeLikelihood();
            
            minimizeA();
            testPhotonNumber();
            
            minimizeB();
            testPhotonNumber();
            
            this.computeLikelihood();
            this.sortLikelihood();
            IJ.log("best lik ("+t+")=   "+this.bestLikelihood);
        }
        
        copyWinnerAtstartingPosition();
        this.computeLikelihood();
        this.sortLikelihood();
        
        
        
        
        
        
        this.finish();
        
        
        return true;
    }
    
    
    public boolean localizeGeneticAlgorithm(){
        
        
        
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeB();
        testPhotonNumber();
        minimizeA();
        testPhotonNumber();
        minimizeB();
        testPhotonNumber();
        
        this.computeLikelihood();
        this.sortLikelihood();
        
        for (int i=0;i<this.popsize;i++){
            //IJ.log("value pop("+i+")  =  "+thetaA[0][(int)sortedLikelihood[i][1]][0]+"  "+thetaB[0][(int)sortedLikelihood[i][1]][0]+"  "+thetaX[(int)sortedLikelihood[i][1]][0]+"  "+thetaY[(int)sortedLikelihood[i][1]][0]+"  "+thetaZ[(int)sortedLikelihood[i][1]][0]+"  "+likelihood[0][(int)sortedLikelihood[i][1]]);
        }
        
        //for (int t=0;t<iterMax;t++){
        for (int t=0;t<20;t++){
            
            this.makeCouples(50);
            
            this.crossOver();
            
            this.mutation(10);
            
            minimizeA();
            testPhotonNumber();
            minimizeB();
            testPhotonNumber();
            
            this.computeLikelihood();
            this.sortLikelihood();
            IJ.log("best lik ("+t+")=   "+this.bestLikelihood);
        }
        
        copyWinnerAtstartingPosition();
        this.computeLikelihood();
        this.sortLikelihood();
        
        
        
        
        
        
        this.finish();
        
        
        return true;
    }
    
    
    
    public void sortLikelihood(){
        
        
        //then copy likelihood to sortedLikelihood (add for each camera)
        for (int pop=0;pop<popsize;pop++){
            sortedLikelihood[pop][0]=0;
            sortedLikelihood[pop][1]=pop;
            for (int p=0;p<nbProcess;p++){
                sortedLikelihood[pop][0]+=likelihood[p][pop];
            }
        }
        //finally: sort likelihood
        Arrays.sort(sortedLikelihood, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o1[0]).compareTo(o2[0]);
            }
        });
        
        if (bestLikelihood>sortedLikelihood[0][0]){
            bestLikelihood=sortedLikelihood[0][0];
            for (int i=0;i<this.numberPSFperModel;i++){
                for (int p=0;p<this.nbProcess;p++){
                    bestA[p][i]=this.thetaA[p][(int)sortedLikelihood[0][1]][i];
                    bestB[p][i]=this.thetaB[p][(int)sortedLikelihood[0][1]][i];
                }
                bestX[i]=this.thetaX[(int)sortedLikelihood[0][1]][i];
                bestY[i]=this.thetaY[(int)sortedLikelihood[0][1]][i];
                bestZ[i]=this.thetaZ[(int)sortedLikelihood[0][1]][i];
            }
        }
        
        
        //for (int pop=0;pop<popsize;pop++){
        //    IJ.log("lik "+sortedLikelihood[pop][0]+"  "+sortedLikelihood[pop][1]+"  "+likelihood[0][(int)sortedLikelihood[pop][1]]);
        //}
    }
    
    
    
    
    
    
    public void computeLikelihood(){
        
        
        //first compute likelihood:
        for (int pop=0;pop<popsize;pop++){
            for (int i=0;i<this.nbProcess;i++){
                
                
                dparam[0].modelMany.setParameters(id[i][pop], thetaX[pop], thetaY[pop], thetaZfocus[i], thetaZ[pop], thetaA[i][pop], thetaB[i][pop]);
                
            }
        }
            

        cpt=0;
        synchronized(monitor0){
            for (int k=0;k<nbProcess;k++){
                for (int pop=0;pop<popsize;pop++){
                
                    synchronized(monitor1[k][pop]){
                        //toBeblocked[k][pop]=true;
                        monitor1[k][pop].notify();
                    }
                }
                
            }
            try{
                monitor0.wait();
            }catch(Exception ee){IJ.log("error wait function "+ee);}
        }
        
    }
    
    
    
    
    
    
    
    public void computeLikelihoodWithoutPSFupdating(){
        
        
        //first compute likelihood:
        for (int pop=0;pop<popsize;pop++){
            for (int i=0;i<this.nbProcess;i++){
                
                
                dparam[0].modelMany.setParameters(id[i][pop], thetaX[pop], thetaY[pop], thetaZfocus[i], thetaZ[pop], thetaA[i][pop], thetaB[i][pop]);
                
            }
        }
            

        cpt=0;
        synchronized(monitor02){
            for (int k=0;k<nbProcess;k++){
                for (int pop=0;pop<popsize;pop++){
                
                    synchronized(monitor2[k][pop]){
                        //toBeblocked[k][pop]=true;
                        monitor2[k][pop].notify();
                    }
                }
                
            }
            try{
                monitor02.wait();
            }catch(Exception ee){IJ.log("error wait function "+ee);}
        }
        
    }
    
    
    void copyWinnerAtstartingPosition(){
        if (bestLikelihood<sortedLikelihood[0][0]){
            sortedLikelihood[0][0]=bestLikelihood;
            for (int i=0;i<this.numberPSFperModel;i++){
                for (int p=0;p<this.nbProcess;p++){
                    this.thetaA[p][(int)sortedLikelihood[0][1]][i]=bestA[p][i];
                    this.thetaB[p][(int)sortedLikelihood[0][1]][i]=bestB[p][i];
                }
                this.thetaX[(int)sortedLikelihood[0][1]][i]=bestX[i];
                this.thetaY[(int)sortedLikelihood[0][1]][i]=bestY[i];
                this.thetaZ[(int)sortedLikelihood[0][1]][i]=bestZ[i];
            }
        }
    }
    
    
    
    public void computeMeanStdOverPopulation(){
        for (int i=0;i<this.numberPSFperModel;i++){
            
            //init mean & std to 0
            for (int p=0;p<this.nbProcess;p++){
                this.meanA[p][i]=0;
                this.meanB[p][i]=0;
                this.stdA[p][i]=0;
                this.stdB[p][i]=0;
            }
            this.meanX[i]=0;
            this.stdX[i]=0;
            this.meanY[i]=0;
            this.stdY[i]=0;
            this.meanZ[i]=0;
            this.stdZ[i]=0;
            
            //compute mean & std
            for (int pop=0;pop<popsizeSelection;pop++){
                for (int p=0;p<this.nbProcess;p++){
                    meanA[p][i]+=this.thetaA[p][(int)sortedLikelihood[pop][1]][i];
                    meanB[p][i]+=this.thetaB[p][(int)sortedLikelihood[pop][1]][i];
                }
                meanX[i]+=this.thetaX[(int)sortedLikelihood[pop][1]][i];
                meanY[i]+=this.thetaY[(int)sortedLikelihood[pop][1]][i];
                meanZ[i]+=this.thetaZ[(int)sortedLikelihood[pop][1]][i];
                
            }
            for (int p=0;p<this.nbProcess;p++){
                meanA[p][i]/=popsizeSelection;
                meanB[p][i]/=popsizeSelection;
            }
            meanX[i]/=popsizeSelection;
            meanY[i]/=popsizeSelection;
            meanZ[i]/=popsizeSelection;
            
            
            
            
            for (int pop=0;pop<popsizeSelection;pop++){
                for (int p=0;p<this.nbProcess;p++){
                    stdA[p][i]+=(this.thetaA[p][(int)sortedLikelihood[pop][1]][i]-meanA[p][i])*(this.thetaA[p][(int)sortedLikelihood[pop][1]][i]-meanA[p][i]);
                    stdB[p][i]+=(this.thetaB[p][(int)sortedLikelihood[pop][1]][i]-meanB[p][i])*(this.thetaB[p][(int)sortedLikelihood[pop][1]][i]-meanB[p][i]);
                }
                stdX[i]+=(this.thetaX[(int)sortedLikelihood[pop][1]][i]-meanX[i])*(this.thetaX[(int)sortedLikelihood[pop][1]][i]-meanX[i]);
                stdY[i]+=(this.thetaY[(int)sortedLikelihood[pop][1]][i]-meanY[i])*(this.thetaY[(int)sortedLikelihood[pop][1]][i]-meanY[i]);
                stdZ[i]+=(this.thetaZ[(int)sortedLikelihood[pop][1]][i]-meanZ[i])*(this.thetaZ[(int)sortedLikelihood[pop][1]][i]-meanZ[i]);
            }
            for (int p=0;p<this.nbProcess;p++){
                stdA[p][i]/=popsizeSelection;
                stdB[p][i]/=popsizeSelection;
                stdA[p][i]=Math.sqrt(stdA[p][i]);
                stdB[p][i]=Math.sqrt(stdB[p][i]);
                //stdA[p][i]=.0;
                //stdB[p][i]=.0;
                
                        
            }
            stdX[i]/=popsizeSelection;
            stdY[i]/=popsizeSelection;
            stdZ[i]/=popsizeSelection;
            stdX[i]=Math.sqrt(stdX[i]);
            stdY[i]=Math.sqrt(stdY[i]);
            stdZ[i]=Math.sqrt(stdZ[i]);
            
            
            
            //IJ.log("WARNING STD fixed");
            
            
            
        }
        
    }
    
    
    
    
    
    
    
    //cross couples
    public void crossOver(){
        
        crossOver(70);
        
    }
    
    public void crossOver(double percentage){
        
        for (int pop=0;pop<popsize/2;pop++){
            
            if (Math.random()*100<percentage){//crossover
                
                int randomPosition=1+(int)(Math.random()*(numberPSFperModel-1));
                
                for (int i=0;i<randomPosition;i++){
                    for (int p=0;p<this.nbProcess;p++){
                        thetaA[p][pop*2][i]=this.thetaAhusband[p][pop][i];
                        thetaB[p][pop*2][i]=this.thetaBhusband[p][pop][i];
                        thetaA[p][pop*2+1][i]=this.thetaAwife[p][pop][i];
                        thetaB[p][pop*2+1][i]=this.thetaBwife[p][pop][i];
                    }
                    thetaX[pop*2][i]=this.thetaXhusband[pop][i];
                    thetaY[pop*2][i]=this.thetaYhusband[pop][i];
                    thetaZ[pop*2][i]=this.thetaZhusband[pop][i];
                    thetaX[pop*2+1][i]=this.thetaXwife[pop][i];
                    thetaY[pop*2+1][i]=this.thetaYwife[pop][i];
                    thetaZ[pop*2+1][i]=this.thetaZwife[pop][i];
                }
                for (int i=0;i<this.numberPSFperModel;i++){
                    for (int p=0;p<this.nbProcess;p++){
                        thetaA[p][pop*2+1][i]=this.thetaAhusband[p][pop][i];
                        thetaB[p][pop*2+1][i]=this.thetaBhusband[p][pop][i];
                        thetaA[p][pop*2][i]=this.thetaAwife[p][pop][i];
                        thetaB[p][pop*2][i]=this.thetaBwife[p][pop][i];
                    }
                    thetaX[pop*2+1][i]=this.thetaXhusband[pop][i];
                    thetaY[pop*2+1][i]=this.thetaYhusband[pop][i];
                    thetaZ[pop*2+1][i]=this.thetaZhusband[pop][i];
                    thetaX[pop*2][i]=this.thetaXwife[pop][i];
                    thetaY[pop*2][i]=this.thetaYwife[pop][i];
                    thetaZ[pop*2][i]=this.thetaZwife[pop][i];
                }
                
                
            }
            else{//no crossover
                for (int i=0;i<this.numberPSFperModel;i++){
                    for (int p=0;p<this.nbProcess;p++){
                        thetaA[p][pop*2][i]=this.thetaAhusband[p][pop][i];
                        thetaB[p][pop*2][i]=this.thetaBhusband[p][pop][i];
                        thetaA[p][pop*2+1][i]=this.thetaAwife[p][pop][i];
                        thetaB[p][pop*2+1][i]=this.thetaBwife[p][pop][i];
                    }
                    thetaX[pop*2][i]=this.thetaXhusband[pop][i];
                    thetaY[pop*2][i]=this.thetaYhusband[pop][i];
                    thetaZ[pop*2][i]=this.thetaZhusband[pop][i];
                    thetaX[pop*2+1][i]=this.thetaXwife[pop][i];
                    thetaY[pop*2+1][i]=this.thetaYwife[pop][i];
                    thetaZ[pop*2+1][i]=this.thetaZwife[pop][i];
                }
            }
            
            
        }
        
        
        
    }
    
    
    //gradient descent for photons numbers
    void minimizeA(){
        
        double h=.01;
        
        double [][] lik=new double [nbProcess][popsize];
        double [][] lik1=new double [nbProcess][popsize];
        double [][] lik2=new double [nbProcess][popsize];
        double [][] lik3=new double [nbProcess][popsize];
        double [][] grad=new double [nbProcess][popsize];
        for (int p=0;p<this.nbProcess;p++){
            for (int pop=0;pop<this.popsize;pop++){
                lik[p][pop]=0;
                lik1[p][pop]=0;
                lik2[p][pop]=0;
                lik3[p][pop]=0;
            }
        }
        
        
        for (int i=0;i<this.numberPSFperModel;i++){
            
            //lik2
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    save[p][pop][i]=thetaA[p][pop][i];
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik2[p][pop]=this.likelihood[p][pop];
                }
            }
            
            //lik1
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    thetaA[p][pop][i]=save[p][pop][i]-h;
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik1[p][pop]=this.likelihood[p][pop];
                }
            }
            
            //lik3
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    thetaA[p][pop][i]=save[p][pop][i]+h;
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik3[p][pop]=this.likelihood[p][pop];
                }
            }
            
            
            //grad
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    if (Math.abs((lik3[p][pop]+lik1[p][pop]-2*lik2[p][pop])/(h*h))==0){
                        grad[p][pop]=((lik3[p][pop]-lik1[p][pop])*(2*h));
                    }
                    else{
                        grad[p][pop]=((lik3[p][pop]-lik1[p][pop])/(2*h))/Math.abs((lik3[p][pop]+lik1[p][pop]-2*lik2[p][pop])/(h*h));
                    }
                    
                }
            }
            
            
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    thetaA[p][pop][i]=save[p][pop][i]-grad[p][pop];
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik[p][pop]=this.likelihood[p][pop];
                }
            }
            
            //go back to previous value if it is not worse
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    if (lik[p][pop]>lik2[p][pop]){
                        thetaA[p][pop][i]=save[p][pop][i];
                    }
                }
            }
            
            
        }
    }
    
    
    
    
    
    
    //gradient descent for photons numbers
    void minimizeB(){
        
        double h=.01;
        
        double [][] lik=new double [nbProcess][popsize];
        double [][] lik1=new double [nbProcess][popsize];
        double [][] lik2=new double [nbProcess][popsize];
        double [][] lik3=new double [nbProcess][popsize];
        double [][] grad=new double [nbProcess][popsize];
        for (int p=0;p<this.nbProcess;p++){
            for (int pop=0;pop<this.popsize;pop++){
                lik[p][pop]=0;
                lik1[p][pop]=0;
                lik2[p][pop]=0;
                lik3[p][pop]=0;
            }
        }
        
        //here, it should be enough to update only the first background (except if it has to be negative  due to large other values)
        for (int i=0;i<1;i++){
            
            //lik2
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    save[p][pop][i]=thetaB[p][pop][i];
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik2[p][pop]=this.likelihood[p][pop];
                }
            }
            
            //lik1
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    thetaB[p][pop][i]=save[p][pop][i]-h;
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik1[p][pop]=this.likelihood[p][pop];
                }
            }
            
            //lik3
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    thetaB[p][pop][i]=save[p][pop][i]+h;
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik3[p][pop]=this.likelihood[p][pop];
                }
            }
            
            
            //grad
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    if (Math.abs((lik3[p][pop]+lik1[p][pop]-2*lik2[p][pop])/(h*h))==0){
                        grad[p][pop]=((lik3[p][pop]-lik1[p][pop])*(2*h));
                    }
                    else{
                        grad[p][pop]=((lik3[p][pop]-lik1[p][pop])/(2*h))/Math.abs((lik3[p][pop]+lik1[p][pop]-2*lik2[p][pop])/(h*h));
                    }
                    
                }
            }
            
            
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    thetaB[p][pop][i]=save[p][pop][i]-grad[p][pop];
                }
            }
            this.computeLikelihoodWithoutPSFupdating();
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    lik[p][pop]=this.likelihood[p][pop];
                }
            }
            
            //go back to previous value if it is not worse
            for (int p=0;p<this.nbProcess;p++){
                for (int pop=0;pop<popsize;pop++){
                    if (lik[p][pop]>lik2[p][pop]){
                        thetaB[p][pop][i]=save[p][pop][i];
                    }
                }
            }
            
            
        }
    }
    
    
    
    public void mutation(){
        mutation(1);
    }
    public void mutation(double percentage){
        
        for (int pop=0;pop<popsize;pop++){
            
            if (Math.random()*100<percentage){//mutation
                
                for (int i=0;i<this.numberPSFperModel;i++){
                    for (int p=0;p<this.nbProcess;p++){
                        double rA=this.thetaA[p][pop][i]*probaPhotonChange;
                        this.thetaA[p][pop][i]=this.thetaA[p][pop][i]+Math.random()*rA-rA/2;
                        double rB=this.thetaB[p][pop][i]*1;
                        this.thetaB[p][pop][i]=meanB[p][i]+Math.random()*rB-rB/2;
                    }
                    
                    this.thetaX[pop][i]=this.thetaX[pop][i]+(Math.random()*amplitudeMutation)-amplitudeMutation/2;
                    this.thetaY[pop][i]=this.thetaY[pop][i]+(Math.random()*amplitudeMutation)-amplitudeMutation/2;
                    this.thetaZ[pop][i]=this.thetaZ[pop][i]+(Math.random()*amplitudeMutation)-amplitudeMutation/2;
                }
                
                
                
            }
            
            
            
        }
        
        
        
    }
        
    
    
    public void makeCouples(){
        makeCouples(20);//make couple from 20 best % of the population
    }
    public void makeCouples(double percentage){
        
        int maxNumber=(int)(popsize*percentage/100);
        maxNumber=Math.max(maxNumber, 5);//at least 5
        maxNumber=Math.min(maxNumber, popsize-1);//at max (popsize-1)
        
        sortedLikelihood[0][2]=1;
        for (int pop=1;pop<popsize;pop++){
            //reversed order// +normalization (x-min/(max-min)) in [0,1] //+cumul
            sortedLikelihood[pop][2]=sortedLikelihood[pop-1][2]+(-sortedLikelihood[pop][0]-(-sortedLikelihood[maxNumber][0]))/(-sortedLikelihood[0][0]-(-sortedLikelihood[maxNumber][0]));
        }
        
        for (int pop=0;pop<popsize/2;pop++){
            
            double randHusband=Math.random()*sortedLikelihood[maxNumber][2];
            int indexHusband=0;
            
            while (randHusband>sortedLikelihood[indexHusband][2]){
                indexHusband++;
            }
            
            double randWife=Math.random()*sortedLikelihood[maxNumber][2];
            int indexWife=0;
            while (randWife>sortedLikelihood[indexWife][2]){
                indexWife++;
            }
            
            if (indexWife==indexHusband){
                looper:for (int u=0;u<10;u++){//try 10 times to take a wife that is different than husband
                    randWife=Math.random()*sortedLikelihood[maxNumber][2];
                    indexWife=0;
                    while (randWife>sortedLikelihood[indexWife][2]){
                        indexWife++;
                    }
                    if (indexWife!=indexHusband){
                        break looper;
                    }
                }
            }
            
            
            
            for (int i=0;i<this.numberPSFperModel;i++){
                for (int p=0;p<this.nbProcess;p++){
                    thetaAhusband[p][pop][i]=this.thetaA[p][(int)sortedLikelihood[indexHusband][1]][i];
                    thetaBhusband[p][pop][i]=this.thetaB[p][(int)sortedLikelihood[indexHusband][1]][i];
                    thetaAwife[p][pop][i]=this.thetaA[p][(int)sortedLikelihood[indexWife][1]][i];
                    thetaBwife[p][pop][i]=this.thetaB[p][(int)sortedLikelihood[indexWife][1]][i];
                }
                thetaXhusband[pop][i]=this.thetaX[(int)sortedLikelihood[indexHusband][1]][i];
                thetaYhusband[pop][i]=this.thetaY[(int)sortedLikelihood[indexHusband][1]][i];
                thetaZhusband[pop][i]=this.thetaZ[(int)sortedLikelihood[indexHusband][1]][i];
                thetaXwife[pop][i]=this.thetaX[(int)sortedLikelihood[indexWife][1]][i];
                thetaYwife[pop][i]=this.thetaY[(int)sortedLikelihood[indexWife][1]][i];
                thetaZwife[pop][i]=this.thetaZ[(int)sortedLikelihood[indexWife][1]][i];
            }
        }
        
            
            
        
    }
    
    
    
    
    
    
    public void makeNewPopulationSelection(){
        for (int i=0;i<this.numberPSFperModel;i++){
            for (int pop=0;pop<popsize;pop++){
                for (int p=0;p<this.nbProcess;p++){
                    this.thetaA[p][pop][i]=meanA[p][i]+random.nextGaussian()*stdA[p][i];
                    this.thetaB[p][pop][i]=meanB[p][i]+random.nextGaussian()*stdB[p][i];
                }
                this.thetaX[pop][i]=meanX[i]+random.nextGaussian()*stdX[i];
                this.thetaY[pop][i]=meanY[i]+random.nextGaussian()*stdY[i];
                this.thetaZ[pop][i]=meanZ[i]+random.nextGaussian()*stdZ[i];
                
            }
        }
        this.testPhotonNumber();
    }
    
    
    
    
    public void setParameters(int idParam){
        
        if (this.withregistration){
            if (idParam==0){
                for (int pop=0;pop<popsize;pop++){
                    dparam[idParam].modelMany.setParameters(id[idParam][pop], thetaX[pop], thetaY[pop],this.thetaZfocus[idParam], thetaZ[pop], thetaA[idParam][pop], thetaB[idParam][pop]);
                }
            }
            else{
                
                
                
                for (int pop=0;pop<popsize;pop++){
                    for (int j=0;j<this.numberPSFperModel;j++){
                        //first, reverse x & y and apply the drift correction of cam 1 (necessary for registration)
                        x2[j]=(-(thetaX[pop][j])-this.dx1);
                        y2[j]=(-(thetaY[pop][j])-this.dy1);
                        z2[j]=thetaZ[pop][j]-this.dz1;

                        //registration: cam1 -> cam2

                        v[0]=positionRefX+x2[j];
                        v[1]=positionRefY+y2[j];
                        v[2]=z2[j];
                        double [] res=pf.transform(v);
                        x2[j]=res[0]-positionRefX;
                        y2[j]=res[1]-positionRefY;
                        z2[j]=res[2];


                        double xs=x2[j]+this.dx2;


                        //remove drift correction of second cam
                        x2[j]=x2[j]+this.dx2;
                        y2[j]=y2[j]+this.dy2;
                        z2[j]=z2[j]+this.dz2;



                        //finally, shift according to the window position and reverse x & y
                        x2[j]=-(x2[j]-this.decX);
                        y2[j]=-(y2[j]-this.decY);
                    }
                    
                    dparam[idParam].modelMany.setParameters(id[idParam][pop], x2, y2,this.thetaZfocus[idParam], z2, thetaA[idParam][pop], thetaB[idParam][pop]);
                    
                }
            }


        }
        else{
            for (int pop=0;pop<popsize;pop++){
                dparam[idParam].modelMany.setParameters(id[idParam][pop], thetaX[pop], thetaY[pop],this.thetaZfocus[idParam], thetaZ[pop], thetaA[idParam][pop], thetaB[idParam][pop]);
            }
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        
    //necessary to launch each cam in parallel
    //obligatoire car la minimization de X, Y, et Z pour chaque camera n'est pas appelé au meme moment et on a des pb de synchro sans thread
    public class ComputeLikelihood extends Thread{
        
        
        int idProcess;
        int idpop;
        
        Object monitor1;
        
        boolean killer=false;
        
        ComputeLikelihood(Object monitor1,int idProcess,int idpop){
            this.idpop=idpop;
            this.idProcess=idProcess;
            this.monitor1=monitor1;
            
        }
        
        
        
        
        public void kill(){
            killer=true;
        }
        
        
        
        
        public void run(){
            synchronized(monitor1) {
                monitor1.notify();
                
                loop:while (true){
                
                    
                    try{ 
                        monitor1.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                    
                    if (killer){
                        break loop;
                    }
                    
                    
                    
                    double lik=dparam[idProcess].modelMany.getLikelihood(id[idProcess][idpop]);
                    likelihood[idProcess][idpop]=lik;
                    
//                    try{
//                        Thread.sleep(30);//To optimize
//                    }catch(Exception ee){}
//                    monitor1.notify();
                    
                    synchronized(monitor0) {
                        cpt++;
                        //IJ.log("cpt"+cpt);
                        if (cpt==nbProcess*popsize){
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
    
    
    
    
    
        
        
    //necessary to launch each cam in parallel
    //obligatoire car la minimization de X, Y, et Z pour chaque camera n'est pas appelé au meme moment et on a des pb de synchro sans thread
    public class ComputeLikelihoodWithoutPSFupdating extends Thread{
        
        
        int idProcess;
        int idpop;
        
        Object monitor2;
        
        boolean killer=false;
        
        ComputeLikelihoodWithoutPSFupdating(Object monitor2,int idProcess,int idpop){
            this.idpop=idpop;
            this.idProcess=idProcess;
            this.monitor2=monitor2;
            
        }
        
        
        
        
        public void kill(){
            killer=true;
        }
        
        
        
        
        public void run(){
            synchronized(monitor2) {
                monitor2.notify();
                
                loop:while (true){
                
                    
                    try{ 
                        monitor2.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                    
                    if (killer){
                        break loop;
                    }
                    
                    
                    
                    double lik=dparam[idProcess].modelMany.getLikelihoodWithoutPSFupdating(id[idProcess][idpop]);
                    likelihood[idProcess][idpop]=lik;
                    
                    
                    synchronized(monitor02) {
                        cpt++;
                        //IJ.log("cpt"+cpt);
                        if (cpt==nbProcess*popsize){
                            monitor02.notify();
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
