/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.utils.*;
import java.util.ArrayList;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import java.util.Arrays;
/**
 *
 * @author benoit
 */
public class DriftCorrectionFiducial {
    
    Object monitor0;
    int cpt=0;
        
    float epsilon=(float)0.0001;
        
    StackLocalization sl1;
    
    int iterations;
    
    double maxDistanceMergingXY;
    double maxDistanceMergingZ;
    
    int [][] idLoc;
    
    int smoothFrameNumber=100;
    
    double sizePix=100;//nm
    
    double maxX=Double.MIN_VALUE;
    double minX=Double.MAX_VALUE;
    double maxY=Double.MIN_VALUE;
    double minY=Double.MAX_VALUE;
       
    int offFrame=0;
    
    int thresholdConcecutiveFrame;
    
    //boolean toBeblocked=true;
    
    
    double initDriftX=0;
    double initDriftY=0;
    double initDriftZ=0;
    
    
    String pathDrift="";
    
    public DriftCorrectionFiducial(StackLocalization sl1,double maxDistanceMergingXY,double maxDistanceMergingZ,int thresholdConcecutiveFrame,int iterations,int smoothingFrameNumber,String path,String pathDrift,int offFrame){
        this.offFrame=offFrame;
        this.pathDrift=pathDrift;
        this.smoothFrameNumber=smoothingFrameNumber;
        this.iterations=iterations;
        this.sl1=sl1;
        
        this.thresholdConcecutiveFrame=thresholdConcecutiveFrame;
        this.maxDistanceMergingXY=maxDistanceMergingXY;
        this.maxDistanceMergingZ=maxDistanceMergingZ;
        
        
        if (path.length()>1){
            try{
                StackLocalization sl_attach=new StackLocalization(path);
                int maxFrame=-1;
                int idF=-1;
                int idPL=-1;
                for (int i=0;i<sl_attach.fl.size();i++){
                    FrameLocalization ffl=sl_attach.fl.get(i);
                    for (int ii=0;ii<ffl.loc.size();ii++){
                        PLocalization pl=ffl.loc.get(ii);
                        if (pl.frame>maxFrame){
                            maxFrame=pl.frame;
                            idPL=ii;
                            idF=i;
                        }
                    }
                }
                if (idPL>=0){
                    initDriftX=sl_attach.fl.get(idF).loc.get(idPL).drift_X;
                    initDriftY=sl_attach.fl.get(idF).loc.get(idPL).drift_Y;
                    initDriftZ=sl_attach.fl.get(idF).loc.get(idPL).drift_Z;
                }
            }
            catch(Exception eeee){}
        }
        
        
        
    }
    
    
    public StackLocalization run(){
        ChainList [][] chain = ChainList.getChainList(sl1,maxDistanceMergingXY,maxDistanceMergingZ,offFrame);
        
        driftCorrection(chain);
        
        return sl1;
    }
    
    
    public void driftCorrection(ChainList [][] chain){
        
        
        
        
        int maxFrame=Integer.MIN_VALUE;
        int minFrame=Integer.MAX_VALUE;
        for (int i=0;i<sl1.fl.size();i++){
            if (sl1.fl.get(i).numFrame>maxFrame){
                maxFrame=sl1.fl.get(i).numFrame;
            }
            if (sl1.fl.get(i).numFrame<minFrame){
                minFrame=sl1.fl.get(i).numFrame;
            }
        }
        int frameNumber=maxFrame-minFrame+1;
        if (frameNumber<thresholdConcecutiveFrame){
            IJ.log("WARNING: impossible to find beads, 'frame number' ("+frameNumber+") should be larger than 'min. consecutive frame' parameter ("+thresholdConcecutiveFrame+")");
            IJ.log("WARNING: Here, 'min. consecutive frame' is set to frame number automatically");
            IJ.log("WARNING: Please, reduce 'min. consecutive frame' parameter ("+thresholdConcecutiveFrame+") if the result is not good");
            thresholdConcecutiveFrame=frameNumber;
        }
        
        int curveNumber=0;
        
        for (int j=0;j<sl1.fl.size();j++){
            
            FrameLocalization fl=new FrameLocalization(sl1.fl.get(j).numFrame);
            for (int ii=0;ii<sl1.fl.get(j).loc.size();ii++){
                if (!chain[j][ii].foundPrevious){
                    double number=0;
                    boolean is=chain[j][ii].foundNext;
                    int idF=chain[j][ii].frameNext;
                    int idL=chain[j][ii].locNext;
                    while (is){
                        number++;
                        is=chain[idF][idL].foundNext;
                        if (is){
                            ChainList cl = chain[idF][idL];
                            idF=cl.frameNext;
                            idL=cl.locNext;
                        }
                    }
                    if (((number+1)>=thresholdConcecutiveFrame)&&(number>0)){
                        curveNumber++;
                    }
                }
            }
        }
        
        
        
        int [] idFrameReverse=new int[maxFrame+1];//this variable returns the position of frame with frame number as input... (necessary because not ordered)
        for (int i=0;i<maxFrame+1;i++){
            idFrameReverse[i]=-1;
        }
        for (int i=0;i<sl1.fl.size();i++){
            idFrameReverse[sl1.fl.get(i).numFrame]=i;
        }
        
//        ArrayList<Double> [] X = new ArrayList[sl1.fl.size()];
//        
//        for (int i=0;i<sl1.fl.size();i++){
//            X[i]=new ArrayList<Double>();
//        }
        
        
        
//        IJ.log("curve number "+curveNumber);
//        IJ.log("frame number "+frameNumber);
        
//        Boolean [][] kro = new Boolean[curveNumber][frameNumber];
        float [] d = new float [frameNumber];//drift at each frame
        float [] r = new float [curveNumber];//curve registration
        
        
        
        float [][] C = new float[curveNumber][];
        int [] Cfold = new int [curveNumber];
        int [][] Cf = new int [curveNumber][];
        
//        ArrayList<Double> [] CX = new ArrayList[curveNumber];
//        ArrayList<Double> [] CY = new ArrayList[curveNumber];
//        ArrayList<Double> [] CZ = new ArrayList[curveNumber];
        
        
        //correctDrift(chain,C,Cf,d,r,minFrame,maxFrame,"X");
        IJ.log("drift correction started ; frame number "+frameNumber+" ; curve number "+curveNumber);
        
        
        
        
        correctDrift(chain,C,Cf,d,r,minFrame,maxFrame,"X",iterations);
        
        
        correctDrift(chain,C,Cf,d,r,minFrame,maxFrame,"Y",iterations);
        
        
        correctDrift(chain,C,Cf,d,r,minFrame,maxFrame,"Z",iterations);
        
        
        
        if (this.pathDrift.length()>1){
            double [][] tableDrift= new double [sl1.fl.size()][4];
            for (int i=0;i<sl1.fl.size();i++){
                if (sl1.fl.get(i).loc.size()>0){
                    PLocalization pp = sl1.fl.get(i).loc.get(0);
                    tableDrift[i][0]=pp.frame;
                    tableDrift[i][1]=pp.drift_X;
                    tableDrift[i][2]=pp.drift_Y;
                    tableDrift[i][3]=pp.drift_Z;
                }
                
                
            }
            FileVectorLoader.saveTableInFile(pathDrift, tableDrift, ",");
        }
        
        IJ.log("drift correction finished");
        
        
        
        
        //correctDrift(chain,C,Cf,d,r,minFrame,maxFrame,"Y",iterations);
        
        //correctDrift(chain,C,Cf,d,r,minFrame,maxFrame,"Z");
        
        
        /*plot(ffx,fdx,"Z drift","frame","X drift");
        
        plot(ffy,fdy,"Z drift","frame","Y drift");
        
        plot(ffz,fdz,"Z drift","frame","Z drift");
        */
        
        
        IJ.showProgress(0);
        
    }
    
    
    
    
    
    private void correctDrift(ChainList [][] chain,float [][] C,int [][] Cf,float [] d, float [] r,int minFrame,int maxFrame,String variable,int iter){
        int curveNumber=C.length;
        int frameNumber=maxFrame-minFrame+1;
        
        int [] positOfFrame=new int[frameNumber];
        
        int idcurveNumber=0;
        
        ArrayList<PLocalization> p_next = new ArrayList<PLocalization>();
        
        for (int j=0;j<sl1.fl.size();j++){
            FrameLocalization fl=new FrameLocalization(sl1.fl.get(j).numFrame);
            positOfFrame[sl1.fl.get(j).numFrame-minFrame]=j;
            for (int ii=0;ii<sl1.fl.get(j).loc.size();ii++){
                if (!chain[j][ii].foundPrevious){
                    PLocalization p_current=sl1.fl.get(j).loc.get(ii).copy();
                    
                    int number=0;
                    boolean is=chain[j][ii].foundNext;
                    int idF=chain[j][ii].frameNext;
                    int idL=chain[j][ii].locNext;
                    p_next.clear();
                    while (is){
                        
                        p_next.add(sl1.fl.get(idF).loc.get(idL));
                        number++;
                        is=chain[idF][idL].foundNext;
                        if (is){
                            ChainList cl = chain[idF][idL];
                            idF=cl.frameNext;
                            idL=cl.locNext;
                        }
                    }
                    if (((number+1)>=thresholdConcecutiveFrame)&&(number>0)){
                        C[idcurveNumber]=new float[number+1];
                        Cf[idcurveNumber]=new int[number+1];
                        Cf[idcurveNumber][0]=p_current.frame;
                        //Cf[idcurveNumber]=p_current.frame;
                        if (variable.startsWith("X")){
                            C[idcurveNumber][0]=new Float((float)p_current.X);
                        }
                        else if (variable.startsWith("Y")){
                            C[idcurveNumber][0]=new Float((float)p_current.Y);
                        }
                        else if (variable.startsWith("Z")){
                            C[idcurveNumber][0]=new Float((float)p_current.Z);
                        }
                        
                        
                        for (int u=0;u<p_next.size();u++){
                            Cf[idcurveNumber][u+1]=p_next.get(u).frame;
                            if (variable.startsWith("X")){
                                C[idcurveNumber][u+1]=new Float((float)p_next.get(u).X);
                            }
                            else if (variable.startsWith("Y")){
                                C[idcurveNumber][u+1]=new Float((float)p_next.get(u).Y);
                            }
                            else if (variable.startsWith("Z")){
                                C[idcurveNumber][u+1]=new Float((float)p_next.get(u).Z);
                            }
                            
                        }
                        idcurveNumber++;
                    }
                }
            }
        }
        
        IJ.log("curve numb "+curveNumber+"  "+frameNumber);
        
        float [] curveNumberI=new float[frameNumber];
        float [] curveLengthJ=new float[curveNumber];
        boolean [] curveAccounted=new boolean[curveNumber];
        
        for (int i=0;i<frameNumber;i++){
            curveNumberI[i]=0;
            d[i]=0;
        }
        double sumtot=0;
        for (int j=0;j<curveNumber;j++){
            curveAccounted[j]=true;
            curveLengthJ[j]=C[j].length;
            sumtot+=curveLengthJ[j];
            r[j]=0;
        }
        for (int j=0;j<curveNumber;j++){
            for (int i=0;i<C[j].length;i++){
                
                //curveNumberI[Cf[j]+i-minFrame]++;
                curveNumberI[Cf[j][i]-minFrame]++;
            }
            
        }
        
        
        double Jprev=Double.POSITIVE_INFINITY;
        double [] copyd = new double [frameNumber];
        int iterInitR=5;
        int nbThread=16;
        
        ParallelComputing [] pc = new ParallelComputing[nbThread];
        Object [] monitor1 = new Object[nbThread];
        monitor0 = new Object();
        int curveNumberSplit=(int)Math.ceil((double)curveNumber/(double)(nbThread));
        int frameNumberSplit=(int)Math.ceil((double)frameNumber/(double)(nbThread));
        for (int k=0;k<nbThread;k++){
            int startI=frameNumberSplit*k;
            int stopI=frameNumberSplit*(k+1);
            if (startI>=frameNumber){
                startI=-1;
                stopI=-1;
            }
            int startJ=curveNumberSplit*k;
            int stopJ=curveNumberSplit*(k+1);
            if (startJ>=curveNumber){
                startJ=-1;
                stopJ=-1;
            }
            if (k==nbThread-1){
                if (stopI!=-1){
                    stopI=frameNumber;
                }
                if (stopJ!=-1){
                    stopJ=curveNumber;
                }
            }
            if (startI>=frameNumber){
                startI=-1;
                stopI=-1;
            }
            if (startJ>=curveNumber){
                startJ=-1;
                stopJ=-1;
            }
            if (stopI>frameNumber){
                stopI=frameNumber;
            }
            if (stopJ>curveNumber){
                stopJ=curveNumber;
            }
            
            
            monitor1[k] = new Object();
            pc[k] = new ParallelComputing(C,Cf,d,r,curveAccounted,minFrame,curveNumberI,curveLengthJ,startJ,stopJ,startI,stopI,monitor1[k],nbThread);
        
        }
        
        
        
        for (int k=0;k<nbThread;k++){
            synchronized(monitor1[k]){
                pc[k].start();
                try{
                    monitor1[k].wait();//to be sure the threads are launched at this stage and that wait() in the thread is called
                }catch(Exception er){IJ.log("oops: problem wait function for sync monitor 1");}
            }
        }
        
        
        boolean firstPass=true;
        loop:for (int t=0;t<iter;t++){
            
            
            if (firstPass){
                if (variable.startsWith("X")){
                    IJ.showProgress((float)t*.5*.33/(float)iter);
                }
                else if (variable.startsWith("Y")){
                    IJ.showProgress(.33+(float)t*.5*.33/(float)iter);
                }
                else if (variable.startsWith("Z")){
                    IJ.showProgress(.66+(float)t*.5*.33/(float)iter);
                }
                
            }
            else{
                if (variable.startsWith("X")){
                    IJ.showProgress(.166+(float)t*.5*.33/(float)iter);
                }
                else if (variable.startsWith("Y")){
                    IJ.showProgress(.33+.166+(float)t*.5*.33/(float)iter);
                }
                else if (variable.startsWith("Z")){
                    IJ.showProgress(.66+.166+(float)t*.5*.33/(float)iter);
                }
            }
            
            
            
            
        
            double [] ff = new double [frameNumber];
            int nbok=0;
            boolean [] ok = new boolean [frameNumber];
            for (int i=0;i<frameNumber;i++){
                ff[i]=i+minFrame;
                if (curveNumberI[i]>0){
    //                IJ.log(""+ff[i]+","+d[i]);
                    ok[i]=true;
                    nbok++;
                }
                else{
                    ok[i]=false;
                }
            }
            
            //registration curve update://///////////////////////////***********************************
            
            int loopNumber=1;
            if (t==0){
                loopNumber=iterInitR;
            }
            
            for (int tt=0;tt<loopNumber;tt++){
                
                for (int j=0;j<curveNumber;j++){
                    //IJ.log("reg "+tt+"  "+j+"  "+r[j]);
                }

                synchronized(monitor0){
                    cpt=0;
                
                
                    for (int k=0;k<nbThread;k++){
                        synchronized(monitor1[k]){
                            //toBeblocked=true;
                            pc[k].switchToR();
                            monitor1[k].notify();
                        }
                    }
//                    for (int k=0;k<nbThread;k++){
//                        synchronized(monitor1[k]){
//                            try{
//                                if (toBeblocked){
//                                    monitor1[k].wait();
//                                }
//                            }catch(Exception ee){IJ.log("error wait function "+ee);}
//                        }
//                    }


                
                    try{
                        monitor0.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                }
                
            }
            
            
            //drift update://///////////////////////////***********************************
            synchronized(monitor0){
                cpt=0;
                for (int k=0;k<nbThread;k++){
                    synchronized(monitor1[k]){
                        //toBeblocked=true;
                        pc[k].switchToD();
                        monitor1[k].notify();
                    }
                }
//            for (int k=0;k<nbThread;k++){
//                synchronized(monitor1[k]){
//                    try{
//                        if (toBeblocked){
//                            monitor1[k].wait();
//                        }
//                    }catch(Exception ee){IJ.log("error wait function "+ee);}
//                }
//            }
                try{
                    monitor0.wait();
                }catch(Exception ee){IJ.log("error wait function "+ee);}
            }
            
            
            
        
            


            //smoothing of D:
            
            

            for (int i=0;i<frameNumber;i++){
                int count=0;
                copyd[i]=0;
                for (int a=-smoothFrameNumber/2;a<=+smoothFrameNumber/2;a++){
                    if ((a+i>=0)&&(a+i<frameNumber)){
                        if (ok[a+i]){
                            copyd[i]+=d[i+a];
                            count++;
                        }

                    }
                }
                if (count>0){
                    copyd[i]/=count;
                }
            }
            for (int i=0;i<frameNumber;i++){
                d[i]=(float)copyd[i];
            }
            
            
            
            
            //criterionComputation/////////////////////////////***********************************
            synchronized(monitor0){
                cpt=0;
                for (int k=0;k<nbThread;k++){
                    synchronized(monitor1[k]){
                        //toBeblocked=true;
                        pc[k].switchToJ();
                        monitor1[k].notify();
                    }
                }
//                for (int k=0;k<nbThread;k++){
//                    synchronized(monitor1[k]){
//                        try{
//                            if (toBeblocked){
//                                monitor1[k].wait();
//                            }
//                        }catch(Exception ee){IJ.log("error wait function "+ee);}
//                    }
//                }
                try{
                    monitor0.wait();
                }catch(Exception ee){IJ.log("error wait function "+ee);}
            }
            double J=0;
            for (int k=0;k<nbThread;k++){
                J+=pc[k].getJ();
            }
            J/=sumtot;
            J=Math.sqrt(J);
            //IJ.log("iteration "+(iter-t)+"  ; drift of "+variable+"    mean distance:"+J);
            
            
            
            
            
            
            
            
            
            
            if ((Math.abs(Jprev-J)<epsilon)||(t==(iter-1))){
                if (firstPass){
                    t=0;
                    double meanError;
                    for (int j=0;j<curveNumber;j++){
                        meanError=0;
                        for (int i=0;i<C[j].length;i++){
                            //meanError+=(double)(d[i+Cf[j]-minFrame]-(C[j][i]-r[j]))*(double)(d[i+Cf[j]-minFrame]-(C[j][i]-r[j]));
                            meanError+=(double)(d[Cf[j][i]-minFrame]-(C[j][i]-r[j]))*(double)(d[Cf[j][i]-minFrame]-(C[j][i]-r[j]));
                        }
                        meanError/=(double)C[j].length;
                        meanError=Math.sqrt(meanError);
                        if (meanError>3*J){
                            curveAccounted[j]=false;
                        }
                        else{
                            curveAccounted[j]=true;
                        }
                    }
                    Jprev=Double.POSITIVE_INFINITY;
                    firstPass=false;
                }
                else{
                    break loop;
                }
            }
            
            
            Jprev=J;
            
        }
        
        for (int k=0;k<nbThread;k++){
            pc[k].kill();
            synchronized(monitor1[k]){
                monitor1[k].notify();
            }
            try{
                pc[k].join();
            }catch(Exception eee){IJ.log("Thread impossible to join() "+eee);}
        }
        
        
        
        
        
        ////////////////////////////////////////////////////////////////////////TO REMOVE*********************************************
        
        /*loop:for (int t=0;t<iter;t++){
            IJ.log("iteration "+(iter-t)+"  ; drift of "+variable);
            
            
            
            
            //registration curve update://///////////////////////////***********************************
            int loopNumber=1;
            if (t==0){
                loopNumber=iterInitR;
            }
            for (int tt=0;tt<loopNumber;tt++){    
                
                for (int j=0;j<curveNumber;j++){
                    float sum=0;
                    for (int i=0;i<C[j].length;i++){
                        sum+=(C[j][i]-d[i+Cf[j]-minFrame]);
                    }
                    if (curveLengthJ[j]>0){
                        r[j]=sum/curveLengthJ[j];
                    }
                    else{
                        r[j]=-1;
                    }
                }
                

            }
            
            
            //drift update://///////////////////////////***********************************
            for (int i=0;i<frameNumber;i++){
                float sum=0;
                for (int j=0;j<curveNumber;j++){
                    if ((i>=Cf[j]-minFrame)&&(i<Cf[j]+C[j].length-minFrame)){
                        sum+=(C[j][i+minFrame-Cf[j]]-r[j]);
                    }
                }
                if (curveNumberI[i]>0){
                    d[i]=sum/curveNumberI[i];
                }
                else{
                    d[i]=-1;
                }
            }

            
            //criterionComputation/////////////////////////////***********************************
            double J=0;
            for (int j=0;j<curveNumber;j++){
                for (int i=0;i<C[j].length;i++){
                    J+=(double)(d[i+Cf[j]-minFrame]-(C[j][i]-r[j]))*(double)(d[i+Cf[j]-minFrame]-(C[j][i]-r[j]));
                }
            }
            IJ.log("                     J "+J);
            if ((Math.abs(Jprev-J)<epsilon)){
                break loop;
            }
            Jprev=(float)J;
            
        }*/
        
        
        
        
        
        //////////////////////////////////////////////////////////POST PROCESSING - filtering bad data
        
        
        
        
        
        
        
            
        
        
        
        
        
        
        
        
        
        
        
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        ArrayList<Integer> alstarthole = new ArrayList();
        ArrayList<Integer> alstophole = new ArrayList();
        ArrayList<Integer> alstart = new ArrayList();
        ArrayList<Integer> alstop = new ArrayList();
        for (int i=0;i<frameNumber;i++){
            curveNumberI[i]=0;//recompute curve Number dealing with curve removed
        }
        for (int j=0;j<curveNumber;j++){
            if (curveAccounted[j]){
                for (int i=0;i<C[j].length;i++){
                    curveNumberI[Cf[j][i]-minFrame]++;
                }
            }
            else{
                IJ.log("some curve removed... "+j);
            }
        }
        
        
        
        boolean nope=true;
        for (int i=0;i<frameNumber;i++){
            if (nope){
                if (curveNumberI[i]>0){
                    nope=false;
                    alstart.add(i);
                    
                }
            }
            else{
                if (curveNumberI[i]<=0){
                    nope=true;
                    alstop.add(i);
                }
            }
        }
        if (!nope){
            alstop.add(frameNumber);
        }
        
        int index=-1;
        if (curveNumberI[0]<=0){
            alstarthole.add(0);
            index++;
        }
        for (int i=1;i<frameNumber;i++){
            if (curveNumberI[i]<=0){
                if (curveNumberI[i-1]>0){
                    alstarthole.add(i);
                    index++;
                }
            }
            else{
                if (curveNumberI[i-1]<=0){
                    alstophole.add(i);
                }
            }
        }
        if (curveNumberI[frameNumber-1]<=0){
            alstophole.add(frameNumber);
        }
        
        
        
        
        double [] ff = new double [frameNumber];
        int nbok=0;
        boolean [] ok = new boolean [frameNumber];
        for (int i=0;i<frameNumber;i++){
            ff[i]=i+minFrame;
            if (curveNumberI[i]>0){
//                IJ.log(""+ff[i]+","+d[i]);
                ok[i]=true;
                nbok++;
            }
            else{
                ok[i]=false;
            }
        }
            
            
        
        
        double [] newD = new double [d.length];
        
        
        
        //we should deal with holes here.....
        if (1==2){
                                                            for (int i=0;i<alstarthole.size();i++){


                                                                int sizevect1=smoothFrameNumber/2;
                                                                sizevect1=Math.min(sizevect1,alstarthole.get(i)-1);

                                                                if (sizevect1>0){
                                                                    int startvect1=Math.max(0, alstarthole.get(i)-sizevect1);
                                                                    //IJ.log("Vect1 "+startvect1+"  "+sizevect1);
                                                                    double [][] vect1=new double[1][sizevect1];
                                                                    double [][] frame1=new double[1][sizevect1];
                                                                    //IJ.log("sizevect1  "+sizevect1+"  "+startvect1+"  "+alstart.get(i));
                                                                    for (int j=startvect1,k=0;k<sizevect1;j++,k++){
                                                                        vect1[0][k]=d[j];
                                                                        frame1[0][k]=ff[j];
                                                                        //IJ.log("vect 1 : "+frame1[0][k]+"  "+vect1[0][k]);
                                                                    }

                                                                    double a1=0;
                                                                    double b1=vect1[0][vect1[0].length-1];
                                                                    PolynomialFit pf1 = new PolynomialFit(1,vect1,frame1);
                                                                    if (sizevect1>1){
                                                                        boolean bool=pf1.run();
                                                                        if (bool){
                                                                            b1=pf1.a[0][0];
                                                                            a1=pf1.a[0][1];

                                                                        }
                                                                        else{
                                                                            //IJ.log("vect 1 size "+vect1.length);
                                                                            pf1.a=new double [vect1.length][2];
                                                                            pf1.a[0][0]=b1;
                                                                            pf1.a[0][1]=a1;
                                                                        }
                                                                    }
                                                                    else{
                                                                        pf1.a=new double [vect1.length][2];
                                                                        pf1.a[0][0]=0;
                                                                        pf1.a[0][1]=0;
                                                                        IJ.log("oops... size vect1 should be at least 2 length");
                                                                    }

                                                                    //IJ.log("1-ax+b "+variable+"  "+pf1.a[0][1]+"  "+pf1.a[0][0]);


                                                                    int sizevect2=smoothFrameNumber/2;
                                                                    sizevect2=Math.min(sizevect2,frameNumber-alstophole.get(i));
                                                                    if (i!=alstarthole.size()-1){//if other hole just after
                                                                        sizevect2=Math.min(sizevect2,alstarthole.get(i+1)-alstophole.get(i));
                                                                    }
                                                                    if (sizevect2==smoothFrameNumber/2){//if hole has full information after also (no hole after)
                                                                        int startvect2=alstophole.get(i);
                                                                        //IJ.log("Vect2 "+startvect2+"  "+sizevect2);
                                                                        double [][] vect2=new double[1][sizevect2];
                                                                        double [][] frame2=new double[1][sizevect2];
                                                                        for (int j=startvect2,k=0;k<sizevect2;j++,k++){
                                                                            vect2[0][k]=d[j];
                                                                            frame2[0][k]=ff[j];

                                                                        }
                                                                        double a2=0;
                                                                        double b2=vect2[0][0];
                                                                        if (sizevect2>1){
                                                                            PolynomialFit pf = new PolynomialFit(1,vect2,frame2);
                                                                            boolean bobool=pf.run();
                                                                            if (bobool){
                                                                                b2=pf.a[0][0];
                                                                                a2=pf.a[0][1];
                                                                            }
                                                                            else{
                                                                                //IJ.log("vect 2 size "+vect2.length);
                                                                            }
                                                                        }
                                                                        else{
                                                                            IJ.log("oops... size vect2 should be at least 2 length");
                                                                        }

                                                                        //IJ.log("2-ax+b "+variable+"  "+a2+"  "+b2);


                                                                        double meanSlope=(a1+a2)/2;
                                                                        double newB1=a1*frame1[0][vect1[0].length-1]+b1-(meanSlope*frame1[0][vect1[0].length-1]);

                                                                        pf1.a[0][0]=newB1;
                                                                        pf1.a[0][1]=meanSlope;
                                                    //                    IJ.log("truc : "+(a1*frame1[0][vect1[0].length-1]+b1));
                                                    //                    IJ.log("truc : "+(meanSlope*frame1[0][vect1[0].length-1]+newB1));
                                                                        //IJ.log("correction between "+alstart.get(i)+" and "+alstop.get(i));
                                                                        for (int k=alstarthole.get(i);k<alstophole.get(i);k++){
                                                                            double [] v = new double [1];
                                                                            v[0]=ff[k];
                                                                            v=pf1.transform(v);

                                                                            d[k]=(float)v[0];
                                                                            //IJ.log(""+ff[k]+"  "+d[k]);
                                                                        }
                                                                        double difference=(meanSlope*frame2[0][0]+newB1)-(a2*frame2[0][0]+b2);
                                                    //                    IJ.log("al stop "+alstop.get(i));
                                                                        for (int k=alstophole.get(i);k<frameNumber;k++){
                                                                            d[k]+=difference;
                                                                        }



                                                                    }
                                                                    else{//if hole depends on previous only
                                                                        for (int k=alstarthole.get(i);k<alstophole.get(i);k++){
                                                                            double [] v = new double [1];
                                                                            v[0]=ff[k];
                                                                            v=pf1.transform(v);

                                                                            d[k]=(float)v[0];
                                                                        }
                                                                    }
                                                                }
                                                                else{//if hole is at the begining of frames
                                                                    int sizevect2=smoothFrameNumber/2;
                                                                    sizevect2=Math.min(sizevect2,frameNumber-alstophole.get(i));
                                                                    if (i!=alstarthole.size()-1){//if other hole just after
                                                                        sizevect2=Math.min(sizevect2,alstarthole.get(i+1)-alstophole.get(i));
                                                                    }
                                                                    if (sizevect2>0){
                                                                        int startvect2=alstophole.get(i);
                                                                        //IJ.log("Vect2 "+startvect2+"  "+sizevect2);
                                                                        double [][] vect2=new double[1][sizevect2];
                                                                        double [][] frame2=new double[1][sizevect2];
                                                                        for (int j=startvect2,k=0;k<sizevect2;j++,k++){
                                                                            vect2[0][k]=d[j];
                                                                            frame2[0][k]=ff[j];
                                                                        }

                                                                        PolynomialFit pf2=null;

                                                                        if (sizevect2>1){
                                                                            pf2 = new PolynomialFit(1,vect2,frame2);
                                                                            boolean bool=pf2.run();
                                                                            if (bool){
                                                                            }
                                                                            else{
                                                                                pf2.a=new double [vect2.length][2];
                                                                                pf2.a[0][0]=0;
                                                                                pf2.a[0][1]=vect2[0][0];
                                                                            }
                                                                            //pf2.log();
                                                                        }
                                                                        else{
                                                                            IJ.log("oops... size vect2 should be at least 2 length");
                                                                        }



                                                                        for (int k=alstarthole.get(i);k<alstophole.get(i);k++){
                                                                            double [] v = new double [1];
                                                                            v[0]=ff[k];
                                                                            v=pf2.transform(v);

                                                                            d[k]=(float)v[0];
                                                                        }


                                                                    }
                                                                    else{
                                                                        IJ.log("WARNING: WRONG "+variable+" drift because 0 bead is found between frames "+(alstart.get(i)+minFrame)+" and "+(alstop.get(i)+minFrame));
                                                                    }

                                                                }





                                                                IJ.log("WARNING: "+variable+" drift is approximative between frames "+(alstart.get(i)+minFrame)+" and "+(alstop.get(i)+minFrame));



                                                            }
        }
        else if (1==2){
            
            double meanPrev=0;
            for (int i=0;i<alstart.size();i++){
                
                //IJ.log("alalalalalalalalal "+alstart.get(i)+"  "+alstop.get(i));
                
                int sizeVect=alstop.get(i)-alstart.get(i);
                sizeVect=Math.min(sizeVect, Math.max(smoothFrameNumber/2,1));
                
                if ((i==0)&&(alstart.get(i)>0)){
                    double mean=0;
                    double count=0;
                    double [] vv=new double [sizeVect-alstart.get(i)];
                    for (int u=alstart.get(i),k=0;u<sizeVect;u++){
                        vv[k++]=d[u];
                        mean+=d[u];
                        count++;
                    }
                    Arrays.sort(vv);
                    mean/=count;
                    for (int u=0;u<alstart.get(i)-1;u++){
                        //d[u]=(float)mean;
                        newD[u]=(float)vv[vv.length/2];
                    }
                    IJ.log("WARNING: "+variable+" drift is approximative between frames "+(0)+" and "+(alstart.get(i)));

                }
                
                
                if (i>=1){//it means that we have holes
                    for (int u=alstop.get(i-1);u<alstart.get(i);u++){
                        d[u]=(float)meanPrev;
                    }
                    double [] vv=new double [alstart.get(i)+sizeVect-alstart.get(i)];
                    double mean=0;
                    double count=0;
                    for (int u=alstart.get(i),k=0;u<alstart.get(i)+sizeVect;u++){
                        vv[k++]=d[u];
                        mean+=d[u];
                        count++;
                    }
                    mean/=count;
                    Arrays.sort(vv);
                    for (int u=alstart.get(i);u<alstop.get(i);u++){
                        //d[u]-=(float)mean-meanPrev;
                        newD[u]=(float)(vv[vv.length/2]-meanPrev);
                    }
                    //IJ.log("mean "+meanPrev+"  "+mean);
                    IJ.log("WARNING: "+variable+" drift is approximative between frames "+(alstop.get(i-1))+" and "+(alstart.get(i)));
                }
                
                
                
                meanPrev=0;
                double count=0;
                double [] vv=new double [alstop.get(i)-(alstop.get(i)-sizeVect)];
                for (int u=alstop.get(i)-sizeVect,k=0;u<alstop.get(i);u++){
                    vv[k++]=d[u];
                    meanPrev+=d[u];
                    count++;
                }
                meanPrev/=count;
                Arrays.sort(vv);
                meanPrev=vv[vv.length/2];
                if (i==alstart.size()-1){
                    for (int u=alstop.get(i);u<d.length;u++){
                        newD[u]=(float)meanPrev;
                    }
                    IJ.log("WARNING: "+variable+" drift is approximative between frames "+(alstop.get(i))+" and "+(d.length));
                }
                
                
                
                
                
            }
            
            
            
                    
        }
        
        
        
        
        
        
        
        
        
        double [] fd = new double [frameNumber];
        
        
        //smoothing:
        for (int i=0;i<frameNumber;i++){
            
            
            
            int count=0;
            
            for (int a=-smoothFrameNumber/2;a<=+smoothFrameNumber/2;a++){
                if ((a+i>=0)&&(a+i<frameNumber)){
                    //if (ok[a+i]){
                        //fd[i]+=d[i+a];
                        count++;
                    //}

                }
            }
            double [] vv=new double [count];
            for (int a=-smoothFrameNumber/2,k=0;a<=+smoothFrameNumber/2;a++){
                if ((a+i>=0)&&(a+i<frameNumber)){
                    //if (ok[a+i]){
                        //fd[i]+=d[i+a];
                        vv[k++]=d[i+a];
                        
                    //}

                }
            }
            //fd[i]/=count;
            Arrays.sort(vv);
            fd[i]=vv[vv.length/2];
            fd[i]=d[i];
//            IJ.log(""+ff[i]+"  "+d[i]+"  "+fd[i]);
        }
        
        
        
        
        int nbPrint = Math.min(50000000, frameNumber);
        if (nbPrint>0){
            double [] fff=new double [nbPrint];
            double [] ddd=new double [nbPrint];
            boolean [] okok=new boolean [nbPrint];

            for (int i=0,k=0;i<frameNumber;i+=(frameNumber/nbPrint),k++){
                if (k<nbPrint){
                    fff[k]=ff[i];
                    if (variable.startsWith("X")){
                        ddd[k]=fd[i]+this.initDriftX;
                    }
                    else if (variable.startsWith("Y")){
                        ddd[k]=fd[i]+this.initDriftY;
                    }
                    else if (variable.startsWith("Z")){
                        ddd[k]=fd[i]+this.initDriftZ;
                    }
                    
                    okok[k]=ok[i];
                }
            }

            this.plot(fff, ddd,okok, ""+variable+" Drift", "frame number", ""+variable+" Drift (nm)");
        }
        else{
            IJ.log("WARNING: wrong "+variable+" drift because no bead found. Maybe, try again after reducing minimum frame number");
            
        }
        
        
        
        for (int i=0;i<frameNumber;i++){
            int posit=positOfFrame[i];
            for (int ii=0;ii<sl1.fl.get(posit).loc.size();ii++){
                PLocalization p=sl1.fl.get(posit).loc.get(ii);
                if (variable.startsWith("X")){
                    p.setDrift_X(fd[i]+this.initDriftX);
                }
                else if (variable.startsWith("Y")){
                    p.setDrift_Y(fd[i]+this.initDriftY);
                }
                else if (variable.startsWith("Z")){
                    p.setDrift_Z(fd[i]+this.initDriftZ);
                }
            }
        }
        
        
        
        
        
        
        
        
        
    }
    
    
    
    
    
        
        
        
        
        
        
        
        
    
    
    
    
    
    
    public void plot(double [] x, double [] y,boolean [] ok,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel);
        
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
        
        p.setLimits(xMin-200, xMax+200, yMin-200, yMax+200);
        p.setLineWidth(2);
        p.setFont(0, 18);
        
        
        int oknumber=0;
        int notoknumber=0;
        for (int i=0;i<ok.length;i++){
            if (ok[i]){
                oknumber++;
            }
            else{
                notoknumber++;
            }
        }
        
        if (oknumber>0){
            double [] xok = new double [oknumber];
            double [] yok = new double [oknumber];
            for (int i=0,k=0;i<ok.length;i++){
                if (ok[i]){
                    xok[k]=x[i];
                    yok[k]=y[i];
                    k++;
                }
            }
            p.setColor(Color.red);
            p.addPoints(xok, yok, Plot.CROSS);
        
        }
        if (notoknumber>0){
            double [] xnotok = new double [notoknumber];
            double [] ynotok = new double [notoknumber];
            for (int i=0,k=0;i<ok.length;i++){
                if (!ok[i]){
                    xnotok[k]=x[i];
                    ynotok[k]=y[i];
                    k++;
                }
            }
            p.setColor(Color.blue);
            p.addPoints(xnotok, ynotok, Plot.CROSS);
        }
        
        
        
        
        
        p.show();
        
    }
    
    
    
    public void plot(double [] x, double [] y,String title,String xlabel,String ylabel){
        
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
        p.setLimits(xMin-200, xMax+200, yMin-200, yMax+200);
        p.setLineWidth(2);
        p.setFont(0, 18);
        for (int ii=0;ii<x.length;ii++){
            
            p.setColor(Color.blue);
            p.add("CIRCLE",  x,y);
            //p.show();


            
        }
        
        p.show();
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class ParallelComputing extends Thread{
        
        int nbThread;
        float [][] C;
        //int [] Cf;
        int [][] Cf;
        float [] d; 
        float [] r;
        float [] curveNumberI;
        float [] curveLengthJ;
        int minFrame;
        boolean [] curveAccounted;
        int startJ;
        int stopJ;
        int startI;
        int stopI;
        
        double J=0;
                
        int computation=0;
        Object monitor1;
        boolean killer=false;
        ParallelComputing(float [][] C,int [][] Cf,float [] d, float [] r,boolean [] curveAccounted,int minFrame,float [] curveNumberI,float [] curveLengthJ,int startJ,int stopJ,int startI,int stopI,Object monitor1,int nbThread){
            this.C=C;
            this.Cf=Cf;
            this.d=d;
            this.r=r;
            this.minFrame=minFrame;
            this.curveNumberI=curveNumberI;
            this.curveLengthJ=curveLengthJ;
            this.startJ=startJ;
            this.stopJ=stopJ;
            this.startI=startI;
            this.stopI=stopI;
            this.monitor1=monitor1;
            this.curveAccounted=curveAccounted;
            this.nbThread=nbThread;
        }
        
        
        public double getJ(){
            return J;
        }
        
        public void switchToR(){
            computation=0;
        }
        public void switchToD(){
            computation=1;
        }
        public void switchToJ(){
            computation=2;
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
                    if (computation==0){
                        runR();
                    }
                    else if (computation==1){
                        runD();
                    }
                    else if (computation==2){
                        runJ();
                    }
                    try{
                        Thread.sleep(100);
                    }catch(Exception ee){}
                    
                    synchronized(monitor0) {
                        cpt++;
                        if (cpt==nbThread){
                            monitor0.notify();
                        }
                    }
//                    monitor1.notify();
//                    toBeblocked=false;
                
                    if (killer){
                        break loop;
                    }
                }


            }
        }
        
        
        void runR(){
            if (startJ!=-1){
                for (int j=startJ;j<stopJ;j++){
                    
                    float sum=0;
//                    if (j>=C.length){
//                        IJ.log("ERROR out of bounds "+j+"  "+stopJ+"  "+C.length);
//                    }
                    for (int i=0;i<C[j].length;i++){
                        //sum+=(C[j][i]-d[i+Cf[j]-minFrame]);
                        sum+=(C[j][i]-d[Cf[j][i]-minFrame]);
                        
                    }
                    
                    if (C[j].length>0){
                        r[j]=sum/C[j].length;
                    }
                    else{
                        r[j]=-1;
                    }
                }
            }
        }
        
        
        void runD(){
            if (startI!=-1){
                for (int i=startI;i<stopI;i++){
                    float sum=0;
                    float count=0;
                    for (int j=0;j<curveLengthJ.length;j++){
//                        if ((i>=Cf[j]-minFrame)&&(i<Cf[j]+C[j].length-minFrame)){
//                            if (curveAccounted[j]){
//                                sum+=(C[j][i+minFrame-Cf[j]]-r[j]);
//                                count++;
//                            }
//                        }
                        boolean ok=false;
                        int posit=0;
                        loloop:for (int ii=0;ii<Cf[j].length;ii++){
                            if (Cf[j][ii]-minFrame==i){
                                ok=true;
                                posit=ii;
                                break loloop;
                            }
                        }
                        if (ok&&curveAccounted[j]){
                            sum+=(C[j][posit]-r[j]);
                            count++;
                        }
                    }
                    
                    
                    
                    if (count>0){
                        d[i]=sum/count;
                    }
                    else{
                        d[i]=-1;
                    }
//                    if (curveNumberI[i]>0){
//                        d[i]=sum/curveNumberI[i];
//                    }
//                    else{
//                        d[i]=-1;
//                    }
                }
            }
        }
        
        
        
        void runJ(){
            J=0;
            
            if (startJ!=-1){
                for (int j=startJ;j<stopJ;j++){
                    if (curveAccounted[j]){
                        for (int i=0;i<C[j].length;i++){
                            //J+=(double)((double)(d[i+Cf[j]-minFrame]-(C[j][i]-r[j]))*(double)(d[i+Cf[j]-minFrame]-(C[j][i]-r[j])));
                            J+=(double)((double)(d[Cf[j][i]-minFrame]-(C[j][i]-r[j]))*(double)(d[Cf[j][i]-minFrame]-(C[j][i]-r[j])));
                        }
                    }
                }
            }
            
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}
