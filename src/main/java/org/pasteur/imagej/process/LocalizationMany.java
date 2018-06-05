/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;


import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.utils.Matrixe;
import jcuda.Pointer;
import ij.IJ;
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

/** 
 *
 * @author benoit
 * One LocalizationMany per thread
 * LocalizationMany could potentially deal with multi cam and multi frame
 */
public class LocalizationMany {
    
    
    double [][][][] fisherMatrix;
    
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
    
    double maxShift=50;//At one step, the maximum gradient shift is 50nm only
    
    MinimizeThread [] mA ;
    MinimizeThread [] mB ;
    ComputeLikelihoodThread [] cl;
    
    double thetaXtmp=0;
    double thetaYtmp=0;
    double thetaZtmp=0;
    
    
    double likelihood=0;
    double thetaX=0;
    double thetaY=0;
    double thetaZ=0;
    double [][] thetaA;
    double [] thetaB;
    
    
    
    double h=.0005;//.5 nm ; do not reduce it because we use float instead of double for PSF generator
    
    double epsilon=.0001; // distance in micrometer (stop criterion)
    
    DataPhase [] dparam;
    int nbParam;//number of phase (1 per camera)
    int size;//third dimension of image : multiframes
    int width;//width of image
    int height;//height of image
    double [] subwindow;// 4 dimensions in the vect: []nbParam ; size;width;height
    int totalsize;
    int iterMax;
    //for one cam
    double photonThreshold=0;
    
    double minZ;
    double maxZ;
    
    int [] id;//vector for each camera...WARNING... each size should also correspond to a different model
    
    
    double [][] modelPSF;
    
    double [] resCRLB;
            
    //in case of multiframe fitting, the likelihood is just the sum of likelihoods according to each image...
    
    
    public LocalizationMany(int size,DataPhase dparam,int iterMax,double minZ, double maxZ){
        
        
        id = new int [1];
        
        this.minZ=minZ;
        this.maxZ=maxZ;
        nbParam=1;
        this.iterMax=iterMax;
        this.dparam = new DataPhase[1];
        this.dparam[0] = dparam;
        this.width=dparam.param.sizeoutput;
        this.height=dparam.param.sizeoutput;
        this.size=size;//z size (number of slices)
        modelPSF=new double [nbParam][width*height];
        fisherMatrix=new double [nbParam][size][5][5];
        resCRLB=new double [3];
        for (int i=0;i<nbParam;i++){
            id[i]=0;
        }
        
        maxShift=dparam.param.xystep/2.;
                
        thetaA=new double[this.nbParam][this.size];
        thetaB=new double[this.nbParam];
        
        
            
        
    }
    
    
    
    // for multiple cam
    public LocalizationMany(int size,DataPhase [] dparam, int iterMax,double minZ, double maxZ){
        
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.iterMax=iterMax;
        this.dparam = dparam;
        nbParam=dparam.length;
        this.width=dparam[0].param.sizeoutput;//width & height are supposed to be the same for each source
        this.height=dparam[0].param.sizeoutput;
        maxShift=dparam[0].param.xystep/2.;
        this.size=size;//z size (number of slices)
        
        id = new int [nbParam];
        resCRLB=new double [3];
        for (int i=0;i<nbParam;i++){
            id[i]=0;
        }
        fisherMatrix=new double [nbParam][size][5][5];
        modelPSF=new double [nbParam][width*height];
        thetaA=new double[this.nbParam][this.size];
        thetaB=new double[this.nbParam];
        if (nbParam>1){
            mA = new MinimizeThread[nbParam-1];
            mB = new MinimizeThread[nbParam-1];
            cl = new ComputeLikelihoodThread[nbParam-1];
        }
    }
    
    
    
    
    
    
    
    
    
    public void setMinimumPhotonNumber(double value){
        this.photonThreshold=value;
        
        
        
    }
    
    
    public double [][] getA(){
        return thetaA;
    }
    public double [] getB(){
        return thetaB;
    }
    public double getX(){
        return thetaX;
    }
    public double getY(){
        return thetaY;
    }
    public double getZ(){
        return (thetaZ);
    }
    
    
    
    public double getGlobalLikelihood(){
        return (likelihood);
    }
    
    
    
    
    
    public void init(double a, double b, double x, double y, double z){
        
        for (int i=0;i<this.nbParam;i++){
            for (int j=0;j<this.size;j++){
                this.thetaA[i][j]=a;
            }
            this.thetaB[i]=b;
        }
        this.thetaX=x;
        this.thetaY=y;
        this.thetaZ=z;
        
        
    }
    
    public void init(double [][] a, double [] b, double x, double y, double z){
        for (int i=0;i<this.nbParam;i++){
            for (int j=0;j<this.size;j++){
                this.thetaA[i][j]=a[i][j];
            }
            this.thetaB[i]=b[i];
        }
        this.thetaX=x;
        this.thetaY=y;
        this.thetaZ=z;
    }
    
    
    
    public void setTmpCam2(double x, double y, double z){
        
        this.thetaXtmp=x;
        this.thetaYtmp=y;
        this.thetaZtmp=z;
        
        //IJ.log("init "+a+"  "+b+"  "+x+"  "+y+"  "+z+"  ");
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
        
        id[0]=dparam[0].modelMany.setSubWindow(subwindow);
        
        
    }
    
    
    
    
    
    public void setSubWindow(double [][] subwindow){
        
        for (int u=0;u<subwindow.length;u++){
            //IJ.log("subWin "+id.length+"  "+subwindow.length+"  "+dparam.length+"  "+dparam[u].modelMany);
            id[u]=dparam[u].modelMany.setSubWindow(subwindow[u]);
        }
        
        
    }
    
    
    
    
    public void setSubWindowScmos(double [][] subwindow,double [][] subwindowSCMOS){
        
        for (int u=0;u<subwindow.length;u++){
            //IJ.log("subWin "+id.length+"  "+subwindow.length+"  "+dparam.length+"  "+dparam[u].modelMany);
            id[u]=dparam[u].modelMany.setSubWindowScmos(subwindow[u],subwindowSCMOS[u]);
        }
        
        
    }
    
    
    public void finish(){
        
        double lik=0;
        for (int s=0;s<this.size;s++){
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p,thetaX,thetaY,thetaZ,thetaA[p][s],thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik+=getLikelihood(0,thetaX,thetaY,thetaZ,thetaA[0][s],thetaB[0]);
            dparam[0].modelMany.getModel(id[0],modelPSF[0]);
            this.computeFisher(0,s, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0], h);
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik+=cl[p-1].getGlobalLikelihood();
                dparam[p].modelMany.getModel(id[p],modelPSF[p]);
                this.computeFisher(p,s, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p], h);
            }
        }
            
        this.computeCRLB();

        dparam[0].modelMany.freePosit(id[0]);
        for (int p=1;p<this.nbParam;p++){
            dparam[p].modelMany.freePosit(id[p]);
            cl[p-1]=null;
        }
            
        
        
    }
    
    
    
    
    public void finishWrongPSF(){
        
        for (int s=0;s<this.size;s++){
            
            
            dparam[0].modelMany.freePosit(id[0]);
            for (int p=1;p<this.nbParam;p++){
                dparam[p].modelMany.freePosit(id[p]);
            }
            
        }
        
    }
    
    
    public double [] mytest(){
        
        this.getLikelihood(0, thetaX, thetaY, thetaZ, thetaA[0][0], thetaB[0]);
        dparam[0].modelMany.getModel(id[0],modelPSF[0]);
        return modelPSF[0];
    }
    
    
    
    
    public void finishold(){
        for (int i=0;i<this.nbParam;i++){
            int maxPhotonid=0;
            double meanBackg=0;
            for (int u=0;u<thetaA[i].length;u++){
                if (thetaA[i][u] > thetaA[i][maxPhotonid]){
                    maxPhotonid=u;
                }
                meanBackg+=thetaB[u];
            }
            meanBackg/=(double)thetaA[i].length;
            
            
            getLikelihood(i, thetaX, thetaY, thetaZ, thetaA[i][maxPhotonid], meanBackg);
            dparam[i].modelMany.getModel(id[i],modelPSF[i]);
            
            this.computeFisher(i,0, thetaX, thetaY, thetaZ, thetaA[i][maxPhotonid], meanBackg, h);
            
            dparam[i].modelMany.freePosit(id[i]);
        }
    }
    
    
    /*
    public void finishtest(){
        for (int i=0;i<this.nbParam;i++){
            
            double maxPhoton=0;
            double meanBackg=0;
            for (int u=0;u<thetaA[i].length;u++){
                if (thetaA[i][u]> maxPhoton){
                    maxPhoton=thetaA[i][u];
                }
                meanBackg+=thetaB[u];
            }
            meanBackg/=(double)thetaA[i].length;
            
            dparam[i].modelMany.getLikelihood(id[i], thetaX, thetaY, thetaZ, maxPhoton, meanBackg);
            
            dparam[i].modelMany.getModel(id[i],modelPSF[i]);
            
            dparam[i].modelMany.getCRLB(id[i], thetaX, thetaY, thetaZ, maxPhoton, meanBackg, h,resCRLB[i]);
            
            dparam[i].modelMany.freePosit(id[i]);
            
        }
    }*/
    
    
    
    
    
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
        for (int i=0;i<this.nbParam;i++){
            for (int j=0;j<this.size;j++){
                this.thetaA[i][j]=Math.max(1, this.thetaA[i][j]);
            }
            this.thetaB[i]=Math.max(.01, thetaB[i]);
        }
    }
    
    
    public boolean localize(){
        
                
        testPhotonNumber();
        
        for (int uu=0;uu<this.size;uu++){
            for (int u=0;u<nbParam;u++){
                likelihood+=getLikelihood(u, thetaX, thetaY, thetaZ, thetaA[u][uu], thetaB[u]);
            }
        }
        
        
        double oldX=Double.MAX_VALUE;
        double oldY=Double.MAX_VALUE;
        double oldZ=Double.MAX_VALUE;
        double distance=Double.MAX_VALUE;
        int nbIterDone=0;
        
        //IJ.log("theta A init "+thetaA[0][0]);
        //improve initialization of photon number
        
        likelihood=this.minimizeA();
        testPhotonNumber();
        likelihood=this.minimizeB();
        testPhotonNumber();
        
        
        
        //IJ.log("theta A first "+thetaA[0][0]);
        
        //for (int u=0;u<nbParam;u++){
        //    IJ.log("start loop    cam:"+u+"  id:"+id[u]);
        //}
        
        boolean break4phCount=false;
        loop:for (int i=0;i<iterMax;i++){
            
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("begin loop ok cam:"+u+"  id:"+id[u]);
            //}
            
            
            
            likelihood=this.minimizeZ(likelihood); 
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("min Z ok cam:"+u+"  id:"+id[u]);
            //}
            
            if ((thetaZ<minZ)){
                thetaZ=minZ;
            }
            if ((thetaZ>maxZ)){
                thetaZ=maxZ;
            }
            
            
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("comp thetaZ ok cam:"+u+"  id:"+id[u]);
            //}
            
            likelihood=this.minimizeX(likelihood);
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("min X ok cam:"+u+"  id:"+id[u]);
            //}
            
            
            likelihood=this.minimizeY(likelihood);
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("min Y ok cam:"+u+"  id:"+id[u]);
            //}
            
            //IJ.log("likelihood: "+i+"  "+likelihood+"    X:"+thetaX+"    Y:"+thetaY+"    Z:"+(thetaZ)+"    A:"+thetaA[0][0]+"    B:"+thetaB[0]);
            //try{Thread.sleep(80);}catch(Exception uut){}
            likelihood=this.minimizeA();
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("min A ok cam:"+u+"  id:"+id[u]);
            //}
            
            //thetaA[0][0]=3300;
            //IJ.log("likelihood: "+i+"  "+likelihood+"    X:"+thetaX+"    Y:"+thetaY+"    Z:"+(thetaZ)+"    A:"+thetaA[0][0]+"    B:"+thetaB[0]);
            //try{Thread.sleep(90);}catch(Exception uut){}
            testPhotonNumber();
            
            //IJ.log("test photon numb ok ");
            
            likelihood=this.minimizeB();
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("min B ok cam:"+u+"  id:"+id[u]);
            //}
            
            testPhotonNumber();
            //IJ.log("likelihood: "+i+"  "+likelihood+"    X:"+thetaX+"    Y:"+thetaY+"    Z:"+(thetaZ)+"    A:"+thetaA[0][0]+"    B:"+thetaB[0]);
            
            
            
            
            distance=Math.sqrt((this.thetaX-oldX)*(this.thetaX-oldX)+(this.thetaY-oldY)*(this.thetaY-oldY)+(this.thetaZ-oldZ)*(this.thetaZ-oldZ));

            
            boolean breaking=true;
            if (i>4){
                for (int a=0;a<thetaA.length;a++){
                    for (int aa=0;aa<thetaA[0].length;aa++){
                        if ((thetaA[a][aa]>photonThreshold)){//si tous ont un faible nombre de photons -> break
                            breaking=false;
                        }
                    }
                }
                if (breaking){
                    break4phCount=true;
                    IJ.log("breaked");
                }
            }
            else{
                breaking=false;
            }
        
            if ((distance<epsilon)||(breaking)){
                nbIterDone=i;
                break loop;
            }
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log("end ok cam:"+u+"  id:"+id[u]);
            //}
            
            
            
            if (Double.isNaN(likelihood)){
                this.finishWrongPSF();
                return false;
            }
            
            
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log(""+i+"  end sec ok cam:"+u+"  id:"+id[u]);
            //}
            
            oldX=thetaX;
            oldY=thetaY;
            oldZ=thetaZ;
            
        }
        
        
        this.finish();
        
        
        return !break4phCount;
    }
    
    
    
    
    
    
    
    
    
    
    public double getLikelihood(int idParam,double x,double y,double z,double a,double b){
        double lik=0;
         
        if (this.withregistration){
            
            if (idParam==0){
                lik=dparam[idParam].modelMany.getLikelihood(id[idParam], x, y,dparam[idParam].param.Zfocus, z, a, b);
            }
            else{
                double x2, y2, z2;
                
                //first, reverse x & y and apply the drift correction of cam 1 (necessary for registration)
                x2=(-(x)-this.dx1);
                y2=(-(y)-this.dy1);
                z2=z-this.dz1;
                
                //registration: cam1 -> cam2
                double [] v = new double [3];
                v[0]=positionRefX+x2;
                v[1]=positionRefY+y2;
                v[2]=z2;
                double [] res=pf.transform(v);
                x2=res[0]-positionRefX;
                y2=res[1]-positionRefY;
                z2=res[2];
                
                
                double xs=x2+this.dx2;
                
                
                //remove drift correction of second cam
                x2=x2+this.dx2;
                y2=y2+this.dy2;
                z2=z2+this.dz2;
                
                
                
                //finally, shift according to the window position and reverse x & y
                x2=-(x2-this.decX);
                y2=-(y2-this.decY);
                
                lik=dparam[idParam].modelMany.getLikelihood(id[idParam], x2, y2,dparam[idParam].param.Zfocus, z2, a, b);
            }
            
            
        }
        else{
            lik=dparam[idParam].modelMany.getLikelihood(id[idParam], x, y,dparam[idParam].param.Zfocus, z, a, b);
        }
        
        
        
        
        return lik;
    }
    
    public void computeCRLB(){
        
        double [][] mat = new double [5][5];
        for (int i=0;i<5;i++){
            for (int ii=0;ii<5;ii++){
                mat[i][ii]=0;
                for (int p=0;p<nbParam;p++){
                    for (int s=0;s<size;s++){
                        mat[i][ii]+=fisherMatrix[p][s][i][ii];
                    }
                }
            }
        }
        //inverse matrix
        Matrixe mm = new Matrixe(mat);
        Matrixe mmres = new Matrixe(mm.getNrows(), mm.getNcols());
        try{
            Matrixe.inverse(mm,mmres);
            mmres.getMatrixe(mat);
            for (int i=0;i<resCRLB.length;i++){
                resCRLB[i]=Math.sqrt(mat[i][i]);
            }
        }catch(Exception ee){

            //dont take into account covar if non inversible
            //IJ.log("fisher matrix non inversible at z="+z);
            for (int i=0;i<resCRLB.length;i++){
                if (mat[i][i]!=0){
                    resCRLB[i]=Math.sqrt(1/mat[i][i]);
                }
                else{
                    resCRLB[i]=Double.MAX_VALUE;
                }
            }

        }
        mmres.free();
        mm.free();
        mat=null;
            
        
    }
     
    public void computeFisher(int idParam,int frame,double x,double y,double z,double a,double b,double h){
        
        
        
        dparam[idParam].modelMany.computeFisherMatrix(id[idParam], x, y, dparam[idParam].param.Zfocus,z, a, b, h,this.fisherMatrix[idParam][frame]);
        
        
        
        
    }
    
    
    double minimizeB(){
        
        
        
        for (int p=1;p<this.nbParam;p++){
            mB[p-1]=new MinimizeThread(h,p);
            mB[p-1].minimizeB();
            mB[p-1].start();
        }
        
        double likelihood=minimizeB(h,0);
        
        for (int p=1;p<this.nbParam;p++){
            try{
                mB[p-1].join();
            }
            catch(Exception e){
                IJ.log("join LocalizationMany impossible minimizationB");
            }
            likelihood+=mB[p-1].getGlobalLikelihood();
            mB[p-1]=null;
        }
        
        return likelihood;
    }
    
    
    
    double minimizeBold(){
        double h=this.h*10;
        for (int p=0;p<this.nbParam;p++){
            double save=thetaB[p];
            
            
            double lik1=0;
            for (int s=0;s<size;s++){
                lik1+=getLikelihood(p,thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]-h);
            }
            
            double lik2=0;
            for (int s=0;s<size;s++){
                lik2+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
            
            double lik3=0;
            for (int s=0;s<size;s++){
                lik3+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]+h);
            }
                    
                
            
            
            this.thetaB[p]=save;
            
            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)*(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }
            loop:for (double gamma=1;gamma>.02;gamma/=10){
                thetaB[p]-=gamma*grad;
                
                double lik=0;
                for (int s=0;s<size;s++){
                    lik+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
            
                if (lik<lik2){
                    break loop;
                }
                else{
                    thetaB[p]=save;
                }
            }
        }
        double lik=0;
        for (int p=0;p<this.nbParam;p++){
            for (int i=0;i<this.size;i++){
                lik+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][i], thetaB[p]);
            }
        }
        return lik;
    }
    
    
    
    double minimizeA(){
        
        
        for (int p=1;p<this.nbParam;p++){
            mA[p-1]=new MinimizeThread(h,p);
            mA[p-1].minimizeA();
            mA[p-1].start();
        }
        
        double likelihood=minimizeA(h,0);
        
        for (int p=1;p<this.nbParam;p++){
            try{
                mA[p-1].join();
            }
            catch(Exception e){
                IJ.log("join LocalizationMany impossible minimizationA");
            }
            likelihood+=mA[p-1].getGlobalLikelihood();
            mA[p-1]=null;
        }
        
        return likelihood;
    }
    
    
    
    double minimizeAold(){
        double h=this.h*1000;
        
        for (int p=0;p<this.nbParam;p++){
            
            
            for (int i=0;i<this.size;i++){
                double save=thetaA[p][i];
                
                
                double lik1=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][i]-h, thetaB[p]);
                
                
                double lik2=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][i], thetaB[p]);
                
                
                double lik3=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][i]+h, thetaB[p]);
                

                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)*(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    thetaA[p][i]-=gamma*grad;

                    double lik=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][i], thetaB[p]);
                    //IJ.log("grad "+gamma+"  "+grad+"  "+thetaA[p][i]+"  "+lik);
                    if (lik<lik2){
                        break loop;
                    }
                    else{
                        thetaA[p][i]=save;
                    }
                }
            }
        }
        
        double lik=0;
        for (int p=0;p<this.nbParam;p++){
            for (int i=0;i<this.size;i++){
                lik+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][i], thetaB[p]);
            }
        }
        return lik;
    }
    
    
    
    
    
    
    double minimizeZ(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p,thetaX,thetaY,thetaZ,thetaA[p][s],thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik2+=getLikelihood(0,thetaX,thetaY,thetaZ,thetaA[0][s],thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik2+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
        }
        return minimizeZ(lik2);
    }
    double minimizeZ(double prevLik){
        
        double save=thetaZ;
        double lik1=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ-h, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik1+=getLikelihood(0, thetaX, thetaY, thetaZ-h, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik1+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
        }
        
        
        
        double lik2=prevLik;
        
        double lik3=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ+h, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik3+=getLikelihood(0, thetaX, thetaY, thetaZ+h, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik3+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
        }
        
        
        
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        if (grad>maxShift){
            grad=maxShift;
        }
        else if (grad<-maxShift){
            grad=-maxShift;
        }
        
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            
            
            thetaZ-=gamma*grad;
            
            double lik=0;
            for (int s=0;s<this.size;s++){
                
                for (int p=1;p<this.nbParam;p++){
                    cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
                for (int p=1;p<this.nbParam;p++){
                    cl[p-1].start();
                }
                
                lik+=getLikelihood(0, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
                
                for (int p=1;p<this.nbParam;p++){
                    try{
                        cl[p-1].join();
                    }
                    catch(Exception e){
                        IJ.log("join LocalizationMany impossible minimizeZ");
                    }
                    lik+=cl[p-1].getGlobalLikelihood();
                    cl[p-1]=null;
                }
            
            }
            
            if (lik<lik2){
                prevLik=lik;
                break loop;
            }
            else{
                thetaZ=save;
            }
        }
        
        return prevLik;
    }
    
    
    
    
    
    double minimizeZold(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                
            
                lik2+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        return minimizeZ(lik2);
    }
    double minimizeZold(double prevLik){
        
        double save=thetaZ;
        double lik1=0;
        IJ.log("Z... "+id[0]);
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik1+=getLikelihood(p, thetaX, thetaY, thetaZ-h, thetaA[p][s], thetaB[p]);
                IJ.log("Z... lik 1 ... "+p+"   "+id[0]);
            }
        }
        
        IJ.log("Z... lik1 ok "+id[0]);
        
        double lik2=prevLik;
        
        double lik3=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik3+=getLikelihood(p, thetaX, thetaY, thetaZ+h, thetaA[p][s], thetaB[p]);
            }
        }
        
        IJ.log("Z... lik3 ok "+id[0]);
        
        
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            thetaZ-=gamma*grad;
            
            double lik=0;
            for (int s=0;s<this.size;s++){
                for (int p=0;p<this.nbParam;p++){
                    lik+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
            }
            
            if (lik<lik2){
                prevLik=lik;
                break loop;
            }
            else{
                thetaZ=save;
            }
        }
        return prevLik;
    }
    
    
    
    
    
    
    
    
    double minimizeX(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik2+=getLikelihood(0, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik2+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
            
        }
        return minimizeX(lik2);
    }
    double minimizeX(double prevLik){
        
        double save=thetaX;
        
        double lik1=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX-h, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik1+=getLikelihood(0, thetaX-h, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik1+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
            
        }
        
        
        
        double lik2=prevLik;
        
        double lik3=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX+h, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik3+=getLikelihood(0, thetaX+h, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeX");
                }
                lik3+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
            
        }
        
        
        
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        if (grad>maxShift){
            grad=maxShift;
        }
        else if (grad<-maxShift){
            grad=-maxShift;
        }
        
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            
            thetaX-=gamma*grad;
            
            
            double lik=0;
            for (int s=0;s<this.size;s++){
                
                for (int p=1;p<this.nbParam;p++){
                    cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
                for (int p=1;p<this.nbParam;p++){
                    cl[p-1].start();
                }
                
                lik+=getLikelihood(0, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
                
                for (int p=1;p<this.nbParam;p++){
                    try{
                        cl[p-1].join();
                    }
                    catch(Exception e){
                        IJ.log("join LocalizationMany impossible minimizeZ");
                    }
                    lik+=cl[p-1].getGlobalLikelihood();
                    cl[p-1]=null;
                }
                
            }
            if (lik<lik2){
                prevLik=lik;
                break loop;
            }
            else{
                
                thetaX=save;
            }
        }
        return prevLik;
    }
    
    
    
    
    
    
    
    double minimizeXold(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik2+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        return minimizeX(lik2);
    }
    double minimizeXold(double prevLik){
        
        double save=thetaX;
        
        double lik1=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik1+=getLikelihood(p, thetaX-h, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        
        
        
        double lik2=prevLik;
        
        double lik3=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik3+=getLikelihood(p, thetaX+h, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        
        
        
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            thetaX-=gamma*grad;

            double lik=0;
            for (int s=0;s<this.size;s++){
                for (int p=0;p<this.nbParam;p++){
                    lik+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
            }
            if (lik<lik2){
                prevLik=lik;
                break loop;
            }
            else{
                thetaX=save;
            }
        }
        return prevLik;
    }
    
    
    
    
    
    double minimizeY(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik2+=getLikelihood(0, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik2+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
            
        }
        return minimizeY(lik2);
    }
    double minimizeY(double prevLik){
        double save=thetaY;
        
        double lik1=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY-h, thetaZ, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik1+=getLikelihood(0, thetaX, thetaY-h, thetaZ, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik1+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
            
        }
        
        double lik2=prevLik;
        double lik3=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=1;p<this.nbParam;p++){
                cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY+h, thetaZ, thetaA[p][s], thetaB[p]);
            }
            for (int p=1;p<this.nbParam;p++){
                cl[p-1].start();
            }
            
            lik3+=getLikelihood(0, thetaX, thetaY+h, thetaZ, thetaA[0][s], thetaB[0]);
            
            for (int p=1;p<this.nbParam;p++){
                try{
                    cl[p-1].join();
                }
                catch(Exception e){
                    IJ.log("join LocalizationMany impossible minimizeZ");
                }
                lik3+=cl[p-1].getGlobalLikelihood();
                cl[p-1]=null;
            }
            
        }
        
        
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        if (grad>maxShift){
            grad=maxShift;
        }
        else if (grad<-maxShift){
            grad=-maxShift;
        }
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            
            thetaY-=gamma*grad;
            
            
            double lik=0;
            for (int s=0;s<this.size;s++){
                
                for (int p=1;p<this.nbParam;p++){
                    cl[p-1]=new ComputeLikelihoodThread(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
                for (int p=1;p<this.nbParam;p++){
                    cl[p-1].start();
                }
                
                lik+=getLikelihood(0, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
                
                for (int p=1;p<this.nbParam;p++){
                    try{
                        cl[p-1].join();
                    }
                    catch(Exception e){
                        IJ.log("join LocalizationMany impossible minimizeZ");
                    }
                    lik+=cl[p-1].getGlobalLikelihood();
                    cl[p-1]=null;
                }
                
            }
            if (lik<lik2){
                prevLik=lik;
                break loop;
            }
            else{
                thetaY=save;
            }
        }
        return prevLik;
    }
    
    
    
    
    
    
    
    double minimizeYold(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik2+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        return minimizeY(lik2);
    }
    double minimizeYold(double prevLik){
        double save=thetaY;
        
        double lik1=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik1+=getLikelihood(p, thetaX, thetaY-h, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        
        double lik2=prevLik;
        double lik3=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                lik3+=getLikelihood(p, thetaX, thetaY+h, thetaZ, thetaA[p][s], thetaB[p]);
            }
        }
        
        
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            thetaY-=gamma*grad;

            double lik=0;
            for (int s=0;s<this.size;s++){
                for (int p=0;p<this.nbParam;p++){
                    lik+=getLikelihood(p, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p]);
                }
            }
            if (lik<lik2){
                prevLik=lik;
                break loop;
            }
            else{
                thetaY=save;
            }
        }
        return prevLik;
    }
    
    
    
    
    double minimizeA(double h,int idParam){
        h=h*1000;
        


        for (int i=0;i<size;i++){
            double save=thetaA[idParam][i];


            double lik1=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i]-h, thetaB[idParam]);


            double lik2=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);


            double lik3=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i]+h, thetaB[idParam]);


            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)*(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }

            loop:for (double gamma=1;gamma>.02;gamma/=10){
                thetaA[idParam][i]-=gamma*grad;

                double lik=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
                //IJ.log("grad "+gamma+"  "+grad+"  "+thetaA[p][i]+"  "+lik);
                if (lik<lik2){
                    lik2=lik;
                    save=thetaA[idParam][i];
                    //I can not break here because getLikelihood has to be called the same numb of time for each camera (synchro problem otherwise)
                    //best stategy: I dont break and I change thetaB[idParam] at each loop if condition is accepted
                }
                else{
                    thetaA[idParam][i]=save;
                }
            }
        }

        double likGlobal=0;
        for (int i=0;i<size;i++){
            likGlobal+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
        }
        return likGlobal;
    }
    
    
    
    
    double minimizeB(double h,int idParam){
        h=h*10;


        double save=thetaB[idParam];

        double lik1=0;
        for (int s=0;s<size;s++){
            lik1+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]-h);
        }


        double lik2=0;
        for (int s=0;s<size;s++){
            lik2+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
        }


        double lik3=0;
        for (int s=0;s<size;s++){
            lik3+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]+h);
        }



        thetaB[idParam]=save;

        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)*(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            thetaB[idParam]-=gamma*grad;
            
            double lik=0;
            for (int s=0;s<size;s++){
                lik+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
            }


            if (lik<lik2){
                lik2=lik;
                save=thetaB[idParam];
                //I can not break here because getLikelihood has to be called the same numb of time for each camera (synchro problem otherwise)
                //best stategy: I dont break and I change thetaB[idParam] at each loop if condition is accepted
            }
            else{
                thetaB[idParam]=save;
            }

        }


        double likGlobal=0;
        for (int i=0;i<size;i++){
            likGlobal+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
        }
        return likGlobal;
    }
    
    //necessary to launch each cam in parallel
    //obligatoire car la minimization de A et B pour chaque camera n'est pas appelé au meme moment et on a des pb de synchro sans thread
    public class MinimizeThread extends Thread{
        
        int launcher=0;
        
        int idParam=0;
        double h;
        
        double likGlobal=0;
        
        MinimizeThread(double h,int idParam){
            this.h=h;
            this.idParam=idParam;
            
        }
        
        
        
        
        void minimizeA(){
            
            launcher=0;
        }
        
        void minimizeB(){
            
            launcher=1;
        }
        
        
        public double getGlobalLikelihood(){
            return likGlobal;
        }
        
        
        public void run(){
            
            /*
            if (launcher==0){
                double h=this.h*1000;
                //IJ.log("this get lik bef "+this.getLikelihood());



                for (int i=0;i<size;i++){
                    double save=thetaA[idParam][i];


                    double lik1=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i]-h, thetaB[idParam]);


                    double lik2=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);


                    double lik3=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i]+h, thetaB[idParam]);


                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)*(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                    }

                    double saveprec=save;
                    loop:for (double gamma=1;gamma>.02;gamma/=10){
                        thetaA[idParam][i]-=gamma*grad;

                        double lik=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
                        //IJ.log("grad "+gamma+"  "+grad+"  "+thetaA[p][i]+"  "+lik);
                        if (lik<lik2){
                            lik2=lik;
                            saveprec=thetaA[idParam][i];
                            //I can not break here because getLikelihood has to be called the same numb of time for each camera (synchro problem otherwise)
                            //best stategy: I dont break and I change thetaB[idParam] at each loop if condition is accepted
                        }
                        else{
                            thetaA[idParam][i]=saveprec;
                        }
                    }
                }

                likGlobal=0;
                for (int i=0;i<size;i++){
                    likGlobal+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
                }
                
                
            }*/
            
            
            if (launcher==0){
                double h=this.h*10;
                
                
                

                
                for (int s=0;s<size;s++){
                    
                    double save=thetaA[idParam][s];
                    double lik1=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s]-h, thetaB[idParam]);
                    double lik2=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
                    double lik3=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s]+h, thetaB[idParam]);
                    
                    thetaA[idParam][s]=save;
                    
                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)*(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                    }
                    double saveprec=save;
                    
                    
                    loop:for (double gamma=1;gamma>.02;gamma/=10){
                        thetaA[idParam][s]-=gamma*grad;

                        double lik=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
                        


                        if (lik<lik2){
                            saveprec=thetaA[idParam][s];
                            lik2=lik;
                            //I can not break here because getLikelihood has to be called the same numb of time for each camera (synchro problem otherwise)
                            //best stategy: I dont break and I change thetaB[idParam] at each loop if condition is accepted
                        }
                        else{
                            thetaA[idParam][s]=saveprec;
                        }

                    }
                }
                




                
            }
            
            else if (launcher==1){
                double h=this.h*10;
                
                
                double save=thetaB[idParam];

                double lik1=0;
                for (int s=0;s<size;s++){
                    lik1+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]-h);
                }
                

                double lik2=0;
                for (int s=0;s<size;s++){
                    lik2+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
                }
                

                double lik3=0;
                for (int s=0;s<size;s++){
                    lik3+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]+h);
                }



                thetaB[idParam]=save;

                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)*(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                double saveprec=save;
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    thetaB[idParam]-=gamma*grad;
                    
                    double lik=0;
                    for (int s=0;s<size;s++){
                        lik+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
                    }
                    
                    
                    if (lik<lik2){
                        saveprec=thetaB[idParam];
                        lik2=lik;
                        //I can not break here because getLikelihood has to be called the same numb of time for each camera (synchro problem otherwise)
                        //best stategy: I dont break and I change thetaB[idParam] at each loop if condition is accepted
                    }
                    else{
                        thetaB[idParam]=saveprec;
                    }
                    
                }
            }
            
            likGlobal=0;
            for (int i=0;i<size;i++){
                likGlobal+=getLikelihood(idParam, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
            }
        }
        
    }
        
    
    
    
    
    
    
        
        
        
    //necessary to launch each cam in parallel
    //obligatoire car la minimization de X, Y, et Z pour chaque camera n'est pas appelé au meme moment et on a des pb de synchro sans thread
    public class ComputeLikelihoodThread extends Thread{
        
        
        
        int id;
        double a,b,x,y,z;
        
        double likGlobal;
        boolean executed;
        
        ComputeLikelihoodThread(int idParam,double x, double y,double z, double a,double b){
            
            this.id=idParam;
            this.x=x;
            this.y=y;
            this.z=z;
            this.a=a;
            this.b=b;
            executed=false;
            
        }
        
        
        public double getGlobalLikelihood(){
            if (executed){
                return likGlobal;
            }
            else{
                IJ.log("oops: problem in LocalizationMany, you should call join before getting likelihood");
                return 0;
            }
        }
        
        
        
        
        public void run(){
            

            likGlobal=getLikelihood(id, x, y, z, a, b);
            executed=true;

        }
            
            
        
        
        
        
    }
    
    
    
    
    
}
