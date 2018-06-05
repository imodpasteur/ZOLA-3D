/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;

/**
 *
 * @author benoit
 */


import org.pasteur.imagej.utils.Matrixe;
import org.pasteur.imagej.utils.PolynomialFit;

import ij.IJ;


/** 
 *
 * @author benoit
 * One LocalizationMany per thread
 * LocalizationMany could potentially deal with multi cam and multi frame
 */
public class Localization {
    
    
    double h=.0005;//.5 nm ; 
    
    double epsilon=.0001; // distance in micrometer (stop criterion)
    
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
    
    
    
    int size;
    int id;
    int width;
    int height;
    int nbParam;//camera number
    DataPhase [] dparam;//multi dim to allow dual camera fit
    int iterMax;
    double minZ;
    double maxZ;
    
    double maxShift=50;//At one step, the maximum gradient shift is 50nm only
    
    double [][][][] subwindow;// 4 dimensions in the vect: nbParam;size;width;height
    
    double [][][][] subwindowSCMOS;// 4 dimensions in the vect: nbParam;size;width;height
    int totalsize;
    
    double photonThreshold=0;
    
    double [][][][] fisherMatrix;
    double [][] f1;
    double [][][] f0;
    double [][][] f2;
    
    double [][][][] modelPSF;
    
    double [] resCRLB;
    
    double [][] thetaA;
    double [] thetaB;
    
    double thetaX;
    double thetaY;
    double thetaZ;
    
    double likelihood;
    boolean isSCMOS=false;
    //size is the number of successive frames that are considered for the same fit
    //dparam.length is the number of cameras that are considered for the same fit
    public Localization(int size,DataPhase dparam,int iterMax,double minZ, double maxZ,int id){
        this.size=size;
        nbParam=1;
        subwindow=new double[1][size][][];
        subwindowSCMOS=new double[1][size][][];
        this.id=id;
        
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.iterMax=iterMax;
        
        this.dparam = new DataPhase[1];
        this.dparam[0] = dparam;
        
        this.width=dparam.param.sizeoutput;
        this.height=dparam.param.sizeoutput;
        
        modelPSF=new double [nbParam][size][width][height];
        fisherMatrix=new double [nbParam][size][5][5];
        f1=new double [width][height];
        f0=new double [5][width][height];
        f2=new double [5][width][height];
        resCRLB=new double [3];
        
        
        
        maxShift=dparam.param.xystep/2.;
                
        thetaA=new double[this.nbParam][size];
        thetaB=new double[this.nbParam];
        
        
            
        
    }
    
    
    //size is the number of successive frames that are considered for the same fit
    //dparam.length is the number of cameras that are considered for the same fit
    public Localization(int size,DataPhase [] dparam,int iterMax,double minZ, double maxZ,int id){
        
        this.size=size;
        nbParam=dparam.length;
        subwindow=new double[nbParam][size][][];
        subwindowSCMOS=new double[nbParam][size][][];
        this.id=id;
        
        this.minZ=minZ;
        this.maxZ=maxZ;
        this.iterMax=iterMax;
        
        this.dparam = new DataPhase[nbParam];
        for (int i=0;i<nbParam;i++){
            this.dparam[i] = dparam[i];
        }
        
        
        this.width=dparam[0].param.sizeoutput;
        this.height=dparam[0].param.sizeoutput;
        
        modelPSF=new double [nbParam][size][width][height];
        fisherMatrix=new double [nbParam][size][5][5];
        resCRLB=new double [3];
        
        
        
        maxShift=dparam[0].param.xystep/2.;
                
        thetaA=new double[this.nbParam][size];
        thetaB=new double[this.nbParam];
        
        
            
        
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
    
    
    
    
    
    public void setSubWindow(int idCam,int idFrame,double [][] subwindow){
        
        
        this.subwindow[idCam][idFrame]=subwindow;
        
        
    }
    
    public void setSubWindowScmos(int idCam,int idFrame,double [][] scmos){
        
        isSCMOS=true;
        this.subwindowSCMOS[idCam][idFrame]=scmos;
        
        
    }
    
    
    
    
    
    
    public double [][] getPSF(int idCamera,int idFrame){
        
        return modelPSF[idCamera][idFrame];
            
    }
    
    
    
    
    public void finish(){
        
        double lik=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                
                this.computeFisher(p,s, thetaX, thetaY, thetaZ, thetaA[p][s], thetaB[p], h);
                lik+=getLikelihood(p,s,thetaX,thetaY,thetaZ,thetaA[p][s],thetaB[p]);
                
            }
        }
            
        
        this.computeCRLB();
            
            
            
            
        
        
    }
    
    
    void computeFisher(int idCam,int idFrame,double x,double y,double z,double a,double b,double hdec){
        
        dparam[idCam].psf_many.computePSF(id, x, y,dparam[idCam].param.Zfocus, z);
        double [][] f =dparam[idCam].psf_many.getPSFPointer(id);

        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                f1[i][ii]=f[i][ii]*a+b;
            }
        }
        
        dparam[idCam].psf_many.computePSF(id, x+hdec, y, dparam[idCam].param.Zfocus,z);
        double [][] x2 =dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<x2.length;i++){
            for (int ii=0;ii<x2[0].length;ii++){
                f2[0][i][ii]=x2[i][ii]*a+b;
            }
        }
        dparam[idCam].psf_many.computePSF(id, x-hdec, y, dparam[idCam].param.Zfocus,z);
        double [][] x0 =dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<x0.length;i++){
            for (int ii=0;ii<x0[0].length;ii++){
                f0[0][i][ii]=x0[i][ii]*a+b;
            }
        }
        dparam[idCam].psf_many.computePSF(id, x, y+hdec, dparam[idCam].param.Zfocus,z);
        double [][] y2 =dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<y2.length;i++){
            for (int ii=0;ii<y2[0].length;ii++){
                f2[1][i][ii]=y2[i][ii]*a+b;
            }
        }
        dparam[idCam].psf_many.computePSF(id, x, y-hdec, dparam[idCam].param.Zfocus,z);
        double [][] y0 =dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<y0.length;i++){
            for (int ii=0;ii<y0[0].length;ii++){
                f0[1][i][ii]=y0[i][ii]*a+b;
            }
        }

        dparam[idCam].psf_many.computePSF(id, x, y, dparam[idCam].param.Zfocus,z+hdec);
        double [][] z2 =dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<z2.length;i++){
            for (int ii=0;ii<z2[0].length;ii++){
                f2[2][i][ii]=z2[i][ii]*a+b;
            }
        }
        dparam[idCam].psf_many.computePSF(id,x, y, dparam[idCam].param.Zfocus,z-hdec);
        double [][] z0 =dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<z0.length;i++){
            for (int ii=0;ii<z0[0].length;ii++){
                f0[2][i][ii]=z0[i][ii]*a+b;
            }
        }
        
        for (int i=0;i<f1.length;i++){
            for (int ii=0;ii<f1[0].length;ii++){
                f0[3][i][ii]=f1[i][ii]*(a-hdec)+b;
                f2[3][i][ii]=f1[i][ii]*(a+hdec)+b;
            }
        }
        
        for (int i=0;i<f1.length;i++){
            for (int ii=0;ii<f1[0].length;ii++){
                f0[4][i][ii]=f1[i][ii]*a+(b-hdec);
                f2[4][i][ii]=f1[i][ii]*a+(b+hdec);
            }
        }

        //compute fisher
        double d1,d2;
        for (int p=0;p<5;p++){
            for (int pp=0;pp<5;pp++){
                this.fisherMatrix[idCam][idFrame][p][pp]=0;
                for (int i=0;i<f1.length;i++){
                    for (int ii=0;ii<f1[0].length;ii++){
                        d1=(f2[p][i][ii]-f0[p][i][ii])/(2*Math.abs(hdec));
                        d2=(f2[pp][i][ii]-f0[pp][i][ii])/(2*Math.abs(hdec));
                        this.fisherMatrix[idCam][idFrame][p][pp]+=(1/f1[i][ii])*d1*d2;
                    }
                }
                
            }
        }
        
        
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
    
    
    
    
    public double getCRLBX(){
        
        return resCRLB[0];
            
    }
    
    
    
    
    public double getCRLBY(){
        
        return resCRLB[1];
            
    }
    
    
    
    
    public double getCRLBZ(){
        
        return resCRLB[2];
            
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
        likelihood=0;
        for (int uu=0;uu<this.size;uu++){
            for (int u=0;u<nbParam;u++){
                likelihood+=getLikelihood(u,uu, thetaX, thetaY, thetaZ, thetaA[u][uu], thetaB[u]);
            }
        }
        
        double oldX=Double.MAX_VALUE;
        double oldY=Double.MAX_VALUE;
        double oldZ=Double.MAX_VALUE;
        double distance=Double.MAX_VALUE;
        int nbIterDone=0;
        boolean break4phCount=false;
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
                
                return false;
            }
            
            
            
            //for (int u=0;u<nbParam;u++){
            //    IJ.log(""+i+"  end sec ok cam:"+u+"  id:"+id[u]);
            //}
            
            oldX=thetaX;
            oldY=thetaY;
            oldZ=thetaZ;
            
        }
//        if (id==0){
//            ImageShow.imshow(modelPSF[0][0],"mod1");
//            ImageShow.imshow(subwindow[0][0],"1");
//        }
        
        this.finish();
//        if (id==0){
//            ImageShow.imshow(modelPSF[0][0],"mod2");
//            ImageShow.imshow(subwindow[0][0],"2");
//        }
        return !break4phCount;
    }
    
    
    
    
    
    double minimizeZ(){
        double lik2=0;
        for (int s=0;s<this.size;s++){
            for (int p=0;p<this.nbParam;p++){
                
                lik2+=getLikelihood(p,s,thetaX,thetaY,thetaZ,thetaA[0][s],thetaB[0]);
            }
            
        }
        
        return minimizeZ(lik2);
    }
        
    double minimizeZ(double prevLik){
        
        double save=thetaZ;
        double lik1=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                lik1+=getLikelihood(p,s, thetaX, thetaY, thetaZ-h, thetaA[0][s], thetaB[0]);
                
            }
        }
        
        
        double lik2=prevLik;
        double lik3=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                lik3+=getLikelihood(p,s, thetaX, thetaY, thetaZ+h, thetaA[0][s], thetaB[0]);
            
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
                
                for (int p=0;p<this.nbParam;p++){
                    lik+=getLikelihood(p,s, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
                
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
            
            for (int p=0;p<this.nbParam;p++){
                lik2+=getLikelihood(p,s, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            }
            
        }
        return minimizeX(lik2);
    }
    double minimizeX(double prevLik){
        
        double save=thetaX;
        
        double lik1=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                lik1+=getLikelihood(p,s, thetaX-h, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            }
            
        }
        
        
        
        double lik2=prevLik;
        
        double lik3=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                lik3+=getLikelihood(p,s, thetaX+h, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            
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
                
                for (int p=0;p<this.nbParam;p++){
                    lik+=getLikelihood(p,s, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
                
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
            
            for (int p=0;p<this.nbParam;p++){
                lik2+=getLikelihood(p,s, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
            
            
            }
            
        }
        return minimizeY(lik2);
    }
    double minimizeY(double prevLik){
        double save=thetaY;
        
        double lik1=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                lik1+=getLikelihood(p,s, thetaX, thetaY-h, thetaZ, thetaA[0][s], thetaB[0]);
            
            
            }
            
        }
        
        double lik2=prevLik;
        double lik3=0;
        for (int s=0;s<this.size;s++){
            
            for (int p=0;p<this.nbParam;p++){
                lik3+=getLikelihood(p,s, thetaX, thetaY+h, thetaZ, thetaA[0][s], thetaB[0]);
            
            
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
                
                for (int p=0;p<this.nbParam;p++){
                    lik+=getLikelihood(p,s, thetaX, thetaY, thetaZ, thetaA[0][s], thetaB[0]);
                
                
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
    
    
    
    
    
    
    
    double minimizeA(){
        
        
        double likelihood=0;
        
        for (int p=0;p<this.nbParam;p++){
            likelihood+=minimizeA(h,0);
        }
        
        return likelihood;
    }
    
    
    
    
    double minimizeA(double h,int idParam){
        h=h*1000;
        


        for (int i=0;i<size;i++){
            double save=thetaA[idParam][i];


            double lik1=getLikelihood(idParam,i, thetaX, thetaY, thetaZ, thetaA[idParam][i]-h, thetaB[idParam]);


            double lik2=getLikelihood(idParam,i, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);


            double lik3=getLikelihood(idParam,i, thetaX, thetaY, thetaZ, thetaA[idParam][i]+h, thetaB[idParam]);


            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)*(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }

            loop:for (double gamma=1;gamma>.02;gamma/=10){
                thetaA[idParam][i]-=gamma*grad;

                double lik=getLikelihood(idParam,i, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
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
            likGlobal+=getLikelihood(idParam,i, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
        }
        return likGlobal;
    }
    
    
    double minimizeB(){
        
        
        
        double likelihood=0;
        
        for (int p=0;p<this.nbParam;p++){
            
            likelihood+=minimizeB(h,p);
        }
        
        return likelihood;
    }
    
    
    
    
    double minimizeB(double h,int idParam){
        h=h*10;


        double save=thetaB[idParam];

        double lik1=0;
        for (int s=0;s<size;s++){
            lik1+=getLikelihood(idParam,s, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]-h);
        }


        double lik2=0;
        for (int s=0;s<size;s++){
            lik2+=getLikelihood(idParam,s, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
        }


        double lik3=0;
        for (int s=0;s<size;s++){
            lik3+=getLikelihood(idParam,s, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]+h);
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
                lik+=getLikelihood(idParam,s, thetaX, thetaY, thetaZ, thetaA[idParam][s], thetaB[idParam]);
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
            likGlobal+=getLikelihood(idParam,i, thetaX, thetaY, thetaZ, thetaA[idParam][i], thetaB[idParam]);
        }
        return likGlobal;
    }
    
    
    
    public double getLikelihood(int idCam,int idFrame,double x,double y,double z,double a,double b){
        double lik=0;
         
        if (this.withregistration){
            
            if (idCam==0){
                lik=computeLikelihood(idCam,idFrame, x, y,dparam[idCam].param.Zfocus, z, a, b);
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
                
                
                lik=computeLikelihood(idCam,idFrame, x2, y2,dparam[idCam].param.Zfocus, z2, a, b);
            }
            
            
        }
        else{
            //lik=dparam[idParam].modelMany.getLikelihood(id[idParam], x, y,dparam[idParam].param.Zfocus, z, a, b);
            lik=this.computeLikelihood(idCam,idFrame, x, y,dparam[idCam].param.Zfocus, z, a, b);
        }
        
        
        
        
        return lik;
    }
    
    
    
    
    
    
    double computeLikelihood(int idCam,int idFrame,double x,double y, double zfocus,double z, double a, double b){
        dparam[idCam].psf_many.computePSF(id, x, y, zfocus, z);
        double lik=0;
        double [][] psf=dparam[idCam].psf_many.getPSFPointer(id);
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                if (isSCMOS){
                    modelPSF[idCam][idFrame][i][ii]=psf[i][ii]*(a)+b+subwindowSCMOS[idCam][idFrame][i][ii];
                }
                else{
                    modelPSF[idCam][idFrame][i][ii]=psf[i][ii]*(a)+b;
                }
                lik+=modelPSF[idCam][idFrame][i][ii]-subwindow[idCam][idFrame][i][ii]*Math.log(modelPSF[idCam][idFrame][i][ii]);
                    
            }
        }
        
        return lik;
    }
    
}