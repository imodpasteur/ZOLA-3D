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
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.utils.Zernike;
import ij.IJ;


/**
 *
 * @author benoit
 */
public class ZernikePhase {
    
    int nbMax=5;
    public int computeAll=nbMax+1;//start by computing everything
    int [] posit = new int[nbMax];//at max nbMax changes without recomputing everything
    double [] diff = new double[nbMax];//at max nbMax changes without recomputing everything
    
    double [] ones;
    
    
    
    double [] weightedZernikePhase;
    
    double [] Zvect;
    public int [] coef;
    double [] a;
    //double [][] A;
    double [] tmp;
    public int nbDataPerImage;
    PhaseParameters param;
    int m;
    public int numCoef;
    int incr;
    int lda;
    int numberOfCoef;
    public int [] complexity;//complexity of zernike polynomial
    boolean multiTrainingWithParabola=false;
    
    public ZernikePhase(PhaseParameters param,int zernikePolyNumber,int method){
        this.multiTrainingWithParabola=false;
        this.param=param;
        Zernike z = new Zernike(param.size_cpu,param.sizeRadiusRingPixel_cpu,zernikePolyNumber);
        //ImageShow.imshow(z.Z,"Z");
        
        numberOfCoef=zernikePolyNumber;
        
        
        
        if (method==1){
            numberOfCoef=0;
            for (int i=0;i<zernikePolyNumber;i++){
                if (z.poly[i][1]==0){
                    numberOfCoef++;
                }
            }
        }
        else if (method==2){
            numberOfCoef=0;
            for (int i=0;i<zernikePolyNumber;i++){
                if (z.poly[i][1]%2==0){
                    numberOfCoef++;
                }
            }
        }
        
        nbDataPerImage=param.sizeDisk_cpu;
        weightedZernikePhase=new double[nbDataPerImage];
        Zvect=new double [nbDataPerImage*numberOfCoef];
        this.a=new double [numberOfCoef];
        //this.A=new double [dim][numberOfCoef];
        this.coef=new int [numberOfCoef];
        numCoef=coef.length;
        this.complexity=new int [numberOfCoef];
        for (int rz=0,id=0;rz<zernikePolyNumber;rz++){
            if (method==1){
                if (z.poly[rz][1]==0){
//                    for (int or=0;or<dim;or++){
//                        A[or][id]=0;
//                    }
                    a[id]=0;
                    this.complexity[id]=z.complexity[rz];
                    this.coef[id]=rz;
                    for (int i=0;i<nbDataPerImage;i++){
                        Zvect[id*nbDataPerImage+i]=z.Z[rz][param.disk2D_cpu[i][0]][param.disk2D_cpu[i][1]];
                    }
                    id++;
                }
            }
            else if (method==2){
                if (z.poly[rz][1]%2==0){
//                    for (int or=0;or<dim;or++){
//                        A[or][id]=0;
//                    }
                    a[id]=0;
                    this.complexity[id]=z.complexity[rz];
                    this.coef[id]=rz;
                    for (int i=0;i<nbDataPerImage;i++){
                        Zvect[id*nbDataPerImage+i]=z.Z[rz][param.disk2D_cpu[i][0]][param.disk2D_cpu[i][1]];
                    }
                    id++;
                }
            }
            else{
//                for (int or=0;or<dim;or++){
//                    A[or][rz]=0;
//                }
                a[rz]=0;
                this.coef[rz]=rz;
                this.complexity[rz]=z.complexity[rz];
                for (int i=0;i<nbDataPerImage;i++){
                    Zvect[rz*nbDataPerImage+i]=z.Z[rz][param.disk2D_cpu[i][0]][param.disk2D_cpu[i][1]];
                }
            }
        }
        
        
        
    }
    
    
    
    
    
    public ZernikePhase(PhaseParameters param,int [] coef){
        this.multiTrainingWithParabola=false;
        int maxiCoef=0;
        for (int i=0;i<coef.length;i++){
            if (maxiCoef<coef[i]){
                maxiCoef=coef[i];
            }
        }
        numCoef=coef.length;
        this.coef=new int[coef.length];
        
        this.param=param;
        Zernike z = new Zernike(param.size_cpu,param.sizeRadiusRingPixel_cpu,maxiCoef+1);
        
        numberOfCoef=coef.length;
        
            
        
        nbDataPerImage=param.sizeDisk_cpu;
        weightedZernikePhase=new double[nbDataPerImage];
        Zvect=new double [nbDataPerImage*numberOfCoef];
        //this.A=new double [dim][numberOfCoef];
        this.a=new double [numberOfCoef];
        this.coef=new int [numberOfCoef];
        this.complexity=new int [numberOfCoef];
        for (int rz=0,id=0;id<numberOfCoef;rz++){
            if (rz==coef[id]){
//                for (int or=0;or<dim;or++){
//                    A[or][id]=0;
//                }
                a[id]=0;
                this.coef[id]=coef[id];
                this.complexity[id]=z.complexity[rz];
                for (int i=0;i<nbDataPerImage;i++){
                    Zvect[id*nbDataPerImage+i]=z.Z[rz][param.disk2D_cpu[i][0]][param.disk2D_cpu[i][1]];
                }
                id++;
            }
            
            
        }
        
        
        
    }
    
    
    
    
    public void computeCombination(){
        
        if ((computeAll>=1)&&(computeAll<=nbMax)){//to compute faster...if only 1 element changed in a (as it is often the case)
            for (int i=0;i<nbMax;i++){
                for (int ii=0;ii<this.nbDataPerImage;ii++){
                    this.weightedZernikePhase[ii]+=Zvect[posit[i]*nbDataPerImage+ii]*this.diff[i];
                }
            }
        }
        else {
            for (int ii=0;ii<this.nbDataPerImage;ii++){
                this.weightedZernikePhase[ii]=0;
            }

            for (int i=0;i<this.numberOfCoef;i++){
                for (int ii=0;ii<this.nbDataPerImage;ii++){
                    this.weightedZernikePhase[ii]+=Zvect[i*nbDataPerImage+ii]*this.a[i];
                }

            }
        }
        computeAll=0;
        
        
    }
    
    
    
    
    public void computeCombinationPlusOtherPhase(double [] otherPhase){
        computeAll=nbMax+1;
        this.computeCombination();
        for (int i=0;i<weightedZernikePhase.length;i++){
            weightedZernikePhase[i]+=otherPhase[i];
        }
        computeAll=nbMax+1;
    }
    
    
    
    
    
    
    public void setA(int posit,double value){
        
        if (computeAll<nbMax){
            boolean ok=true;
            for (int i=0;i<computeAll;i++){
                if (this.posit[i]==posit){
                    ok=false;
                    this.diff[i]+=value-a[posit];
                }
            }
            if (ok){
                this.diff[computeAll]=value-a[posit];
                this.posit[computeAll]=posit;
                computeAll++;
            }
            //IJ.log("A to modif "+this.posit[computeAll-1]+"  "+this.diff[computeAll-1]+"  "+value+"  "+a[posit]);
        }
        
        a[posit]=value;
        
    }
    
    public double getA(int posit){
        return a[posit];
    }
    
    
    public void setA(double [] a){
        for (int i=0;i<a.length;i++){
            this.a[i]=a[i];
        }
        computeAll=nbMax+1;
    }
    
    public double [] getAPointer(){
        return a;
    }
    
    public double [] getPhasePointer(){
        return this.weightedZernikePhase;
    }
    
    
    
    public void setMatAtPosit(double [] k,int posit){
        
        for (int i=0;i<nbDataPerImage;i++){
            Zvect[posit*nbDataPerImage+i]=k[i];
        }
        computeAll=nbMax+1;
        
    }
    
    
    
}

