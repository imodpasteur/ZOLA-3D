/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;
import org.pasteur.imagej.process.*;

import org.pasteur.imagej.cuda.*;

import org.pasteur.imagej.utils.ImageShow;
import ij.IJ;
import java.util.ArrayList;
import java.util.Collections;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.cublasHandle;
import jcuda.jcusparse.cusparseHandle;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import jcuda.runtime.cudaError;
import jcuda.runtime.cudaMemcpyKind;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;

/**
 *
 * @author benoit
 */
public class GenericPhase_ {
    
    
    private double modulo=100*2*Math.PI;//100* : like no modulo operator
    
    public int angleNumber=0;
    public int distNumber=0;
    
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;
    
    CUstream custream;
    
    
    int cudaResult;
    
    double [] phase;
    Pointer host_phase;
    Pointer device_phase;
    
    
    //Pointer device_id;
    //Pointer host_id;
    int [] id;
    int [][] id_table;//  first dim: idNumber   ; second dim: the corresponding positions
    double [][] phase_list;//  first dim: idNumber   ; second dim: the corresponding positions
    
    PhaseParameters param;
    int nbDataPerImage;
    
    
    Pointer  host_tmp;
    double [] tmp;
    
    
    int maxId=0;
    
    GenericPhase_(PhaseParameters param){
        
        
        this.param=param;
        
        
        
        nbDataPerImage=param.sizeDisk;
        
        id=new int [nbDataPerImage];
        //host_id=Pointer.to(id);
        
        
        
        phase=new double [nbDataPerImage]; 
        host_phase=Pointer.to(phase);
        
        device_phase=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_phase, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 8");}
        
        
        //device_id=new Pointer();
        //cudaResult = JCuda.cudaMalloc(device_id, nbDataPerImage * Sizeof.DOUBLE);
        //if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 8");}
        
        tmp=new double [1];
        host_tmp=Pointer.to(tmp);
        
    }
    
    
    
    public void free(){
        
    
        
        JCuda.cudaFree(device_phase);
        
        //JCuda.cudaFree(device_id);
        
        
    }
    
    
    
    public void showPhase(){
        
        double [][] mat = new double [param.size][param.size];
        for (int i=0;i<param.sizeDisk;i++){
            mat[param.disk2D[i][0]][param.disk2D[i][1]]=phase[i];
        }
        ImageShow.imshow(mat,"generic phase");
    }
    
    public int createSplits(int rad_number, int rho_number){
        this.angleNumber=rad_number;
        this.distNumber=rho_number;
        
        maxId=rad_number*rho_number;
        
        //double shiftAngle=3.141592/(double)rad_number;//half of angle
        double shiftAngle=Math.random()*Math.PI/(double)rad_number;//half of angle
        
        id_table= new int [maxId][];
        phase_list= new double [maxId][];
        double maxDist=param.sizeRadiusRingPixel;
        
        //int [][] mat = new int [param.size][param.size];
        double angle,dist,x,y;
        
        for (int i=0;i<nbDataPerImage;i++){
            x=param.disk2D[i][0]-param.centerFourierImage;
            y=param.disk2D[i][1]-param.centerFourierImage;
            if (x!=0){
                angle=Math.atan2(y, x)+Math.PI+shiftAngle;
                angle=angle%(2*Math.PI);
            }
            else{
                angle=0;
            }
            dist=Math.sqrt(x*x+y*y);
            
            double maxArea=Math.sqrt(maxDist)*maxDist*3.141592;
            double theArea=Math.sqrt(dist)*dist*3.141592;
            
            //id[i]=rad_number*(int)(dist/(maxDist/(double)rho_number))+(int)(angle/(2*Math.PI/rad_number));
            
            id[i]=rad_number*(int)(theArea/(maxArea/(double)rho_number))+(int)(angle/(2*Math.PI/rad_number));
            
            //mat[param.disk2D[i][0]][param.disk2D[i][1]]=id;
        }
        
        int [] nb = new int[maxId];
        for (int i=0;i<nbDataPerImage;i++){
            nb[id[i]]++;
        }
        for (int i=0;i<maxId;i++){
            if (nb[i]!=0){
                id_table[i]=new int[nb[i]];
                phase_list[i]=new double[nb[i]];
            }
            else{
                id_table[i]=null;
                phase_list[i]=null;
            }
            nb[i]=0;
        }
        for (int i=0;i<nbDataPerImage;i++){
            id_table[id[i]][nb[id[i]]]=i;
            phase_list[id[i]][nb[id[i]]]=phase[i];
            nb[id[i]]++;
        }
        
        
        return maxId;
        
    }
    
    
    
    //rad: [0:2*PI]      rho: [0:1]      radius: [0:1]
    public int createSplitsAs2Rings(double rad, double rho,double radius){
        
        maxId=3;
        
        
        
        id_table= new int [maxId][];
        phase_list= new double [maxId][];
        double maxDist=param.sizeRadiusRingPixel;
        
        //int [][] mat = new int [param.size][param.size];
        double angle,dist1,dist2,x,y;
        
        for (int i=0;i<nbDataPerImage;i++){
            
            double xp1=((rho*maxDist)*Math.cos(rad))+param.centerFourierImage;
            double yp1=((rho*maxDist)*Math.sin(rad))+param.centerFourierImage;
            
            double xp2=((rho*maxDist)*Math.cos(rad+Math.PI))+param.centerFourierImage;
            double yp2=((rho*maxDist)*Math.sin(rad+Math.PI))+param.centerFourierImage;
            
            x=param.disk2D[i][0];
            y=param.disk2D[i][1];
            
            
            dist1=Math.sqrt((x-xp1)*(x-xp1)+(y-yp1)*(y-yp1));
            dist2=Math.sqrt((x-xp2)*(x-xp2)+(y-yp2)*(y-yp2));
            
            
            if (dist1<radius*maxDist){
                id[i]=0;
            }
            else if (dist1<radius*maxDist){
                id[i]=1;
            }
            else{
                id[i]=2;
            }
            //id[i]=rad_number*(int)(dist/(maxDist/(double)rho_number))+(int)(angle/(2*Math.PI/rad_number));
            
            
            //mat[param.disk2D[i][0]][param.disk2D[i][1]]=id;
        }
        
        int [] nb = new int[maxId];
        for (int i=0;i<nbDataPerImage;i++){
            nb[id[i]]++;
        }
        for (int i=0;i<maxId;i++){
            if (nb[i]!=0){
                id_table[i]=new int[nb[i]];
                phase_list[i]=new double[nb[i]];
            }
            else{
                id_table[i]=null;
                phase_list[i]=null;
            }
            nb[i]=0;
        }
        for (int i=0;i<nbDataPerImage;i++){
            id_table[id[i]][nb[id[i]]]=i;
            phase_list[id[i]][nb[id[i]]]=phase[i];
            nb[id[i]]++;
        }
        
        
        return maxId;
        
    }
    
    
    
    
    
    
    
    
    
    
    
    public double [] getValuesPhase(int id){
        if ((id>=0)&&(id<maxId)){
            if (phase_list[id]!=null){
                double [] v = new double [phase_list[id].length] ;
                for (int i=0;i<phase_list[id].length;i++){
                    v[i]=phase_list[id][i];
                }
                return v;
            }
            else{
                return null;
            }
        }
        else{
            IJ.log("error position index in getValuesPhasePointer of PhaseJcudaFastDouble class");
            return null;
        }
    }
    
    
    public int [] getPositSplitPointer(int id){
        if ((id>=0)&&(id<maxId)){
            return id_table[id];
        }
        else{
            IJ.log("error position index in getFirstPositSplit of PhaseJcudaFastDouble class");
            return null;
        }
    }
    
    
    public void setValuesPhase(int id,double [] v){
        if ((id>=0)&&(id<maxId)){
            if (v.length==phase_list[id].length){
                for (int i=0;i<v.length;i++){
                    phase_list[id][i]=v[i]%modulo;
                    phase[id_table[id][i]]=v[i]%modulo;
                    
                }
            }
            else{
                IJ.log("error position second index in setValuesPhase of PhaseJcudaFastDouble class");
            }
        }
        else{
            IJ.log("error position index in setValuesPhase of PhaseJcudaFastDouble class");
            
        }
    }
    
    
    public double [] getValuesPhase(){
        double [] v = new double [phase.length] ;
        for (int i=0;i<phase.length;i++){
            v[i]=phase[i];
        }
        return v;
    }
    
    
    public void setValuesPhase(double [] v){
        
        for (int i=0;i<v.length;i++){
            phase[i]=v[i]%modulo;
        }
            
    }
    
    
    
    
    public Pointer getPointerPhase(){
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        cudaMemcpyAsync( device_phase,host_phase,nbDataPerImage*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro modif piston  "+cudaResult);}
        return device_phase;
    }
    
    
    
    public void setValuePixel(int posit,double value){
        
        if ((posit>=0)&&(posit<nbDataPerImage)){
            phase[posit]=value%modulo;
           
        }   
        else{
            IJ.log("error wrong posit setA Zernike function");
        }
    }
    
    
    
    
    public double getValuePixel(int posit){
        if ((posit>=0)&&(posit<nbDataPerImage)){
            return phase[posit];
        }
        IJ.log("error position index in getValuePixel of PhaseJcudaFastDouble class");
        return 0;
    }
    
}
