/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;
import org.pasteur.imagej.process.*;


import ij.IJ;

/**
 *
 * @author benoit
 */
public class Phase {
    

    
    
    
    
    private double modulo=2*Math.PI;
    public int angleNumber=0;
    public int distNumber=0;
    
    
    
    int cudaResult;
    
    double [] phase;
    
    
    
    //Pointer device_id;
    //Pointer host_id;
    int [] id;
    int [][] id_table;//  first dim: idNumber   ; second dim: the corresponding positions
    double [][] phase_list;//  first dim: idNumber   ; second dim: the corresponding positions
    
    PhaseParameters param;
    int nbDataPerImage;
    
    
    
    double [] tmp;
    
    
    int maxId=0;
    
    Phase(PhaseParameters param){
        
        
        this.param=param;
        
        
        
        nbDataPerImage=param.sizeDisk_cpu;
        
        id=new int [nbDataPerImage];
        //host_id=Pointer.to(id);
        
        
        
        phase=new double [nbDataPerImage]; 
        
        
        
        //device_id=new Pointer();
        //cudaResult = JCuda.cudaMalloc(device_id, nbDataPerImage * Sizeof.DOUBLE);
        //if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 8");}
        
        tmp=new double [1];
        
        
    }
    
    
    
    
    public int createSplits(int rad_number, int rho_number){
        this.angleNumber=rad_number;
        this.distNumber=rho_number;
        
        maxId=rad_number*rho_number;
        
        id_table= new int [maxId][];
        phase_list= new double [maxId][];
        double maxDist=param.sizeRadiusRingPixel_cpu;
        
        //int [][] mat = new int [param.size][param.size];
        double angle,dist,x,y;
        
        for (int i=0;i<nbDataPerImage;i++){
            x=param.disk2D_cpu[i][0]-param.centerFourierImage;
            y=param.disk2D_cpu[i][1]-param.centerFourierImage;
            if (x!=0){
                angle=Math.atan2(y, x)+Math.PI;
            }
            else{
                angle=0;
            }
            dist=Math.sqrt(x*x+y*y);
            
            id[i]=rad_number*(int)(dist/(maxDist/(double)rho_number))+(int)(angle/(2*Math.PI/rad_number));
            
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
    
    public void setValuesPhase(int id,double [] v){
        if ((id>=0)&&(id<maxId)){
            if (v.length==phase_list[id].length){
                for (int i=0;i<v.length;i++){
                    phase_list[id][i]=v[i];
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
    
    
    public void setValuesPhase(double [] v){
        
        for (int i=0;i<v.length;i++){
            phase[i]=v[i]%modulo;
        }
            
    }
    
    
    
    
    public double [] getPhasePointer(){
        return phase;
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
