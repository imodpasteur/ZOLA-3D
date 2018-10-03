/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;



import ij.IJ;
import jcuda.jcublas.JCublas2;
import static jcuda.jcublas.JCublas2.cublasCreate;
import static jcuda.jcublas.JCublas2.cublasSetPointerMode;
import jcuda.jcublas.cublasHandle;
import static jcuda.jcublas.cublasPointerMode.CUBLAS_POINTER_MODE_DEVICE;
import jcuda.jcufft.JCufft;
import jcuda.jcufft.cufftHandle;
import jcuda.jcufft.cufftType;
import jcuda.jcusparse.JCusparse;
import static jcuda.jcusparse.JCusparse.cusparseCreate;
import jcuda.jcusparse.cusparseHandle;

/**
 *
 * @author benoit
 */
public class PhaseParameters {
    
    
    
    public int PSF_number=1;
    public String pathcalib="";
    public double sigmaGaussianKernel;
    
    public double centerX;
    public double centerY;
    public double centerZ;
    
    public double centerFourierImage;
    public double centerFourierImage_cpu;
    
    public int size;//taille de l'image Fourier domain
    public int size_cpu;//taille de l'image Fourier domain
    public double sizeRadiusRingPixel;//taille du rayon en pixels
    public double sizeRadiusRingPixel_cpu;//taille du rayon en pixels
    public int sizeDisk;//nombre de pixels dans le disque de phase
    public int sizeDisk_cpu;//nombre de pixels dans le disque de phase
    public int [][] disk2D;//coordonnées du disque 2D
    public int [][] disk2D_cpu;//coordonnées du disque 2D
    public int sizeoutput;//size de l'image espace réel
    
    public double weightZ=1;
    public double weightTMP=1;
    public double xystep;
    public double zstep;
    public double wavelength;
    public double noil;
    public double nwat=1.33;
    public double na;
    public double ringsize;
    
    public double A;//mean numb of photons in phase retrieval
    public double B;//background intensity in phase retrieval
    public int orderGaussianPupil=0;
    public int stream=0;
    public double Zfocus=0;
    public double ZfocusCenter=0;
    
    public boolean zernikedPSF=true;
    
    
    public PhaseParameters(int size,int sizeoutput,int orderGaussianPupil,double xystep,double zstep,double wavelength,double noil,double na,double wz,double sigmaGaussianKernel){
        
        this.sigmaGaussianKernel=sigmaGaussianKernel;
        
        this.orderGaussianPupil=orderGaussianPupil;
        this.weightZ=wz;
        this.xystep=xystep; 
        this.zstep=zstep; 
        this.stream=0;
        this.wavelength=wavelength;
        this.noil=noil;
        
        this.na=na;
        if (size%2!=0){
            size--;
        }
        if (size<=0){
            IJ.log("problem size initialization at 0");
        }
        if (sizeoutput%2!=0){
            sizeoutput--;
        }
    
        if (sizeoutput>size){
            sizeoutput=size;
        }
        
        if (sizeoutput<=0){
            sizeoutput=size;
        }
        this.size=size;
        size_cpu=(int)Math.pow(2, (int)Math.ceil(Math.log(size)/Math.log(2)));
        this.sizeoutput=sizeoutput;
        
        
        this.centerFourierImage=(this.size/2)-.5;
        this.centerFourierImage_cpu=(this.size_cpu/2)-.5;

        ringsize=(na/wavelength)*(na/wavelength);
        sizeRadiusRingPixel=Math.sqrt(ringsize)*(size*xystep);
        sizeRadiusRingPixel_cpu=Math.sqrt(ringsize)*(size_cpu*xystep);
        double center=size/2;
        double center_cpu=size_cpu/2;
        center=centerFourierImage;
        center_cpu=centerFourierImage_cpu;
        sizeDisk=0;
        int id=0;
        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                double disk=((((i-center)/(size*xystep))*((i-center)/(size*xystep))+((ii-center)/(size*xystep))*((ii-center)/(size*xystep))));
                if (disk<=ringsize){
                    sizeDisk++;
                    id++;
                }
            }
        }
        id=0;
        sizeDisk_cpu=0;
        for (int i=0;i<size_cpu;i++){
            for (int ii=0;ii<size_cpu;ii++){
                double disk=((((i-center_cpu)/(size_cpu*xystep))*((i-center_cpu)/(size_cpu*xystep))+((ii-center_cpu)/(size_cpu*xystep))*((ii-center_cpu)/(size_cpu*xystep))));
                if (disk<=ringsize){
                    sizeDisk_cpu++;
                    id++;
                }
            }
        }
        //IJ.log("sizeDisk="+sizeDisk);
        disk2D=new int[sizeDisk][2];
        disk2D_cpu=new int[sizeDisk_cpu][2];
        id=0;
        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                double disk=((((i-center)/(size*xystep))*((i-center)/(size*xystep))+((ii-center)/(size*xystep))*((ii-center)/(size*xystep))));
                
                if (disk<=ringsize){
                    
                    this.disk2D[id][0]=i;
                    this.disk2D[id][1]=ii;
                    id++;
                }
            }
        }
        
        
        id=0;
        for (int i=0;i<size_cpu;i++){
            for (int ii=0;ii<size_cpu;ii++){
                double disk=((((i-center_cpu)/(size_cpu*xystep))*((i-center_cpu)/(size_cpu*xystep))+((ii-center_cpu)/(size_cpu*xystep))*((ii-center_cpu)/(size_cpu*xystep))));
                
                if (disk<=ringsize){
                    
                    this.disk2D_cpu[id][0]=i;
                    this.disk2D_cpu[id][1]=ii;
                    id++;
                }
            }
        }
        
        
        
    }
    
    
    
    
    public PhaseParameters(PhaseParameters paramInit){
        this.pathcalib=paramInit.pathcalib;
        this.sigmaGaussianKernel=paramInit.sigmaGaussianKernel;
        this.centerX=paramInit.centerX;
        this.centerY=paramInit.centerY;
        this.centerZ=paramInit.centerZ;
        this.orderGaussianPupil=paramInit.orderGaussianPupil;
        
        this.A=paramInit.A;
        this.B=paramInit.B;
        
        this.weightZ=paramInit.weightZ;
        this.weightTMP=paramInit.weightTMP;
        this.xystep=paramInit.xystep; 
        this.zstep=paramInit.zstep; 
        
        this.wavelength=paramInit.wavelength;
        this.noil=paramInit.noil;
        this.nwat=paramInit.nwat;
        this.Zfocus=paramInit.Zfocus;
        this.ZfocusCenter=paramInit.ZfocusCenter;
        this.na=paramInit.na;
        
        this.stream=paramInit.stream;
        this.size=paramInit.size;
        this.size_cpu=paramInit.size_cpu;
        
        this.centerFourierImage=(this.size/2)-.5;
        this.centerFourierImage_cpu=(this.size_cpu/2)-.5;
        
        this.sizeoutput=paramInit.sizeoutput;
        
        

        ringsize=paramInit.ringsize;
        sizeRadiusRingPixel=paramInit.sizeRadiusRingPixel;
        
        sizeDisk=paramInit.sizeDisk;
        
        disk2D=new int[sizeDisk][2];
        for (int i=0;i<sizeDisk;i++){
            disk2D[i][0]=paramInit.disk2D[i][0];
            disk2D[i][1]=paramInit.disk2D[i][1];
        }
        
        
        sizeRadiusRingPixel_cpu=paramInit.sizeRadiusRingPixel_cpu;
        
        sizeDisk_cpu=paramInit.sizeDisk_cpu;
        
        disk2D_cpu=new int[sizeDisk_cpu][2];
        for (int i=0;i<sizeDisk_cpu;i++){
            disk2D_cpu[i][0]=paramInit.disk2D_cpu[i][0];
            disk2D_cpu[i][1]=paramInit.disk2D_cpu[i][1];
        }
        
        
    }
    
    
    
    public PhaseParameters(int size, int sizeoutput,PhaseParameters paramInit){
        this.pathcalib=paramInit.pathcalib;
        this.sigmaGaussianKernel=paramInit.sigmaGaussianKernel;
        this.size=size;
        this.size_cpu=(int)Math.pow(2, (int)Math.ceil(Math.log(size)/Math.log(2)));
        
        this.centerFourierImage_cpu=(this.size_cpu/2)-.5;
        
        this.centerFourierImage=(this.size/2)-.5;
        
        this.sizeoutput=sizeoutput;
        this.orderGaussianPupil=paramInit.orderGaussianPupil;
        this.centerX=paramInit.centerX;
        this.centerY=paramInit.centerY;
        this.centerZ=paramInit.centerZ;
        
        this.stream=paramInit.stream;
        this.A=paramInit.A;
        this.B=paramInit.B;
        
        this.weightZ=paramInit.weightZ;
        this.weightTMP=paramInit.weightTMP;
        this.xystep=paramInit.xystep; 
        this.zstep=paramInit.zstep; 
        
        this.wavelength=paramInit.wavelength;
        this.noil=paramInit.noil;
        this.nwat=paramInit.nwat;
        this.Zfocus=paramInit.Zfocus;
        this.ZfocusCenter=paramInit.ZfocusCenter;
        this.na=paramInit.na;
        
        
        //IJ.log("nwat copy "+nwat+"  "+paramInit.nwat);
        
        

        ringsize=paramInit.ringsize;
        sizeRadiusRingPixel=Math.sqrt(ringsize)*(size*xystep);
        double center=size/2;
        center=this.centerFourierImage;
        sizeDisk=0;
        int id=0;
        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                double disk=((((i-center)/(size*xystep))*((i-center)/(size*xystep))+((ii-center)/(size*xystep))*((ii-center)/(size*xystep))));
                if (disk<=ringsize){
                    sizeDisk++;
                    id++;
                }
            }
        }
        //IJ.log("sizeDisk="+sizeDisk);
        disk2D=new int[sizeDisk][2];
        id=0;
        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                double disk=((((i-center)/(size*xystep))*((i-center)/(size*xystep))+((ii-center)/(size*xystep))*((ii-center)/(size*xystep))));
                
                if (disk<=ringsize){
                    
                    this.disk2D[id][0]=i;
                    this.disk2D[id][1]=ii;
                    id++;
                }
            }
        }
        
        
        
        
        
        
        
        sizeRadiusRingPixel_cpu=Math.sqrt(ringsize)*(size_cpu*xystep);
        double center_cpu=size_cpu/2;
        center_cpu=this.centerFourierImage_cpu;
        sizeDisk_cpu=0;
        id=0;
        for (int i=0;i<size_cpu;i++){
            for (int ii=0;ii<size_cpu;ii++){
                double disk=((((i-center_cpu)/(size_cpu*xystep))*((i-center_cpu)/(size_cpu*xystep))+((ii-center_cpu)/(size_cpu*xystep))*((ii-center_cpu)/(size_cpu*xystep))));
                if (disk<=ringsize){
                    sizeDisk_cpu++;
                    id++;
                }
            }
        }
        //IJ.log("sizeDisk="+sizeDisk);
        disk2D_cpu=new int[sizeDisk_cpu][2];
        id=0;
        for (int i=0;i<size_cpu;i++){
            for (int ii=0;ii<size_cpu;ii++){
                double disk=((((i-center_cpu)/(size_cpu*xystep))*((i-center_cpu)/(size_cpu*xystep))+((ii-center_cpu)/(size_cpu*xystep))*((ii-center_cpu)/(size_cpu*xystep))));
                
                if (disk<=ringsize){
                    
                    this.disk2D_cpu[id][0]=i;
                    this.disk2D_cpu[id][1]=ii;
                    id++;
                }
            }
        }
        
        
    }
    
    
    
    
    public void updateweightZ(double weight){
        weightZ=weight;
    }
    
    
    
    public double getweightZ(){
        return weightZ;
    }
    
    
    
    public void updateSigmaGaussianKernel(double sigma){
        this.sigmaGaussianKernel=sigma;
    }
    
    
    
    public double getSigmaGaussianKernel(){
        return this.sigmaGaussianKernel;
    }
    
    
    public void updateweightTMP(double weight){
        weightTMP=weight;
    }
    
    
    
    public double getweightTMP(){
        return weightTMP;
    }
    
    
}
