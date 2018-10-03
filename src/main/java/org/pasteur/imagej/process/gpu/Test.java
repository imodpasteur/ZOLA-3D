/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.cublasHandle;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_N;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_C;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import org.pasteur.imagej.cuda.MyCudaStream;

/**
 *
 * @author benoit
 */
public class Test {
    
    Pointer device_tmp;
    Pointer device_alpha;
    Pointer device_beta;
    Pointer device_res;
    Pointer device_ones;
    cublasHandle  handlecublas;
    
    public Test(){
        double [] res=new double [6];
        double [] value=new double[]{1,2,3, 4,5,6, 7,8,9, 10,11,12};//res should be 5,7,9, 17,19,21
        int size=value.length;
        
        double [] ones=new double[]{1,1,1};
        
        handlecublas=MyCudaStream.getHandleCublas(0);
        
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.DOUBLE);
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));
        
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.DOUBLE);
        //int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(dparam[0].param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda PSF localize() 0 "+cudaResult);}
        cudaMemcpyAsync(device_beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));
        
        device_tmp = new Pointer();
        cudaMalloc(device_tmp, size* Sizeof.DOUBLE);
        cudaMemcpyAsync(device_tmp, Pointer.to(value), size*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));
        
        
        device_res = new Pointer();
        cudaMalloc(device_res, 6* Sizeof.DOUBLE);
        
        device_ones = new Pointer();
        cudaMalloc(device_ones, 3* Sizeof.DOUBLE);
        cudaMemcpyAsync(device_ones, Pointer.to(ones), 3*Sizeof.DOUBLE, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));
        
        int row=6;
        int col=2;
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,row,col,device_alpha,device_tmp,row,device_ones,1,device_beta,device_res,1);
        
        cudaMemcpyAsync(Pointer.to(res), device_res, 6*Sizeof.DOUBLE, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(0));
        
        for (int i=0;i<6;i++){
            IJ.log("res:"+res[i]);
        }
        
    }
}
