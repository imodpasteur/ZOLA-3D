/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.cuda;


import ij.IJ;
import jcuda.driver.CUstream;
import jcuda.jcublas.JCublas;
import jcuda.jcublas.JCublas2;
import static jcuda.jcublas.JCublas2.cublasCreate;
import static jcuda.jcublas.JCublas2.cublasSetPointerMode;
import static jcuda.jcublas.JCublas2.cublasGetPointerMode;
import static jcuda.jcublas.JCublas2.cublasSetStream;
import jcuda.jcublas.cublasHandle;
import static jcuda.jcublas.cublasPointerMode.CUBLAS_POINTER_MODE_DEVICE;
import jcuda.jcusparse.JCusparse;
import static jcuda.jcusparse.JCusparse.cusparseCreate;
import jcuda.jcusparse.cusparseHandle;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaStream_t;
/**
 *
 * @author benoit
 */
public class MyCudaStream {
    private static CUstream [] mycustream;
    private static cudaStream_t [] mystream;
    private static cublasHandle [] handlecublas;
    private static cusparseHandle [] handlecusparse;
    
    public static int numberStream;
    Boolean ok=null;
    public static void init(int numberOfStream){
        //IJ.log("init cuda "+numberOfStream);
        //numberOfStream*=2;
        
        
        
        if (numberOfStream>=1){
            

            numberStream=numberOfStream;
            mycustream=new CUstream[numberOfStream];
            mystream=new cudaStream_t[numberOfStream];
            handlecublas=new cublasHandle[numberOfStream];
            handlecusparse=new cusparseHandle[numberOfStream];
            
            for (int i=0;i<numberOfStream;i++){
                //mystream[i];
                mystream[i]= new cudaStream_t();



                //JCuda.cudaStreamCreate(mystream[i]);
                JCuda.cudaStreamCreateWithFlags(mystream[i],JCuda.cudaStreamNonBlocking);
                //JCuda.cudaStreamCreateWithFlags(mystream[i],JCuda.cudaEventDisableTiming);
                
                
                mycustream[i]=new CUstream(mystream[i]);
                
                
                handlecublas[i] = new cublasHandle();
                
                
                // Initialize JCusblas library
                cublasCreate(handlecublas[i]);//in cublas2
                cublasSetPointerMode(handlecublas[i], CUBLAS_POINTER_MODE_DEVICE);
                JCublas.cublasSetKernelStream(mystream[i]);
                cublasSetStream(handlecublas[i], mystream[i]);
                //jcuda.jcublas.JCublas2.setExceptionsEnabled(true);
                
                handlecusparse[i] = new cusparseHandle();
                // Initialize JCusparse library
                cusparseCreate(handlecusparse[i]);
                JCusparse.cusparseSetStream(handlecusparse[i], mystream[i]);
                
                
                
            }
            MyVecDouble.init();
            
        }
        
    }
    
    
    
    public static CUstream getCUstream(int id){
        
        return mycustream[id];
    }
    
    
    public static cudaStream_t getCudaStream_t(int id){
        
        return mystream[id];
    }
    
    
    
    public static cublasHandle getHandleCublas(int id){
        
        return handlecublas[id];
    }
    
    public static cusparseHandle getHandleCuSparse(int id){
        
        return handlecusparse[id];
    }
    
    
    
    
    public static void destroy(){
        //JCufft.cufftDestroy(plan);
        
        for (int i=0;i<numberStream;i++){
            if (handlecusparse[i]!=null){
                JCusparse.cusparseDestroy(handlecusparse[i]);
                JCublas2.cublasDestroy(handlecublas[i]);
                handlecusparse[i]=null;
                //plan=null;
                handlecublas[i]=null;
                
            }
        }
        
        MyVecDouble.shutdown();
    }
    
    
}
