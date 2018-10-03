/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;


import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.utils.Zernike;
import org.pasteur.imagej.utils.GaussianElimination;
import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas2;
import static  jcuda.jcublas.JCublas2.cublasCreate;
import jcuda.jcublas.cublasHandle;
import jcuda.jcublas.cublasOperation;
import static jcuda.jcusparse.JCusparse.cusparseCreate;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaError;




 
 
import static jcuda.jcublas.JCublas2.cublasCreate;
import static jcuda.jcublas.JCublas2.cublasDestroy;
import static jcuda.jcublas.JCublas2.cublasDgemm;
import static jcuda.jcublas.JCublas2.cublasGetVector;
import static jcuda.jcublas.JCublas2.cublasSetPointerMode;
import static jcuda.jcublas.JCublas2.cublasSetVector;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_N;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.jcublas.cublasPointerMode.CUBLAS_POINTER_MODE_DEVICE;
import static jcuda.runtime.JCuda.cudaFree;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import static jcuda.runtime.cudaMemcpyKind.*;
 
import java.util.Random;
 
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUstream;
import jcuda.jcublas.JCublas;
import jcuda.jcublas.JCublas2;
import jcuda.jcublas.cublasHandle;
import static jcuda.jcusparse.JCusparse.cusparseDgthr;
import jcuda.jcusparse.cusparseHandle;
import static jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaMemcpyKind;

/**
 *
 * @author benoit
 */
public class ZernikePhase_{
     
    cublasHandle  handlecublas;
    cusparseHandle handlecusparse;
    
    CUstream custream;
            
    int cudaResult;
    
    Pointer device_combination;
    Pointer  device_zernike;
    Pointer  device_tmp;
    Pointer  device_tmp2;
//    Pointer  device_zer1;
//    Pointer  device_zer2;
//    Pointer  device_zer4;
    Pointer  device_a;
    //Pointer [] device_A;
    Pointer  device_beta;
    Pointer  device_alpha;
    
    
    Pointer  host_tmp;
    Pointer  host_zernike;
    Pointer  host_a;
    //Pointer [] host_A;
    Pointer hostT;
    double [] ones;
    
    
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
    //method
    //=0 -> use all the coefs
    //=1 -> use only symetric coefs
    //=2 -> use only even coefs
    
    
    
    public ZernikePhase_(PhaseParameters param,int zernikePolyNumber,int method){
        this.multiTrainingWithParabola=false;
        this.param=param;
        Zernike z = new Zernike(param.size,param.sizeRadiusRingPixel,zernikePolyNumber);
        //ImageShow.imshow(z.Z,"Z");
        
        numberOfCoef=zernikePolyNumber;
//        for (int i=0;i<z.Z[0].length;i++){
//            for (int ii=0;ii<z.Z[0][0].length;ii++){
//                if (z.Z.length>2){
//                IJ.log("Z "+i+"  "+ii+"  "+z.Z[2][i][ii]+"  ");
//                }
//            }
//        }
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
        
        nbDataPerImage=param.sizeDisk;
        
        Zvect=new double [nbDataPerImage*numberOfCoef];
        this.a=new double [numberOfCoef];
        //this.A=new double [dim][numberOfCoef];
        this.coef=new int [numberOfCoef];
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
                        Zvect[id*nbDataPerImage+i]=z.Z[rz][param.disk2D[i][0]][param.disk2D[i][1]];
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
                        Zvect[id*nbDataPerImage+i]=z.Z[rz][param.disk2D[i][0]][param.disk2D[i][1]];
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
                    Zvect[rz*nbDataPerImage+i]=z.Z[rz][param.disk2D[i][0]][param.disk2D[i][1]];
                }
            }
        }
        
        
        
        builder();
        
    }
    
    
    public void setMatAtPosit(Pointer k,int posit){
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda ZER setkZ 1 "+cudaResult+"   "+param.stream);}
        
        cudaMemcpyAsync(device_zernike.withByteOffset(posit*nbDataPerImage*Sizeof.DOUBLE), k, nbDataPerImage*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda ZER setkZ 2 "+cudaResult+"   "+param.stream);}
    }
    
    
    
    
    
    
    public ZernikePhase_(PhaseParameters param,int [] coef){
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
        Zernike z = new Zernike(param.size,param.sizeRadiusRingPixel,maxiCoef+1);
        
        numberOfCoef=coef.length;
        
            
        
        nbDataPerImage=param.sizeDisk;
        
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
                    Zvect[id*nbDataPerImage+i]=z.Z[rz][param.disk2D[i][0]][param.disk2D[i][1]];
                }
                id++;
            }
            
            
        }
        
        
        
        /*host_zernike=Pointer.to(Zvect);
        host_a=Pointer.to(a);
        
        tmp=new double [1];
        host_tmp=Pointer.to(tmp);
        device_zernike = new Pointer();
        int cudaResult = JCuda.cudaMalloc(device_zernike, nbDataPerImage*numberOfCoef * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        
        
        
        
        device_zer1 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_zer1, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        device_zer2 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_zer2, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        device_zer4 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_zer4, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        device_tmp = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        device_tmp2 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp2, nbDataPerImage*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        
        double [] one = new double [1];
        one[0]=1;
        Pointer hostOne=Pointer.to(one);
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.DOUBLE);
        cudaMemcpyAsync(device_alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.DOUBLE);
        cudaMemcpyAsync(device_beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice);
        
        
        
        device_a = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_a, numberOfCoef * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        
        device_combination= new Pointer();
        cudaResult = JCuda.cudaMalloc(device_combination, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda");return ;
        }
        
        
        JCuda.cudaMemcpyAsync(device_a, host_a, numberOfCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice);
        JCuda.cudaMemcpyAsync(device_zernike, host_zernike, nbDataPerImage*numberOfCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice);
        
        
        
        JCuda.cudaMemcpyAsync(device_zer1, device_zernike.withByteOffset(1*nbDataPerImage*Sizeof.DOUBLE), nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
        JCuda.cudaMemcpyAsync(device_zer2, device_zernike.withByteOffset(2*nbDataPerImage*Sizeof.DOUBLE), nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
        JCuda.cudaMemcpyAsync(device_zer4, device_zernike.withByteOffset(4*nbDataPerImage*Sizeof.DOUBLE), nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
        
        
        double [] ones = new double [nbDataPerImage];
        for (int i=0;i<ones.length;i++){
            ones[i]=1;
        }
        Pointer hostT = Pointer.to(ones);
        JCuda.cudaMemcpyAsync(device_tmp2, hostT, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice);//half first part is 1
        
        
        
        m=nbDataPerImage;
        numCoef=numberOfCoef;
        incr=1;
        lda=m;*/
        
        builder();
        
        
    }
    
    
    
    
    
    void reLaunchCudaTool(){
        
        
        
        handlecublas=MyCudaStream.getHandleCublas(param.stream);
        handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
        custream=MyCudaStream.getCUstream(param.stream);
        
    }
    
    
    
    private void builder(){
        
        
        host_zernike=Pointer.to(Zvect);
        host_a=Pointer.to(a);
//        host_A=new Pointer[dim];
//        for (int or=0;or<dim;or++){
//            host_A[or]=Pointer.to(A[or]);
//        }
        
        tmp=new double [1];
        host_tmp=Pointer.to(tmp);
        device_zernike = new Pointer();
        int cudaResult = JCuda.cudaMalloc(device_zernike, nbDataPerImage*numberOfCoef * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda build 1");return ;
        }
        
        
        
        
        
//        device_zer1 = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_zer1, nbDataPerImage * Sizeof.DOUBLE);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            IJ.log("ERROR malloc cuda  build 2");return ;
//        }
//        device_zer2 = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_zer2, nbDataPerImage * Sizeof.DOUBLE);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            IJ.log("ERROR malloc cuda build 3");return ;
//        }
//        device_zer4 = new Pointer();
//        cudaResult = JCuda.cudaMalloc(device_zer4, nbDataPerImage * Sizeof.DOUBLE);
//        if (cudaResult != cudaError.cudaSuccess)
//        {
//            IJ.log("ERROR malloc cuda build 4");return ;
//        }
        device_tmp = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda build 5");return ;
        }
        device_tmp2 = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_tmp2, nbDataPerImage*2 * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda build 6");return ;
        }
        
        double [] one = new double [1];
        one[0]=1;
        Pointer hostOne=Pointer.to(one);
        device_alpha = new Pointer();
        cudaMalloc(device_alpha, 1 * Sizeof.DOUBLE);
        cudaMemcpyAsync(device_alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        
        device_beta = new Pointer();
        cudaMalloc(device_beta, 1 * Sizeof.DOUBLE);
        cudaMemcpyAsync(device_beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        
        
////////////////        device_A=new Pointer[dim];
////////////////        for (int or=0;or<dim;or++){
////////////////            device_A[or] = new Pointer();
////////////////            cudaResult = JCuda.cudaMalloc(device_A[or], numberOfCoef * Sizeof.DOUBLE);
////////////////            if (cudaResult != cudaError.cudaSuccess)
////////////////            {
////////////////                IJ.log("ERROR malloc cuda");return ;
////////////////            }
////////////////        }
        
        device_a = new Pointer();
        cudaResult = JCuda.cudaMalloc(device_a, numberOfCoef * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess)
        {
            IJ.log("ERROR malloc cuda build 7");return ;
        }
        
        device_combination=new Pointer();
        cudaResult = JCuda.cudaMalloc(device_combination, nbDataPerImage * Sizeof.DOUBLE);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 8");}
//        for (int or=0;or<dim;or++){
//            device_combination[or]= new Pointer();
//            cudaResult = JCuda.cudaMalloc(device_combination[or], nbDataPerImage * Sizeof.DOUBLE);
//            if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda");}
//        }
        
//        for (int or=0;or<dim;or++){
//            JCuda.cudaMemcpyAsync(device_A[or], host_A[or], numberOfCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
//        }
        JCuda.cudaMemcpyAsync(device_a, host_a, numberOfCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult = JCuda.cudaMemcpyAsync(device_zernike, host_zernike, nbDataPerImage*numberOfCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 9");}
        
        
//        cudaResult = JCuda.cudaMemcpyAsync(device_zer1, device_zernike.withByteOffset(1*nbDataPerImage*Sizeof.DOUBLE), nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 10");}
//        cudaResult = JCuda.cudaMemcpyAsync(device_zer2, device_zernike.withByteOffset(2*nbDataPerImage*Sizeof.DOUBLE), nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 11");}
//        cudaResult = JCuda.cudaMemcpyAsync(device_zer4, device_zernike.withByteOffset(4*nbDataPerImage*Sizeof.DOUBLE), nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 12");}
        
        
        ones = new double [nbDataPerImage];
        for (int i=0;i<ones.length;i++){
            ones[i]=1;
        }
        hostT = Pointer.to(ones);
        cudaResult = JCuda.cudaMemcpyAsync(device_tmp2, hostT, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda build 13");}
        
        
        
        m=nbDataPerImage;
        numCoef=numberOfCoef;
        incr=1;
        lda=m;
        
        
        
        
        handlecublas=MyCudaStream.getHandleCublas(param.stream);
        handlecusparse=MyCudaStream.getHandleCuSparse(param.stream);
        custream=MyCudaStream.getCUstream(param.stream);
        
    }
    
    
    
    public void init(Pointer device_kz){
    
        
        
        double [] y = new double [this.nbDataPerImage];
        Pointer hostKz = Pointer.to(y);
        cudaResult = JCuda.cudaMemcpyAsync(hostKz,device_kz, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda!");}
        
        
        double [][] x = new double [this.nbDataPerImage][this.numCoef];
        
        
        for (int id=0;id<numCoef;id++){
            for (int i=0;i<nbDataPerImage;i++){
                x[i][id]=Zvect[id*nbDataPerImage+i];
            }
        }
        
        double [][] xx = new double [this.numCoef][this.numCoef];
        
        for (int id=0;id<numCoef;id++){
            for (int iid=0;iid<numCoef;iid++){
                xx[id][iid]=0;
                for (int i=0;i<nbDataPerImage;i++){
                    xx[id][iid]+=x[i][id]*x[i][iid];
                }
            }
        }
        
        double [] xy = new double [this.numCoef];
        
        for (int id=0;id<numCoef;id++){
            
            xy[id]=0;
            for (int i=0;i<nbDataPerImage;i++){
                xy[id]+=x[i][id]*y[i];
            }
            
        }
        
        for (int i=0;i<xy.length;i++){
        }
        
        double [] res=GaussianElimination.lsolve(xx, xy);
        //String s="";
        //for (int i=0;i<res.length;i++){
        //    s+=res[i]+"  ";
        //}
        //IJ.log("A init "+s);
        this.setA(res);
    }
    
    
    public void free(){
        
    
        
        JCuda.cudaFree(device_zernike);
        //for (int or=0;or<dim;or++){
            JCuda.cudaFree(device_combination);
            JCuda.cudaFree(device_a);
        //}
        ////////////////////////////JCuda.cudaFree(device_a);
        JCuda.cudaFree(device_alpha);
        JCuda.cudaFree(device_beta);
        
        JCuda.cudaFree(device_tmp);
        JCuda.cudaFree(device_tmp2);
//        JCuda.cudaFree(device_zer1);
//        JCuda.cudaFree(device_zer2);
//        JCuda.cudaFree(device_zer4);
        
    }
    
    
    
    //matrix order with CUBLAS_OP_N;
    // 0 3 6
    // 1 4 7
    // 2 5 8
    
    public Pointer computeCombination(){
        //IJ.log("truc "+MyCudaStream.numberStream);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_a,incr,device_beta,device_combination,incr);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 1 phase "+cudaResult+"   "+param.stream);}
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        
        return device_combination;
    }
    
    
    
    
    
    public Pointer computeCombinationPlusOtherPhase(Pointer otherPhase){
        //IJ.log("truc "+MyCudaStream.numberStream);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_a,incr,device_beta,device_combination,incr);
        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 1 phase "+cudaResult+"   "+param.stream);}
        
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        MyVecDouble.add(custream, this.nbDataPerImage, device_combination, device_combination, otherPhase);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        
        return device_combination;
    }
    
    
    
    
    
////////////    public Pointer [] computeCombination(){
////////////        
////////////        for (int or=0;or<dim;or++){
////////////            
////////////            //device_alpha instead of device_beta (0->1)
////////////            
////////////            cudaResult=jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_A[or],incr,device_beta,device_combination[or],incr);
////////////            if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR combination cuda");}
////////////        }
////////////        return device_combination;
////////////    }
    
    
    /*
    public Pointer [] computeCombinationPositive(){
        for (int or=0;or<order;or++){
            jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_A[or],incr,device_beta,device_combination[or],incr);

            VecDouble.gtScalar(nbDataPerImage, device_tmp, device_combination[or], 0);
            VecDouble.mul(nbDataPerImage, device_combination[or], device_combination[or], device_tmp);


            //new part
            VecDouble.mul(nbDataPerImage, device_combination[or], device_combination[or], device_combination[or]);

            
            Pointer p = new Pointer();
            JCuda.cudaMalloc(p, 1 * Sizeof.DOUBLE);
            double [] som = new double [1];
            JCublas2.cublasDasum(handlecublas,nbDataPerImage,device_combination[or],1,p);//norm sum=1
            JCuda.cudaMemcpyAsync(Pointer.to(som), p, 1 * Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);

            VecDouble.divScalar(nbDataPerImage, device_combination[or], device_combination[or],som);//square
            VecDouble.sqrt(nbDataPerImage, device_combination[or], device_combination[or]);
        }
        return device_combination;
    }
    */
    
    //used for derivative according to posit
////////////////////    public Pointer [] computeCombination(int orderposit,int posit,double shift){
////////////////////        cudaMemcpyAsync( host_tmp,device_A[orderposit].withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
////////////////////        tmp[0]+=shift;
////////////////////        cudaResult=cudaMemcpyAsync( device_A[orderposit].withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
////////////////////        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
////////////////////        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));
////////////////////        device_combination=this.computeCombination();
////////////////////        //jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_A[0],incr,device_beta,device_combination[0],incr);
////////////////////        
////////////////////        
////////////////////        tmp[0]-=shift;
////////////////////        cudaResult=cudaMemcpyAsync( device_A[orderposit].withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
////////////////////        if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
////////////////////        JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));
////////////////////        return device_combination;
////////////////////    }
    
    
    
    //used for derivative according to posit
    public Pointer computeCombination(int posit,double shift){
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        cudaMemcpyAsync( host_tmp,device_a.withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        tmp[0]+=shift;
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        //IJ.log("tmp "+tmp[0]+"  "+posit);
        cudaMemcpyAsync( device_a.withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_a,incr,device_beta,device_combination,incr);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        tmp[0]-=shift;
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        cudaMemcpyAsync( device_a.withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        return device_combination;
    }
    
    
    //used for derivative according to posit
    public Pointer computeCombinationPlusOtherPhase(Pointer otherPhase, int posit,double shift){
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        cudaMemcpyAsync( host_tmp,device_a.withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
        int cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        tmp[0]+=shift;
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        //IJ.log("tmp "+tmp[0]+"  "+posit);
        cudaMemcpyAsync( device_a.withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_a,incr,device_beta,device_combination,incr);
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        tmp[0]-=shift;
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        MyVecDouble.add(custream, this.nbDataPerImage, device_combination, device_combination, otherPhase);
        cudaMemcpyAsync( device_a.withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda zernike  "+cudaResult);}
        
        return device_combination;
    }
    
    
    
    
    public void modifyOnePixelOfPiston(int posit,double newvalue){
        
        
        if ((posit>=0)&&(posit<this.nbDataPerImage)){
            this.Zvect[posit]=newvalue;
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro modif piston  "+cudaResult);}
            
            
            cudaResult = JCuda.cudaMemcpyAsync(device_zernike.withByteOffset(posit* Sizeof.DOUBLE), host_zernike.withByteOffset(posit* Sizeof.DOUBLE), 1* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR malloc cuda modif piston 9");}
            
            
            
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro modif piston  "+cudaResult);}
        }
        else{
            IJ.log("oops: seg fault getPistonValue ZernikePhaseJcuda");
        }
        
        
        
        
        
        
    }
    
    public double getPistonValue(int posit){
        if ((posit>=0)&&(posit<this.nbDataPerImage)){
            return this.Zvect[posit];
        }
        else{
            IJ.log("oops: seg fault getPistonValue ZernikePhaseJcuda");
            return 0;
        }
    }
    
    
   /* //used for derivative according to posit
    not possible because division by zero if all elements equal 0 -> NaN values
    public Pointer [] computeCombinationPositive(int orderposit,int posit,double shift){
        
        cudaMemcpyAsync( host_tmp,device_A[orderposit].withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost);
        tmp[0]+=shift;
        cudaMemcpyAsync( device_A[orderposit].withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice);
        
        
        //jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_N,m,numCoef,device_alpha,device_zernike,m,device_A[orderposit],incr,device_beta,device_combination,incr);
        device_combination=computeCombinationPositive();
        
        tmp[0]-=shift;
        cudaMemcpyAsync( device_A[orderposit].withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice);
        
        VecDouble.gtScalar(nbDataPerImage, device_tmp, device_combination[orderposit], 0);
        VecDouble.mul(nbDataPerImage, device_combination[orderposit], device_combination[orderposit], device_tmp);
        
        return device_combination;
    }*/
    
    
//////////////    public void setA(int positorder,int posit,double value){
//////////////        if ((posit>=0)&&(posit<numCoef)){
//////////////            cudaResult =cudaMemcpyAsync( host_tmp,device_A[positorder].withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}//unuseful
//////////////            tmp[0]=value;
//////////////            cudaResult =cudaMemcpyAsync( device_A[positorder].withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
//////////////        }
//////////////        else{
//////////////            IJ.log("error wrong posit setA Zernike function");
//////////////        }
//////////////    }
    
    
    
    public void setA(int posit,double value){
        if ((posit>=0)&&(posit<numCoef)){
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
            cudaMemcpyAsync( host_tmp,device_a.withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
            tmp[0]=value;
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
            cudaMemcpyAsync( device_a.withByteOffset(posit*Sizeof.DOUBLE),host_tmp,1*Sizeof.DOUBLE,cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
            cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        
        }   
        else{
            IJ.log("error wrong posit setA Zernike function");
        }
    }
    
//////////////    public void setApointer(Pointer [] device){
//////////////        
//////////////        for (int or=0;or<dim;or++){
//////////////            cudaResult=JCuda.cudaMemcpyAsync(device_A[or], device[or], numCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
//////////////        }
//////////////        
//////////////        
//////////////    }
    
    
    public void setApointer(Pointer device){
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        JCuda.cudaMemcpyAsync(device_a, device, numCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));
        cudaResult=JCuda.cudaStreamSynchronize(MyCudaStream.getCudaStream_t(param.stream));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR synchro cuda Zerknike 2 phase "+cudaResult+"   "+param.stream);}
        
    }
    
////////////    public Pointer [] getApointer(){
////////////        
////////////        return device_A;
////////////    }
    
    
    public Pointer getApointer(){
        
        return device_a;
    }
    
//////////    public void setA(double [] value){
//////////        if (value.length!=this.numCoef){
//////////            IJ.log("error size setA Zernike function");
//////////        }
//////////        else{
//////////            //String aa="";
//////////            for (int i=0;i<numCoef;i++){
//////////                
//////////                this.A[0][i]=value[i];
//////////                //aa+="  "+value[i];
//////////            }
//////////            //IJ.log("aa "+aa);
//////////            host_A[0]=Pointer.to(A[0]);
//////////            cudaResult=JCuda.cudaMemcpyAsync(device_A[0], host_A[0], numCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
//////////        }
//////////        
//////////    }
//////////    
//////////    
//////////    public void setA(double [][] value){
//////////        if ((value[0].length!=this.numCoef)&&(value.length!=this.dim)){
//////////            IJ.log("error size setA Zernike function");
//////////        }
//////////        else{
//////////            //String aa="";
//////////            for (int i=0;i<numCoef;i++){
//////////                for (int or=0;or<dim;or++){
//////////                    this.A[or][i]=value[or][i];
//////////                }
//////////                //aa+="  "+value[i];
//////////            }
//////////            //IJ.log("aa "+aa);
//////////            for (int or=0;or<dim;or++){
//////////                host_A[or]=Pointer.to(A[or]);
//////////                cudaResult=JCuda.cudaMemcpyAsync(device_A[or], host_A[or], numCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
//////////            }
//////////        }
//////////        
//////////    }
    
    
    
    public void setA(double [] value){
        if (value.length!=this.numCoef){
            IJ.log("error size setA Zernike function");
        }
        else{
            //String aa="";
            for (int i=0;i<numCoef;i++){
                
                this.a[i]=value[i];
                //aa+="  "+value[i];
            }
            //IJ.log("aa "+aa);
            host_a=Pointer.to(a);
            JCuda.cudaMemcpyAsync(device_a, host_a, numCoef* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
        }
        
    }
    
    
    
    
    
//////////    public double [][] getA(){
//////////        for (int or=0;or<dim;or++){
//////////            cudaResult=cudaMemcpyAsync( host_A[or],device_A[or],numCoef*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
//////////        }
//////////        return A;
//////////    }
    
    
    public double [] getA(){
        cudaMemcpyAsync( host_a,device_a,numCoef*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
        return a;
    }
    
    
    public double getA(int posit){
        cudaResult=cudaMemcpyAsync( host_tmp,device_a.withByteOffset(posit*Sizeof.DOUBLE),1*Sizeof.DOUBLE,cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream)); if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memcpy cuda");}
        return tmp[0];
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//    
//    
//    
//    
//    
//    
//    
//    public double computeCenterX(Pointer device_kx){
//        
//        JCuda.cudaMemcpyAsync(device_tmp2, hostT, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        //coef zernike 2
//        //JCuda.cudaMemsetAsync(device_tmp2, 0, this.nbDataPerImage * Sizeof.DOUBLE);
//        JCuda.cudaMemcpyAsync(device_tmp2.withByteOffset(nbDataPerImage*Sizeof.DOUBLE), device_kx, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        
//        
//        Pointer sqMat = new Pointer();
//        JCuda.cudaMalloc(sqMat, 4 * Sizeof.DOUBLE);
//        
//        //dfghdrfeghjcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T, CUBLAS_OP_N,2,2,m,device_alpha,device_tmp2,m,device_tmp2,m,device_beta,sqMat,2);
//                           
//        //jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
//        
//        double [] mat = new double [4];
//        Pointer hmat=Pointer.to(mat);
//        JCuda.cudaMemcpyAsync(hmat, sqMat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
//        
//        double det=(mat[0]*mat[3]-mat[1]*mat[2]);
//        if (det!=0){
//            det=1/det;
//        }
//        else{
//            IJ.log("error inversion matrix det=0");
//            return 0;
//        }
//        double tmp=mat[3];
//        mat[3]=det*mat[0];
//        mat[0]=det*tmp;
//        mat[1]=-det*mat[1];
//        mat[2]=-det*mat[2];
//        //IJ.log("mat inv "+mat[0]+"  "+mat[1]+"  "+mat[2]+"  "+mat[3]+"  size : "+this.nbDataPerImage);
//        JCuda.cudaMemcpyAsync(sqMat,hmat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_T,2,m,2,device_alpha,sqMat,2,device_tmp2,m,device_beta,device_tmp2,2);
//                          
//        
//        VecDouble.mulScalar(this.nbDataPerImage, device_tmp, device_zer2, this.getA(2));
//        
//        
//        Pointer dbMat = new Pointer();
//        JCuda.cudaMalloc(dbMat, 2 * Sizeof.DOUBLE);
//        
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_N,2,1,m,device_alpha,device_tmp2,2,device_tmp,m,device_beta,dbMat,2);
//         
//        
//        
//        double [] res =new double [2];
//        Pointer pres= Pointer.to(res);
//        JCuda.cudaMemcpyAsync(pres, dbMat,  2*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
////        for (int i=0;i<res.length;i++){
////            IJ.log("res "+res[i]+"  ");
////        }
//        
//        //add or substract to X res[1]
//        if (res[1]!=0)
//            IJ.log("shift X "+res[1]+"   "+(1/res[1]));
//        else
//            IJ.log("shift X "+res[1]);
//        
//        JCuda.cudaFree(sqMat);
//        return res[1];
//    }
//    
//    
//    
//    public double computeCenterY(Pointer device_ky){
//        
//        JCuda.cudaMemcpyAsync(device_tmp2, hostT, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        //coef zernike 1
//        
//        JCuda.cudaMemcpyAsync(device_tmp2.withByteOffset(nbDataPerImage*Sizeof.DOUBLE), device_ky, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        
//        
//        Pointer sqMat = new Pointer();
//        JCuda.cudaMalloc(sqMat, 4 * Sizeof.DOUBLE);
//        
//        //dfghdrfeghjcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T, CUBLAS_OP_N,2,2,m,device_alpha,device_tmp2,m,device_tmp2,m,device_beta,sqMat,2);
//                           
//        //jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
//        
//        double [] mat = new double [4];
//        Pointer hmat=Pointer.to(mat);
//        JCuda.cudaMemcpyAsync(hmat, sqMat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
//        
//        double det=(mat[0]*mat[3]-mat[1]*mat[2]);
//        if (det!=0){
//            det=1/det;
//        }
//        else{
//            IJ.log("error inversion matrix det=0");
//            return 0;
//        }
//        double tmp=mat[3];
//        mat[3]=det*mat[0];
//        mat[0]=det*tmp;
//        mat[1]=-det*mat[1];
//        mat[2]=-det*mat[2];
//        //IJ.log("mat inv "+mat[0]+"  "+mat[1]+"  "+mat[2]+"  "+mat[3]+"  size : "+this.nbDataPerImage);
//        JCuda.cudaMemcpyAsync(sqMat,hmat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_T,2,m,2,device_alpha,sqMat,2,device_tmp2,m,device_beta,device_tmp2,2);
//                          
//        
//        VecDouble.mulScalar(this.nbDataPerImage, device_tmp, device_zer1, this.getA(1));
//        
//        
//        Pointer dbMat = new Pointer();
//        JCuda.cudaMalloc(dbMat, 2 * Sizeof.DOUBLE);
//        
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_N,2,1,m,device_alpha,device_tmp2,2,device_tmp,m,device_beta,dbMat,2);
//         
//        
//        
//        double [] res =new double [2];
//        Pointer pres= Pointer.to(res);
//        JCuda.cudaMemcpyAsync(pres, dbMat,  2*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
////        for (int i=0;i<res.length;i++){
////            IJ.log("res "+res[i]+"  ");
////        }
//        
//        //add or substract to X res[1]
//        if (res[1]!=0)
//            IJ.log("shift Y "+res[1]+"   "+(1/res[1]));
//        else
//            IJ.log("shift Y "+res[1]);
//        JCuda.cudaFree(sqMat);
//        return res[1];
//        
//    }
//    
//    
//    
//    //fit kz to zernike coef number 4 (parabola)
//    public double computeCenterZ(Pointer device_kz){
//        
//        //coef zernike 4
//        JCuda.cudaMemcpyAsync(device_tmp2, hostT, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        JCuda.cudaMemcpyAsync(device_tmp2.withByteOffset(nbDataPerImage*Sizeof.DOUBLE), device_kz, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        
//        
//        Pointer sqMat = new Pointer();
//        JCuda.cudaMalloc(sqMat, 4 * Sizeof.DOUBLE);
//        
//        //dfghdrfeghjcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T, CUBLAS_OP_N,2,2,m,device_alpha,device_tmp2,m,device_tmp2,m,device_beta,sqMat,2);
//                           
//        //jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
//        
//        double [] mat = new double [4];
//        Pointer hmat=Pointer.to(mat);
//        JCuda.cudaMemcpyAsync(hmat, sqMat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));
//        
//        double det=(mat[0]*mat[3]-mat[1]*mat[2]);
//        if (det!=0){
//            det=1/det;
//        }
//        else{
//            IJ.log("error inversion matrix det=0");
//            return 0;
//        }
//        double tmp=mat[3];
//        mat[3]=det*mat[0];
//        mat[0]=det*tmp;
//        mat[1]=-det*mat[1];
//        mat[2]=-det*mat[2];
//        //IJ.log("mat inv "+mat[0]+"  "+mat[1]+"  "+mat[2]+"  "+mat[3]+"  size : "+this.nbDataPerImage);
//        JCuda.cudaMemcpyAsync(sqMat,hmat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
//        
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_T,2,m,2,device_alpha,sqMat,2,device_tmp2,m,device_beta,device_tmp2,2);
//                          
//        
//        VecDouble.mulScalar(this.nbDataPerImage, device_tmp, device_zer4, this.getA(4));//4e zernike coef
//        
//        
//        Pointer dbMat = new Pointer();
//        JCuda.cudaMalloc(dbMat, 2 * Sizeof.DOUBLE);
//        
//        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_N,2,1,m,device_alpha,device_tmp2,2,device_tmp,m,device_beta,dbMat,2);
//         
//        
//        
//        double [] res =new double [2];
//        Pointer pres= Pointer.to(res);
//        JCuda.cudaMemcpyAsync(pres, dbMat,  2*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost, MyCudaStream.getCudaStream_t(param.stream));//half first part is 1
////        for (int i=0;i<res.length;i++){
////            IJ.log("res "+res[i]+"  ");
////        }
//        
//        //add or substract to X res[1]
//        if (res[1]!=0)
//            IJ.log("shift Z "+res[1]+"   "+(1/res[1]));
//        else
//            IJ.log("shift Z "+res[1]);
//        JCuda.cudaFree(sqMat);
//        return res[1];
//    }
//    
//    
//    
    
    
    
    
    
    
    
    
    //fit coef to kz (parabola)
    /*public double initCoefZ(Pointer device_kz){
        
        
        /////////to do///
        
        JCuda.cudaMemcpyAsync(device_tmp2.withByteOffset(nbDataPerImage*Sizeof.DOUBLE), device_kz, nbDataPerImage* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice);//half first part is 1
        
        
        
        Pointer sqMat = new Pointer();
        JCuda.cudaMalloc(sqMat, 4 * Sizeof.DOUBLE);
        
        //dfghdrfeghjcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_T, CUBLAS_OP_N,2,2,m,device_alpha,device_tmp2,m,device_tmp2,m,device_beta,sqMat,2);
                           
        //jcuda.jcublas.JCublas2.cublasDgemv(handlecublas,CUBLAS_OP_T,m,2,device_alpha,device_tmp2,m,device_tmp2,incr,device_beta,sqMat,incr);
        
        double [] mat = new double [4];
        Pointer hmat=Pointer.to(mat);
        JCuda.cudaMemcpyAsync(hmat, sqMat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);
        
        double det=(mat[0]*mat[3]-mat[1]*mat[2]);
        if (det!=0){
            det=1/det;
        }
        else{
            IJ.log("error inversion matrix det=0");
            return 0;
        }
        double tmp=mat[3];
        mat[3]=det*mat[0];
        mat[0]=det*tmp;
        mat[1]=-det*mat[1];
        mat[2]=-det*mat[2];
        //IJ.log("mat inv "+mat[0]+"  "+mat[1]+"  "+mat[2]+"  "+mat[3]+"  size : "+this.nbDataPerImage);
        JCuda.cudaMemcpyAsync(sqMat,hmat, 4* Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice);//half first part is 1
        
        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_T,2,m,2,device_alpha,sqMat,2,device_tmp2,m,device_beta,device_tmp2,2);
                          
        
        VecDouble.mulScalar(this.nbDataPerImage, device_tmp, device_zer4, this.getA(4));//4e zernike coef
        
        
        Pointer dbMat = new Pointer();
        JCuda.cudaMalloc(dbMat, 2 * Sizeof.DOUBLE);
        
        jcuda.jcublas.JCublas2.cublasDgemm(handlecublas,CUBLAS_OP_N, CUBLAS_OP_N,2,1,m,device_alpha,device_tmp2,2,device_tmp,m,device_beta,dbMat,2);
         
        
        
        double [] res =new double [2];
        Pointer pres= Pointer.to(res);
        JCuda.cudaMemcpyAsync(pres, dbMat,  2*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);//half first part is 1
//        for (int i=0;i<res.length;i++){
//            IJ.log("res "+res[i]+"  ");
//        }
        
        //add or substract to X res[1]
        if (res[1]!=0)
            IJ.log("shift Z "+res[1]+"   "+(1/res[1]));
        else
            IJ.log("shift Z "+res[1]);
        
        return res[1];
    }*/
    
    
    
    
//    
//    
//    
//    
//    
//    boolean test2(int N)
//    {
//        JCublas2.setExceptionsEnabled(true);
//        JCuda.setExceptionsEnabled(true);
//       
//        double h_A[];
//        double h_B[];
//        double h_C[];
//        Pointer d_A = new Pointer();
//        Pointer d_B = new Pointer();
//        Pointer d_C = new Pointer();
//        Pointer alpha = new Pointer();
//        Pointer beta = new Pointer();
//        int n2 = N * N;
//        int i;
//        cublasHandle handle = new cublasHandle();
// 
//        System.out.printf("CUBLAS test %d running...\n", N);
// 
//        cublasCreate(handle);
//       
//        h_A = new double[n2];
//        h_B = new double[n2];
//        h_C = new double[n2];
//        Random random = new Random(0);
//        for (i = 0; i < n2; i++) {
//            h_A[i] = random.nextDouble();
//            h_B[i] = random.nextDouble();
//            h_C[i] = random.nextDouble();
//        }
// 
//        cudaMalloc(d_A, n2 * Sizeof.DOUBLE);
//        cudaMalloc(d_B, n2 * Sizeof.DOUBLE);
//        cudaMalloc(d_C, n2 * Sizeof.DOUBLE);
//       
//        cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(h_A), 1, d_A, 1);
//        cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(h_B), 1, d_B, 1);
//        cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(h_C), 1, d_C, 1);
// 
//        cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
// 
//        cudaMalloc(alpha, 1 * Sizeof.DOUBLE);
//        cudaMemcpyAsync(alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
// 
//        cudaMalloc(beta, 1 * Sizeof.DOUBLE);
//        cudaMemcpyAsync(beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
// 
//        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, alpha, d_A, N, d_B, N, beta, d_C, N);
//       
//        h_C = new double[n2];
//        cublasGetVector(n2, Sizeof.DOUBLE, d_C, 1, Pointer.to(h_C), 1);
// 
//        boolean result = true;
//        for (i=0; i<N*N; i++)
//        {
//            if (Double.isNaN(h_C[i]))
//            {
//                System.out.printf("Error for size %d at %d\n", N, i);
//                result = false;
//                break;
//            }
//        }
// 
//        cudaFree(d_A);
//        cudaFree(d_B);
//        cudaFree(d_C);
//        cudaFree(alpha);
//        cudaFree(beta);
// 
//        cublasDestroy(handle);
//        return result;
//    }
//    
//    
//    
//    
//    
//    boolean runTest(int N)
//    {
//        JCublas2.setExceptionsEnabled(true);
//        JCuda.setExceptionsEnabled(true);
//       
//        double h_A[];
//        double h_B[];
//        double h_C[];
//        Pointer d_A = new Pointer();
//        Pointer d_B = new Pointer();
//        Pointer d_C = new Pointer();
//        Pointer alpha = new Pointer();
//        Pointer beta = new Pointer();
//        int n2 = N * N;
//        int i;
//        cublasHandle handle = new cublasHandle();
// 
//        System.out.printf("CUBLAS test %d running...\n", N);
// 
//        cublasCreate(handle);
//       
//        h_A = new double[n2];
//        h_B = new double[n2];
//        h_C = new double[n2];
//        Random random = new Random(0);
//        for (i = 0; i < n2; i++) {
//            h_A[i] = random.nextDouble();
//            h_B[i] = random.nextDouble();
//            h_C[i] = random.nextDouble();
//        }
// 
//        cudaMalloc(d_A, n2 * Sizeof.DOUBLE);
//        cudaMalloc(d_B, n2 * Sizeof.DOUBLE);
//        cudaMalloc(d_C, n2 * Sizeof.DOUBLE);
//       
//        cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(h_A), 1, d_A, 1);
//        cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(h_B), 1, d_B, 1);
//        cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(h_C), 1, d_C, 1);
// 
//        cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
// 
//        cudaMalloc(alpha, 1 * Sizeof.DOUBLE);
//        cudaMemcpyAsync(alpha, Pointer.to(new double[]{1.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
// 
//        cudaMalloc(beta, 1 * Sizeof.DOUBLE);
//        cudaMemcpyAsync(beta, Pointer.to(new double[]{0.0}), 1*Sizeof.DOUBLE, cudaMemcpyHostToDevice, MyCudaStream.getCudaStream_t(param.stream));
// 
//        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, alpha, d_A, N, d_B, N, beta, d_C, N);
//       
//        h_C = new double[n2];
//        cublasGetVector(n2, Sizeof.DOUBLE, d_C, 1, Pointer.to(h_C), 1);
// 
//        boolean result = true;
//        for (i=0; i<N*N; i++)
//        {
//            if (Double.isNaN(h_C[i]))
//            {
//                System.out.printf("Error for size %d at %d\n", N, i);
//                result = false;
//                break;
//            }
//        }
// 
//        cudaFree(d_A);
//        cudaFree(d_B);
//        cudaFree(d_C);
//        cudaFree(alpha);
//        cudaFree(beta);
// 
//        cublasDestroy(handle);
//        return result;
//    }
    
}
