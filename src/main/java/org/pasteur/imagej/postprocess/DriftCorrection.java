/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import org.pasteur.imagej.cuda.*;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.utils.GaussianElimination;
import org.pasteur.imagej.utils.SplineFit;
import org.pasteur.imagej.utils.Matrixe;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import java.util.List;


import org.pasteur.imagej.utils.ImageShow;

import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import java.util.ArrayList;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcufft.JCufft;
import jcuda.jcufft.cufftHandle;
import jcuda.jcufft.cufftType;
import static jcuda.jcusparse.JCusparse.cusparseSgthr;
import static jcuda.jcusparse.JCusparse.cusparseSsctr;
import static jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import jcuda.runtime.JCuda;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpyAsync;
import jcuda.runtime.cudaError;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import jcuda.runtime.JCuda;

/**
 *
 * @author benoit
 */
public class DriftCorrection {
    
    
    
    double precisionX;
    double precisionY;
    double precisionZ;
    
    
    
    double [] driftIndex;
    double [] driftX;
    double [] driftY;
    double [] driftZ;
    
    cufftHandle plan;
    
    Pointer device_input;
    Pointer device_step;
    Pointer device_gaussian;
    
    Pointer host_input;
    Pointer host_step;
    
    float [] input;
    float [] step;
    float [] gaussian;
    float [][][][] result;
    
    
    int [][] framePosition = new int[2][];
    
    StackLocalization stackloc;
    
    int bin;
    
    int sizeFFT=512;
    int sizeFFTPadded;
    
    
    int sizeFFTx=512;
    int sizeFFTPaddedx;
    int sizeFFTy=512;
    int sizeFFTPaddedy;
    int sizeFFTz=512;
    int sizeFFTPaddedz;
    
    double [][][][] gridPosition;//X;Y;Z;debut,fin
    
    int nbX,nbY,nbZ;
    
    
    int cudaResult;
    //In this class, the rendered 3D image is splitted into many 3D images of size sizeFFT to deal with low memory on GPU
    
    
    double minX,maxX,minY,maxY,minZ,maxZ;

    int minF,maxF;
               
    int width,height,depth;
            
    double pixelsizeNM;
    
    double initDriftX=0;
    double initDriftY=0;
    double initDriftZ=0;
    
    //WARNING : [x][y][z] instead of [z][x][y]
    public DriftCorrection(StackLocalization stackloc,double pixelsizeNM,int bin, double subImageSizeUm,String path){
        
        
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
            catch(Exception eeee){IJ.log("WARNING: file to attach not found");}
            if ((initDriftX==0)&&(initDriftX==0)&&(initDriftX==0)){
                IJ.log("WARNING: file attached not found or not drift corrected");
            }
        }
        
        
        
        
        
        driftIndex=new double[bin];
        driftX=new double[bin];
        driftY=new double[bin];
        driftZ=new double[bin];
        
        
        /*long [] memtot=new long[1];
        long [] memfree=new long[1];
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        
        //IJ.log("mem available "+memfree[0]+" / "+memtot[0]);
       
        
        long neededMemory=((long)(3*2))*((long)sizeFFT*(long)sizeFFT*(long)sizeFFT*(long)8*(long)2 * (long)Sizeof.FLOAT);//+200 mo free
        //IJ.log("need memory   "+neededMemory+"  "+3*(sizeFFT*sizeFFT*sizeFFT*8*2 * Sizeof.FLOAT));
        
        while (neededMemory>memfree[0]){
            sizeFFT*=.9;
            if (sizeFFT%2==1){
                sizeFFT++;
            }
            neededMemory=((long)(3*2))*((long)sizeFFT*(long)sizeFFT*(long)sizeFFT*(long)8*(long)2 * (long)Sizeof.FLOAT)+(long)100000000;//+200 mo free
        }
        
        IJ.log("size FFT according to available memory: "+sizeFFT+"  need mem:"+neededMemory);*/
        
        
        
        
        this.bin=bin;
        
        this.stackloc=stackloc;
        
        minX=Double.MAX_VALUE;
        maxX=Double.NEGATIVE_INFINITY;
        
        
        minY=Double.MAX_VALUE;
        maxY=Double.NEGATIVE_INFINITY;
        
        minZ=Double.MAX_VALUE;
        maxZ=Double.NEGATIVE_INFINITY;
        
        minF=Integer.MAX_VALUE;
        maxF=Integer.MIN_VALUE;
        
        
        
        double x;
        double y;
        double z;
        int f;
        framePosition = new int[stackloc.fl.size()][2];
        
        precisionX=0;
        precisionY=0;
        precisionZ=0;
        int countPrecX=0;
        int countPrecY=0;
        int countPrecZ=0;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<stackloc.fl.size();i++){
            framePosition[i][0]=stackloc.fl.get(i).numFrame;
            framePosition[i][1]=i;
            for (int j=0;j<stackloc.fl.get(i).loc.size();j++){
                if (stackloc.fl.get(i).loc.get(j).exists){
                    
                    x=stackloc.fl.get(i).loc.get(j).X;
                    y=stackloc.fl.get(i).loc.get(j).Y;
                    z=stackloc.fl.get(i).loc.get(j).Z;
                    f=stackloc.fl.get(i).loc.get(j).frame;
                    
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                    if (f<minF){
                        minF=f;
                    }
                    if (f>maxF){
                        maxF=f;
                    }
                    double crlbx=stackloc.fl.get(i).loc.get(j).crlb_X;
                    double crlby=stackloc.fl.get(i).loc.get(j).crlb_Y;
                    double crlbz=stackloc.fl.get(i).loc.get(j).crlb_Z;
                    if (crlbx!=-1){
                        precisionX+=crlbx;
                        countPrecX++;
                    }
                    if (crlby!=-1){
                        precisionY+=crlby;
                        countPrecY++;
                    }
                    if (crlbz!=-1){
                        precisionZ+=crlbz;
                        countPrecZ++;
                    }
                    
                }
            }
        }
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        
        
        if (countPrecX==0){
            //IJ.log("no CRLB available, the precision threshold is set to 50 nm");
            precisionX=50;//default 50 nm
        }
        else{
            precisionX/=countPrecX;
        }
        if (countPrecY==0){
            precisionY=50;
        }
        else{
            precisionY/=countPrecY;
        }
        if (countPrecZ==0){
            precisionZ=50;
        }
        else{
            precisionZ/=countPrecZ;
        }
        
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        
        double nbMem=Double.POSITIVE_INFINITY;
        
        
        
        long [] memtot=new long[1];
        long [] memfree=new long[1];
        
        JCuda.cudaMemGetInfo(memfree, memtot);
        //IJ.log("needed GPU memory >"+Math.round(nbMem/1000000.)+" MB");
        //IJ.log("mem available  "+Math.round(memfree[0]/1000000)+" / "+Math.round(memtot[0]/1000000)+"  MB");
        
        
        
        int nbloop=0;
        while (nbMem>=memfree[0]){
            
            
            if (nbloop!=0){
                pixelsizeNM+=5;
            }
            
            sizeFFT=(int)(subImageSizeUm/(pixelsizeNM/1000.));
            
            sizeFFTPadded=sizeFFT*2;
            
            

            width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
            height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
            depth=(int)Math.ceil((maxZ-minZ)/pixelsizeNM);
            width=Math.max(width, 1);
            height=Math.max(height, 1);
            depth=Math.max(depth, 1);

            sizeFFTx=Math.min(sizeFFT, width);
            sizeFFTy=Math.min(sizeFFT, height);
            sizeFFTz=Math.min(sizeFFT, depth);

            //sizeFFTx+=sizeFFTx%2;//becomes even
            //sizeFFTy+=sizeFFTy%2;//becomes even
            //sizeFFTz+=sizeFFTz%2;//becomes even


            
            sizeFFTPaddedx=sizeFFTx*2;
            sizeFFTPaddedy=sizeFFTy*2;
            sizeFFTPaddedz=sizeFFTz*2;
            
            nbMem=3.5*((double)sizeFFTPaddedx*(double)sizeFFTPaddedy*(double)sizeFFTPaddedz*2. * (double)Sizeof.FLOAT*3);
            nbloop++;
            if (nbloop==20){
                //drift not possible... wrong parameters ???
                break;
            }
        }
        
        this.pixelsizeNM=pixelsizeNM;
        
        
        
        if(nbloop>1){
            IJ.log("Pixel size increased to "+pixelsizeNM+" nm to deal with available GPU memory");
            IJ.log("You can also try do reduce the maximum_drift parameter if you need better precision");
        }
        
//        IJ.log("width "+width+"  "+sizeFFTx);
//        IJ.log("width "+depth+"  "+sizeFFTz);
//        
//        nbX=(int)Math.ceil((double)(width-this.sizeFFTx/2)/(double)(this.sizeFFTx/2));//   /2 car chevauchement
//        nbY=(int)Math.ceil((double)(height-this.sizeFFTy/2)/(double)(this.sizeFFTy/2));//   /2 car chevauchement
//        nbZ=(int)Math.ceil((double)(depth-this.sizeFFTz/2)/(double)(this.sizeFFTz/2));//   /2 car chevauchement
//        
//        IJ.log("nbXYZ "+nbX+"  "+nbY+"  "+nbZ);
        
        nbX=(int)(Math.ceil(2.*(width/this.sizeFFTx))-1);
        nbY=(int)(Math.ceil(2.*(height/this.sizeFFTy))-1);
        nbZ=(int)(Math.ceil(2.*(depth/this.sizeFFTz))-1);

//        IJ.log("nbXYZ "+nbX+"  "+nbY+"  "+nbZ);
        
        gridPosition=new double[nbX][nbY][nbZ][6];

        for (int u=0;u<nbX;u++){
            for (int uu=0;uu<nbY;uu++){
                for (int uuu=0;uuu<nbZ;uuu++){
                    //avec chevauchement de taille (sizeFFT/2)
                    gridPosition[u][uu][uuu][0]=(u/2)*sizeFFTx*pixelsizeNM+(u%2)*(sizeFFTx/2)*pixelsizeNM;
                    gridPosition[u][uu][uuu][1]=(u/2+1)*sizeFFTx*pixelsizeNM+(u%2)*(sizeFFTx/2)*pixelsizeNM;
                    gridPosition[u][uu][uuu][2]=(uu/2)*sizeFFTy*pixelsizeNM+(uu%2)*(sizeFFTy/2)*pixelsizeNM;
                    gridPosition[u][uu][uuu][3]=(uu/2+1)*sizeFFTy*pixelsizeNM+(uu%2)*(sizeFFTy/2)*pixelsizeNM;
                    gridPosition[u][uu][uuu][4]=(uuu/2)*sizeFFTz*pixelsizeNM+(uuu%2)*(sizeFFTz/2)*pixelsizeNM;
                    gridPosition[u][uu][uuu][5]=(uuu/2+1)*sizeFFTz*pixelsizeNM+(uuu%2)*(sizeFFTz/2)*pixelsizeNM;
                    //IJ.log("grid x "+gridPosition[u][uu][uuu][0]+"   "+gridPosition[u][uu][uuu][1]+"  "+width*pixelsizeNM);
                    //IJ.log("grid y "+gridPosition[u][uu][uuu][2]+"   "+gridPosition[u][uu][uuu][3]);
                    //IJ.log("grid z "+gridPosition[u][uu][uuu][4]+"   "+gridPosition[u][uu][uuu][5]);
                    //IJ.log("");
                }
            }
        }
        
        
        
        
        
        Arrays.sort(framePosition, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                return ((Integer) o1[0]).compareTo(o2[0]);
            }
        });
        
        int xx;
        int yy;
        int zz;
        
        
        
        input=new float[sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2];
        
        
        
        step=new float[sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2];
        
        gaussian=new float[sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2];
        
        
        result=new float[1][sizeFFTz][sizeFFTy][sizeFFTx];
        
        
        
        if (nbMem<memfree[0]){
            device_input=new Pointer();
            cudaResult =cudaMalloc(device_input, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2 * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 1 "+cudaResult);}
            
            
            
            
            device_step=new Pointer();
            cudaResult =cudaMalloc(device_step, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2 * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 2 "+cudaResult);}


            device_gaussian=new Pointer();
            cudaResult =cudaMalloc(device_gaussian, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2 * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 2 "+cudaResult);}
            double sigmaNM=3*Math.max(Math.max(precisionX, precisionY), precisionZ);
            double sigma=sigmaNM/pixelsizeNM;
            double sigmaFourierPowx=((1./(sigma))*(double)sizeFFTPaddedx/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTPaddedx/(Math.PI*2.));

            double sigmaFourierPowy=((1./(sigma))*(double)sizeFFTPaddedy/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTPaddedy/(Math.PI*2.));

            double sigmaFourierPowz=((1./(sigma))*(double)sizeFFTPaddedz/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTPaddedz/(Math.PI*2.));



            int xxshift;
            int yyshift;
            int zzshift;

            double dx,dy,dz;
            for (int xxx=0;xxx<sizeFFTPaddedx;xxx++){
                    for (int yyy=0;yyy<sizeFFTPaddedy;yyy++){
                        for (int zzz=0;zzz<sizeFFTPaddedz;zzz++){
                        xxshift=((xxx+sizeFFTPaddedx/2)%sizeFFTPaddedx);
                        yyshift=((yyy+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                        zzshift=((zzz+sizeFFTPaddedz/2)%sizeFFTPaddedz);

                        dx=((double)xxx-(double)sizeFFTPaddedx/2.)*((double)xxx-(double)sizeFFTPaddedx/2.);
                        dy=((double)yyy-(double)sizeFFTPaddedy/2.)*((double)yyy-(double)sizeFFTPaddedy/2.);
                        dz=((double)zzz-(double)sizeFFTPaddedz/2.)*((double)zzz-(double)sizeFFTPaddedz/2.);

                        gaussian[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2]=(float)Math.exp(-.5*(dx/sigmaFourierPowx+dy/sigmaFourierPowy+dz/sigmaFourierPowz));
                        gaussian[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2+1]=(float)Math.exp(-.5*(dx/sigmaFourierPowx+dy/sigmaFourierPowy+dz/sigmaFourierPowz));


                        //inputTMP[zzshift][xxshift][yyshift]=(float)Math.exp(-.5*(dx+dy+dz)/sigmaFourierPow);
                    }
                }
            }

            //ImageShow.imshow(inputTMP,"gauss");

            cudaResult =cudaMemcpyAsync(device_gaussian, Pointer.to(gaussian), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}



            plan = new cufftHandle();
            //IJ.log("WARNING order padd x y z may be WRONG");
            cudaResult =JCufft.cufftPlan3d(plan, sizeFFTPaddedz,sizeFFTPaddedy,sizeFFTPaddedx, cufftType.CUFFT_C2C);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 11 "+cudaResult);}

            cudaResult =JCufft.cufftSetStream(plan, MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 12 "+cudaResult);}



            /*

            int id=0;
            int posit;

            loopinput:for (;true;id++){

                posit=framePosition[id][1];

                int idBin=bin*(stackloc.fl.get(posit).numFrame-minF)/(1+maxF-minF);

                if (idBin>0){
                    break loopinput;
                }

                for (int j=0;j<stackloc.fl.get(posit).loc.size();j++){
                    if (stackloc.fl.get(posit).loc.get(j).exists){
                        x=stackloc.fl.get(posit).loc.get(j).X;
                        y=stackloc.fl.get(posit).loc.get(j).Y;
                        z=stackloc.fl.get(posit).loc.get(j).Z;

                        xx=Math.min(Math.max((int)((x-minX)/pixelsizeNM),0),width-1);
                        yy=Math.min(Math.max((int)((y-minY)/pixelsizeNM),0),height-1);
                        zz=Math.min(Math.max((int)((z-minZ)/pixelsizeNM),0),depth-1);


                        if (idBin==0){
                            if (input[zz*width*height+yy*width+xx]<Short.MAX_VALUE){
                                input[zz*width*height+yy*width+xx]++;
                            }
                        }


                    }
                }
            }


            //ImageShow.imshow(input,"input");

            for (int s=1;s<bin;s++){


                loopstep:for (;true;id++){

                    posit=framePosition[id][1];

                    int idBin=bin*(stackloc.fl.get(posit).numFrame-minF)/(1+maxF-minF);

                    if (idBin>s){
                        break loopstep;
                    }

                    for (int j=0;j<stackloc.fl.get(posit).loc.size();j++){
                        if (stackloc.fl.get(posit).loc.get(j).exists){
                            x=stackloc.fl.get(posit).loc.get(j).X;
                            y=stackloc.fl.get(posit).loc.get(j).Y;
                            z=stackloc.fl.get(posit).loc.get(j).Z;

                            xx=Math.min(Math.max((int)((x-minX)/pixelsizeNM),0),width-1);
                            yy=Math.min(Math.max((int)((y-minY)/pixelsizeNM),0),height-1);
                            zz=Math.min(Math.max((int)((z-minZ)/pixelsizeNM),0),depth-1);


                            if (idBin==s){
                                if (step[zz*width*height+yy*width+xx]<Short.MAX_VALUE){
                                    step[zz*width*height+yy*width+xx]++;
                                }
                            }


                        }
                    }
                }


                //do cross correl here





                step=null;
            }
            */
            runFast();
            //run();
            free();
        }
        else{
            IJ.log("DRIFT CORRECTION IMPOSSIBLE, not enough memory in GPU.");
            if (subImageSizeUm>5){
                IJ.log("If possible, try to increase the rendering pixel size or reduce the maximum drift bound (~ 5 µm)");
            }
            else{
                IJ.log("If possible, try to increase the rendering pixel size or reduce the maximum drift bound");
            }
        }
        
    }
    
    
    
    
    
    
    void run(){
        
        long time0=System.currentTimeMillis();
        double x;
        double y;
        double z;
        
        
        
        int xx,xxx;
        int yy,yyy;
        int zz,zzz;
        
        int xxshift;
        int yyshift;
        int zzshift;
        
        //int xxshiftEven;
        //int yyshiftEven;
        //int zzshiftEven;
        
        int posit;
        
        
        long t1=System.currentTimeMillis();
        
        
        
        //step 0 : no drift
        double frameNumberPerBin=((double)(maxF-minF))/((double)(bin));
        driftIndex[0]=((double)minF)+(frameNumberPerBin/2.);//mean of frame number
        driftX[0]=0;
        driftY[0]=0;
        driftZ[0]=0;
        
        
        
        //for each step
        int passNumber=bin-1;
        //int passNumber=1;
        int numbId=0;
        
        int totNumb=0;
        for (int times=0;times<passNumber;times++){//times is the reference for cross correlation ... 0 at first... 1 after... etc
            for (int s=times+1;s<bin;s++){
                totNumb++;
            }
        }
        
        double [][] A=new double[totNumb][bin-1];
        for (int i=0;i<A.length;i++){
            for (int ii=0;ii<A[0].length;ii++){
                A[i][ii]=0;
            }
        }
        double [][] R=new double[3][totNumb];
        
        
        for (int times=0;times<passNumber;times++){//times is the reference for cross correlation ... 0 at first... 1 after... etc
            for (int s=times+1;s<bin;s++){


                for (xxx=0;xxx<sizeFFTz;xxx++){
                    for (yyy=0;yyy<sizeFFTy;yyy++){
                        for (zzz=0;zzz<sizeFFTx;zzz++){
                            result[0][xxx][yyy][zzz]=0;
                        }
                    }
                }
boolean mybool=true;
if (mybool){
                loopi:for (int u=0;u<nbX;u++){

                    IJ.showProgress((float)numbId/(float)totNumb);

                    for (int uu=0;uu<nbY;uu++){
                        for (int uuu=0;uuu<nbZ;uuu++){
                            //int [][][] inputTMP=new int[sizeFFTPadded][sizeFFTPadded][sizeFFTPadded];
                            for (int i=0;i<input.length;i++){
                                input[i]=0;
                                step[i]=0;
                            }

                            int id=0;

                            int count=0;

                            
                            //manage frame shift when times!=0
                            for (int up=0;up<times;up++){
                                loopinput0:for (;id<framePosition.length;id++){

                                    posit=framePosition[id][1];

                                    int idBin=bin*(stackloc.fl.get(posit).numFrame-minF)/(1+maxF-minF);
                                    //IJ.log("posit "+posit+"   "+stackloc.fl.get(posit).numFrame+"   "+idBin+"    "+bin);
                                    if (idBin>=times){
                                        break loopinput0;
                                    }
                                }
                            }
                            
                            
                            
                            
                            
                            //copy input 
                            loopinput:for (;id<framePosition.length;id++){

                                posit=framePosition[id][1];

                                int idBin=bin*(stackloc.fl.get(posit).numFrame-minF)/(1+maxF-minF);
                                //IJ.log("posit "+posit+"   "+stackloc.fl.get(posit).numFrame+"   "+idBin+"    "+bin);
                                if (idBin>times){
                                    break loopinput;
                                }

                                for (int j=0;j<stackloc.fl.get(posit).loc.size();j++){
                                    if (stackloc.fl.get(posit).loc.get(j).exists){
                                        x=stackloc.fl.get(posit).loc.get(j).X-minX;
                                        y=stackloc.fl.get(posit).loc.get(j).Y-minY;
                                        z=stackloc.fl.get(posit).loc.get(j).Z-minZ;




                                        if (idBin==times){


                                            if ((x>=this.gridPosition[u][uu][uuu][0])&&(x<this.gridPosition[u][uu][uuu][1])&&(y>=this.gridPosition[u][uu][uuu][2])&&(y<this.gridPosition[u][uu][uuu][3])&&(z>=this.gridPosition[u][uu][uuu][4])&&(z<this.gridPosition[u][uu][uuu][5])){

                                                //make histogram padded:
                                                xx=(sizeFFTx/2)+Math.min(Math.max((int)(((x)%(sizeFFTx*pixelsizeNM))/pixelsizeNM),0),sizeFFTx-1);
                                                yy=(sizeFFTy/2)+Math.min(Math.max((int)(((y)%(sizeFFTy*pixelsizeNM))/pixelsizeNM),0),sizeFFTy-1);
                                                zz=(sizeFFTz/2)+Math.min(Math.max((int)(((z)%(sizeFFTz*pixelsizeNM))/pixelsizeNM),0),sizeFFTz-1);

                                                //xxshiftEven=((xx+sizeFFTPadded/2)%sizeFFTPadded)*2;//*2 because even
                                                //yyshiftEven=((yy+sizeFFTPadded/2)%sizeFFTPadded)*2;
                                                //zzshiftEven=((zz+sizeFFTPadded/2)%sizeFFTPadded)*2;

                                                xxshift=((xx+sizeFFTPaddedx/2)%sizeFFTPaddedx);//*2 because even
                                                yyshift=((yy+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                                                zzshift=((zz+sizeFFTPaddedz/2)%sizeFFTPaddedz);


                                                if (input[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2]<Float.MAX_VALUE){
                                                    input[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2]++;
                                                    //inputTMP[zz][xx][yy]=zz;
                                                    count++;
                                                }
                                            }


                                        }


                                    }
                                }
                            }

                            //ImageShow.imshow(inputTMP,"input");

                            cudaResult =cudaMemcpyAsync(device_input, Pointer.to(input), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}

                            cudaResult =JCufft.cufftExecC2C(plan, device_input,device_input,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("drift ERROR cuda fft1 "+cudaResult);}



                            
                            
                            

                            //inputTMP=new int[sizeFFTPadded][sizeFFTPadded][sizeFFTPadded];

                            loopstep:for (;id<framePosition.length;id++){

                                posit=framePosition[id][1];

                                int idBin=bin*(stackloc.fl.get(posit).numFrame-minF)/(1+maxF-minF);

                                if (idBin>s){
                                    break loopstep;
                                }

                                for (int j=0;j<stackloc.fl.get(posit).loc.size();j++){
                                    if (stackloc.fl.get(posit).loc.get(j).exists){
                                        x=stackloc.fl.get(posit).loc.get(j).X-minX;
                                        y=stackloc.fl.get(posit).loc.get(j).Y-minY;
                                        z=stackloc.fl.get(posit).loc.get(j).Z-minZ;

                                        if (idBin==s){
                                            
                                            
                                            if ((x>=this.gridPosition[u][uu][uuu][0])&&(x<this.gridPosition[u][uu][uuu][1])&&(y>=this.gridPosition[u][uu][uuu][2])&&(y<this.gridPosition[u][uu][uuu][3])&&(z>=this.gridPosition[u][uu][uuu][4])&&(z<this.gridPosition[u][uu][uuu][5])){

                                                //make histogram padded:
                                                xx=(sizeFFTx/2)+Math.min(Math.max((int)(((x)%(sizeFFTx*pixelsizeNM))/pixelsizeNM),0),sizeFFTx-1);
                                                yy=(sizeFFTy/2)+Math.min(Math.max((int)(((y)%(sizeFFTy*pixelsizeNM))/pixelsizeNM),0),sizeFFTy-1);
                                                zz=(sizeFFTz/2)+Math.min(Math.max((int)(((z)%(sizeFFTz*pixelsizeNM))/pixelsizeNM),0),sizeFFTz-1);

                                                //xxshiftEven=((xx+sizeFFTPadded/2)%sizeFFTPadded)*2;//*2 because even
                                                //yyshiftEven=((yy+sizeFFTPadded/2)%sizeFFTPadded)*2;
                                                //zzshiftEven=((zz+sizeFFTPadded/2)%sizeFFTPadded)*2;

                                                xxshift=((xx+sizeFFTPaddedx/2)%sizeFFTPaddedx);//*2 because even
                                                yyshift=((yy+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                                                zzshift=((zz+sizeFFTPaddedz/2)%sizeFFTPaddedz);


                                                if (step[(zzshift*sizeFFTPaddedy*sizeFFTPaddedx+yyshift*sizeFFTPaddedx+xxshift)*2]<Float.MAX_VALUE){
                                                    step[(zzshift*sizeFFTPaddedy*sizeFFTPaddedx+yyshift*sizeFFTPaddedx+xxshift)*2]++;
                                                    //inputTMP[zz][xx][yy]=zz;
                                                    count++;
                                                }
                                            }


                                        }


                                    }
                                }
                            }

                            //ImageShow.imshow(inputTMP,"step"+s);

                            cudaResult =cudaMemcpyAsync(device_step, Pointer.to(step), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
                            cudaResult =JCufft.cufftExecC2C(plan, device_step,device_step,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("drift ERROR cuda fft1 "+cudaResult);}




                            //do cross correl here




                            MyVecDouble.complexeConjugateKernel(MyCudaStream.getCUstream(0), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz, device_input, device_input, device_step);

                            MyVecDouble.mul_fl(MyCudaStream.getCUstream(0), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2, device_input, device_input, device_gaussian);



                            JCufft.cufftExecC2C(plan, device_input,device_input,JCufft.CUFFT_INVERSE);

                            //could be removed to win few time
                            MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(0), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2, device_input, device_input, (float)Math.sqrt(sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz));

                            cudaResult =cudaMemcpyAsync(Pointer.to(input), device_input, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(0));



                            for (xxx=0;xxx<sizeFFTx;xxx++){
                                for (yyy=0;yyy<sizeFFTy;yyy++){
                                    for ( zzz=0;zzz<sizeFFTz;zzz++){
                                        //xxshiftEven=((xxx+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;//*2 because even
                                        //yyshiftEven=((yyy+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;
                                        //zzshiftEven=((zzz+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;

                                        xxshift=((xxx+(sizeFFTx/2)+sizeFFTPaddedx/2)%sizeFFTPaddedx);//*2 because even
                                        yyshift=((yyy+(sizeFFTy/2)+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                                        zzshift=((zzz+(sizeFFTz/2)+sizeFFTPaddedz/2)%sizeFFTPaddedz);

                                        result[0][zzz][yyy][xxx]+=input[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2];
                                    }
                                }
                            }




                            //break loopi;

                        }
                    }
                }


                /*
                int sizesub=50;
                float [][][] subResult=new float [sizesub][sizesub][sizesub];
                for (int i=0;i<sizesub;i++){
                    for (int ii=0;ii<sizesub;ii++){
                        for (int iii=0;iii<sizesub;iii++){
                            subResult[iii][i][ii]=result[0][i+sizeFFT/2-sizesub/2][ii+sizeFFT/2-sizesub/2][iii+sizeFFT/2-sizesub/2];
                        }
                    }
                }*/
                //ImageShow.imshow(subResult,"subResult ");
                

                //ImageShow.imshow(result,"result "+s);

}

                LocalMaxima locmax = new LocalMaxima(result[0]);
                locmax.run(.95);
                if (times==0){
                    driftIndex[s]=(((double)s)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.);//mean of frame number
                    if (mybool){
                    driftX[s]=locmax.getDriftXinPix()*pixelsizeNM;
                    driftY[s]=locmax.getDriftYinPix()*pixelsizeNM;
                    driftZ[s]=locmax.getDriftZinPix()*pixelsizeNM;
                    }
//                    IJ.log("drift X "+driftX[numbId]+"  nm");
//                    IJ.log("drift Y "+driftY[numbId]+"  nm");
//                    IJ.log("drift Z "+driftZ[numbId]+"  nm");
//                    
//                    IJ.log("");
                }
                
                R[0][numbId]=locmax.getDriftXinPix()*pixelsizeNM;
                R[1][numbId]=locmax.getDriftYinPix()*pixelsizeNM;
                R[2][numbId]=locmax.getDriftZinPix()*pixelsizeNM;
                
                
                
                //A[numbId][times]=-((((double)times)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.));
                
                //A[numbId][s]=((((double)s)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.));
                if (times!=0){
                    A[numbId][times-1]=-1;
                }
                
                A[numbId][s-1]=1;
                
                
                numbId++;
            }
        }
        
        
        //double [][] AA={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0}};;
        //double [][] RR={{-47.77584571501592,-106.8553677959379,-0.7675404540464115},{-133.4594686872343,-132.7645498087783,-4.543981918508422},{-182.92575929587613,-101.1243637780396,-22.19715951653427},{-186.17080657785223,-82.66789094763425,-38.72344841612474},{-75.51451595803726,-16.538492083225975,-7.230904698883478},{-132.69461898940094,6.694369064555872,-19.976627901110078},{-130.36682727684124,34.06332981424143,-37.18433264789738},{-44.83839636955267,26.70459137640364,-10.569133284470666},{-56.44697715252249,48.07490244759762,-29.148866554700703},{-6.520611747652794,23.982255741189817,-15.682598027674999}};
        
        
        
        
        double [] threshold = new double [3];
        
        
        threshold[0]=precisionX;
        threshold[1]=precisionY;
        threshold[2]=precisionZ;
        
        //IJ.log("threshold X:"+threshold[0]+"  Y:"+threshold[1]+"  Z:"+threshold[2]);
                
//        IJ.log("RR :");
//        for (int u=0;u<totNumb;u++){
            /*R[0][u]=RR[u][0];
            R[1][u]=RR[u][1];
            R[2][u]=RR[u][2];*/
//            IJ.log(""+R[0][u]+"  "+R[1][u]+"  "+R[2][u]+"  ");
//        }
//        IJ.log("");
//        IJ.log("AA");
//        for (int u=0;u<A.length;u++){
//            String s="";
//            for (int uu=0;uu<A[0].length;uu++){
//                //A[u][uu]=AA[u][uu];
//                s+=A[u][uu]+" ";
//            }
//            IJ.log(""+s);
//        }
        
//        IJ.log("drift before :");
//        for (int u=0;u<driftZ.length;u++){
//            IJ.log(""+driftX[u]+"  "+driftY[u]+"  "+driftZ[u]+"  ");
//        }
//        IJ.log("");
        
        
        if ((bin>2)&&(A.length>bin)){

            ArrayList<Integer> [] toDelete = new ArrayList[3];
            for (int i=0;i<3;i++){
                toDelete[i]= new ArrayList<Integer>();
            }





            
            for (int check=0;check<2;check++){//2 passes to remove outliers with precision < threshold
                
                
                
                double [][][] Anew=new double[3][][];
                double [][] Rnew=new double [3][];

                for (int dim=0;dim<3;dim++){
                    
                    
//                    IJ.log("dim "+dim);
                    int thelength=totNumb-toDelete[dim].size();
                    
                    Anew[dim]=new double[thelength][bin-1];
                    
                    if (check==1){
//                        IJ.log("reject number="+toDelete[dim].size());
                        
                        //toDelete[dim].clear();
                        //IJ.log("WARNING WARNING WARNING WARNING: no rejection instead reject number replaced by 0");
                    }

                    Rnew[dim]=new double [thelength];


                    for (int i=0,it=0;i<A.length;i++){
                        boolean reject=false;
                        for (int u=0;u<toDelete[dim].size();u++){
                            if (i==(int)toDelete[dim].get(u)){
                                reject=true;
                            }
                        }
                        if (!reject){
                            for (int ii=0;ii<A[0].length;ii++){
                                Anew[dim][it][ii]=A[i][ii];
                            }
                            Rnew[dim][it]=R[dim][i];
                            it++;
                        }
                    }
                    
                    //IJ.log("todelok ");
                    double [] residual = new double [thelength];


                    if (thelength>=bin-1){
                        try{
                            

                            Matrixe m = new Matrixe(Anew[dim]);
                            Matrixe r = new Matrixe(Rnew[dim]);
                            Matrixe mt = Matrixe.transpose(m);
                            Matrixe mm = mt.times(m);

                            Matrixe minv = Matrixe.inverse(mm);
                            Matrixe M = minv.times(mt);


                            Matrixe d = M.times(r);

                            Matrixe res = m.times(d);

                            double [][] res_=res.getMatrixe();

                            for (int i=0;i<residual.length;i++){
                                
                                
                                residual[i]= Math.abs(res_[i][0]-Rnew[dim][i]);

                                if (residual[i]>threshold[dim]){
                                    toDelete[dim].add(i);
                                }
//                                IJ.log("residualOld ("+dim+") "+residual[i]);
                            }
                            double [][] dd=d.getMatrixe();
                            for (int u=0;u<dd.length;u++){
                                
                                if (dim==0){
                                    driftX[u+1]=dd[u][0];
                                }
                                else if (dim==1){
                                    driftY[u+1]=dd[u][0];
                                }
                                else if (dim==2){
                                    driftZ[u+1]=dd[u][0];
                                }
                            }
                            
                            
                            /*
                            //IJ.log("solve... ");
                            double [] D=MatrixColt.solve(Anew[dim], Rnew[dim]);
                            //IJ.log("solve ok ");
                            
                            
                            
                            double [] Ad=new double[Rnew[dim].length];
                            for (int i=0;i<residual.length;i++){
                                Ad[i]=0;
                                for (int ii=0;ii<Anew[dim][i].length;ii++){
                                    Ad[i]+=Anew[dim][i][ii]*D[ii];
                                }
                                residual[i]= Math.abs(Ad[i]-Rnew[dim][i]);

                                if (residual[i]>threshold[dim]){
                                    toDelete[dim].add(i);
                                }
                                //IJ.log("residualNew ("+dim+") "+residual[i]);
                            }
                            
                            
                            for (int u=0;u<D.length;u++){
                                
                                if (dim==0){
                                    driftX[u+1]=D[u];
                                }
                                else if (dim==1){
                                    driftY[u+1]=D[u];
                                }
                                else if (dim==2){
                                    driftZ[u+1]=D[u];
                                }
                            }*/



                        }
                        catch(Exception eee){IJ.log("oops matrix non inversible -> maybe result will not be super precise");}
                    }



                }

//                IJ.log("drift after :");
//                for (int u=0;u<driftZ.length;u++){
//                    IJ.log(""+driftX[u]+"  "+driftY[u]+"  "+driftZ[u]+"  ");
//                }
//                IJ.log("");


            }
            
            
        }

        
        
        
        
        
        //plot(driftIndex,driftX,"drift of X","frame number","drift (nm)");
        
        //plot(driftIndex,driftY,"drift of Y","frame number","drift (nm)");
        
        //plot(driftIndex,driftZ,"drift of Z","frame number","drift (nm)");
        
        
        
        
        
        double [][] di=new double [1][bin];
        double [][] dx=new double [1][bin];
        double [][] dy=new double [1][bin];
        double [][] dz=new double [1][bin];
        
        
        
        for (int i=0;i<bin;i++){
            dx[0][i]=driftX[i];
            dy[0][i]=driftY[i];
            dz[0][i]=driftZ[i];
            di[0][i]=driftIndex[i];
            
        }
        
        SplineFit pfX = new SplineFit(dx,di);
        pfX.run();
        SplineFit pfY = new SplineFit(dy,di);
        pfY.run();
        SplineFit pfZ = new SplineFit(dz,di);
        pfZ.run();
        
        double [] li=new double[maxF-minF];
        double [] lx=new double[maxF-minF];
        double [] ly=new double[maxF-minF];
        double [] lz=new double[maxF-minF];
        
        double [] v=new double[1];//for poly fit result
        double [] vv=new double[1];//for poly fit result
        for (int i=0;i<maxF-minF;i++){
            li[i]=i+minF;
            v[0]=li[i];
            vv=pfX.transform(v);
            lx[i]=vv[0];
            vv=pfY.transform(v);
            ly[i]=vv[0];
            vv=pfZ.transform(v);
            lz[i]=vv[0];
        }
        
        
        
        for (int i=0;i<stackloc.fl.size();i++){
            for (int ii=0;ii<stackloc.fl.get(i).loc.size();ii++){
                v[0]=stackloc.fl.get(i).loc.get(ii).frame;
                vv=pfX.transform(v);
                stackloc.fl.get(i).loc.get(ii).setDrift_X(vv[0]);
                vv=pfY.transform(v);
                stackloc.fl.get(i).loc.get(ii).setDrift_Y(vv[0]);
                vv=pfZ.transform(v);
                stackloc.fl.get(i).loc.get(ii).setDrift_Z(vv[0]);
            }
        }
        
        
        plot(driftIndex,driftX,li,lx,"X_drift","frame number","X drift (nm)");
        
        plot(driftIndex,driftY,li,ly,"Y_drift","frame number","Y drift (nm)");
        
        plot(driftIndex,driftZ,li,lz,"Z_drift","frame number","Z drift (nm)");
        
        
        long t2=System.currentTimeMillis();

        //IJ.log("t2-t1 "+(t2-t1)); 
        
        
        int minut=(int)(3600000.*(((double)(System.currentTimeMillis()-time0)/3600000.)-(double)((System.currentTimeMillis()-time0)/3600000)))/60000;
        IJ.log("elapsed time "+((System.currentTimeMillis()-time0)/3600000)+" h "+minut+" m");
        IJ.showProgress(0);
        IJ.log("Drift correction finished");
        //IJ.log("");
        //IJ.log("In case of extremely large drift, a residual drift can happend using cross-correlation due to .");
        //IJ.log("Please, apply drift correction a second time if you think the correction is not perfect");
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    void runFast(){
        
        long time0=System.currentTimeMillis();
        double x;
        double y;
        double z;
        
        
        
        int xx,xxx;
        int yy,yyy;
        int zz,zzz;
        
        int xxshift;
        int yyshift;
        int zzshift;
        
        //int xxshiftEven;
        //int yyshiftEven;
        //int zzshiftEven;
        
        int posit;
        
        
        long t1=System.currentTimeMillis();
        
        
        
        //step 0 : no drift
        double frameNumberPerBin=((double)(maxF-minF))/((double)(bin));
        driftIndex[0]=((double)minF)+(frameNumberPerBin/2.);//mean of frame number
        driftX[0]=0;
        driftY[0]=0;
        driftZ[0]=0;
        
        
        
        //for each step
        
        int passNumber=bin-1;
        //int passNumber=1;
        
        
        int totNumb=0;
        for (int times=0;times<passNumber;times++){//times is the reference for cross correlation ... 0 at first... 1 after... etc
            for (int s=times+1;s<bin;s++){
                totNumb++;
            }
        }
        
        //IJ.log("result needs "+Math.round((totNumb*sizeFFTx*sizeFFTy*sizeFFTz)*4./1000000.)+" mB ");
        result = new float[totNumb][sizeFFTz][sizeFFTy][sizeFFTx];
        
        for (int i=0;i<result.length;i++){
            for (int ii=0;ii<result[0].length;ii++){
                for (int iii=0;iii<result[0][0].length;iii++){
                    for (int iiii=0;iiii<result[0][0][0].length;iiii++){
                        result[i][ii][iii][iiii]=0;
                    }
                }
            }
        }
        
        double [][] A=new double[totNumb][bin-1];
        for (int i=0;i<A.length;i++){
            for (int ii=0;ii<A[0].length;ii++){
                A[i][ii]=0;
            }
        }
        double [][] R=new double[3][totNumb];
        
        
        
        
        //float [][] stepSaved = new float[bin][sizeFFTPadded*sizeFFTPadded*sizeFFTPadded*2];
        float [][] stepSavedSmall = new float[bin][sizeFFTx*sizeFFTy*sizeFFTz];
        
        for (int i=0;i<stepSavedSmall.length;i++){
            for (int ii=0;ii<stepSavedSmall[0].length;ii++){
                stepSavedSmall[i][ii]=0;
            }
        }
        
        
        for (int u=0,uprocess=0;u<nbX;u++){
            for (int uu=0;uu<nbY;uu++){
                for (int uuu=0;uuu<nbZ;uuu++,uprocess++){
                    
                    
                    IJ.showProgress(((float)uprocess)/(float)(nbX*nbY*nbZ));
                    
                    //init images
                    //for (int s=0;s<bin;s++){

                        
                        int count=0;
                        loopi:for (int id=0;id<framePosition.length;id++){

                            posit=framePosition[id][1];

                            int idBin=bin*(stackloc.fl.get(posit).numFrame-minF)/(1+maxF-minF);
                            
                            

                            for (int j=0;j<stackloc.fl.get(posit).loc.size();j++){
                                if (stackloc.fl.get(posit).loc.get(j).exists){
                                    x=stackloc.fl.get(posit).loc.get(j).X-minX;
                                    y=stackloc.fl.get(posit).loc.get(j).Y-minY;
                                    z=stackloc.fl.get(posit).loc.get(j).Z-minZ;





                                    if ((x>=this.gridPosition[u][uu][uuu][0])&&(x<this.gridPosition[u][uu][uuu][1])&&(y>=this.gridPosition[u][uu][uuu][2])&&(y<this.gridPosition[u][uu][uuu][3])&&(z>=this.gridPosition[u][uu][uuu][4])&&(z<this.gridPosition[u][uu][uuu][5])){

                                        //make histogram padded:
                                        xx=Math.min(Math.max((int)(((x)%(sizeFFTx*pixelsizeNM))/pixelsizeNM),0),sizeFFTx-1);
                                        yy=Math.min(Math.max((int)(((y)%(sizeFFTy*pixelsizeNM))/pixelsizeNM),0),sizeFFTy-1);
                                        zz=Math.min(Math.max((int)(((z)%(sizeFFTz*pixelsizeNM))/pixelsizeNM),0),sizeFFTz-1);


                                        //xxshift=((((sizeFFT/2)+xx)+sizeFFTPadded/2)%sizeFFTPadded);//*2 because even
                                        //yyshift=((((sizeFFT/2)+yy)+sizeFFTPadded/2)%sizeFFTPadded);
                                        //zzshift=((((sizeFFT/2)+zz)+sizeFFTPadded/2)%sizeFFTPadded);

                                        
                                        if (stepSavedSmall[idBin][(zz*sizeFFTx*sizeFFTy+yy*sizeFFTx+xx)]<Float.MAX_VALUE){
                                            //stepSaved[idBin][(zzshift*sizeFFTPadded*sizeFFTPadded+yyshift*sizeFFTPadded+xxshift)*2]++;
                                            
                                            
                                            stepSavedSmall[idBin][(zz*sizeFFTx*sizeFFTy+yy*sizeFFTx+xx)]++;
                                            
                                            count++;
                                        }
                                        
                                    }


                                    


                                }
                            }
                            
                        }

                    //}
                    
                    
                    
                    
                    
                    //follow
                    for (int times=0,numbId=0;times<passNumber;times++){//times is the reference for cross correlation ... 0 at first... 1 after... etc
                        for (int i=0;i<input.length;i++){
                            input[i]=0;
                        }
                        for (int i=0;i<stepSavedSmall[times].length;i++){
                            
                            //xx=i%sizeFFT;
                            //zz=i/(sizeFFT*sizeFFT);
                            //yy=(i-xx-sizeFFT*sizeFFT*zz)/sizeFFT;

                            zz=i/(sizeFFTx*sizeFFTy);
                            yy=(i%(sizeFFTx*sizeFFTy))/sizeFFTx;
                            xx=(i%(sizeFFTx*sizeFFTy))%sizeFFTx;
                            
                            xxshift=((((sizeFFTx/2)+xx)+sizeFFTPaddedx/2)%sizeFFTPaddedx);//*2 because even
                            yyshift=((((sizeFFTy/2)+yy)+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                            zzshift=((((sizeFFTz/2)+zz)+sizeFFTPaddedz/2)%sizeFFTPaddedz);
                            input[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2]=stepSavedSmall[times][(zz*sizeFFTx*sizeFFTy+yy*sizeFFTx+xx)];
                            
                            
                            
                        }
                        
                        cudaResult =cudaMemcpyAsync(device_input, Pointer.to(input), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}

                        cudaResult =JCufft.cufftExecC2C(plan, device_input,device_input,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("drift ERROR cuda fft1 "+cudaResult);}
                        
                        MyVecDouble.mul_fl(MyCudaStream.getCUstream(0), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2, device_input, device_input, device_gaussian);

                        
                        for (int s=times+1;s<bin;s++,numbId++){
                            for (int i=0;i<step.length;i++){
                                step[i]=0;
                            }
                            for (int i=0;i<stepSavedSmall[s].length;i++){
                                //xx=i%sizeFFT;
                                //zz=i/(sizeFFT*sizeFFT);
                                //yy=(i-xx-sizeFFT*sizeFFT*zz)/sizeFFT;
                                
                                zz=i/(sizeFFTx*sizeFFTy);
                                yy=(i%(sizeFFTx*sizeFFTy))/sizeFFTx;
                                xx=(i%(sizeFFTx*sizeFFTy))%sizeFFTx;
                                
                                xxshift=((((sizeFFTx/2)+xx)+sizeFFTPaddedx/2)%sizeFFTPaddedx);//*2 because even
                                yyshift=((((sizeFFTy/2)+yy)+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                                zzshift=((((sizeFFTz/2)+zz)+sizeFFTPaddedz/2)%sizeFFTPaddedz);
                                step[(zzshift*sizeFFTPaddedy*sizeFFTPaddedx+yyshift*sizeFFTPaddedx+xxshift)*2]=stepSavedSmall[s][(zz*sizeFFTx*sizeFFTy+yy*sizeFFTx+xx)];
                            }
                            
                            cudaResult =cudaMemcpyAsync(device_step, Pointer.to(step), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(0));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
                            cudaResult =JCufft.cufftExecC2C(plan, device_step,device_step,JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("drift ERROR cuda fft1 "+cudaResult);}




                            //do cross correl here



                            MyVecDouble.complexeConjugateKernel(MyCudaStream.getCUstream(0), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz, device_step, device_input, device_step);



                            JCufft.cufftExecC2C(plan, device_step,device_step,JCufft.CUFFT_INVERSE);

                            //could be removed to win few time
                            MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(0), sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2, device_step, device_step, (float)Math.sqrt(sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz));

                            cudaResult =cudaMemcpyAsync(Pointer.to(step), device_step, sizeFFTPaddedx*sizeFFTPaddedy*sizeFFTPaddedz*2*Sizeof.FLOAT, cudaMemcpyDeviceToHost,MyCudaStream.getCudaStream_t(0));

                            for (xxx=0;xxx<sizeFFTx;xxx++){
                                for (yyy=0;yyy<sizeFFTy;yyy++){
                                    for ( zzz=0;zzz<sizeFFTz;zzz++){
                                        //xxshiftEven=((xxx+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;//*2 because even
                                        //yyshiftEven=((yyy+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;
                                        //zzshiftEven=((zzz+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;

                                        xxshift=((xxx+(sizeFFTx/2)+sizeFFTPaddedx/2)%sizeFFTPaddedx);//*2 because even
                                        yyshift=((yyy+(sizeFFTy/2)+sizeFFTPaddedy/2)%sizeFFTPaddedy);
                                        zzshift=((zzz+(sizeFFTz/2)+sizeFFTPaddedz/2)%sizeFFTPaddedz);

                                        result[numbId][zzz][yyy][xxx]+=step[(zzshift*sizeFFTPaddedx*sizeFFTPaddedy+yyshift*sizeFFTPaddedx+xxshift)*2];
                                    }
                                }
                            }
                            
                            
                        }
                    }
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                }
            }
        }
        
        
//        for (int i=0;i<result.length;i++){
//            ImageShow.imshow(result[i],"show "+i);
//        }
        
        //follow
        for (int times=0,numbId=0;times<passNumber;times++){//times is the reference for cross correlation ... 0 at first... 1 after... etc
            for (int s=times+1;s<bin;s++,numbId++){
                
                
                
                

                //ImageShow.imshow(result,"result "+s);



                LocalMaxima locmax = new LocalMaxima(result[numbId]);
                locmax.run(.95);
                if (times==0){
                    driftIndex[s]=(((double)s)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.);//mean of frame number
                    
                    driftX[s]=locmax.getDriftXinPix()*pixelsizeNM;
                    driftY[s]=locmax.getDriftYinPix()*pixelsizeNM;
                    driftZ[s]=locmax.getDriftZinPix()*pixelsizeNM;
                    
//                    IJ.log("drift X "+driftX[numbId]+"  nm");
//                    IJ.log("drift Y "+driftY[numbId]+"  nm");
//                    IJ.log("drift Z "+driftZ[numbId]+"  nm");
//                    
//                    IJ.log("");
                }
                
                R[0][numbId]=locmax.getDriftXinPix()*pixelsizeNM;
                R[1][numbId]=locmax.getDriftYinPix()*pixelsizeNM;
                R[2][numbId]=locmax.getDriftZinPix()*pixelsizeNM;
                
                
                
                //A[numbId][times]=-((((double)times)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.));
                
                //A[numbId][s]=((((double)s)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.));
                if (times!=0){
                    A[numbId][times-1]=-1;
                }
                
                A[numbId][s-1]=1;
                
                
            }
        }
        
        
        //double [][] AA={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0}};;
        //double [][] RR={{-47.77584571501592,-106.8553677959379,-0.7675404540464115},{-133.4594686872343,-132.7645498087783,-4.543981918508422},{-182.92575929587613,-101.1243637780396,-22.19715951653427},{-186.17080657785223,-82.66789094763425,-38.72344841612474},{-75.51451595803726,-16.538492083225975,-7.230904698883478},{-132.69461898940094,6.694369064555872,-19.976627901110078},{-130.36682727684124,34.06332981424143,-37.18433264789738},{-44.83839636955267,26.70459137640364,-10.569133284470666},{-56.44697715252249,48.07490244759762,-29.148866554700703},{-6.520611747652794,23.982255741189817,-15.682598027674999}};
        
        
        
        
        double [] threshold = new double [3];
        
        
        threshold[0]=precisionX;
        threshold[1]=precisionY;
        threshold[2]=precisionZ;
        
        //IJ.log("threshold X:"+threshold[0]+"  Y:"+threshold[1]+"  Z:"+threshold[2]);
                
//        IJ.log("RR :");
//        for (int u=0;u<totNumb;u++){
//            /*R[0][u]=RR[u][0];
//            R[1][u]=RR[u][1];
//            R[2][u]=RR[u][2];*/
//            IJ.log(""+R[0][u]+"  "+R[1][u]+"  "+R[2][u]+"  ");
//        }
//        IJ.log("");
//        IJ.log("AA");
//        for (int u=0;u<A.length;u++){
//            String s="";
//            for (int uu=0;uu<A[0].length;uu++){
//                //A[u][uu]=AA[u][uu];
//                s+=A[u][uu]+" ";
//            }
//            IJ.log(""+s);
//        }
        
//        IJ.log("drift before :");
//        for (int u=0;u<driftZ.length;u++){
//            IJ.log(""+driftX[u]+"  "+driftY[u]+"  "+driftZ[u]+"  ");
//        }
//        IJ.log("");
        
        
        if ((bin>2)&&(A.length>bin)){

            ArrayList<Integer> [] toDelete = new ArrayList[3];
            for (int i=0;i<3;i++){
                toDelete[i]= new ArrayList<Integer>();
            }






            for (int check=0;check<1;check++){//2 passes to remove outliers with precision < threshold
                
                
                
                double [][][] Anew=new double[3][][];
                double [][] Rnew=new double [3][];

                for (int dim=0;dim<3;dim++){
                    
                    
//                    IJ.log("dim "+dim);
                    int thelength=totNumb-toDelete[dim].size();
                    
                    Anew[dim]=new double[thelength][bin-1];
                    
                    if (check==1){
//                        IJ.log("reject number="+toDelete[dim].size());
                    }

                    Rnew[dim]=new double [thelength];


                    for (int i=0,it=0;i<A.length;i++){
                        boolean reject=false;
                        for (int u=0;u<toDelete[dim].size();u++){
                            if (i==(int)toDelete[dim].get(u)){
                                reject=true;
                            }
                        }
                        if (!reject){
                            for (int ii=0;ii<A[0].length;ii++){
                                Anew[dim][it][ii]=A[i][ii];
                            }
                            Rnew[dim][it]=R[dim][i];
                            it++;
                        }
                    }
                    
                    //IJ.log("todelok ");
                    double [] residual = new double [thelength];


                    if (thelength>=bin-1){
                        //try{
                            
                            Matrixe m = new Matrixe(Anew[dim]);
                            Matrixe r = new Matrixe(Rnew[dim]);
                            Matrixe mt = Matrixe.transpose(m);
                            Matrixe mm = mt.times(m);
                            
                            //too long... better with Gaussian elimination (pivot gauss)
                            //Matrixe minv;
                            //try{
                            //     minv = Matrixe.inverse(mm);
                            //}catch(Exception eee){IJ.log("oops matrix non inversible -> maybe result will not be super precise");continue;}
                            
                            //Matrixe M = minv.times(mt);
                            

                            //Matrixe d = M.times(r);

                            //
                            
                            
                            Matrixe mr = mt.times(r);
                            
                            double [][] myMatA=mm.getMatrixe();
                            double [][] myMatR=mr.getMatrixe();
                            
                            
                            
                            double [][] dresGE=GaussianElimination.lsolve(myMatA, myMatR);
                            
                            Matrixe d= new Matrixe(dresGE);
                            Matrixe res = m.times(d);
                            double [][] res_=res.getMatrixe();
                            
                            for (int i=0;i<residual.length;i++){
                                
                                
                                residual[i]= Math.abs(res_[i][0]-Rnew[dim][i]);

                                if (residual[i]>threshold[dim]){
                                    toDelete[dim].add(i);
                                }
//                                IJ.log("residualOld ("+dim+") "+residual[i]);
                            }
                            
                            
                            for (int u=0;u<dresGE.length;u++){
                                
                                if (dim==0){
                                    driftX[u+1]=dresGE[u][0];
                                }
                                else if (dim==1){
                                    driftY[u+1]=dresGE[u][0];
                                }
                                else if (dim==2){
                                    driftZ[u+1]=dresGE[u][0];
                                }
                            }
                            
                            
                            /*
                            //IJ.log("solve... ");
                            double [] D=MatrixColt.solve(Anew[dim], Rnew[dim]);
                            //IJ.log("solve ok ");
                            
                            
                            
                            double [] Ad=new double[Rnew[dim].length];
                            for (int i=0;i<residual.length;i++){
                                Ad[i]=0;
                                for (int ii=0;ii<Anew[dim][i].length;ii++){
                                    Ad[i]+=Anew[dim][i][ii]*D[ii];
                                }
                                residual[i]= Math.abs(Ad[i]-Rnew[dim][i]);

                                if (residual[i]>threshold[dim]){
                                    toDelete[dim].add(i);
                                }
                                //IJ.log("residualNew ("+dim+") "+residual[i]);
                            }
                            
                            
                            for (int u=0;u<D.length;u++){
                                
                                if (dim==0){
                                    driftX[u+1]=D[u];
                                }
                                else if (dim==1){
                                    driftY[u+1]=D[u];
                                }
                                else if (dim==2){
                                    driftZ[u+1]=D[u];
                                }
                            }*/



                        //}
                        //catch(Exception eee){IJ.log("oops matrix non inversible -> maybe result will not be super precise");}
                    }



                }

//                IJ.log("drift after :");
//                for (int u=0;u<driftZ.length;u++){
//                    IJ.log(""+driftX[u]+"  "+driftY[u]+"  "+driftZ[u]+"  ");
//                }
//                IJ.log("");


            }
            
            
        }

        
        
        
        
        
        //plot(driftIndex,driftX,"drift of X","frame number","drift (nm)");
        
        //plot(driftIndex,driftY,"drift of Y","frame number","drift (nm)");
        
        //plot(driftIndex,driftZ,"drift of Z","frame number","drift (nm)");
        
        
        
        
        
        double [][] di=new double [1][bin];
        double [][] dx=new double [1][bin];
        double [][] dy=new double [1][bin];
        double [][] dz=new double [1][bin];
        
        
        
        for (int i=0;i<bin;i++){
            driftX[i]+=this.initDriftX;
            driftY[i]+=this.initDriftY;
            driftZ[i]+=this.initDriftZ;
            dx[0][i]=driftX[i];
            dy[0][i]=driftY[i];
            dz[0][i]=driftZ[i];
            di[0][i]=driftIndex[i];
            
        }
        
        SplineFit pfX = new SplineFit(dx,di);
        pfX.run();
        SplineFit pfY = new SplineFit(dy,di);
        pfY.run();
        SplineFit pfZ = new SplineFit(dz,di);
        pfZ.run();
        
        double [] li=new double[maxF-minF];
        double [] lx=new double[maxF-minF];
        double [] ly=new double[maxF-minF];
        double [] lz=new double[maxF-minF];
        
        double [] v=new double[1];//for poly fit result
        double [] vv=new double[1];//for poly fit result
        for (int i=0;i<maxF-minF;i++){
            li[i]=i+minF;
            v[0]=li[i];
            vv=pfX.transform(v);
            lx[i]=vv[0];
            vv=pfY.transform(v);
            ly[i]=vv[0];
            vv=pfZ.transform(v);
            lz[i]=vv[0];
        }
        
        
        
        for (int i=0;i<stackloc.fl.size();i++){
            for (int ii=0;ii<stackloc.fl.get(i).loc.size();ii++){
                v[0]=stackloc.fl.get(i).loc.get(ii).frame;
                vv=pfX.transform(v);
                stackloc.fl.get(i).loc.get(ii).setDrift_X(vv[0]);
                vv=pfY.transform(v);
                stackloc.fl.get(i).loc.get(ii).setDrift_Y(vv[0]);
                vv=pfZ.transform(v);
                stackloc.fl.get(i).loc.get(ii).setDrift_Z(vv[0]);
            }
        }
        
        
        plot(driftIndex,driftX,li,lx,"X_drift","frame number","X drift (nm)");
        
        plot(driftIndex,driftY,li,ly,"Y_drift","frame number","Y drift (nm)");
        
        plot(driftIndex,driftZ,li,lz,"Z_drift","frame number","Z drift (nm)");
        
        long t2=System.currentTimeMillis();

        //IJ.log("t2-t1 "+(t2-t1)); 
        
        IJ.showProgress(1);
        int minut=(int)(3600000.*(((double)(System.currentTimeMillis()-time0)/3600000.)-(double)((System.currentTimeMillis()-time0)/3600000)))/60000;
        int hour=(int)((System.currentTimeMillis()-time0)/3600000);
        if (hour+minut>0){
            IJ.log("elapsed time "+hour+" h "+minut+" m");
        }
        else{
            IJ.log("elapsed time "+((System.currentTimeMillis()-time0)/1000)+" second");
        }
        
        IJ.showProgress(0);
        IJ.log("Drift correction finished");
        //IJ.log("");
        //IJ.log("In case of extremely large drift, a residual drift can happend using cross-correlation due to .");
        //IJ.log("Please, apply drift correction a second time if you think the correction is not perfect");
    }
    
    
    
    
    
    
    
    
    
    
    public void plot(double [] x, double [] y,double [] xline,double [] yline,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,xline,yline);
        
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
        for (int i=0;i<yline.length;i++){
            if (xline[i]<xMin){
                xMin=xline[i];
            }
            if (xline[i]>xMax){
                xMax=xline[i];
            }
            if (yline[i]<yMin){
                yMin=yline[i];
            }
            if (yline[i]>yMax){
                yMax=yline[i];
            }
        }
        p.setLimits(xMin-200, xMax+200, yMin-200, yMax+200);
        p.setLineWidth(2);
        p.setFont(0, 18);
        for (int ii=0;ii<x.length;ii++){
            
            p.setColor(Color.blue);
            p.add("CIRCLE",  x,y);
            p.setColor(Color.red);
            //p.show();


            
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
    
    
    
    
    void free(){
        
        
        JCuda.cudaFree(device_input) ;
        JCuda.cudaFree(device_step) ;
        JCuda.cudaFree(device_gaussian) ;
        JCufft.cufftDestroy(plan);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    class LocalMaxima{
        float [][][] image;
        double x,y,z;//position to estimate -> drift
        int sizePatchx=11;//odd value;
        int sizePatchy=11;//odd value;
        int sizePatchz=11;//odd value;
        LocalMaxima(float [][][] image){
            this.image=image;
            width=image.length;
            height=image[0].length;
            depth=image[0][0].length;
            
            sizePatchx=Math.min(sizePatchx, width);
            sizePatchy=Math.min(sizePatchy, height);
            sizePatchz=Math.min(sizePatchz, depth);
            
            if ((width%2==1)||(height%2==1)||(depth%2==1)){
                //IJ.log("WARNING, drift image computation should have even size");
            }
            
            if (sizePatchx%2==0){
                sizePatchx=sizePatchx-1;
            }
            if (sizePatchy%2==0){
                sizePatchy=sizePatchy-1;
            }
            if (sizePatchz%2==0){
                sizePatchz=sizePatchz-1;
            }
            sizePatchx=Math.max(sizePatchx, 1);
            sizePatchy=Math.max(sizePatchy, 1);
            sizePatchz=Math.max(sizePatchz, 1);
        }
        
        
        public void run(double threshold){//threshold=.75 -> 75% of data used to compute center
            
            
            double max=Double.NEGATIVE_INFINITY;
            
            int xx=0,yy=0,zz=0;
            
            for (int i=0,t=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        if (max<image[i][ii][iii]){
                            max=image[i][ii][iii];
                            xx=i;
                            yy=ii;
                            zz=iii;
                        }
                    }
                }
            }
            
            
            float [] vect=new float[sizePatchx*sizePatchy*sizePatchz];
            int number=0;
            for (int i=0;i<sizePatchx;i++){
                for (int ii=0;ii<sizePatchy;ii++){
                    for (int iii=0;iii<sizePatchz;iii++){
                        if ((i+xx-sizePatchx/2>=0)&&(ii+yy-sizePatchy/2>=0)&&(iii+zz-sizePatchz/2>=0)&&(i+xx-sizePatchx/2<width)&&(ii+yy-sizePatchy/2<height)&&(iii+zz-sizePatchz/2<depth)){
                            vect[number++]=image[i+xx-sizePatchx/2][ii+yy-sizePatchy/2][iii+zz-sizePatchz/2];
                        }
                    }
                }
            }
            //IJ.log("number "+number+"  "+xx+"  "+yy+"  "+zz);
            Arrays.sort(vect);
            
            //IJ.log("threshold "+vect[(int)((number-1)*threshold)]+"   "+max);
            
            x=0;
            y=0;
            z=0;
            double sum=0;
            for (int i=0;i<sizePatchx;i++){
                for (int ii=0;ii<sizePatchy;ii++){
                    for (int iii=0;iii<sizePatchz;iii++){
                        if ((i+xx-sizePatchx/2>=0)&&(ii+yy-sizePatchy/2>=0)&&(iii+zz-sizePatchz/2>=0)&&(i+xx-sizePatchx/2<width)&&(ii+yy-sizePatchy/2<height)&&(iii+zz-sizePatchz/2<depth)){
                            if (image[i+xx-sizePatchx/2][ii+yy-sizePatchy/2][iii+zz-sizePatchz/2]>vect[(int)((number-1)*threshold)]){
                                sum+=image[i+xx-sizePatchx/2][ii+yy-sizePatchy/2][iii+zz-sizePatchz/2];
                                x+=(i+xx-sizePatchx/2)*image[i+xx-sizePatchx/2][ii+yy-sizePatchy/2][iii+zz-sizePatchz/2];
                                y+=(ii+yy-sizePatchy/2)*image[i+xx-sizePatchx/2][ii+yy-sizePatchy/2][iii+zz-sizePatchz/2];
                                z+=(iii+zz-sizePatchz/2)*image[i+xx-sizePatchx/2][ii+yy-sizePatchy/2][iii+zz-sizePatchz/2];
                            }
                        }
                    }
                }
            }
            //IJ.log("posit "+x+"  "+y+"  "+z);
            x/=sum;
            y/=sum;
            z/=sum;
            //IJ.log("posit "+x+"  "+y+"  "+z);
        }
        
        public double getDriftZinPix(){
            return -(x-width/2);
        }
        
        public double getDriftYinPix(){
            return -(y-height/2);
        }
        
        public double getDriftXinPix(){
            return -(z-depth/2);
        }
        
    }
    
    class Fit{
        
        float [][][] image;
        double [][][] model;
        double [][][] partial;
        double [][] matVarCovar;
        double [][] sig;
        double sigA, sigB, sigC, sigE, sigF, sigI;//Mat var/covar in 3D -> 9 elements
        double x,y,z;//position to estimate -> drift
        double A; //photon number
        double likelihood;
        int width,height,depth;
        
        Fit(float [][][] image){
            this.image=image;
            width=image.length;
            height=image[0].length;
            depth=image[0][0].length;
            
            //search max      
            
            
            model=new double [width][height][depth];
            partial=new double [width][height][depth];
            init();
        }
        
        
        void init(){
            
            sigA=.1;
            sigB=0;
            sigC=0;
            sigE=.1;
            sigF=0;
            sigI=.1;
            A=0;
            double max=Double.NEGATIVE_INFINITY;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        if (image[i][ii][iii]>max){
                            max=image[i][ii][iii];
                            x=i;
                            y=ii;
                            z=iii;
                        }
                        
                    }
                }
            }
            A=max;
            matVarCovar=new double[3][3];
            for (int i=0;i<3;i++){//if non inversible -> inverse diagonale
                for (int ii=0;ii<3;ii++){
                    if (i!=ii){
                        matVarCovar[i][ii]=0;
                    }
                }
            }
            matVarCovar[0][0]=width;
            matVarCovar[1][1]=height;
            matVarCovar[2][2]=depth;
            sig=new double[3][3];//matVarCovar inverse
            computeSigmaInv();
            
        }
        
        void computeSigmaInv(){
            Matrixe m = new Matrixe(matVarCovar);
            try{
                m=Matrixe.inverse(m);
                sig=m.getMatrixe();
            }catch(Exception ee){
                for (int i=0;i<3;i++){//if non inversible -> inverse diagonale
                    for (int ii=0;ii<3;ii++){
                        if (i!=ii){
                            sig[i][ii]=0;
                        }
                    }
                }
                sig[0][0]=1/matVarCovar[0][0];
                sig[1][1]=1/matVarCovar[1][1];
                sig[2][2]=1/matVarCovar[2][2];
            }
            sigA=sig[0][0];
            sigE=sig[1][1];
            sigI=sig[2][2];
            sigB=sig[1][0];
            sigC=sig[2][0];
            sigF=sig[1][2];
        }
        
        public double [][][] getModel(){
            computeModel();
            return model;
        }
        public void run(int iter){
            x-=2;
            for (int t=0;t<iter;t++){
                computeLikelihood();
                IJ.log("lik:"+likelihood+"      A:"+A+"       x:"+x+"   y:"+y+"   z:"+z+"      sigA:"+sigA+"      sigE:"+sigE+"      sigI:"+sigI+"   ");
                
                IJ.log("sigA:"+sigA+"      sigE:"+sigE+"      sigI:"+sigI+"   "+"      sigB:"+sigB+"   "+"      sigC:"+sigC+"   "+"      sigF:"+sigF+"   ");
                
                A-=computePartialModel_A();
                x-=computePartialModel_X();
                y-=computePartialModel_Y();
                z-=computePartialModel_Z();
                
                
                sigA-=.01*computePartialModel_sigA();
                
                
                //sigA-=.01*computePartialModel_sigAOffset(.0000001);
                //sigE-=.01*computePartialModel_sigE();
                //sigI-=.01*computePartialModel_sigI();
                
            }
            
        }
        
        
        void computeModel(){
            double xx,yy,zz,d;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        model[i][ii][iii]=A*Math.exp(-.5*d);
                    }
                }
            }
        }
        
        
        double computeLikelihood(){
            double xx,yy,zz,d;
            likelihood=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        model[i][ii][iii]=A*Math.exp(-.5*d);
                        likelihood+=(model[i][ii][iii]-image[i][ii][iii])*(model[i][ii][iii]-image[i][ii][iii]);
                    }
                }
            }
            return likelihood;
        }
        
        double computePartialModel_X(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);//distance
                        p=-.5*(   (xx*sigA+yy*sigB+zz*sigC)+sigA*xx +yy*sigB+ zz*sigC   );//partial of distance
                        M=A*Math.exp(-.5*d);//model
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2=2*(-.5*(   sigA*2   ) * M  +  p*p*M);//partial of 2*p*M
                        som2+=( (M-image[i][ii][iii])*p2 + p*M*2*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        
        
        
        double computePartialModel_Y(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   xx*sigB+ (xx*sigB+yy*sigE+zz*sigF)+sigE*yy +zz*sigF   );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2=2*(-.5*(   sigE*2   ) * M  +  p*p*M);//partial of 2*p*M
                        som2+=( (M-image[i][ii][iii])*p2 + p*M*2*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        double computePartialModel_Z(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   xx*sigC +yy*sigF + (xx*sigC+yy*sigF+zz*sigI)+sigI*zz   );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2=2*(-.5*(   sigI*2   ) * M  +  p*p*M);//partial of 2*p*M
                        som2+=( (M-image[i][ii][iii])*p2 + p*M*2*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        
        
        
        double computePartialModel_sigA(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   xx*xx  );
                        M=A*Math.exp(-.5*d);
                        
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2= 2*p*p*M;
                        som2+=( (M-image[i][ii][iii])*p2 + 2*p*M*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        
        double computePartialModel_sigAOffset(double hdec){
            
            double xx,yy,zz,d0,d1,d2,p,C0,C1,C2,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d0=xx*(xx*(sigA-hdec)+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        C0=(A*Math.exp(-.5*d0)-image[i][ii][iii])*(A*Math.exp(-.5*d0)-image[i][ii][iii]);
                        d1=xx*(xx*(sigA)+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        C1=(A*Math.exp(-.5*d1)-image[i][ii][iii])*(A*Math.exp(-.5*d1)-image[i][ii][iii]);
                        d2=xx*(xx*(sigA+hdec)+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        C2=(A*Math.exp(-.5*d2)-image[i][ii][iii])*(A*Math.exp(-.5*d2)-image[i][ii][iii]);
                        
                        som+=(C2-C0)/(2*hdec);
                        
                        som2+=(C2+C0-2*C1)/(hdec*hdec);
                    }
                }
            }
            return som/som2;
        }
        
        double computePartialModel_sigB(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   xx*yy*2 );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2= 2*p*p*M;
                        som2+=( (M-image[i][ii][iii])*p2 + 2*p*M*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        double computePartialModel_sigC(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   xx*zz*2  );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2= 2*p*p*M;
                        som2+=( (M-image[i][ii][iii])*p2 + 2*p*M*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        
        double computePartialModel_sigE(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   yy*yy  );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2= 2*p*p*M;
                        som2+=( (M-image[i][ii][iii])*p2 + 2*p*M*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        
        double computePartialModel_sigF(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   yy*zz*2  );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2= 2*p*p*M;
                        som2+=( (M-image[i][ii][iii])*p2 + 2*p*M*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        double computePartialModel_sigI(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        p=-.5*(   zz*zz  );
                        M=A*Math.exp(-.5*d);
                        
                        pp=2*p*M*(M-image[i][ii][iii]);//partial1
                        som+=pp;
                        p2= 2*p*p*M;
                        som2+=( (M-image[i][ii][iii])*p2 + 2*p*M*p*M );//partial 2
                    }
                }
            }
            return som/som2;
        }
        
        
        double computePartialModel_A(){
            
            double xx,yy,zz,d,p,M,p2,pp;
            double som=0;
            double som2=0;
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        xx=x-i;
                        yy=y-ii;
                        zz=z-iii;
                        d=xx*(xx*sigA+yy*sigB+zz*sigC)+yy*(xx*sigB+yy*sigE+zz*sigF)+zz*(xx*sigC+yy*sigF+zz*sigI);
                        M=A*Math.exp(-.5*d);
                        som+=Math.exp(-.5*d)*(M-image[i][ii][iii]);
                        som2+=Math.exp(-.5*d)*Math.exp(-.5*d);
                    }
                }
            }
            return som/som2;
        }
        
        
    }
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    



}
