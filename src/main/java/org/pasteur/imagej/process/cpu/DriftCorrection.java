/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;
import org.pasteur.imagej.utils.Matrixe;
import org.pasteur.imagej.utils.FourierTransform;
import org.pasteur.imagej.utils.SplineFit;
import org.pasteur.imagej.utils.GaussianElimination;
import org.pasteur.imagej.data.*;
/**
 *
 * @author benoit
 */
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import java.util.List;


//import MatrixColt;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import java.util.ArrayList;


/**
 *
 * @author benoit
 */
public class DriftCorrection {
    
    double nbMem=0;
    
    private  Object monitor = new Object();
    
    int thread=Runtime.getRuntime().availableProcessors();
    FourierTransform ft = new FourierTransform();
    
    double precisionX;
    double precisionY;
    double precisionZ;
    
    
    
    double [] driftIndex;
    double [] driftX;
    double [] driftY;
    double [] driftZ;
    
    
    /*float [][][]  device_input_real;
    float [][][]  device_step_real;
    float [][][]  device_input_imag;
    float [][][]  device_step_imag;
    float [][][]  device_gaussian;*/
    
    
    
    float [][][]  host_input;
    float [][][]  host_step;
    
    
    float [][][][] image;
    
    float [][][][] output_real;
    float [][][][] step_real;
    float [][][][] output_imag;
    float [][][][] step_imag;
    
    
    float [][][] gaussian;
    
    float [][][][] result;
    
    
    int [][] framePosition = new int[2][];
    
    StackLocalization stackloc;
    
    int bin;
    
    
    
    
    int sizeFFTz=512;
    int sizeFFTy=512;
    int sizeFFTx=512;
    
    
    
    
    
    int sizeFFTPaddedx;
    int sizeFFTPaddedy;
    int sizeFFTPaddedz;
    
    double [][][][] gridPosition;//X;Y;Z;debut,fin
    
    int nbX,nbY,nbZ;
    
    
    int cudaResult;
    //In this class, the rendered 3D image is splitted into many 3D images of size sizeFFT to deal with low memory on GPU
    
    
    double minX,maxX,minY,maxY,minZ,maxZ;

    int minF,maxF;
               
    int width,height,depth;
            
    double pixelsizeNM;

    
    //WARNING : [x][y][z] instead of [z][x][y]
    public DriftCorrection(StackLocalization stackloc,double pixelsizeNM,int bin, double subImageSizeUm){
    
        
        
        
        int sizeFFT=(int)Math.pow(2,(int)(Math.ceil(Math.log((subImageSizeUm/(pixelsizeNM/1000.)))/Math.log(2))));
        
        int nbloop=0;
        while (sizeFFT>128){
            pixelsizeNM*=1.01;
            sizeFFT=(int)Math.pow(2,(int)(Math.ceil(Math.log((subImageSizeUm/(Math.round(pixelsizeNM)/1000.)))/Math.log(2))));
            nbloop++;
        }
        pixelsizeNM=Math.round(pixelsizeNM);
        
        
        int sizeFFTPadded;
        
        
        if(nbloop>1){
            IJ.log("Pixel size increased to "+pixelsizeNM+" nm to deal with CPU memory");
            IJ.log("You can also try do reduce the maximum_drift parameter if you need better precision");
        }
        
        driftIndex=new double[bin];
        driftX=new double[bin];
        driftY=new double[bin];
        driftZ=new double[bin];
        
        
        
        /*long memfree = Runtime.getRuntime().freeMemory();
        //IJ.log("mem available "+memfree[0]+" / "+memtot[0]);
       
        
        long neededMemory=((long)(3*2))*((long)sizeFFT*(long)sizeFFT*(long)sizeFFT*(long)8*(long)2 * (long)Sizeof.FLOAT)+(long)100000000;//+200 mo free
        //IJ.log("need memory   "+neededMemory+"  "+3*(sizeFFT*sizeFFT*sizeFFT*8*2 * Sizeof.FLOAT));
        
        while (neededMemory>memfree){
            sizeFFT*=.9;
            if (sizeFFT%2==1){
                sizeFFT++;
            }
            neededMemory=((long)(3*2))*((long)sizeFFT*(long)sizeFFT*(long)sizeFFT*(long)8*(long)2 * (long)Sizeof.FLOAT)+(long)100000000;//+200 mo free
        }
        
        IJ.log("size FFT according to available memory: "+sizeFFT+"  need mem:"+neededMemory);*/
        
        /*sizeFFT=256;
        sizeFFTx=256;
        sizeFFTy=256;
        sizeFFTz=128;*/
        
        
        sizeFFTPadded=sizeFFT*2;
        sizeFFTPaddedx=sizeFFT*2;
        sizeFFTPaddedy=sizeFFT*2;
        sizeFFTPaddedz=sizeFFT*2;
        
        this.pixelsizeNM=pixelsizeNM;
        this.bin=bin;
        
        this.stackloc=stackloc;
        
        minX=Double.MAX_VALUE;
        maxX=Double.MIN_VALUE;
        
        minY=Double.MAX_VALUE;
        maxY=Double.MIN_VALUE;
        
        minZ=Double.MAX_VALUE;
        maxZ=Double.MIN_VALUE;
        
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
        
        width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        depth=(int)Math.ceil((maxZ-minZ)/pixelsizeNM);
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        depth=Math.max(depth, 1);
        
        sizeFFTx=Math.min(sizeFFT, width);
        sizeFFTy=Math.min(sizeFFT, height);
        sizeFFTz=Math.min(sizeFFT, depth);
        
        sizeFFTx=Math.max(sizeFFTx, (int)Math.pow(2, Math.ceil(Math.log(sizeFFTx)/Math.log(2))));
        sizeFFTy=Math.max(sizeFFTy, (int)Math.pow(2, Math.ceil(Math.log(sizeFFTy)/Math.log(2))));
        sizeFFTz=Math.max(sizeFFTz, (int)Math.pow(2, Math.ceil(Math.log(sizeFFTz)/Math.log(2))));
        
        
        nbX=(int)Math.ceil((double)(width-this.sizeFFTx/2)/(double)(this.sizeFFTx/2));//   /2 car chevauchement
        nbY=(int)Math.ceil((double)(height-this.sizeFFTy/2)/(double)(this.sizeFFTy/2));//   /2 car chevauchement
        nbZ=(int)Math.ceil((double)(depth-this.sizeFFTz/2)/(double)(this.sizeFFTz/2));//   /2 car chevauchement
        
        //IJ.log("grid size "+nbX+"   "+nbY+"   "+nbZ);
        
        
        
        //nbX=nbX*2-1;//*2-1 car chevauchement
        //nbY=nbY*2-1;//*2-1 car chevauchement
        //nbZ=nbZ*2-1;//*2-1 car chevauchement
        
        //IJ.log("grid size "+nbX+"   "+nbY+"   "+nbZ);
        
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
        
        sizeFFTPaddedx=sizeFFTx*2;
        sizeFFTPaddedy=sizeFFTy*2;
        sizeFFTPaddedz=sizeFFTz*2;
        
        Arrays.sort(framePosition, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                return ((Integer) o1[0]).compareTo(o2[0]);
            }
        });
        
        int xx;
        int yy;
        int zz;
        
        
        nbMem+=(5.*(double)bin*(double)sizeFFTPaddedz*(double)sizeFFTPaddedy*(double)sizeFFTPaddedx*4.)/1000000.+(double)sizeFFTPaddedz*(double)sizeFFTPaddedy*(double)sizeFFTPaddedx*4./1000000.;
        
        //IJ.log("needed CPU memory > "+(double)Math.round(nbMem/100.)/10.+" GB");
        //IJ.log("if this value is too big, try to increase the rendering pixel size or reduce the maximum drift bound");
        
        image=new float[bin][sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        
        step_real=new float[bin][sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        step_imag=new float[bin][sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        
        output_real=new float[bin][sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        output_imag=new float[bin][sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        
        gaussian=new float[sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        
        
        double sigmaNM=3*Math.max(Math.max(precisionX, precisionY), precisionZ);
        double sigma=sigmaNM/pixelsizeNM;
        double sigmaFourierPowx=((1./(sigma))*(double)sizeFFTPaddedx/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTPaddedx/(Math.PI*2.));
        
        double sigmaFourierPowy=((1./(sigma))*(double)sizeFFTPaddedy/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTPaddedy/(Math.PI*2.));
        
        double sigmaFourierPowz=((1./(sigma))*(double)sizeFFTPaddedz/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTPaddedz/(Math.PI*2.));
        
        
        //float [][][] inputTMP=new float[sizeFFTPadded][sizeFFTPadded][sizeFFTPadded];
        
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
                    
                    gaussian[zzshift][yyshift][xxshift]=(float)Math.exp(-.5* (dx/sigmaFourierPowx + dy/sigmaFourierPowy +dz/sigmaFourierPowz));
                    
                    
                }
            }
        }
        
        runFast();
        
        
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
        //int [][] processPosition = new int [passNumber][];
        for (int times=0;times<passNumber;times++){//times is the reference for cross correlation ... 0 at first... 1 after... etc
            //processPosition[times]=new int[bin-times-1];
            for (int s=times+1;s<bin;s++){
                //processPosition[times][s]=totNumb;
                totNumb++;
            }
        }
        
        int numberProcess=Math.min(thread, totNumb);
        
        int [][] processPosition = new int [totNumb][2];
        for (int times=0,ip=0;times<passNumber;times++){
            for (int s=times+1;s<bin;s++){
                processPosition[ip][0]=times;
                processPosition[ip][1]=s;
                ip++;
            }
        }
        
        
        
        
        double truc1=(double)totNumb*(double)sizeFFTz*(double)sizeFFTy*(double)sizeFFTx*4./1000000.;
        double truc2=(double)numberProcess*(double)sizeFFTPaddedz*(double)sizeFFTPaddedy*(double)sizeFFTPaddedx*2.*4./1000000.;
        
        nbMem+=truc1;
        
        nbMem+=truc2;
        
        //IJ.log("needed CPU memory > "+(float)Math.round(nbMem/100.)/10.+" GB");
        
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
        
        
        
        
        
        
        for (int u=0,uprocess=0;u<nbX;u++){
            for (int uu=0;uu<nbY;uu++){
                for (int uuu=0;uuu<nbZ;uuu++,uprocess++){
                    
                    
                    
                    
                    IJ.showProgress(((float)uprocess)/(float)(nbX*nbY*nbZ));
                    
                    
                    //init
                    for (int i=0;i<step_real.length;i++){
                        for (int ii=0;ii<step_real[0].length;ii++){
                            for (int iii=0;iii<step_real[0][0].length;iii++){
                                for (int iiii=0;iiii<step_real[0][0][0].length;iiii++){
                                    step_real[i][ii][iii][iiii]=0;
                                    step_imag[i][ii][iii][iiii]=0;
                                    image[i][ii][iii][iiii]=0;
                                }
                            }
                        }
                    }


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
                                    xx=sizeFFTx/2+Math.min(Math.max((int)(((x)%(sizeFFTx*pixelsizeNM))/pixelsizeNM),0),sizeFFTx-1);
                                    yy=sizeFFTy/2+Math.min(Math.max((int)(((y)%(sizeFFTy*pixelsizeNM))/pixelsizeNM),0),sizeFFTy-1);
                                    zz=sizeFFTz/2+Math.min(Math.max((int)(((z)%(sizeFFTz*pixelsizeNM))/pixelsizeNM),0),sizeFFTz-1);


                                    //xxshift=((((sizeFFT/2)+xx)+sizeFFTPadded/2)%sizeFFTPadded);//*2 because even
                                    //yyshift=((((sizeFFT/2)+yy)+sizeFFTPadded/2)%sizeFFTPadded);
                                    //zzshift=((((sizeFFT/2)+zz)+sizeFFTPadded/2)%sizeFFTPadded);

                                    if (image[idBin][zz][yy][xx]<Float.MAX_VALUE){
                                            //stepSaved[idBin][(zzshift*sizeFFTPadded*sizeFFTPadded+yyshift*sizeFFTPadded+xxshift)*2]++;


                                        image[idBin][zz][yy][xx]++;

                                        count++;

                                    }





                                }
                            }

                        }
                    }
                    
                    Init [] init = new Init[bin];
                    
                    for (int i=0;i<bin;i++){
                        init[i]=new Init(i);
                        init[i].start();
                    }
                    for (int i=0;i<bin;i++){
                        try{
                            init[i].join();
                        }catch(Exception e){IJ.log("join init impossible "+e);}
                    }
                    
                    
                    //here, step contains all fft images for one position
                    
                    
                    int nbBlock=processPosition.length/numberProcess;
                    int sizeLastBlock=processPosition.length%numberProcess;
                    if (sizeLastBlock!=0){
                        nbBlock++;
                    }
                    
                    
                    Process [] process = new Process [numberProcess];
                    
                    for (int block=0,id=0;block<nbBlock;block++){
                        int processNumber=numberProcess;
                        if ((sizeLastBlock!=0)&&(block==nbBlock-1)){
                            processNumber=sizeLastBlock;
                        }
                        
                        
                        
                        for (int ip=0;ip<processNumber;ip++){
                            
                            process[ip]=new Process(ip,id,processPosition[id][0],processPosition[id][1]);
                            process[ip].start();
                            
                            
                            
                            id++;
                        }
                        
                        for (int ip=0;ip<processNumber;ip++){
                            try{
                                process[ip].join();
                            }catch(Exception e){IJ.log("join process impossible "+e);}
                        }
                    }
                    
                    process=null;
                    
                    for (int i=0;i<totNumb;i++){
                        //ImageShow.imshow(result[i],"res "+processPosition[i][0]+" "+processPosition[i][1]);
                    }
                    
                            
                          
                }
            }
        }
        
                    

        
        
        //follow
        for (int numbId=0;numbId<totNumb;numbId++){
            int times=processPosition[numbId][0];
            int s=processPosition[numbId][1];
            
                
                
            


            LocalMaxima locmax = new LocalMaxima(result[numbId]);
            locmax.run(.95);
            if (times==0){
                driftIndex[s]=(((double)s)/((double)bin))*((double)(1.+maxF-minF))+((double)minF)+(frameNumberPerBin/2.);//mean of frame number

                driftX[s]=locmax.getDriftXinPix()*pixelsizeNM;
                driftY[s]=locmax.getDriftYinPix()*pixelsizeNM;
                driftZ[s]=locmax.getDriftZinPix()*pixelsizeNM;

                //IJ.log("drift X "+driftX[numbId]+"  nm");
                //IJ.log("drift Y "+driftY[numbId]+"  nm");
                //IJ.log("drift Z "+driftZ[numbId]+"  nm");

                //IJ.log("");
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
        
        
        //double [][] AA={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0}};;
        //double [][] RR={{-47.77584571501592,-106.8553677959379,-0.7675404540464115},{-133.4594686872343,-132.7645498087783,-4.543981918508422},{-182.92575929587613,-101.1243637780396,-22.19715951653427},{-186.17080657785223,-82.66789094763425,-38.72344841612474},{-75.51451595803726,-16.538492083225975,-7.230904698883478},{-132.69461898940094,6.694369064555872,-19.976627901110078},{-130.36682727684124,34.06332981424143,-37.18433264789738},{-44.83839636955267,26.70459137640364,-10.569133284470666},{-56.44697715252249,48.07490244759762,-29.148866554700703},{-6.520611747652794,23.982255741189817,-15.682598027674999}};
        
        
        
        
        double [] threshold = new double [3];
        
        
        threshold[0]=precisionX;
        threshold[1]=precisionY;
        threshold[2]=precisionZ;
        
        //IJ.log("threshold X:"+threshold[0]+"  Y:"+threshold[1]+"  Z:"+threshold[2]);
                
        
        if ((bin>2)&&(A.length>bin)){

            ArrayList<Integer> [] toDelete = new ArrayList[3];
            for (int i=0;i<3;i++){
                toDelete[i]= new ArrayList<Integer>();
            }






            for (int check=0;check<2;check++){//2 passes to remove outliers with precision < threshold
                
                
                
                double [][][] Anew=new double[3][][];
                double [][] Rnew=new double [3][];

                for (int dim=0;dim<3;dim++){
                    
                    
                    //IJ.log("dim "+dim);
                    int thelength=totNumb-toDelete[dim].size();
                    
                    Anew[dim]=new double[thelength][bin-1];
                    
                    if (check==1){
                        //IJ.log("reject number="+toDelete[dim].size());
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
                                //IJ.log("residualOld ("+dim+") "+residual[i]);
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



                        }
                        catch(Exception eee){IJ.log("oops matrix non inversible -> maybe result will not be super precise");}
                    }



                }

                //IJ.log("drift after :");
                for (int u=0;u<driftZ.length;u++){
                    //IJ.log(""+driftX[u]+"  "+driftY[u]+"  "+driftZ[u]+"  ");
                }
                //IJ.log("");


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
        
        
        plot(driftIndex,driftX,li,lx,"drift of X","frame number","drift (nm)");
        
        plot(driftIndex,driftY,li,ly,"drift of Y","frame number","drift (nm)");
        
        plot(driftIndex,driftZ,li,lz,"drift of Z","frame number","drift (nm)");
        
        
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
    
    
    
    
    
    
    
    
    
    
    class LocalMaxima{
        float [][][] image;
        double x,y,z;//position to estimate -> drift
        int sizePatchX=11;//odd value;
        int sizePatchY=11;//odd value;
        int sizePatchZ=11;//odd value;
        LocalMaxima(float [][][] image){
            this.image=image;
            width=image.length;
            height=image[0].length;
            depth=image[0][0].length;
            
            sizePatchX=Math.min(sizePatchX, width-1);
            sizePatchY=Math.min(sizePatchY, height-1);
            sizePatchZ=Math.min(sizePatchZ, depth-1);
            if ((width%2==1)||(height%2==1)||(depth%2==1)){
                IJ.log("WARNING, drift image computation should have even size");
            }
            
            if (sizePatchX%2==0){
                sizePatchX=sizePatchX-1;
            }
            if (sizePatchY%2==0){
                sizePatchY=sizePatchY-1;
            }
            if (sizePatchZ%2==0){
                sizePatchZ=sizePatchZ-1;
            }
            
        }
        
        
        public void run(double threshold){//threshold=.75 -> 75% of data used to compute center
            
            
            double max=Double.MIN_VALUE;
            
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
            
            
            float [] vect=new float[sizePatchX*sizePatchY*sizePatchZ];
            int number=0;
            for (int i=0;i<sizePatchX;i++){
                for (int ii=0;ii<sizePatchY;ii++){
                    for (int iii=0;iii<sizePatchZ;iii++){
                        if ((i+xx-sizePatchX/2>=0)&&(ii+yy-sizePatchY/2>=0)&&(iii+zz-sizePatchZ/2>=0)&&(i+xx-sizePatchX/2<width)&&(ii+yy-sizePatchY/2<height)&&(iii+zz-sizePatchZ/2<depth)){
                            vect[number++]=image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
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
            for (int i=0;i<sizePatchX;i++){
                for (int ii=0;ii<sizePatchY;ii++){
                    for (int iii=0;iii<sizePatchZ;iii++){
                        if ((i+xx-sizePatchX/2>=0)&&(ii+yy-sizePatchY/2>=0)&&(iii+zz-sizePatchZ/2>=0)&&(i+xx-sizePatchX/2<width)&&(ii+yy-sizePatchY/2<height)&&(iii+zz-sizePatchZ/2<depth)){
                            if (image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2]>vect[(int)((number-1)*threshold)]){
                                sum+=image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                                x+=(i+xx-sizePatchX/2)*image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                                y+=(ii+yy-sizePatchY/2)*image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                                z+=(iii+zz-sizePatchZ/2)*image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
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
        
        
        //WARNING reversed because matrix order...
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
            double max=Double.MIN_VALUE;
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
            org.pasteur.imagej.utils.Matrixe m = new org.pasteur.imagej.utils.Matrixe(matVarCovar);
            try{
                m=org.pasteur.imagej.utils.Matrixe.inverse(m);
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
    

    class Init extends Thread{
        
        int idBin;
        Init(int idBin){
            this.idBin=idBin;
        }
        
        
        
        public void run(){
            
            ft.shift3D(image[idBin], step_real[idBin]);

            ft.fft(step_real[idBin], step_imag[idBin],output_real[idBin], output_imag[idBin]);



            for (int ii=0;ii<output_real[0].length;ii++){
                for (int iii=0;iii<output_real[0][0].length;iii++){
                    for (int iiii=0;iiii<output_real[0][0][0].length;iiii++){
                        step_real[idBin][ii][iii][iiii]=output_real[idBin][ii][iii][iiii]*gaussian[ii][iii][iiii];
                        step_imag[idBin][ii][iii][iiii]=output_imag[idBin][ii][iii][iiii]*gaussian[ii][iii][iiii];
                    }
                }
            }
        }
        
        
    }
    
    
    
    
    
    class Process extends Thread{
        int ip;
        int position1;
        int position2;
        int globalPosition;
        float [][][] res_A;
        float [][][] res_B;
        Process(int ip, int globalPosition,int position1, int position2){
            this.ip=ip;
            this.position1=position1;
            this.position2=position2;
            this.globalPosition=globalPosition;
            res_A = new float [sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
            res_B = new float [sizeFFTPaddedz][sizeFFTPaddedy][sizeFFTPaddedx];
        }
    
    
        public void run(){



            for (int ii=0;ii<output_real[0].length;ii++){
                for (int iii=0;iii<output_real[0][0].length;iii++){
                    for (int iiii=0;iiii<output_real[0][0][0].length;iiii++){
                        res_A[ii][iii][iiii]=step_real[position1][ii][iii][iiii]*step_real[position2][ii][iii][iiii]+step_imag[position1][ii][iii][iiii]*step_imag[position2][ii][iii][iiii];
                        res_B[ii][iii][iiii]=step_imag[position1][ii][iii][iiii]*step_real[position2][ii][iii][iiii]-step_real[position1][ii][iii][iiii]*step_imag[position2][ii][iii][iiii];
                    }
                }
            }
            
            
            ft.ifft(res_A, res_B,res_A, res_B);
            
            ft.shift3D(res_A, res_B);
            
            //here, result is in res_B
            
            
            //IJ.log("WARNING: needs synchronization of result... to do...!!!");
            int xxx,yyy,zzz,xxshift,yyshift,zzshift;
            synchronized(monitor) {
                for (xxx=0;xxx<sizeFFTx;xxx++){
                    for (yyy=0;yyy<sizeFFTy;yyy++){
                        for ( zzz=0;zzz<sizeFFTz;zzz++){
                            //xxshiftEven=((xxx+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;//*2 because even
                            //yyshiftEven=((yyy+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;
                            //zzshiftEven=((zzz+(sizeFFT/2)+sizeFFTPadded/2)%sizeFFTPadded)*2;

                            xxshift=(xxx+(sizeFFTx/2));//*2 because even
                            yyshift=(yyy+(sizeFFTy/2));
                            zzshift=(zzz+(sizeFFTz/2));

                            result[globalPosition][zzz][yyy][xxx]+=res_B[zzshift][yyshift][xxshift];
                        }
                    }
                }
            }
            
        }
        
    }

}
