/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process;


import java.util.ArrayList;
import org.pasteur.imagej.utils.SCMOScamera;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.utils.PolynomialFit;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import java.util.Arrays;

/**
 *
 * @author benoit
 */
public class InitBackgroundAndPhotonNumber {
    
    //compute photon number and apply 0 padding
    
    double [][] polyBackground;
    
    int orderPolyBackground=2;//plan
    
    public double [] B;
    public double [][] A;
    
    double exclusionRadius=20;//in pixels
    
    double backgroundRadius=50;//in pixels
    
    int sizeSubImage=64;
    
    public double [][][][] image;
    public double [][][] scmos;
    //private double [][] scmosImage;
    int width;
    int height;
    
    double filtSizeMicron=0.5;
    
    int sizeFilt;
    
    
    public InitBackgroundAndPhotonNumber(double [][][] image,float [] xpoints, float [] ypoints,int sizePatch,double adu, double gain,double offset,double zstep){
        
        
        width=image[0].length;
        height=image[0][0].length;
        
        double slope=adu/gain;
        double intercept=-offset*adu/gain;
        for (int z=0;z<image.length;z++){
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    
                    image[z][i][ii]=(image[z][i][ii]*slope+intercept)+1;
                    
                    if (image[z][i][ii]<0){
                        image[z][i][ii]=0;
                    }
                    
                    
                    
                }
            }
            
            
        }
        initBackgroundAndPhotonNumberMany(image,xpoints,ypoints,sizePatch,zstep,null);
        
    }
    
    public InitBackgroundAndPhotonNumber(double [][][] image,float [] xpoints, float [] ypoints,int sizePatch,SCMOScamera scmoscam,double zstep){
        
        
        width=image[0].length;
        height=image[0][0].length;
        
        for (int z=0;z<image.length;z++){
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    
                    image[z][i][ii]=((image[z][i][ii]-scmoscam.scmosoffset[i][ii])/scmoscam.scmosgain[i][ii]);
                    
                    if (image[z][i][ii]<0){
                        image[z][i][ii]=0;
                    }
                }
            }
            
        }
        initBackgroundAndPhotonNumberMany(image,xpoints,ypoints,sizePatch,zstep,scmoscam.scmosvargain);
    }
    
    void initBackgroundAndPhotonNumberMany(double [][][] image,float [] xpoints, float [] ypoints,int sizePatch,double zstep,double [][] scmos){
        
        
        
        //sizeFilt is 0.5um :
        sizeFilt=Math.max((int)(filtSizeMicron/zstep)/2,1);
            
            
        //IJ.log("slope inter "+slope+"  "+intercept);
        
        
        refinePoints(image,xpoints,ypoints,sizePatch);
        
        
        double tailleFiltre=1;
        
        
        sizeSubImage=sizePatch;
        
        exclusionRadius=3*sizePatch/4;
        backgroundRadius=sizePatch;
        
        width=image[0].length;
        height=image[0][0].length;
        int sizeImage=Math.max(width, height);
        //GaussianKernel gk = new GaussianKernel(sizeImage,tailleFiltre,0);
        
        
        
        //ImageShow.imshow(image,"imimim");
        
        computeBackground(image);
        
        
        
        
        int nbpoints=xpoints.length;
        
        if (nbpoints<=0){
            IJ.log("no points...Please add seeds on PSFs");
        }
        else{
            //IJ.log("nb points : "+nbpoints);
        }
        
        
        
        
        /*
        ArrayList<Integer> [][] backgroundPosition = new ArrayList[nbpoints][2];
        
        
        for (int p=0;p<nbpoints;p++){
            backgroundPosition[p][0]=new ArrayList<Integer>();
            backgroundPosition[p][1]=new ArrayList<Integer>();
        }
        
        
        
        int [] count = new int [nbpoints];
        
        for (int ii=0;ii<width;ii++){
            for (int iii=0;iii<height;iii++){
                for (int p=0;p<nbpoints;p++){
                    double dist=Math.sqrt((ii-xpoints[p])*(ii-xpoints[p])+(iii-ypoints[p])*(iii-ypoints[p]));
                    if ((dist<backgroundRadius)&&(dist>exclusionRadius)){
                        backgroundPosition[p][0].add(ii);
                        backgroundPosition[p][1].add(iii);
                        count[p]++;
                        
                    }
//                    if ((dist<backgroundRadius)){
//                        backgroundAndPSFPosition[p][0].add(ii);
//                        backgroundAndPSFPosition[p][1].add(iii);
//                    }
                }
            }
        }
        double [][] vectBackground = new double[nbpoints][];
        for (int p=0;p<nbpoints;p++){
            vectBackground[p] = new double[count[p]*image.length];
        }
        
        for (int p=0;p<nbpoints;p++){
            int k=0;
            for (int z=0;z<image.length;z++){
                for (int i=0;i<backgroundPosition[p][0].size();i++){
                    vectBackground[p][k++]=image[z][(int)backgroundPosition[p][0].get(i)][(int)backgroundPosition[p][1].get(i)];
                }
            }
        }
        
        double [] medianBackground = new double[nbpoints];
                
                
        for (int p=0;p<nbpoints;p++){
            Arrays.sort(vectBackground[p]);
            medianBackground[p]=2+vectBackground[p][vectBackground[p].length/2];//+2 to deal with low values with asymetric function where median is not good // and also gaussian additive noise
        
        }
        
        
        
        
        
        //compute mean according to a threshold based on median valu : non noisy solution
        double [] mean = new double [nbpoints];
        double [] count2 = new double [nbpoints];
        for (int p=0;p<nbpoints;p++){
            for (int u=0;u<vectBackground[p].length;u++){
                if ((vectBackground[p][u]>medianBackground[p]-3*Math.sqrt(medianBackground[p]))&&(vectBackground[p][u]<medianBackground[p]+3*Math.sqrt(medianBackground[p]))){
                    mean[p]+=vectBackground[p][u];
                    count2[p]++;
                }
            }
        }
        
        //IJ.log("meanback+std: "+meanback+"       median: "+medianBackground);
        
        B = new double [nbpoints];
        for (int p=0;p<nbpoints;p++){
            this.B[p]=mean[p]/count2[p];
            IJ.log("background: "+B[p]);
        }
        */
        
        //IJ.log("mean Back="+newMeanBcg);
        
        //ImageShow.imshow(image,"image meanbckg"+newMeanBcg);
        
        //IJ.log("mean background : "+B);
        
        
        
        
        
        
        
        
        
        
        
        
        this.B=new double [nbpoints];
        for (int p=0;p<nbpoints;p++){
            this.B[p]=0;
            double countB=0;
            //use the square
            for (int u=0;u<sizeSubImage;u++){
                for (int uu=0;uu<sizeSubImage;uu++){
                    if (((int)(xpoints[p])-(sizeSubImage/2)+u>=0)&&((int)(xpoints[p])-(sizeSubImage/2)+u<image[0].length)&&((int)(ypoints[p])-(sizeSubImage/2)+uu>=0)&&((int)(ypoints[p])-(sizeSubImage/2)+uu<image[0][0].length)){
                        B[p]+=this.polyBackground[(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu];
                        countB++;
                    }
                }
            }
            
            if (countB>=1){
                B[p]/=countB;
            }
            //IJ.log("Background "+B[p]);
        }
        
        
        
        
        
        
        
        double [][] A ;
        double [][] xAxis;
        A = new double [1][image.length*nbpoints];
        double [] res = new double [image.length*nbpoints];
        xAxis = new double [1][image.length*nbpoints];
        this.A=new double [nbpoints][image.length];
        
        double [] sumAPerPoint= new double[nbpoints];
        double maxphoton=0;
        for (int p=0;p<nbpoints;p++){
            sumAPerPoint[p]=0;
            for (int z=0;z<image.length;z++){
                A[0][z*nbpoints+p]=0;
                xAxis[0][z*nbpoints+p]=z;
                
                //use the square
                for (int u=0;u<sizeSubImage;u++){
                    for (int uu=0;uu<sizeSubImage;uu++){
                        if (((int)(xpoints[p])-(sizeSubImage/2)+u>=0)&&((int)(xpoints[p])-(sizeSubImage/2)+u<image[z].length)&&((int)(ypoints[p])-(sizeSubImage/2)+uu>=0)&&((int)(ypoints[p])-(sizeSubImage/2)+uu<image[z][0].length)){
                            A[0][z*nbpoints+p]+=image[z][(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu]-this.B[p];
                            sumAPerPoint[p]+=image[z][(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu]-this.B[p];
                            
                            //we do not remove background yet
                            //A[0][z*nbpoints+p]+=image[z][(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu];
                            //sumAPerPoint[p]+=image[z][(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu];
                            
                            
                        }
                    }
                }
                
//                for (int n=0;n<backgroundAndPSFPosition[p][0].size();n++){
//                    A[0][z*nbpoints+p]+=image[z][backgroundAndPSFPosition[p][0].get(n)][backgroundAndPSFPosition[p][1].get(n)]-this.B;
//                    sumAPerPoint[p]+=image[z][backgroundAndPSFPosition[p][0].get(n)][backgroundAndPSFPosition[p][1].get(n)]-this.B;
//                    
//                }
                //IJ.log("A "+p+"  "+z+"  "+A[0][z*nbpoints+p]+"   "+backgroundAndPSFPosition[p][0].size());
            }
            
            sumAPerPoint[p]/=(double)image.length;
            if (sumAPerPoint[p]>maxphoton){
                maxphoton=sumAPerPoint[p];
            }
        }
        
        
        
        
        if (true){
            //////////********************************************************BEGIN
            //mean filter of curves
            
            
            

            for (int p=0;p<nbpoints;p++){
                for (int z=0;z<image.length;z++){
                    double meanA=0;
                    double countma=0;
                    for (int f=-sizeFilt;f<=sizeFilt;f++){
                        if ((z+f>=0)&&(z+f<image.length)){
                            meanA+=A[0][(z+f)*nbpoints+p];
                            countma++;
                        }
                    }
                    
                    if (countma>=1){
                        meanA/=countma;
                    }
                    res[z*nbpoints+p]=meanA;
                    this.A[p][z]=meanA;
                    
                }
            }


            //////////********************************************************END
        
        
        
        }
        else{


            


            //////////********************************************************BEGIN
            //polynomial fit
            int order=10;
            //IJ.log("order="+order);
            //normalize
            for (int p=0;p<nbpoints;p++){
                for (int z=0;z<image.length;z++){
                    A[0][z*nbpoints+p]-=sumAPerPoint[p];
                }
            }





            PolynomialFit pf = new PolynomialFit(order,A,xAxis);
            pf.run();
            double [] v = new double [1];

            for (int p=0;p<nbpoints;p++){
                for (int z=0;z<image.length;z++){

                    v[0]=xAxis[0][z*nbpoints+p];
                    v=pf.transform(v);

                    this.A[p][z]=v[0]+sumAPerPoint[p]-B[p];


                    //same photon number for plot
                    res[z*nbpoints+p]=v[0]+maxphoton;
                    A[0][z*nbpoints+p]+=maxphoton;
                }
            }
            //////////********************************************************END
            
        }
        
        //this.plot(xAxis[0], A[0], res,"Photon number estimation","z position (nm)","photon number");
        
//        for (int i=0;i<this.A.length;i++){
//            for (int ii=0;ii<this.A[0].length;ii++){
//                IJ.log("A "+this.A[i][ii]);
//            }
//        }
        
        this.image=new double[nbpoints][image.length][this.sizeSubImage][this.sizeSubImage];
        
        double [][] sumB = new double [nbpoints][image.length];
        for (int p=0;p<nbpoints;p++){
            for (int z=0;z<image.length;z++){
                for (int u=0;u<sizeSubImage;u++){
                    for (int uu=0;uu<sizeSubImage;uu++){
                        if (((int)(xpoints[p])-(sizeSubImage/2)+u>=0)&&((int)(xpoints[p])-(sizeSubImage/2)+u<image[z].length)&&((int)(ypoints[p])-(sizeSubImage/2)+uu>=0)&&((int)(ypoints[p])-(sizeSubImage/2)+uu<image[z][0].length)){
                            if (scmos!=null){
                                this.image[p][z][u][uu]=image[z][(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu]+scmos[(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu];
                            }
                            else{
                                this.image[p][z][u][uu]=image[z][(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu];

                            }
                        }
                        else{
                            this.image[p][z][u][uu]=B[p];
                        }
                    }
                }
            }
            
            if (scmos!=null){
                this.scmos=new double[nbpoints][this.sizeSubImage][this.sizeSubImage];
                for (int u=0;u<sizeSubImage;u++){
                    for (int uu=0;uu<sizeSubImage;uu++){
                        if (((int)(xpoints[p])-(sizeSubImage/2)+u>=0)&&((int)(xpoints[p])-(sizeSubImage/2)+u<image[0].length)&&((int)(ypoints[p])-(sizeSubImage/2)+uu>=0)&&((int)(ypoints[p])-(sizeSubImage/2)+uu<image[0][0].length)){
                            this.scmos[p][u][uu]=scmos[(int)(xpoints[p])-(sizeSubImage/2)+u][(int)(ypoints[p])-(sizeSubImage/2)+uu];

                        }
                        else{
                            this.scmos[p][u][uu]=0;
                        }
                    }
                }
            }
            else{
                this.scmos=new double[nbpoints][][];
                for (int i=0;i<nbpoints;i++){
                    this.scmos[i]=null;
                }
            }
            double themin=Double.POSITIVE_INFINITY;
            double themax=Double.NEGATIVE_INFINITY;
            for (int i=0;i<this.image[p].length;i++){
                for (int ii=0;ii<this.image[p][i].length;ii++){
                    for (int iii=0;iii<this.image[p][i][ii].length;iii++){
                        if (themin>this.image[p][i][ii][iii]){
                            themin=this.image[p][i][ii][iii];
                        }
                        if (themax<this.image[p][i][ii][iii]){
                            themax=this.image[p][i][ii][iii];
                        }
                    }
                }
            }
            
        }
        
        
        
        
    }
    
    
    
    
    
    public void plot(double [] x, double [] y,double [] yline,String title,String xlabel,String ylabel){
        if (x.length>1){
            Plot p = new Plot(""+title,xlabel,ylabel);
            p.setFont(0, 18);
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
                if (yline[i]<yMin){
                    yMin=yline[i];
                }
                if (yline[i]>yMax){
                    yMax=yline[i];
                }
            }
            p.setLimits(xMin, xMax, yMin, yMax);

            for (int ii=0;ii<x.length;ii++){

                p.setColor(Color.blue);
                p.setLineWidth(2);
                p.add("CIRCLE",  x,y);
                p.setColor(Color.red);
                p.setLineWidth(2);
                p.add("CIRCLE",  x,yline);
                //p.show();



            }
            p.show();
        }
        
    }
    
    
    public void plot(double [] x, double [] y,String title,String xlabel,String ylabel){
        if (x.length>1){
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
            p.setLimits(xMin, xMax, yMin, yMax);
            p.setFont(0, 18);

            for (int ii=0;ii<x.length;ii++){

                p.setColor(Color.blue);
                p.setLineWidth(2);
                p.add("CIRCLE",  x,y);
                //p.show();



            }

            p.show();
        }
        
    }
    
    
    
    
    
    
    void refinePoints(double [][][] image,float [] xpoints, float [] ypoints,int sizePatch){
        
        for (int p=0;p<xpoints.length;p++){
            double [][] subImage=new double [sizePatch][sizePatch];
            double [] vect=new double [sizePatch*sizePatch];
            
            for (int i=0;i<image.length;i++){
                for (int u=0;u<sizePatch;u++){
                    for (int uu=0;uu<sizePatch;uu++){
                        if (((int)(xpoints[p])-(sizePatch/2)+u>=0)&&((int)(xpoints[p])-(sizePatch/2)+u<image[i].length)&&((int)(ypoints[p])-(sizePatch/2)+uu>=0)&&((int)(ypoints[p])-(sizePatch/2)+uu<image[i][0].length)){
                            subImage[u][uu]+=image[i][(int)(xpoints[p])-(sizePatch/2)+u][(int)(ypoints[p])-(sizePatch/2)+uu];
                            
                        }
                    }
                }
            }
            
            for (int u=0,t=0;u<sizePatch;u++){
                for (int uu=0;uu<sizePatch;uu++){
                    vect[t++]=subImage[u][uu];
                }
            }
            Arrays.sort(vect);
            
            double meanX=0;
            double meanY=0;
            double sum=0;
            int count=0;
            for (int u=sizePatch/4;u<3*sizePatch/4;u++){
                for (int uu=sizePatch/4;uu<3*sizePatch/4;uu++){
                    if (subImage[u][uu]>vect[(int)((vect.length-1)*.9)]){
                        meanX+=((int)(xpoints[p])-(sizePatch/2)+u)*subImage[u][uu];
                        meanY+=((int)(ypoints[p])-(sizePatch/2)+uu)*subImage[u][uu];
                        sum+=subImage[u][uu];
                        count++;
                    }
                }
            }
            //ImageShow.imshow(subImage);
            //IJ.log("before "+xpoints[p]+"   "+ypoints[p]);
            if (count>=1){
                
                if (Math.sqrt((xpoints[p]-(float)(meanX/sum))*(xpoints[p]-(float)(meanX/sum)))<sizePatch/4){
                    
                    xpoints[p]=(float)(meanX/sum);
                }
                if (Math.sqrt((ypoints[p]-(float)(meanY/sum))*(ypoints[p]-(float)(meanY/sum)))<sizePatch/4){
                    ypoints[p]=(float)(meanY/sum);
                }
            }
            //IJ.log("after  "+xpoints[p]+"   "+ypoints[p]);
        }
    }
    
    
    
    
    
    public void computeBackground(double [][][] image){
        
        polyBackground=new double [image[0].length][image[0][0].length];
        
        
        //double [][][] imtmp=new double  [image.length][image[0].length][image[0][0].length];
        
        boolean [][][] is2=new boolean  [image.length][image[0].length][image[0][0].length];
        boolean [][][] is3=new boolean  [image.length][image[0].length][image[0][0].length];
        
        
        //first pass...suppress high intensity
        
        
        double mean0=0;
        double std0=0;
        int i,ii,iii,a,aa,aaa,t;
        
        for (i=0,t=0;i<image.length;i++){
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    
                    
                    
                    mean0+=image[i][ii][iii];
                }
            }
        }
        
        
        mean0/=image.length*image[0].length*image[0][0].length;
        for (i=0;i<image.length;i++){
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    std0+=(image[i][ii][iii]-mean0)*(image[i][ii][iii]-mean0);
                }
            }
        }
        std0/=image.length*image[0].length*image[0][0].length;
        std0=Math.sqrt(std0);
        
        
        //second pass...compute mean background
        
        double mean=0;
        double std=0;
        int countMean0=0;
        for (i=0,t=0;i<image.length;i++){
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    if (image[i][ii][iii]<=mean0+2*std0){
                        mean+=image[i][ii][iii];
                        countMean0++;
                    }
                }
            }
        }
        
        
        mean/=countMean0;
        for (i=0;i<image.length;i++){
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    if (image[i][ii][iii]<=mean0+3*std0){
                        std+=(image[i][ii][iii]-mean)*(image[i][ii][iii]-mean);
                    }
                }
            }
        }
        std/=countMean0;
        std=Math.sqrt(std);
        
        
        for (i=0;i<image.length;i++){
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    if (image[i][ii][iii]>mean+3*std){
                        is2[i][ii][iii]=false;
                        //imtmp[i][ii][iii]=-1;
                    }
                    else{
                        //imtmp[i][ii][iii]=image[i][ii][iii];
                        is2[i][ii][iii]=true;
                    }
                }
            }
        }
        
        //ImageShow.imshow(imtmp,"image");
        
        
        double [] meanStack = new double [image.length];
        double [] meanStackXAxis = new double [image.length];
        
        for (i=0;i<image.length;i++){
            meanStack[i]=0;
            meanStackXAxis[i]=i;
            int count=0;
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    if (is2[i][ii][iii]){
                        meanStack[i]+=image[i][ii][iii];
                        count++;
                    }
                }
            }
            if (count!=0){
                meanStack[i]/=count;
            }
        }
        
        
        double [] res = new double [image.length];
        double [] vu = new double [1];
        
        
        
        //this.plot(meanStackXAxis, meanStack,"mean background","z position (nm)","photon number");
        
        int miniposition=0;
        
        double [] filtMean=new double [meanStack.length];
        
        for (i=0;i<filtMean.length;i++){
            filtMean[i]=0;
            int count=0;
            for (int u=-sizeFilt;u<=sizeFilt;u++){
                if ((i+u>=0)&&(i+u<filtMean.length)){
                    filtMean[i]+=meanStack[i+u];
                    count++;
                }
            }
            filtMean[i]/=count;
            if (filtMean[i]<filtMean[miniposition]){
                miniposition=i;
            }
        }
        
        
        //this.plot(meanStackXAxis, meanStack,filtMean,"mean background","z position (nm)","photon number");
        
        
        
        //IJ.log(""+(-p_b/(2*p_a))+"  "+p_b+"  "+p_a);
        int numberImageNeighbor=sizeFilt;
        
        int totalImageNumberProcess=Math.min(numberImageNeighbor*2+1, image.length);
        
        int sizeDilate=2;
        
        if ((miniposition<numberImageNeighbor)){
            miniposition=numberImageNeighbor;
        }
        if ((miniposition>=image.length-numberImageNeighbor)){
            miniposition=image.length-numberImageNeighbor-1;
        }
        
        //IJ.log("miniposition "+miniposition);
        
        for (i=0;i<image.length;i++){
            for (ii=0;ii<image[0].length;ii++){
                for (iii=0;iii<image[0][0].length;iii++){
                    is3[i][ii][iii]=true;
                    if (i<miniposition-(numberImageNeighbor)){
                        is3[i][ii][iii]=false;
                        //imtmp[i][ii][iii]=-1;
                    }
                    else if (i>miniposition+(numberImageNeighbor)){
                        is3[i][ii][iii]=false;
                        //imtmp[i][ii][iii]=-1;
                    }
                    else{

                        theloop:for (int u=-sizeDilate;u<=sizeDilate;u++){
                            for (int uu=-sizeDilate;uu<=sizeDilate;uu++){
                                for (int uuu=-sizeDilate;uuu<=sizeDilate;uuu++){
                                    if ((i+u>=0)&&(i+u<image.length)&&(ii+uu>=0)&&(ii+uu<image[0].length)&&(iii+uuu>=0)&&(iii+uuu<image[0][0].length)){
                                        if (!is2[i+u][ii+uu][iii+uuu]){
                                            is3[i][ii][iii]=false;
                                            //imtmp[i][ii][iii]=-1;
                                            break theloop;
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }
        
        
        //ImageShow.imshow(imtmp,"image");
        
        
        
        
        
        double meantmp=0;
        int count=0;
        //double [] vect= new double [image.length];
        ArrayList<double []> al = new ArrayList<double []>();
        for (ii=0;ii<image[0].length;ii++){
            for (iii=0;iii<image[0][0].length;iii++){
                meantmp=0;
                count=0;
                for (i=0;i<image.length;i++){
                    if (is3[i][ii][iii]){
                        count++;
                    }
                    
                }
                if (count>=1){
                    
                    for (i=0;i<image.length;i++){
                        if (is3[i][ii][iii]){
                            meantmp+=image[i][ii][iii];
                            //vect[i]=image[i][ii][iii];
                        }
                        else{
                            //vect[i]=-1;
                        }
                    }
                    //Arrays.sort(vect);
                    
                    meantmp/=count;
                    double [] data=new double [3];
                    data[0]=ii;
                    data[1]=iii;
                    data[2]=meantmp;
                    //data[2]=vect[vect.length/2+(vect.length-count)/2];
                    al.add(data);
                }
            }
        }
        
        
        
        
        
        int dim=2;
        
        double [][] x=new double [dim][al.size()];
        double [][] y=new double [dim][al.size()];
        double [] truc;
        
        
        for (t=0;t<al.size();t++){
            truc=al.get(t);
            x[0][t]=truc[0];
            x[1][t]=truc[1];
            y[0][t]=truc[2];
            y[1][t]=truc[2];
        }
        
        
        
        PolynomialFit pf = new PolynomialFit(orderPolyBackground,y,x);
        
        
        pf.run();
        
        
        double [] v = new double [dim];
        for (ii=0;ii<image[0].length;ii++){
            for (iii=0;iii<image[0][0].length;iii++){
                v[0]=ii;
                v[1]=iii;
                v=pf.transform(v);
                polyBackground[ii][iii]=v[0];
            }
        }
        //ImageShow.imshow(polyBackground,"background");
    }
    
    
    
}
