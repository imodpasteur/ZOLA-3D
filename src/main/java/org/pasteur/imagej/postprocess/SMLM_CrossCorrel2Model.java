/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;
import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.data.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.process.ImageProcessor;
import java.util.Arrays;
import java.util.Comparator;
import java.util.ArrayList;
/**
 *
 * @author benoit
 */
public class SMLM_CrossCorrel2Model {


    StackLocalization sl;
    
    
    double sizePix=100;//nm
    
    
    double maxX=Double.NEGATIVE_INFINITY;
    double minX=Double.POSITIVE_INFINITY;
    double maxY=Double.NEGATIVE_INFINITY;
    double minY=Double.POSITIVE_INFINITY;
       
    
    FastFourierTransform fft;
    
    double [][] imageRendered;
    ArrayList<PLocalization> [][] imageLocalization;
    int width;//input image size
    int height;
    
    double [][] image;//image for fft
    double [][] imagevar;//image for fft
    double [][] model;//model for fft
    double modelvar;//model for fft
    int sizeX;//fft image size
    int sizeY;
    int shiftX;//because of 0 padding
    int shiftY;//because of 0 padding
    
    
    int patchsize;
    double modelradius;//radius of the model
    double modelprecision;//sigma gaussian
    double modelcenter;
    
    
    
    
    //Post processing that performs 2D cross correlation to detect objects
    public SMLM_CrossCorrel2Model(StackLocalization sl,double sizePix){
        
        this.sl=sl;
        
        this.sizePix=sizePix;
        
        
        maxX=Double.NEGATIVE_INFINITY;
        minX=Double.POSITIVE_INFINITY;
        maxY=Double.NEGATIVE_INFINITY;
        minY=Double.POSITIVE_INFINITY;
        for (int i=0;i<sl.fl.size();i++){
            for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                PLocalization p=sl.fl.get(i).loc.get(ii);
                if (p.X>maxX){
                    maxX=p.X;
                }
                if (p.X<minX){
                    minX=p.X;
                }
                if (p.Y>maxY){
                    maxY=p.Y;
                }
                if (p.Y<minY){
                    minY=p.Y;
                }
            }
        }
        
        
        this.imageRendering();
        
        
        //in case of 0 padding
        //sizeX=2*(int)Math.pow(2,(int)(Math.ceil(Math.log(width)/Math.log(2))));
        //sizeY=2*(int)Math.pow(2,(int)(Math.ceil(Math.log(height)/Math.log(2))));
        //shiftX=sizeX/4;
        //shiftY=sizeY/4;
        
        
        //NO padding
        sizeX=(int)Math.pow(2,(int)(Math.ceil(Math.log(width)/Math.log(2))));
        sizeY=(int)Math.pow(2,(int)(Math.ceil(Math.log(height)/Math.log(2))));
        shiftX=0;
        shiftY=0;
        
        IJ.log("size of FFT:  width:"+sizeX+"  height:"+sizeY);
        image = new double [sizeX][sizeY];
        for (int i=0,j=shiftX;i<width;i++,j++){
            for (int ii=0,jj=shiftY;ii<height;ii++,jj++){
                image[j][jj]=imageRendered[i][ii];
                
            }
        }
        
        fft=new FastFourierTransform(sizeX,sizeY);
        
        
    }
    
    
    
    
    
    
    
    //sigma corresponds to gaussian blur amount
    public void setModelAsRing(double diameter,double sigma){
        modelradius=diameter/2;//radius
        patchsize=(int)((diameter+3*sigma)/this.sizePix);
        this.modelprecision=sigma;
        modelcenter=(double)patchsize/2.;
        
        IJ.log("modelcenter "+modelcenter+"  ");
        model = new double[sizeX][sizeY];
        for (int i=shiftX;i<Math.min(patchsize+shiftX,sizeX);i++){
            for (int ii=shiftY;ii<Math.min(patchsize+shiftY,sizeY);ii++){
                //dist to ring center
                double d=sizePix*Math.sqrt((double)((i-modelcenter-shiftX)*(i-modelcenter-shiftX)+(ii-modelcenter-shiftX)*(ii-modelcenter-shiftX)));
                //dist to ring
                double d2=Math.abs(d-modelradius);
                model[i][ii]=Math.exp(-.5*d2*d2/(sigma*sigma));
                
            }
        }
        
        
        //normalize model (remove mean)
        double avg=0;
        for (int i=shiftX;i<Math.min(patchsize+shiftX,sizeX);i++){
            for (int ii=shiftY;ii<Math.min(patchsize+shiftY,sizeY);ii++){
                avg+=model[i][ii];
                
            }
        }
        avg/=(double)(Math.min(patchsize+shiftX,sizeX)*Math.min(patchsize+shiftY,sizeY));
        modelvar=0;
        for (int i=shiftX;i<Math.min(patchsize+shiftX,sizeX);i++){
            for (int ii=shiftY;ii<Math.min(patchsize+shiftY,sizeY);ii++){
                model[i][ii]-=avg;
                modelvar+=model[i][ii]*model[i][ii];
            }
        }
        modelvar/=(double)(Math.min(patchsize+shiftX,sizeX)*Math.min(patchsize+shiftY,sizeY));
    }
    
    
    //input: threshold cross correl amplitude
    public void detect(double threshold){
        
        
        
        
        double [][] crosscorr=crosscorrel();
        
        ArrayList<double []> al = computeLocalMax(crosscorr);
        
        
        //Lateral threshold
        double shiftThreshod=modelradius+modelprecision*2;
        
        
        
        
        //add new variable
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                PLocalization p = sl.fl.get(i).loc.get(j);
                ArrayList<String> variableNew = new ArrayList<String>();
                variableNew.add("pattern_index");
                variableNew.add("pattern_X");
                variableNew.add("pattern_Y");
                p.addListOfVariables(variableNew);
                p.otherVariable[0]=0;
            }
        }
        
        
        double [][] detections=new double [width][height];
        double maxValue=Double.NEGATIVE_INFINITY;
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                detections[i][ii]=imageRendered[i][ii];
                if (imageRendered[i][ii]>maxValue){
                    maxValue=imageRendered[i][ii];
                }
            }
        }
        
        
        int shiftsize=(int)Math.ceil(shiftThreshod/sizePix)+1;
        int indexDetection=1;
        for (int i=0;i<al.size();i++){
            double [] detect=al.get(i);
            if (detect[0]>threshold){
                int xx=(int)(detect[1]+detect[3]+modelcenter);
                int yy=(int)(detect[2]+detect[4]+modelcenter);
                //int xx=(int)(detect[1]+modelcenter);
                //int yy=(int)(detect[2]+modelcenter);
                double xxsub=(detect[1]+detect[3]+modelcenter)*sizePix+minX;
                double yysub=(detect[2]+detect[4]+modelcenter)*sizePix+minY;
                if ((xx>0)&&(xx<width)&&(yy>0)&&(yy<height)){
                    detections[xx][yy]=maxValue;
                    
                    for (int a =-shiftsize;a<shiftsize;a++){
                        for (int aa =-shiftsize;aa<shiftsize;aa++){
                            if ((xx+a>=0)&&(xx+a<width)&&(yy+aa>=0)&&(yy+aa<height)&&(imageLocalization[xx+a][yy+aa]!=null)){
                                for (int j=0;j<this.imageLocalization[xx+a][yy+aa].size();j++){
                                    
                                    PLocalization p =this.imageLocalization[xx+a][yy+aa].get(j);
                                    double d=Math.sqrt((p.X-xxsub)*(p.X-xxsub)+(p.Y-yysub)*(p.Y-yysub));
                                    if (d<shiftThreshod){
                                        p.otherVariable[0]=indexDetection;
                                        p.otherVariable[1]=xxsub;
                                        p.otherVariable[2]=yysub;
                                    }
                                }
                                
                            }
                        }
                    }
                    indexDetection++;
                }
            }
        }
        
        ImageShow.imshow(detections,"detections");
        
        IJ.log("detection finished");
        IJ.log("The localization table has been updated:");
        IJ.log("There are 3 new columns indicating the index and central position:  pattern_index ; pattern_X ; pattern_Y");

        
        
        
    }
    
    
    
    
    private ArrayList<double []> computeLocalMax(double [][] crosscorr){
        int decalFilter=3;
        //get local max:
        ArrayList<double []> al = new ArrayList<double []>();
        double [][] tmp = new double [width][height];
        for (int i=decalFilter;i<width-decalFilter;i++){
            for (int ii=decalFilter;ii<height-decalFilter;ii++){
                if ((crosscorr[i][ii]>0)){
                    tmp[i][ii]=Double.NEGATIVE_INFINITY;
                    for (int a=-decalFilter;a<=decalFilter;a++){
                        tmp[i][ii]=Math.max(crosscorr[i][ii+a],tmp[i][ii]);
                    }
                }
            }
        }
        double max=0;
        for (int i=decalFilter;i<width-decalFilter;i++){
            for (int ii=decalFilter;ii<height-decalFilter;ii++){
                if (crosscorr[i][ii]==tmp[i][ii]){
                    max=Double.NEGATIVE_INFINITY;
                    for (int a=-decalFilter;a<=decalFilter;a++){
                        max=Math.max(tmp[i+a][ii],max);
                    }
                    if (max==crosscorr[i][ii]){
                            
                            
                        double [] machin = new double[5];

                        double topX=0;
                        if ((i>1)&&(i<crosscorr.length-1)){
                            topX=topOfParabola(-1,crosscorr[i-1][ii],0,crosscorr[i][ii],1,crosscorr[i+1][ii]);
                        }

                        double topY=0;
                        if ((ii>1)&&(ii<crosscorr[0].length-1)){
                            topY=topOfParabola(-1,crosscorr[i][ii-1],0,crosscorr[i][ii],1,crosscorr[i][ii+1]);
                        }

                        machin[0]=crosscorr[i][ii];
                        machin[1]=i;
                        machin[2]=ii;
                        machin[3]=topX;
                        machin[4]=topY;
                        
                        al.add(machin);
                    }
                }
            }
        }
        return al;
    }
    
    
    
    
    
    
    
    private void imageRendering(){
        double binColor=255;
        width=(int)Math.ceil((maxX-minX)/sizePix);
        height=(int)Math.ceil((maxY-minY)/sizePix);
        
        
        double x;
        double y;
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        imageLocalization = new ArrayList[width][height];
        imageRendered = new double [width][height];
        int xx,yy;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    xx=(int)((x-minX)/sizePix);
                    yy=(int)((y-minY)/sizePix);
                    if (true){
                        //IJ.log("z "+z);
                        if ((xx>=0)&&(yy>=0)&&(xx<width)&&(yy<height)){
                            imageRendered[xx][yy]+= 1;
                            if (imageLocalization[xx][yy]==null){
                                imageLocalization[xx][yy]=new ArrayList<PLocalization>();
                            }
                            imageLocalization[xx][yy].add(sl.fl.get(i).loc.get(j));
                        }
                    }
                    
                }
            }
        }
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    private double topOfParabola(double xa,double ya,double xb,double yb,double xc,double yc){
        double a=(yc-ya)/((xc-xa)*(xc-xb))-(yb-ya)/((xb-xa)*(xc-xb));
        double b=((yb-ya)/(xb-xa))-a*(xb+xa);
        return(-b/(2*a));
    }
    
    
    
    
    private void computeNormalizedVarImage(){
        
        double [][] meanImage=new double [sizeX][sizeY];
        double [][] sumImage=new double [sizeX][sizeY];
        double [][] sumSquareImage=new double [sizeX][sizeY];
        double [][] sumImageSecond=new double [sizeX][sizeY];
        double [][] normImage=new double [sizeX][sizeY];
        double [][] sumSquareImageSecond=new double [sizeX][sizeY];
        
        
        
        
        
        for (int i=0;i<sizeX;i++){
            for (int ii=0;ii<sizeY;ii++){
                int valx1=patchsize/2;
                int valx2=patchsize/2;
                if (i<patchsize/2-1){
                    valx1=i+1;
                }
                if (i>sizeX-patchsize/2-2){
                    valx2=sizeY-i-1;
                }
                int valy1=patchsize/2;
                int valy2=patchsize/2;
                if (ii<patchsize/2-1){
                    valy1=ii+1;
                }
                if (ii>sizeX-patchsize/2-2){
                    valy2=sizeY-ii-2;
                }
                normImage[i][ii]=(valx1+valx2)*(valy1+valy2);
            }
        }
        
        for (int i=0;i<sizeX;i++){
            sumImage[i][0]=0;
            sumSquareImage[i][0]=0;
            for (int t=0;t<patchsize/2+1;t++){//before 0 it is outside
                if (t<sizeY){
                    sumImage[i][0]+=image[i][t];
                    sumSquareImage[i][0]+=image[i][t]*image[i][t];
                }
            }
        
            for (int ii=1;ii<sizeY;ii++){
                sumImage[i][ii]=sumImage[i][ii-1];
                if (ii-patchsize/2>=0){
                    sumImage[i][ii]-=image[i][ii-patchsize/2];
                }
                else{
                    //outside
                }
                if (ii+patchsize/2<sizeY){
                    sumImage[i][ii]+=image[i][ii+patchsize/2];
                }
                else{
                    //outside
                }
                sumSquareImage[i][ii]=sumSquareImage[i][ii-1];
                if (ii-patchsize/2>=0){
                    sumSquareImage[i][ii]-=image[i][ii-patchsize/2]*image[i][ii-patchsize/2];
                }
                if (ii+patchsize/2<sizeY){
                    sumSquareImage[i][ii]+=image[i][ii+patchsize/2]*image[i][ii+patchsize/2];
                }
                
            }
        }
        
        for (int ii=0;ii<sizeY;ii++){
            sumImageSecond[0][ii]=0;
            sumSquareImageSecond[0][ii]=0;
            for (int t=0;t<patchsize/2+1;t++){
                if (t<sizeX){
                    sumImageSecond[0][ii]+=sumImage[t][ii];
                    sumSquareImageSecond[0][ii]+=sumSquareImage[t][ii];
                }
            }
            for (int i=1;i<sizeX;i++){
                sumImageSecond[i][ii]=sumImageSecond[(i-1)][ii];
                if (i-patchsize/2>=0){
                    sumImageSecond[i][ii]-=sumImage[(i-patchsize/2)][ii];
                }
                if (i+patchsize/2<sizeX){
                    sumImageSecond[i][ii]+=sumImage[(i+patchsize/2)][ii];
                }
                
                sumSquareImageSecond[i][ii]=sumSquareImageSecond[(i-1)][ii];
                if (i-patchsize/2>=0){
                    sumSquareImageSecond[i][ii]-=sumSquareImage[(i-patchsize/2)][ii];
                }
                if (i+patchsize/2<sizeX){
                    sumSquareImageSecond[i][ii]+=sumSquareImage[(i+patchsize/2)][ii];
                }
                
            }
            for (int i=0;i<sizeX;i++){
                meanImage[i][ii]=sumImageSecond[i][ii]/normImage[i][ii];
            }
        }
        
        
        
        for (int i=0;i<sizeX;i++){
            for (int ii=0;ii<sizeY;ii++){
                imagevar[i][ii]=(sumSquareImageSecond[i][ii]/normImage[i][ii])-meanImage[i][ii]*meanImage[i][ii];
            }
        }
        
        //ImageShow.imshow(imagevar,"var");
        
        //ImageShow.imshow(meanImage,"mean");
        
        
    }
    
    
    
    
    private double [][] crosscorrel(){
        imagevar=new double [sizeX][sizeY];
        for (int i=0;i<sizeX;i++){
            for (int ii=0;ii<sizeY;ii++){
                imagevar[i][ii]=1;
            }
        }
        //computeNormalizedVarImage();
        
        double [][] tmp=new double [sizeX][sizeY];
        double [][] areal;
        double [][] aimag;
        double [][] breal;
        double [][] bimag;
        double [][] tmpreal=new double [sizeX][sizeY];
        double [][] tmpimag=new double [sizeX][sizeY];
        //fft input image
        fft.setReal(image);
        fft.setImag(tmp);
        fft.fft2D();
        areal=fft.getPointerRealOut2D();
        aimag=fft.getPointerImagOut2D();
        
        for (int i=0;i<sizeX;i++){
            for (int ii=0;ii<sizeY;ii++){
                //complexeConjugate
                tmpreal[i][ii]=areal[i][ii];
                tmpimag[i][ii]=aimag[i][ii];
            }
        }
        
        
        //fft model
        fft.setReal(model);
        fft.setImag(tmp);
        fft.fft2D();
        breal=fft.getPointerRealOut2D();
        bimag=fft.getPointerImagOut2D();
        
        
        for (int i=0;i<sizeX;i++){
            for (int ii=0;ii<sizeY;ii++){
                //complexeConjugate
                double re=tmpreal[i][ii];
                double im=tmpimag[i][ii];
                tmpreal[i][ii]=im*bimag[i][ii]+re*breal[i][ii];
                tmpimag[i][ii]=im*breal[i][ii]-re*bimag[i][ii];
            }
        }
        
        
        
        fft.setReal(tmpreal);
        fft.setImag(tmpimag);
        fft.ifft2D();
        areal=fft.getPointerRealOut2D();
        aimag=fft.getPointerImagOut2D();
        
        
       
        
        
        
        for (int i=0;i<sizeX;i++){
            for (int ii=0;ii<sizeY;ii++){
                //tmpreal[i][ii]=areal[i][ii];
                try{
                    if (modelvar*imagevar[i+shiftX-patchsize/2+1][ii+shiftY-patchsize/2+1]>0.0000000){

                        tmpreal[i][ii]=areal[i][ii]/Math.sqrt(modelvar*imagevar[i+shiftX-patchsize/2+1][ii+shiftY-patchsize/2+1]);

                    }
                    else {
                        tmpreal[i][ii]=-1;
                    }
                    
                
                }catch(Exception et){}
            }
        }
        
        
        return tmpreal;
        
    }
    
    
    
}
