/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;
import org.pasteur.imagej.process.gpu.DataPhase_;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Font;
import org.pasteur.imagej.cuda.MyCudaStream;
import ij.io.FileSaver;
import org.pasteur.imagej.data.PLocalization;
import org.pasteur.imagej.data.StackLocalization;
/**
 *
 * @author benoit
 */
public  class Model_generator_ {
    int width;
    int height;
    int model_number;
    DataPhase_ dp;
    double bckg=0;
    
    public Model_generator_(DataPhase_ dp){
        
        this.dp=dp;
        model_number=dp.param.PSF_number;
        
    }
    
    
    
    
    public void setImageSize(int imagesize){
        this.width=imagesize;
        this.height=imagesize;
    }
    
    public void setBackground(double bckg){
        this.bckg=bckg;
    }
    
    
    public void computePSF(double [][] list,String path_result){
        
        
        
        
        double [] x = new double [model_number];
        double [] y = new double [model_number];
        double [] z = new double [model_number];
        double [] f = new double [model_number];
        int [] frame = new int [model_number];
        double [] photon = new double [model_number];
        
        int maxFrame=0;
        for (int i=0;i<list.length;i++){
            int currentframe=(int)list[i][4];
            if (currentframe>maxFrame){
                maxFrame=currentframe;
            }
        }
        
        if (dp.param.sizeoutput>width){
            width=dp.param.sizeoutput;
        }
        if (dp.param.sizeoutput>height){
            height=dp.param.sizeoutput;
        }
        double [][][] res= new double [maxFrame+1][width][height];
        
        int dec=(dp.param.sizeoutput)/2;
        
        for (int i=0;i<=list.length/model_number;i++){
            
            for (int ii=0;ii<model_number;ii++){
                if (i*model_number+ii<list.length){
                    x[ii]=.001*(-(list[i*model_number+ii][0]%(1000.*dp.param.xystep)));
                    y[ii]=.001*(-(list[i*model_number+ii][1]%(1000.*dp.param.xystep)));
                    z[ii]=.001*(list[i*model_number+ii][2]);
                    f[ii]=.001*(list[i*model_number+ii][3]);
                    frame[ii]=(int)(list[i*model_number+ii][4]);
                    photon[ii]=(int)(list[i*model_number+ii][5]);
                }
            }
            dp.psf_fMany.computePSF(x, y, f, z);
            float [] im=dp.psf_fMany.getAllPSF();
            int xp,yp;
            for (int u=0;u<model_number;u++){
                if (i*model_number+u<list.length){
                    for (int uu=0;uu<dp.param.sizeoutput;uu++){
                        for (int uuu=0;uuu<dp.param.sizeoutput;uuu++){
                            xp=uu+(int)(list[i*model_number+u][0]/(1000.*dp.param.xystep))-dec;
                            yp=uuu+(int)(list[i*model_number+u][1]/(1000.*dp.param.xystep))-dec;
                            if ((xp>=0)&&(xp<this.width)&&(yp>=0)&&(yp<this.height)){
                                res[frame[u]][xp][yp]+=photon[u]*im[u*dp.param.sizeoutput*dp.param.sizeoutput+uu*dp.param.sizeoutput+uuu];
                            }
                        }
                    }
                }
            } 
        }
        
        for (int i=0;i<res.length;i++){
            for (int ii=0;ii<res[i].length;ii++){
                for (int iii=0;iii<res[i][ii].length;iii++){
                    res[i][ii][iii]+=bckg;
                }
            }
        }
        
        //path_result
        
        //org.pasteur.imagej.utils.ImageShow.imshow(res,"PSF generated");
         
        
        ImageStack ims = new ImageStack(res[0].length,res[0][0].length);
        for (int k=0;k<res.length;k++){
            ImageProcessor ip = new FloatProcessor(res[0].length,res[0][0].length);
            for (int i=0;i<res[0].length;i++){
                for (int ii=0;ii<res[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)res[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus("PSF generated",ims);
        
        try{
            FileSaver fs = new FileSaver(imp);
            fs.saveAsTiff(path_result);
            IJ.log("result image saved");
        }
        catch(Exception e){
            IJ.log("ERROR: impossible to save file");
        }
        //imp.show();
    }
    
    
    
    public void computePSFfromTable(StackLocalization sl,String path_result,boolean startAt0){
        computePSFfromTable(sl,path_result,startAt0,-1);
    }
    
    public void computePSFfromTable(StackLocalization sl,String path_result,boolean startAt0,int fixedimagesize){
        
        
        
        double minX=Double.POSITIVE_INFINITY,maxX=Double.NEGATIVE_INFINITY,minY=Double.POSITIVE_INFINITY,maxY=Double.NEGATIVE_INFINITY,minZ=Double.POSITIVE_INFINITY,maxZ=Double.NEGATIVE_INFINITY;
            
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                PLocalization p = sl.fl.get(i).loc.get(j);


                if (p.X<minX){
                    minX=p.X;
                }
                if (p.Y<minY){
                    minY=p.Y;
                }
                if (p.Z<minZ){
                    minZ=p.Z;
                }
                if (p.X>maxX){
                    maxX=p.X;
                }
                if (p.Y>maxY){
                    maxY=p.Y;
                }
                if (p.Z>maxZ){
                    maxZ=p.Z;
                }

            }
        }
        
        if (startAt0){
            IJ.log("image origin is (0,0)");
            minX=0;
            minY=0;
        }
        if (fixedimagesize==-1){
            width=(int)Math.ceil(((maxX-minX)/1000)/dp.param.xystep)+dp.param.sizeoutput;
            height=(int)Math.ceil(((maxY-minY)/1000)/dp.param.xystep)+dp.param.sizeoutput;
        }
        else{
            width=fixedimagesize;
            height=fixedimagesize;
        }
        
        IJ.log("min max "+minX+"  "+maxX+"  "+Math.ceil(((maxX-minX)/1000)/dp.param.xystep)+"  "+dp.param.xystep);
        IJ.log("wh "+width+"  "+height+"  "+dp.param.sizeoutput);
        
        
        
        
        
        
        int locNumber=0;
        int maxFrame=0;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                PLocalization p = sl.fl.get(i).loc.get(j);
                if (p.exists){
                    locNumber++;
                    if (p.frame>maxFrame){
                        maxFrame=p.frame;
                    }
                }
            }
        }
        
        
        
        int [][] index = new int[locNumber][];//index of localizations
        for (int i=0,u=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                PLocalization p = sl.fl.get(i).loc.get(j);
                if (p.exists){
                    int  [] v = new int [2];
                    v[0]=i;
                    v[1]=j;
                    index[u++]=v;
                }
            }
        }
        
        double [] x = new double [model_number];
        double [] y = new double [model_number];
        double [] z = new double [model_number];
        double [] f = new double [model_number];
        int [] frame = new int [model_number];
        double [] photon = new double [model_number];
        
        if (fixedimagesize>=0){
            if (dp.param.sizeoutput>fixedimagesize){

                IJ.log("WARNING: image size extended to psf output size");

            }
        }
        
        
        if (dp.param.sizeoutput>width){
            width=dp.param.sizeoutput;
        }
        if (dp.param.sizeoutput>height){
            height=dp.param.sizeoutput;
        }
        double [][][] res= new double [maxFrame+1][width][height];
        
        //int dec=(dp.param.sizeoutput)/2;
        
        for (int i=0;i<=locNumber/model_number;i++){
            IJ.log(""+(locNumber/model_number-i));
            for (int ii=0;ii<model_number;ii++){
                if (i*model_number+ii<locNumber){
                    PLocalization p=sl.fl.get(index[i*model_number+ii][0]).loc.get(index[i*model_number+ii][1]);
                    x[ii]=.001*(-((p.X-minX)%(1000.*dp.param.xystep)));
                    y[ii]=.001*(-((p.Y-minY)%(1000.*dp.param.xystep)));
                    z[ii]=.001*(p.Z);
                    f[ii]=(dp.param.Zfocus);
                    frame[ii]=(int)(p.frame);
                    photon[ii]=(int)(p.I);
                    
                }
            }
            dp.psf_fMany.computePSF(x, y, f, z);
            float [] im=dp.psf_fMany.getAllPSF();
            //dp.psf_fMany.imshow(height, im, "toto");
            int xp,yp;
            for (int u=0;u<model_number;u++){
                if (i*model_number+u<locNumber){
                    for (int uu=0;uu<dp.param.sizeoutput;uu++){
                        for (int uuu=0;uuu<dp.param.sizeoutput;uuu++){
                            PLocalization p=sl.fl.get(index[i*model_number+u][0]).loc.get(index[i*model_number+u][1]);
                            if (startAt0){
                                xp=uu+(int)((p.X-minX)/(1000.*dp.param.xystep)) - (dp.param.sizeoutput/2);
                                yp=uuu+(int)((p.Y-minY)/(1000.*dp.param.xystep)) - (dp.param.sizeoutput/2);
                            }
                            else{
                                xp=uu+(int)((p.X-minX)/(1000.*dp.param.xystep));
                                yp=uuu+(int)((p.Y-minY)/(1000.*dp.param.xystep));
                            }
                            if ((xp>=0)&&(xp<this.width)&&(yp>=0)&&(yp<this.height)){
                                res[frame[u]][xp][yp]+=photon[u]*im[u*dp.param.sizeoutput*dp.param.sizeoutput+uu*dp.param.sizeoutput+uuu];
                            }
                        }
                    }
                }
            } 
        }
        
        
        
        for (int i=0;i<res.length;i++){
            for (int ii=0;ii<res[i].length;ii++){
                for (int iii=0;iii<res[i][ii].length;iii++){
                    res[i][ii][iii]+=bckg;
                }
            }
        }
        
        //path_result
        
        //org.pasteur.imagej.utils.ImageShow.imshow(res,"PSF generated");
         
        
        ImageStack ims = new ImageStack(res[0].length,res[0][0].length);
        for (int k=0;k<res.length;k++){
            ImageProcessor ip = new FloatProcessor(res[0].length,res[0][0].length);
            for (int i=0;i<res[0].length;i++){
                for (int ii=0;ii<res[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)res[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus("PSF generated",ims);
        
        try{
            FileSaver fs = new FileSaver(imp);
            fs.saveAsTiff(path_result);
            IJ.log("result image saved");
        }
        catch(Exception e){
            IJ.log("ERROR: impossible to save file");
        }
        imp.show();


    }
        
        
        
        
    
    
    
    public void free(){
        this.dp.free();
    }
    
    
    
}
