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
/**
 *
 * @author benoit
 */
public  class Model_generator_ {
    int imageSize;
    int model_number;
    
    DataPhase_ dp;
    
    
    public Model_generator_(DataPhase_ dp,int imageSize){
        
        this.dp=dp;
        model_number=dp.param.PSF_number;
        this.imageSize=imageSize;
    }
    
    
    
    
    public void computePSF(double [][] list,String path_result){
        
        double [] x = new double [model_number];
        double [] y = new double [model_number];
        double [] z = new double [model_number];
        double [] f = new double [model_number];
        
        if (dp.param.sizeoutput>imageSize){
            imageSize=dp.param.sizeoutput;
        }
        double [][][] res= new double [list.length][imageSize][imageSize];
        
        int dec=(dp.param.sizeoutput)/2;
        
        for (int i=0;i<=list.length/model_number;i++){
            
            for (int ii=0;ii<model_number;ii++){
                if (i*model_number+ii<list.length){
                    x[ii]=.001*(-(list[i*model_number+ii][0]%(1000.*dp.param.xystep)));
                    y[ii]=.001*(-(list[i*model_number+ii][1]%(1000.*dp.param.xystep)));
                    z[ii]=.001*(list[i*model_number+ii][2]);
                    f[ii]=.001*(list[i*model_number+ii][3]);
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
                            if ((xp>=0)&&(xp<this.imageSize)&&(yp>=0)&&(yp<this.imageSize)){
                                res[i*model_number+u][xp][yp]=im[u*dp.param.sizeoutput*dp.param.sizeoutput+uu*dp.param.sizeoutput+uuu];
                            }
                        }
                    }
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
        
        
        
        
    
    
    
    public void free(){
        this.dp.free();
    }
    
    
    
}
