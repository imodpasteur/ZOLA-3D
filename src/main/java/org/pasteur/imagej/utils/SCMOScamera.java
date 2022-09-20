/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ImageProcessor;
import java.awt.Rectangle;

/**
 *
 * @author benoit
 */
public class SCMOScamera {
    
    
    public double [][] scmosvargain=null;
    public double [][] scmosvariance=null;
    public double [][] scmosgain=null;
    public double [][] scmosoffset=null;
    
    public String path_SCMOSvariance="";
    public String path_SCMOSoffset="";
    public String path_SCMOSgain="";
    
    public SCMOScamera(String path_SCMOSvariance,String path_SCMOSoffset,String path_SCMOSgain){
        this.path_SCMOSvariance=path_SCMOSvariance;
        this.path_SCMOSoffset=path_SCMOSoffset;
        this.path_SCMOSgain=path_SCMOSgain;
    }
    
    public boolean loadSCMOScameraFiles(int width, int height,Roi r){
        
        //SCMOS MAP loading
        ImagePlus impvar=null;
        ImagePlus impoff=null;
        ImagePlus impgain=null;
        try{
            impvar=new ImagePlus(path_SCMOSvariance);
            impoff=new ImagePlus(path_SCMOSoffset);
            impgain=new ImagePlus(path_SCMOSgain);
        }catch(Exception e){}

        double [][] scmos=null;
        if ((impvar!=null)&&(impoff!=null)&&(impgain!=null)){
            
            int w1=impvar.getWidth();
            int w2=impoff.getWidth();
            int w3=impgain.getWidth();
            int h1=impvar.getHeight();
            int h2=impoff.getHeight();
            int h3=impgain.getHeight();
            
            if ((w1!=w2)||(w1!=w3)||(w2!=w3)||(h1!=h2)||(h1!=h3)||(h2!=h3)){
                IJ.log("ERROR : SCMOS images (variance/gain/offset) should have the same size");
                IJ.log("widths "+w1+"  "+w2+"  "+w3);
                IJ.log("heights "+h1+"  "+h2+"  "+h3);
                return false;
            }
            if ((w1!=width)||(h1!=height)){
                IJ.log("ERROR : SCMOS images parameters (variance/gain/offset) and opened image should have the same size");
                return false;
            }
            int startX=0;
            int startY=0;
            if (r!=null){
                Rectangle bound=r.getBounds();
                width=bound.width;
                height=bound.height;
                startX=bound.x;
                startY=bound.y;
            }
                
                
                
                
            
            scmosvargain=new double [width][height];
            scmosvariance=new double [width][height];
            scmosgain=new double [width][height];
            scmosoffset=new double [width][height];
            
            ImageProcessor ipvar=impvar.getProcessor();
            ImageProcessor ipoff=impoff.getProcessor();
            ImageProcessor ipgain=impgain.getProcessor();
            for (int p=0;p<width;p++){
                for (int pp=0;pp<height;pp++){
                    scmosvariance[p][pp]=ipvar.getPixelValue(startX+p, startY+pp);
                    scmosoffset[p][pp]=ipoff.getPixelValue(startX+p, startY+pp);
                    scmosgain[p][pp]=ipgain.getPixelValue(startX+p, startY+pp);
                    scmosvargain[p][pp]=scmosvariance[p][pp]/(scmosgain[p][pp]*scmosgain[p][pp]);
                }
            }
            return true;
        }
        else{
            IJ.log("ERROR, problem loading SCMOS images (variance/gain/offset)");
            return false;
        }
                
    }
    
    
}
