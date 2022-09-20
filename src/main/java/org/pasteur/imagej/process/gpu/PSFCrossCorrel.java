/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import ij.IJ;
import org.pasteur.imagej.utils.FileVectorLoader;
import org.pasteur.imagej.utils.ImageShow;

/**
 *
 * @author benoit
 */
public class PSFCrossCorrel {
    
    
    
    
    double axialRange; 
    double stepZ;
    
    DataPhase_ dp;
    
    
    
    
    public PSFCrossCorrel(DataPhase_ dp,double axialRange, double stepZ){
        this.dp=dp;
        this.axialRange=axialRange;
        this.stepZ=stepZ;
        
        
    }
    
    
    
    
    
    
    public void run(){
        run(null);
    }
    
    
    public void run(String path){
        
        dp.psf.resetKz();
        
        SearchPSFcenter_ spsfc= new SearchPSFcenter_(dp,axialRange);
        double position=spsfc.getPosition();
        
        double minZ=position-axialRange/2;
        double maxZ=position+axialRange/2;
        
        
        Cross predGPU=new Cross(dp.param.sizeoutput,dp,minZ, maxZ, stepZ,.3);
        double [][][] psf=predGPU.getPSFNonNormalized3D();
        double []range=predGPU.getRange();
        
        
        ImageShow.imshow(psf);
        
        
        
        for (int i=0;i<psf.length;i++){
            predGPU.setImage(psf[i]);
            float [][][] corr=new float[psf.length][psf[0].length][psf[0][0].length];
            predGPU.convolveNormalizedInFourierDomain(corr);
            ImageShow.imshow(corr,"corr "+i+"  "+range[i]);
        }
        
        
        
        
        
        
        
        
        
        
        predGPU.free();
        
        
        
    }
    
    
    
    
    
    
}
