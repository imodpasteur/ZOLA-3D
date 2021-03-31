/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import org.pasteur.imagej.process.InitBackgroundAndPhotonNumber;
import org.pasteur.imagej.process.gpu.DataPhase_;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.utils.ImageShow;
import org.pasteur.imagej.cuda.*;
import jcuda.runtime.JCuda;

import jcuda.runtime.cudaError;
import ij.IJ;
import java.awt.Color;
import ij.gui.Plot;
import jcuda.Pointer;

/**
 *
 * @author benoit
 */

public class ZernikePhaseRetrieval_ {
    //int stream=0;//not parallel process: iterative
    double [] deltaZ;
    double [][][][] image;
    //double [][][][] imageBlur;
    int methodLikelihood=0;
    int center;
    InitBackgroundAndPhotonNumber paramImage;
    DataPhase_ dparam;
    public Model3D_ [] model3D;
    
    double [][] x;
    double [][] y;
    double [][] z;
    int nbSlice;
    int nbstack;
    
    double imageLength;
    
    //double [] fit_b_0;
    //double [] fit_b_1;
    double [] fit_a;
    double [][] fit_a_each;
    double [][][] background;
    
    double epsilon=.00001;//stop criterion
    
    double [][] registrationStack;//z,x,y registration for each stack
    
    String path_calibration;
    
    int fitOrderBackground=2;
    PolynomialBackground [] pbg;
    double sigma;
    double [] zwat ;
    
    double dx=0;//drift X supposed linear
    double dy=0;//drift Y supposed linear
    double dz=0;//drift Z supposed linear
    
    double maxDrift=10;//max drift between each frame (nm)
    
    public ZernikePhaseRetrieval_(int sizeFFT,double xystep,double zstep,double wavelength,double noil,double na,int zernikeCoefNumber,InitBackgroundAndPhotonNumber paramImage,String path_calibration,double sigma,int axialside,boolean withApoFactor){
        this.sigma=sigma;
        this.image=paramImage.image;
        
        nbstack=image.length;
        nbSlice=image[0].length;
        imageLength=image[0][0].length*image[0][0][0].length;
        dparam = new DataPhase_(sizeFFT,image[0][0].length,0,xystep,zstep,wavelength,noil,na,1.0,zernikeCoefNumber,withApoFactor);
        dparam.setNwat(dparam.param.noil);//ca ne change rien normalement car bille collée à lamelle
        
        dparam.param.Zfocus=0;
        //astuce for registration
        dparam.phaseZer.setMatAtPosit(dparam.psf.getPointerkx(),0);
        dparam.phaseZer.setMatAtPosit(dparam.psf.getPointerky(),1);
        dparam.phaseZer.setMatAtPosit(dparam.psf.getPointerkzOil(),2);
            
        
        registrationStack=new double [3][nbstack];
        this.paramImage=paramImage;
        center=image[0].length/2;
        
        deltaZ=new double [nbSlice];
        //IJ.log("center is : "+center+"  !!!!!!!!!!§§§§§§§§§§§§§§§!!!!!!!");
        for (int s=0;s<nbSlice;s++){
            if (axialside==0){
                deltaZ[s]=(((double)s)-center)*dparam.param.zstep;
            }
            else{
                deltaZ[s]=-(((double)s)-center)*dparam.param.zstep;
            }
            
        }
        
        zwat = new double [nbSlice];
        for (int i=0;i<nbSlice;i++){
            zwat[i]=0;
        }
        
        x=new double[nbstack][nbSlice];//registration of each stack
        y=new double[nbstack][nbSlice];//registration of each stack
        z=new double[nbstack][nbSlice];//registration of each stack
        
        dparam.setMany(image[0].length);//images of 1 stack are computed in parallel
        
        fit_a=new double [nbstack];
        fit_a_each=new double [nbstack][nbSlice];
        //fit_b_0=new double [nbstack];
        //fit_b_1=new double [nbstack];
        
        for (int i=0;i<nbstack;i++){
            fit_a[i]=0;
            for (int ii=0;ii<nbSlice;ii++){
                fit_a_each[i][ii]=0;
            }
            //fit_b_0[i]=0;
            //fit_b_1[i]=0;
        }
        
        
        
        
        
        model3D=new Model3D_[image.length];
        
        
        background= new double[nbstack][image[0][0].length][image[0][0][0].length];
        //imageBlur=new double[nbstack][nbSlice][image[0][0].length][image[0][0][0].length];
        pbg=new PolynomialBackground[nbstack];
        for (int z=0;z<nbstack;z++){
            /*for (int ii=0;ii<nbSlice;ii++){
                for (int iii=0;iii<image[0][0].length;iii++){
                    for (int iiii=0;iiii<image[0][0][0].length;iiii++){
                        imageBlur[z][ii][iii][iiii] = 0;
                        double nb=0;
                        for (int a=-1;a<=1;a++){
                            for (int aa=-1;aa<=1;aa++){
                                if (((iii+a>=0)&&(iii+a<image[0][0].length))&&((iiii+aa>=0)&&(iiii+aa<image[0][0][0].length))){
                                    imageBlur[z][ii][iii][iiii] += image[z][ii][iii+a][iiii+aa];
                                    nb++;
                                }
                            }
                        }
                        imageBlur[z][ii][iii][iiii]/=imageBlur[z][ii][iii][iiii];
                    }
                }
            }*/
            model3D[z]=new Model3D_(dparam.param,this.image[z],paramImage.scmos[z]);
            pbg[z]= new PolynomialBackground();
            
            //init background fixed value
            int [] ind=pbg[z].getIndexOfOrder(0);
            pbg[z].a[ind[0]]=paramImage.B[z];
            
            computeBackground2D(z);
            
            
        }
        
        //dparam.setSizeoutput(image[0].length);
        
        for (int i=0;i<image.length;i++){
            
            double [] theA = new double [nbSlice];
            for (int ii=0;ii<theA.length;ii++){
                //theA[ii]=paramImage.A[i][ii]-imageLength*paramImage.B[i];
                theA[ii]=paramImage.A[i][ii];
            }
            model3D[i].setAmplitude(theA);
            /*double [] theB = new double [nbSlice];
            for (int ii=0;ii<theB.length;ii++){
                theB[ii]=paramImage.B[i];
            }
            model3D[i].setBackground(theB);*/
            
            
            
            
            
        }
        
        
        
        
        
        this.path_calibration=path_calibration;
    }
    
    
    
    
    
    public void run(int nbIter){
        
        //IJ.log("zstep "+dparam.param.zstep);
        
        //this.phase_retrieve_zernike_cross_validation(nbIter);
        
        this.phase_retrieve_zernike(nbIter);
        IJ.log("");
        for (int z=0;z<this.nbstack;z++){
            
            IJ.log("registration bead_"+z+"(x) = "+this.registrationStack[0][z]+"  (µm)");
            IJ.log("registration bead_"+z+"(y) = "+this.registrationStack[1][z]+"  (µm)");
            IJ.log("registration bead_"+z+"(z) = "+this.registrationStack[2][z]+"  (µm)");
            
        }
        IJ.log("");
        
        //at last stage, it updates tip, tilt and defocus according to registration of 1st bead
        dparam.phaseZer.setA(0, 0);//tip set to 0
        dparam.phaseZer.setA(1, 0);//tilt Y
        dparam.phaseZer.setA(2, 0);//tilt X
        
        if (path_calibration.length()>2){
            dparam.save(path_calibration);
            
            
            
            IJ.log("calibration file saved");
        }
        else{
            IJ.log("calibration file not saved...you should select a path to save it");
        }
        
        
        
    }
    
    
    
    public void free(){ 
        for (int i=0;i<image.length;i++){
            model3D[i].free();
        }
        dparam.free();
    }
    
    
    
    private void phase_retrieve_zernike_cross_validation(int iterations){
        
        IJ.log("phase retrieval cross validation");
        
        this.dparam.param.sigmaGaussianKernel=this.sigma;///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
        
        
        
        double likelihood=initialization(iterations);
        
        IJ.log("calibration started");
        
        
        int nbStackUse4PR=1;
        
        nbStackUse4PR=Math.min(nbStackUse4PR, nbstack);
                
        if (true){

            
            
            loop1:for (int t=0;t<iterations;t++){
                
                IJ.log("remaining iterations: "+(t)+"/"+iterations);
                //IJ.log("pass 2 ; remaining iterations: "+(iterations-t));
                IJ.showProgress((1.-(double)(iterations-t)/iterations));
                long tim=System.currentTimeMillis();

                

                this.updateRegistrationStacks(nbstack);
                
                
                this.updateSigmaGaussianKernel(nbStackUse4PR);///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                //IJ.log("sigma "+dparam.param.sigmaGaussianKernel);
                
                //this.updatePhotonB_0(nbstack);
                
                
                this.updatePhotonBpoly(nbstack,2);
                
                
                
                this.updatePhotonA(nbstack);
                
                
                
                
                this.updatePhotonAeach(nbstack);
                
                for (int k=0;k<paramImage.A[0].length;k++){

                    //IJ.log("A "+(paramImage.A[0][k]+this.fit_a[0]+this.fit_a_each[0][k]));
                }
                
                
                
                this.updatePhaseZernike(nbStackUse4PR);
                
                
                //this.updateDrift(nbstack);
                
                
                //IJ.log("drift "+dx+"  "+dy+"  "+dz);
                
                //this.updateWeightZ(nbstack);
                
                for (int ss=0;ss<nbstack;ss++){
                    //IJ.log("registration  x:"+this.registrationStack[0][ss]+"  y;"+this.registrationStack[1][ss]+"  z:"+this.registrationStack[2][ss]);
                }
                
                //IJ.log("lik : "+this.getLikelihood(nbstack));
                double lik=this.getLikelihood(nbstack);
                
                if (Math.abs(likelihood-lik)<epsilon){
                    break loop1;
                }
                
                likelihood=lik;
                
                

                //IJ.log("sig "+dparam.param.sigmaGaussianKernel+"   wz "+dparam.param.getweightZ());
                

            }
            
            
            
            
            
        }
        
        
        IJ.log("total number of image for training "+nbStackUse4PR);
        IJ.log("total number of image "+nbstack);
        
        showPhase();
        showImageAndModel();
        double r=0;
        for (int i=0;i<image.length;i++){
            double d=computeResidual(i);
            IJ.log("residual ("+i+") : "+d);
            r+=d;
        }
        IJ.log("global residual : "+r);
        
        
        IJ.showProgress(0);
        
    }
        
    
    private void phase_retrieve_zernike(int iterations){
        
        this.dparam.param.sigmaGaussianKernel=this.sigma;///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
        
        
        
        
        double likelihood=initialization(iterations);
        
        //showImageAndModel();
        
        IJ.log("calibration started");
        
        //second pass...all the stacks
        
        if (true){

            
            
            loop1:for (int t=0;t<iterations;t++){
                
                IJ.log("iteration: "+(t)+"/"+iterations);
                //IJ.log("pass 2 ; remaining iterations: "+(iterations-t));
                IJ.showProgress((1.-(double)(iterations-t)/iterations));
                

                

                this.updateRegistrationStacks(nbstack);
                
                this.updatePhaseZernike(nbstack);
                
                
                this.updateSigmaGaussianKernel(nbstack);///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                //IJ.log("sigma "+dparam.param.sigmaGaussianKernel);
                
                //this.updatePhotonB_0(nbstack);
                
                
                this.updatePhotonBpoly(nbstack,2);
                
                
                
                this.updatePhotonA(nbstack);
                
                
                
                
                this.updatePhotonAeach(nbstack);
                
                for (int k=0;k<paramImage.A[0].length;k++){

                    //IJ.log("A "+(paramImage.A[0][k]+this.fit_a[0]+this.fit_a_each[0][k]));
                }
                
                
                
                
                
                //this.updateDrift(nbstack);
                
                
                //IJ.log("drift "+dx+"  "+dy+"  "+dz);
                
                //this.updateWeightZ(nbstack);
                
                
                for (int ss=0;ss<nbstack;ss++){
                    //IJ.log("registration  x:"+this.registrationStack[0][ss]+"  y;"+this.registrationStack[1][ss]+"  z:"+this.registrationStack[2][ss]);
                }
                
                //IJ.log("lik : "+this.getLikelihood(nbstack));
                double lik=this.getLikelihood(nbstack);
                
                if (Math.abs(likelihood-lik)<epsilon){
                    break loop1;
                }
                
                likelihood=lik;
                
                

                //IJ.log("sig "+dparam.param.sigmaGaussianKernel+"   wz "+dparam.param.getweightZ());
                

            }
            
            
            
            
            
        }
        
        
        
        IJ.log("calibration finished");
        IJ.showProgress(0);
        
        showPhase();
        showImageAndModel();
        double r=0;
        for (int i=0;i<image.length;i++){
            double d=computeResidual(i);
            IJ.log("residual ("+i+") : "+d);
            r+=d;
        }
        IJ.log("global residual : "+r);
        
    }
    
    
    void showImageAndModel(){
        showImageAndModel(-1);
    }
    void showImageAndModel(int id){
        
        
        double [][][] modconcat = new double[image.length*image[0].length][image[0][0].length][image[0][0][0].length];
        double [][][] imageconcat = new double[image.length*image[0].length][image[0][0].length][image[0][0][0].length];
        double themin=Double.POSITIVE_INFINITY;
        double themax=Double.NEGATIVE_INFINITY;
        
        
            
        
        
        themin=Double.POSITIVE_INFINITY;
        themax=Double.NEGATIVE_INFINITY;
        double [][][][] mod = new double [model3D.length][][][];
        
        
         
            
            
        
        for (int i=0;i<model3D.length;i++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][i]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            compute2DModelMany_f(i);
            
            this.computeModel3D(i);
            
            mod[i]=model3D[i].getModel();
        }
        for (int u=0;u<image.length;u++){
            for (int i=0;i<image[u].length;i++){
                for (int ii=0;ii<image[u][i].length;ii++){
                    for (int iii=0;iii<image[u][i][ii].length;iii++){
                        if (themin>this.image[u][i][ii][iii]){
                            themin=this.image[u][i][ii][iii];
                        }
                        if (themax<this.image[u][i][ii][iii]){
                            themax=this.image[u][i][ii][iii];
                        }
                        imageconcat[u*image[0].length+i][ii][iii]=this.image[u][i][ii][iii];
                        modconcat[u*image[0].length+i][ii][iii]=mod[u][i][ii][iii];
                    }
                }
            }
        }
        
        
        
        //computeMeanAndVar(this.image,mod);
        ///////////////////////////////////////////////////////////////////////////
        if (id==-1){
            ImageShow.imshow(imageconcat,modconcat,"input-model image");
        }
        else{
            ImageShow.imshow(imageconcat,modconcat,"input-model image "+id);
        }
        IJ.setMinAndMax(themin, themax);
        
        
    }
    
    void showPhase(){
        
            
        float [][] ph=dparam.psf_fMany.getPhase();
        double themin=Double.POSITIVE_INFINITY;
        double themax=Double.NEGATIVE_INFINITY;
        for (int i=0;i<ph.length;i++){
            for (int ii=0;ii<ph[i].length;ii++){
                if (themin>ph[i][ii]){
                    themin=ph[i][ii];
                }
                if (themax<ph[i][ii]){
                    themax=ph[i][ii];
                }
            }
        }
        ImageShow.imshow(ph,"phase");
        IJ.setMinAndMax(themin, themax);
    }
    
    
    double computeResidual(int image_id){
        double [][][] mod=model3D[image_id].getModel();
        
        double som=0;
        int nbData=0;
        
            for (int ii=0;ii<image[image_id].length;ii++){
                for (int iii=0;iii<image[image_id][ii].length;iii++){
                    for (int iiii=0;iiii<image[image_id][ii][iii].length;iiii++){
                        som+=(image[image_id][ii][iii][iiii]-mod[ii][iii][iiii])*(image[image_id][ii][iii][iiii]-mod[ii][iii][iiii])/mod[ii][iii][iiii];
                        nbData++;
                    }
                }
            }
        som/=nbData;
        
        return som;
    }
    
    
    double initialization(int iterations){
        
        IJ.log("initialization started");
        
        
        double [][][][] stack;
        
        stack=this.image;
        
        
        
        
        int stack_size = stack[0].length;
        
        double stack_center=center;
        
        
        double epsilon=.0001;
        
        
        
        
        int maxorder=1;
        
        //first pass...ont stack
        double likelihood = Double.MAX_VALUE;
        
        
        IJ.log("initialization 1/6");
        
        likelihood=init(likelihood,1, 3,6, -3, 3,1);
        likelihood=init(likelihood,1, 6,10, -1, 1,1);
        IJ.log("initialization 2/6");
        likelihood=init(likelihood,1, 10,15, -1, 1,1);
        likelihood=init(likelihood,1, 15,21, -1, 1,1);
        
        //remove parabola in order to fit the model center using updateRegistrationStacks method
        dparam.phaseZer.setA(4, 0);
        
        IJ.log("initialization 3/6");
        
        loop1_1:for (int t=0;t<Math.min(3,iterations/2+1);t++){
            
            
            //IJ.log("pass 1 ; remaining iterations: "+(iterations-t));
            //IJ.showProgress(.05*(1.-(double)(iterations-t)/iterations));
            
            
            //this.updatePhaseZernike(1,6);
            
            this.updateRegistrationStacks(1);
            
            double lik=this.getLikelihood(1);
            
            if (Math.abs(likelihood-lik)<epsilon){
                break loop1_1;
            }
            
            likelihood=lik;
            
            //IJ.log("likelihood "+likelihood);
            
        }
        
        
        
        
        loop1_2:for (int t=0;t<Math.min(3,iterations/2+1);t++){
            
            
            //IJ.log("pass 1 ; remaining iterations: "+(iterations-t));
            //IJ.showProgress(.1*(1.-(double)(iterations-t)/iterations)+.05);
            
            
            this.updatePhaseZernike(1,10,false);
            this.updateRegistrationStacks(nbstack);
            this.updatePhotonBpoly(nbstack,0);
            
            
            
            double lik=this.getLikelihood(1);
            
            if (Math.abs(likelihood-lik)<epsilon){
                break loop1_2;
            }
            
            likelihood=lik;
            
            //IJ.log("likelihood "+likelihood);
            
        }
        
        
        
        IJ.log("initialization 4/6");
        
        //this.computeModel3D(0);
        //sresd=model3D[0].getModel();
        //ImageShow.imshow(sresd,"PSF_model_pass1");
        
        
        loop1_3:for (int t=0;t<Math.min(3,iterations/2+1);t++){
            
            
            //IJ.log("pass 1 ; remaining iterations: "+(iterations-t));
            //IJ.showProgress(.15*(1.-(double)(iterations-t)/iterations)+.15);
            
            
            this.updatePhaseZernike(1,15,false);
            this.updatePhotonBpoly(nbstack,1);
            this.updatePhotonA(nbstack);
            
            this.updateRegistrationStacks(nbstack);
            
            double lik=this.getLikelihood(1);
            if (Math.abs(likelihood-lik)<epsilon){
                break loop1_3;
            }
            
            likelihood=lik;
            
            //IJ.log("likelihood "+likelihood);
            
        }
        
        IJ.log("initialization 5/6");
        
        //this.computeModel3D(0);
        //sresd=model3D[0].getModel();
        //ImageShow.imshow(sresd,"PSF_model_pass2");
        
        for (int t=0;t<iterations/2;t++){
            this.updateRegistrationStacks(nbstack);
        }
        
        
        
        loop1_4:for (int t=0;t<Math.min(6,iterations/2+1);t++){
            
            
            //IJ.log("pass 1 ; remaining iterations: "+(iterations-t));
            
            
            this.updatePhaseZernike(1,15,false);
            
            this.updatePhotonBpoly(nbstack,2);
            
            
            this.updatePhotonAeach(nbstack);
            
            
            this.updateRegistrationStacks(nbstack);
            
            //this.updateSigmaGaussianKernel(nbstack);
            
            //pbg[0].log();
            //ImageShow.imshow(background,"bckg");

            
            double lik=this.getLikelihood(1);
            
            if (Math.abs(likelihood-lik)<epsilon){
                break loop1_4;
            }
            
            likelihood=lik;
            
            //IJ.log("likelihood "+likelihood);
            
        }
        
        
        
        
        IJ.log("initialization 6/6");
        
        this.updatePhaseZernike(1,28,false);
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        this.updateRegistrationStacks(nbstack);
        this.updatePhaseZernike(1,45,false);
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        this.updateRegistrationStacks(nbstack);
        this.updatePhaseZernike(1,66,false);
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        this.updateRegistrationStacks(nbstack);
        
        
        

        this.updateRegistrationStacks(nbstack);
        this.updatePhaseZernike(nbstack,false);
        this.updateRegistrationStacks(nbstack);
        this.updatePhaseZernike(nbstack,false);
        
        
        //this.computeModel3D(0);
        //sresd=model3D[0].getModel();
        //ImageShow.imshow(sresd,"PSF_model_pass3");
        
        likelihood = this.getLikelihood(nbstack);
        
        
        return likelihood;
    }
    
    
    
    
    
    
    
    
    double initializationNew(int iterations){
        
        
        
        
        double [][][][] stack;
        
        stack=this.image;
        
        
        
        
        int stack_size = stack[0].length;
        
        double stack_center=center;
        
        
        double epsilon=.0001;
        
        
        
        
        int maxorder=1;
        
        //first pass...ont stack
        double likelihood = Double.MAX_VALUE;
        
        
        
        for (int u=3;u<21;u++){dparam.phaseZer.setA(u, Math.random()*10);}
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
        for (int t=0;t<iterations;t++){
            this.updatePhotonBpoly(1,0);
            this.updateRegistrationStacks(1);
            this.updatePhotonA(1);
            double lik=this.getLikelihood(1);
            likelihood=lik;
        }
        this.updatePhaseZernike(1,21);
        double lik0=this.getLikelihood(1);
        showImageAndModel();
        IJ.log("lik0 "+lik0);
        
        
        for (int u=3;u<21;u++){dparam.phaseZer.setA(u, Math.random()*10);}
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
        for (int t=0;t<iterations;t++){
            this.updatePhotonBpoly(1,0);
            this.updateRegistrationStacks(1);
            this.updatePhotonA(1);
            double lik=this.getLikelihood(1);
            likelihood=lik;
        }
        this.updatePhaseZernike(1,21);
        double lik1=this.getLikelihood(1);
        showImageAndModel();
        IJ.log("lik1 "+lik1);
        
        
        for (int u=3;u<21;u++){dparam.phaseZer.setA(u, Math.random()*10);}
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
        for (int t=0;t<iterations;t++){
            this.updatePhotonBpoly(1,0);
            this.updateRegistrationStacks(1);
            this.updatePhotonA(1);
            double lik=this.getLikelihood(1);
            likelihood=lik;
        }
        this.updatePhaseZernike(1,21);
        double lik2=this.getLikelihood(1);
        showImageAndModel();
        IJ.log("lik2 "+lik2);
        
        
        for (int u=3;u<21;u++){dparam.phaseZer.setA(u, Math.random()*10);}
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
        for (int t=0;t<iterations;t++){
            this.updatePhotonBpoly(1,0);
            this.updateRegistrationStacks(1);
            this.updatePhotonA(1);
            double lik=this.getLikelihood(1);
            likelihood=lik;
        }
        this.updatePhaseZernike(1,21);
        double lik3=this.getLikelihood(1);
        showImageAndModel();
        IJ.log("lik3 "+lik3);
        
        
        
        
        
        
        
        
        
        likelihood = this.getLikelihood(nbstack);
        
        return likelihood;
    }
    
    
    
    private void computeMeanAndVar(double [][][][] im1,double [][][][] ref){
        //double [] mean = new double [im1.length*im1[0].length*im1[0][0].length];
        //double [] var = new double [im1.length*im1[0].length*im1[0][0].length];
        int k=0;
        for (int i=0;i<im1.length;i++){
            for (int ii=0;ii<im1[i].length;ii++){
                for (int iii=0;iii<im1[i][ii].length;iii++){
                    for (int iiii=0;iiii<im1[i][ii][iii].length;iiii++){
                        double mean=ref[i][ii][iii][iiii];
                        double var=(ref[i][ii][iii][iiii]-im1[i][ii][iii][iiii])*(ref[i][ii][iii][iiii]-im1[i][ii][iii][iiii]);
                        IJ.write(""+mean+","+var);
                    }
                }
            }
        }
        //this.plot(mean,var,"photon curve","mean","var");
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
    
    
    
    
    
    
    
    
    private double init(double likelihood,int stackNumber,int startcoefNumber,int endcoefNumber,double min,double max,double step){
        
        
        if (startcoefNumber<3){
            startcoefNumber=3;
        }
        
        
        if (startcoefNumber>=dparam.phaseZer.numCoef){
            return likelihood;
        }
        
        
        if (endcoefNumber>dparam.phaseZer.numCoef){
            endcoefNumber=dparam.phaseZer.numCoef;
        }
        
        
        double [] init=dparam.phaseZer.getA();
        
        
        
        int t=0;
        for (double val=min;val<=max;val+=step){
            t++;
        }
        double [] valstep = new double [t];
        t=0;
        for (double val=min;val<=max;val+=step){
            //IJ.log("val "+val+"  "+min+"  "+max+"  "+step);
            valstep[t]=val;
            t++;
        }
        int numberOfCoef=endcoefNumber-startcoefNumber;
        int numberOfVal=valstep.length;
        
        
        double [][] coco = new double [(int)Math.round(Math.pow(numberOfVal,numberOfCoef))][dparam.phaseZer.numCoef];
        
        int [] id = new int [numberOfCoef];
        for (int i=0;i<id.length;i++){
            id[i]=0;
        }
        for (int i=0;i<coco.length;i++){
            //on compte en base numberOfVal
            int number=i;
            for (int b=0;b<id.length;b++){
                id[b]=number%numberOfVal;
                number=number/numberOfVal;
            }
            for (int b=0;b<3;b++){
                coco[i][b]=0;
            }
            for (int b=3;b<startcoefNumber;b++){
                coco[i][b]=init[b];
            }
            for (int b=startcoefNumber,bb=0;b<endcoefNumber;b++,bb++){
                coco[i][b]=valstep[id[bb]]+init[b];
            }
            for (int b=endcoefNumber;b<dparam.phaseZer.numCoef;b++){
                coco[i][b]=init[b];
            }
        }
        
        int bestId=0;
        
        for (int p=0;p<coco.length;p++){
            
            
                
            double lik=0;
            dparam.phaseZer.setA(coco[p]);
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
        
            //showImageAndModel(p);
            for (int z=0;z<stackNumber;z++){

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model1");
                lik+=model3D[z].getLikelihood(methodLikelihood);
            }
            //IJ.log("lik("+p+") = "+lik);
            if (lik<=likelihood){
                likelihood=lik;
                bestId=p;
            }
            //IJ.log("lik "+p+"  "+val+"  "+lik);
            //if (p%100==0){
                //IJ.log("... "+(coco.length-p)/100);
            //}

        }
        
        dparam.phaseZer.setA(coco[bestId]);
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
        
        
        
        return likelihood;
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    private void updatePhaseZernike(int stackNumber){
        updatePhaseZernike(stackNumber,-1,true);
    }
    
    private void updatePhaseZernike(int stackNumber,int coefNumber){
        updatePhaseZernike(stackNumber,coefNumber,true);
    }
    
    
    
    private void updatePhaseZernike(int stackNumber,boolean withParabola){
        updatePhaseZernike(stackNumber,-1,withParabola);
    }
    
    
    private void updatePhaseZernike(int stackNumber,int coefNumber,boolean withParabola){
        
        if (coefNumber<3){
            coefNumber=dparam.phaseZer.numCoef;
        }
        
        if (coefNumber>dparam.phaseZer.numCoef){
            coefNumber=dparam.phaseZer.numCoef;
        }
        
        double h=.01;
        
        
        double lastLikelihood=-1;
        
       
        for (int p=3;p<coefNumber;p++){
        
            if ((withParabola)||(p!=4)){//do not update parabola if withParabola=false
                
                
                double lik1=0;
                for (int z=0;z<stackNumber;z++){
                    
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, -h));
                    
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");
                            
                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model1");
                    lik1+=model3D[z].getLikelihood(methodLikelihood);
                }
                
                

                
                
                
                
                double lik2=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");
                    
                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model2");
                    lik2+=model3D[z].getLikelihood(methodLikelihood);
                }
                
                
                
                double lik3=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, h));
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                    compute2DModelMany_f(z);
                    
                    this.computeModel3D(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                    //ImageShow.imshow(model3D[z].getModel(),"model3");
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                }
                
                
                
                double save=dparam.phaseZer.getA(p);
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                
                
                boolean found=false;
                loop:for (double gamma=100;gamma>.02;gamma/=10){
                    
                    double lik=0;
                    for (int z=0;z<stackNumber;z++){
                        for (int pp=0;pp<3;pp++){
                            dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                        }
                        dparam.phaseZer.setA(p, save-gamma*grad);
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                        compute2DModelMany_f(z);
                        this.computeModel3D(z);
                        lik+=model3D[z].getLikelihood(methodLikelihood);
                    }
                    
                    
                    if (lik<lik2){
                        lastLikelihood=lik;
                        found=true;
                        break loop;
                    }
                    else{
                        dparam.phaseZer.setA(p, save);
                    }
                }
                if (!found){
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    
                    
                }
            }
        }
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    private void updateDrift(int stackNumber){
        
        
        double h=.00001;
        
        
        double lastLikelihood=-1;
        
       
            
        {
            double lik1=0;
            for (int z=0;z<stackNumber;z++){

                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.x[z][u]=(this.dx-h)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model1");
                lik1+=model3D[z].getLikelihood(methodLikelihood);
            }







            double lik2=0;
            for (int z=0;z<stackNumber;z++){
                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.x[z][u]=(this.dx)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model2");
                lik2+=model3D[z].getLikelihood(methodLikelihood);
            }



            double lik3=0;
            for (int z=0;z<stackNumber;z++){
                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.x[z][u]=(this.dx+h)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                compute2DModelMany_f(z);

                this.computeModel3D(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                //ImageShow.imshow(model3D[z].getModel(),"model3");
                lik3+=model3D[z].getLikelihood(methodLikelihood);
            }

            double save=dx;
            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)/(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }

            //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad+"   "+dparam.param.stream);

            boolean found=false;
            loop:for (double gamma=100;gamma>.02;gamma/=10){
                this.dx-=grad*gamma;
                double lik=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }

                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                    for (int u=0;u<nbSlice;u++){
                        this.x[z][u]=(dx)*this.deltaZ[u]/dparam.param.zstep;
                    }

                    compute2DModelMany_f(z);
                    this.computeModel3D(z);
                    lik+=model3D[z].getLikelihood(methodLikelihood);
                }


                if (lik<lik2){
                    lastLikelihood=lik;
                    found=true;
                    break loop;
                }
                else{
                    this.dx=save;
                }
            }

            if (!found){
                
                for (int z=0;z<stackNumber;z++){
                    for (int u=0;u<nbSlice;u++){
                        x[z][u]=(dx)*this.deltaZ[u]/dparam.param.zstep;
                    }
                }

            }
            else{
                if (maxDrift<dx*1000.){
                    dx=maxDrift/1000.;
                    for (int z=0;z<stackNumber;z++){
                        for (int u=0;u<nbSlice;u++){
                            this.x[z][u]=(dx)*this.deltaZ[u]/dparam.param.zstep;
                        }
                    }
                }
            }
        }
        
        
        
        {
            double lik1=0;
            for (int z=0;z<stackNumber;z++){

                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.y[z][u]=(this.dy-h)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model1");
                lik1+=model3D[z].getLikelihood(methodLikelihood);
            }







            double lik2=0;
            for (int z=0;z<stackNumber;z++){
                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.y[z][u]=(this.dy)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model2");
                lik2+=model3D[z].getLikelihood(methodLikelihood);
            }



            double lik3=0;
            for (int z=0;z<stackNumber;z++){
                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.y[z][u]=(this.dy+h)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                compute2DModelMany_f(z);

                this.computeModel3D(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                //ImageShow.imshow(model3D[z].getModel(),"model3");
                lik3+=model3D[z].getLikelihood(methodLikelihood);
            }

            double save=dy;
            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)/(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }

            //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad+"   "+dparam.param.stream);

            boolean found=false;
            loop:for (double gamma=100;gamma>.02;gamma/=10){
                this.dy-=grad*gamma;
                double lik=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }

                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                    for (int u=0;u<nbSlice;u++){
                        this.y[z][u]=(dy)*this.deltaZ[u]/dparam.param.zstep;
                    }

                    compute2DModelMany_f(z);
                    this.computeModel3D(z);
                    lik+=model3D[z].getLikelihood(methodLikelihood);
                }


                if (lik<lik2){
                    lastLikelihood=lik;
                    found=true;
                    break loop;
                }
                else{
                    this.dy=save;
                }
            }

            if (!found){

                for (int z=0;z<stackNumber;z++){
                    for (int u=0;u<nbSlice;u++){
                        this.y[z][u]=(dy)*this.deltaZ[u]/dparam.param.zstep;
                    }
                }

            }
            else{
                if (maxDrift<dy*1000.){
                    dy=maxDrift/1000.;
                    for (int z=0;z<stackNumber;z++){
                        for (int u=0;u<nbSlice;u++){
                            this.y[z][u]=(dy)*this.deltaZ[u]/dparam.param.zstep;
                        }
                    }
                }
            }
        }
        
        
        
        if (false){
            double lik1=0;
            for (int z=0;z<stackNumber;z++){

                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.z[z][u]=(this.dz-h)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model1");
                lik1+=model3D[z].getLikelihood(methodLikelihood);
            }







            double lik2=0;
            for (int z=0;z<stackNumber;z++){
                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.z[z][u]=(this.dz)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                compute2DModelMany_f(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");

                this.computeModel3D(z);
                //ImageShow.imshow(model3D[z].getModel(),"model2");
                lik2+=model3D[z].getLikelihood(methodLikelihood);
            }



            double lik3=0;
            for (int z=0;z<stackNumber;z++){
                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                for (int u=0;u<nbSlice;u++){
                    this.z[z][u]=(this.dz+h)*this.deltaZ[u]/dparam.param.zstep;
                }

                //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                compute2DModelMany_f(z);

                this.computeModel3D(z);
                //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                //ImageShow.imshow(model3D[z].getModel(),"model3");
                lik3+=model3D[z].getLikelihood(methodLikelihood);
            }

            double save=dz;
            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)/(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }

            //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad+"   "+dparam.param.stream);

            boolean found=false;
            loop:for (double gamma=100;gamma>.02;gamma/=10){
                this.dz-=grad*gamma;
                double lik=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }

                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());

                    for (int u=0;u<nbSlice;u++){
                        this.z[z][u]=(dz)*this.deltaZ[u]/dparam.param.zstep;
                    }

                    compute2DModelMany_f(z);
                    this.computeModel3D(z);
                    lik+=model3D[z].getLikelihood(methodLikelihood);
                }


                if (lik<lik2){
                    lastLikelihood=lik;
                    found=true;
                    break loop;
                }
                else{
                    this.dz=save;
                }
            }

            if (!found){

                for (int z=0;z<stackNumber;z++){
                    for (int u=0;u<nbSlice;u++){
                        this.z[z][u]=(dz)*this.deltaZ[u]/dparam.param.zstep;
                    }
                }

            }
            else{
                if (maxDrift<dz*1000.){
                    dz=maxDrift/1000.;
                    for (int z=0;z<stackNumber;z++){
                        for (int u=0;u<nbSlice;u++){
                            this.z[z][u]=(dz)*this.deltaZ[u]/dparam.param.zstep;
                        }
                    }
                }
            }
        }
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    private void updateRegistrationStacks(int stackNumber){
        
        if (stackNumber>=1){
            double h=.01;


            


            for (int p=0;p<3;p++){//first 3 coefs: used for X, Y, Z shifts
                
                
                for (int z=0;z<stackNumber;z++){
                    
                    //dparam.phaseZer.setA(p, this.registrationStack[p][z]);
                    for (int pp=0;pp<3;pp++){
                            dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    
                    double lik1=0;
                    
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, -h));

                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");

                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model1");
                    lik1+=model3D[z].getLikelihood(methodLikelihood);
                    


                    




                    double lik2=0;
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");

                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model2");
                    lik2+=model3D[z].getLikelihood(methodLikelihood);
                    
                    


                    double lik3=0;
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination(p, h));
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                    compute2DModelMany_f(z);

                    this.computeModel3D(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                    //ImageShow.imshow(model3D[z].getModel(),"model3");
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                    

                    //IJ.log("save "+dparam.phaseZer.getA(p)+"  "+this.registrationStack[p][z]);
                    double save=this.registrationStack[p][z];
                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad+"   "+dparam.param.stream);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.02;gamma/=10){

                        double lik=0;
                        dparam.phaseZer.setA(p, save-gamma*grad);
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                        compute2DModelMany_f(z);
                        this.computeModel3D(z);
                        lik+=model3D[z].getLikelihood(methodLikelihood);

                        
                        if (lik<lik2){
                            
                            this.registrationStack[p][z]=dparam.phaseZer.getA(p);
                            found=true;
                            break loop;
                        }
                        else{
                            dparam.phaseZer.setA(p, save);
                        }
                    }
                    if (!found){
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());


                    }
                }
            }
            
        }
        
    }
    
    
    
    
    
    
    
    private void updateSigmaGaussianKernel(int stackNumber){
        
        
        
        double h=.001;
        
        
        double lastLikelihood=-1;
        
       
            
        
        double lik1=0;
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel-h);
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik1+=model3D[z].getLikelihood(methodLikelihood);
        }
        
        
        double lik2=0;
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik2+=model3D[z].getLikelihood(methodLikelihood);
        }
        

        double lik3=0;
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel+h);
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik3+=model3D[z].getLikelihood(methodLikelihood);
        }
        
        

        double save=dparam.param.sigmaGaussianKernel;
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)/(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }



        boolean found=false;
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            
            dparam.param.sigmaGaussianKernel=save-gamma*grad;
            dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
            double lik=0;
            for (int z=0;z<stackNumber;z++){
                
                for (int pp=0;pp<3;pp++){//update due to registration of stacks
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
                compute2DModelMany_f(z);
                this.computeModel3D(z);
                lik+=model3D[z].getLikelihood(methodLikelihood);
            }
            
            
            if (lik<lik2){
                lastLikelihood=lik;
                found=true;
                break loop;
            }
            else{
                dparam.param.sigmaGaussianKernel=save;
            }
        }
        if (!found){
            dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
            
        }
            
        
    }
    
    
    
    
    
    
    
    
    private void updateWeightZ(int stackNumber){
        
        
        
        double h=.001;
        
        
        double lastLikelihood=-1;
        
       
            
        
        double lik1=0;
        double wSave=dparam.param.getweightZ();
       
        dparam.param.updateweightZ(wSave-h);
        
        
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik1+=model3D[z].getLikelihood(methodLikelihood);
        }
        
        
        double lik2=0;
        dparam.param.updateweightZ(wSave);
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik2+=model3D[z].getLikelihood(methodLikelihood);
        }
        

        double lik3=0;
        dparam.param.updateweightZ(wSave+h);
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik3+=model3D[z].getLikelihood(methodLikelihood);
        }
        
        

        double save=wSave;
        double grad=0;
        if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
            grad=((lik3-lik1)/(2*h));
        }
        else{
            grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
        }



        boolean found=false;
        loop:for (double gamma=1;gamma>.02;gamma/=10){
            dparam.param.updateweightZ(save-gamma*grad);
            double lik=0;
            for (int z=0;z<stackNumber;z++){
                
                for (int pp=0;pp<3;pp++){//update due to registration of stacks
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
                compute2DModelMany_f(z);
                this.computeModel3D(z);
                lik+=model3D[z].getLikelihood(methodLikelihood);
            }
            
            
            if (lik<lik2){
                lastLikelihood=lik;
                found=true;
                break loop;
            }
            else{
                dparam.param.updateweightZ(save);
            }
        }
        if (!found){
            dparam.param.updateweightZ(save);
            
        }
            
        
    }
    
    
    
    
    
    
    
    
    private void updatePhotonA(int stackNumber){
        
            double h=1;

            
            
            
                for (int z=0;z<stackNumber;z++){
                    
                    double [] theA = new double [paramImage.A[z].length];
                    
                    
                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    compute2DModelMany_f(z);
                    
                    
                    
                    
                    double lik1=0;
                    
                    
                    for (int i=0;i<theA.length;i++){
                        //theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+fit_b_1[z]*(double)i)-h;
                        theA[i]=paramImage.A[z][i]+fit_a[z]+fit_a_each[z][i]-h;
                        
                    }
                    
                    model3D[z].setAmplitude(theA);
                    
                    this.computeModel3D(z);
                    
                    lik1+=model3D[z].getLikelihood(methodLikelihood);
                    







                    double lik2=0;
                    
                    for (int i=0;i<theA.length;i++){
                        //theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+fit_b_1[z]*(double)i);
                        theA[i]=paramImage.A[z][i]+fit_a[z]+fit_a_each[z][i];
                    }
                    
                    model3D[z].setAmplitude(theA);
                    
                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model2");
                    lik2+=model3D[z].getLikelihood(methodLikelihood);
                    
                    


                    double lik3=0;
                    
                    
                    for (int i=0;i<theA.length;i++){
                        //theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+fit_b_1[z]*(double)i)+h;
                        theA[i]=paramImage.A[z][i]+fit_a[z]+fit_a_each[z][i]+h;
                    }
                    
                    model3D[z].setAmplitude(theA);
                    
                    this.computeModel3D(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                    //ImageShow.imshow(model3D[z].getModel(),"model3");
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                    


                    
                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.002;gamma/=10){

                        double lik=0;
                        for (int i=0;i<theA.length;i++){
                            //theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+fit_b_1[z]*(double)i)-gamma*grad;
                            theA[i]=paramImage.A[z][i]+fit_a[z]+fit_a_each[z][i]-gamma*grad;
                        }

                        model3D[z].setAmplitude(theA);
                        this.computeModel3D(z);
                        lik+=model3D[z].getLikelihood(methodLikelihood);


                        if (lik<lik2){
                            
                            fit_a[z]=fit_a[z]-gamma*grad;
                            
                            found=true;
                            break loop;
                        }
                        else{
                            //nothing to replace
                        }
                    }
                    if (!found){
                        for (int i=0;i<theA.length;i++){
                            //theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+fit_b_1[z]*(double)i);
                            theA[i]=paramImage.A[z][i]+fit_a[z]+fit_a_each[z][i];
                        }
                        model3D[z].setAmplitude(theA);

                    }
                }
            
            
        
        
    }
    
    
    
    
    
    private void updatePhotonAeach(int stackNumber){
        
        double h=1;

        for (int z=0;z<stackNumber;z++){
            
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            compute2DModelMany_f(z);
            
            double [] theA = new double [paramImage.A[z].length];
            double [] saveA = new double [paramImage.A[z].length];
            
            for (int i=0;i<theA.length;i++){
                //theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+fit_b_1[z]*(double)i)-h;
                theA[i]=paramImage.A[z][i]+fit_a[z]+fit_a_each[z][i];
                saveA[i]=theA[i];
            }
                   
            double [] lik1;
            double [] lik2;
            double [] lik3;
            double [] lik;
            
            for (int a=0;a<theA.length;a++){
                
                theA[a]=saveA[a]-h;
            }
                
            model3D[z].setAmplitude(theA);

            this.computeModel3D(z);

            lik1=model3D[z].getLikelihoodStack(methodLikelihood);
            //double llik1=model3D[z].getLikelihood(methodLikelihood);
            
            
            for (int a=0;a<theA.length;a++){
                
                theA[a]=saveA[a];
            }

            model3D[z].setAmplitude(theA);

            this.computeModel3D(z);
            
            lik2=model3D[z].getLikelihoodStack(methodLikelihood);
            //double llik2=model3D[z].getLikelihood(methodLikelihood);
            
            

            for (int a=0;a<theA.length;a++){
                theA[a]=saveA[a]+h;
            }
            
            model3D[z].setAmplitude(theA);

            this.computeModel3D(z);

            lik3=model3D[z].getLikelihoodStack(methodLikelihood);
            //double llik3=model3D[z].getLikelihood(methodLikelihood);

            


            double [] grad = new double [theA.length];
            for (int a=0;a<theA.length;a++){
                if (Math.abs((lik3[a]+lik1[a]-2*lik2[a])/(h*h))==0){
                    grad[a]=((lik3[a]-lik1[a])/(2*h));
                }
                else{
                    grad[a]=((lik3[a]-lik1[a])/(2*h))/Math.abs((lik3[a]+lik1[a]-2*lik2[a])/(h*h));
                }
                //IJ.log("lik(a)  "+a+"  "+lik1[a]+"  "+lik2[a]+"  "+lik3[a]+"  "+grad[a]);
            }

                
            
            boolean [] found= new boolean [theA.length];
            for (int a=0;a<theA.length;a++){
                found[a]=false;
            }
            loop:for (double gamma=1;gamma>.002;gamma/=10){
                
                for (int a=0;a<theA.length;a++){
                    
                    if (!found[a]){

                        theA[a]=saveA[a]-gamma*grad[a];
                    }
                }

                model3D[z].setAmplitude(theA);
                this.computeModel3D(z);
                lik=model3D[z].getLikelihoodStack(methodLikelihood);
                
                for (int a=0;a<theA.length;a++){
                    if (!found[a]){
                        if (lik[a]<lik2[a]){

                            fit_a_each[z][a]=fit_a_each[z][a]-gamma*grad[a];

                            found[a]=true;
                        }
                    }
                }
                int numberNotFount=0;
                for (int a=0;a<theA.length;a++){
                    if (!found[a]){
                        theA[a]=saveA[a];
                        numberNotFount++;
                    }
                }
                model3D[z].setAmplitude(theA);
                //IJ.log("numberNotFount "+gamma+"  "+numberNotFount);
                if (numberNotFount==0){
                    break loop;

                }
            }
        }
            
            
        
        
    }
    
    
    
    
    
    
    
//    private void updatePhotonB_0(int stackNumber){
//        
//            double h=0.001;
//
//            
//            
//            
//                for (int z=0;z<stackNumber;z++){
//                    
//                    double [] theA = new double [this.nbSlice];
//                    double [] theB = new double [this.nbSlice];
//                    
//                    
//                    
//                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
//                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
//                    }
//                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
//                    compute2DModelMany_f(z);
//                    
//                    
//                    
//                    
//                    double lik1=0;
//                    
//                    
//                    for (int i=0;i<theB.length;i++){
//                        theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-h;
//                    }
//                    model3D[z].setBackground(theB);
//                    
//                    for (int i=0;i<theA.length;i++){
//                        theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-h);
//                    }
//                    model3D[z].setAmplitude(theA);
//                    
//                    this.computeModel3D(z);
//                    
//                    lik1+=model3D[z].getLikelihood(methodLikelihood);
//                    
//
//
//
//
//
//
//
//                    double lik2=0;
//                    
//                    for (int i=0;i<theB.length;i++){
//                        theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i);
//                    }
//                    model3D[z].setBackground(theB);
//                    
//                    for (int i=0;i<theA.length;i++){
//                        theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i));
//                    }
//                    model3D[z].setAmplitude(theA);
//                    
//                    this.computeModel3D(z);
//                    //ImageShow.imshow(model3D[z].getModel(),"model2");
//                    lik2+=model3D[z].getLikelihood(methodLikelihood);
//                    
//                    
//
//
//                    double lik3=0;
//                    
//                    
//                    for (int i=0;i<theB.length;i++){
//                        theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)+h;
//                    }
//                    model3D[z].setBackground(theB);
//                    
//                    for (int i=0;i<theA.length;i++){
//                        theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)+h);
//                    }
//                    model3D[z].setAmplitude(theA);
//                    
//                    this.computeModel3D(z);
//                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
//                    //ImageShow.imshow(model3D[z].getModel(),"model3");
//                    lik3+=model3D[z].getLikelihood(methodLikelihood);
//                    
//
//
//                    
//                    double grad=0;
//                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
//                        grad=((lik3-lik1)/(2*h));
//                    }
//                    else{
//                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
//                    }
//
//                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);
//
//                    boolean found=false;
//                    loop:for (double gamma=1;gamma>.02;gamma/=10){
//
//                        double lik=0;
//                        
//                        for (int i=0;i<theB.length;i++){
//                            theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-gamma*grad;
//                        }
//                        model3D[z].setBackground(theB);
//                        
//                        for (int i=0;i<theA.length;i++){
//                            theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-gamma*grad);
//                        }
//                        model3D[z].setAmplitude(theA);
//                    
//                        this.computeModel3D(z);
//                        lik+=model3D[z].getLikelihood(methodLikelihood);
//
//
//                        if (lik<lik2){
//                            
//                            fit_b_0[z]=fit_b_0[z]-gamma*grad;
//                            
//                            found=true;
//                            break loop;
//                        }
//                        else{
//                            //nothing to replace
//                        }
//                    }
//                    if (!found){
//                        
//                        for (int i=0;i<theB.length;i++){
//                            theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i);
//                        }
//                        model3D[z].setBackground(theB);
//                        
//                        for (int i=0;i<theA.length;i++){
//                            theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i));
//                        }
//                        model3D[z].setAmplitude(theA);
//                    }
//                }
//            
//            
//        
//        
//    }
//    
//    
//    
//    
//    
//    
//    
//    private void updatePhotonB_1(int stackNumber){
//        
//            double h=0.001;
//
//            
//            
//            
//                for (int z=0;z<stackNumber;z++){
//                    
//                    double [] theA = new double [this.nbSlice];
//                    double [] theB = new double [this.nbSlice];
//                    
//                    
//                    
//                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
//                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
//                    }
//                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
//                    compute2DModelMany_f(z);
//                    
//                    
//                    
//                    
//                    double lik1=0;
//                    
//                    
//                    for (int i=0;i<theB.length;i++){
//                        theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-(h*(double)i);
//                    }
//                    model3D[z].setBackground(theB);
//                    
//                    for (int i=0;i<theA.length;i++){
//                        theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-(h*(double)i));
//                    }
//                    model3D[z].setAmplitude(theA);
//                    
//                    this.computeModel3D(z);
//                    
//                    lik1+=model3D[z].getLikelihood(methodLikelihood);
//                    
//
//
//
//
//
//
//
//                    double lik2=0;
//                    
//                    for (int i=0;i<theB.length;i++){
//                        theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i);
//                    }
//                    model3D[z].setBackground(theB);
//                    
//                    for (int i=0;i<theA.length;i++){
//                        theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i));
//                    }
//                    model3D[z].setAmplitude(theA);
//                    
//                    this.computeModel3D(z);
//                    //ImageShow.imshow(model3D[z].getModel(),"model2");
//                    lik2+=model3D[z].getLikelihood(methodLikelihood);
//                    
//                    
//
//
//                    double lik3=0;
//                    
//                    
//                    for (int i=0;i<theB.length;i++){
//                        theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)+(h*(double)i);
//                    }
//                    model3D[z].setBackground(theB);
//                    
//                    for (int i=0;i<theA.length;i++){
//                        theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)+(h*(double)i));
//                    }
//                    model3D[z].setAmplitude(theA);
//                    
//                    this.computeModel3D(z);
//                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
//                    //ImageShow.imshow(model3D[z].getModel(),"model3");
//                    lik3+=model3D[z].getLikelihood(methodLikelihood);
//                    
//
//
//                    
//                    double grad=0;
//                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
//                        grad=((lik3-lik1)/(2*h));
//                    }
//                    else{
//                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
//                    }
//
//                    IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);
//
//                    boolean found=false;
//                    loop:for (double gamma=1;gamma>.02;gamma/=10){
//
//                        double lik=0;
//                        
//                        for (int i=0;i<theB.length;i++){
//                            theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-(gamma*grad*(double)i);
//                        }
//                        model3D[z].setBackground(theB);
//                        
//                        for (int i=0;i<theA.length;i++){
//                            theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i)-(gamma*grad*(double)i));
//                        }
//                        model3D[z].setAmplitude(theA);
//                    
//                        this.computeModel3D(z);
//                        lik+=model3D[z].getLikelihood(methodLikelihood);
//
//
//                        if (lik<lik2){
//                            
//                            fit_b_1[z]=fit_b_1[z]-(gamma*grad);
//                            
//                            found=true;
//                            break loop;
//                        }
//                        else{
//                            //nothing to replace
//                        }
//                    }
//                    if (!found){
//                        IJ.log("not found");
//                        for (int i=0;i<theB.length;i++){
//                            theB[i]=paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i);
//                        }
//                        model3D[z].setBackground(theB);
//                        
//                        for (int i=0;i<theA.length;i++){
//                            theA[i]=paramImage.A[z][i]+fit_a[z]-imageLength*(paramImage.B[z]+fit_b_0[z]+(fit_b_1[z]*(double)i));
//                        }
//                        model3D[z].setAmplitude(theA);
//                    }
//                }
//            
//            
//        
//        
//    }
    
    
    
    
    
    
    /*
    
    private void updatePhotonB(int stackNumber){
        
            double h=.001;

            
            
            
                for (int z=0;z<stackNumber;z++){
                    
                    double theB;
                    
                    
                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    compute2DModelMany_f(z);
                    
                    
                    
                    
                    double lik1=0;
                    
                    
                    theB = paramImage.B[z]+fit_b_0[z]-h;
                    
                    
                    
                    this.addPhotons(z, theB);
                    
                    lik1+=model3D[z].getLikelihood(methodLikelihood);
                    







                    double lik2=0;
                    
                    theB = paramImage.B[z]+fit_b_0[z];
                    
                    
                    this.addPhotons(z, theB);
                    
                    
                    lik2+=model3D[z].getLikelihood(methodLikelihood);
                    
                    


                    double lik3=0;
                    
                    
                    theB = paramImage.B[z]+fit_b_0[z]+h;
                    
                    
                    this.addPhotons(z, theB);
                    
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                    


                    
                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.02;gamma/=10){

                        double lik=0;
                        

                        theB = paramImage.B[z]+fit_b_0[z]-gamma*grad;
                    
                        this.addPhotons(z, theB);
                        
                        lik+=model3D[z].getLikelihood(methodLikelihood);


                        if (lik<lik2){
                            
                            fit_b_0[z]=fit_b_0[z]-gamma*grad;
                            
                            found=true;
                            break loop;
                        }
                        else{
                            //nothing to replace
                        }
                    }
                    if (!found){
                        
                        //nothing to replace 

                    }
                }
            
            
        
        
    }
    
    
    */
    
    
    //order c
    private void updatePhotonBpoly(int stackNumber,int order){
        
        
        
        double h=.0001;



        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            compute2DModelMany_f(z);
            
            for (int o=0;o<=order;o++){
                
                int [] index=this.pbg[z].getIndexOfOrder(o);
                for (int ind=0;ind<index.length;ind++){
                    
                    
                    
                    
                    double saveB=pbg[z].a[index[ind]];

                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
                    compute2DModelMany_f(z);




                    double lik1=0;


                    pbg[z].a[index[ind]] = saveB-h;
                    this.computeBackground2D(z);
                    this.computeModel3D(z);

                    lik1+=model3D[z].getLikelihood(methodLikelihood);








                    double lik2=0;

                    pbg[z].a[index[ind]] = saveB;
                    

                    this.computeBackground2D(z);
                    this.computeModel3D(z);

                    lik2+=model3D[z].getLikelihood(methodLikelihood);




                    double lik3=0;


                    pbg[z].a[index[ind]] = saveB+h;


                    this.computeBackground2D(z);
                    this.computeModel3D(z);
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                    



                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.02;gamma/=10){

                        double lik=0;


                        pbg[z].a[index[ind]] = saveB-gamma*grad;

                        this.computeBackground2D(z);
                        
                        this.computeModel3D(z);
                        lik+=model3D[z].getLikelihood(methodLikelihood);


                        if (lik<lik2){

                            //IJ.log("found "+lik+"  "+lik2+"  "+pbg[z].a[index[ind]]);

                            found=true;
                            break loop;
                        }
                        else{
                            
                            //nothing to replace
                        }
                    }
                    if (!found){
                        pbg[z].a[index[ind]]= saveB;
                        this.computeBackground2D(z);
                        
                        //nothing to replace 

                    }
                }
            }
        }


        
        
    }
    
    
    
    
    
    private double  getLikelihood(int stackNumber){
        
        
        
        
        double likelihood=0;
        
       
            
        
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombination());
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            likelihood+=model3D[z].getLikelihood(methodLikelihood);
        }
        
        
            
        return likelihood;
        
    }
    
    
    void compute2DModelMany_f(int id){
        
        //dparam.psf.computePSF(dparam.param.centerX,dparam.param.centerY,z);.
        double [] zz = new double [nbSlice];
        for (int i=0;i<nbSlice;i++){
            zz[i]=deltaZ[i]+z[id][i];
        }
        
        dparam.psf_fMany.computePSF(x[id], y[id], zz,zwat);
        //dparam.psf.computePSF(0,0,z);
        
        model3D[id].setPSFMany(dparam.psf_fMany.getPointerPSF());
        
    }
    
    
    
    
    
    
    
    
    void computeModel3D(int id){
        //model3D[id].computeModel();
        model3D[id].computeModel_background2D();
    }
    
    
    void computeBackground2D(int idZ){
        
        for (int i=0;i<background[idZ].length;i++){
            for (int ii=0;ii<background[idZ][0].length;ii++){

                    double [] v = new double [2];
                    v[0]=i;
                    v[1]=ii;
                    background[idZ][i][ii]=pbg[idZ].transform(v);

                    
            }
        }
        
        model3D[idZ].setBackground2D(background[idZ]);
        
    }
    
    class PolynomialBackground{
        
        PolynomialFit pf;
        public double [] a;//parameters
        int order;
        int [][] indices;//[size of a][dimension] : index ef exposant
        int dim=2;
        int base;
        int size;
        PolynomialBackground(){
            this.order=2;
            base=this.order+1;
            
            constructIndices();
            size=this.indices.length;
            a=new double [size];
        }
        
        
        
        int [] getIndexOfOrder(int ord){
            
            int number=0;
            for (int i=0;i<indices.length;i++){
                int sum=0;
                for (int j=0;j<indices[i].length;j++){
                    sum+=indices[i][j];
                }
                if (sum==ord){
                    number++;
                }
            }
            
            int [] r = new int[number];
            for (int i=0,k=0;i<indices.length;i++){
                int sum=0;
                for (int j=0;j<indices[i].length;j++){
                    sum+=indices[i][j];
                }
                if (sum==ord){
                    r[k++]=i;
                }
            }
            return r;
            
        }
    
        //return position from i, j, k, ...
        int getIndice(int [] i){
            int posit=0;
            for (int t=0;t<i.length;t++){
                posit+=i[t]*((int)Math.pow(base,t));
            }
            return posit;
        }
        
        
    
        public double transform(double [] X){//transform one position of N dimensions
            if (X.length!=dim){
                System.out.println("problem dimension registration transform");
            }
            
                double res=0;
                for (int i=0;i<a.length;i++){
                    double prod=a[i];
                    //IJ.log("a "+a[d][i]);
                    int [] indiceLigne = this.indices[i];
                    for (int ii=0;ii<X.length;ii++){
    //                    if (indiceLigne[ii]==0){
    //                        
    //                    }
    //                    else if (indiceLigne[ii]==1){
    //                        prod*=X[ii];
    //                    }
    //                    else{
                            prod*=Math.pow(X[ii], indiceLigne[ii]);
    //                    }
                        //IJ.log("* "+X[ii]+" ^ "+indiceLigne[ii]);
                    }
                    res+=prod;
                    //IJ.log("res "+res[d]+"   "+prod);
                }
            
            return res;
        }
    



        //inverse of getIndice : return i j k... from position
        int [] getIndices(int val){
            int [] i = new int[dim];
            for (int t=dim-1;t>=0;t--){
                i[t]=val/((int)Math.pow(base,t));
                val=val%((int)Math.pow(base,t));
            }
            return i;
        }



        void constructIndices(){
            int nbMajore=(int)Math.pow(base,dim);
            boolean [] ok = new boolean[nbMajore];
            int count=0;
            for (int u=0;u<nbMajore;u++){
                int [] p=this.getIndices(u);
                int valid=0;
                for (int uu=0;uu<p.length;uu++){
                    valid+=p[uu];
                }
                if (valid<base){
                    count++;
                    ok[u]=true;
                }
                else{
                    ok[u]=false;
                }
            }
            this.indices=new int[count][dim];
            int j=0;
            for (int u=0;u<nbMajore;u++){
                if (ok[u]){
                    int [] p=this.getIndices(u);
                    this.indices[j++]=p;
                }
            }
        }


        public void log(){
            String [] vect = new String [7];
            vect[0]="x^";
            vect[1]="y^";
            vect[2]="z^";
            vect[3]="t^";
            vect[4]="u^";
            vect[5]="v^";
            vect[6]="w^";

            if (dim<7){
                
                    String s="";
                    for (int u=0;u<size;u++){
                        s+=a[u];
                        for (int uu=0;uu<dim;uu++){
                            if (indices[u][uu]!=0){
                                if (indices[u][uu]!=1){
                                    s+="*("+vect[uu]+indices[u][uu]+")";
                                }
                                else{
                                    s+="*"+vect[uu].substring(0, 1);
                                }
                            }
                        }
                        if ((u<size-1))
                            s+=" + ";
                    }
                    IJ.log("p ="+s);

                
            }
            else{
                IJ.log("print yet not implemented with dimension > 7");
            }


        }



    
    
    }
    
    
}
 