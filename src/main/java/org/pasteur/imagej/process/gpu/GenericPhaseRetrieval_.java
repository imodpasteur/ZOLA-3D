/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;

import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.process.*;
import jcuda.runtime.JCuda;

import jcuda.runtime.cudaError;
import ij.IJ;
import java.awt.Color;
import ij.gui.Plot;
import jcuda.Pointer;
/*
 *
 * @author benoit
 */
public class GenericPhaseRetrieval_ {
    
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
    
    double lambda=0;//regularization
    
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
    
    int [][] phasePosit;// from x,y -> return position in phase vector
    
    
    GS_ gs;
    
    public GenericPhaseRetrieval_(int sizeFFT,double xystep,double zstep,double wavelength,double noil,double na,InitBackgroundAndPhotonNumber paramImage,String path_calibration,double sigma,int axialside,boolean withApoFactor){
        this.sigma=sigma;
        this.image=paramImage.image;
        
        nbstack=image.length;
        nbSlice=image[0].length;
        imageLength=image[0][0].length*image[0][0][0].length;
        //dparam = new classes.cudaProcess.DataPhase(sizeFFT,image[0][0].length,0,xystep,zstep,wavelength,noil,na,1.0,zernikeCoefNumber);
        dparam = new DataPhase_(sizeFFT,image[0][0].length,0,xystep,zstep,wavelength,noil,na,1.0,withApoFactor);
        dparam.setNwat(dparam.param.noil);//ca ne change rien normalement car bille collée à lamelle
        dparam.param.zernikedPSF=false;
        dparam.param.Zfocus=0;
        //astuce for registration
        dparam.phaseZer.setMatAtPosit(dparam.psf.getPointerkx(),0);
        dparam.phaseZer.setMatAtPosit(dparam.psf.getPointerky(),1);
        dparam.phaseZer.setMatAtPosit(dparam.psf.getPointerkzOil(),2);
            
        phasePosit=new int[dparam.param.size][dparam.param.size];
        for (int i=0;i<dparam.param.size;i++){
            for (int ii=0;ii<dparam.param.size;ii++){
                phasePosit[i][ii]=-1;//outside
            }
        }
        
        for (int i=0;i<dparam.param.sizeDisk;i++){
            phasePosit[dparam.param.disk2D[i][0]][dparam.param.disk2D[i][1]]=i;
        }
        
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
        
        double [][][] stackGS = new double [image[0].length][image[0][0].length][image[0][0][0].length];
        for (int i=0;i<image[0].length;i++){
            for (int ii=0;ii<image[0][0].length;ii++){
                for (int iii=0;iii<image[0][0][0].length;iii++){
                    stackGS[i][ii][iii]=(image[0][i][ii][iii]-paramImage.B[0])/paramImage.A[0][i];
                }
            }
        }
        gs=new GS_(stackGS,dparam.param,axialside);
        
        this.path_calibration=path_calibration;
    }
    
    
    
    
    
    public void run(int nbIter){
        
        //IJ.log("zstep "+dparam.param.zstep);
        
        //this.phase_retrieve_zernike_cross_validation(nbIter);
        
        IJ.log("size "+dparam.param.sizeDisk+"  "+dparam.param.size+"  "+dparam.param.sizeRadiusRingPixel);
        
        this.phase_retrieve(nbIter);
        
        if (path_calibration.length()>2){
            dparam.save(path_calibration);
            dparam.saveJSON(path_calibration);
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
    
    
    
//    private void phase_retrieve_cross_validation(int iterations){
//        
//        IJ.log("phase retrieval cross validation");
//        
//        this.dparam.param.sigmaGaussianKernel=this.sigma;///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
//        
//        
//        
//        double likelihood=initialization(iterations);
//        
//        IJ.log("calibration started");
//        
//        
//        int nbStackUse4PR=1;
//        
//        nbStackUse4PR=Math.min(nbStackUse4PR, nbstack);
//                
//        if (true){
//
//            
//            
//            loop1:for (int t=0;t<iterations;t++){
//                
//                IJ.log("remaining iterations: "+(t)+"/"+iterations);
//                //IJ.log("pass 2 ; remaining iterations: "+(iterations-t));
//                IJ.showProgress((1.-(double)(iterations-t)/iterations));
//                long tim=System.currentTimeMillis();
//
//                
//
//                this.updateRegistrationStacks(nbstack);
//                
//                
//                this.updateSigmaGaussianKernel(nbStackUse4PR);///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                
//                //IJ.log("sigma "+dparam.param.sigmaGaussianKernel);
//                
//                //this.updatePhotonB_0(nbstack);
//                
//                
//                this.updatePhotonBpoly(nbstack,2);
//                
//                
//                
//                this.updatePhotonA(nbstack);
//                
//                
//                
//                
//                this.updatePhotonAeach(nbstack);
//                
//                for (int k=0;k<paramImage.A[0].length;k++){
//
//                    //IJ.log("A "+(paramImage.A[0][k]+this.fit_a[0]+this.fit_a_each[0][k]));
//                }
//                
//                
//                
//                this.updatePhaseNonZernike(nbStackUse4PR);
//                
//                
//                //this.updateDrift(nbstack);
//                
//                
//                //IJ.log("drift "+dx+"  "+dy+"  "+dz);
//                
//                //this.updateWeightZ(nbstack);
//                
//                for (int ss=0;ss<nbstack;ss++){
//                    //IJ.log("registration  x:"+this.registrationStack[0][ss]+"  y;"+this.registrationStack[1][ss]+"  z:"+this.registrationStack[2][ss]);
//                }
//                
//                //IJ.log("lik : "+this.getLikelihood(nbstack));
//                double lik=this.getLikelihood(nbstack);
//                
//                if (Math.abs(likelihood-lik)<epsilon){
//                    break loop1;
//                }
//                
//                likelihood=lik;
//                
//                
//
//                //IJ.log("sig "+dparam.param.sigmaGaussianKernel+"   wz "+dparam.param.getweightZ());
//                
//
//            }
//            
//            
//            
//            
//            
//        }
//        
//        
//        IJ.log("total number of image for training "+nbStackUse4PR);
//        IJ.log("total number of image "+nbstack);
//        
//        showPhase();
//        showImageAndModel();
//        double r=0;
//        for (int i=0;i<image.length;i++){
//            double d=computeResidual(i);
//            IJ.log("residual ("+i+") : "+d);
//            r+=d;
//        }
//        IJ.log("global residual : "+r);
//        
//        
//        IJ.showProgress(0);
//        
//    }
        
    
    private void phase_retrieve(int iterations){
        //this.sigma=0.;///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        this.dparam.param.sigmaGaussianKernel=this.sigma;///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
        
        
        
        
        initializationGS(100);
        double likelihood=Double.POSITIVE_INFINITY;
        //double likelihood=initialization2(10);
        
        
        IJ.log("calibration started");
        
        //second pass...all the stacks
        
        if (true){

            
            loop1:for (int t=0;t<iterations;t++){
                
                IJ.log("iteration: "+(t)+"/"+iterations);
                //IJ.log("pass 2 ; remaining iterations: "+(iterations-t));
                IJ.showProgress((1.-(double)(iterations-t)/iterations));
                


                this.updateRegistrationStacks(nbstack);
                
                this.updatePhaseNonZernike(nbstack);
                
                
                //this.updatePhotonB_0(nbstack);
                
                
                this.updatePhotonBpoly(nbstack,2);
                
                
                
                this.updatePhotonA(nbstack);
                
                
                
                
                this.updatePhotonAeach(nbstack);
                
                for (int k=0;k<paramImage.A[0].length;k++){

                    //IJ.log("A "+(paramImage.A[0][k]+this.fit_a[0]+this.fit_a_each[0][k]));
                }
                
                
                
                
                this.updateSigmaGaussianKernel(nbstack);///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                IJ.log("sigma "+dparam.param.sigmaGaussianKernel);
                
                
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
                
                if (t>0&&(t%10==0||t==iterations-1)){
                    smoothPhase();
                    showPhase("phase smoothed iter "+t);
                    showImageAndModel("psf with smoothed phase iter "+t);
                }
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
        showImageAndModel("input-model image");
    }
    void showImageAndModel(String figname){
        
        
        double [][][] modconcat = new double[image.length*image[0].length][image[0][0].length][image[0][0][0].length];
        double [][][] imageconcat = new double[image.length*image[0].length][image[0][0].length][image[0][0][0].length];
        double themin=Double.POSITIVE_INFINITY;
        double themax=Double.NEGATIVE_INFINITY;
        this.computeModel3D(0);
        
            
        
        themin=Double.POSITIVE_INFINITY;
        themax=Double.NEGATIVE_INFINITY;
        double [][][][] mod = new double [model3D.length][][][];
        for (int i=0;i<model3D.length;i++){
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
        
        ImageShow.imshow(imageconcat,modconcat,""+figname);
        IJ.setMinAndMax(themin, themax);
        
        
    }
    
    
    void showPhase(){
        showPhase("Phase");
    }
    void showPhase(String namefig){
        
        double [][] ph=new double[dparam.param.size][dparam.param.size];
        for (int p=0;p<dparam.param.sizeDisk;p++){
            double v=dparam.phaseNonZer.getValuePixel(p);
            ph[dparam.param.disk2D[p][0]][dparam.param.disk2D[p][1]]=v;
        }
                
        double themin=Double.POSITIVE_INFINITY;
        double themax=Double.NEGATIVE_INFINITY;
        for (int i=0;i<ph.length;i++){
            for (int ii=0;ii<ph[i].length;ii++){
                //ph[i][ii]=(float)((double)ph[i][ii]%(2.*Math.PI));
                if (themin>ph[i][ii]){
                    themin=ph[i][ii];
                }
                if (themax<ph[i][ii]){
                    themax=ph[i][ii];
                }
            }
        }
        ImageShow.imshow(ph,""+namefig);
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
    
    
    void initializationGS(int iterations){
        for (int z=0;z<nbstack;z++){
            IJ.log("registration("+z+")  "+this.registrationStack[0][z]+"  "+this.registrationStack[1][z]+"  "+this.registrationStack[2][z]);
        }
        gs.run(iterations);
        this.updatePhotonBpoly(nbstack,2);
        this.updatePhotonAeach(nbstack);
        for (int z=0;z<nbstack;z++){
            IJ.log("registration("+z+")  "+this.registrationStack[0][z]+"  "+this.registrationStack[1][z]+"  "+this.registrationStack[2][z]);
        }
        this.updateRegistrationStacks(nbstack);
        for (int z=0;z<nbstack;z++){
            IJ.log("registration("+z+")  "+this.registrationStack[0][z]+"  "+this.registrationStack[1][z]+"  "+this.registrationStack[2][z]);
        }
        this.updatePhotonBpoly(nbstack,2);
        this.updatePhotonAeach(nbstack);
        this.updatePhotonBpoly(nbstack,2);
        this.updatePhotonAeach(nbstack);
        
        double [][] phase= gs.getPhase();
        for (int i=0;i<dparam.param.disk2D.length;i++){
            int j=dparam.param.disk2D[i][0];
            int jj=dparam.param.disk2D[i][1];
            phase[j][jj]=phase[j][jj]%(2*Math.PI);
            
            dparam.phaseNonZer.setValuePixel(i, phase[j][jj]);
        }
        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
        gs.free();
        
        
        for (int u=0;u<iterations;u++){
            this.updatePhotonBpoly(nbstack,2);
            this.updatePhotonAeach(nbstack);
            this.updateRegistrationStacks(nbstack);
            for (int z=0;z<nbstack;z++){
            IJ.log("registration("+z+")  "+this.registrationStack[0][z]+"  "+this.registrationStack[1][z]+"  "+this.registrationStack[2][z]);
        }
        }
        showPhase("Gerchberg-Saxton initialization");
        this.showImageAndModel("Gerchberg-Saxton initialization");
        
        
        
        
        
    }
    
    
    double initialization2(int iterations){
        
        for (int u=0;u<iterations;u++){
            int angleSplit=4+(int)(Math.random()*58);
            int distSplit=20+(int)(Math.random()*10);
            this.updatePhaseNonZernike(nbstack,angleSplit,distSplit);
            this.updatePhotonBpoly(nbstack,2);
            this.updatePhotonAeach(nbstack);
            this.updateRegistrationStacks(nbstack);
        }
        showPhase("second step initialization");
        this.showImageAndModel("second step initialization");
        double likelihood = this.getLikelihood(nbstack);
        
        IJ.log("initialization 7/7");
        
        
        
        return likelihood;
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
        
        
        IJ.log("initialization 1/7");
        
        this.updatePhaseNonZernike(1,8,2);
        for (int u=0;u<10;u++){
            this.updateRegistrationStacks(1);
        }
        for (int u=0;u<100;u++){
            int angleSplit=4+(int)(Math.random()*3);
            int distSplit=1;//2+(int)(Math.random()*3);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
            
        }
        
        IJ.log("z"+registrationStack[2][0]);
        
        
        IJ.log("initialization 2/7");
        showPhase();
        this.showImageAndModel();
        
        for (int u=0;u<20;u++){
            int angleSplit=4+(int)(Math.random()*12);
            int distSplit=1;//2+(int)(Math.random()*12);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        
        IJ.log("initialization 3/7");
        showPhase();
        this.showImageAndModel();
        
        for (int u=0;u<20;u++){
            int angleSplit=4+(int)(Math.random()*15);
            int distSplit=1;//2+(int)(Math.random()*15);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        showPhase();
        this.showImageAndModel();
        
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        
        
        for (int u=0;u<15;u++){
            int angleSplit=4+(int)(Math.random()*20);
            int distSplit=2+(int)(Math.random()*20);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        IJ.log("z"+registrationStack[2][0]);
        
        showPhase();
        this.showImageAndModel();
        
        for (int u=0;u<10;u++){
            int angleSplit=10+(int)(Math.random()*25);
            int distSplit=8+(int)(Math.random()*14);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        
        
        showPhase();
        this.showImageAndModel();
        
        this.updatePhotonBpoly(nbstack,2);
        this.updatePhotonAeach(nbstack);
        this.updateRegistrationStacks(nbstack);
        
        IJ.log("z"+registrationStack[2][0]);
        for (int u=0;u<10;u++){
            int angleSplit=4+(int)(Math.random()*58);
            int distSplit=20+(int)(Math.random()*12);
            this.updatePhaseNonZernike(nbstack,angleSplit,distSplit);
            this.updatePhotonBpoly(nbstack,2);
            this.updatePhotonAeach(nbstack);
            this.updateRegistrationStacks(nbstack);
        }
        showPhase();
        this.showImageAndModel();
        
        
        likelihood = this.getLikelihood(nbstack);
        
        IJ.log("initialization 7/7");
        
        
        
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
    
    
    
    
    
    
    
    
    
    
    private void smoothPhase(){
        double [] phase= dparam.phaseNonZer.getValuesPhase();
        double [] phaseMean= dparam.phaseNonZer.getValuesPhase();
        int size=dparam.param.size;
        int [][] ph=new int[size][size];
        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                ph[i][ii]=-1;
            }
        }
        for (int p=0;p<dparam.param.sizeDisk;p++){
            ph[dparam.param.disk2D[p][0]][dparam.param.disk2D[p][1]]=p;
        }
        for (int i=0;i<size;i++){
            for (int ii=0;ii<size;ii++){
                if (ph[i][ii]>=0){
                    phaseMean[ph[i][ii]]=0;
                    double count=0;
                    for (int j=-1;j<=1;j++){
                        for (int jj=-1;jj<=1;jj++){
                            if (((i+j)>=0)&&((i+j)<size)&&((ii+jj)>=0)&&((ii+jj)<size)){
                                if (ph[i+j][ii+jj]>=0){
                                    double shift=0;
                                    double minDist1=Math.abs(phase[ph[i][ii]]-phase[ph[i+j][ii+jj]]);
                                    double minDist2=Math.abs(phase[ph[i][ii]]-phase[ph[i+j][ii+jj]]-Math.PI*2);
                                    double minDist3=Math.abs(phase[ph[i][ii]]-phase[ph[i+j][ii+jj]]+Math.PI*2);
                                    if (minDist2<minDist1){
                                        shift=Math.PI*2;
                                    }
                                    if (minDist3<minDist1){
                                        shift=-Math.PI*2;
                                    }
                                    
                                    phaseMean[ph[i][ii]]+=phase[ph[i+j][ii+jj]]+shift;
                                    count++;
                                }
                            }
                        }
                    }
                    phaseMean[ph[i][ii]]/=count;
                }
            }
        }
        dparam.phaseNonZer.setValuesPhase(phaseMean);
    }
    
    
    
    
    
    
    
    
    //compute each pixel of the phase
    private void updatePhaseNonZernike(int stackNumber){
        
        
        
        double h=.001;
        
        
        double lastLikelihood=-1;
        
       
        for (int p=0;p<dparam.param.sizeDisk;p++){
            
            if (true){
                
                double save=dparam.phaseNonZer.getValuePixel(p);
                
                
                double lik1=0;
                for (int z=0;z<stackNumber;z++){
                    
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    
                    
                    dparam.phaseNonZer.setValuePixel(p, save-h);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");
                            
                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model1");
                    lik1+=model3D[z].getLikelihood(methodLikelihood);
                }
                
                //IJ.log("lik1 "+lik1+"  "+getPhaseRegularization(p,save-h));
                lik1+=lambda*getPhaseRegularization(p,save-h);
                

                
                
                
                
                double lik2=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.phaseNonZer.setValuePixel(p, save);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");
                    
                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model2");
                    lik2+=model3D[z].getLikelihood(methodLikelihood);
                }
                lik2+=lambda*getPhaseRegularization(p,save);
                
                
                double lik3=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.phaseNonZer.setValuePixel(p, save+h);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                    compute2DModelMany_f(z);
                    
                    this.computeModel3D(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                    //ImageShow.imshow(model3D[z].getModel(),"model3");
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                }
                lik3+=lambda*getPhaseRegularization(p,save+h);
                
                
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                if (grad>Math.PI/20.){
                    grad=Math.PI/20.;
                }
                if (grad<-Math.PI/20.){
                    grad=-Math.PI/20.;
                }
                
                boolean found=false;
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    
                    double lik=0;
                    for (int z=0;z<stackNumber;z++){
                        for (int pp=0;pp<3;pp++){
                            dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                        }
                        
                        dparam.phaseNonZer.setValuePixel(p, save-gamma*grad);
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                        compute2DModelMany_f(z);
                        this.computeModel3D(z);
                        lik+=model3D[z].getLikelihood(methodLikelihood);
                    }
                    lik+=lambda*getPhaseRegularization(p,save-gamma*grad);
                    
                    if (lik<lik2){
                        lastLikelihood=lik;
                        found=true;
                        break loop;
                    }
                    else{
                        dparam.phaseNonZer.setValuePixel(p, save);
                        
                    }
                }
                if (!found){
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    
                    
                }
            }
        }
        
        
    }
    
    
    
    //compute small portions of the phase
    private void updatePhaseNonZernike(int stackNumber,int angleSplit,int distSplit){
        
        int splitNumber=angleSplit*distSplit;
        
        if ((dparam.phaseNonZer.angleNumber!=angleSplit)||(dparam.phaseNonZer.distNumber!=distSplit)){
            dparam.phaseNonZer.createSplits(angleSplit, distSplit);
        }
        
        
        double h=.001;
        
        
        double lastLikelihood=-1;
        
       
        for (int p=0;p<splitNumber;p++){
            
            double [] save=dparam.phaseNonZer.getValuesPhase(p);
            
            if (save!=null){
                
                
                
                double [] tmp=new double [save.length];
                
                double lik1=0;
                for (int z=0;z<stackNumber;z++){
                    
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    
                    for (int i=0;i<tmp.length;i++){
                        tmp[i]=save[i]-h;
                    }
                    dparam.phaseNonZer.setValuesPhase(p, tmp);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    
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
                    dparam.phaseNonZer.setValuesPhase(p, save);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
                    for (int i=0;i<tmp.length;i++){
                        tmp[i]=save[i]+h;
                    }
                    dparam.phaseNonZer.setValuesPhase(p, tmp);
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph3"+p); 
                    compute2DModelMany_f(z);
                    
                    this.computeModel3D(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image3","DOUBLE");
                    //ImageShow.imshow(model3D[z].getModel(),"model3");
                    lik3+=model3D[z].getLikelihood(methodLikelihood);
                }
                
                
                
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                if (grad>Math.PI/20){
                    grad=Math.PI/20;
                }
                if (grad<-Math.PI/20.){
                    grad=-Math.PI/20.;
                }
                
                boolean found=false;
                loop:for (double gamma=1;gamma>.00002;gamma/=10){
                    
                    double lik=0;
                    for (int z=0;z<stackNumber;z++){
                        for (int pp=0;pp<3;pp++){
                            dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                        }
                        for (int i=0;i<tmp.length;i++){
                            tmp[i]=save[i]-gamma*grad;
                        }
                        dparam.phaseNonZer.setValuesPhase(p, tmp);
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
                        //IJ.write("old / new "+lik2+"   "+lik+"   "+grad);
                        dparam.phaseNonZer.setValuesPhase(p, save);
                    }
                }
                if (!found){
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    
                    
                }
                else{
                    
                }
            }
        }
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    /*
    
    private void updatePhaseZernike(int stackNumber){
        updatePhaseZernike(stackNumber,-1);
    }
    
    
    private void updatePhaseZernike(int stackNumber,int coefNumber){
        
        if (coefNumber<3){
            coefNumber=dparam.phaseZer.numCoef;
        }
        
        if (coefNumber>dparam.phaseZer.numCoef){
            coefNumber=dparam.phaseZer.numCoef;
        }
        
        double h=.01;
        
        
        double lastLikelihood=-1;
        
       
        for (int p=3;p<coefNumber;p++){
        
            if (true){
            
                
                double lik1=0;
                for (int z=0;z<stackNumber;z++){
                    
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase(),p, -h));
                    
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
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
                    
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase(),p, h));
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
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    
                    
                }
            }
        }
        
        
    }
    
    */
    
    
    
    
    
    
    
    
    
    
    private void updateDrift(int stackNumber){
        
        
        double h=.00001;
        
        
        double lastLikelihood=-1;
        
       
            
        {
            double lik1=0;
            for (int z=0;z<stackNumber;z++){

                for (int pp=0;pp<3;pp++){
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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

                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));

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
                    
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase(),p, -h));

                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph1 "+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image1","DOUBLE");

                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model1");
                    lik1+=model3D[z].getLikelihood(methodLikelihood);
                    


                    




                    double lik2=0;
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
                    //dparam.psf_fMany.imshowFloat(dparam.psf_fMany.device_phase, "ph2"+p); 
                    compute2DModelMany_f(z);
                    //dparam.psf_fMany.imshow(dparam.param.sizeoutput*dparam.param.sizeoutput*nbSlice, dparam.param.sizeoutput,dparam.psf_fMany.getPointerPSF(),"image2","DOUBLE");

                    this.computeModel3D(z);
                    //ImageShow.imshow(model3D[z].getModel(),"model2");
                    lik2+=model3D[z].getLikelihood(methodLikelihood);
                    
                    


                    double lik3=0;
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase(),p, h));
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
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
                        dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));


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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
            compute2DModelMany_f(z);
            this.computeModel3D(z);
            lik2+=model3D[z].getLikelihood(methodLikelihood);
        }
        //IJ.log("lik "+lik2);
        

        double lik3=0;
        dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel+h);
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
            double shift=0;
            if (gamma*grad>0){
                shift=Math.min(gamma*grad,0.05);
            }
            else{
                shift=Math.max(gamma*grad,-0.05);
            }
            dparam.param.sigmaGaussianKernel=save-gamma*grad;
            dparam.psf_fMany.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
            double lik=0;
            for (int z=0;z<stackNumber;z++){
                
                for (int pp=0;pp<3;pp++){//update due to registration of stacks
                    dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                }
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
        //IJ.log("lik "+lastLikelihood);
            
        
    }
    
    
    
    
    
    
    
    /*
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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
                dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
    */
    
    
    
    
    
    
    
    private void updatePhotonA(int stackNumber){
        
            double h=1;

            
            
            
                for (int z=0;z<stackNumber;z++){
                    
                    double [] theA = new double [paramImage.A[z].length];
                    
                    
                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
//                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
//                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            compute2DModelMany_f(z);
            
            for (int o=0;o<=order;o++){
                
                int [] index=this.pbg[z].getIndexOfOrder(o);
                for (int ind=0;ind<index.length;ind++){
                    
                    
                    
                    
                    double saveB=pbg[z].a[index[ind]];

                    for (int pp=0;pp<3;pp++){//update due to registration of stacks
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
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
    
    
    
    
    
    private double getPhaseRegularization(int posit,double currentValue){
        
        int x=this.dparam.param.disk2D[posit][0];
        int y=this.dparam.param.disk2D[posit][1];
        
        
        
        //deal with edge effet --> symmetry if outsize of the disk
        int shiftx=1;
        if (this.phasePosit[x+shiftx][y]<0){
            shiftx=-1;
        }
        int shifty=1;
        if (this.phasePosit[x][y+shifty]<0){
            shifty=-1;
        }
        
        
        double X=this.dparam.phaseNonZer.getValuePixel(this.phasePosit[x+shiftx][y]);
        double Y=this.dparam.phaseNonZer.getValuePixel(this.phasePosit[x][y+shifty]);
        
        double reg=Math.sqrt((X-currentValue)*(X-currentValue)+(Y-currentValue)*(Y-currentValue));
        
        return reg;
        
    }
    
    
    
    private double  getLikelihood(int stackNumber){
        
        
        
        
        double likelihood=0;
        
       
            
        
        for (int z=0;z<stackNumber;z++){
            
            for (int pp=0;pp<3;pp++){//update due to registration of stacks
                dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
            }
            dparam.psf_fMany.updatePhase(dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPointerPhase()));
            
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
 