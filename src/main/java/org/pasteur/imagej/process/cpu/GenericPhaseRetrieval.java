/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;
import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.process.*;


import ij.IJ;
import java.awt.Color;
import ij.gui.Plot;


/**
 *
 * @author benoit
 */
public class GenericPhaseRetrieval {
    
    int cpt=0;
    int idStack=0;
    
    
    double [] deltaZ;
    
    double [][][][] image;
    int methodLikelihood=0;
    int center;
    InitBackgroundAndPhotonNumber paramImage;
    DataPhase dparam;
    //public Model3DJCudaFastDouble [] model3D;
    
    int nbProcess;
    
    double [][] x;
    double [][] y;
    double [][] z;
    int nbSlice;
    int nbstack;
    
    double imageLength;
    
    Model [] model;
    Object [] monitor1;
    Object monitor0;
    boolean [] toBeblocked;
    
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
    double zwat=0 ;
    
    double dx=0;//drift X supposed linear
    double dy=0;//drift Y supposed linear
    double dz=0;//drift Z supposed linear
    
    double maxDrift=10;//max drift between each frame (nm)
    
    public GenericPhaseRetrieval(int sizeFFT,double xystep,double zstep,double wavelength,double noil,double na,InitBackgroundAndPhotonNumber paramImage,String path_calibration,double sigma,int axialside){
        this.sigma=sigma;
        this.image=paramImage.image;
        
        nbstack=image.length;
        nbSlice=image[0].length;
        imageLength=image[0][0].length*image[0][0][0].length;
        
        nbProcess=nbSlice;
        
        dparam = new DataPhase(sizeFFT,image[0][0].length,0,xystep,zstep,wavelength,noil,na,1.0);
        dparam.param.zernikedPSF=false;
        dparam.param.Zfocus=0;
        dparam.phaseZer.setMatAtPosit(dparam.psf.getKxPointer(),0);
        dparam.phaseZer.setMatAtPosit(dparam.psf.getKyPointer(),1);
        dparam.phaseZer.setMatAtPosit(dparam.psf.getKzPointer(),2);
        dparam.setNwat(dparam.param.noil);//ca ne change rien normalement car bille collée à lamelle
        dparam.setMany(nbProcess);
        monitor1 = new Object[nbProcess];
        monitor0 = new Object();
        toBeblocked=new boolean [nbProcess];
        model= new Model[nbProcess];
        for (int i=0;i<nbProcess;i++){
            monitor1[i] = new Object();
            toBeblocked[i]=true;
            model[i]=new Model(monitor1[i],i);
        }
        
        
        
        for (int k=0;k<nbProcess;k++){
            synchronized(monitor1[k]){
                model[k].start();
                try{
                    monitor1[k].wait();//to be sure the threads are launched at this stage and that wait() in the thread is called
                }catch(Exception er){IJ.log("oops: problem wait function for sync monitor 1");}
            }
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
        
        
        
        x=new double[nbstack][nbSlice];//registration of each stack
        y=new double[nbstack][nbSlice];//registration of each stack
        z=new double[nbstack][nbSlice];//registration of each stack
        
        
        
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
        
        //model3D=new Model3DJCudaFastDouble[image.length];
        
        
        background= new double[nbstack][image[0][0].length][image[0][0][0].length];
        pbg=new PolynomialBackground[nbstack];
        for (int z=0;z<nbstack;z++){
            //model3D[z]=new Model3DJCudaFastDouble(dparam.param,this.image[z]);
            pbg[z]= new PolynomialBackground();
            
            //init background fixed value
            int [] ind=pbg[z].getIndexOfOrder(0);
            pbg[z].a[ind[0]]=paramImage.B[z];
            
            computeBackground2D(z);
            
            
        }
        
        

        
        this.path_calibration=path_calibration;
    }
    
    
    
    public void run(int nbIter){
        
        //IJ.log("zstep "+dparam.param.zstep);
        
        this.phase_retrieve(nbIter);
        
        if (path_calibration.length()>2){
            dparam.save(path_calibration);
            IJ.log("calibration file saved");
        }
        else{
            IJ.log("calibration file not saved...you should select a path to save it");
        }
        
        
        
    }
    
    
    
    
    
    
    
    
    
    
    private void phase_retrieve(int iterations){
        
        
        
        this.dparam.param.sigmaGaussianKernel=this.sigma;///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dparam.psf_many.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
        
        
        
        
        double likelihood=initialization(iterations);
        
        
        
        
        IJ.log("calibration started");
        
        if (true){
            
            
            likelihood=0;
                
                
            for (int i=0;i<image.length;i++){
                likelihood+=this.computePSFLik(i);
            }


            
            
            loop1:for (int t=0;t<iterations;t++){
                
                IJ.log("iteration: "+(t)+"/"+iterations);
                //IJ.log("pass 2 ; remaining iterations: "+(iterations-t));
                IJ.showProgress((1.-(double)(iterations-t)/iterations));
                

                

                this.updateRegistrationStacks(nbstack);
                
                
                this.updateSigmaGaussianKernel(nbstack);///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                //this.updatePhotonB_0(nbstack);
                
                
                this.updatePhotonBpoly(nbstack,2);
                
                
                
                this.updatePhotonA(nbstack);
                
                
                
                
                this.updatePhotonAeach(nbstack);
                
                for (int k=0;k<paramImage.A[0].length;k++){

                    //IJ.log("A "+(paramImage.A[0][k]+this.fit_a[0]+this.fit_a_each[0][k]));
                }
                
                
                
                this.updatePhaseNonZernike(nbstack);
                
                
                //this.updateDrift(nbstack);
                
                
                //IJ.log("drift "+dx+"  "+dy+"  "+dz);
                
                //this.updateWeightZ(nbstack);
                
                for (int ss=0;ss<nbstack;ss++){
                    //IJ.log("registration  x:"+this.registrationStack[0][ss]+"  y;"+this.registrationStack[1][ss]+"  z:"+this.registrationStack[2][ss]);
                }
                
                //IJ.log("lik : "+this.getLikelihood(nbstack));
                double lik=0;
                
                
                for (int i=0;i<image.length;i++){
                    lik+=this.computePSFLik(i);
                }
                
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
        
        
        for (int k=0;k<nbProcess;k++){
            model[k].kill();
            synchronized(monitor1[k]){
                monitor1[k].notify();
            }
            try{
                model[k].join();
            }catch(Exception eee){IJ.log("Thread impossible to join() "+eee);}
        }
        
    }
    
   
    
    
    void showPhase(){
        
            
        double [][] ph=dparam.psf_many.getPhase();
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
        this.computePSFLik(image_id);
        double [][][] mod=new double [image[image_id].length][][];
        for (int i=0;i<image[image_id].length;i++){
            mod[i]=model[i].getModelPointer();
        }
        
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
    
    
    
    void showImageAndModel(){
        showImageAndModel("");
    }
    void showImageAndModel(String title){
        
        
        double [][][] modconcat = new double[image.length*image[0].length][image[0][0].length][image[0][0][0].length];
        double [][][] imageconcat = new double[image.length*image[0].length][image[0][0].length][image[0][0][0].length];
        double themin=Double.POSITIVE_INFINITY;
        double themax=Double.NEGATIVE_INFINITY;
        this.computePSFLik(0);
        
        
        
        
        themin=Double.POSITIVE_INFINITY;
        themax=Double.NEGATIVE_INFINITY;
        double [][][][] mod = new double [image.length][][][];
        for (int i=0;i<image.length;i++){
            this.computePSFLik(i);
            mod[i] = new double [image[i].length][][];
            for (int ii=0;ii<image[i].length;ii++){
                
                mod[i][ii]=model[ii].getModelPointer();
            }
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
        
        ImageShow.imshow(imageconcat,modconcat,"input/model image "+title);
        IJ.setMinAndMax(themin, themax);
        
        
    }
    
    
    
    
    double initialization(int iterations){
        
        IJ.log("initialization started");
        
        IJ.log("This program does not work...");
        
        double [][][][] stack;
        
        stack=this.image;
        
        
        
        
        int stack_size = stack[0].length;
        
        double stack_center=center;
        
        
        double epsilon=.0001;
        
        
        
        
        int maxorder=1;
        
        //first pass...ont stack
        double likelihood = Double.MAX_VALUE;
        
        
        IJ.log("initialization 1/6");
        
        this.updatePhaseNonZernike(1,8,2);
        for (int u=0;u<10;u++){
            this.updateRegistrationStacks(1);
        }
        for (int u=0;u<20;u++){
            int angleSplit=4+(int)(Math.random()*3);
            int distSplit=2+(int)(Math.random()*3);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
            
        }
        
        
        
        IJ.log("initialization 2/6");
        showPhase();
        this.showImageAndModel();
        
        for (int u=0;u<20;u++){
            int angleSplit=4+(int)(Math.random()*12);
            int distSplit=2+(int)(Math.random()*12);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        
        
        IJ.log("initialization 3/6");
        showPhase();
        this.showImageAndModel();
        
        for (int u=0;u<20;u++){
            int angleSplit=4+(int)(Math.random()*15);
            int distSplit=2+(int)(Math.random()*15);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        showPhase();
        this.showImageAndModel();
        
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        
        
        IJ.log("initialization 4/6");
        
        for (int u=0;u<15;u++){
            int angleSplit=4+(int)(Math.random()*20);
            int distSplit=2+(int)(Math.random()*20);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
            this.updateRegistrationStacks(1);
        }
        
        this.updatePhotonBpoly(1,2);
        this.updatePhotonAeach(1);
        
        
        IJ.log("initialization 5/6");
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
        
        
        IJ.log("initialization 6/6");
        
        for (int u=0;u<5;u++){
            int angleSplit=4+(int)(Math.random()*58);
            int distSplit=20+(int)(Math.random()*12);
            this.updatePhaseNonZernike(1,angleSplit,distSplit);
        }
        
        showPhase();
        this.showImageAndModel();
        
        this.updatePhotonBpoly(nbstack,2);
        this.updatePhotonAeach(nbstack);
        this.updateRegistrationStacks(nbstack);
        
        
        
        
        
        
        double lik=0;
        for (int i=0;i<nbstack;i++){
            lik+=this.computePSFLik(i);
        }
        
        
        showPhase();
        
        
        return likelihood;
    }
    
    
    
    private void updatePhotonAeach(int stackNumber){
        
        double h=1;

        for (int z=0;z<stackNumber;z++){
            
            computePSFOnly(z);
            
            double [] saveA=new double[paramImage.A[z].length];
            
            for (int i=0;i<paramImage.A[z].length;i++){
                
                saveA[i]=fit_a_each[z][i];
                
            }
                   
            double [] lik1;
            double [] lik2;
            double [] lik3;
            double [] lik;
            
            for (int a=0;a<saveA.length;a++){
                fit_a_each[z][a]=saveA[a]-h;
            }
                
            lik1=this.computeLikOnlyStack(z);

            
            for (int a=0;a<saveA.length;a++){
                fit_a_each[z][a]=saveA[a];
            }
                
            lik2=this.computeLikOnlyStack(z);

            
            
            for (int a=0;a<saveA.length;a++){
                fit_a_each[z][a]=saveA[a]+h;
            }
                
            lik3=this.computeLikOnlyStack(z);



            double [] grad = new double [saveA.length];
            for (int a=0;a<saveA.length;a++){
                if (Math.abs((lik3[a]+lik1[a]-2*lik2[a])/(h*h))==0){
                    grad[a]=((lik3[a]-lik1[a])/(2*h));
                }
                else{
                    grad[a]=((lik3[a]-lik1[a])/(2*h))/(Math.abs((lik3[a]+lik1[a]-2*lik2[a])/(h*h)));
                }
                //IJ.log("lik(a)  "+a+"  "+lik1[a]+"  "+lik2[a]+"  "+lik3[a]+"  "+grad[a]);
            }

                
            
            boolean [] found= new boolean [saveA.length];
            for (int a=0;a<saveA.length;a++){
                found[a]=false;
            }
            loop:for (double gamma=1;gamma>.002;gamma/=10){
                
                for (int a=0;a<saveA.length;a++){
                    
                    if (!found[a]){

                        fit_a_each[z][a]=saveA[a]-gamma*grad[a];
                    }
                }
                
                lik=this.computeLikOnlyStack(z);
                
                int numberNotFount=0;
                for (int a=0;a<saveA.length;a++){
                    if (!found[a]){
                        if (lik[a]<lik2[a]){

                            

                            found[a]=true;
                        }
                        else{
                            fit_a_each[z][a]=saveA[a];
                            numberNotFount++;
                        }
                    }
                }
                
                
                //IJ.log("numberNotFount "+gamma+"  "+numberNotFount);
                if (numberNotFount==0){
                    
                    break loop;

                }
            }
            this.computeLikOnlyStack(z);
        }
            
            
        
        
    }
    
    
    
    private void updatePhotonA(int stackNumber){
        
            double h=1;

            
            
            
                for (int z=0;z<stackNumber;z++){
                    
                    computePSFOnly(z);
                    
                    
                    
                    double saveA=fit_a[z];
                    
                    
                    
                    
                    fit_a[z]=saveA-h;
                    
                    double lik1=computeLikOnly(z);
                    
                    
                    fit_a[z]=saveA;
                    
                    double lik2=computeLikOnly(z);


                    fit_a[z]=saveA+h;
                    
                    double lik3=computeLikOnly(z);

                    

                    
                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/(Math.abs((lik3+lik1-2*lik2)/(h*h)));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.002;gamma/=10){
                        
                        
                        fit_a[z]=saveA-gamma*grad;
                    
                        double lik=computeLikOnly(z);
                        

                        if (lik<lik2){
                            
                            fit_a[z]=fit_a[z]-gamma*grad;
                            
                            found=true;
                            break loop;
                        }
                    }
                    if (!found){
                        fit_a[z]=saveA;
                        computeLikOnly(z);
                    }
                }
            
            
        
        
    }
    
    
    
    
    
    
    
    
    
    private void updatePhotonBpoly(int stackNumber,int order){
        
        
        
        double h=.0001;



        for (int z=0;z<stackNumber;z++){
            
            computePSFOnly(z);
            
            
            for (int o=0;o<=order;o++){
                
                int [] index=this.pbg[z].getIndexOfOrder(o);
                for (int ind=0;ind<index.length;ind++){
                    
                    
                    
                    
                    double saveB=pbg[z].a[index[ind]];

                    




                    double lik1=0;


                    pbg[z].a[index[ind]] = saveB-h;
                    this.computeBackground2D(z);

                    lik1=computeLikOnly(z);








                    double lik2=0;

                    pbg[z].a[index[ind]] = saveB;
                    

                    this.computeBackground2D(z);

                    lik2=computeLikOnly(z);




                    double lik3=0;


                    pbg[z].a[index[ind]] = saveB+h;


                    this.computeBackground2D(z);
                    
                    lik3=computeLikOnly(z);
                    



                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/(Math.abs((lik3+lik1-2*lik2)/(h*h)));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.02;gamma/=10){

                        double lik=0;


                        pbg[z].a[index[ind]] = saveB-gamma*grad;

                        this.computeBackground2D(z);
                        
                        lik=computeLikOnly(z);


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
    
    
    
    
    
    
    
    
    
    
    
    private void updateSigmaGaussianKernel(int stackNumber){
        
        
        
        double h=.001;
        
        
        double lastLikelihood=-1;
        
       
        dparam.phaseZer.computeCombination();
        dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
        
        double lik1=0;
        dparam.psf_many.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel-h);
        
        for (int z=0;z<stackNumber;z++){
            
            lik1+=computePSFLik(z);
        }
        
        
        double lik2=0;
        dparam.psf_many.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
        
        for (int z=0;z<stackNumber;z++){
            
            lik2+=computePSFLik(z);
        }
        //IJ.log("lik "+lik2);
        

        double lik3=0;
        dparam.psf_many.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel+h);
        for (int z=0;z<stackNumber;z++){
            
            lik3+=computePSFLik(z);
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
            dparam.psf_many.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
            double lik=0;
            for (int z=0;z<stackNumber;z++){
                
                lik+=computePSFLik(z);
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
            dparam.psf_many.updateSigmaGaussianKernel(dparam.param.sigmaGaussianKernel);
            
        }
        //IJ.log("lik "+lastLikelihood);
            
        
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
                    
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    lik1+=computePSFLik(z);
                    
                    
                }
                
                

                
                
                
                
                double lik2=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.phaseNonZer.setValuePixel(p, save);
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    lik2+=computePSFLik(z);
                }
                
                
                
                double lik3=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.phaseNonZer.setValuePixel(p, save+h);
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    lik3+=computePSFLik(z);
                }
                
                
                
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                if (grad>Math.PI/10.){
                    grad=Math.PI/10.;
                }
                if (grad<-Math.PI/10.){
                    grad=-Math.PI/10.;
                }
                
                boolean found=false;
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    
                    double lik=0;
                    for (int z=0;z<stackNumber;z++){
                        for (int pp=0;pp<3;pp++){
                            dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                        }
                        
                        dparam.phaseNonZer.setValuePixel(p, save-gamma*grad);
                        dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                        dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                        lik+=computePSFLik(z);
                    }
                    
                    
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
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    
                            
                    
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
                    
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    lik1+=computePSFLik(z);
                    
                }
                
                

                
                
                
                
                double lik2=0;
                for (int z=0;z<stackNumber;z++){
                    for (int pp=0;pp<3;pp++){
                        dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                    }
                    dparam.phaseNonZer.setValuesPhase(p, save);
                    
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    lik2+=computePSFLik(z);
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
                    
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    lik3+=computePSFLik(z);
                }
                
                
                
                double grad=0;
                if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                    grad=((lik3-lik1)/(2*h));
                }
                else{
                    grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
                }
                
                if (grad>Math.PI/10){
                    grad=Math.PI/10;
                }
                if (grad<-Math.PI/10.){
                    grad=-Math.PI/10.;
                }
                
                boolean found=false;
                loop:for (double gamma=1;gamma>.02;gamma/=10){
                    
                    double lik=0;
                    for (int z=0;z<stackNumber;z++){
                        for (int pp=0;pp<3;pp++){
                            dparam.phaseZer.setA(pp, this.registrationStack[pp][z]);
                        }
                        for (int i=0;i<tmp.length;i++){
                            tmp[i]=save[i]-gamma*grad;
                        }
                        dparam.phaseNonZer.setValuesPhase(p, tmp);
                        dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                        dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                        lik+=computePSFLik(z);
                    }
                    
                    
                    if (lik<lik2){
                        lastLikelihood=lik;
                        found=true;
                        break loop;
                    }
                    else{
                        dparam.phaseNonZer.setValuesPhase(p, save);
                    }
                }
                if (!found){
                    dparam.phaseZer.computeCombinationPlusOtherPhase(dparam.phaseNonZer.getPhasePointer());
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                    
                }
            }
        }
        
        
    }
    
    
    
    
    
    private void updateRegistrationStacks(int stackNumber){
        
        if (stackNumber>=1){
            double h=.01;


            


            for (int p=0;p<3;p++){//first 3 coefs: used for X, Y, Z shifts
                
                
                
                for (int z=0;z<stackNumber;z++){
                    
                    double save=this.registrationStack[p][z];
                    
                    this.registrationStack[p][z]=save-h;
                    
                    double lik1=computePSFLik(z);
                    
                    
                    


                    


                    this.registrationStack[p][z]=save;
                    
                    double lik2=computePSFLik(z);
                    
                    



                    this.registrationStack[p][z]=save+h;
                    
                    double lik3=computePSFLik(z);
                    
                    
                    

                    //IJ.log("save "+dparam.phaseZer.getA(p)+"  "+this.registrationStack[p][z]);
                    
                    double grad=0;
                    if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                        grad=((lik3-lik1)/(2*h));
                    }
                    else{
                        grad=((lik3-lik1)/(2*h))/(Math.abs((lik3+lik1-2*lik2)/(h*h)));
                    }

                    //IJ.log("lik  "+lik1+"  "+lik2+"  "+lik3+"  "+grad);

                    boolean found=false;
                    loop:for (double gamma=1;gamma>.02;gamma/=10){

                        
                        this.registrationStack[p][z]=save-grad;
                        double lik=computePSFLik(z);
                        
                        if (lik<lik2){
                            
                            found=true;
                            break loop;
                        }
                    }
                    if (!found){
                        this.registrationStack[p][z]=save;
                        computePSFLik(z);


                    }
                }
            }
            
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
        
        
        double [] init=dparam.phaseZer.getAPointer();
        
        
        
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
                coco[i][b]=valstep[id[bb]];
            }
            for (int b=endcoefNumber;b<dparam.phaseZer.numCoef;b++){
                coco[i][b]=init[b];
            }
        }
        
        
        int bestId=0;
        
        for (int p=0;p<coco.length;p++){
            
            
                
            double lik=0;
            if (p==0){
                dparam.phaseZer.setA(coco[p]);
                dparam.phaseZer.computeCombination();
                dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
            }
            else{
                int nbChange=0;
                int idChange=0;
                for (int i=0;i<coco[p].length;i++){
                    if (coco[p][i]!=coco[p-1][i]){
                        nbChange++;
                        idChange=i;
                    }
                }
                if (nbChange==1){
                    
                    dparam.phaseZer.setA(idChange, coco[p][idChange]);
                    dparam.phaseZer.computeCombination();
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                }
                else{
                    dparam.phaseZer.setA(coco[p]);
                    dparam.phaseZer.computeCombination();
                    dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
                }
            }
            
            for (int z=0;z<stackNumber;z++){
                lik+=computePSFLik(z);
            }
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
        dparam.phaseZer.computeCombination();
        dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
        
        
        
        return likelihood;
    }
    
    
    
    double computePSFLik(int stack){
        
        for (int k=0;k<nbProcess;k++){
            model[k].runWithPSF();
            model[k].runWithLikelihood();
        }
        return compute(stack);
    }
    
    
    
    double computeLikOnly(int stack){
        
        for (int k=0;k<nbProcess;k++){
            model[k].runWithoutPSF();
            model[k].runWithLikelihood();
        }
        return compute(stack);
    }
    
    double [] computeLikOnlyStack(int stack){
        
        for (int k=0;k<nbProcess;k++){
            model[k].runWithoutPSF();
            model[k].runWithLikelihood();
        }
        return computeStack(stack);
    }
    
    double computePSFOnly(int stack){
        
        for (int k=0;k<nbProcess;k++){
            model[k].runWithPSF();
            model[k].runWithoutLikelihood();
        }
        return compute(stack);
    }
    
    
    
    double compute(int stack){
        double lik=0;
        switchToStack(stack);
        cpt=0;
        synchronized(monitor0){
            for (int k=0;k<nbProcess;k++){
                synchronized(monitor1[k]){
                    toBeblocked[k]=true;
                    monitor1[k].notify();
                }
            }
//            for (int k=0;k<nbProcess;k++){
//                synchronized(monitor1[k]){
//                    try{
//                        if (toBeblocked[k]){
//                            monitor1[k].wait(10000);
//                        }
//                    }catch(Exception ee){IJ.log("error wait function "+ee);}
//                }
//            }
            try{
                monitor0.wait();
            }catch(Exception ee){IJ.log("error wait function "+ee);}
        }
        for (int k=0;k<nbProcess;k++){
            lik+=model[k].getLikelihood();
        }
        return lik;
    }
    
    double [] computeStack(int stack){
        double [] lik=new double [nbProcess];
        switchToStack(stack);
        cpt=0;
        synchronized(monitor0){
            for (int k=0;k<nbProcess;k++){
                synchronized(monitor1[k]){
                    toBeblocked[k]=true;
                    monitor1[k].notify();
                }
            }
    //        for (int k=0;k<nbProcess;k++){
    //            synchronized(monitor1[k]){
    //                try{
    //                    if (toBeblocked[k]){
    //                        monitor1[k].wait(10000);
    //                    }
    //                }catch(Exception ee){IJ.log("error wait function "+ee);}
    //            }
    //        }
            try{
                monitor0.wait();
            }catch(Exception ee){IJ.log("error wait function "+ee);}
        }
        for (int k=0;k<nbProcess;k++){
            lik[k]=model[k].getLikelihood();
        }
        return lik;
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
    
    
    
    
    
    void switchToStack(int id){
        for (int pp=0;pp<3;pp++){
            dparam.phaseZer.setA(pp, registrationStack[pp][id]);
        }
        dparam.phaseZer.computeCombination();
        dparam.psf_many.updatePhasePointer(dparam.phaseZer.getPhasePointer());
        idStack=id;
    }
    
    
    
    
    
    class Model extends Thread{
        
        boolean withPSF=true;
        boolean withLik=true;
        int idSlice;
        double [][] mod;
        Object monitor1;
        boolean killer=false;
        double variable;
        int computation=0;
        double likelihood;
        Model(Object monitor1,int idSlice){
            this.monitor1=monitor1;
            this.mod=new double [dparam.param.sizeoutput][dparam.param.sizeoutput];
            this.idSlice=idSlice;
        }
        
        
        
        double getLikelihood(){
            return likelihood;
        }
        
        double [][] getModelPointer(){
            return mod;
        }
        
        
        
        public void kill(){
            killer=true;
        }
        
        
        void runWithoutPSF(){
            withPSF=false;
        }
        
        void runWithPSF(){
            withPSF=true;
        }
        
        void runWithoutLikelihood(){
            withLik=false;
        }
        
        void runWithLikelihood(){
            withLik=true;
        }
        
        public void run(){
            synchronized(monitor1) {
                monitor1.notify();
                
                loop:while (true){
                
                    
                    try{ 
                        monitor1.wait();
                    }catch(Exception ee){IJ.log("error wait function "+ee);}
                    
                    
                    if (killer){
                        break loop;
                    }
                    
                    if (withPSF){
                        computePSF();
                    }
                    if (withLik){
                        computeLikelihood();
                    }
                    
                    
//                    toBeblocked[idSlice]=false;
//                    try{
//                        Thread.sleep(30);//To optimize
//                    }catch(Exception ee){}
//                    monitor1.notify();
                    
                    synchronized(monitor0) {
                        cpt++;
                        if (cpt==nbProcess){
                            monitor0.notify();
                        }
                    }
                    
                
                    if (killer){
                        break loop;
                    }
                }


            }
        }
        
        
        
        
        void computePSF(){
            

            dparam.psf_many.computePSF(idSlice,x[idStack][idSlice], y[idStack][idSlice], deltaZ[idSlice]+z[idStack][idSlice],zwat);
            
            
        }
        
        void computeLikelihood(){
            

            double [][] psf = dparam.psf_many.getPSFPointer(idSlice);
            likelihood=0;
            for (int i=0;i<dparam.param.sizeoutput;i++){
                for (int ii=0;ii<dparam.param.sizeoutput;ii++){
                    if (paramImage.scmos[idStack]!=null){
                        mod[i][ii]=psf[i][ii]*(paramImage.A[idStack][idSlice]+fit_a[idStack]+fit_a_each[idStack][idSlice])+background[idStack][i][ii]+paramImage.scmos[idStack][i][ii];

                    }
                    else{
                        mod[i][ii]=psf[i][ii]*(paramImage.A[idStack][idSlice]+fit_a[idStack]+fit_a_each[idStack][idSlice])+background[idStack][i][ii];
                    }
                    likelihood+=mod[i][ii]-image[idStack][idSlice][i][ii]*Math.log(mod[i][ii]);
                    
                    
                }
            }
            
        }
            
            
        
        
    }


    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
}