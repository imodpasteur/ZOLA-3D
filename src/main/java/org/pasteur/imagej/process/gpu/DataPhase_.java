/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;


import  org.pasteur.imagej.process.*;
import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.utils.*;
import ij.IJ;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaError;



import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;


/**
 *
 * @author benoit
 */
public class DataPhase_ {
    
    
    public PhaseParameters param=null;
    public PSF_double_ psf=null;
    
    public PSFmany_float_ psf_fMany=null;
    
    
    public PSFmany_float_ psf_f_crlbMany=null;
    
    
    //public PSFphaseJCudaFastDoubleMany psfMany=null;
    public Modelmany_double_ modelMany=null;
    
    
    //public ZernikePhaseJCudaFastDouble kzZer;
    
    public ZernikePhase_ phaseZer=null;
    public GenericPhase_ phaseNonZer=null;
    
    
    public boolean loading;
    
//    public ZernikePhaseJCudaFastDouble pupilZer;
    //public GaussianPupil pupilGauss;
    
    
    public DataPhase_(int sizeFFT,int sizeoutput,int orderGaussianPupil,double xystep,double zstep,double wavelength,double noil,double na,double sigmaGaussianKernel,int zernikeCoefNumber){
        
        
        
        param=new PhaseParameters(sizeFFT,sizeoutput,orderGaussianPupil,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSF_double_(param);
        
        phaseZer = new ZernikePhase_(param,zernikeCoefNumber,0);
        
        
        //pupilGauss = new GaussianPupil(param);
        
        
//        pupilZer = new ZernikePhaseJCudaFastDouble(param,1,1);
//        double value=1./(double)Math.sqrt(pupilZer.nbDataPerImage);
//        pupilZer.setA(0,0, value);
        
        //int [] coefZ = new int [4];
        //coefZ[0]=0;coefZ[1]=1;coefZ[2]=2;coefZ[3]=4;
        //kzZer = new ZernikePhaseJCudaFastDouble(param,coefZ);
        
    }
    
    
    
    
    public DataPhase_(int sizeFFT,int sizeoutput,int orderGaussianPupil,double xystep,double zstep,double wavelength,double noil,double na,double sigmaGaussianKernel){
        
        
        
        param=new PhaseParameters(sizeFFT,sizeoutput,orderGaussianPupil,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSF_double_(param);
        
        phaseZer = new ZernikePhase_(param,3,0);//for piston(replaced by defocus) and tilt (useful for bead registration)
        
        phaseNonZer = new GenericPhase_(param);
        
        
        
    }
    
    
    
    public DataPhase_(int sizeFFT){
        
        String path=FileVectorLoader.chooseLoadingPath("Please select the calibration file (retrieved phase)");
        loading=this.load(sizeFFT,path);
        
        
    }
    
    
    
    public DataPhase_(int sizeFFT,String path){
        if (path!=null){
            loading=load(sizeFFT,path);
        }
        
    }
    
    
    public DataPhase_(){
        
        
    }
    
    
    
    //clone dparam but sizeFFT change
    public DataPhase_(int sizeFFT,int sizeoutput,DataPhase_ dphase){
        if (dphase.param.zernikedPSF){
            param=new PhaseParameters(sizeFFT,sizeoutput,dphase.param);
        }
        else{
            param=new PhaseParameters(dphase.param.size,sizeoutput,dphase.param);
        }
        
        psf=new PSF_double_(param);
        
        if (dphase.param.zernikedPSF){
            phaseZer = new ZernikePhase_(param,dphase.phaseZer.coef);
            phaseZer.setApointer(dphase.phaseZer.getApointer());
            psf.updatePhase(phaseZer.computeCombination());
        }
        else{
            phaseNonZer  = new GenericPhase_(param);
            phaseNonZer.setValuesPhase(dphase.phaseNonZer.phase);
            psf.updatePhase(phaseNonZer.getPointerPhase());
        }
        
        
        //pupilGauss = new GaussianPupil(param);
        
//        pupilZer = new ZernikePhaseJCudaFastDouble(param,dphase.pupilZer.coef);
//        pupilZer.setApointer(dphase.pupilZer.getApointer());
        
        
        //kzZer = new ZernikePhaseJCudaFastDouble(param,dphase.kzZer.coef);
        //kzZer.setApointer(dphase.kzZer.getApointer());
        
        
        
        
        //psf.updatePupil(pupilGauss.getPupil());
       // psf.updateKz(kzZer.computeCombination());
    }
    
    
    //clone dparam
    public DataPhase_(DataPhase_ dphase,int stream){
        
        
        
        //param=new PhaseParameters(dphase.param.size,dphase.param.sizeoutput,dphase.param.xystep,dphase.param.zstep,dphase.param.wavelength,dphase.param.noil,dphase.param.na);
        
        
        param=new PhaseParameters(dphase.param);
        this.param.stream=stream;
        
        psf=new PSF_double_(param);
        
//        phaseZer = new ZernikePhaseJCudaFastDouble(param,dphase.phaseZer.coef);
//        phaseZer.setApointer(dphase.phaseZer.getApointer());
        
        if (dphase.param.zernikedPSF){
            phaseZer = new ZernikePhase_(param,dphase.phaseZer.coef);
            phaseZer.setApointer(dphase.phaseZer.getApointer());
            psf.updatePhase(phaseZer.computeCombination());
        }
        else{
            phaseNonZer  = new GenericPhase_(param);
            phaseNonZer.setValuesPhase(dphase.phaseNonZer.phase);
            psf.updatePhase(phaseNonZer.getPointerPhase());
        }
        
        //pupilGauss = new GaussianPupil(param);
        
//        pupilZer = new ZernikePhaseJCudaFastDouble(param,dphase.pupilZer.coef);
//        pupilZer.setApointer(dphase.pupilZer.getApointer());
        
        
        //kzZer = new ZernikePhaseJCudaFastDouble(param,dphase.kzZer.coef);
        //kzZer.setApointer(dphase.kzZer.getApointer());
        
        
        
        
        //psf.updatePupil(pupilGauss.getPupil());
//        psf.updatePupil(pupilZer.computeCombination());
        //psf.updateKz(kzZer.computeCombination());
    }
    
    
    public DataPhase_(PhaseParameters p,int zerNumber){
        
        
        
        param=new PhaseParameters(p);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSF_double_(p);
        
        phaseZer = new ZernikePhase_(param,zerNumber,0);
        
        
        
        
    }
    
    
    //useful for single emitter fitting where numberPSF=numberModel
    public void setMany(int numberModel){
        setMany(numberModel,false);
    }
    public void setMany(int numberModel,boolean scmosCamera){
        //psfMany=new PSFphaseJCudaFastDoubleMany(param,numberPSF);
        //psfMany.updatePhase(phaseZer.computeCombination());
        
        psf_fMany=new PSFmany_float_(param,numberModel);
        psf_fMany.updatePhase(phaseZer.computeCombination());
        
        
        psf_f_crlbMany=new PSFmany_float_(param,11);//for CRLB computation: 5*5 matrix -> 11 elem to compute derivatives
        psf_f_crlbMany.updatePhase(phaseZer.computeCombination());
        
        
        
        modelMany=new Modelmany_double_(psf_fMany,psf_f_crlbMany,param,numberModel,scmosCamera);
        
    }
    
    
    
    
    //useful for multi emitter fitting where numberPSF>numberModel
    public void setMany(int numberModel,int numberPSFperModel){
        setMany(numberModel,numberPSFperModel,false);
    }
    public void setMany(int numberModel,int numberPSFperModel,boolean scmosCamera){
        //psfMany=new PSFphaseJCudaFastDoubleMany(param,numberPSF);
        //psfMany.updatePhase(phaseZer.computeCombination());
        
        psf_fMany=new PSFmany_float_(param,numberPSFperModel*numberModel);
        psf_fMany.updatePhase(phaseZer.computeCombination());
        
        
        psf_f_crlbMany=new PSFmany_float_(param,11);//for CRLB computation: 5*5 matrix -> 11 elem to compute derivatives
        psf_f_crlbMany.updatePhase(phaseZer.computeCombination());
        
        
        
        modelMany=new Modelmany_double_(psf_fMany,psf_f_crlbMany,param,numberModel,numberPSFperModel,scmosCamera);
        
    }
    
    
    
    public void setStream(int stream){
        this.param.stream=stream;
    }
    
    public int getStream(){
        return this.param.stream;
    }
    
    public void setSizeoutput(int sizeoutput){
        
        psf.setSizeoutput(sizeoutput);
        
        if (this.psf_fMany!=null){
            this.psf_fMany.setSizeoutput(sizeoutput);
        }
        
        if (this.psf_f_crlbMany!=null){
            this.psf_f_crlbMany.setSizeoutput(sizeoutput);
        }
    }
    
    public void setNwat(double nwat){
        param.nwat=nwat;
        psf.resetKz();
        
        if (this.psf_fMany!=null){
            this.psf_fMany.resetKz();
        }
        
        if (this.psf_f_crlbMany!=null){
            this.psf_f_crlbMany.resetKz();
        }
    }
    
    
    public void free(){
        
        
        if (this.modelMany!=null){
            this.modelMany.free();
            modelMany=null;
        }
        
        if (this.psf_fMany!=null){
            this.psf_fMany.free();
            psf_fMany=null;
            
        }
        
        
        
        
        if (this.psf_f_crlbMany!=null){
            this.psf_f_crlbMany.free();
            psf_f_crlbMany=null;
        }
        if (phaseZer!=null)
            this.phaseZer.free();
        if (phaseNonZer!=null)
            this.phaseNonZer.free();
        //this.pupilGauss.free();
        //this.pupilZer.free();
        
        //this.kzZer.free();
        this.psf.free();
        //if (this.psfMany!=null){
        //    this.psfMany.free();
        //}
        
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    public boolean load(int sizeFFT,String path){
        
        
        if (path.endsWith(".csv")){
            loadCSV(sizeFFT,path);
        }
        else if (path.endsWith(".json")){
            loadJSON(sizeFFT,path);
        }
        else{
            IJ.log("ERROR: file extension unknown");
            return false;
        }
        return true;
    }
       
    public boolean loadCSV(int sizeFFT,String path){
        
        
        double na=0;
        double noil=0;
        double wavelength=0;
        int size=0;
        double xystep=0;
        double zstep=0;
        int order=0;
        int dimph=1;
        
        String regex=",";
        
        int nbLine=0;
        try{
            InputStream ips=new FileInputStream(path); 
            InputStreamReader ipsr=new InputStreamReader(ips);
            BufferedReader br=new BufferedReader(ipsr);
            
            String ligne;
            String [] lin;
            ligne=br.readLine();

            ligne=br.readLine();//nb phase (should be 1 here)
            lin=ligne.split(regex);
            
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            zstep=Double.parseDouble(lin[1]);
            
            
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            dimph=Integer.parseInt(lin[1]);
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            order=Integer.parseInt(lin[1]);
            
            
            
            
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            size=Integer.parseInt(lin[1]);
            

            //IJ.log("lin "+lin[1]);
            ligne=br.readLine();
            lin=ligne.split(regex);
            xystep=Double.parseDouble(lin[1]);

            ligne=br.readLine();
            lin=ligne.split(regex);
            na=Double.parseDouble(lin[1]);

            ligne=br.readLine();
            lin=ligne.split(regex);
            noil=Double.parseDouble(lin[1]);
            /*if (lin.length>=3){
                nwat=Double.parseDouble(lin[2]);
            }
            else{
                nwat=noil;
            }
            if (lin.length>=4){
                zfocus=Double.parseDouble(lin[2]);
            }
            else{
                zfocus=0;
            }*/
            ligne=br.readLine();
            lin=ligne.split(regex);
            wavelength=Double.parseDouble(lin[1]);
            
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            //IJ.log("lin "+lin[0]+"  "+lin[1]);
            
            if (lin[0].startsWith("sigma")){
                double sigmaGaussianKernel=Double.parseDouble(lin[1]);
                
                param = new PhaseParameters(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
                ligne=br.readLine();
            }
            else{
                param = new PhaseParameters(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,1);
            }
            //IJ.log("sigGauss "+param.sigmaGaussianKernel);
            param.pathcalib=path;
            
            ligne=br.readLine();

            ligne=br.readLine();//position
            lin=ligne.split(regex);
            
            
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            
            param.A=Double.parseDouble(lin[1]);

            ligne=br.readLine();
            lin=ligne.split(regex);
            param.B=Double.parseDouble(lin[1]);

            ligne=br.readLine();
            lin=ligne.split(regex);
            //param.updateweightZ(Double.parseDouble(lin[1]));
            
            ligne=br.readLine();
            lin=ligne.split(regex);
            int nbCoefa=lin.length-1;
            int [] coefa = new int [nbCoefa] ;
            for (int i=1;i<lin.length;i++){
                coefa[i-1]=Integer.parseInt(lin[i]);
                
            }
            
            this.phaseZer=new ZernikePhase_(param,coefa);
            
            
            double [] a = new double [nbCoefa];
            //for (int or=0;or<dimph;or++){
                ligne=br.readLine();
                lin=ligne.split(regex);
                for (int i=1;i<lin.length;i++){
                    
                    a[i-1]=Double.parseDouble(lin[i]);
                }
            //}
                
            phaseZer.setA(a);
            
            
            
            
            
            
            psf=new PSF_double_(param);
            
            psf.updatePhase(phaseZer.computeCombination());
            
            
            
            
            br.close();
            
            
            IJ.log("Calibration file loaded");
            return true;
        }		
        catch (Exception e){
            IJ.log("UNABLE to load calibration file "+e);
                System.out.println(e.toString());
        }
        return false;
        
    }
    
    
    
    
    
    
    
    public boolean save(String path){
        
        
        if (path.endsWith(".csv")){
            saveCSV(path);
        }
        else if (path.endsWith(".json")){
            saveJSON(path);
        }
        else{
            IJ.log("ERROR: file extension unknown");
            return false;
        }
        return true;
    }
    
    
    
    public void saveCSV(String path){
        if (path!=null){
            

            try {


                PrintWriter sortie;
                
                IJ.log("path : "+path);
                    sortie = new PrintWriter(new FileWriter(path, false));
                    sortie.println("Parameters");
                    sortie.println("phaseNumber,"+1);
                    sortie.println("zstep,"+param.zstep);
                    sortie.println("dimPhase,"+0);
                    sortie.println("orderPupil,"+0);
                    sortie.println("size,"+param.size);
                    sortie.println("sizePix,"+param.xystep);
                    sortie.println("NA,"+param.na);
                    sortie.println("index,"+param.noil);
                    sortie.println("wavelength,"+param.wavelength);
                    sortie.println("sigmaGaussianKernel,"+param.sigmaGaussianKernel);
                    sortie.println("sizeRingInPixel,"+param.sizeRadiusRingPixel);
                    sortie.println("Result Phase Retrieval Zernike Decomposition");
                    sortie.println("position,"+0);
                    sortie.println("A,"+param.A);
                    sortie.println("B,"+param.B);
                    sortie.println("Wz,"+1);
                    
                    //param = new PhaseParameters(param.size,param.size,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
                
                    
                    String sca="ZernikeCoef";
                    
//                    String scb="ZernikeCoef";
//                    String scp="ZernikeCoef";
                    
                    
//                    String sp="parabola";
                    for (int i=0;i<phaseZer.numCoef;i++){
                        sca+=","+phaseZer.coef[i];
                    }
                    sortie.println(""+sca);
                    //for (int or=0;or<param.dimPhase;or++){
                        String sa="phase";
                        for (int i=0;i<phaseZer.numCoef;i++){
                            if (i<3){
                                sa+=","+0;
                            }
                            else{
                                sa+=","+phaseZer.getA(i);
                            }
                            
                        }
                        sortie.println(""+sa);
                    //}
                    
                    
//                    String sgX="GaussianPupVarX";
//                    for (int i=0;i<pupilGauss.polyNomialSigmaX.length;i++){
//                        sgX+=","+pupilGauss.polyNomialSigmaX[i];
//                    }
//                    sortie.println(""+sgX);
//                    /*String sgY="GaussianPupVarY";
//                    for (int i=0;i<pupilGauss.polyNomialSigmaY.length;i++){
//                        sgY+=","+pupilGauss.polyNomialSigmaY[i];
//                    }
//                    sortie.println(""+sgY);*/
//                    sortie.println("GaussianPupSuperPower,"+pupilGauss.powerFlatTop);
                    
                    
                    
//                    for (int i=0;i<pupilZer.numCoef;i++){
//                        scb+=","+pupilZer.coef[i];
//                    }
//                    sortie.println(""+scb);
//                    for (int or=0;or<param.order;or++){
//                        String sb="pupil("+or+")";
//                        for (int i=0;i<pupilZer.numCoef;i++){
//                            sb+=","+pupilZer.getA(or,i);
//                        }
//                        sortie.println(""+sb);
//                    }
                    
                    
                    
                    //for (int i=0;i<kzZer.numCoef;i++){
                    //    scp+=","+kzZer.coef[i];
                    //}
                    //for (int i=0;i<kzZer.numCoef;i++){
                    //    sp+=","+kzZer.getA(i);
                    //}
                    
                    
                    
                    
                    //sortie.println(""+scp);
                    //sortie.println(""+sp);
                    IJ.log("Calibration file saved");
                    sortie.close();
                    
                    
                    
                    
            } catch (Exception e) {
                    e.printStackTrace();
            }
            
            
            
            
            //psf.save(path);
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    public boolean loadJSON(int sizeFFT,String path){
        
        //JSONParser parser = new JSONParser();

        //try{
            double na=0;
            double noil=0;
            double wavelength=0;
            int size=0;
            double xystep=0;
            double zstep=0;
            int order=0;
            int dimph=1;
            double sigmaGaussianKernel=0;
            boolean zernikedPSF;
            //Object obj = parser.parse(new FileReader(path));


            JSONformat json= new JSONformat();

            json.load(path);
            
            //JSONObject jsonObject = (JSONObject) obj;

            zernikedPSF = Boolean.parseBoolean(json.get("zernikedPSF").toString());

            zstep = Double.parseDouble(json.get("zstep").toString());


            size = Integer.parseInt(json.get("sizeFFT").toString());


            xystep = Double.parseDouble(json.get("sizePix").toString());


            na = Double.parseDouble(json.get("NA").toString());


            noil = Double.parseDouble(json.get("index").toString());


            wavelength = Double.parseDouble(json.get("wavelength").toString());


            sigmaGaussianKernel = Double.parseDouble(json.get("sigmaGaussianKernel").toString());



            if (zernikedPSF){
                param = new PhaseParameters(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
            }
            else{
                param = new PhaseParameters(size,size,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
            }
            param.zernikedPSF=zernikedPSF;

            psf=new PSF_double_(param);

            if (zernikedPSF){
                // loop array
                //JSONArray msg = (JSONArray) jsonObject.get("zernike");

                String [] msg = json.getVect("zernike");

                int nbCoefa=msg.length;
                int [] coefa = new int [nbCoefa] ;
                double [] a = new double [nbCoefa];
                for (int i=0;i<msg.length;i++){
                    coefa[i]=i;
                    a[i]=Double.parseDouble((msg[i]));
                    //IJ.log("load "+i+"  "+a[i]);
                }

                this.phaseZer=new ZernikePhase_(param,coefa);

                phaseZer.setA(a);

                psf.updatePhase(phaseZer.computeCombination());
            }
            else{

                String [] msg = json.getVect("pixelPhaseValue");
                //JSONArray msg = (JSONArray) jsonObject.get("pixelPhaseValue");


                IJ.log("len "+msg.length+"  "+msg[0]+"  "+msg[msg.length-1]);
                int nbCoefa=msg.length;
                double [] a = new double [nbCoefa];
                for (int i=0;i<msg.length;i++){
                    a[i]=Double.parseDouble((msg[i]));
                    //IJ.log("load "+i+"  "+a[i]);
                }

                this.phaseNonZer=new GenericPhase_(param);
                phaseNonZer.setValuesPhase(a);

                psf.updatePhase(phaseNonZer.getPointerPhase());
            }
        /*}
        catch(Exception e){
            IJ.log("exception loading "+e);
        }*/
        
        
        return true;
    }
    
    
    
    
    public void saveJSON(String path){
        if (path!=null){
            
            
            
            //JSONObject obj = new JSONObject();
            
            JSONformat obj = new JSONformat();
            
            obj.put("zstep",param.zstep);
            //obj.put("orderPupil",param.orderGaussianPupil);
            obj.put("sizeFFT",param.size);
            obj.put("sizePix",param.xystep);
            obj.put("NA",param.na);
            obj.put("index",param.noil);
            obj.put("wavelength",param.wavelength);
            obj.put("sigmaGaussianKernel",param.sigmaGaussianKernel);
            obj.put("zernikedPSF",param.zernikedPSF);
            
            if (param.zernikedPSF){
                //JSONArray zernikelist = new JSONArray();
                double [] zer = new double [phaseZer.numCoef];
                for (int i=0;i<phaseZer.numCoef;i++){
                    //zernikelist.add(phaseZer.getA(i));
                    zer[i]=phaseZer.getA(i);
                    //IJ.log("save "+i+"  "+phaseZer.getA(i));

                }
                obj.put("zernike", zer);
            }
            else{
                //JSONArray list = new JSONArray();
                double [] list = new double [phaseNonZer.phase.length];
                for (int i=0;i<phaseNonZer.phase.length;i++){
                    //list.add(phaseNonZer.phase[i]);
                    list[i]=phaseNonZer.phase[i];

                }
                obj.put("pixelPhaseValue", list);
            }
            
            
            obj.save(path);
            
            
            
            
        }
    }
    
    
    
    
    
}