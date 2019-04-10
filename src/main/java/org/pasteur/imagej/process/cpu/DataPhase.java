/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;

/**
 *
 * @author benoit
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import org.pasteur.imagej.process.PhaseParameters;
import org.pasteur.imagej.utils.*;
import ij.IJ;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;


import ij.IJ;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;

/**
 *
 * @author benoit
 */
public class DataPhase {
    
    public PhaseParameters param=null;
    
    public ZernikePhase phaseZer=null;
    
    public Phase phaseNonZer=null;
    
    public PSFPhase psf=null;
    
    public PSFPhaseMany psf_many=null;
    
    public boolean loading;
    
    
    public DataPhase(int sizeFFT,int sizeoutput,int orderGaussianPupil,double xystep,double zstep,double wavelength,double noil,double na,double sigmaGaussianKernel,int zernikeCoefNumber,boolean withApoFactor){
        
        
        
        param=new PhaseParameters(sizeFFT,sizeoutput,orderGaussianPupil,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel,withApoFactor);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSFPhase(param);
        
        phaseZer = new ZernikePhase(param,zernikeCoefNumber,0);
        
        //phaseNonZer = new Phase(param);
        
        
        
    }
    
    
    
    public DataPhase(int sizeFFT,int sizeoutput,int orderGaussianPupil,double xystep,double zstep,double wavelength,double noil,double na,double sigmaGaussianKernel,boolean withApoFactor){
        
        
        
        param=new PhaseParameters(sizeFFT,sizeoutput,orderGaussianPupil,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel,withApoFactor);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSFPhase(param);
        
        phaseZer = new ZernikePhase(param,3,0);
        
        phaseNonZer = new Phase(param);
        
        
        
    }
    
    
    public DataPhase(PhaseParameters p,int zerNumber){
        
        
        
        param=new PhaseParameters(p);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSFPhase(param);
        
        phaseZer = new ZernikePhase(param,zerNumber,0);
        
        
        
        
    }
    
    
    
    public DataPhase(int sizeFFT){
        
        String path=FileVectorLoader.chooseLoadingPath("Please select the calibration file (retrieved phase)");
        loading=this.load(sizeFFT,path);
        
        
    }
    
    
    
    public DataPhase(int sizeFFT,String path){
        if (path!=null){
            loading=load(sizeFFT,path);
        }
        
    }
    
    
    public DataPhase(){
        
        
    }
    
    
    
    //clone dparam but sizeFFT change
    public DataPhase(int sizeFFT,int sizeoutput,DataPhase dphase){
        
        param=new PhaseParameters(sizeFFT,sizeoutput,dphase.param);
        
        psf=new PSFPhase(param);
        
        
        if (dphase.param.zernikedPSF){
            phaseZer = new ZernikePhase(param,dphase.phaseZer.coef);
            phaseZer.setA(dphase.phaseZer.getAPointer());
            phaseZer.computeCombination();
            psf.updatePhase(phaseZer.getPhasePointer());
        }
        else{
            phaseNonZer  = new Phase(param);
            phaseNonZer.setValuesPhase(dphase.phaseNonZer.phase);
            psf.updatePhase(phaseNonZer.getPhasePointer());
        }
        
        
        
        
        //psf.updatePupil(pupilGauss.getPupil());
       // psf.updateKz(kzZer.computeCombination());
    }
    
    
    //clone dparam
    public DataPhase(DataPhase dphase){
        
        
        
        
        param=new PhaseParameters(dphase.param);
        
        
        psf=new PSFPhase(param);
        
        if (dphase.param.zernikedPSF){
            phaseZer = new ZernikePhase(param,dphase.phaseZer.coef);
            phaseZer.setA(dphase.phaseZer.getAPointer());
            phaseZer.computeCombination();
            psf.updatePhase(phaseZer.getPhasePointer());
        }
        else{
            phaseNonZer  = new Phase(param);
            phaseNonZer.setValuesPhase(dphase.phaseNonZer.phase);
            psf.updatePhase(phaseNonZer.getPhasePointer());
        }
        
    }
    
    
    
    
    
    
    
    public void setSizeoutput(int sizeoutput){
        
        psf.setSizeoutput(sizeoutput);
        
        
    }
    
    public void setNwat(double nwat){
        param.nwat=nwat;
        psf.resetKz();
        
        
    }
    
    
    
    
    public void setMany(int numberPSF){
        //psfMany=new PSFphaseJCudaFastDoubleMany(param,numberPSF);
        //psfMany.updatePhase(phaseZer.computeCombination());
        psf_many=new PSFPhaseMany(param,numberPSF);
        
        if (param.zernikedPSF){
            phaseZer.computeCombination();
            psf_many.updatePhasePointer(phaseZer.getPhasePointer());
        }
        else{
            psf_many.updatePhasePointer(phaseNonZer.getPhasePointer());
        }
        
        
        
        
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
        
        int nn=(path.lastIndexOf("."));
        if (nn>path.length()-5){
            path=path.substring(0,nn);
        }
        
       
        
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
            InputStream ips=new FileInputStream(path+".csv"); 
            
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
                
                param = new PhaseParameters(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel,false);
                ligne=br.readLine();
            }
            else{
                param = new PhaseParameters(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,1,false);
            }
            param.pathcalib=path+".csv";
            //IJ.log("sigGauss "+param.sigmaGaussianKernel);
            
            
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
            
            this.phaseZer=new ZernikePhase(param,coefa);
            double [] a = new double [nbCoefa];
            //for (int or=0;or<dimph;or++){
                ligne=br.readLine();
                lin=ligne.split(regex);
                for (int i=1;i<lin.length;i++){
                    a[i-1]=Double.parseDouble(lin[i]);
                }
            //}
            phaseZer.setA(a);
            
            
            
            
            psf=new PSFPhase(param);
            
            phaseZer.computeCombination();
            
            psf.updatePhase(phaseZer.getPhasePointer());
            
            
            
            
            //IJ.log("get pup...");
            
            //psf.updatePupil(pupilGauss.getPupil());
            //IJ.log("get pup ok");
            //ImageShow.imshow(psf.getPhase(),"phase");
            
            //ImageShow.imshow(psf.getPupil(),"pupil");
            //psf.updatePupil(pupilZer.computeCombination());
            //psf.updateKz(kzZer.computeCombination());
            
            br.close();
            //ImageShow.imshow(psf.getPhase(),"ph");
            //ImageShow.imshow(psf.getPupil(),"pup");
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
            int nn=(path.lastIndexOf("."));
            if (nn>path.length()-5){
                path=path.substring(0,nn);
            }

            try {


                PrintWriter sortie;

                //IJ.log("path : "+path+".csv");
                    sortie = new PrintWriter(new FileWriter(path+".csv", false));
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
            boolean withApoFactor;
            JSONformat json = new JSONformat();
            
            json.load(path);
            
            //Object obj = parser.parse(new FileReader(path));

            //JSONObject jsonObject = (JSONObject) obj;
            
            zernikedPSF = Boolean.parseBoolean(json.get("zernikedPSF").toString());
            
            zstep = Double.parseDouble(json.get("zstep").toString());
            
            
            size = Integer.parseInt(json.get("sizeFFT").toString());
            
            
            xystep = Double.parseDouble(json.get("sizePix").toString());
            
            
            na = Double.parseDouble(json.get("NA").toString());
            
            
            noil = Double.parseDouble(json.get("index").toString());
            
            
            wavelength = Double.parseDouble(json.get("wavelength").toString());
            
            
            sigmaGaussianKernel = Double.parseDouble(json.get("sigmaGaussianKernel").toString());
            
            
            String tmp=json.get("withApoFactor");
            if (tmp!=null){
                withApoFactor = Boolean.parseBoolean(tmp);
            }
            else{
                withApoFactor = false;
            }
            
            
            
            if (zernikedPSF){
                param = new PhaseParameters(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel,withApoFactor);
            }
            else{
                param = new PhaseParameters(size,size,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel,withApoFactor);
            }
            
            
            param.zernikedPSF=zernikedPSF;
            
            
            
            
            psf=new PSFPhase(param);
            
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
                this.phaseZer=new ZernikePhase(param,coefa);
                phaseZer.setA(a);
                

                phaseZer.setA(a);
                phaseZer.computeCombination();
                psf.updatePhase(phaseZer.getPhasePointer());
            }
            else{
                //JSONArray msg = (JSONArray) jsonObject.get("pixelPhaseValue");

                String [] msg = json.getVect("pixelPhaseValue");

                int nbCoefa=msg.length;
                
                double [] a = new double [nbCoefa];
                for (int i=0;i<msg.length;i++){
                    a[i]=Double.parseDouble((msg[i]));
                    //IJ.log("load "+i+"  "+a[i]);
                }
                
                this.phaseNonZer=new Phase(param);
                IJ.log("nbCoef "+nbCoefa+"    "+sizeFFT+"  "+param.size+"  "+param.size_cpu);
                IJ.log("nbCoefphaseNonZer "+phaseNonZer.phase.length+"  ");
                phaseNonZer.setValuesPhase(a);
                
                psf.updatePhase(phaseNonZer.getPhasePointer());
            }
            
            
            
            
            
        return true;
    }
    
    
    
    
    public void saveJSON(String path){
        if (path!=null){
            
            
            JSONformat obj = new JSONformat();
            //JSONObject obj = new JSONObject();
            obj.put("zstep",param.zstep);
            //obj.put("orderPupil",param.orderGaussianPupil);
            obj.put("sizeFFT",param.size);
            obj.put("sizePix",param.xystep);
            obj.put("NA",param.na);
            obj.put("index",param.noil);
            obj.put("wavelength",param.wavelength);
            obj.put("sigmaGaussianKernel",param.sigmaGaussianKernel);
            obj.put("zernikedPSF",param.zernikedPSF);
            obj.put("withApoFactor",param.withApoFactor);
            
            if (param.zernikedPSF){
                //JSONArray zernikelist = new JSONArray();
                double [] zer = new double[phaseZer.numCoef];
                for (int i=0;i<phaseZer.numCoef;i++){
                    //zernikelist.add(phaseZer.getA(i));
                    zer[i]=phaseZer.getA(i);
                    //IJ.log("save "+i+"  "+phaseZer.getA(i));

                }
                obj.put("zernike", zer);
            }
            else{
                //JSONArray list = new JSONArray();
                double [] list = new double[phaseNonZer.phase.length];
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
