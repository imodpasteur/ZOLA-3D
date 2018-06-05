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
import org.pasteur.imagej.process.PhaseRetrievalParametersDouble;
import org.pasteur.imagej.utils.FileVectorLoader;
import ij.IJ;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;

/**
 *
 * @author benoit
 */
public class DataPhase {
    
    public PhaseRetrievalParametersDouble param=null;
    
    public ZernikePhase phaseZer=null;
    
    public PSFPhase psf=null;
    
    public PSFPhaseMany psf_many=null;
    
    public boolean loading;
    
    
    public DataPhase(int sizeFFT,int sizeoutput,int orderGaussianPupil,double xystep,double zstep,double wavelength,double noil,double na,double sigmaGaussianKernel,int zernikeCoefNumber){
        
        
        
        param=new PhaseRetrievalParametersDouble(sizeFFT,sizeoutput,orderGaussianPupil,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
        
        
        //IJ.log("sigmaGaussianKernel "+sigmaGaussianKernel+"  "+param.sigmaGaussianKernel);
        
        
        psf=new PSFPhase(param);
        
        phaseZer = new ZernikePhase(param,zernikeCoefNumber,0);
        
        
        
    }
    
    
    public DataPhase(PhaseRetrievalParametersDouble p,int zerNumber){
        
        
        
        param=new PhaseRetrievalParametersDouble(p);
        
        
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
        
        param=new PhaseRetrievalParametersDouble(sizeFFT,sizeoutput,dphase.param);
        
        psf=new PSFPhase(param);
        
        phaseZer = new ZernikePhase(param,dphase.phaseZer.coef);
        phaseZer.setA(dphase.phaseZer.getAPointer());
        
        
        
        phaseZer.computeCombination();
        
        psf.updatePhase(phaseZer.getPhasePointer());
        //psf.updatePupil(pupilGauss.getPupil());
       // psf.updateKz(kzZer.computeCombination());
    }
    
    
    //clone dparam
    public DataPhase(DataPhase dphase){
        
        
        
        
        param=new PhaseRetrievalParametersDouble(dphase.param);
        
        
        psf=new PSFPhase(param);
        
        phaseZer = new ZernikePhase(param,dphase.phaseZer.coef);
        phaseZer.setA(dphase.phaseZer.getAPointer());
        
        
        
        
        
        psf.updatePhase(phaseZer.getPhasePointer());
        
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
        phaseZer.computeCombination();
        psf_many.updatePhasePointer(phaseZer.getPhasePointer());
        
        
        
        
    }
    
    
    
    
    public boolean load(int sizeFFT,String path){
        
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
                
                param = new PhaseRetrievalParametersDouble(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,sigmaGaussianKernel);
                ligne=br.readLine();
            }
            else{
                param = new PhaseRetrievalParametersDouble(sizeFFT,sizeFFT,order,xystep,zstep,wavelength,noil,na,1,1);
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
            param.updateweightZ(Double.parseDouble(lin[1]));
            
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
            
            
//            ligne=br.readLine();
//            lin=ligne.split(regex);
//            int nbCoefb=lin.length-1;
//            int [] coefb = new int [nbCoefb] ;
//            for (int i=1;i<lin.length;i++){
//                coefb[i-1]=Integer.parseInt(lin[i]);
//            }
//            
//            this.pupilZer=new ZernikePhaseJCudaFastDouble(param,coefb);
//            double [][] b = new double [order][nbCoefb];
//            for (int or=0;or<order;or++){
//                ligne=br.readLine();
//                lin=ligne.split(regex);
//
//                for (int i=1;i<lin.length;i++){
//                    b[or][i-1]=Double.parseDouble(lin[i]);
//                }
//                pupilZer.setA(b);
//            }
            
            
//            double [] varX = new double [order+1];
//            ligne=br.readLine();
//            lin=ligne.split(regex);
//            for (int or=1;or<lin.length;or++){
//                IJ.log(""+lin[or]);
//                varX[or-1]=Double.parseDouble(lin[or]);
//            }
//            
//            
//            
//            ligne=br.readLine();
//            lin=ligne.split(regex);
//            IJ.log(""+lin[1]);
//            double powerFlatTop=Double.parseDouble(lin[1]);
//            
//            IJ.log("varX "+varX[0]);
//            IJ.log("powerFlatTop "+powerFlatTop);
//           /this.pupilGauss=new GaussianPupil(param,varX,powerFlatTop);
            
            //ligne=br.readLine();
            //lin=ligne.split(regex);
            //int nbCoefp=lin.length-1;
            //int [] coefp = new int [nbCoefp] ;
            //for (int i=1;i<lin.length;i++){
            //    coefp[i-1]=Integer.parseInt(lin[i]);
            //}
            //
            //ligne=br.readLine();
            //lin=ligne.split(regex);
            //this.kzZer=new ZernikePhaseJCudaFastDouble(param,coefp);
            //double [] p = new double [nbCoefp];
            //for (int i=1;i<lin.length;i++){
            //    p[i-1]=Double.parseDouble(lin[i]);
            //    
            //}
            //kzZer.setA(p);
            
            
            
            
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
    
    public void save(String path){
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
                    sortie.println("orderPupil,"+param.orderGaussianPupil);
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
                    sortie.println("Wz,"+param.getweightZ());
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
    
}
