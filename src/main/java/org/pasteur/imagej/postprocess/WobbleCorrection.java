/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;

import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import org.pasteur.imagej.data.StackLocalization;
import org.pasteur.imagej.utils.FileVectorLoader;
import org.pasteur.imagej.utils.PolynomialFit;
import org.pasteur.imagej.utils.SplineFit;
import org.pasteur.imagej.utils.JSONformat;

/**
 *
 * @author benoit
 */
public class WobbleCorrection {
    
    
    
    
    
    public static void calibration(String wobbleFilePath,StackLocalization sl,double z_step,double xy_step){
        
        
        
        
        
        
        
        
        
        int maxframe=0;
        Condensation cond=new Condensation(sl,xy_step*3, 10000000,10000000);//zstep=infinity ; off_frame=infinity
        StackLocalization smerged=cond.run();
        //after that, merged loc. share the same "id"
        
        for (int i=0;i<sl.fl.size();i++){
            for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                //IJ.log(""+sl.fl.get(i).loc.get(ii).id+"  "+sl.fl.get(i).loc.get(ii).frame+"  "+sl.fl.get(i).loc.get(ii).X+"  "+sl.fl.get(i).loc.get(ii).Y+"  "+sl.fl.get(i).loc.get(ii).Z);
                if (sl.fl.get(i).loc.get(ii).frame>maxframe){
                    maxframe=sl.fl.get(i).loc.get(ii).frame;
                }
            }
        }
        ArrayList<Integer> ids2keep = new ArrayList<Integer>();//keep all ids that have occ number > frame nb / 2 to remove noise
        for (int i=0;i<smerged.fl.size();i++){
            for (int ii=0;ii<smerged.fl.get(i).loc.size();ii++){
                if (smerged.fl.get(i).loc.get(ii).occurrence>maxframe/2){
                    ids2keep.add(smerged.fl.get(i).loc.get(ii).id);
                }
            }
        }
        
        
        
        int zeroFrameReference=maxframe/2;
        
        if (ids2keep.size()>0){
            
            //expected X=0
            //expected Y=0
            //expected Z=slope
            double [][][][] curve= new double [ids2keep.size()][3][2][];//4 dim: X,Y,Z
            double [][] frame = new double [ids2keep.size()][];
            int counttot=0;
            for (int c=0;c<ids2keep.size();c++){
                
                int count=0;
                for (int i=0;i<sl.fl.size();i++){
                    for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                        if (sl.fl.get(i).loc.get(ii).id==ids2keep.get(c)){
                            count++;
                        }
                    }
                }
                counttot+=count;
                curve[c][0][0]=new double[count];
                curve[c][1][0]=new double[count];
                curve[c][2][0]=new double[count];
                curve[c][0][1]=new double[count];
                curve[c][1][1]=new double[count];
                curve[c][2][1]=new double[count];
                frame[c]=new double[count];
                int posit=0;
                for (int i=0;i<sl.fl.size();i++){
                    for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                        if (sl.fl.get(i).loc.get(ii).id==ids2keep.get(c)){
                            curve[c][0][0][posit]=sl.fl.get(i).loc.get(ii).Z;
                            curve[c][0][1][posit]=sl.fl.get(i).loc.get(ii).X;
                            curve[c][1][0][posit]=sl.fl.get(i).loc.get(ii).Z;
                            curve[c][1][1][posit]=sl.fl.get(i).loc.get(ii).Y;
                            curve[c][2][0][posit]=sl.fl.get(i).loc.get(ii).Z;
                            curve[c][2][1][posit]=sl.fl.get(i).loc.get(ii).Z;
                            frame[c][posit]=sl.fl.get(i).loc.get(ii).frame;
                            posit++;
                        }
                    }
                }
                
                //X/Y centering (pre-registration)
                double mean=0;
                for (int i=0;i<curve[c][0][1].length;i++){
                    mean+=curve[c][0][1][i];
                }
                mean/=curve[c][0][1].length;
                for (int i=0;i<curve[c][0][1].length;i++){
                    curve[c][0][1][i]-=mean;
                }
                mean=0;
                for (int i=0;i<curve[c][1][1].length;i++){
                    mean+=curve[c][1][1][i];
                }
                mean/=curve[c][1][1].length;
                for (int i=0;i<curve[c][1][1].length;i++){
                    curve[c][1][1][i]-=mean;
                }
                
                //center Z
                for (int i=0;i<curve[c][2][1].length;i++){
                    curve[c][2][1][i]-=(frame[c][i]-zeroFrameReference)*z_step;//here we remve slope
                }
                mean=0;
                for (int i=0;i<curve[c][2][1].length;i++){
                    mean+=curve[c][2][1][i];
                }
                mean/=curve[c][2][1].length;
                for (int i=0;i<curve[c][2][1].length;i++){
                    curve[c][2][1][i]-=mean;
                }
                
                
            }
            
            
            
            
            //X
            
            double [] trueZ = new double [maxframe*2];
            for (int i=0;i<trueZ.length;i++){
                trueZ[i]=((double)i-(double)zeroFrameReference)*z_step;
            }
            //register all curve to curve 0
            
            for (int c=1;c<ids2keep.size();c++){
                
                MSE(curve[0][0][0],curve[0][0][1],curve[c][0][0],curve[c][0][1]);
                MSE(curve[0][1][0],curve[0][1][1],curve[c][1][0],curve[c][1][1]);
                MSE(curve[0][2][0],curve[0][2][1],curve[c][2][0],curve[c][2][1]);
                
            }
//            for (int c=0;c<ids2keep.size();c++){
//                
//                MSE(trueZ,trueZ,curve[c][2][0],curve[c][2][1]);
//            }
            //plot(curve[0][0][0],curve[0][0][1],curve[1][0][0],curve[1][0][1],"curveX","z","x");
            //plot(curve[0][1][0],curve[0][1][1],curve[1][1][0],curve[1][1][1],"curveY","z","y");
            //plot(trueZ,trueZ,curve[0][2][0],curve[0][2][1],"curveZt","z","z");
            
            double [][][] mergecurve= new double [3][2][counttot];//4 dim: X,Y,Z
            int nbpos=0;
            for (int c=0;c<ids2keep.size();c++){
                for (int i=0;i<curve[c][0][0].length;i++){
                    
                    mergecurve[0][0][nbpos]=curve[c][0][0][i];
                    mergecurve[0][1][nbpos]=curve[c][0][1][i];
                    mergecurve[1][0][nbpos]=curve[c][1][0][i];
                    mergecurve[1][1][nbpos]=curve[c][1][1][i];
                    mergecurve[2][0][nbpos]=curve[c][2][0][i];
                    mergecurve[2][1][nbpos]=curve[c][2][1][i];
                    nbpos++;
                }
            }
            //plot(curve[0][0][0],curve[0][0][1],mergecurve[0][0],mergecurve[0][1],"curveXmerge","z","x");
            //plot(curve[0][1][0],curve[0][1][1],mergecurve[1][0],mergecurve[1][1],"curveYmerge","z","x");
            //plot(curve[0][2][0],curve[0][2][1],mergecurve[2][0],mergecurve[2][1],"curveZmerge","z","x");
            
            
            
            
            
            
            
            
//            double [][] refX = new double [1][counttot];//x axis
//            double [][] refY = new double [1][counttot];//x axis
//            double [][] regX = new double [1][counttot];//x merge
//            double [][] regY = new double [1][counttot];//y merge
//            
//            double [][] refZ = new double [1][counttot];
//            double [][] regZ = new double [1][counttot];
//            
//            for (int i=0;i<counttot;i++){
//                refX[0][i]=mergecurve[0][0][i];
//                regX[0][i]=mergecurve[0][1][i];
//                refY[0][i]=mergecurve[1][0][i];
//                regY[0][i]=mergecurve[1][1][i];
//                refZ[0][i]=mergecurve[2][0][i];
//                regZ[0][i]=mergecurve[2][1][i];
//            }
//            
//            PolynomialFit pfX = new PolynomialFit(40,regX,refX);
//            PolynomialFit pfY = new PolynomialFit(40,regY,refY);
//            PolynomialFit pfZ = new PolynomialFit(6,regZ,refZ);
//            pfX.run();
//            pfY.run();
//            pfZ.run();
//            double [] xcurve = new double [trueZ.length];
//            double [] ycurve = new double [trueZ.length];
//            double [] zcurve = new double [trueZ.length];
//            double [] v = new double[1];
//            for (int i=0;i<trueZ.length;i++){
//                v[0]=trueZ[i];
//                v=pfX.transform(v);
//                xcurve[i]=v[0];
//                v[0]=trueZ[i];
//                v=pfY.transform(v);
//                ycurve[i]=v[0];
//                v[0]=trueZ[i];
//                v=pfZ.transform(v);
//                zcurve[i]=v[0];
//            }
//            plot(refX[0],regX[0],trueZ,xcurve,"polyfit","z","x");
//            plot(refY[0],regY[0],trueZ,ycurve,"polyfit","z","y");
//            plot(refZ[0],regZ[0],trueZ,zcurve,"polyfit","z","z");
            
            
                    
            int numberSplinePoints=counttot/10;
            
            
            double [][] refX = new double [1][numberSplinePoints];//x axis
            double [][] regX = new double [1][numberSplinePoints];//x merge
            smooth(mergecurve[0][0],mergecurve[0][1],refX[0],regX[0],numberSplinePoints);
            SplineFit sfX = new SplineFit(regX,refX);
            sfX.run();
            double [][] refY = new double [1][numberSplinePoints];//x axis
            double [][] regY = new double [1][numberSplinePoints];//x merge
            smooth(mergecurve[1][0],mergecurve[1][1],refY[0],regY[0],numberSplinePoints);
            SplineFit sfY = new SplineFit(regY,refY);
            sfY.run();
            double [][] refZ = new double [1][numberSplinePoints];//x axis
            double [][] regZ = new double [1][numberSplinePoints];//x merge
            smooth(mergecurve[2][0],mergecurve[2][1],refZ[0],regZ[0],numberSplinePoints);
            SplineFit sfZ = new SplineFit(regZ,refZ);
            sfZ.run();
            double [] xcurve = new double [trueZ.length];
            double [] ycurve = new double [trueZ.length];
            double [] zcurve = new double [trueZ.length];
            double [] v = new double[1];
            for (int i=0;i<trueZ.length;i++){
                v[0]=trueZ[i];
                v=sfX.transform(v);
                xcurve[i]=v[0];
                
                v[0]=trueZ[i];
                v=sfY.transform(v);
                ycurve[i]=v[0];
                
                v[0]=trueZ[i];
                v=sfZ.transform(v);
                zcurve[i]=v[0];
            }
            plot(mergecurve[0][0],mergecurve[0][1],trueZ,xcurve,"polyfitX","z","x bias");
            plot(mergecurve[1][0],mergecurve[1][1],trueZ,ycurve,"polyfitY","z","y bias");
            plot(mergecurve[2][0],mergecurve[2][1],trueZ,zcurve,"polyfitZ","z","z bias");
            
            
            save(wobbleFilePath,sfX,sfY,sfZ);
            
            /*SplineFit [] sfnew =load("/home/benoit/data/Challenge/testfit");
            
            
            
            for (int i=0;i<trueZ.length;i++){
                v[0]=trueZ[i];
                v=sfnew[0].transform(v);
                xcurve[i]=v[0];
                
                v[0]=trueZ[i];
                v=sfnew[1].transform(v);
                ycurve[i]=v[0];
                
                v[0]=trueZ[i];
                v=sfnew[2].transform(v);
                zcurve[i]=v[0];
            }
            plot(mergecurve[0][0],mergecurve[0][1],trueZ,xcurve,"polyfitXn","z","x");
            plot(mergecurve[1][0],mergecurve[1][1],trueZ,ycurve,"polyfitYn","z","x");
            plot(mergecurve[2][0],mergecurve[2][1],trueZ,zcurve,"polyfitZn","z","x");*/
            
            
        }
        else{
            IJ.log("ERROR: wobble calibration impossible because there is no bead detected over full axial range");
        }
        
        
        
    }
    
    
    
    
    
    public static StackLocalization correction(String pathcalibwobble,StackLocalization sl){
        
        
        SplineFit [] sf =load(pathcalibwobble);
            
        double [] v = new double [1];
        double [] v2 = new double [1];
        
        for (int i=0;i<sl.fl.size();i++){
            for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                v[0]=sl.fl.get(i).loc.get(ii).Z;
                
                v2=sf[0].transform(v);
                
                
                
                
                sl.fl.get(i).loc.get(ii).X-=v2[0];
                
                
                v2=sf[1].transform(v);
                sl.fl.get(i).loc.get(ii).Y-=v2[0];
                
                v2=sf[2].transform(v);
                sl.fl.get(i).loc.get(ii).Z-=v2[0];
            }
        }
        
        return sl;
    }
    
    private static void save(String path,SplineFit sfX,SplineFit sfY,SplineFit sfZ){
        
        
        JSONformat json = new JSONformat();
        json.put("sfX", sfX.toString());
        json.put("sfY", sfY.toString());
        json.put("sfZ", sfZ.toString());
        json.save(path);
        IJ.log("wobble calibration file saved");
    }
    
    
    private static SplineFit [] load(String path){
        
        JSONformat json = new JSONformat();
        json.load(path);
        SplineFit [] sf = new SplineFit[3];
        String strX = json.get("sfX");
        sf[0] = new SplineFit(strX);
        
        String strY = json.get("sfY");
        
        sf[1] = new SplineFit(strY);
        
        String strZ = json.get("sfZ");
        sf[2] = new SplineFit(strZ);
        
        
        
        
        return sf;

        
    }
    
    
    
    
    //this function provides sizeoutput averaging from input
    private static void smooth(double [] xinput,double [] yinput,double [] xoutput,double [] youtput,int sizeoutput){
        
        
        
        
        double [][] vect= new double[xinput.length][2];
        for (int i=0;i<xinput.length;i++){
            vect[i][0]=xinput[i];
            vect[i][1]=yinput[i];
        }
        
        Arrays.sort(vect, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o2[0]).compareTo(o1[0]);
            }
        });
        
        //if xinputs are stricly the same -> merge them (it prevents inversion matrix error in polynomial fitting)
        int count=1;
        for (int i=1;i<vect.length;i++){
            if (vect[i][0]!=vect[i-1][0]){
                count++;
            }
        }
        double [][] vectCopy = new double[count][2];
        int idcurrent=0;
        int idnew=0;
        for (int i=1;i<vect.length;i++){
            if (vect[i][0]!=vect[i-1][0]){
                if (idcurrent!=i-1){//mean from idcurrent -> i
                    double mean=0;
                    double num=0;
                    for (int j=idcurrent;j<i;j++){
                        mean+=vect[j][1];
                        num++;
                    }
                    mean/=num;
                    vectCopy[idnew][0]=vect[idcurrent][0];
                    vectCopy[idnew][1]=mean;
                    idnew++;
                }
                else{
                    vectCopy[idnew][0]=vect[idcurrent][0];
                    vectCopy[idnew][1]=vect[idcurrent][1];
                    idnew++;
                }
                idcurrent=i;
            }
            //treat last case
            if (i==vect.length){
                if (idcurrent!=i-1){//mean from idcurrent -> i+1
                    double mean=0;
                    double num=0;
                    for (int j=idcurrent;j<=i;j++){
                        mean+=vect[j][1];
                        num++;
                    }
                    mean/=num;
                    vectCopy[idnew][0]=vect[i][0];
                    vectCopy[idnew][1]=mean;
                    idnew++;
                }
                else{
                    vectCopy[idnew][0]=vect[i][0];
                    vectCopy[idnew][1]=vect[i][1];
                    idnew++;
                }
                idcurrent=i;
            }
            
            
        }
        
        
        int nb=(vectCopy.length/sizeoutput);
        double nbf=(double)vectCopy.length/(double)sizeoutput;
        for (int i=0; i<sizeoutput;i++){
            xoutput[i]=0;
            youtput[i]=0;
            double countn=0;
            for (int j=(int)((double)i*nbf-nbf/2.),k=0;k<nb*2;k++,j++){
                
                if ((j<vectCopy.length)&&(j>=0)){
                    xoutput[i]+=vectCopy[j][0];
                    youtput[i]+=vectCopy[j][1];
                    countn++;
                }
            }
            xoutput[i]/=countn;
            youtput[i]/=countn;
            
        }
        
        
    }
    
    private static void MSE(double [] zref,double [] ref,double [] zreg,double [] reg){
        
        
        double shift_parameter=0;
        double h=0.0001;
        
        
        
        double likelihood=Double.POSITIVE_INFINITY;
        
        iter:for (int it=0;it<1;it++){

            double lastLikelihood=likelihood;

            


            double lik1=error(zref,ref,zreg,reg,shift_parameter-h);
            double lik2=error(zref,ref,zreg,reg,shift_parameter);
            double lik3=error(zref,ref,zreg,reg,shift_parameter+h);

            double grad=0;
            if (Math.abs((lik3+lik1-2*lik2)/(h*h))==0){
                grad=((lik3-lik1)/(2*h));
            }
            else{
                grad=((lik3-lik1)/(2*h))/Math.abs((lik3+lik1-2*lik2)/(h*h));
            }
            

            boolean found=false;
            loop:for (double gamma=1;gamma>.001;gamma/=10){

                double lik=error(zref,ref,zreg,reg,shift_parameter-gamma*grad);



                if (lik<lik2){
                    shift_parameter-=gamma*grad;
                    likelihood=lik;
                    found=true;
                    break loop;
                }

            }
            
            
            //IJ.log("lik "+likelihood+"  "+shift_parameter);
            if (Math.abs(likelihood-lastLikelihood)<.0001){
                break iter;
            }
            
            
        }
        
        for (int c=0;c<reg.length;c++){
            reg[c]-=shift_parameter;
        }
        
        
        
    }
    
    
    
    
    
    private static double error(double [] zref,double [] ref,double [] zreg,double [] reg,double shift_parameter){
        double error=0;
        for (int i=0;i<ref.length;i++){
            double [] r=interpolate(zref[i],zreg,reg);
            
            if (r[1]>0){
                error+=(ref[i]-r[0]+shift_parameter)*(ref[i]-r[0]+shift_parameter);
                
            }
        }
        
        return error;
    }
    
    
    
    
    private static double [] interpolate(double z_expected,double [] z,double [] reg){
        
        double [] res = new double [2];
        
        int closest=-1;
        boolean isabove=false;
        double d=Double.POSITIVE_INFINITY;
        for (int i=0;i<z.length;i++){
            if (Math.abs(z[i]-z_expected)<d){
                d=Math.abs(z[i]-z_expected);
                closest=i;
                if (z[i]>z_expected){
                    isabove=true;
                }
                else{
                    isabove=false;
                }
            }
        }
        
        int secondclosest=-1;
        d=Double.POSITIVE_INFINITY;
        for (int i=0;i<z.length;i++){
            if ((i!=closest) && (( (isabove) && (z[i]<z_expected))) || ( (!isabove) && (z[i]>z_expected))){
                
                if (Math.abs(z[i]-z_expected)<d){
                    d=Math.abs(z[i]-z_expected);
                    secondclosest=i;
                }
            }
        }
        
        if (secondclosest==-1){//second closest not found --> edges --> no interpolation
            res[0]=reg[closest];
            res[1]=0;
            return res;
        }
        
        res[0]=((reg[closest]-reg[secondclosest])/(z[closest]-z[secondclosest]))*(z_expected-z[secondclosest])+reg[secondclosest];
        res[1]=1;
        return res;
        
        
    }
    
    
    
    private static void plot(double [] x, double [] y, double [] xfit, double [] yfit,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel);
        
        //Plot p = new Plot(""+title,xlabel,ylabel,x,yfit,3);
        
        //p.setFont(0, 18);
        
        double xmin=Double.POSITIVE_INFINITY;
        double xmax=0;
        double ymin=Double.POSITIVE_INFINITY;
        double ymax=0;
        /*for (int i=0;i<xfit.length;i++){
            if (xfit[i]<xmin){
                xmin=xfit[i];
            }
            if (xfit[i]>xmax){
                xmax=xfit[i];
            }
            if (yfit[i]<ymin){
                ymin=yfit[i];
            }
            if (yfit[i]>ymax){
                ymax=yfit[i];
            }
            
            
        }*/
        
        for (int i=0;i<x.length;i++){
            
            
            if (x[i]<xmin){
                xmin=x[i];
            }
            if (x[i]>xmax){
                xmax=x[i];
            }
            if (y[i]<ymin){
                ymin=y[i];
            }
            if (y[i]>ymax){
                ymax=y[i];
            }
        }
        p.setLimits(xmin-1, xmax+1, ymin-(ymax*.1), ymax+(ymax*.1));
        
        p.setLineWidth(2);
        p.setFont(0, 18);
        
        p.setColor(Color.red);
        p.addPoints(xfit, yfit, Plot.CROSS);
        
        p.setColor(Color.blue);
        p.addPoints(x, y, Plot.CROSS);
        
        
        p.show();
        
    }
    
    
    
    
    
}
