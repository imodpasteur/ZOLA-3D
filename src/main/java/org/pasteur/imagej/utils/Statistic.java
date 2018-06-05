/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;
import org.pasteur.imagej.data.*;
/**
 *
 * @author benoit
 */
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
/**
 *
 * @author benoit
 */
public class Statistic {
    
    public static void computeHist(StackLocalization sl,int variable,int bin){
        
        double mini=Double.POSITIVE_INFINITY;
        double maxi=Double.NEGATIVE_INFINITY;
        
        PLocalization p=null;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                p=sl.fl.get(i).loc.get(j);
                double v=p.getValueVariable(variable);
                if (maxi<v){
                    maxi=v;
                }
                if (mini>v){
                    mini=v;
                }
            }
        }
        if (p==null){
            IJ.log("Impossible to compute histogram: no localization");
            return;
        }
        if (maxi==mini){
            IJ.log("Impossible to compute histogram: all values = "+mini);
            return;
        }
        String nameVariable=p.getLabel(variable);
        double step=(maxi-mini)/(double)bin;
        
        double [] axis=new double [bin];
        double [] hist =new double [bin];
        for (int i=0;i<bin;i++){
            hist[i]=0;
            axis[i]=mini+(double)i*step;
        }
        
        
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                p=sl.fl.get(i).loc.get(j);
                double v=p.getValueVariable(variable);
                int posit=(int)((v-mini)/step);
                if ((posit>=0)&&(posit<bin)){
                    hist[posit]++;
                }
            }
        }
        
        plotHist(axis, hist,"histogram_"+nameVariable,nameVariable,"occurrence number");
        
    }
    
    
    
    
    
    
    public static void computeMean(StackLocalization sl){
        
        
        PLocalization p=null;
            
        looper:for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                p=sl.fl.get(i).loc.get(j);
                break looper;
            }
        }
        if (p!=null){
            
            int number=p.getNumberVariable();
            double [] mean = new double [number];
            double [] std = new double [number];
            String [] name = new String [number];
            for (int i=0;i<number;i++){
                mean[i]=0;
                std[i]=0;
                name[i]=p.getLabel(i);
            }
            double count=0;
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    p=sl.fl.get(i).loc.get(j);
                    for (int k=0;k<number;k++){
                        double v = p.getValueVariable(k);
                        mean[k]+=v;
                    }
                    count++;
                }
            }
            for (int i=0;i<number;i++){
                mean[i]/=count;
            }
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    p=sl.fl.get(i).loc.get(j);
                    for (int k=0;k<number;k++){
                        double v = p.getValueVariable(k);
                        std[k]+=(v-mean[k])*(v-mean[k]);
                    }
                }
            }
            
            
            
            IJ.log("statistics:");
            for (int i=0;i<number;i++){
                std[i]/=count;
                std[i]=Math.sqrt(std[i]);
                if (mean[i]>=0){
                    IJ.log("Mean("+name[i]+") = "+mean[i]+" +- "+std[i]);
                }
            }
            
        }
        
        
        
    }
    
    
    
    
    
    
     static void plotHist(double [] x, double [] y,String title,String xlabel,String ylabel){
        
         Plot p = new Plot(""+title,xlabel,ylabel,x,y,3);
        
        //Plot p = new Plot(""+title,xlabel,ylabel,x,yfit,3);
        
        //p.setFont(0, 18);
        
        double xmin=Double.POSITIVE_INFINITY;
        double xmax=0;
        double ymin=Double.POSITIVE_INFINITY;
        double ymax=0;
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
        
        
        
        p.show();
        
    }
    
    
    
}

