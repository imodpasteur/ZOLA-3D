/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

import ij.IJ;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import java.awt.Color;

/**
 *
 * @author benoit
 */
public class Regression {
    
    //compute linear regression from 2 vectors (y=r[0]x+r[1])
    public static double [] linear(double [] x, double [] y,boolean show){
        double [] r= new double [2];
        
        CurveFitter cf = new CurveFitter(x,y);
        
        cf.doFit(CurveFitter.STRAIGHT_LINE);
        double b=0;
        double a=0;
        double [] p = cf.getParams();
        if (p!=null){
            r[1]=p[0];
            r[0]=p[1];
            IJ.log("regression result "+r[0]+"  "+r[1]);
        }
        else{
            IJ.log("Sorry, curve impossible to fit using CurveFitter tool");
            
        }
        if (show){
            double [] xxx = new double [500];
            double [] yyy = new double [500];
            double [] xx = new double [10];
            double [] yy = new double [10];
            double minx=Double.POSITIVE_INFINITY;
            double maxx=Double.NEGATIVE_INFINITY;
            int maxPosit=0;
            int minPosit=0;
            for (int i=0;i<x.length;i++){
                if (minx>x[i]){
                    minx=x[i];
                    minPosit=i;
                }
                if (maxx<x[i]){
                    maxx=x[i];
                    maxPosit=i;
                }
                
            }
            
            for (int i=0;i<10;i++){
                xx[i]=minx+(double)i*(maxx-minx)/10.;
                yy[i]=r[0]*xx[i]+r[1];
            }
            xxx[0]=x[minPosit];
            yyy[0]=y[minPosit];
            xxx[1]=x[maxPosit];
            yyy[1]=y[maxPosit];
            for (int i=2;i<500;i++){
                int rand=(int)(Math.random()*x.length);
                
                xxx[i]=x[rand];
                yyy[i]=y[rand];
            }
            Regression.plot(xxx,yyy,xx,yy,"linear regression result (500 samples plotted)","X","Y");
        }
        return r;
    }
    
    
    
    
    
    public static void plot(double [] x, double [] y,double [] xline,double [] yline,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,xline,yline);
        
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
        for (int i=0;i<yline.length;i++){
            if (xline[i]<xMin){
                xMin=xline[i];
            }
            if (xline[i]>xMax){
                xMax=xline[i];
            }
            if (yline[i]<yMin){
                yMin=yline[i];
            }
            if (yline[i]>yMax){
                yMax=yline[i];
            }
        }
        p.setLimits(xMin-200, xMax+200, yMin-200, yMax+200);
        p.setLineWidth(2);
        p.setFont(0, 18);
        for (int ii=0;ii<x.length;ii++){
            
            p.setColor(Color.blue);
            p.add("CIRCLE",  x,y);
            p.setColor(Color.red);
            //p.show();


            
        }
        p.show();
        
    }
}
