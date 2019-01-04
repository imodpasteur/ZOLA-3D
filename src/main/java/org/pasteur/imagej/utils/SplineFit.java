/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

import ij.IJ;
import java.io.FileWriter;
import java.io.PrintWriter;
/**
 *
 * @author benoit
 */
public class SplineFit {
    
    public double [] x;
    public double [] y;
    
    public int nbFit;
    
    public PolynomialFit [] pf;
    
    public SplineFit(double [] y,double [] x){
        this.x=new double[x.length];
        this.y=new double[y.length];
        nbFit=Math.max(1,x.length+1-4);
        pf=new PolynomialFit[nbFit];
        boolean croissant=true;
        for (int i=0;i<x.length;i++){
            this.x[i]=x[i];
            this.y[i]=y[i];
        }
        for (int i=1;i<x.length;i++){
            if (this.x[i]<this.x[i-1]){
                croissant=false;
            }
        }
        if (!croissant){//tri bulle
            System.out.println("WARNING... spline fit... data should be sorted... slow sort of data");
            double tmp;
            for (int i=0;i<this.x.length;i++){
                for (int ii=1;ii<this.x.length-i;ii++){
                    if (this.x[ii]<this.x[ii-1]){
                        tmp=this.x[ii];
                        this.x[ii]=this.x[ii-1];
                        this.x[ii-1]=tmp;
                        tmp=this.y[ii];
                        this.y[ii]=this.y[ii-1];
                        this.y[ii-1]=tmp;
                    }
                }

            }
        }
        
    }
    
    
    //take only first dimension
    public SplineFit(double [][] y,double [][] x){
        if (x.length>1){
            IJ.log("WARNING... spline fit only for one dimension");
        }
        this.x=new double[x[0].length];
        this.y=new double[y[0].length];
        nbFit=Math.max(1,x[0].length+1-4);
        pf=new PolynomialFit[nbFit];
        
        boolean croissant=true;
        for (int i=0;i<x[0].length;i++){
            this.x[i]=x[0][i];
            this.y[i]=y[0][i];
        }
        for (int i=1;i<x[0].length;i++){
            if (this.x[i]<this.x[i-1]){
                croissant=false;
            }
        }
        if (!croissant){//tri bulle
            //IJ.log("WARNING... spline fit... data should be sorted... slow sort of data");
            double tmp;
            for (int i=0;i<this.x.length;i++){
                for (int ii=1;ii<this.x.length-i;ii++){
                    if (this.x[ii]<this.x[ii-1]){
                        tmp=this.x[ii];
                        this.x[ii]=this.x[ii-1];
                        this.x[ii-1]=tmp;
                        tmp=this.y[ii];
                        this.y[ii]=this.y[ii-1];
                        this.y[ii-1]=tmp;
                    }
                }

            }
        }
        
    }
    
    
    public SplineFit(String jsonContent){
    
        this.fromString(jsonContent);
        
    }
    
    public void run(){
        
        if (x.length>1){
            
            if (x.length==1){
                double [][] x0 = new double [1][1];//before first position
                double [][] y0 = new double [1][1];
                x0[0][0]=x[0];
                y0[0][0]=y[0];
                pf[0]=new PolynomialFit(0,y0,x0);
                pf[0].run();
            }
            else if (x.length==2){
                double [][] x0 = new double [1][2];//before first position
                double [][] y0 = new double [1][2];
                x0[0][0]=x[0];
                x0[0][1]=x[1];
                y0[0][0]=y[0];
                y0[0][1]=y[1];
                pf[0]=new PolynomialFit(1,y0,x0);
                pf[0].run();
            }
            else if (x.length==3){
                double [][] x0 = new double [1][3];//before first position
                double [][] y0 = new double [1][3];
                x0[0][0]=x[0];
                x0[0][1]=x[1];
                x0[0][2]=x[2];
                y0[0][0]=y[0];
                y0[0][1]=y[1];
                y0[0][2]=y[2];
                pf[0]=new PolynomialFit(2,y0,x0);
                pf[0].run();
            }
            else if (x.length>=4){
                double [][] x0 = new double [1][4];//middle points : cubic polynomial fit with 4 points
                double [][] y0 = new double [1][4];

                for (int id=0;id<nbFit;id++){
                    x0[0][0]=x[id];//between last and previous of last position
                    x0[0][1]=x[id+1];
                    x0[0][2]=x[id+2];
                    x0[0][3]=x[id+3];
                    y0[0][0]=y[id];
                    y0[0][1]=y[id+1];
                    y0[0][2]=y[id+2];
                    y0[0][3]=y[id+3];
                    //IJ.log("fit "+id+"  "+x[id]+"  "+y[id]+"     "+x[id+1]+"  "+y[id+1]+"     "+x[id+2]+"  "+y[id+2]+"     "+x[id+3]+"  "+y[id+3]);
                    pf[id]=new PolynomialFit(3,y0,x0);
                    
                    pf[id].run();
                }
            }
            
            
            
        }
        
        
        
    }
    
    
    
    
    public double  transform(double X){//transform one position (length is 1 only !!!)
        
        
        if (x.length>1){
            double [] r =new double [1];
            r[0]=X;
            
            if (x.length<4){
                r=pf[0].transform(r);
                return r[0];
            }
            else{
                int id=0;
        
                loopSearch:for (int i=0;i<x.length;i++){
                    if (X<x[i]){
                        id=i;
                        break loopSearch;
                    }
                }
                if(X>x[x.length-1]){
                    id=x.length;
                }



                if (id==1){
                    id=2;
                }

                if (id==x.length-1){
                    id=x.length-2;
                }

                if (id==0){
                    double [] rr =new double [1];
                    rr[0]=x[0];
                    rr=pf[id].transformPartial1(rr, 0);
                    r[0]=rr[0]*(X-x[0])+y[0];
                }
                else if (id==x.length){
                    double [] rr =new double [1];
                    rr[0]=x[x.length-1];
                    rr=pf[id-4].transformPartial1(rr, 0);
                    r[0]=rr[0]*(X-x[x.length-1])+y[x.length-1];
                }
                else{
                    
                    r=pf[id-2].transform(r);
                }

                return r[0];

            }
            
        }
        
        
        
        
        
        return -1;
    }
    
    
    
    public double [] transform(double [] X){//transform one position (length is 1 only !!!)
        if (X.length>1){
            IJ.log("WARNING... spline fit only for one dimension");
        }
        
        if (x.length>1){
            double [] r =new double [1];
            r[0]=X[0];
            
            if (x.length<4){
                r=pf[0].transform(r);
                return r;
            }
            else{
                int id=0;
        
                loopSearch:for (int i=0;i<x.length;i++){
                    id=i;
                    if (X[0]<x[i]){
                        id=i;
                        break loopSearch;
                    }
                }
                
                
                
                if(X[0]>x[x.length-1]){
                    id=x.length;
                    
                }



                if (id==1){
                    id=2;
                    
                }

                if (id==x.length-1){
                    id=x.length-2;
                    
                }

                if (id==0){
                    
                    double [] rr =new double [1];
                    rr[0]=x[0];
                    rr=pf[id].transformPartial1(rr, 0);
                    r[0]=rr[0]*(X[0]-x[0])+y[0];
                }
                else if (id==x.length){
                    double [] rr =new double [1];
                    rr[0]=x[x.length-1];
                    rr=pf[id-4].transformPartial1(rr, 0);
                    r[0]=rr[0]*(X[0]-x[x.length-1])+y[x.length-1];
                }
                else{
                    
                    r=pf[id-2].transform(r);
                }
                
                

                return r;

            }
            
        }
        
        
        
        
        
        return null;
    }
    
    
    
    
    
    
    
    
    
    
    
    public String toString(){
        
        JSONformat json = new JSONformat();
        
        json.put("x", x);
        
        json.put("y", x);
        
        String [] s = new String[pf.length];
        
        for (int i=0;i<pf.length;i++){
            s[i]=pf[i].toString();
        }
        json.put("pf", s);

        return json.toString();
        
            
        
    }
    
    public void fromString(String s){
        
        JSONformat json = new JSONformat();
        json.fromString(s);
        String [] sx = json.getVect("x");
        String [] sy = json.getVect("y");
        String [] spf = json.getVect("pf");
        
        this.nbFit=spf.length;
        this.pf = new PolynomialFit[nbFit];
        this.x = new double [sx.length];
        this.y = new double [sy.length];
        
        
        if (x.length!=y.length){
            IJ.log("error vector size splineFit class string loading");
        }
        
        for (int i=0;i<x.length;i++){
            this.x[i]=Double.parseDouble(sx[i]);
            this.y[i]=Double.parseDouble(sy[i]);
        }
        
        
        
        for (int i=0;i<nbFit;i++){
            pf[i]=new PolynomialFit();
            pf[i].fromString(spf[i]);
            
        }
        
        
        
        
        
    }
    
    
    
}
