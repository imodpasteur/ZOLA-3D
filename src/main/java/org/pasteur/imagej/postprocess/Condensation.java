/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;

import org.pasteur.imagej.data.*;
import ij.IJ;
import ij.gui.Plot;
import java.util.ArrayList;
import ij.measure.CurveFitter;
import java.awt.Color;
/**
 *
 * @author benoit
 */
public class Condensation {
    
    StackLocalization sl1;
    
    
    
    double maxDistanceMergingXY;
    double maxDistanceMergingZ;
    
    int [][] idLoc;
    
    
    double sizePix=100;//nm
    
    double maxX=Double.NEGATIVE_INFINITY;
    double minX=Double.MAX_VALUE;
    double maxY=Double.NEGATIVE_INFINITY;
    double minY=Double.MAX_VALUE;
       
    int offFrame;
    
    int bin=100;
    double [] histaxis = new double[bin];
    double [] hist = new double[bin];
    
    public Condensation(StackLocalization sl1,double maxDistanceMergingXY,double maxDistanceMergingZ,int offFrame){
        
        this.sl1=sl1;
        
        this.offFrame=offFrame;
        
        this.maxDistanceMergingXY=maxDistanceMergingXY;
        this.maxDistanceMergingZ=maxDistanceMergingZ;
        
        
        for (int i=0;i<bin;i++){
            histaxis[i]=i+1;
        }
        
        
    }
    
    
    public StackLocalization run(){
        ChainList [][] chain = ChainList.getChainList(sl1,maxDistanceMergingXY,maxDistanceMergingZ,offFrame);
        StackLocalization slnew=merge(chain);
        return slnew;
    }
    
    
    
    
    
    
    
    public StackLocalization merge(ChainList [][] chain){
        int condenseNumber=0;
        int totNumber=0;
        StackLocalization slnew=new StackLocalization();
        double mean=0;
        double countmean=0;
        
        
        //forward merging
        for (int i=0;i<sl1.fl.size();i++){
            FrameLocalization fl=new FrameLocalization(sl1.fl.get(i).numFrame);
            totNumber+=sl1.fl.get(i).loc.size();
            for (int ii=0;ii<sl1.fl.get(i).loc.size();ii++){
                if (!chain[i][ii].foundPrevious){//dont merge because if previous true...already merged (forward merging)
                    
                    if (!chain[i][ii].foundNext){//no merging
                        fl.loc.add(sl1.fl.get(i).loc.get(ii).copy());
                        hist[0]++;
                        mean+=1;
                        countmean++;
                    }
                    
                    else{//merge
                        PLocalization p_current=sl1.fl.get(i).loc.get(ii).copy();
                        //IJ.log("");
                        //IJ.log("idF "+sl1.fl.get(i).numFrame+"  "+0+"   "+sl1.fl.get(i).loc.get(ii).X);
                        double number=0;
                        boolean is=chain[i][ii].foundNext;
                        int idF=chain[i][ii].frameNext;
                        int idL=chain[i][ii].locNext;
                        ArrayList<PLocalization> p_next = new ArrayList<PLocalization>();
                        while (is){
                            p_next.add(sl1.fl.get(idF).loc.get(idL));
                            number++;
                            is=chain[idF][idL].foundNext;
                            if (is){
                                ChainList cl = chain[idF][idL];
                                idF=cl.frameNext;
                                idL=cl.locNext;
                            }
                        }
                        
                        
                        int ppphist=p_next.get(p_next.size()-1).frame-p_current.frame+1;
                        if (ppphist<bin){
                            hist[ppphist]++;
                        }
                        mean+=p_next.size();
                        countmean++;
                        double  x =0;
                        double  y =0;
                        double  z =0;
                        double [] cx = new double [p_next.size()+1];
                        double [] cy = new double [p_next.size()+1];
                        double [] cz = new double [p_next.size()+1];
                        condenseNumber+=p_next.size();
                        double sumcx=0;
                        double sumcy=0;
                        double sumcz=0;
                        boolean withCRLBX=true;
                        boolean withCRLBY=true;
                        boolean withCRLBZ=true;
                        boolean withPhoton=true;
                        boolean withBackground=true;
                        if (p_current.I<=0){
                            withPhoton=false;
                        }
                        if (p_current.BG<=0){
                            withBackground=false;
                        }
                        if (p_current.crlb_X<=0){
                            withCRLBX=false;
                        }
                        if (p_current.crlb_Y<=0){
                            withCRLBY=false;
                        }
                        if (p_current.crlb_Z<=0){
                            withCRLBZ=false;
                        }
                        for (int u=0;u<p_next.size();u++){
                            if (p_next.get(u).I<=0){
                                withPhoton=false;
                            }
                            if (p_next.get(u).BG<=0){
                                withBackground=false;
                            }
                            if (p_next.get(u).crlb_X<=0){
                                withCRLBX=false;
                            }
                            if (p_next.get(u).crlb_Y<=0){
                                withCRLBY=false;
                            }
                            if (p_next.get(u).crlb_Z<=0){
                                withCRLBZ=false;
                            }
                        }
                        if (withCRLBX){
                            cx[0]=1/p_current.crlb_X;
                            cx[0]=cx[0]*cx[0];
                            sumcx+=cx[0];
                            for (int u=0;u<p_next.size();u++){
                                cx[u+1]=1/p_next.get(u).crlb_X;
                                cx[u+1]=cx[u+1]*cx[u+1];
                                sumcx+=cx[u+1];
                            }
                        }
                        else if (withPhoton){
                            cx[0]=p_current.I;
                            sumcx+=cx[0];
                            for (int u=0;u<p_next.size();u++){
                                cx[u+1]=p_next.get(u).I;
                                sumcx+=cx[u+1];
                            }
                        }
                        else{
                            cx[0]=1./(double)(p_next.size()+1);
                            sumcx+=cx[0];
                            for (int u=0;u<p_next.size();u++){
                                cx[u+1]=1./(double)(p_next.size()+1);
                                sumcx+=cx[u+1];
                            }
                        }
                        
                        if (withCRLBY){
                            cy[0]=1/p_current.crlb_Y;
                            cy[0]=cy[0]*cy[0];
                            sumcy+=cy[0];
                            for (int u=0;u<p_next.size();u++){
                                cy[u+1]=1/p_next.get(u).crlb_Y;
                                cy[u+1]=cy[u+1]*cy[u+1];
                                sumcy+=cy[u+1];
                            }
                        }
                        else if (withPhoton){
                            cy[0]=p_current.I;
                            sumcy+=cy[0];
                            for (int u=0;u<p_next.size();u++){
                                cy[u+1]=p_next.get(u).I;
                                sumcy+=cy[u+1];
                            }
                        }
                        else{
                            cy[0]=1./(double)(p_next.size()+1);
                            sumcy+=cy[0];
                            for (int u=0;u<p_next.size();u++){
                                cy[u+1]=1./(double)(p_next.size()+1);
                                sumcy+=cy[u+1];
                            }
                        }
                        
                        
                        if (withCRLBZ){
                            cz[0]=1/p_current.crlb_Z;
                            cz[0]=cz[0]*cz[0];
                            sumcz+=cz[0];
                            for (int u=0;u<p_next.size();u++){
                                cz[u+1]=1/p_next.get(u).crlb_Z;
                                cz[u+1]=cz[u+1]*cz[u+1];
                                sumcz+=cz[u+1];
                            }
                        }
                        else if (withPhoton){
                            cz[0]=p_current.I;
                            sumcz+=cz[0];
                            for (int u=0;u<p_next.size();u++){
                                cz[u+1]=p_next.get(u).I;
                                sumcz+=cz[u+1];
                            }
                        }
                        else{
                            cz[0]=1./(double)(p_next.size()+1);
                            sumcz+=cz[0];
                            for (int u=0;u<p_next.size();u++){
                                cz[u+1]=1./(double)(p_next.size()+1);
                                sumcz+=cz[u+1];
                            }
                        }
                        
                        x=p_current.X*cx[0];
                        y=p_current.Y*cy[0];
                        z=p_current.Z*cz[0];
                        for (int u=0;u<p_next.size();u++){
                            x+=p_next.get(u).X*cx[u+1];
                            y+=p_next.get(u).Y*cy[u+1];
                            z+=p_next.get(u).Z*cz[u+1];
                        }
                        x/=sumcx;
                        y/=sumcy;
                        z/=sumcz;
                        p_current.X=x;
                        p_current.Y=y;
                        p_current.Z=z;
                        if (withCRLBX){
                            p_current.crlb_X=Math.sqrt(1/sumcx);
                        }
                        if (withCRLBY){
                            p_current.crlb_Y=Math.sqrt(1/sumcy);
                        }
                        if (withCRLBZ){
                            p_current.crlb_Z=Math.sqrt(1/sumcz);
                        }
                        if (withPhoton){
                            for (int u=0;u<p_next.size();u++){
                                p_current.I+=p_next.get(u).I;
                            }
                        }
                        if (withBackground){
                            for (int u=0;u<p_next.size();u++){
                                p_current.BG+=p_next.get(u).BG;
                            }
                            p_current.BG/=(double)(p_next.size()+1);
                        }
                        for (int u=0;u<p_next.size();u++){
                            p_next.get(u).id=p_current.id;
                        }
                        p_current.occurrence=p_next.get(p_next.size()-1).frame-p_current.frame+1;
                        fl.loc.add(p_current);
                 
                    }
                }
            }
            slnew.fl.add(fl);
        }
        
        
        IJ.log(""+condenseNumber+" among "+totNumber+"localizations are merged ("+Math.round((double)condenseNumber*100./(double)totNumber)+"% condensed)");
        IJ.log("It remains "+(totNumber-condenseNumber)+" localizations after merging");
        
        mean/=countmean;
        

        
        //IJ.log("persistence mean "+mean);
        
        
    
        return slnew;
    }
    
    
    
    
    
    public double [] fitLogExponential1(){
        int posit=0;
        loop:for (int i=0;i<hist.length-1;i++){
            if (hist[i]<1){
                posit=i;
                break loop;
            }
        }
        int len=posit;
        double [] x= new double [len];
        double [] y= new double [len];
        double [] yplot= new double [len];
        for (int i=0;i<len;i++){
            x[i]=histaxis[i];
            y[i]=Math.log(hist[i]);
            yplot[i]=(hist[i]);
        }
        CurveFitter cf = new CurveFitter(x,y);
        
        cf.doFit(CurveFitter.STRAIGHT_LINE);
        double b=0;
        double a=0;
        double [] p = cf.getParams();
        if (p!=null){
            a=p[0];
            b=p[1];
            IJ.log("fitted curve:  y = "+Math.exp(a)+" exp(-t/ "+(-1/b)+" )");
        }
        else{
            IJ.log("curve impossible to fit with exponential");
            
        }
        
        int count=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            count++;
        }
        double [] xf= new double [count];
        double [] yf= new double [count];
        int id=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            xf[id]=v;
            yf[id]=Math.exp(p[0])*Math.exp(p[1]*v);
            id++;
        }
        plotHist(x, yplot,xf,yf,"histogram of condensation","number of frames","occurrence number");
        double [] v = new double[2];
        v[0]=a;
        v[1]=(-1/b);
        return v;
    }
    
    
    
    public double [] fitExponential1(int start,int end){
        int posit=hist.length;
        loop:for (int i=start;i<hist.length-1;i++){
            if (hist[i]<1){
                posit=i;
                break loop;
            }
        }
        posit=Math.min(end, posit);
        if (posit<3){
            IJ.log("ERROR: not enough points for fitting");
            return null;
        }
        int len=posit-start;
        double [] x= new double [len];
        double [] y= new double [len];
        double [] xfull= new double [posit];
        double [] yfull= new double [posit];
        
        for (int i=0;i<posit;i++){
            xfull[i]=histaxis[i];
            IJ.write(""+hist[i]);
            yfull[i]=Math.log(hist[i]);
        }
        for (int i=0;i<len;i++){
            x[i]=histaxis[i+start];
            y[i]=Math.log(hist[i+start]);
        }
        CurveFitter cf = new CurveFitter(x,y);
        
        cf.doFit(CurveFitter.STRAIGHT_LINE);
        double b=0;
        double a=0;
        double [] p = cf.getParams();
        if (p!=null){
            a=p[0];
            b=p[1];
            IJ.log("fitted curve:  y = log["+a+" exp(-t/ "+(-1/b)+" )]");
        }
        else{
            IJ.log("curve impossible to fit with exponential");
            
        }
        
        int count=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            count++;
        }
        double [] xf= new double [count];
        double [] yf= new double [count];
        int id=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            xf[id]=v;
            yf[id]=p[0]+p[1]*v;
            id++;
        }
        
        plotHist(xfull, yfull,xf,yf,"histogram of condensation","number of frames","occurrence number");
        double [] v = new double[2];
        v[0]=a;
        v[1]=(-1/b);
        return v;
    }
    
    
    
    public double fitExponential2(double tau,double amplitude){
        int posit=0;
        loop:for (int i=0;i<hist.length-1;i++){
            if (hist[i]<5){
                posit=i;
                break loop;
            }
        }
        int len=posit;
        double [] x= new double [len];
        double [] y= new double [len];
        for (int i=0;i<len;i++){
            x[i]=histaxis[i];
            y[i]=hist[i];
        }
        
        CurveFitter cf1 = new CurveFitter(x,y);
        
        cf1.doFit(CurveFitter.EXPONENTIAL);
        
        double [] p = cf1.getParams();
        double ainit=1;
        double binit=1;
        if (p!=null){
            ainit=p[0];
            binit=p[1];
        }
        else{
        }
        
        
        CurveFitter cf = new CurveFitter(x,y);
        String function="y = "+amplitude+"*exp(x*"+(-1/tau)+")+a*exp(x*b)";
        
        double [] xp = {Math.abs(ainit-amplitude),binit};
        IJ.log("function "+function);
        cf.doCustomFit(function,xp,false);
        double b=0;
        p = cf.getParams();
        if (p!=null){
            double a=p[0];
            b=p[1];
            IJ.log("fitted curve:  y = "+amplitude+" exp(-t/ "+tau+" ) + "+a+" exp(-t/ "+(-1/b)+" )");
            IJ.log("non specific attach time= "+tau);
            IJ.log("specific attach time= "+(-1/b));
        }
        else{
            IJ.log("curve impossible to fit with exponential");
        }
        
        int count=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            count++;
        }
        double [] xf= new double [count];
        double [] yf= new double [count];
        int id=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            xf[id]=v;
            yf[id]=amplitude*Math.exp(-v/tau)+p[0]*Math.exp(v*p[1]);
            id++;
        }
        plotHist(x, y,xf,yf,"histogram of condensation","number of frames","occurrence number");
        
        return (1/b);
    }
    
    
    
    
    
    
    public void fitExponential2(){
        int posit=0;
        loop:for (int i=0;i<hist.length-1;i++){
            if (hist[i]<5){
                posit=i;
                break loop;
            }
        }
        int len=posit;
        double [] x= new double [len];
        double [] y= new double [len];
        for (int i=0;i<len;i++){
            x[i]=histaxis[i];
            y[i]=hist[i];
        }
        
        CurveFitter cf1 = new CurveFitter(x,y);
        
        cf1.doFit(CurveFitter.EXPONENTIAL);
        
        double [] p = cf1.getParams();
        double ainit=1;
        double binit=1;
        if (p!=null){
            ainit=p[0];
            binit=p[1];
        }
        else{
        }
        
        CurveFitter cf = new CurveFitter(x,y);
        String function="y = a*exp(x*b)+c*exp(x*d)";
        
        double [] xp = {ainit/2,binit/2,ainit/2,binit*2};
        //IJ.log("function "+function);
        cf.doCustomFit(function,xp,false);
        double a=0;
        double b=0;
        double c=0;
        double d=0;
        p = cf.getParams();
        if (p!=null){
            a=p[0];
            b=p[1];
            c=p[2];
            d=p[3];
            IJ.log("fitted curve:  y = "+a+" exp(-t/ "+(-1/b)+" ) + "+c+" exp(-t/ "+(-1/d)+" )");
            IJ.log("attach time 1= "+(-1/b));
            IJ.log("attach time 2= "+(-1/d));
        }
        else{
            IJ.log("curve impossible to fit with exponential");
        }
        
        int count=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            count++;
        }
        double [] xf= new double [count];
        double [] yf= new double [count];
        int id=0;
        for (double v=(double)x[0];v<(double)x[len-1];v+=((double)x[len-1]-(double)x[0])/(len*5)){
            xf[id]=v;
            yf[id]=a*Math.exp(v*b)+c*Math.exp(v*d);
            id++;
        }
        plotHist(x, y,xf,yf,"histogram of condensation","number of frames","occurrence number");
        
    }
    
    
    
    void plotHist(double [] x, double [] y, double [] xfit, double [] yfit,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel);
        
        //Plot p = new Plot(""+title,xlabel,ylabel,x,yfit,3);
        
        //p.setFont(0, 18);
        
        double xmin=Double.POSITIVE_INFINITY;
        double xmax=0;
        double ymin=Double.POSITIVE_INFINITY;
        double ymax=0;
        for (int i=0;i<xfit.length;i++){
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
            
            
        }
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
        p.addPoints(xfit, yfit, Plot.LINE);
        
        p.setColor(Color.blue);
        p.addPoints(x, y, Plot.CROSS);
        
        
        p.show();
        
    }
    
    
        
    
    
}
