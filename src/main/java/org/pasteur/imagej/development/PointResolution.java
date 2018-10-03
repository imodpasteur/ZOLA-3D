/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.development;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.postprocess.ZRendering;
import ij.IJ;
import ij.gui.Plot;
import java.util.ArrayList;
import java.awt.Color;
import java.util.Arrays;
/**
 *
 * @author benoit
 */
public class PointResolution{

    double precisionSamplingPolyCurve=0.1;//(nm)

    public PointResolution(StackLocalization sl,double sizePix,double [] xpoints, double [] ypoints, double [] zpoints,double maximumDistance,double binSize,String pathRes,String pathResexclude){
        
        
        
        StackLocalization slnew=new StackLocalization();
        StackLocalization slnewexclude=new StackLocalization();
        
        int nbpoints=xpoints.length;
        IJ.log("nb pt "+nbpoints);
        double minx=Double.POSITIVE_INFINITY;
        double miny=Double.POSITIVE_INFINITY;
        double minz=Double.POSITIVE_INFINITY;
        double maxx=Double.NEGATIVE_INFINITY;
        double maxy=Double.NEGATIVE_INFINITY;
        double maxz=Double.NEGATIVE_INFINITY;
        
        for (int p=0;p<nbpoints;p++){
            xpoints[p]*=sizePix;
            ypoints[p]*=sizePix;
            zpoints[p]*=sizePix;
            
        }
        
        
        
        
        for (int p=0;p<nbpoints;p++){
            if (minx>xpoints[p]){
                minx=xpoints[p];
            }
            if (miny>ypoints[p]){
                miny=ypoints[p];
            }
            if (minz>zpoints[p]){
                minz=zpoints[p];
            }
            if (maxx<xpoints[p]){
                maxx=xpoints[p];
            }
            if (maxy<ypoints[p]){
                maxy=ypoints[p];
            }
            if (maxz<zpoints[p]){
                maxz=zpoints[p];
            }
        }
        minx-=maximumDistance;
        miny-=maximumDistance;
        minz-=maximumDistance;
        maxx+=maximumDistance;
        maxy+=maximumDistance;
        maxz+=maximumDistance;
        
        //IJ.log("min/max "+minx+"  "+maxx+"  "+miny+"  "+maxy+"  "+minz+"  "+maxz);
        
        int nbloc2=0;
        double x,y,z,d,dmin;
        double a,b,c;
        double minXlocalization=Double.POSITIVE_INFINITY;
        double minYlocalization=Double.POSITIVE_INFINITY;
        double minZlocalization=Double.POSITIVE_INFINITY;
        for (int f=0;f<sl.fl.size();f++){
            for (int p=0;p<sl.fl.get(f).loc.size();p++){
                x=sl.fl.get(f).loc.get(p).X;
                if (x<minXlocalization){
                    minXlocalization=x;
                }
                y=sl.fl.get(f).loc.get(p).Y;
                if (y<minYlocalization){
                    minYlocalization=y;
                }
                z=sl.fl.get(f).loc.get(p).Z;
                if (z<minZlocalization){
                    minZlocalization=z;
                }
            }
        }
        for (int f=0;f<sl.fl.size();f++){
            FrameLocalization flnew = new FrameLocalization(sl.fl.get(f).numFrame);
            FrameLocalization flnewexclude = new FrameLocalization(sl.fl.get(f).numFrame);
            for (int p=0;p<sl.fl.get(f).loc.size();p++){
                x=sl.fl.get(f).loc.get(p).X-minXlocalization;
                y=sl.fl.get(f).loc.get(p).Y-minYlocalization;
                z=sl.fl.get(f).loc.get(p).Z-minZlocalization;
                
                //IJ.log("posit "+x+"  "+y+"  "+z);
                boolean found=false;
                if ((x>=minx)&&(x<=maxx)&&(y>=miny)&&(y<=maxy)&&(z>=minz)&&(z<=maxz)){
                    
                    dmin=Double.POSITIVE_INFINITY;
                    int idpt=0;
                    for (int pt=0;pt<nbpoints;pt++){
                        a=(x-xpoints[pt])*(x-xpoints[pt]);
                        b=(y-ypoints[pt])*(y-ypoints[pt]);
                        c=(z-zpoints[pt])*(z-zpoints[pt]);
                        d=a+b+c;//it does not work when i put directly the result in d. i need to use a, b and c... surprising
                        if (dmin>d){
                            dmin=d;
                            idpt=pt;
                        }
                    }
                    dmin=Math.sqrt(dmin);
                    
                    
                    
                    if (dmin<=maximumDistance){
                   // if (true){
                        //localization is ok
                        PLocalization ploc=sl.fl.get(f).loc.get(p).copy();
                        //ploc.frame=0;//because frame number does not matter
                        flnew.loc.add(ploc);
                        found=true;
                        nbloc2++;
                    }
                }
                if (!found){
                    PLocalization ploc=sl.fl.get(f).loc.get(p).copy();
                    flnewexclude.loc.add(ploc);
                }
            }
            
            slnew.fl.add(flnew);
            slnewexclude.fl.add(flnewexclude);
        }
        
        
        
        if (nbloc2>0){
            ZRendering.hist3D(slnew, sizePix,0);
            
        }
        else{
            IJ.log("oops, no localization found in the region");
            return;
        }
        
        if ((pathRes!=null)&&(pathRes.length()>1)){
            slnew.save(pathRes);
            IJ.log("localization table with points saved");
        }
        else{
            IJ.log("localization table with points not saved");
        }
        
        if ((pathResexclude!=null)&&(pathResexclude.length()>1)){
            slnewexclude.save(pathResexclude);
            IJ.log("localization table without points saved");
        }
        else{
            IJ.log("localization table without points not saved");
        }
        
        IJ.log("nbpoints "+nbpoints);
        
        
        
        ArrayList<Double> [] dx=new ArrayList[nbpoints];
        ArrayList<Double> [] dy=new ArrayList[nbpoints];
        ArrayList<Double> [] dz=new ArrayList[nbpoints];
        ArrayList<Double> [] dphoton=new ArrayList[nbpoints];
        ArrayList<Double> [] dbckg=new ArrayList[nbpoints];
        
        for (int i=0;i<nbpoints;i++){
            dx[i]=new ArrayList<Double>();
            dy[i]=new ArrayList<Double>();
            dz[i]=new ArrayList<Double>();
            dphoton[i]=new ArrayList<Double>();
            dbckg[i]=new ArrayList<Double>();
        }
        
        int t=0;
        for (int f=0;f<sl.fl.size();f++){
            for (int p=0;p<sl.fl.get(f).loc.size();p++){
                
                x=sl.fl.get(f).loc.get(p).X-minXlocalization;
                y=sl.fl.get(f).loc.get(p).Y-minYlocalization;
                z=sl.fl.get(f).loc.get(p).Z-minZlocalization;
                
                //IJ.log("posit "+x+"  "+y+"  "+z);
                if ((x>=minx)&&(x<=maxx)&&(y>=miny)&&(y<=maxy)&&(z>=minz)&&(z<=maxz)){
                    
                    dmin=Double.POSITIVE_INFINITY;
                    int idpt=0;
                    for (int pt=0;pt<nbpoints;pt++){
                        a=(x-xpoints[pt])*(x-xpoints[pt]);
                        b=(y-ypoints[pt])*(y-ypoints[pt]);
                        c=(z-zpoints[pt])*(z-zpoints[pt]);
                        d=a+b+c;//it does not work when i put directly the result in d. i need to use a, b and c... surprising
                        if (dmin>d){
                            dmin=d;
                            idpt=pt;
                        }
                    }
                    dmin=Math.sqrt(dmin);
                    if (dmin<=maximumDistance){
                        dx[idpt].add(sl.fl.get(f).loc.get(p).X);
                        dy[idpt].add(sl.fl.get(f).loc.get(p).Y);
                        dz[idpt].add(sl.fl.get(f).loc.get(p).Z);
                        dphoton[idpt].add(sl.fl.get(f).loc.get(p).I);
                        dbckg[idpt].add(sl.fl.get(f).loc.get(p).BG);
                    }
                    
                }
                
                
                
                
                t++;
            }
        }
        
        
        
        double [] meanX = new double [nbpoints];
        double [] meanY = new double [nbpoints];
        double [] meanZ = new double [nbpoints];
        double [] meanPhoton = new double [nbpoints];
        double [] meanBCKG = new double [nbpoints];
        double meanI=0;
        double meanBG=0;
        
        for (int i=0;i<nbpoints;i++){
            meanX[i]=0;
            meanY[i]=0;
            meanZ[i]=0;
            meanPhoton[i]=0;
            meanBCKG[i]=0;
            for (int u=0;u<dx[i].size();u++){
                meanX[i]+=(double)dx[i].get(u);
                meanY[i]+=(double)dy[i].get(u);
                meanZ[i]+=(double)dz[i].get(u);
                meanPhoton[i]+=(double)dphoton[i].get(u);
                meanBCKG[i]+=(double)dbckg[i].get(u);
            }
            meanX[i]/=(double)dx[i].size();
            meanY[i]/=(double)dy[i].size();
            meanZ[i]/=(double)dz[i].size();
            meanPhoton[i]/=(double)dphoton[i].size();
            meanBCKG[i]/=(double)dbckg[i].size();
            meanI+=meanPhoton[i];
            meanBG+=meanBCKG[i];
            
        }
        meanI/=(double)nbpoints;
        meanBG/=(double)nbpoints;
        
        
        double [] stdX = new double [nbpoints];
        double [] stdY = new double [nbpoints];
        double [] stdXY = new double [nbpoints];
        double [] stdZ = new double [nbpoints];
        double [] std3D = new double [nbpoints];
        double meanSTDX=0;
        double meanSTDY=0;
        double meanSTDXY=0;
        double meanSTDZ=0;
        double meanSTD3d=0;
        
        for (int i=0;i<nbpoints;i++){
            stdX[i]=0;
            stdY[i]=0;
            stdZ[i]=0;
            std3D[i]=0;
            for (int u=0;u<dx[i].size();u++){
                stdX[i]+=((double)dx[i].get(u)-meanX[i])*((double)dx[i].get(u)-meanX[i]);
                stdY[i]+=((double)dy[i].get(u)-meanY[i])*((double)dy[i].get(u)-meanY[i]);
                stdZ[i]+=((double)dz[i].get(u)-meanZ[i])*((double)dz[i].get(u)-meanZ[i]);
                
                
            }
            stdX[i]/=(double)dx[i].size();
            stdY[i]/=(double)dy[i].size();
            stdZ[i]/=(double)dz[i].size();
            stdXY[i]=Math.sqrt(stdX[i]+stdY[i]);
            std3D[i]=Math.sqrt(stdX[i]+stdY[i]+stdZ[i]);
            
            
            meanSTD3d+=std3D[i];
            meanSTDXY+=stdXY[i];
            stdX[i]=Math.sqrt(stdX[i]);
            stdY[i]=Math.sqrt(stdY[i]);
            stdZ[i]=Math.sqrt(stdZ[i]);
            meanSTDX+=stdX[i];
            meanSTDY+=stdY[i];
            meanSTDZ+=stdZ[i];
        }
        
        meanSTDX/=(double)nbpoints;
        meanSTDY/=(double)nbpoints;
        meanSTDZ/=(double)nbpoints;
        meanSTD3d/=(double)nbpoints;
        meanSTDXY/=(double)nbpoints;
        
        IJ.log("mean photon "+meanI);
        IJ.log("mean background "+meanBG);
        IJ.log("mean X precision "+meanSTDX);
        IJ.log("mean Y precision "+meanSTDY);
        IJ.log("mean Lateral precision "+meanSTDXY);
        IJ.log("mean Z precision "+meanSTDZ);
        IJ.log("mean 3D precision "+meanSTD3d);
        IJ.log("");IJ.log("");IJ.log("");
        IJ.log("Each point study:");
        IJ.log("loc number , X precision , Y precision , Z precision , 3D precision , mean photon , mean background");
        
        for (int i=0;i<nbpoints;i++){
            IJ.log(""+dx[i].size()+" , "+stdX[i]+" , "+stdY[i]+" , "+stdZ[i]+" , "+std3D[i]+" , "+meanPhoton[i]+" , "+meanBCKG[i]);
        }
        
        IJ.log("");IJ.log("");IJ.log("");
        IJ.log("Each Localization study (for box plot):");
        IJ.log("idPoint , X precision , Y precision , Z precision , photon number, background number");
        
        int nbLoc=0;
        for (int i=0;i<nbpoints;i++){
            nbLoc+=dx[i].size();
        }
        double [] xarray = new double [nbLoc];
        double [] yarray = new double [nbLoc];
        double [] zarray = new double [nbLoc];
        double [] parray = new double [nbLoc];
        double [] barray = new double [nbLoc];
        for (int i=0,k=0;i<nbpoints;i++){
            for (int ii=0;ii<dx[i].size();ii++){
                xarray[k]=Math.abs(dx[i].get(ii)-meanX[i]);
                yarray[k]=Math.abs(dy[i].get(ii)-meanY[i]);
                zarray[k]=Math.abs(dz[i].get(ii)-meanZ[i]);
                parray[k]=Math.abs(dphoton[i].get(ii));
                barray[k]=Math.abs(dbckg[i].get(ii));
                IJ.log(""+i+" , "+xarray[k]+" , "+yarray[k]+" , "+zarray[k]+" , "+parray[k]+" , "+barray[k]);
                k++;
            }
            
        }
        Arrays.sort(xarray);
        Arrays.sort(yarray);
        Arrays.sort(zarray);
        Arrays.sort(parray);
        Arrays.sort(barray);
        IJ.log("");IJ.log("");
        
        IJ.log("stat : median, quart1, quart3, min, max");
        IJ.log("X stat: "+xarray[xarray.length/2]+", "+xarray[xarray.length/4]+", "+xarray[3*xarray.length/4]+", "+xarray[0]+", "+xarray[xarray.length-1]);
        IJ.log("Y stat: "+yarray[xarray.length/2]+", "+yarray[xarray.length/4]+", "+yarray[3*xarray.length/4]+", "+yarray[0]+", "+yarray[xarray.length-1]);
        IJ.log("Z stat: "+zarray[xarray.length/2]+", "+zarray[xarray.length/4]+", "+zarray[3*xarray.length/4]+", "+zarray[0]+", "+zarray[xarray.length-1]);
        IJ.log("P stat: "+parray[xarray.length/2]+", "+parray[xarray.length/4]+", "+parray[3*xarray.length/4]+", "+parray[0]+", "+parray[xarray.length-1]);
        IJ.log("B stat: "+barray[xarray.length/2]+", "+barray[xarray.length/4]+", "+barray[3*xarray.length/4]+", "+barray[0]+", "+barray[xarray.length-1]);
        
        IJ.log("N="+nbLoc+" localizations");


        /////////////////////NOW COMPUTE HISTOGRAMS
        int binHistogram=(int)(4*meanSTD3d/binSize);
        double [] histbin = new double [binHistogram];
        for (int i=0;i<binHistogram;i++){
            histbin[i]=binSize*(double)i-binSize*(double)(binHistogram/2);
        }
        double [] histX = new double [binHistogram];
        double [] histY = new double [binHistogram];
        double [] histZ = new double [binHistogram];
        int idpt=0;
        t=0;
        
        for (int ii=0;ii<nbpoints;ii++){
            for (int i=0;i<dx[ii].size();i++){
            
            
            
                a=(meanX[ii]-(double)dx[ii].get(i))*(meanX[ii]-(double)dx[ii].get(i));
                b=(meanY[ii]-(double)dy[ii].get(i))*(meanY[ii]-(double)dy[ii].get(i));
                c=(meanZ[ii]-(double)dz[ii].get(i))*(meanZ[ii]-(double)dz[ii].get(i));
                d=a+b+c;


                boolean keep=true;
                if (Math.sqrt(a)>2*meanSTD3d){
                    keep=false;
                }
                if (Math.sqrt(b)>2*meanSTD3d){
                    keep=false;
                }
                if (Math.sqrt(c)>2*meanSTD3d){
                    keep=false;
                }
                if (keep){

                    int idX=(binHistogram/2)+(int)Math.round(Math.signum((meanX[ii]-(double)dx[ii].get(i)))*Math.sqrt(a)/binSize);
                    int idY=(binHistogram/2)+(int)Math.round(Math.signum((meanY[ii]-(double)dy[ii].get(i)))*Math.sqrt(b)/binSize);
                    int idZ=(binHistogram/2)+(int)Math.round(Math.signum((meanZ[ii]-(double)dz[ii].get(i)))*Math.sqrt(c)/binSize);

                    if ((idX>=0)&&(idX<binHistogram)){
                        histX[idX]++;
                    }
                    
                    if ((idY>=0)&&(idY<binHistogram)){
                        histY[idY]++;
                    }

                    if ((idZ>=0)&&(idZ<binHistogram)){
                        histZ[idZ]++;
                    }


                }
            }
            
        }
        
        plotHist(histbin,histX,"Histogram", "X axis (nm)", "occurrence #");
        plotHist(histbin,histY,"Histogram", "Y axis (nm)", "occurrence #");
        plotHist(histbin,histZ,"Histogram", "Z axis (nm)", "occurrence #");
        
    }
    
    
    
    
    
    
    
    
    
    public void plot(double [] x, double [] y,double [] xline,double [] yline,String title,String xlabel,String ylabel){
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
        
        for (int ii=0;ii<Math.min(x.length,10000);ii++){
            
            p.setColor(Color.blue);
            p.add("CIRCLE",  x,y);
            p.setColor(Color.red);
            //p.show();


            
        }
        p.setFont(0, 18);
        p.show();
        
    }
    
    
    
    
    
    
    
    
    public void plotHist(double [] x, double [] y,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,x,y,3);
        p.setFont(0, 18);
        //p.setAxisLabelFont(0, 18);
        p.show();
        
    }
    

      
}
