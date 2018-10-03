/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.development;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.postprocess.ZRendering;
import org.pasteur.imagej.utils.PolynomialFit;
import ij.IJ;
import ij.gui.Plot;
import java.awt.Color;
/**
 *
 * @author benoit
 */
public class FilamentResolution{

    double precisionSamplingPolyCurve=0.01;//(nm)

    public FilamentResolution(StackLocalization sl,double sizePix,double [] xpoints, double [] ypoints, double [] zpoints,int order,double maximumDistance,double binSize,String pathRes,String pathResexclude){
        
        if (order<1){
            order=1;
            IJ.log("WARNING: minimum order is 1");
        }
        
        
        StackLocalization slnew=new StackLocalization();
        StackLocalization slnewexclude=new StackLocalization();
        
        StackLocalization slnewMerge2D=new StackLocalization();
        
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
        
        
        boolean xDirection=true;
        if ((maxy-miny)>(maxx-minx)){
            xDirection=false;
        }
        
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
        
        if ((pathRes!=null)&&(pathRes.length()>1)){
            slnew.save(pathRes);
            IJ.log("localization table with filament saved");
        }
        else{
            IJ.log("localization table with filament not saved");
        }
        
        if ((pathResexclude!=null)&&(pathResexclude.length()>1)){
            slnewexclude.save(pathResexclude);
            IJ.log("localization table without filament saved");
        }
        else{
            IJ.log("localization table without filament not saved");
        }
        
        
        if (nbloc2>0){
            ZRendering.hist3D(slnew, sizePix,0);
            
        }
        else{
            IJ.log("oops, no localization found in the region");
            return;
        }
        
        
        
        
        
        double [][] dx=new double [1][nbloc2];
        double [][] dy=new double [1][nbloc2];
        double [][] dz=new double [1][nbloc2];
        
        double [] photon=new double[nbloc2];
        double [] bckg=new double[nbloc2];
        
        double [] locx=new double [nbloc2];
        double [] locy=new double [nbloc2];
        double [] locz=new double [nbloc2];
        
        double [] locCRLBx=new double [nbloc2];
        double [] locCRLBy=new double [nbloc2];
        double [] locCRLBz=new double [nbloc2];
        
        double minX=Double.POSITIVE_INFINITY;
        double minY=Double.POSITIVE_INFINITY;
        double minZ=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        int t=0;
        for (int i=0;i<slnew.fl.size();i++){
            for (int ii=0;ii<slnew.fl.get(i).loc.size();ii++){
                dx[0][t]=slnew.fl.get(i).loc.get(ii).X;
                dy[0][t]=slnew.fl.get(i).loc.get(ii).Y;
                dz[0][t]=slnew.fl.get(i).loc.get(ii).Z;
                locCRLBx[t]=slnew.fl.get(i).loc.get(ii).crlb_X;
                locCRLBy[t]=slnew.fl.get(i).loc.get(ii).crlb_Y;
                locCRLBz[t]=slnew.fl.get(i).loc.get(ii).crlb_Z;
                if (!xDirection){
                    double tmp=dx[0][t];
                    dx[0][t]=dy[0][t];
                    dy[0][t]=tmp;
                }
                photon[t]=slnew.fl.get(i).loc.get(ii).I;
                bckg[t]=slnew.fl.get(i).loc.get(ii).BG;
                locx[t]=dx[0][t];
                locy[t]=dy[0][t];
                locz[t]=dz[0][t];
                
                if (minX>dx[0][t]){
                    minX=dx[0][t];
                }
                if (minY>dy[0][t]){
                    minY=dy[0][t];
                }
                if (minZ>dz[0][t]){
                    minZ=dz[0][t];
                }
                if (maxX<dx[0][t]){
                    maxX=dx[0][t];
                }
                if (maxY<dy[0][t]){
                    maxY=dy[0][t];
                }
                if (maxZ<dz[0][t]){
                    maxZ=dz[0][t];
                }
                t++;
            }
        }
        
        
        
        int nbValueInPoly=2*(int)((maxX-minX)/this.precisionSamplingPolyCurve);
        
        
        PolynomialFit pfY=null;
        
        
        pfY = new PolynomialFit(order,dy,dx);
        
        pfY.run();
        
        PolynomialFit pfXYLinear= new PolynomialFit(1,dy,dx);
        pfXYLinear.run();
        
        PolynomialFit pfZ = new PolynomialFit(order,dz,dx);
        
        pfZ.run();
        double [] polx=new double[nbValueInPoly];
        double [] poly=new double[nbValueInPoly];
        double [] polz=new double[nbValueInPoly];
        
        
        
        int nbValueInPolyprint=50;
        
        double [] polxprint=new double[nbValueInPolyprint];
        double [] polyprint=new double[nbValueInPolyprint];
        double [] polzprint=new double[nbValueInPolyprint];
        
        double [] v=new double[1];//for poly fit result
        double [] vv=new double[1];//for poly fit result
        for (int i=0;i<nbValueInPoly;i++){
            polx[i]=minX+((double)i)*((maxX-minX)/(((double)nbValueInPoly)/2.))-(maxX-minX)/2;
            
            v[0]=polx[i];
            vv=pfY.transform(v);
            poly[i]=vv[0];
            vv=pfZ.transform(v);
            polz[i]=vv[0];
        }
        
        
        for (int i=0;i<nbValueInPolyprint;i++){
            polxprint[i]=minX+((double)i)*((maxX-minX)/(((double)nbValueInPolyprint)));
            
            v[0]=polxprint[i];
            vv=pfY.transform(v);
            polyprint[i]=vv[0];
            vv=pfZ.transform(v);
            polzprint[i]=vv[0];
        }
        
        
        
        
        
        double precisionX=0;
        double precisionY=0;
        double precisionZ=0;
        double precisionPlan=0;
        double meanI=0;
        double meanBG=0;
        
        
        //compare loc and pol
        int idpt=0,idptz=0;
        double amin, bmin, cmin, aminz;
        double dminz;
        double ddz;
        for (int i=0;i<nbloc2;i++){
            dmin=Double.POSITIVE_INFINITY;
            dminz=Double.POSITIVE_INFINITY;
            amin=Double.POSITIVE_INFINITY;
            bmin=Double.POSITIVE_INFINITY;
            cmin=Double.POSITIVE_INFINITY;
            idpt=0;
            idptz=0;

            int begin=0;
            int end=nbValueInPoly;
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    if (dmin>d){
                        dmin=d;
                        amin=a;
                        bmin=b;
                        idpt=ii;
                    }
                    
                }
                begin=Math.max(0,idpt-step*2);
                end=Math.min(nbValueInPoly,idpt+step*2);
            }
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    ddz=a;
                    
                    if (dminz>ddz){
                        dminz=ddz;
                        cmin=c;
                        aminz=a;
                        idptz=ii;
                    }
                }
                begin=Math.max(0,idptz-step*2);
                end=Math.min(nbValueInPoly,idptz+step*2);
            }
            meanI+=photon[i];
            meanBG+=bckg[i];
            precisionX+=amin;
            precisionY+=bmin;
            precisionZ+=cmin;
            precisionPlan+=amin+bmin;
        }
        
        meanI/=((double)nbloc2);
        meanBG/=((double)nbloc2);
        precisionX=Math.sqrt(precisionX/((double)nbloc2));
        precisionY=Math.sqrt(precisionY/((double)nbloc2));
        precisionZ=Math.sqrt(precisionZ/((double)nbloc2));
        precisionPlan=Math.sqrt(precisionPlan/((double)nbloc2));
        IJ.log("precision using all localizations whose distance is < "+maximumDistance);
        if (xDirection){
            IJ.log("precision X : "+precisionX);
            IJ.log("precision Y : "+precisionY);
        }
        else{
            IJ.log("precision Y : "+precisionX);
            IJ.log("precision X : "+precisionY);
        }
        IJ.log("precision Z : "+precisionZ);
        IJ.log("Plan precision : "+precisionPlan);
        IJ.log("localization number : "+nbloc2);
        IJ.log("mean photon number : "+meanI);
        IJ.log("mean background : "+meanBG);
        
        IJ.log("");
        
        
        
        
        
        int nbloc3=0;
        
        //recompute to reject outliers
        idpt=0;
        
        for (int i=0;i<nbloc2;i++){
            dmin=Double.POSITIVE_INFINITY;
            amin=Double.POSITIVE_INFINITY;
            bmin=Double.POSITIVE_INFINITY;
            cmin=Double.POSITIVE_INFINITY;
            idpt=0;
            
            int begin=0;
            int end=nbValueInPoly;
            idptz=0;
            dminz=Double.POSITIVE_INFINITY;
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    if (dmin>d){
                        dmin=d;
                        amin=a;
                        bmin=b;
                        idpt=ii;
                    }
                    
                }
                begin=Math.max(0,idpt-step*2);
                end=Math.min(nbValueInPoly,idpt+step*2);
            }
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    ddz=a;
                    
                    if (dminz>ddz){
                        dminz=ddz;
                        cmin=c;
                        aminz=a;
                        idptz=ii;
                    }
                }
                begin=Math.max(0,idptz-step*2);
                end=Math.min(nbValueInPoly,idptz+step*2);
            }
            boolean keep=true;
            if (Math.sqrt(amin)>3*precisionX){
                keep=false;
            }
            if (Math.sqrt(bmin)>3*precisionY){
                keep=false;
            }
            if (Math.sqrt(cmin)>3*precisionZ){
                keep=false;
            }
            if (keep){
                nbloc3++;
            }
        }
        
        
        
        
        
        
        double [] locxnew=new double [nbloc3];
        double [] locynew=new double [nbloc3];
        double [] locznew=new double [nbloc3];
        
        double precisionXnew=0;
        double precisionYnew=0;
        double precisionZnew=0;
        double precisionPlannew=0;
        
        double crlbXnew=0;
        double crlbYnew=0;
        double crlbZnew=0;
        double crlbPlannew=0;
        
        double crlbXnewstd=0;
        double crlbYnewstd=0;
        double crlbZnewstd=0;
        double crlbPlannewstd=0;
        meanI=0;
        meanBG=0;
        double meanZ=0;
        
        double maxDmin=Double.NEGATIVE_INFINITY;
        //compare loc and pol
        idpt=0;
        t=0;
        for (int i=0;i<nbloc2;i++){
            dmin=Double.POSITIVE_INFINITY;
            amin=Double.POSITIVE_INFINITY;
            bmin=Double.POSITIVE_INFINITY;
            cmin=Double.POSITIVE_INFINITY;
            idpt=0;

            int begin=0;
            int end=nbValueInPoly;
            idptz=0;
            dminz=Double.POSITIVE_INFINITY;
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    if (dmin>d){
                        dmin=d;
                        amin=a;
                        bmin=b;
                        idpt=ii;
                    }
                    
                }
                begin=Math.max(0,idpt-step*2);
                end=Math.min(nbValueInPoly,idpt+step*2);
            }
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    ddz=a;
                    
                    if (dminz>ddz){
                        dminz=ddz;
                        cmin=c;
                        aminz=a;
                        idptz=ii;
                    }
                }
                begin=Math.max(0,idptz-step*2);
                end=Math.min(nbValueInPoly,idptz+step*2);
            }
            boolean keep=true;
            if (Math.sqrt(amin)>3*precisionX){
                keep=false;
            }
            if (Math.sqrt(bmin)>3*precisionY){
                keep=false;
            }
            if (Math.sqrt(cmin)>3*precisionZ){
                keep=false;
            }
            if (keep){
                meanI+=photon[i];
                meanBG+=bckg[i];
                meanZ+=locz[i];
                locxnew[t]=locx[i];
                locynew[t]=locy[i];
                locznew[t]=locz[i];
                t++;
                precisionXnew+=amin;
                precisionYnew+=bmin;
                precisionZnew+=cmin;
                precisionPlannew+=amin+bmin;
                
                crlbXnew+=locCRLBx[i];
                crlbYnew+=locCRLBy[i];
                crlbZnew+=locCRLBz[i];
                crlbPlannew+=Math.sqrt(locCRLBx[i]*locCRLBx[i]+locCRLBy[i]*locCRLBy[i]);
                
                if (Math.sqrt(dmin)>maxDmin){
                    maxDmin=Math.sqrt(dmin);
                }
                
                
            }
            
        }
        crlbXnew/=((double)nbloc3);
        crlbYnew/=((double)nbloc3);
        crlbZnew/=((double)nbloc3);
        crlbPlannew/=((double)nbloc3);
        
        //std CRLB
        idpt=0;
        t=0;
        for (int i=0;i<nbloc2;i++){
            dmin=Double.POSITIVE_INFINITY;
            amin=Double.POSITIVE_INFINITY;
            bmin=Double.POSITIVE_INFINITY;
            cmin=Double.POSITIVE_INFINITY;
            idpt=0;

            int begin=0;
            int end=nbValueInPoly;
            idptz=0;
            dminz=Double.POSITIVE_INFINITY;
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    if (dmin>d){
                        dmin=d;
                        amin=a;
                        bmin=b;
                        idpt=ii;
                    }
                    
                }
                begin=Math.max(0,idpt-step*2);
                end=Math.min(nbValueInPoly,idpt+step*2);
            }
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    ddz=a;
                    
                    if (dminz>ddz){
                        dminz=ddz;
                        cmin=c;
                        aminz=a;
                        idptz=ii;
                    }
                }
                begin=Math.max(0,idptz-step*2);
                end=Math.min(nbValueInPoly,idptz+step*2);
            }
            boolean keep=true;
            if (Math.sqrt(amin)>3*precisionX){
                keep=false;
            }
            if (Math.sqrt(bmin)>3*precisionY){
                keep=false;
            }
            if (Math.sqrt(cmin)>3*precisionZ){
                keep=false;
            }
            if (keep){
                
                crlbXnewstd+=(crlbXnew-locCRLBx[i])*(crlbXnew-locCRLBx[i]);
                crlbYnewstd+=(crlbYnew-locCRLBy[i])*(crlbYnew-locCRLBy[i]);
                crlbZnewstd+=(crlbZnew-locCRLBz[i])*(crlbZnew-locCRLBz[i]);
                crlbPlannewstd+=(Math.sqrt(locCRLBx[i]*locCRLBx[i]+locCRLBy[i]*locCRLBy[i])-crlbPlannew)*(Math.sqrt(locCRLBx[i]*locCRLBx[i]+locCRLBy[i]*locCRLBy[i])-crlbPlannew);
                
            }
            
        }
        
        crlbXnewstd=Math.sqrt(crlbXnewstd/((double)nbloc3));
        crlbYnewstd=Math.sqrt(crlbYnewstd/((double)nbloc3));
        crlbZnewstd=Math.sqrt(crlbZnewstd/((double)nbloc3));
        crlbPlannewstd=Math.sqrt(crlbPlannewstd/((double)nbloc3));
        
        
        
        
        
        
        
        meanI/=((double)nbloc3);
        meanBG/=((double)nbloc3);
        meanZ/=((double)nbloc3);
        precisionXnew=Math.sqrt(precisionXnew/((double)nbloc3));
        precisionYnew=Math.sqrt(precisionYnew/((double)nbloc3));
        precisionZnew=Math.sqrt(precisionZnew/((double)nbloc3));
        precisionPlannew=Math.sqrt(precisionPlannew/((double)nbloc3));
        
        
        IJ.log("precision after rejecting 3 standars deviation");
        if (xDirection){
            IJ.log("precision X : "+precisionXnew);
            IJ.log("precision Y : "+precisionYnew);
        }
        else{
            IJ.log("precision Y : "+precisionXnew);
            IJ.log("precision X : "+precisionYnew);
        }
        IJ.log("Plan precision : "+precisionPlannew);
        IJ.log("precision Z : "+precisionZnew);
        
        IJ.log("localization number : "+nbloc3);
        IJ.log("mean photon number : "+meanI);
        IJ.log("mean background : "+meanBG);
        IJ.log("mean Z position: "+meanZ);
        
        
        IJ.log("CRLB X : "+crlbXnew+"  +- "+crlbXnewstd);
        IJ.log("CRLB Y : "+crlbYnew+"  +- "+crlbYnewstd);
        //moyenne weighted according to slope
        double slope=Math.abs(pfXYLinear.a[0][1]);
        if (!xDirection){
            if (slope!=0){
                slope=1/slope;
            }
        }
        IJ.log("slope XY "+slope);
        if (slope==0){
            IJ.log("CRLB MeanWeighted(X,Y) : "+crlbYnew);
        }
        else if (slope==1){
            IJ.log("CRLB MeanWeighted(X,Y) : "+crlbXnew);
        }
        else{
            IJ.log("CRLB MeanWeighted(X,Y) : "+(1/(slope+(1/slope)))*(slope*crlbXnew+(crlbYnew/slope)));
        }
        IJ.log("CRLB Plan X*Y : "+crlbPlannew+"  +- "+crlbPlannewstd);
        IJ.log("CRLB Z : "+crlbZnew+"  +- "+crlbZnewstd);
        //IJ.log("Plan CRLB (unreal:crlbX+crlbY...): "+Math.sqrt(crlbPlannew/((double)nbloc3)));
        //IJ.log("Global CRLB (unreal:crlbX+crlbY+crlbZ...): "+Math.sqrt(crlbGlobalnew/((double)nbloc3)));
        
        IJ.log("plot... please wait");
        
//        if (xDirection){
//            plot(locxnew,locynew,polxprint,polyprint,"polynomial fit Y=f(X)", "X localization", "Y localization");
//            plot(locxnew,locznew,polxprint,polzprint,"polynomial fit Z=f(X)", "X localization", "Z localization");
//        }
//        else{
//            plot(locxnew,locynew,polxprint,polyprint,"polynomial fit X=f(Y)", "Y localization", "X localization");
//            plot(locxnew,locznew,polxprint,polzprint,"polynomial fit Z=f(Y)", "Y localization", "Z localization");
//        }
        
        
        
        
        FrameLocalization flnewMerge2D = new FrameLocalization(0);
        
        
        //for histogram along filament
        double [] distBinned=new double[polx.length-1];
        double sumDist=0;
        for (int i=0;i<polx.length-1;i++){
            distBinned[i]=Math.sqrt((polx[i]-polx[i+1])*(polx[i]-polx[i+1])+(poly[i]-poly[i+1])*(poly[i]-poly[i+1])+(polz[i]-polz[i+1])*(polz[i]-polz[i+1]));
            sumDist+=distBinned[i];
        }
        int binNumber=(int)Math.ceil(sumDist/binSize);
        int [] idAlongFil=new int[polx.length];
        double [] histAlongFil=new double[binNumber];
        double [] histAlongFilAxis=new double[binNumber];
        for (int i=0;i<histAlongFilAxis.length;i++){
            histAlongFilAxis[i]=(double)i*binSize;
        }
        sumDist=0;
        for (int i=0;i<polx.length-1;i++){
            double distBin=Math.sqrt((polx[i]-polx[i+1])*(polx[i]-polx[i+1])+(poly[i]-poly[i+1])*(poly[i]-poly[i+1])+(polz[i]-polz[i+1])*(polz[i]-polz[i+1]));
            sumDist+=distBin;
            idAlongFil[i]=(int)(sumDist/binSize);
        }
        
        /////////////////////NOW COMPUTE HISTOGRAMS
        int binHistogram=(int)(6*precisionPlannew/binSize);
        double [] histbin = new double [binHistogram];
        for (int i=0;i<binHistogram;i++){
            histbin[i]=binSize*(double)i-binSize*(double)(binHistogram/2);
        }
        //double [] histXY = new double [binHistogram];
        double [] histY = new double [binHistogram];
        double [] histPlanar = new double [binHistogram];
        double [] histZ = new double [binHistogram];
        idpt=0;
        idptz=0;
        t=0;
        for (int i=0;i<nbloc2;i++){
            dmin=Double.POSITIVE_INFINITY;
            amin=Double.POSITIVE_INFINITY;
            bmin=Double.POSITIVE_INFINITY;
            cmin=Double.POSITIVE_INFINITY;
            
            double signXY=1;
            double signZ=1;
            idpt=0;

            int begin=1;
            int end=nbValueInPoly-1;
            idptz=0;
            dminz=Double.POSITIVE_INFINITY;
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    if (dmin>d){
                        signXY=Math.signum((poly[ii]-locy[i]));
                        dmin=d;
                        amin=a;
                        bmin=b;
                        idpt=ii;
                    }
                    
                }
                begin=Math.max(0,idpt-step*2);
                end=Math.min(nbValueInPoly,idpt+step*2);
            }
            for (int step=10000;step>=1;step/=10){
                
                for (int ii=begin;ii<end;ii+=step){
                    a=(polx[ii]-locx[i])*(polx[ii]-locx[i]);
                    b=(poly[ii]-locy[i])*(poly[ii]-locy[i]);
                    c=(polz[ii]-locz[i])*(polz[ii]-locz[i]);
                    d=a+b;
                    ddz=a;
                    
                    if (dminz>ddz){
                        signZ=Math.signum((polz[ii]-locz[i]));
                        dminz=ddz;
                        cmin=c;
                        aminz=a;
                        idptz=ii;
                    }
                }
                begin=Math.max(0,idptz-step*2);
                end=Math.min(nbValueInPoly,idptz+step*2);
            }
            
            boolean keep=true;
            if (Math.sqrt(amin)>2*precisionX){
                keep=false;
            }
            if (Math.sqrt(bmin)>2*precisionY){
                keep=false;
            }
            if (Math.sqrt(cmin)>2*precisionZ){
                keep=false;
            }
            if (keep){
                //int idXY=(binHistogram/2)+(int)(sign*Math.sqrt(amin+bmin)/binSize);
                int idY=(binHistogram/2)+(int)Math.round(signXY*Math.sqrt(bmin)/binSize);
                int idP=(binHistogram/2)+(int)Math.round(signXY*Math.sqrt(amin+bmin)/binSize);
                
                
                int idZ=(binHistogram/2)+(int)Math.round(signZ*Math.sqrt(cmin)/binSize);
                //if ((idXY>=0)&&(idXY<binHistogram)){
                //    histXY[idXY]++;
                //}
                
                if (idpt<polx.length-1){
                    if ((idAlongFil[idpt]>=0)&&(idAlongFil[idpt]<histAlongFil.length)){
                        
                        histAlongFil[idAlongFil[idpt]]++;
                    }
                }
                
                if ((idY>=0)&&(idY<binHistogram)){
                    histY[idY]++;
                }
                
                if ((idP>=0)&&(idP<binHistogram)){
                    histPlanar[idP]++;
                }
                
                if ((idZ>=0)&&(idZ<binHistogram)){
                    histZ[idZ]++;
                }
                
                //angle computed according to Z axis (because filament supposed orthogonal to axial direction)
                double lx=locx[i]-polx[idpt];
                double ly=locy[i]-poly[idpt];
                double lz=locz[i]-polz[idptz];
                
                double ox=0;
                double oy=0;
                double oz=0;
                
                double norme=Math.sqrt(lx*lx+ly*ly+lz*lz);
                double alpha=Math.acos(lz/norme);
                if (signXY<0){
                    alpha+=Math.PI;
                }
                double dist=Math.sqrt(amin+bmin+cmin);
                double py=dist*Math.cos(alpha);
                double px=dist*Math.sin(alpha);
                
                //IJ.log("alpha "+alpha+"  "+px+"  "+py+"  "+maxDmin+"  "+(maxDmin+px)+"  "+(maxDmin+py));
                PLocalization p= new PLocalization(0,maxDmin+px,maxDmin+py,0);
                flnewMerge2D.loc.add(p);
                
            }
            
        }
        
        slnewMerge2D.fl.add(flnewMerge2D);
        ZRendering.hist2D(slnewMerge2D, binSize, 0);
        
        if (xDirection){
            //plotHist(histbin,histY,"Histogram", "Y axis (nm)", "occurrence #");
            
        }
        else{
            //plotHist(histbin,histY,"Histogram", "X axis (nm)", "occurrence #");
            
        }
        
        plotHist(histAlongFilAxis,histAlongFil,"Histogram along segment", "X/Y segment (nm)", "occurrence #");
        plotHist(histbin,histZ,"Histogram", "Z axis (nm)", "occurrence #");
        plotHist(histbin,histPlanar,"Histogram", "X/Y axis (nm)", "occurrence #");
        
        
        
    }
    
    
    
    
    
    
    
    
    
    public void plot(double [] x, double [] y,double [] xline,double [] yline,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,xline,yline);
        
        p.setFont(0, 18);
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
            p.setLineWidth(1);
            p.setColor(Color.blue);
            p.add("CIRCLE",  x,y);
            p.setColor(Color.red);
            p.setLineWidth(2);
            //p.show();


            
        }
        p.show();
        
    }
    
    
    
    
    
    
    
    
    public void plotHist(double [] x, double [] y,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,x,y,3);
        p.setFont(0, 18);
        p.setLineWidth(2);
        p.show();
        
    }
    

      
}
