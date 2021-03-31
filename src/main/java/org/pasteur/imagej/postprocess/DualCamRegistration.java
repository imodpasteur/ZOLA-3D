/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;
import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.data.*;
import ij.IJ;
import ij.gui.Plot;
import java.util.Arrays;
import java.util.Comparator;
import java.util.ArrayList;
/**
 *
 * @author benoit
 */
public class DualCamRegistration {
    
    StackLocalization sl1;
    StackLocalization sl2;
    
    StackLocalization sl1fuse;
    StackLocalization sl2fuse;
    StackLocalization sl1unfuse;
    StackLocalization sl2unfuse;
    
    String pathRes1;
    String pathRes2;
    
    String pathResFusion;
    
    double maxDistanceMergingXY;
    double maxDistanceMergingZ=Double.MAX_VALUE;
    
    int [][] idLoc;
    
    
    double sizePix=120;//nm
    
    
    double maxX=Double.NEGATIVE_INFINITY;
    double minX=Double.MAX_VALUE;
    double maxY=Double.NEGATIVE_INFINITY;
    double minY=Double.MAX_VALUE;
       
    
    public DualCamRegistration(StackLocalization sl1,StackLocalization sl2,String pathRes1,String pathRes2,String pathResFusion,double maxDistanceMergingXY){
        
        this.sl1=sl1;
        this.sl2=sl2;
        
        
        
        this.maxDistanceMergingXY=maxDistanceMergingXY;
        this.maxDistanceMergingZ=maxDistanceMergingZ;
        
        this.pathRes1=pathRes1;
        this.pathRes2=pathRes2;
        this.pathResFusion=pathResFusion;
    }
    
    
    
    
    public void run(){
        
        
        
        
        IJ.showProgress(0);
        
        
        org.pasteur.imagej.postprocess.RegistrationCrossCorrel reg = new org.pasteur.imagej.postprocess.RegistrationCrossCorrel(sl1,sl2,(int)Math.ceil(sizePix));
        double [] lateralShift=reg.getShift();
        
        sl1fuse=new StackLocalization();
        sl2fuse=new StackLocalization();
        sl1unfuse=new StackLocalization(sl1);
        sl2unfuse=new StackLocalization(sl2);
        
        
        
        
        int [][] idFrameCam1=new int[sl1.fl.size()][2];
        for (int i=0;i<sl1.fl.size();i++){
            idFrameCam1[i][0]=i;
            idFrameCam1[i][1]=sl1.fl.get(i).numFrame;
        }
        Arrays.sort(idFrameCam1, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                return ((Integer) o1[1]).compareTo(o2[1]);
            }
        });
        
        int [][] idFrameCam2=new int[sl2.fl.size()][2];
        for (int i=0;i<sl2.fl.size();i++){
            idFrameCam2[i][0]=i;
            idFrameCam2[i][1]=sl2.fl.get(i).numFrame;
        }
        Arrays.sort(idFrameCam2, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                return ((Integer) o1[1]).compareTo(o2[1]);
            }
        });
        
        
        double x,y;
        
        for (int i=0;i<sl1.fl.size();i++){
            for (int j=0;j<sl1.fl.get(i).loc.size();j++){
                x=sl1.fl.get(i).loc.get(j).X;
                y=sl1.fl.get(i).loc.get(j).Y;
                if ((x>=0)&&(y>=0)){
                    if (maxX<x){
                        maxX=x;
                    }
                    if (maxY<y){
                        maxY=y;
                    }
                    if (minX>x){
                        minX=x;
                    }
                    if (minY>y){
                        minY=y;
                    }
                }
            }
        }
        for (int i=0;i<sl2.fl.size();i++){
            for (int j=0;j<sl2.fl.get(i).loc.size();j++){
                x=sl2.fl.get(i).loc.get(j).X;
                y=sl2.fl.get(i).loc.get(j).Y;
                
                if ((x>=0)&&(y>=0)){
                    if (maxX<x){
                        maxX=x;
                    }
                    if (maxY<y){
                        maxY=y;
                    }
                    if (minX>x){
                        minX=x;
                    }
                    if (minY>y){
                        minY=y;
                    }
                }
            }
        }
        
        int width=(int)Math.ceil((maxX-minX)/sizePix);
        int height=(int)Math.ceil((maxY-minY)/sizePix);
        //IJ.log("width "+width+"  "+height);
        
        //IJ.log("max min "+maxX+"  "+minX+"  "+sizePix);
        idLoc=new int [width][height];
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                idLoc[i][ii]=-1;
            }
        }
        
        ArrayList al = new ArrayList();
        
        int idLocFuse=0;
        int notFuseCam2=0;
        int nbLocCam1=0;
        int fuseCam2=0;
        int i1=0;
        int nbNegative=0;
        for (int i2=0;i2<idFrameCam2.length;i2++){
            IJ.showProgress(.4*(double)i2/(double)idFrameCam2.length);
            int id2=idFrameCam2[i2][0];
            int numframe2=idFrameCam2[i2][1];
            boolean frameFound=false;
            int numframe1=-1;
            int id1=-1;
            searchloop:for (;i1<idFrameCam1.length;i1++){//because sorted, we do not start from the beginning
                id1=idFrameCam1[i1][0];
                numframe1=idFrameCam1[i1][1];
                if (numframe1==numframe2){
                    frameFound=true;
                    break searchloop;
                }
                if (numframe1>numframe2){
                    break searchloop;
                }
            }
            
            
            
            int shift=(int)Math.ceil(maxDistanceMergingXY/sizePix);
            int latshiftX=(int)Math.ceil(lateralShift[0]/sizePix);
            int latshiftY=(int)Math.ceil(lateralShift[1]/sizePix);
            
            int posX1,posY1,posX2,posY2;
            double distance,distZ,xx,yy,zz;
            double prevDistance=Double.MAX_VALUE;
            
            //if frameFound -> id1 for cam 1 corresponds to id2 for cam 2
            if (frameFound){
                FrameLocalization fl1fuse=new FrameLocalization(numframe1);
                FrameLocalization fl2fuse=new FrameLocalization(numframe1);
                FrameLocalization fl1unfuse=new FrameLocalization(numframe1);
                FrameLocalization fl2unfuse=new FrameLocalization(numframe1);
                boolean atLeastOnePartFound=false;
                
                //IJ.log("frame found  "+numframe1+"  "+numframe2);
                
                
                for (int j=0;j<sl1.fl.get(id1).loc.size();j++){
                    if ((sl1.fl.get(id1).loc.get(j).X>=0)&&(sl1.fl.get(id1).loc.get(j).Y>=0)){
                        posX1=(int)((sl1.fl.get(id1).loc.get(j).X-minX)/sizePix);
                        posY1=(int)((sl1.fl.get(id1).loc.get(j).Y-minY)/sizePix);
                        
                        if ((posX1>=0)&&(posY1>=0)&&(posX1<width)&&(posY1<height)){
                            idLoc[posX1][posY1]=j;
                            nbLocCam1++;
                        }
                        else{
                            IJ.log("wrong value "+sl1.fl.get(id1).loc.get(j).X+"  "+sl1.fl.get(id1).loc.get(j).Y+"  "+posX1+"  "+posY1);
                        }
                    }
                    else{
                        nbNegative++;
                    }
                }
                
                for (int j=0;j<sl2.fl.get(id2).loc.size();j++){
                    
                    posX2=(int)((sl2.fl.get(id2).loc.get(j).X-minX)/sizePix);
                    posY2=(int)((sl2.fl.get(id2).loc.get(j).Y-minY)/sizePix);
                    int idPartFound=-1;
                    if ((sl2.fl.get(id2).loc.get(j).X>=0)&&(sl2.fl.get(id2).loc.get(j).Y>=0)){
                        //search in the matrix the closest
                        for (int u=posX2-shift-latshiftX;u<=posX2+shift-latshiftX;u++){
                            for (int v=posY2-shift-latshiftY;v<=posY2+shift-latshiftY;v++){
                                if ((u>=0)&&(v>=0)&&(u<width)&&(v<height)){
                                    if (idLoc[u][v]!=-1){
                                        if (idPartFound==-1){
                                            //new particle


                                            xx=((sl2.fl.get(id2).loc.get(j).X)-lateralShift[0] - (sl1.fl.get(id1).loc.get(idLoc[u][v]).X));
                                            yy=((sl2.fl.get(id2).loc.get(j).Y)-lateralShift[1] - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Y));
                                            zz=((sl2.fl.get(id2).loc.get(j).Z) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Z));
                                            distZ=Math.sqrt(zz*zz);
                                            distance=Math.sqrt(xx*xx+yy*yy);
                                            if ((distance<maxDistanceMergingXY)&&(distZ<maxDistanceMergingZ)){
                                                
                                                idPartFound=idLoc[u][v];
                                                prevDistance=distance;
                                            }
                                        }
                                        else{
                                            //not new -> compare distance
                                            xx=((sl2.fl.get(id2).loc.get(j).X)-lateralShift[0] - (sl1.fl.get(id1).loc.get(idLoc[u][v]).X));
                                            yy=((sl2.fl.get(id2).loc.get(j).Y)-lateralShift[1] - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Y));
                                            zz=((sl2.fl.get(id2).loc.get(j).Z) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Z));
                                            distZ=Math.sqrt(zz*zz);
                                            distance=Math.sqrt(xx*xx+yy*yy);
                                            if ((distance<maxDistanceMergingXY)&&(distZ<maxDistanceMergingZ)&&(distance<prevDistance)){
                                                idPartFound=idLoc[u][v];
                                                prevDistance=distance;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else{
                        nbNegative++;
                    }
                    
                    if (idPartFound!=-1){//it has been found
                        fuseCam2++;
                        PLocalization p1=sl1.fl.get(id1).loc.get(idPartFound).copy();
                        p1.id=idLocFuse;
                        PLocalization p2=sl2.fl.get(id2).loc.get(j).copy();
                        p2.id=idLocFuse;
                        fl1fuse.loc.add(p1);
                        fl2fuse.loc.add(p2);
                        try{
                        sl1unfuse.fl.get(id1).loc.get(idPartFound).exists=false;
                        sl2unfuse.fl.get(id2).loc.get(j).exists=false;
                        }catch(Exception ee){
                            IJ.log("size problem : "+sl1unfuse.fl.size());
                            IJ.log("size problem : "+sl1unfuse.fl.get(id1).loc.size());
                        }
                        
                        idLocFuse++;
                        atLeastOnePartFound=true;
                        //IJ.log("fuse "+(sl2.fl.get(id2).loc.get(j).X)+"  "+(sl1.fl.get(id1).loc.get(idPartFound).X)+"                "+(sl2.fl.get(id2).loc.get(j).Y)+"   "+(sl1.fl.get(id1).loc.get(idPartFound).Y));
                    }
                    else{
                        notFuseCam2++;
                    }
                    
                    
                }
                
                
                
                for (int j=0;j<sl1.fl.get(id1).loc.size();j++){
                    posX1=(int)((sl1.fl.get(id1).loc.get(j).X-minX)/sizePix);
                    posY1=(int)((sl1.fl.get(id1).loc.get(j).Y-minY)/sizePix);
                    if ((posX1>=0)&&(posY1>=0)&&(posX1<width)&&(posY1<height)){
                        idLoc[posX1][posY1]=-1;
                    }
                }
                
                if (atLeastOnePartFound){
                    sl1fuse.fl.add(fl1fuse);
                    sl2fuse.fl.add(fl2fuse);
                }
                
            }
            else{
                //frame not found
            }
            
            
            
            
        }
        
        
        for (int i=0;i<sl1unfuse.fl.size();i++){
            for (int ii=0;ii<sl1unfuse.fl.get(i).loc.size();ii++){
                PLocalization p=sl1unfuse.fl.get(i).loc.get(ii);
                if (!p.exists){
                    sl1unfuse.fl.get(i).loc.remove(ii);
                }
            }
        }
        for (int i=0;i<sl2unfuse.fl.size();i++){
            for (int ii=0;ii<sl2unfuse.fl.get(i).loc.size();ii++){
                PLocalization p=sl2unfuse.fl.get(i).loc.get(ii);
                if (!p.exists){
                    sl2unfuse.fl.get(i).loc.remove(ii);
                }
            }
        }
        
        
        
        
        IJ.log("fuse amount "+fuseCam2+"  /  "+nbLocCam1+" (cam1)  | "+(fuseCam2+notFuseCam2)+" (cam2)");
        
        
        if ((nbNegative/2)>fuseCam2){
            IJ.log("wow, there is a surprisingly high number of negative coordinates: "+(nbNegative/2));
        }
        
        register();
        IJ.showProgress(.6);
        fuseMean();
        
        
        IJ.showProgress(0);
    }
    
    
    public void fuseMean(){
        double x1,y1,z1;
        double cx1,cy1,cz1;
        double x2,y2,z2;
        double cx2,cy2,cz2;
        double cx,cy,cz;
        for (int i=0;i<sl1fuse.fl.size();i++){
            IJ.showProgress(.9+.1*(double)i/(double)sl1fuse.fl.size());
            for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){
                x1=sl1fuse.fl.get(i).loc.get(j).X;
                y1=sl1fuse.fl.get(i).loc.get(j).Y;
                z1=sl1fuse.fl.get(i).loc.get(j).Z;
                x2=sl2fuse.fl.get(i).loc.get(j).X;
                y2=sl2fuse.fl.get(i).loc.get(j).Y;
                z2=sl2fuse.fl.get(i).loc.get(j).Z;
                
                /*if (Math.abs(z1-z2)>100){//xu et al
                    sl1fuse.fl.get(i).loc.get(j).exists=false;
                }
                
                if (true){//xu et al
                    cx1=sl1fuse.fl.get(i).loc.get(j).I;
                    cy1=sl1fuse.fl.get(i).loc.get(j).I;
                    cz1=sl1fuse.fl.get(i).loc.get(j).I;
                    cx2=sl2fuse.fl.get(i).loc.get(j).I;
                    cy2=sl2fuse.fl.get(i).loc.get(j).I;
                    cz2=sl2fuse.fl.get(i).loc.get(j).I;
                }
                else*/
                
                
                if ((sl1fuse.fl.get(i).loc.get(j).crlb_X>0)&&(sl1fuse.fl.get(i).loc.get(j).crlb_Y>0)&&(sl1fuse.fl.get(i).loc.get(j).crlb_Z>0)&&(sl2fuse.fl.get(i).loc.get(j).crlb_X>0)&&(sl2fuse.fl.get(i).loc.get(j).crlb_Y>0)&&(sl2fuse.fl.get(i).loc.get(j).crlb_Z>0)){
                    cx1=1/sl1fuse.fl.get(i).loc.get(j).crlb_X;
                    cx1=cx1*cx1;
                    cy1=1/sl1fuse.fl.get(i).loc.get(j).crlb_Y;
                    cy1=cy1*cy1;
                    cz1=1/sl1fuse.fl.get(i).loc.get(j).crlb_Z;
                    cz1=cz1*cz1;
                    cx2=1/sl2fuse.fl.get(i).loc.get(j).crlb_X;
                    cx2=cx2*cx2;
                    cy2=1/sl2fuse.fl.get(i).loc.get(j).crlb_Y;
                    cy2=cy2*cy2;
                    cz2=1/sl2fuse.fl.get(i).loc.get(j).crlb_Z;
                    cz2=cz2*cz2;
                }
                else if ((sl1fuse.fl.get(i).loc.get(j).I>0)&&(sl2fuse.fl.get(i).loc.get(j).I>0)){
                    cx1=sl1fuse.fl.get(i).loc.get(j).I;
                    cy1=sl1fuse.fl.get(i).loc.get(j).I;
                    cz1=sl1fuse.fl.get(i).loc.get(j).I;
                    cx2=sl2fuse.fl.get(i).loc.get(j).I;
                    cy2=sl2fuse.fl.get(i).loc.get(j).I;
                    cz2=sl2fuse.fl.get(i).loc.get(j).I;
                }
                else{
                    cx1=.5;
                    cy1=.5;
                    cz1=.5;
                    cx2=.5;
                    cy2=.5;
                    cz2=.5;
                }
                cx=cx1+cx2;
                cy=cy1+cy2;
                cz=cz1+cz2;
                sl1fuse.fl.get(i).loc.get(j).X=(x1*cx1+x2*cx2)/cx;
                sl1fuse.fl.get(i).loc.get(j).Y=(y1*cy1+y2*cy2)/cy;
                sl1fuse.fl.get(i).loc.get(j).Z=(z1*cz1+z2*cz2)/cz;
                sl1fuse.fl.get(i).loc.get(j).crlb_X=Math.sqrt(1/cx);
                sl1fuse.fl.get(i).loc.get(j).crlb_Y=Math.sqrt(1/cy);
                sl1fuse.fl.get(i).loc.get(j).crlb_Z=Math.sqrt(1/cz);
                sl1fuse.fl.get(i).loc.get(j).I+=sl2fuse.fl.get(i).loc.get(j).I;
            }
        }
        
        
        if ((pathRes1.length()>2)){
            sl1fuse.save(pathResFusion);
        }
    }
    
    
    
    
    
    
    public void register(){
        
        
        int nbData=0;
        for (int i=0;i<sl1fuse.fl.size();i++){
            nbData+=sl1fuse.fl.get(i).loc.size();
        }
        IJ.log("nbData "+nbData);
        
        double maxdistBeforeXY=Double.NEGATIVE_INFINITY;
        double maxdistAfterXY=Double.NEGATIVE_INFINITY;
        double [] distanceBeforeXY = new double  [nbData];
        double [] distanceAfterXY = new double  [nbData];
        double maxdistBeforeZ=Double.NEGATIVE_INFINITY;
        double maxdistAfterZ=Double.NEGATIVE_INFINITY;
        double [] distanceBeforeZ = new double  [nbData];
        double [] distanceAfterZ = new double  [nbData];
        double [][] data1 = new double [3][nbData];
        double [][] data2 = new double [3][nbData];
        boolean notApparied=false;
        
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        for (int i=0,t=0;i<sl1fuse.fl.size();i++){
            IJ.showProgress(.4+.2*(double)i/(double)sl1fuse.fl.size());
            for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){
                if (((sl1fuse.fl.get(i).loc.get(j).id!=sl2fuse.fl.get(i).loc.get(j).id))||((sl1fuse.fl.get(i).loc.get(j).frame!=sl2fuse.fl.get(i).loc.get(j).frame))){
                    notApparied=true;
                }
                data1[0][t]=sl1fuse.fl.get(i).loc.get(j).X;
                data1[1][t]=sl1fuse.fl.get(i).loc.get(j).Y;
                data1[2][t]=sl1fuse.fl.get(i).loc.get(j).Z;
                data2[0][t]=sl2fuse.fl.get(i).loc.get(j).X;
                data2[1][t]=sl2fuse.fl.get(i).loc.get(j).Y;
                data2[2][t]=sl2fuse.fl.get(i).loc.get(j).Z;
                
                if (Math.min(data1[0][t],data2[0][t])<minX){
                    minX=Math.min(data1[0][t],data2[0][t]);
                }
                if (Math.min(data1[1][t],data2[1][t])<minY){
                    minY=Math.min(data1[1][t],data2[1][t]);
                }
                if (Math.min(data1[2][t],data2[2][t])<minZ){
                    minZ=Math.min(data1[2][t],data2[2][t]);
                }
                if (Math.max(data1[0][t],data2[0][t])>maxX){
                    maxX=Math.max(data1[0][t],data2[0][t]);
                }
                if (Math.max(data1[1][t],data2[1][t])>maxY){
                    maxY=Math.max(data1[1][t],data2[1][t]);
                }
                if (Math.max(data1[2][t],data2[2][t])>maxZ){
                    maxZ=Math.max(data1[2][t],data2[2][t]);
                }
                    
                distanceBeforeXY[t]=Math.sqrt((data1[0][t]-data2[0][t])*(data1[0][t]-data2[0][t])+(data1[1][t]-data2[1][t])*(data1[1][t]-data2[1][t]));
                distanceBeforeZ[t]=Math.sqrt((data1[2][t]-data2[2][t])*(data1[2][t]-data2[2][t]));
                if (distanceBeforeXY[t]>maxdistBeforeXY){
                    maxdistBeforeXY=distanceBeforeXY[t];
                }
                if (distanceBeforeZ[t]>maxdistBeforeZ){
                    maxdistBeforeZ=distanceBeforeZ[t];
                }
                t++;
            }
        }
        if (notApparied){
            IJ.log("oops : the two files are not apparied -> registration not possible");
        }
        else{
            PolynomialFit pf = new PolynomialFit(2,data1,data2);
            pf.removeParameter(9);//remove Z^2 parameter
            pf.removeParameter(8);//remove Y.Z component
            pf.removeParameter(7);//remove X.Z component
            pf.run();
            pf.log();
            double [] v = new double [3];
            for (int i=0,t=0;i<sl1fuse.fl.size();i++){
                IJ.showProgress(.7+.2*(double)i/(double)sl1fuse.fl.size());
                for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){
                    v[0]=sl2fuse.fl.get(i).loc.get(j).X;
                    v[1]=sl2fuse.fl.get(i).loc.get(j).Y;
                    v[2]=sl2fuse.fl.get(i).loc.get(j).Z;
                    double [] res=pf.transform(v);
                    sl2fuse.fl.get(i).loc.get(j).X=res[0];
                    sl2fuse.fl.get(i).loc.get(j).Y=res[1];
                    sl2fuse.fl.get(i).loc.get(j).Z=res[2];
                    distanceAfterXY[t]=Math.sqrt((data1[0][t]-res[0])*(data1[0][t]-res[0])+(data1[1][t]-res[1])*(data1[1][t]-res[1]));
                    //distanceAfterXY[t]=Math.sqrt((data1[1][t]-res[1])*(data1[1][t]-res[1]));
                    if (distanceAfterXY[t]>maxdistAfterXY){
                        maxdistAfterXY=distanceAfterXY[t];
                    }
                    
                    distanceAfterZ[t]=Math.sqrt((data1[2][t]-res[2])*(data1[2][t]-res[2]));
                    if (distanceAfterZ[t]>maxdistAfterZ){
                        maxdistAfterZ=distanceAfterZ[t];
                    }
                    t++;
                }
            }
            
            for (int i=0;i<sl2unfuse.fl.size();i++){
                for (int j=0;j<sl2unfuse.fl.get(i).loc.size();j++){
                    v[0]=sl2unfuse.fl.get(i).loc.get(j).X;
                    v[1]=sl2unfuse.fl.get(i).loc.get(j).Y;
                    v[2]=sl2unfuse.fl.get(i).loc.get(j).Z;
                    double [] res=pf.transform(v);
                    sl2unfuse.fl.get(i).loc.get(j).X=res[0];
                    sl2unfuse.fl.get(i).loc.get(j).Y=res[1];
                    sl2unfuse.fl.get(i).loc.get(j).Z=res[2];
                    
                }
            }
        }
        
        if ((pathRes1.length()>2)&&(pathRes2.length()>2)){
            sl1fuse.save(pathRes1);
            sl2fuse.save(pathRes2);
        }
        
        ZRendering.nameHistPlot="camera1";
        ZRendering.hist2D(sl1fuse, 20, minX, maxX, minY, maxY, minZ, maxZ, 1);
        ZRendering.nameHistPlot="camera2";
        ZRendering.hist2D(sl2fuse, 20, minX, maxX, minY, maxY, minZ, maxZ, 1);
        
        
        double [] axisBeforeXY = new double [(int)maxdistBeforeXY+1];
        double [] histBeforeXY = new double [(int)maxdistBeforeXY+1];
        double [] axisAfterXY = new double [(int)maxdistAfterXY+1];
        double [] histAfterXY = new double [(int)maxdistAfterXY+1];
        
        double [] axisBeforeZ = new double [(int)maxdistBeforeZ+1];
        double [] histBeforeZ = new double [(int)maxdistBeforeZ+1];
        double [] axisAfterZ = new double [(int)maxdistAfterZ+1];
        double [] histAfterZ = new double [(int)maxdistAfterZ+1];
        
        for (int i=0;i<(int)maxdistBeforeXY+1;i++){
            axisBeforeXY[i]=i;
            histBeforeXY[i]=0;
        }
        for (int i=0;i<(int)maxdistAfterXY+1;i++){
            axisAfterXY[i]=i;
            histAfterXY[i]=0;
        }
        
        for (int i=0;i<(int)maxdistBeforeZ+1;i++){
            axisBeforeZ[i]=i;
            histBeforeZ[i]=0;
        }
        for (int i=0;i<(int)maxdistAfterZ+1;i++){
            axisAfterZ[i]=i;
            histAfterZ[i]=0;
        }
        
        
        for (int i=0;i<nbData;i++){
            histBeforeXY[(int)distanceBeforeXY[i]]++;
            histAfterXY[(int)distanceAfterXY[i]]++;
        }
        for (int i=0;i<nbData;i++){
            histBeforeZ[(int)distanceBeforeZ[i]]++;
            histAfterZ[(int)distanceAfterZ[i]]++;
        }
        
        
        plotHist(axisBeforeXY,histBeforeXY,"Histogram of X/Y distances (cam1 vs. cam2) before registration", "distance (nm)", "occurrence #");
        plotHist(axisBeforeZ,histBeforeZ,"Histogram of Z distances (cam1 vs. cam2) before registration", "distance (nm)", "occurrence #");
        plotHist(axisAfterXY,histAfterXY,"Histogram of X/Y distances (cam1 vs. cam2) after registration (avg!=0 because of 2D)", "distance (nm)", "occurrence #");
        plotHist(axisAfterZ,histAfterZ,"Histogram of Z distances (cam1 vs. cam2) after registration", "distance (nm)", "occurrence #");
        
    }
    
    
    
    
    public void plotHist(double [] x, double [] y,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,x,y,3);
        p.setFont(0, 18);
        p.setLineWidth(2);
        p.show();
        
    }
    
}
