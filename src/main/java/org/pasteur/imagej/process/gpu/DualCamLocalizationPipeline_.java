/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.gpu;
import org.pasteur.imagej.utils.*;
import org.pasteur.imagej.data.StackLocalization;
import org.pasteur.imagej.data.FrameLocalization;
import org.pasteur.imagej.data.PLocalization;
import ij.IJ;
import ij.ImageStack;
import java.util.Arrays;
import java.util.Comparator;
import java.util.ArrayList;
import ij.process.ImageProcessor;

import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author benoit
 */
public class DualCamLocalizationPipeline_ {
    
    
    ArrayList<String> otherVariableName = new ArrayList<String>();
    
    StackLocalization stackloc;
    boolean saveOnTheFly=false;
    
    ReentrantLock lock = new ReentrantLock();
    
    int iterMaxLocalization=20;
            
    int nbGlobalMask=100;
    int computedGlobalMask=0;
       
    
    double [][][] globalmaskA ;
    double [][][] globalmaskA2 ;
    
    double [][][] globalmaskB ;
    double [][][] globalmaskB2 ;
    
    int nbStream;
    int nbThread;
    
    int numberFrameAfterMerging=0;
    
    int numberFrame;
    
    Localization_ [][] loc;
    
    int sizePatch;
    
    
    ImageStack ims1;
    ImageStack ims2;
            
        
    StackLocalization sl1;
    StackLocalization sl2;
    
    StackLocalization sl1fuse;
    StackLocalization sl2fuse;
    
    PolynomialFit pf;//registration parameter
    
    
    
    String pathResFusion;
    
    double maxDistanceMergingXY;
    
    
    
    int idLoc=0;
    
    double sizePix;//nm
    
    
    
    
    DataPhase_ [][] dp;//first dim: stream ; second dim:camera number
    DataPhase_ dp1;
    DataPhase_ dp2;
    int width1;
    int height1;
    int width2;
    int height2;
    
    
    double minZ=Double.MAX_VALUE;//in um
    double maxZ=Double.NEGATIVE_INFINITY;//in um
    
    String path_localization=null;
    
    double rescale_slope1;
    double rescale_intercept1;
    double rescale_slope2;
    double rescale_intercept2;
            
    public DualCamLocalizationPipeline_(String imageCam1,StackLocalization sl1,DataPhase_ dp1,String imageCam2,StackLocalization sl2,DataPhase_ dp2,String pathResFusion,double maxDistanceMergingXY,int nbStream, int nbThread,double adu, double gain1,double offset,double gain2){
        
        this.rescale_slope1=adu/gain1;
        this.rescale_intercept1=-offset*adu/gain1;
        this.rescale_slope2=adu/gain2;
        this.rescale_intercept2=-offset*adu/gain2;
        
        
        this.path_localization=pathResFusion;
        
        
        IJ.selectWindow(imageCam1);
        ims1 = (IJ.getImage()).getImageStack();
        
        IJ.selectWindow(imageCam2);
        ims2 = (IJ.getImage()).getImageStack();
        
        numberFrame=ims1.getSize();
        
        
        this.nbStream=nbStream;
        this.nbThread=nbThread;
        
        
        this.dp1=dp1;
        this.dp2=dp2;
        
        this.sizePatch=dp1.param.sizeoutput;
        this.sl1=sl1;
        this.sl2=sl2;
        
        
        
        
        sizePix=dp1.param.xystep*1000;//convert in nm
        
        this.maxDistanceMergingXY=maxDistanceMergingXY;
        
        this.pathResFusion=pathResFusion;
        
        
        
        
        loc=new Localization_[nbStream][nbThread];
        
        otherVariableName.add("intensity2");
        otherVariableName.add("background2");
        
    }
    
    
    
    public StackLocalization run(StackLocalization stacklocinput){
        
        
        long timeBegin=System.currentTimeMillis();
        
        stackloc=stacklocinput;
        
        if (path_localization!=null){
            if (path_localization.length()>2){
                stackloc.saveOnTheFly(path_localization);
                saveOnTheFly=true;
            }
            else{
                IJ.log("Localization file not saved...you should select a path to save it");
            }
        }
        else{
            IJ.log("oops : problem, no path to save the table -> no result");
        }
        
        
        merge();
        
        numberFrame=sl1fuse.fl.size();
        nbThread=Math.min((numberFrame), nbThread);
        nbStream=Math.min(numberFrame/nbThread,nbStream);//manage wrong initialization compared to slice number
        
        
        double axialRange=2;
        IJ.log("WARNING: axial range fixed to 2 µm");
        SearchPSFcenter_ spc1 = new SearchPSFcenter_(dp1,axialRange);
        dp1.param.ZfocusCenter=spc1.getPosition();
        
        SearchPSFcenter_ spc2 = new SearchPSFcenter_(dp2,axialRange);
        dp2.param.ZfocusCenter=spc2.getPosition();
            
        
        
        
        this.dp=new DataPhase_[nbStream][2];
        for (int u=0;u<nbStream;u++){
            dp[u][0]=new DataPhase_(dp1,u*2);
            dp[u][1]=new DataPhase_(dp2,u*2+1);
            dp[u][0].setMany(this.nbThread);
            dp[u][1].setMany(this.nbThread);
        }
        
        
        
        
        IJ.log("dp search "+dp[0][0].param.Zfocus+"  "+dp[0][0].param.ZfocusCenter);
        IJ.log("dp search "+dp[0][1].param.Zfocus+"  "+dp[0][1].param.ZfocusCenter);
        
        
        
        
        
        
        if (register()){
            
            
            
            
            width1=ims1.getWidth();
            height1=ims1.getHeight();
            width2=ims2.getWidth();
            height2=ims2.getHeight();



            computedGlobalMask=0;
            nbGlobalMask=Math.min(numberFrame, nbGlobalMask);

            globalmaskA = new double [nbGlobalMask][width1][height1];
            globalmaskA2 = new double [nbGlobalMask][width1][height1];
            globalmaskB = new double [nbGlobalMask][width2][height2];
            globalmaskB2 = new double [nbGlobalMask][width2][height2];




            for (int up=0;up<nbStream;up++){
                for (int u=0;u<nbThread;u++){
                    loc[up][u]=new Localization_(1,dp[up],iterMaxLocalization,minZ-.5, maxZ+.5);
                }
            }


            LocalizationThread [][] lt = new LocalizationThread[nbStream][nbThread];


            for (int ip=0;ip<lt.length;ip++){
                for (int i=0;i<lt[0].length;i++){
                    lt[ip][i]=new LocalizationThread(ip,i);
                }
            }

            idLoc=0;

            //launch step by step
            for (int ip=0;ip<lt.length;ip++){
                for (int i=0;i<lt[0].length;i++){
                    lt[ip][i].start();
                }
            }




            for (int ip=0;ip<lt.length;ip++){
                for (int i=0;i<lt[0].length;i++){
                    try{
                        lt[ip][i].join();
                    }catch(Exception eeee){System.out.println("join lt impossible");}
                }
            }

            for (int up=0;up<nbStream;up++){
                for (int u=0;u<nbThread;u++){
                    loc[up][u].kill();
                }
            }


            if (saveOnTheFly){
                stackloc.stopSaveOnTheFly();
            }

            long timeEnd=System.currentTimeMillis();

            IJ.showProgress(0);

            IJ.log("localization done");

            IJ.log("elapsed time = "+((double)(timeEnd-timeBegin))/60000.+" min");
            
            for (int u=0;u<nbStream;u++){
                dp[u][0].free();
                dp[u][1].free();
            }
        
            IJ.showStatus("localization finished");
            return stackloc;
        }
        else{
            for (int u=0;u<nbStream;u++){
                dp[u][0].free();
                dp[u][1].free();
                dp1.free();
                dp2.free();
            }
            return null;
        }
        
    }
    
    
    
    
    
    
    
    
    
    public void merge(){
        
        
        int [][] idLoc;
        
        
        
        double maxX=Double.NEGATIVE_INFINITY;
        double minX=Double.MAX_VALUE;
        double maxY=Double.NEGATIVE_INFINITY;
        double minY=Double.MAX_VALUE;
        
    
    
    
        sl1fuse=new StackLocalization();
        sl2fuse=new StackLocalization();
        
        
        
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
        
        
        double x,y,z;
        
        for (int i=0;i<sl1.fl.size();i++){
            for (int j=0;j<sl1.fl.get(i).loc.size();j++){
                x=sl1.fl.get(i).loc.get(j).X;
                y=sl1.fl.get(i).loc.get(j).Y;
                z=sl1.fl.get(i).loc.get(j).Z;
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
                    
                    if (minZ>z/1000.){
                        minZ=z/1000.;
                    }
                    if (maxZ<z/1000.){
                        maxZ=z/1000.;
                    }
                    
                }
            }
        }
        for (int i=0;i<sl2.fl.size();i++){
            for (int j=0;j<sl2.fl.get(i).loc.size();j++){
                x=sl2.fl.get(i).loc.get(j).X;
                y=sl2.fl.get(i).loc.get(j).Y;
                z=sl2.fl.get(i).loc.get(j).Z;
                
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
                    if (minZ>z){
                        minZ=z;
                    }
                    if (maxZ<z/1000.){
                        maxZ=z/1000.;
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
        numberFrameAfterMerging=0;
        for (int i2=0;i2<idFrameCam2.length;i2++){
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
            
            
            
            int shift=(int)Math.ceil(maxDistanceMergingXY/(sizePix));
            
            int posX1,posY1,posX2,posY2;
            double distance,distZ,xx,yy,zz;
            double prevDistance=Double.MAX_VALUE;
            boolean atLeastOneMolFused=false;
            //if frameFound -> id1 for cam 1 corresponds to id2 for cam 2
            if (frameFound){
                FrameLocalization fl1fuse=new FrameLocalization(numframe1);
                FrameLocalization fl2fuse=new FrameLocalization(numframe1);
                
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
                        for (int u=posX2-shift;u<=posX2+shift;u++){
                            for (int v=posY2-shift;v<=posY2+shift;v++){
                                if ((u>=0)&&(v>=0)&&(u<width)&&(v<height)){
                                    if (idLoc[u][v]!=-1){
                                        if (idPartFound==-1){
                                            //new particle


                                            xx=((sl2.fl.get(id2).loc.get(j).X) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).X));
                                            yy=((sl2.fl.get(id2).loc.get(j).Y) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Y));
                                            zz=((sl2.fl.get(id2).loc.get(j).Z) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Z));
                                            distZ=Math.sqrt(zz*zz);
                                            distance=Math.sqrt(xx*xx+yy*yy);
                                            if ((distance<maxDistanceMergingXY)){
                                                idPartFound=idLoc[u][v];
                                                prevDistance=distance;
                                            }
                                        }
                                        else{
                                            //not new -> compare distance
                                            xx=((sl2.fl.get(id2).loc.get(j).X) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).X));
                                            yy=((sl2.fl.get(id2).loc.get(j).Y) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Y));
                                            zz=((sl2.fl.get(id2).loc.get(j).Z) - (sl1.fl.get(id1).loc.get(idLoc[u][v]).Z));
                                            distZ=Math.sqrt(zz*zz);
                                            distance=Math.sqrt(xx*xx+yy*yy);
                                            if ((distance<maxDistanceMergingXY)&&(distance<prevDistance)){
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
                        if (!atLeastOneMolFused){
                            numberFrameAfterMerging++;
                            atLeastOneMolFused=true;
                        }
                        PLocalization p1=sl1.fl.get(id1).loc.get(idPartFound).copy();
                        p1.id=idLocFuse;
                        PLocalization p2=sl2.fl.get(id2).loc.get(j).copy();
                        p2.id=idLocFuse;
                        fl1fuse.loc.add(p1);
                        fl2fuse.loc.add(p2);
                        
                        
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
        
        
        IJ.log("fuse amount "+fuseCam2+"  /  "+nbLocCam1+" (cam1)  | "+(fuseCam2+notFuseCam2)+" (cam2)");
        
        
        if ((nbNegative/2)>fuseCam2){
            IJ.log("note: there is a surprisingly high number of negative coordinates: "+(nbNegative/2));
        }
        
        
        
    }//end merge
    
    
    
    
    
    public boolean register(){
        
        
        int nbData=0;
        for (int i=0;i<sl1fuse.fl.size();i++){
            nbData+=sl1fuse.fl.get(i).loc.size();
        }
        if (nbData>4){
            double [][] data1 = new double [3][nbData];
            double [][] data2 = new double [3][nbData];
            
            /*double [][] data1tmp = new double [1][nbData];
            double [][] data2tmp = new double [1][nbData];*/
            
            double [][] data1nm = new double [3][nbData];
            double [][] data2nm = new double [3][nbData];
            boolean notApparied=false;

            for (int i=0,t=0;i<sl1fuse.fl.size();i++){
                for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){
                    if (((sl1fuse.fl.get(i).loc.get(j).id!=sl2fuse.fl.get(i).loc.get(j).id))||((sl1fuse.fl.get(i).loc.get(j).frame!=sl2fuse.fl.get(i).loc.get(j).frame))){
                        notApparied=true;
                    }
                    data1[0][t]=sl1fuse.fl.get(i).loc.get(j).X/1000.;
                    data1[1][t]=sl1fuse.fl.get(i).loc.get(j).Y/1000.;
                    data1[2][t]=sl1fuse.fl.get(i).loc.get(j).Z/1000.;
                    data2[0][t]=sl2fuse.fl.get(i).loc.get(j).X/1000.;
                    data2[1][t]=sl2fuse.fl.get(i).loc.get(j).Y/1000.;
                    data2[2][t]=sl2fuse.fl.get(i).loc.get(j).Z/1000.;
                    data1nm[0][t]=sl1fuse.fl.get(i).loc.get(j).X;
                    data1nm[1][t]=sl1fuse.fl.get(i).loc.get(j).Y;
                    data1nm[2][t]=sl1fuse.fl.get(i).loc.get(j).Z;
                    data2nm[0][t]=sl2fuse.fl.get(i).loc.get(j).X;
                    data2nm[1][t]=sl2fuse.fl.get(i).loc.get(j).Y;
                    data2nm[2][t]=sl2fuse.fl.get(i).loc.get(j).Z;
                    
                    /*data1tmp[0][t]=sl1fuse.fl.get(i).loc.get(j).Z;
                    data2tmp[0][t]=sl2fuse.fl.get(i).loc.get(j).Z;*/

                    t++;
                }
            }

            if (notApparied){
                IJ.log("oops : the two files are not apparied -> registration not possible");
            }
            else{
                
                //IJ.log("WARNING WARNING WARNING WARNING WARNING poly fit...ORDER "+0);
                int order=2;
                
                //manual registration in Z
                /*PolynomialFit pftmp = new PolynomialFit(1,data2tmp,data1tmp);
                pftmp.run();
                IJ.log("ZZZZZ");
                pftmp.log();*/
                
                pf = new PolynomialFit(order,data2,data1);
                
                pf.removeParameter(9);//remove Z^2 component
                pf.removeParameter(8);//remove Y.Z component
                pf.removeParameter(7);//remove X.Z component
               
                //WARNING: data1 fits data2 here. it is normal for dual localization next
                //in fact, it will be use to shift the coordinate of cam1 to cam2
                pf.run();
                
                IJ.log("size "+pf.a.length+"  "+pf.a[0].length);
                
                //doing the following: no registration
                /*pf.a[0][0]=0;
                pf.a[0][1]=1;
                pf.a[0][2]=0;
                pf.a[0][3]=0;
                pf.a[0][4]=0;
                pf.a[0][5]=0;
                pf.a[0][6]=0;
                
                pf.a[1][0]=0;
                pf.a[1][1]=0;
                pf.a[1][2]=0;
                pf.a[1][3]=1;
                pf.a[1][4]=0;
                pf.a[1][5]=0;
                pf.a[1][6]=0;
                
                pf.a[2][0]=7.34;
                pf.a[2][1]=0;
                pf.a[2][2]=0;
                pf.a[2][3]=0;
                pf.a[2][4]=0;
                pf.a[2][5]=0;
                pf.a[2][6]=-1;*/
                
                pf.log();
                
                double [] v=new double [3];
                
                v[0]=5;
                v[1]=6;
                v[2]=-.6;
                double [] r=pf.transform(v);
                IJ.log("reg: "+r[0]+"  "+r[1]+"  "+r[2]+"  ");
                v[0]=5;
                v[1]=6;
                v[2]=.15;
                r=pf.transform(v);
                IJ.log("reg: "+r[0]+"  "+r[1]+"  "+r[2]+"  ");
                
                
                for (int i=0,t=0;i<sl1fuse.fl.size();i++){
                    for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){
                        v[0]=sl1fuse.fl.get(i).loc.get(j).X/1000;
                        v[1]=sl1fuse.fl.get(i).loc.get(j).Y/1000;
                        v[2]=sl1fuse.fl.get(i).loc.get(j).Z/1000;
                        double [] res=pf.transform(v);
                        //sl1fuse.fl.get(i).loc.get(j).X=res[0];
                        //sl1fuse.fl.get(i).loc.get(j).Y=res[1];
                        //sl1fuse.fl.get(i).loc.get(j).Z=res[2];

                        t++;
                    }
                }



                PolynomialFit pfinit = new PolynomialFit(order,data1nm,data2nm);
                pfinit.removeParameter(9);//remove Z^2 parameter
                pfinit.removeParameter(8);//remove Y.Z component
                pfinit.removeParameter(7);//remove X.Z component
                
                
                pfinit.run();
                
                //doing the following: no registration
                /*pfinit.a[0][0]=0;
                pfinit.a[0][1]=1;
                pfinit.a[0][2]=0;
                pfinit.a[0][3]=0;
                pfinit.a[0][4]=0;
                pfinit.a[0][5]=0;
                pfinit.a[0][6]=0;
                
                pfinit.a[1][0]=0;
                pfinit.a[1][1]=0;
                pfinit.a[1][2]=0;
                pfinit.a[1][3]=1;
                pfinit.a[1][4]=0;
                pfinit.a[1][5]=0;
                pfinit.a[1][6]=0;
                
                pfinit.a[2][0]=7.34*1000;
                pfinit.a[2][1]=0;
                pfinit.a[2][2]=0;
                pfinit.a[2][3]=0;
                pfinit.a[2][4]=0;
                pfinit.a[2][5]=0;
                pfinit.a[2][6]=-1;*/
                
                //IJ.log("poly fit...2");

                //WARNING: now data2 fits data1. it is normal initialization of localization
                //in fact, the initialization will be the mean of cam2 table and cam1 table -> recorded on cam1 table
                
                
                pfinit.log();
                
                //IJ.log("poly fit...3");

                double x1,y1,z1,x2,y2,z2,cx1,cy1,cz1,cx2,cy2,cz2,cx,cy,cz;
                for (int i=0;i<sl1fuse.fl.size();i++){
                    for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){
                        v[0]=sl2fuse.fl.get(i).loc.get(j).X;
                        v[1]=sl2fuse.fl.get(i).loc.get(j).Y;
                        v[2]=sl2fuse.fl.get(i).loc.get(j).Z;
                        double [] res=pfinit.transform(v);
                        x2=res[0];
                        y2=res[1];
                        z2=res[2];



                        x1=sl1fuse.fl.get(i).loc.get(j).X;
                        y1=sl1fuse.fl.get(i).loc.get(j).Y;
                        z1=sl1fuse.fl.get(i).loc.get(j).Z;


                        if ((sl1fuse.fl.get(i).loc.get(j).crlb_X>=0)&&(sl1fuse.fl.get(i).loc.get(j).crlb_Y>=0)&&(sl1fuse.fl.get(i).loc.get(j).crlb_Z>=0)&&(sl2fuse.fl.get(i).loc.get(j).crlb_X>=0)&&(sl2fuse.fl.get(i).loc.get(j).crlb_Y>=0)&&(sl2fuse.fl.get(i).loc.get(j).crlb_Z>=0)){
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
                    else if ((sl1fuse.fl.get(i).loc.get(j).I>=0)&&(sl2fuse.fl.get(i).loc.get(j).I>=0)){
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
                    
                    
                    }
                }


            }
        }
        else{
            IJ.log("Registration impossible: not enough localizations");
            return false;
        }
        return true;
        
        
        
    }
    
    
    
    
    
 




    
    class Killer extends Thread{
        int time_ms;
        Killer(int time_ms){
            this.time_ms=time_ms;
        }
        public void run(){
            try{
                int nb=time_ms/10000;
                while (nb>0){
                    Thread.sleep(10000);
                    IJ.log("killer ... "+nb);
                    nb--;
                }
                
                System.exit(10);
            }catch(Exception e){}
        }
    }
    
    
    
    class LocalizationThread extends Thread{
    
        double [][] patch;
        int idStream, idThread;

        LocalizationThread(int idStream, int idThread){
            this.idStream=idStream;
            this.idThread=idThread;
            patch = new double [2][sizePatch*sizePatch];
        }

        public void run(){


            //IJ.log("id "+idStream+" / "+nbStream+"      "+idThread+" / "+nbThread);




           int width1=ims1.getWidth();
           int height1=ims1.getHeight();
           int width2=ims2.getWidth();
           int height2=ims2.getHeight();
           int depth1=ims1.size();
           int depth2=ims1.size();
           
           

           if (depth1!=depth2){
               IJ.log("WARNING... are you sure the 2 images are compatibles ? the number of frames is different !");
           }

           double x1,x2,y1,y2,z1,z2,dx1,dy1,dz1,dx2,dy2,dz2;
           double [][] a = new double [2][1];
           double [] b = new double [2];
           int xint1,xint2,yint1,yint2,xbeg1,xbeg2,ybeg1,ybeg2;
           int f;

           //double [][][] frame= new double [depth1][width1][height1];

           //double [][][] frame2= new double [depth2][width2][height2];



           ImageProcessor ip1=null;
           ImageProcessor ip2=null;
           
           
           for (int i=idStream*nbThread+idThread,t=0;i<sl1fuse.fl.size();i+=nbStream*nbThread){//Here, each stream process different images
               
               if (i%(nbStream*nbThread)==0){
                   IJ.showProgress(((float)i)/(float)(sl1fuse.fl.size()));
               }
               int theframe=-1;
               FrameLocalization fl=null;
               IJ.log("process image "+i);
               //for (int j=0;j<Math.min(sl1fuse.fl.get(i).loc.size(), 2);j++){IJ.log("WARNING: DualCamLoc: 2 loc per frame only");
                for (int j=0;j<sl1fuse.fl.get(i).loc.size();j++){

                   x1=sl1fuse.fl.get(i).loc.get(j).X;
                   
                   y1=sl1fuse.fl.get(i).loc.get(j).Y;
                   z1=sl1fuse.fl.get(i).loc.get(j).Z;
                   a[0][0]=sl1fuse.fl.get(i).loc.get(j).I;
                   b[0]=sl1fuse.fl.get(i).loc.get(j).BG;
                   dx1=sl1fuse.fl.get(i).loc.get(j).drift_X;
                   dy1=sl1fuse.fl.get(i).loc.get(j).drift_Y;
                   dz1=sl1fuse.fl.get(i).loc.get(j).drift_Z;
                   
                   theframe=sl1fuse.fl.get(i).loc.get(j).frame;
                   
                   if (j==0){
                       ip1=ims1.getProcessor(theframe+1);
                       ip2=ims2.getProcessor(theframe+1);
                       fl=new FrameLocalization(theframe);
                   }

                   x2=sl2fuse.fl.get(i).loc.get(j).X;
                   double xsave=x2;
                   y2=sl2fuse.fl.get(i).loc.get(j).Y;
                   z2=sl2fuse.fl.get(i).loc.get(j).Z;
                   a[1][0]=sl2fuse.fl.get(i).loc.get(j).I;
                   b[1]=sl2fuse.fl.get(i).loc.get(j).BG;
                   dx2=sl2fuse.fl.get(i).loc.get(j).drift_X;
                   dy2=sl2fuse.fl.get(i).loc.get(j).drift_Y;
                   dz2=sl2fuse.fl.get(i).loc.get(j).drift_Z;
                   
                   
                   x1=((x1+dx1)/1000.)/dp[idStream][0].param.xystep;
                   y1=((y1+dy1)/1000.)/dp[idStream][0].param.xystep;
                   z1=(z1+dz1)/1000.;
                           
                   x2=((x2+dx2)/1000.)/dp[idStream][1].param.xystep;
                   y2=((y2+dy2)/1000.)/dp[idStream][1].param.xystep;
                   z2=(z2+dz2)/1000.;
                   
                   xint1 = (int)(x1);
                   yint1 = (int)(y1);

                   xint2 = (int)(x2);
                   yint2 = (int)(y2);

                   x1-=xint1;
                   y1-=yint1;

                   x2-=xint2;
                   y2-=yint2;
                   
                   x1*=dp[idStream][0].param.xystep;
                   y1*=dp[idStream][0].param.xystep;

                   x2*=dp[idStream][1].param.xystep;
                   y2*=dp[idStream][1].param.xystep;

                   xbeg1=xint1-sizePatch/2;
                   xbeg2=xint2-sizePatch/2;

                   ybeg1=yint1-sizePatch/2;
                   ybeg2=yint2-sizePatch/2;
                   
                   dx1/=1000.;
                   dy1/=1000.;
                   dz1/=1000.;
                   
                   dx2/=1000.;
                   dy2/=1000.;
                   dz2/=1000.;
                   
                   

                   if ((xbeg1>=0)&&(xbeg2>=0)&&(ybeg1>=0)&&(ybeg2>=0)&&(xbeg1+sizePatch<width1)&&(xbeg2+sizePatch<width2)&&(ybeg1+sizePatch<height1)&&(ybeg2+sizePatch<height2)){
                       //IJ.log("ok "+id);
                       
                       for (int u=0;u<sizePatch;u++){
                           for (int uu=0;uu<sizePatch;uu++){
                               patch[0][u*sizePatch+uu]=ip1.getPixelValue(xbeg1+u, ybeg1+uu)*rescale_slope1+rescale_intercept1;
                               
                               patch[1][u*sizePatch+uu]=ip2.getPixelValue(xbeg2+u, ybeg2+uu)*rescale_slope2+rescale_intercept2;
                               
                               
                           }
                        }


                        loc[idStream][idThread].setSubWindow(patch);
                        
                        
                        
                        loc[idStream][idThread].setRegistrationParameters(xint1*dp[idStream][1].param.xystep,yint1*dp[idStream][1].param.xystep,xbeg2*dp[idStream][1].param.xystep-xbeg1*dp[idStream][0].param.xystep,ybeg2*dp[idStream][1].param.xystep-ybeg1*dp[idStream][0].param.xystep,dx1,dy1,dz1,dx2,dy2,dz2,pf,xsave);

                        loc[idStream][idThread].init(a, b, -x1, -y1, z1);
                        
                        loc[idStream][idThread].setTmpCam2(-x2, -y2, z2);
                        
                        
                        boolean locate=loc[idStream][idThread].localize();       
                        
                        
                        
                        if (locate){
                            
                            
                            
                            if ((Math.abs(loc[idStream][idThread].getX())<5*dp[idStream][0].param.xystep)&&(Math.abs(loc[idStream][idThread].getY())<5*dp[idStream][0].param.xystep)){
                                int px_pix=(int)(xint1-(loc[idStream][idThread].getX()/dp[idStream][0].param.xystep));
                                int py_pix=(int)(yint1-(loc[idStream][idThread].getY()/dp[idStream][0].param.xystep));


                                if ((px_pix>=0)&&(px_pix<width1)&&(py_pix>=0)&&(py_pix<height1)){


                                    if ((loc[idStream][idThread].getZ()>minZ)&&(loc[idStream][idThread].getZ()<maxZ)){
                                        
                                        double [] modelA=loc[idStream][idThread].getPSF(0);
                                        double [] modelB=loc[idStream][idThread].getPSF(1);
                                        
                                        
                                        
                                        
                                        double my_x=1000*(xint1*dp[idStream][0].param.xystep-loc[idStream][idThread].getX());
                                        double my_y=1000*(yint1*dp[idStream][0].param.xystep-loc[idStream][idThread].getY());
                                        double my_z=(1000*loc[idStream][idThread].getZ());
                                        double my_A=loc[idStream][idThread].getA()[0][0];
                                        double my_A2=loc[idStream][idThread].getA()[1][0];
                                        double my_B=loc[idStream][idThread].getB()[0];
                                        double my_B2=loc[idStream][idThread].getB()[1];
                                        double my_Score=scoreCompute(patch[0],modelA,patch[1],modelB);
                                        double my_crlbx=(1000*loc[idStream][idThread].getCRLBX());
                                        double my_crlby=(1000*loc[idStream][idThread].getCRLBY());
                                        double my_crlbz=(1000*loc[idStream][idThread].getCRLBZ());
                                        
                                        //IJ.log("crlb "+(sl1fuse.fl.get(i).loc.get(j).crlb_X)+"  "+(sl2fuse.fl.get(i).loc.get(j).crlb_X)+"  "+(1000*loc[idStream][idThread].getCRLBX()));
                                        
                                        //IJ.log("position:   X:"+my_x+"    Y:"+my_y+"    Z:"+my_z+"    A:"+my_A+"    B:"+my_B);
                                        
                                        PLocalization p = new PLocalization(idLoc,theframe,my_x,my_y,my_z,my_A,my_B,my_Score,my_crlbx,my_crlby,my_crlbz);
                                        p.addListOfVariables(otherVariableName);
                                        p.setValueOtherVariable(0, my_A2);
                                        p.setValueOtherVariable(1, my_B2);
                                        p.setDrift_X(dx1*1000);
                                        p.setDrift_Y(dy1*1000);
                                        p.setDrift_Z(dz1*1000);
                                        
                                        if (fl!=null){//not necessary to test
                                            fl.loc.add(p); 
                                        }



                                        //replacement in mask of the region corresponting to 20% of max PSF
                                        try{
                                        lock.lock();
                                            idLoc++;
                                        }
                                        finally{
                                            lock.unlock();
                                        }
                                        
                                        for (int ij=0,ji=-sizePatch/2;ij<sizePatch;ij++,ji++){
                                            for (int iij=0,jji=-sizePatch/2;iij<sizePatch;iij++,jji++){
                                                if (i<nbGlobalMask){
                                                    try{
                                                    if ((xint1+ji>=0)&&(xint1+ji<globalmaskA[i].length)&&(yint1+jji>=0)&&(yint1+jji<globalmaskA[i][0].length)){
                                                        
                                                        globalmaskA[i][xint1+ji][yint1+jji]=modelA[ij*sizePatch+iij];
                                                        
                                                        globalmaskA2[i][xint1+ji][yint1+jji]=patch[0][ij*sizePatch+iij];
                                                    
                                                        //globalmaskA[i][xint1][yint1]=loc[idStream][idThread].getZ();
                                                        //WARNING_____________________________________________________________________WARNING... put xybeg2 to cam [1]
                                                        //globalmaskB[i][xint1+ji][yint1+jji]=modelB[ij*sizePatch+iij];
                                                        //globalmaskB2[i][xint1+ji][yint1+jji]=patch[1][ij*sizePatch+iij];
                                                    }
                                                    if ((xint2+ji>=0)&&(xint2+ji<globalmaskB[i].length)&&(yint2+jji>=0)&&(yint2+jji<globalmaskB[i][0].length)){
                                                        globalmaskB[i][xint2+ji][yint2+jji]=modelB[ij*sizePatch+iij];
                                                        
                                                        globalmaskB2[i][xint2+ji][yint2+jji]=patch[1][ij*sizePatch+iij];
                                                    }
                                                    }
                                                    catch(Exception ee){
                                                        
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                }
                                
                            }
                            
                            
                        }
                        
                        
                        
                   }



               }

           
                 try{
                     lock.lock();
                     stackloc.fl.add(fl);
                     if ((sl1fuse.fl.size()-i-1)<nbStream*nbThread){
                         //IJ.log("decrement "+idStream);
                             dp[idStream][0].modelMany.decrementNumberPSFToCompute();
                             dp[idStream][1].modelMany.decrementNumberPSFToCompute();

                     }
                     else{
                         //IJ.log("non decrement "+idStream);
                     }


                     if (i<nbGlobalMask){
                         //maybe starts lock here
                         
                         //IJ.log("comp "+computedGlobalMask);
                         if ( computedGlobalMask == nbGlobalMask-1 ){
                             ImageShow.imshow(globalmaskA,"model_frame1-100");
                             ImageShow.imshow(globalmaskA2,"image_frame1-100");
                             ImageShow.imshow(globalmaskB,"model_frame2-100");
                             ImageShow.imshow(globalmaskB2,"image_frame2-100");
                             globalmaskA=null;
                             globalmaskA2=null;
                             globalmaskB=null;
                             globalmaskB2=null;
                         }
                         computedGlobalMask++;
                     }
                 }
                 finally{
                     lock.unlock();
                 }
           }
           //ImageShow.imshow(frame,"frame");
           //ImageShow.imshow(frame2,"frame2");

       }
   
   
   
   
   
    }
    
    public double scoreCompute(double [] patch1, double [] model1,double [] patch2, double [] model2){
        
        
        
        double sum=0;
        double sumModel=0;
        for (int i=0;i<patch1.length;i++){
            //if (model[i]-background>maxModelThreshold){
                sum+=(patch1[i]-model1[i])*(patch1[i]-model1[i])/model1[i];
                sumModel++;
            //}

        }
        for (int i=0;i<patch2.length;i++){
            //if (model[i]-background>maxModelThreshold){
                sum+=(patch2[i]-model2[i])*(patch2[i]-model2[i])/model2[i];
                sumModel++;
            //}

        }
        sum/=sumModel;
        return sum;
        
        
    }
   
   
}

