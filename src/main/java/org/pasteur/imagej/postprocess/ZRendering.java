/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;

/**
 *
 * @author benoit
 */

import org.pasteur.imagej.data.*;
import org.pasteur.imagej.utils.ImageShow;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.ImagePlus;
import ij.ImageStack;
import ij.IJ;
import ij.process.ColorProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import ij.LookUpTable;

import ij.plugin.LutLoader;
import ij.process.LUT;

/**
 * Simple rendering using scatter plot. If there is any molecule at the pixel
 * location, the pixel value will be a constant positive number and zero
 * otherwise.
 */
public class ZRendering {
    
    public static String nameHistPlot="histogram";
    public static String nameScatterPlot="scatter plot";
    
    private static String pathLutFolder=IJ.getDirectory("startup")+"luts"+File.separator;
    public static String lut="ZOLA.lut";
    private static String lut2="ZOLA2.lut";
    
    
    public static String [] getFijiLUTs(){
        
        String [] liste=null;
        try{
            
            File file = new File(pathLutFolder);
            if (!file.exists()){
                IJ.log("The 'luts' folder not found in the ImageJ folder");
                return null;
            }
            liste = file.list();
            int count=0;
            for (int i=0; i<liste.length; i++) {
                if (liste[i].endsWith(".lut")) {
                    try{
                        LUT lut=LutLoader.openLut(pathLutFolder+liste[i]);
                        count++;
                    }catch(Exception e ){}
                    
                }
            }
            return liste;
        }
        catch(Exception e2){}
        
        return null;
    }
    
    
    
    
    
    
    public static void makeLUT(){
        String pathlutfolder=IJ.getDirectory("startup")+""+"luts";
        File ff = new File(pathlutfolder);
        if (!ff.exists()){
            ff.mkdir();
        }
        String pathlut=IJ.getDirectory("startup")+""+"luts/"+lut;
        String pathlut2=IJ.getDirectory("startup")+""+"luts/"+lut2;
        {
            File f = new File(pathlut);
            String s="0\t0\t0\t0\n";
            if (!f.exists()){
                for (int i=1;i<255;i++){
                    double angle=1.-((double)i)/256.;
                    double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                    if (H<0){
                        H+=360;
                    }
                    double truc=0;
                    double S = 1;
                    double V =1;
                    int [] RGB=getRGB(H,S,V);
                    //IJ.log("H "+H);
                    s+=i+"\t"+RGB[0]+"\t"+RGB[1]+"\t"+RGB[2]+"\n";

                }
                s+="255\t0\t0\t0\n";
                try{
                    PrintWriter sortie;

                    sortie = new PrintWriter(new FileWriter(pathlut, false));

                    sortie.println(""+s);
                    sortie.close();
                }catch(Exception e){
                    IJ.log("lookup table impossible to create. Can we create files in ImageJ ?");
                }
            }
        }
        /*{
            File f2 = new File(pathlut2);

            String s2="";
            if (!f2.exists()){
                for (int i=0;i<=255;i++){
                    double angle=1.-((double)i)/256.;
                    double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                    if (H<0){
                        H+=360;
                    }
                    double truc=0;
                    double S = 1;
                    double V =1;
                    int [] RGB=getRGB(H,S,V);
                    //IJ.log("H "+H);
                    s2+=i+"\t"+RGB[0]+"\t"+RGB[1]+"\t"+RGB[2]+"\n";

                }
                try{
                    PrintWriter sortie;

                    sortie = new PrintWriter(new FileWriter(pathlut2, false));
                    sortie.println(""+s2);
                    sortie.close();
                }catch(Exception e){
                    IJ.log("lookup table impossible to create. Can we create files in ImageJ ?");
                }
            }
        }*/
    }
    
    
    
    
    static int [] getRGB(double H, double S, double V){

        int ti=((int)H/60)%6;
        double fi=H/60.-(double)ti;
        double l=V*(1-S);
        double m=V*(1-fi*S);
        double n=V*(1-(1-fi)*S);
        double Rp=0,Gp=0,Bp=0;
        if (ti==0){
            Rp=V;Gp=n;Bp=l;
        }
        if (ti==1){
            Rp=m;Gp=V;Bp=l;
        }
        if (ti==2){
            Rp=l;Gp=V;Bp=n;
        }
        if (ti==3){
            Rp=l;Gp=m;Bp=V;
        }
        if (ti==4){
            Rp=n;Gp=l;Bp=V;
        }
        if (ti==5){
            Rp=V;Gp=l;Bp=m;
        }

        int [] RGB=new int [3];
        RGB[0]=(int)((Rp)*255);
        RGB[1]=(int)((Gp)*255);
        RGB[2]=(int)((Bp)*255);
        return RGB;
    }
    
    
    
    
    
    public static ImagePlus scatter2D(StackLocalization sl,double pixelsizeNM,boolean color) {
        ZRendering.makeLUT();
        double binColor=256;
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        
        
        return scatter2D(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,color);
        
    }
    
    
    
    
    
    public static ImagePlus scatter2D(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ,boolean color) {
        ZRendering.makeLUT();
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        ZRendering.makeLUT();
        double binColor=256;
        
        
        double x;
        double y;
        double z;
        
        
        
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        
        
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        
        
        
        ImageProcessor ip = new FloatProcessor(width,height);
        for (int u=0;u<width;u++){
//            for (int uu=0;uu<height;uu++){
//                ip.putPixelValue(u, uu, 0);
//            }
        }
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (z<minZ){
                        z=minZ;
                    }
                    if (z>maxZ){
                        z=maxZ;
                    }
                    x=Math.min(Math.max((int)((x-minX)/pixelsizeNM),0),width-1);
                    y=Math.min(Math.max((int)((y-minY)/pixelsizeNM),0),height-1);
                    if ((x>=0)&&(x<width)&&(y>=0)&&(y<height)){
                        if (Math.abs(maxZ-minZ)<.0000000000001){
                            ip.putPixelValue((int)x, (int)y, 1);
                        }
                        else{
                            if (color){
                                ip.putPixelValue((int)x, (int)y, z-minZ+((maxZ-minZ)/binColor));//((maxZ-minZ)/binColor) just to have non black bottom
                            }
                            else{
                                ip.putPixelValue((int)x, (int)y, 1);
                            }
                        }
                    }
                    
                }
            }
        }
        ImagePlus imp = new ImagePlus("2D scatter plot ",ip);
        imp.show();
        if (color){
            try{
                Thread.sleep(200);
                IJ.run("ZOLA");
            }
            catch(Exception e){}
        }
        IJ.setMinAndMax(imp, 0, maxZ-minZ+(1./binColor));
        return imp;
    }
    
    
    
    
    
    
    
    //scatter 3D plot according to selected variable (by default Z position)
    public static ImagePlus scatter3D(StackLocalization sl,double pixelsizeNM,boolean color) {
        ZRendering.makeLUT();
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        
        return scatter3D(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,-1);
        
    }
    
    
    
    
    
    public static ImagePlus scatter3D(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ) {
        return scatter3D(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,-1);
    }
    
    public static ImagePlus scatter3D(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ,int idVariable) {
        ZRendering.makeLUT();
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        
        if (idVariable<0){
            idVariable=4;//Z by default
        }
        
        double minV=Double.POSITIVE_INFINITY;
        double maxV=Double.NEGATIVE_INFINITY;
        
        
        double v;
        
        //get min and max of selected variable
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    v=sl.fl.get(i).loc.get(j).getValueVariable(idVariable);
                    
                    if (v<minV){
                        minV=v;
                    }
                    if (v>maxV){
                        maxV=v;
                    }
                }
            }
        }
        
        ZRendering.makeLUT();
        double binColor=255;
        
        double x;
        double y;
        double z;
        
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        int depth=(int)Math.ceil((maxZ-minZ)/pixelsizeNM);
        
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        depth=Math.max(depth, 1);
        
        
        ImageStack ims = new ImageStack(width,height);
        
        ImageProcessor [] ip = new FloatProcessor[depth];
        for (int t=0;t<depth;t++){
            ip[t]= new FloatProcessor(width,height);
//            for (int u=0;u<width;u++){
//                for (int uu=0;uu<height;uu++){
//                    ip[t].putPixelValue(u, uu, 0);
//                }
//            }
        }
        int zz;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    v=sl.fl.get(i).loc.get(j).getValueVariable(idVariable);
                    if (z<minZ){
                        z=minZ;
                    }
                    if (z>maxZ){
                        z=maxZ;
                    }
                    x=Math.min(Math.max((int)((x-minX)/pixelsizeNM),0),width-1);
                    y=Math.min(Math.max((int)((y-minY)/pixelsizeNM),0),height-1);
                    //z=Math.min(Math.max((int)((z-minZ)/pixelsizeNM),0),depth-1);
                    int zzz=(int)Math.min(Math.max((int)((z-minZ)/pixelsizeNM),0),depth-1);
                    if ((zzz>=0)&&(zzz<ip.length)&&(x>=0)&&(x<width)&&(y>=0)&&(y<height)){
                        ip[zzz].putPixelValue((int)x, (int)y, (v-minV+1));
                            
                    }
                    
                }
            }
        }
        for (int t=0;t<depth;t++){
            ims.addSlice("z="+(((double)t)*pixelsizeNM), ip[t]);
        }
        ImagePlus imp = new ImagePlus(""+nameScatterPlot+" ",ims);
        imp.show();
        try{
            Thread.sleep(200);
            IJ.run("ZOLA");
        }
        catch(Exception e){}
        IJ.setMinAndMax(imp, 0, maxV-minV+1+((maxZ-minZ)/binColor));
        return imp;
    }

    
    
    
    
    
    
    
    
    
        
    
        
    /*public static ImagePlus smartColorRendering(ImagePlus imp,StackLocalization sl,double pixelsizeNM,int shift,boolean printLUT)  {
        
        
        try{
            ZRendering.makeLUT();


            double minX=Double.POSITIVE_INFINITY;
            double maxX=Double.NEGATIVE_INFINITY;

            double minY=Double.POSITIVE_INFINITY;
            double maxY=Double.NEGATIVE_INFINITY;

            double minZ=Double.POSITIVE_INFINITY;
            double maxZ=Double.NEGATIVE_INFINITY;

            double x;
            double y;
            double z;

            //maybe use Arrays.sort to be fast
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        x=sl.fl.get(i).loc.get(j).X;
                        y=sl.fl.get(i).loc.get(j).Y;
                        z=sl.fl.get(i).loc.get(j).Z;
                        if (x<minX){
                            minX=x;
                        }
                        if (y<minY){
                            minY=y;
                        }
                        if (x>maxX){
                            maxX=x;
                        }
                        if (y>maxY){
                            maxY=y;
                        }
                        if (z<minZ){
                            minZ=z;
                        }
                        if (z>maxZ){
                            maxZ=z;
                        }
                    }
                }
            }
            if (maxZ-minZ<pixelsizeNM){
                maxZ=minZ+pixelsizeNM;
            }
            if (maxY-minY<pixelsizeNM){
                maxY=minY+pixelsizeNM;
            }
            if (maxX-minX<pixelsizeNM){
                maxX=minX+pixelsizeNM;
            }
            double crlbXTh=0;
            double crlbYTh=0;
            double crlbZTh=0;
            double khiTh=Double.POSITIVE_INFINITY;
            double meancrlbX=0;
            double meancrlbY=0;
            double meancrlbZ=0;
            double stdcrlbX=0;
            double stdcrlbY=0;
            double stdcrlbZ=0;
            double minZsmart=minZ;
            double maxZsmart=maxZ;
            double count=0;
            
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        meancrlbX+=sl.fl.get(i).loc.get(j).crlb_X;
                        meancrlbY+=sl.fl.get(i).loc.get(j).crlb_Y;
                        meancrlbZ+=sl.fl.get(i).loc.get(j).crlb_Z;
                        count++;
                    }
                }
            }
            if (count==0){
                return null;
            }
            meancrlbX/=count;
            meancrlbY/=count;
            meancrlbZ/=count;
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        stdcrlbX+=(sl.fl.get(i).loc.get(j).crlb_X-meancrlbX)*(sl.fl.get(i).loc.get(j).crlb_X-meancrlbX);
                        stdcrlbY+=(sl.fl.get(i).loc.get(j).crlb_Y-meancrlbY)*(sl.fl.get(i).loc.get(j).crlb_Y-meancrlbY);
                        stdcrlbZ+=(sl.fl.get(i).loc.get(j).crlb_Z-meancrlbZ)*(sl.fl.get(i).loc.get(j).crlb_Z-meancrlbZ);
                    }
                }
            }
            stdcrlbX/=count;
            stdcrlbY/=count;
            stdcrlbZ/=count;

            crlbXTh=meancrlbX+3*Math.sqrt(stdcrlbX);
            crlbYTh=meancrlbY+3*Math.sqrt(stdcrlbY);
            crlbZTh=meancrlbZ+3*Math.sqrt(stdcrlbZ);
            
            
            int bin=(int)Math.ceil(((maxZ-minZ)/pixelsizeNM));
            
            int [] hist = new int[bin];
            for (int i=0;i<bin;i++){
                hist[i]=0;
            }
            for (int i=0;i<sl.fl.size();i++){
                for (int j=0;j<sl.fl.get(i).loc.size();j++){
                    if (sl.fl.get(i).loc.get(j).exists){
                        PLocalization p =sl.fl.get(i).loc.get(j);
                        if ((crlbXTh>0)&&(crlbXTh>0)&&(crlbXTh>0)){
                            if (((p.crlb_X<crlbXTh)||(p.crlb_X<0))&&((p.crlb_Y<crlbYTh)||(p.crlb_Y<0))&&((p.crlb_Z<crlbZTh)||(p.crlb_Z<0))&&((p.score<khiTh)||(p.score<0))){
                                double zi=sl.fl.get(i).loc.get(j).Z;
                                if ((zi>=minZ)&&(zi<=maxZ)){
                                    int b=(int)((double)bin*(zi-minZ)/(maxZ-minZ));
                                    if ((b>=0)&&(b<bin)){
                                     hist[b]++;
                                    }
                                }
                            }
                        }
                        else{
                            double zi=sl.fl.get(i).loc.get(j).Z;
                            if ((zi>=minZ)&&(zi<=maxZ)){
                                int b=(int)((double)bin*(zi-minZ)/(maxZ-minZ));
                                if ((b>=0)&&(b<bin)){
                                 hist[b]++;
                                }
                            }
                        }
                    }
                }
            }
            //threshold hist according to max 10%
            int maxhist=hist[0];
            for (int i=1;i<bin;i++){
                if (hist[i]>maxhist){
                    maxhist=hist[i];
                }
            }
            for (int i=0;i<bin;i++){

                if (hist[i]>.1*(double)maxhist){
                    hist[i]=1;
                }
                else{
                    hist[i]=0;
                }
            }
            ArrayList<Integer> alstart = new ArrayList();
            ArrayList<Integer> alstop = new ArrayList();
            ArrayList<Integer> alsize = new ArrayList();
            int index=-1;
            if (hist[0]==1){
                alstart.add(0);
                alsize.add(1);
                index++;
            }
            for (int i=1;i<bin;i++){
                if (hist[i]==1){
                    if (hist[i-1]==0){
                        alstart.add(i);
                        alsize.add(1);
                        index++;
                    }
                    else{
                        alsize.set(index,((int)alsize.get(index)+1));
                    }
                }
                else{
                    if (hist[i-1]==1){
                        alstop.add(i);
                    }
                }
            }
            if (hist[bin-1]==1){
                alstop.add(bin);
            }

            //IJ.log(""+alstart.size()+"  "+alsize.size()+"  "+alstop.size());


            maxhist=0;
            for (int i=1;i<alsize.size();i++){
                if ((int)alsize.get(i)>(int)alsize.get(maxhist)){
                    maxhist=i;

                }
            }
            minZsmart=((((double)((int)alstart.get(maxhist)))/((double)bin))*(maxZ-minZ))+minZ;
            maxZsmart=((((double)((int)alstop.get(maxhist)))/((double)bin))*(maxZ-minZ))+minZ;
            
            
            if (imp==null){
                imp = new ImagePlus();
                imp.setTitle("2D color histogram (automatic filtering) "+pixelsizeNM+"nm per px");
            }

            return colorRendering(imp,sl,pixelsizeNM,minX,maxX,minY,maxY,minZsmart,maxZsmart,shift,crlbXTh,crlbYTh,crlbZTh,khiTh,printLUT);
        }
        catch(Exception problem){
            //???? should not be any problem... just in case
        }
        return imp;
    }*/
    
    
    
    
    
    
    public static ImagePlus colorRendering(ImagePlus imp,StackLocalization sl,double pixelsizeNM,int shift,boolean printLUT) {
        
        ZRendering.makeLUT();
        
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        
        return colorRendering(imp,sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,shift,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,printLUT,ZRendering.lut);
        
    }
    
    public static ImagePlus colorRendering(ImagePlus imp,StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ,double maxZ,int shift,boolean printLUT) {
        
        ZRendering.makeLUT();
        
        
        return colorRendering(imp,sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,shift,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,printLUT,ZRendering.lut);
        
    }
    
    
    public static ImagePlus colorRendering(ImagePlus imp,StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ,double maxZ,int shift,boolean printLUT,String lut) {
        
        ZRendering.makeLUT();
        
        
        return colorRendering(imp,sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,shift,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,printLUT,lut);
        
    }
    
    
    
    
    public static ImagePlus colorRendering(ImagePlus imp,StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZsmart, double maxZsmart,int shift,double crlbXTh,double crlbYTh,double crlbZTh,double khiTh,boolean printLUT,String lut)  {
        
        
        byte [] rlut=null; 
        byte [] glut=null; 
        byte [] blut=null; 
        
        if (lut!=ZRendering.lut){
            try{
                LUT lut_=LutLoader.openLut(pathLutFolder+lut);
                rlut = new byte[256] ;
                glut = new byte[256] ;
                blut = new byte[256] ;
                
                lut_.getReds(rlut);
                lut_.getBlues(blut);
                lut_.getGreens(glut);
                
                
            }catch(Exception e ){lut=ZRendering.lut;}
            
        }
        
        
        if (maxZsmart-minZsmart<pixelsizeNM){
            maxZsmart=minZsmart+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        int depth=1;
        width+=shift*2;
        height+=shift*2;
        int nbShift=1+shift*2;
        int [][] shifT = new int[1+shift*2][1+shift*2];
        
        //create first element
        shifT[0][0]=1;
        //create first line
        for (int i=1;i<=shift;i++){
            shifT[0][i]=shifT[0][i-1]+1;
        }
        for (int i=shift+1;i<1+shift*2;i++){
            shifT[0][i]=shifT[0][i-1]-1;
        }
        //create each column
        for (int ii=0;ii<nbShift;ii++){
            for (int i=1;i<=shift;i++){
                shifT[i][ii]=shifT[i-1][ii]+shifT[0][ii];
            }
        }
        for (int ii=0;ii<nbShift;ii++){
            for (int i=shift+1;i<1+shift*2;i++){
                shifT[i][ii]=shifT[i-1][ii]-shifT[0][ii];
            }
        }
        double maxShift=0;
        double meanShift=0;
        for (int i=0;i<shifT.length;i++){
            for (int ii=0;ii<shifT[i].length;ii++){
                meanShift+=shifT[i][ii];
                if (shifT[i][ii]>maxShift){
                    maxShift=shifT[i][ii];
                }
            }
        }
        meanShift/=(double)(shifT.length*shifT[0].length);
        
        
        double x;
        double y;
        double z;
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        
        ArrayList<Double> [][] im2Dz = new ArrayList[width][height];
        ArrayList<Integer> [][] im2Dw = new ArrayList[width][height];
        int xx,yy;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    PLocalization p =sl.fl.get(i).loc.get(j);
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    xx=(int)((x-minX)/pixelsizeNM)+shift;
                    yy=(int)((y-minY)/pixelsizeNM)+shift;
                    if ((x>=minX)&&(x<=maxX)&&(y>=minY)&&(y<=maxY)&&(z>=minZsmart)&&(z<=maxZsmart)&&((p.score<=khiTh)||(p.score<0))&&((p.crlb_X<=crlbXTh)||(p.crlb_X<0))&&((p.crlb_Y<=crlbYTh)||(p.crlb_Y<0))&&((p.crlb_Z<=crlbZTh)||(p.crlb_Z<0))){
                        
                        if ((xx-shift>=0)&&(yy-shift>=0)&&(xx+shift<width)&&(yy+shift<height)){
                            for (int a=-shift,aa=0;a<=shift;a++,aa++){
                                for (int b=-shift,bb=0;b<=shift;b++,bb++){
                                    if ((xx+a>=0)&&(xx+a<width)&&(yy+b>=0)&&(yy+b<height)){
                                        if (im2Dz[xx+a][yy+b]==null){
                                            im2Dz[xx+a][yy+b]=new ArrayList<Double>();
                                            im2Dw[xx+a][yy+b]=new ArrayList<Integer>();
                                            
                                        }
                                        im2Dz[xx+a][yy+b].add(z);
                                        im2Dw[xx+a][yy+b].add(shifT[aa][bb]);
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //compute weighted mean and std
        double [][] meanZ = new double[width][height];
        double [][] stdZ = new double[width][height];
        double [][] weight = new double[width][height];
        double meanWeight=0;
        double nbWeight=0;
        double countm=0;
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                meanZ[i][ii]=0;
                stdZ[i][ii]=0;
                countm=0;
                if (im2Dz[i][ii]!=null){
                    for (int u=0;u<im2Dz[i][ii].size();u++){
                        meanZ[i][ii]+=im2Dz[i][ii].get(u)*(double)im2Dw[i][ii].get(u);
                        countm+=(double)im2Dw[i][ii].get(u);
                    }
                    weight[i][ii]=countm;
                    if (countm!=0){
                        meanWeight+=weight[i][ii];
                        nbWeight++;
                        meanZ[i][ii]/=countm;
                        for (int u=0;u<im2Dz[i][ii].size();u++){
                            stdZ[i][ii]+=(double)im2Dw[i][ii].get(u)*(im2Dz[i][ii].get(u)-meanZ[i][ii])*(im2Dz[i][ii].get(u)-meanZ[i][ii]);
                        }
                        stdZ[i][ii]/=countm;
                        stdZ[i][ii]=Math.sqrt(stdZ[i][ii]);
                        
                    }
                }
                else{
                    weight[i][ii]=countm;
                }
            }
        }
        meanWeight/=nbWeight;
        meanWeight*=maxShift/meanShift;
        ImageProcessor ip = new ColorProcessor(width,height);
        
        
        
        double val;
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                val=weight[i][ii];
                if (val!=0){
                    
                    double angle=1.-((double)meanZ[i][ii]-minZsmart)/(maxZsmart-minZsmart);
                    int [] RGB=null;
                    if (lut.equals(ZRendering.lut)){
                        double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                        if (H<0){
                            H+=360;
                        }
                        double V=0,S=0;

                        /*if (val-minrange<((maxrange-minrange)/2)){
                            S=1;
                            V =((val)*2-minrange*2)/(maxrange-minrange);
                        }
                        else{
                            S=1-((val-((maxrange-minrange)/2))*2-minrange)/(maxrange-minrange);
                            V =1;
                        }*/

                        //color bar cut in 3 maincolor (RGB)/ I put a factor 2 in order to have white color when  std is extremely large (blue/red) 
                        //otherwise, a mix color is kept

                        S=1-2*((double)stdZ[i][ii])/(maxZsmart-minZsmart);

                        V =val/(meanWeight*2);//*2 to avoid saturation
                        //V=1;
                        V=Math.max(Math.min(V,1), 0);
                        S=Math.max(Math.min(S,1), 0);
                        RGB=getRGB(H,S,V);
                    }
                    else{
                        RGB=new int[3];
                        double red=rlut[(int)Math.ceil((1-angle+.00000001)*255)-1];
                        double green=glut[(int)Math.ceil((1-angle+.00000001)*255)-1];
                        double blue=blut[(int)Math.ceil((1-angle+.00000001)*255)-1];
                        if (red<0){
                            red+=256;
                        }
                        if (green<0){
                            green+=256;
                        }
                        if (blue<0){
                            blue+=256;
                        }
                        
                        double V =val/(meanWeight*2);//*2 to avoid saturation
                        //V=1;
                        V=Math.max(Math.min(V,1), 0);
                        
                        RGB[0]=(int)(red*V);
                        RGB[1]=(int)(green*V);
                        RGB[2]=(int)(blue*V);
                        
                    }
                    
                    
                    ip.putPixel(i, ii,RGB);
                }
            }
        }
        
        if (imp==null){
            imp = new ImagePlus("2D color histogram ",ip);
            imp.show();
        }
        else{
            imp.setProcessor(ip);
            imp.updateAndDraw();
            if (!imp.isVisible()){
                imp.show();
            }
        }
        
        
        
        if (printLUT){
            if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro())){
                addLUT(imp,minZsmart, maxZsmart,lut);
            }
        }
        
        
        
        IJ.run("Set Scale...", "distance=1 known="+pixelsizeNM+" unit=nm");
        
        return imp;
    }
    
    
    
    
    
    
    
    
    
    public static ImagePlus timeRendering2D(StackLocalization sl,double pixelsizeNM,int frameSplit) {
        double binColor=256;
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        
        
        return timeRendering2D(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,frameSplit);
        
    }
    
    
    
    
    
    
    
    
    
    
    public static ImagePlus timeRendering2D(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ,int frameSplit) {
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        ZRendering.makeLUT();
        double binColor=256;
        
        
        int maxF=Integer.MIN_VALUE;
        int minF=Integer.MAX_VALUE;
        for (int i=0;i<sl.fl.size();i++){
            int f=sl.fl.get(i).numFrame;
            if (f<minF){
                minF=f;
            }
            if (f>maxF){
                maxF=f;
            }
        }
        int frameNumb=maxF-minF+1;
        if ((frameSplit<1)||(frameSplit>frameNumb)){
            frameSplit=frameNumb;
        }
        int bin=frameNumb/frameSplit;
        
        
        
        
        double saturation=2;
        
        
        
        
        
        
        
        
        double x;
        double y;
        double z;
        
        
        
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        
        
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        
        
        
        
        
        
        
        IJ.log("slice number "+bin);
        
        int [][][][] mat = new int[bin][width][height][];
        int [][][] count = new int[bin][width][height];
        
        for (int b=0;b<bin;b++){
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    count[b][i][ii]=0;
                    mat[b][i][ii]=null;
                }
            }
        }
        
        
        for (int i=0;i<sl.fl.size();i++){
            int f=sl.fl.get(i).numFrame;
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    

                    int b=(f-minF)/frameSplit;
                    
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    int xxx=Math.min(Math.max((int)((x-minX)/pixelsizeNM),0),width-1);
                    int yyy=Math.min(Math.max((int)((y-minY)/pixelsizeNM),0),height-1);
                    
                    if ((b>=0)&&(b<bin)&&(xxx>=0)&&(yyy>=0)&&(xxx<width)&&(yyy<height)){
                        count[b][xxx][yyy]++;
                    }
                }
            }
        }
        
        
        for (int i=0;i<sl.fl.size();i++){
            int f=sl.fl.get(i).numFrame;
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){


                    int b=(f-minF)/frameSplit;
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    int xxx=Math.min(Math.max((int)((x-minX)/pixelsizeNM),0),width-1);
                    int yyy=Math.min(Math.max((int)((y-minY)/pixelsizeNM),0),height-1);
                    
                    if ((b>=0)&&(b<bin)&&(xxx>=0)&&(yyy>=0)&&(xxx<width)&&(yyy<height)){
                        if (mat[b][xxx][yyy]==null){
                            mat[b][xxx][yyy]=new int[count[b][xxx][yyy]];
                            count[b][xxx][yyy]=0;//recount for id position
                        }
                        mat[b][xxx][yyy][count[b][xxx][yyy]]=f;
                        count[b][xxx][yyy]++;
                    }
                }
            }
        }
        
        
        
        
        
        IJ.log("NOTE:");
        
        
        IJ.log("The color defines the frame: purple->blue->cian->green->yellow->red->pink in range: [0-"+frameSplit+"]");
        IJ.log("The color corresponds to the mean of frame position in each pixel");
        IJ.log("The saturation corresponds to the stdev of frame position in each pixel:");
        IJ.log("The higher the stdev, the whiter it is ; range(std):[0-"+(frameSplit/(2*saturation))+"]");
        
        
        
        
        
        
        
        
        ImageStack ims = new ImageStack(width,height);
        
        for (int b=0;b<bin;b++){
            ColorProcessor ip = new ColorProcessor(width,height);
            //FloatProcessor ip = new FloatProcessor(imageWidthM,imageHeightM);
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    if (mat[b][i][ii]!=null){
                        //double mx=0;
                        //double my=0;
                        double std=0;
                        double ma=0;
                        for (int k=0;k<mat[b][i][ii].length;k++){
                            double angle=2.*Math.PI-2.*Math.PI*((double)(((mat[b][i][ii][k]-minF)%frameSplit)))/((double)frameSplit);//2*Pi-x : to reverse color range: from blue to red
                            ma+=angle;
                            
                            
                            //mx+=Math.cos(angle);
                            //my+=Math.sin(angle);
                        }
                        //mx/=(double)mat[b][i][ii].length;
                        //my/=(double)mat[b][i][ii].length;
                        ma/=(double)mat[b][i][ii].length;
                        for (int k=0;k<mat[b][i][ii].length;k++){
                            double angle=2.*Math.PI-2*Math.PI*((double)((mat[b][i][ii][k]-minF)%frameSplit))/((double)frameSplit);
                            std+=(angle-ma)*(angle-ma);
                        }
                        std/=(double)mat[b][i][ii].length;
                        std=Math.sqrt(std)/Math.PI;
                        
                        //double anglenew=Math.PI+Math.atan2(my, mx);
                        
                        //double sat=Math.sqrt(mx*mx+my*my);
                        //if ((anglenew<Math.PI/4.)&&(mat[b][i][ii].length>1)){//45 degres -> max saturation
                        //    sat=0;
                        //}
                        double H =-45.+((315.)*ma/(Math.PI*2.)); //convert in degree 315->max color=red     +45 degres: min color blue
                        if (H<0){
                            H+=360;
                        }
                        //double H =360. *ma/(Math.PI*2.); //convert in degree 315->max color=red     +45 degres: min color blue
                        double S = (1-(std*saturation));
                        S=Math.max(S, 0);
                        //double S = 1;
                        double V =1;
                        int [] RGB=getRGB(H,S,V);
                        ip.putPixel(i, ii,RGB);
                        
                        //ip.putPixelValue(i, ii,l);
                        
                    }
                    
                    
                    
                }
            }
            
            //draw calibration bar
            
            int thesizeX=Math.min(width, height)/40;
            int thesizeY=Math.min(width, height)/10;
            
            int sizex=5;
            int sizey=20;
            
            if (thesizeX>sizex){
                sizex=thesizeX;
            }
            if (thesizeY>sizey){
                sizey=thesizeY;
            }
            
            int decalPosit=sizex;
            int x2=width-decalPosit;
            int x1=x2-sizex;
            int y1=decalPosit;
            int y2=y1+sizey;


            if ((y2<height)&&(x1>=0)){
                for (int i=y1,ii=sizey;i<y2;i++,ii--){
                    double angle=1-((double)ii)/(double)sizey;
                    double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                    if (H<0){
                        H+=360;
                    }
                    
                    for (int j=x1,jj=0;j<x2;j++,jj++){
                        double S = 1-(double)jj/(double)sizex;
                        //double S = 1;
                        double V =1;
                        int [] RGB = getRGB(H,S,V);
                        ip.putPixel(j, i,RGB);
                        
                    }
                }

            }
            else{
                IJ.log("no calibration bar: image too small");
            }
            ims.addSlice("frame_"+(int)(b*frameSplit)+"->"+(int)((b+1)*frameSplit),ip);
        }
        
        ImagePlus impp = new ImagePlus("colorStackFrame",ims);
        
        impp.show();
        return impp;
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public static ImagePlus hist2D(StackLocalization sl,double pixelsizeNM,int shift) {
        
        ZRendering.makeLUT();
        
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        return hist2D(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,shift);
        
    }
        
    
        
        
    public static ImagePlus hist2D(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ,int shift)  {
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        double binColor=255;
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        int depth=1;
        width+=shift*2;
        height+=shift*2;
        int nbShift=1+shift*2;
        int [][] shifT = new int[1+shift*2][1+shift*2];
        
        //create first element
        shifT[0][0]=1;
        //create first line
        for (int i=1;i<=shift;i++){
            shifT[0][i]=shifT[0][i-1]+1;
        }
        for (int i=shift+1;i<1+shift*2;i++){
            shifT[0][i]=shifT[0][i-1]-1;
        }
        //create each column
        for (int ii=0;ii<nbShift;ii++){
            for (int i=1;i<=shift;i++){
                shifT[i][ii]=shifT[i-1][ii]+shifT[0][ii];
            }
        }
        for (int ii=0;ii<nbShift;ii++){
            for (int i=shift+1;i<1+shift*2;i++){
                shifT[i][ii]=shifT[i-1][ii]-shifT[0][ii];
            }
        }
        
        
        double x;
        double y;
        double z;
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        
        
        ImageStack ims = new ImageStack(width,height);
        
        ImageProcessor [] ip = new FloatProcessor[depth];
        for (int t=0;t<depth;t++){
            ip[t]= new FloatProcessor(width,height);
//            for (int u=0;u<width;u++){
//                for (int uu=0;uu<height;uu++){
//                    ip[t].putPixelValue(u, uu, 0);
//                }
//            }
        }
        int xx,yy;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    xx=(int)((x-minX)/pixelsizeNM)+shift;
                    yy=(int)((y-minY)/pixelsizeNM)+shift;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if ((z>=minZ)&&(z<maxZ)){
                        //IJ.log("z "+z);
                        if ((xx-shift>=0)&&(yy-shift>=0)&&(xx+shift<width)&&(yy+shift<height)){
                            for (int a=-shift,aa=0;a<=shift;a++,aa++){
                                for (int b=-shift,bb=0;b<=shift;b++,bb++){
                                    if ((xx+a>=0)&&(xx+a<width)&&(yy+b>=0)&&(yy+b<height)){
                                        ip[0].putPixelValue(xx+a, yy+b, ip[0].getPixelValue(xx+a, yy+b)+shifT[aa][bb]);
                                    }
                                }
                            }
                        }
                    }
                    
                }
            }
        }
        for (int t=0;t<depth;t++){
            ims.addSlice(ip[t]);
        }
        ImagePlus imp = new ImagePlus(""+nameHistPlot+"",ims);
        imp.show();
        IJ.run("Set Scale...", "distance=1 known="+pixelsizeNM+" unit=nm");
        return imp;
        
    }
    
    
    
    
    
    
    
    
    
    
    
    public static ImagePlus hist2Dgaussian(StackLocalization sl,double pixelsizeNM) {
        
        ZRendering.makeLUT();
        
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        return hist2Dgaussian(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ);
        
    }
        
    
        
    
    public static ImagePlus hist2Dgaussian(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ)  {
        
        if (maxZ-minZ<pixelsizeNM){
            maxZ=minZ+pixelsizeNM;
        }
        if (maxY-minY<pixelsizeNM){
            maxY=minY+pixelsizeNM;
        }
        if (maxX-minX<pixelsizeNM){
            maxX=minX+pixelsizeNM;
        }
        double binColor=255;
        int width=(int)Math.ceil((maxX-minX)/pixelsizeNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeNM);
        int depth=1;
        
        
        
        double x;
        double y;
        double z;
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        
        
        ImageStack ims = new ImageStack(width,height);
        
        ImageProcessor [] ip = new FloatProcessor[depth];
        for (int t=0;t<depth;t++){
            ip[t]= new FloatProcessor(width,height);
//            for (int u=0;u<width;u++){
//                for (int uu=0;uu<height;uu++){
//                    ip[t].putPixelValue(u, uu, 0);
//                }
//            }
        }
        
        
        
        double cste=1/(2*Math.PI);
        int xx,yy,shiftA,shiftAA;
        double decx,decy;
        double sx,sy;//sigma value
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=((sl.fl.get(i).loc.get(j).X-minX)/pixelsizeNM);
                    y=((sl.fl.get(i).loc.get(j).Y-minY)/pixelsizeNM);
                    xx=(int)x;
                    yy=(int)y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    
                    sx=sl.fl.get(i).loc.get(j).crlb_X;
                    sy=sl.fl.get(i).loc.get(j).crlb_Y;
                    if (Math.random()<.001){
                        IJ.log("z"+z+"  "+sx+"  "+sy);
                    }
                    if (sx<=0){
                        sx=pixelsizeNM;
                    }
                    if (sy<=0){
                        sy=pixelsizeNM;
                    }
                    double leftGauss=cste/(sx*sy);
                    decx=(x-xx)*pixelsizeNM;
                    decy=(y-yy)*pixelsizeNM;
                    
                    shiftA=(int)Math.ceil(2*sx/pixelsizeNM);
                    shiftAA=(int)Math.ceil(2*sy/pixelsizeNM);
                    if ((z>=minZ)&&(z<maxZ)){
                        //IJ.log("z "+z);
                        for (int a=-shiftA;a<=shiftA;a++){
                            for (int aa=-shiftAA;aa<=shiftAA;aa++){
                                if ((xx+a>=0)&&(xx+a<width)&&(yy+aa>=0)&&(yy+aa<height)){
                                    ip[0].putPixelValue(xx+a, yy+aa, ip[0].getPixelValue(xx+a, yy+aa)+leftGauss*Math.exp(-.5*(((a*pixelsizeNM+decx)*(a*pixelsizeNM+decx)/(sx*sx)) + ((aa*pixelsizeNM+decy)*(aa*pixelsizeNM+decy)/(sy*sy)) )));
                                }
                            }

                        }
                    }
                    
                }
            }
        }
        
        
        for (int t=0;t<depth;t++){
            ims.addSlice(ip[t]);
        }
        ImagePlus imp = new ImagePlus(""+ZRendering.nameHistPlot+"",ims);
        imp.show();
        IJ.run("Set Scale...", "distance=1 known="+pixelsizeNM+" unit=nm");
        return imp;
        
    }
    
    
    
    
    
    
    
    
    public static ImagePlus hist3D(StackLocalization sl,double pixelsizeNM,int shift) {
        
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        return hist3D(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,shift);
        
        
    }
        
        
        
        
    public static ImagePlus hist3D(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ,int shift)  {
        return hist3D(sl,pixelsizeNM,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ,shift);
    }
    
    
    
    
        
    public static ImagePlus hist3D(StackLocalization sl,double pixelsizeXYNM,double pixelsizeZNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ,int shift)  {
        if (maxZ-minZ<pixelsizeZNM){
            maxZ=minZ+pixelsizeZNM;
        }
        if (maxY-minY<pixelsizeXYNM){
            maxY=minY+pixelsizeXYNM;
        }
        if (maxX-minX<pixelsizeXYNM){
            maxX=minX+pixelsizeXYNM;
        }
        double binColor=255;
        int width=(int)Math.ceil((maxX-minX)/pixelsizeXYNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeXYNM);
        int depth=(int)Math.ceil((maxZ-minZ)/pixelsizeZNM);
        
        width+=shift*2;
        height+=shift*2;
        depth+=shift*2;
        
        int nbShift=1+shift*2;
        int [][][] shifT = new int[1+shift*2][1+shift*2][1+shift*2];
        
        //create first element
        shifT[0][0][0]=1;
        //create first line
        for (int i=1;i<=shift;i++){
            shifT[0][0][i]=shifT[0][0][i-1]+1;
        }
        for (int i=shift+1;i<1+shift*2;i++){
            shifT[0][0][i]=shifT[0][0][i-1]-1;
        }
        //create each column
        for (int ii=0;ii<nbShift;ii++){
            for (int i=1;i<=shift;i++){
                shifT[0][i][ii]=shifT[0][i-1][ii]+shifT[0][0][ii];
            }
        }
        for (int ii=0;ii<nbShift;ii++){
            for (int i=shift+1;i<1+shift*2;i++){
                shifT[0][i][ii]=shifT[0][i-1][ii]-shifT[0][0][ii];
            }
        }
        //create each Z
        for (int iii=0;iii<nbShift;iii++){
            for (int ii=0;ii<nbShift;ii++){
                for (int i=1;i<=shift;i++){
                    shifT[i][ii][iii]=shifT[i-1][ii][iii]+shifT[0][ii][iii];
                }
            }
        }
        for (int iii=0;iii<nbShift;iii++){
            for (int ii=0;ii<nbShift;ii++){
                for (int i=shift+1;i<1+shift*2;i++){
                    shifT[i][ii][iii]=shifT[i-1][ii][iii]-shifT[0][ii][iii];
                }
            }
        }
        
        
        double x;
        double y;
        double z;
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        depth=Math.max(depth, 1);
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        
        
        
        
        ImageStack ims = new ImageStack(width,height);
        
        ImageProcessor [] ip = new FloatProcessor[depth];
        for (int t=0;t<depth;t++){
            ip[t]= new FloatProcessor(width,height);
//            for (int u=0;u<width;u++){
//                for (int uu=0;uu<height;uu++){
//                    ip[t].putPixelValue(u, uu, 0);
//                }
//            }
        }
        int xx,yy,zz;
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    xx=(int)((x-minX)/pixelsizeXYNM)+shift;
                    yy=(int)((y-minY)/pixelsizeXYNM)+shift;
                    zz=(int)((z-minZ)/pixelsizeZNM)+shift;
                    //IJ.log("z "+z);
                    if ((xx-shift>=0)&&(yy-shift>=0)&&(xx+shift<width)&&(yy+shift<height)){
                        for (int a=-shift,aa=0;a<=shift;a++,aa++){
                            for (int b=-shift,bb=0;b<=shift;b++,bb++){
                                for (int c=-shift,cc=0;c<=shift;c++,cc++){
                                    if ((zz+c>=0)&&(zz+c<ip.length)&&(xx+a>=0)&&(xx+a<width)&&(yy+b>=0)&&(yy+b<height)){
                                        ip[zz+c].putPixelValue(xx+a, yy+b, ip[zz+c].getPixelValue(xx+a, yy+b)+shifT[aa][bb][cc]);
                                    }
                                }
                            }
                        }
                    }
                    
                    
                }
            }
        }
        for (int t=0;t<depth;t++){
            ims.addSlice("z="+(((double)t)*pixelsizeZNM), ip[t]);
        }
        ImagePlus imp = new ImagePlus(""+ZRendering.nameHistPlot+"",ims);
        imp.show();
        IJ.run("Set Scale...", "distance=1 known="+pixelsizeXYNM+" unit=nm");
        if (pixelsizeXYNM!=pixelsizeZNM){
            IJ.log("WARNING:  X/Y pixel size="+pixelsizeXYNM+"  //  Z pixel size="+pixelsizeZNM);
        }
        return imp;
        
    }
    
    
    
    
    
    
    
    
    
    public static ImagePlus hist3Dgaussian(StackLocalization sl,double pixelsizeNM) {
        
        double minX=Double.POSITIVE_INFINITY;
        double maxX=Double.NEGATIVE_INFINITY;
        
        double minY=Double.POSITIVE_INFINITY;
        double maxY=Double.NEGATIVE_INFINITY;
        
        double minZ=Double.POSITIVE_INFINITY;
        double maxZ=Double.NEGATIVE_INFINITY;
        
        double x;
        double y;
        double z;
        
        //maybe use Arrays.sort to be fast
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=sl.fl.get(i).loc.get(j).X;
                    y=sl.fl.get(i).loc.get(j).Y;
                    z=sl.fl.get(i).loc.get(j).Z;
                    if (x<minX){
                        minX=x;
                    }
                    if (y<minY){
                        minY=y;
                    }
                    if (x>maxX){
                        maxX=x;
                    }
                    if (y>maxY){
                        maxY=y;
                    }
                    if (z<minZ){
                        minZ=z;
                    }
                    if (z>maxZ){
                        maxZ=z;
                    }
                }
            }
        }
        return hist3Dgaussian(sl,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ);
        
        
    }
        
        
        
        
    public static ImagePlus hist3Dgaussian(StackLocalization sl,double pixelsizeNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ)  {
        return hist3Dgaussian(sl,pixelsizeNM,pixelsizeNM,minX,maxX,minY,maxY,minZ,maxZ);
    }
    
    
    
    
        
    public static ImagePlus hist3Dgaussian(StackLocalization sl,double pixelsizeXYNM,double pixelsizeZNM,double minX, double maxX, double minY, double maxY,double minZ, double maxZ)  {
        if (maxZ-minZ<pixelsizeZNM){
            maxZ=minZ+pixelsizeZNM;
        }
        if (maxY-minY<pixelsizeXYNM){
            maxY=minY+pixelsizeXYNM;
        }
        if (maxX-minX<pixelsizeXYNM){
            maxX=minX+pixelsizeXYNM;
        }
        double binColor=255;
        int width=(int)Math.ceil((maxX-minX)/pixelsizeXYNM);
        int height=(int)Math.ceil((maxY-minY)/pixelsizeXYNM);
        int depth=(int)Math.ceil((maxZ-minZ)/pixelsizeZNM);
        
        
        
        
        double x;
        double y;
        double z;
        
        width=Math.max(width, 1);
        height=Math.max(height, 1);
        depth=Math.max(depth, 1);
        
        //IJ.log("depth "+depth+"  "+minZ+"  "+maxZ);
        
        
        
        
        ImageStack ims = new ImageStack(width,height);
        
        ImageProcessor [] ip = new FloatProcessor[depth];
        for (int t=0;t<depth;t++){
            ip[t]= new FloatProcessor(width,height);
//            for (int u=0;u<width;u++){
//                for (int uu=0;uu<height;uu++){
//                    ip[t].putPixelValue(u, uu, 0);
//                }
//            }
        }
        double cste=1/Math.pow(2*Math.PI,1.5);
        int xx,yy,zz,shiftA,shiftAA,shiftAAA;
        double decx,decy,decz;
        double sx,sy,sz;//sigma value
        
        for (int i=0;i<sl.fl.size();i++){
            for (int j=0;j<sl.fl.get(i).loc.size();j++){
                if (sl.fl.get(i).loc.get(j).exists){
                    x=((sl.fl.get(i).loc.get(j).X-minX)/pixelsizeXYNM);
                    y=((sl.fl.get(i).loc.get(j).Y-minY)/pixelsizeXYNM);
                    z=((sl.fl.get(i).loc.get(j).Z-minZ)/pixelsizeZNM);
                    
                    xx=(int)x;
                    yy=(int)y;
                    zz=(int)z;
                    sx=sl.fl.get(i).loc.get(j).crlb_X;
                    sy=sl.fl.get(i).loc.get(j).crlb_Y;
                    sz=sl.fl.get(i).loc.get(j).crlb_Z;
                    
                    if (sx<=0){
                        sx=pixelsizeXYNM;
                    }
                    if (sy<=0){
                        sy=pixelsizeXYNM;
                    }
                    if (sz<=0){
                        sz=pixelsizeZNM;
                    }
                    double leftGauss=cste/(sx*sy*sz);
                    decx=(x-xx)*pixelsizeXYNM;
                    decy=(y-yy)*pixelsizeXYNM;
                    decz=(z-zz)*pixelsizeZNM;
                    //if (Math.random()<.01){
                        //IJ.log("z "+z+"  "+x+"  "+y+"  "+sx+"  "+sy+"  "+sz);
                    //}
                    shiftA=(int)Math.ceil(2*sx/pixelsizeXYNM);
                    shiftAA=(int)Math.ceil(2*sy/pixelsizeXYNM);
                    shiftAAA=(int)Math.ceil(2*sz/pixelsizeZNM);
                    for (int a=-shiftA;a<=shiftA;a++){
                        for (int aa=-shiftAA;aa<=shiftAA;aa++){
                            for (int aaa=-shiftAAA;aaa<=shiftAAA;aaa++){
                                if ((zz+aaa>=0)&&(zz+aaa<ip.length)&&(xx+a>=0)&&(xx+a<width)&&(yy+aa>=0)&&(yy+aa<height)){
                                    ip[zz+aaa].putPixelValue(xx+a, yy+aa, ip[zz+aaa].getPixelValue(xx+a, yy+aa)+leftGauss*Math.exp(-.5*(((a*pixelsizeXYNM+decx)*(a*pixelsizeXYNM+decx)/(sx*sx)) + ((aa*pixelsizeXYNM+decy)*(aa*pixelsizeXYNM+decy)/(sy*sy)) + ((aaa*pixelsizeXYNM+decz)*(aaa*pixelsizeXYNM+decz)/(sz*sz)))));
                                }
                            }
                        }
                        
                    }
                    
                    
                    
                }
            }
        }
        for (int t=0;t<depth;t++){
            ims.addSlice("z="+(((double)t)*pixelsizeZNM), ip[t]);
        }
        ImagePlus imp = new ImagePlus(""+ZRendering.nameHistPlot+"",ims);
        imp.show();
        IJ.run("Set Scale...", "distance=1 known="+pixelsizeXYNM+" unit=nm");
        if (pixelsizeXYNM!=pixelsizeZNM){
            IJ.log("WARNING:  X/Y pixel size="+pixelsizeXYNM+"  //  Z pixel size="+pixelsizeZNM);
        }
        return imp;
        
    }
    
    
    
    
    
    
    public static void colorizeHist(ImagePlus imp){
        
        int width=imp.getWidth();
        int height=imp.getHeight();
        int depth=imp.getNSlices();
        double minrange=imp.getDisplayRangeMin();
        double maxrange=imp.getDisplayRangeMax();
        ImageStack ims = new ImageStack(width,height);
        double maxval=100;
        ColorProcessor [] ip = new ColorProcessor[depth];
        for (int t=0;t<depth;t++){
            imp.setSlice(t+1);
            ImageProcessor ipin=imp.getProcessor();
            ip[t]= new ColorProcessor(width,height);
            for (int u=0;u<width;u++){
                for (int uu=0;uu<height;uu++){
                    double val=ipin.getPixelValue(u, uu);
                    if (val!=0){
                        double angle=1.-((double)t)/depth;
                        double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                        if (H<0){
                            H+=360;
                        }
                        double V=0,S=0;
                        
                        /*if (val-minrange<((maxrange-minrange)/2)){
                            S=1;
                            V =((val)*2-minrange*2)/(maxrange-minrange);
                        }
                        else{
                            S=1-((val-((maxrange-minrange)/2))*2-minrange)/(maxrange-minrange);
                            V =1;
                        }*/
                        
                        S=1;
                        V =((val)-minrange)/(maxrange-minrange);
                            
                        V=Math.max(Math.min(V,1), 0);
                        S=Math.max(Math.min(S,1), 0);
                        int [] RGB=getRGB(H,S,V);
                        ip[t].putPixel(u, uu,RGB);
                    }
                    
                }
            }
            ims.addSlice(ip[t]);
        }
        ImagePlus impout = new ImagePlus("colorized 3D histogram",ims);
        impout.show();
        
    }
    
    
    
    
    
    
    
    
    public static void localMaxima(ImagePlus imp){
        
        int width=imp.getWidth();
        int height=imp.getHeight();
        int depth=imp.getNSlices();
        ImageStack ims = new ImageStack(width,height);
        FloatProcessor [] ip = new FloatProcessor[depth];
        double maxi;
        for (int t=0;t<depth;t++){
            imp.setSlice(t+1);
            ImageProcessor ipin=imp.getProcessor();
            ip[t]= new FloatProcessor(width,height);
            for (int u=1;u<width-1;u++){
                for (int uu=1;uu<height-1;uu++){
                    
                    maxi=Double.NEGATIVE_INFINITY;
                    
                    for (int a=-1;a<=1;a++){
                        for (int aa=-1;aa<=1;aa++){
                            if (ipin.getPixelValue(u+a, uu+aa)>maxi){
                                maxi=ipin.getPixelValue(u+a, uu+aa);
                            }
                        }
                    }
                    if (maxi==ipin.getPixelValue(u, uu)){
                        ip[t].putPixelValue(u, uu, maxi);
                    }
                    
                    
                }
            }
            ims.addSlice(ip[t]);
        }
        ImagePlus impout = new ImagePlus("local maxima image",ims);
        impout.show();
        
    }
    
    
    
    
    
    
    
    private static void printLUT(double minZ, double maxZ){
        int dec=5;
        int width=80;
        int height=150+dec*4;
        ImageProcessor ip = new ColorProcessor(width,height);
        for (int i=0;i<width;i++){
            for (int ii=0;ii<height;ii++){
                int [] RGB={0,0,0};
                ip.putPixel(i, ii,RGB);
            }
        }
        for (int i=0;i<width/4;i++){
            for (int ii=10,iii=0;ii<height-10;ii++,iii++){
                double angle=((double)iii)/((double)(height-dec*4));
                double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                if (H<0){
                    H+=360;
                }
                double V=1,S=1;
                
                int [] RGB=getRGB(H,S,V);
                ip.putPixel(i, ii,RGB);
            }
            
        }
        Color c = new Color(255,255,255);
        ip.setColor(c);
        ip.drawString(""+(int)Math.round(maxZ),width/4+2,3*dec);
        ip.drawString(""+(int)Math.round((maxZ+minZ)/2),width/4+2,height/2);
        ip.drawString(""+(int)Math.round(minZ),width/4+2,height-dec);
        ImagePlus imp = new ImagePlus("Calibration bar",ip);
        imp.show();
    }
    
    
    private static void addLUT(ImagePlus imp,double minZ, double maxZ){
        addLUT(imp,minZ, maxZ,null);
    }
    
    private static void addLUT(ImagePlus imp,double minZ, double maxZ,String pathLut){
        
        byte [] r=null; 
        byte [] g=null; 
        byte [] b=null; 
        
        if (pathLut!=ZRendering.lut){
            try{
                LUT lut=LutLoader.openLut(pathLutFolder+pathLut);
                r = new byte[256] ;
                g = new byte[256] ;
                b = new byte[256] ;
                lut.getReds(r);
                lut.getBlues(b);
                lut.getGreens(g);
            }catch(Exception e ){pathLut=ZRendering.lut;}
            
        }
        
        if (minZ<0){
            maxZ-=minZ;
            minZ-=minZ;
        }
        
        ImageProcessor ip = imp.getProcessor();
        int w=ip.getWidth();
        int h=ip.getHeight();
        
        
        int sizeCh=(String.valueOf((int)Math.round(maxZ))).length();
        sizeCh=(int)Math.max(sizeCh,(String.valueOf((int)Math.round(minZ))).length());
        
        
        
        
        int heightbar=h/6;
        int widthbar=heightbar/8;
        int sizeFont=heightbar/6;
        int margin=sizeFont;
        int width=widthbar+(3*sizeFont/4)*sizeCh;  //  *3/4 because width charactere < than height character
        int height=heightbar;
        
        
        for (int i=w-width-margin;i<w-width-margin+widthbar;i++){
            for (int ii=margin,iii=0;ii<height+margin;ii++,iii++){
                
                double angle=((double)iii)/((double)(height));
                int [] RGB=null;
                if (pathLut==ZRendering.lut){
                    double H =-45.+((315.)*angle); //convert in degree 315->max color=red     +45 degres: min color blue
                    if (H<0){
                        H+=360;
                    }
                    double V=1,S=1;

                    RGB=getRGB(H,S,V);
                }
                else{
                    RGB=new int[3];
                    int red=r[(int)Math.ceil((1-angle)*256)-1];
                    int green=g[(int)Math.ceil((1-angle)*256)-1];
                    int blue=b[(int)Math.ceil((1-angle)*256)-1];
                    if (red<0){
                        red+=256;
                    }
                    if (green<0){
                        green+=256;
                    }
                    if (blue<0){
                        blue+=256;
                    }
                    RGB[0]=red;
                    RGB[1]=green;
                    RGB[2]=blue;
                }
                ip.putPixel(i, ii,RGB);
            }
            
        }
        Color c = new Color(255,255,255);
        Color cb = new Color(0,0,0);
        Font f = new Font(Font.SANS_SERIF,Font.PLAIN,sizeFont);
        ip.setFont(f);
        ip.setColor(c);
        ip.drawString(""+(int)Math.round(maxZ),w-width-margin+widthbar+2,margin+sizeFont/2,cb);
        ip.drawString(""+(int)Math.round((maxZ+minZ)/2),w-width-margin+widthbar+2,height/2+margin+sizeFont/2,cb);
        ip.drawString(""+(int)Math.round(minZ),w-width-margin+widthbar+2,height+margin+sizeFont/2,cb);
        
        imp.updateAndDraw();
    }
    
    
    
    
    
    
    
    

}
