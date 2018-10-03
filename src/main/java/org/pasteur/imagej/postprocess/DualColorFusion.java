


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;
import org.pasteur.imagej.data.*;
import org.pasteur.imagej.utils.*;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ImageProcessor;
import java.util.Arrays;
import java.util.Comparator;
import java.util.ArrayList;
/**
 *
 * @author benoit
 */
public class DualColorFusion {
    
    StackLocalization sl1;
    StackLocalization sl2;
    
    StackLocalization sl1fuse;
    StackLocalization sl2fuse;
    
    
    int consecutiveFrameOn=100;
    
    int pixelSize2D=256;
    
    double pixsize;//nm
    int magnification=2;
    int magnificationConvolution=1;
    ImagePlus im1;
    ImagePlus im2;
    ImagePlus im3;
    
    
    StackLocalization sl2Registered;
    
    
    String pathRes;
    String pathLocRes;
    
    double maxDistanceMergingXY;
    
    int [][] idLoc;
    
    
    double sizePix=100;//nm
    double zstep=100;//nm
    
    double maxX=Double.NEGATIVE_INFINITY;
    double minX=Double.MAX_VALUE;
    double maxY=Double.NEGATIVE_INFINITY;
    double minY=Double.MAX_VALUE;
    double maxZ=Double.NEGATIVE_INFINITY;
    double minZ=Double.MAX_VALUE;
       
    double [][] matrix1;
    double [][] matrix2;
    int order;
    
    int way=0;
    double sigmaGaussian;//nm
    
    //constructor to compute registration
    public DualColorFusion(StackLocalization sl1,StackLocalization sl2,String pathRes,double maxDistanceMergingXY, int order){
        
        this.order=order;
        this.sl1=sl1;
        this.sl2=sl2;
        way=1;
        
        this.maxDistanceMergingXY=maxDistanceMergingXY;
        
        this.pathRes=pathRes;
    }
    
    
    
    
    //constructor to compute registration with widefield
    public DualColorFusion(ImagePlus im1,StackLocalization sl2,String pathRes,double pixSize,double zstep, int magnification,double sigma,ImagePlus im2,ImagePlus im3){
        this.sigmaGaussian=sigma;
        this.im1=im1;
        this.im2=im2;
        this.im3=im3;
        
        
        if ((im2!=null)&&(im1.getWidth()!=im2.getWidth())&&(im1.getHeight()!=im2.getHeight())&&(im1.getNSlices()!=im2.getNSlices())&&(im1.getNFrames()!=im2.getNFrames())){
            IJ.log("WARNING: dimensions of widefield image 1 and 2 do not match");
            this.im2=null;
        }
        if ((im3!=null)&&(im1.getWidth()!=im3.getWidth())&&(im1.getHeight()!=im3.getHeight())&&(im1.getNSlices()!=im3.getNSlices())&&(im1.getNFrames()!=im3.getNFrames())){
            IJ.log("WARNING: dimensions of widefield image 1 and 3 do not match");
            this.im3=null;
        }
        
        this.sl2=sl2;
        this.pixsize=pixSize;
        this.magnification=magnification;
        way=2;
        this.zstep= zstep;
        
        this.pathRes=pathRes;
    }
    
    
    //constructor to apply registration
    public DualColorFusion(StackLocalization sl2,String path_registration,String pathLocRes){
        way=0;
        this.pathRes=path_registration;
        this.pathLocRes=pathLocRes;
        this.sl2=sl2;
    }
    
    
    public void run(){
        if (way==1){
            run1();
        }
        else if (way==2){
            run2();
        }
        else{
            run0();
        }
            
    }
    
    
    
    
    //just apply registration
    public void run0(){
        
        
        
        PolynomialFit pf = new PolynomialFit(pathRes);
        double [] v = new double [3];
        for (int i=0;i<sl2.fl.size();i++){
            IJ.showProgress((double)i/(double)sl2.fl.size());
            for (int j=0;j<sl2.fl.get(i).loc.size();j++){
                v[0]=sl2.fl.get(i).loc.get(j).X;
                v[1]=sl2.fl.get(i).loc.get(j).Y;
                v[2]=sl2.fl.get(i).loc.get(j).Z;
                double [] res=pf.transform(v);
                sl2.fl.get(i).loc.get(j).X=res[0];
                sl2.fl.get(i).loc.get(j).Y=res[1];
                sl2.fl.get(i).loc.get(j).Z=res[2];
                
            }
        }
        
        sl2.save(pathLocRes);
        
        IJ.showProgress(1);
        
        
    }
    
    //registration to image
    public void run2(){
        
        
        
        boolean stacked=true;
        int nbImage=im1.getNSlices();
        if (nbImage==1){
            nbImage=im1.getNFrames();
            stacked=false;
        }
        if (nbImage==1){
            IJ.log("WARNING, you should provide a stack with more than 1 image");
        }
        ImageStack ims= im1.getStack();
        
        int width=im1.getWidth();
        int height=im1.getHeight();
        
        int sizeFFTZ=(int)Math.pow(2,(int)(Math.ceil(Math.log(nbImage*magnificationConvolution*2)/Math.log(2))));
        int sizeFFTX=(int)Math.pow(2,(int)(Math.ceil(Math.log(width*magnificationConvolution*2)/Math.log(2))));
        int sizeFFTY=(int)Math.pow(2,(int)(Math.ceil(Math.log(height*magnificationConvolution*2)/Math.log(2))));
        
        float [][][] image =new float [sizeFFTZ][sizeFFTX][sizeFFTY];//*2 for padding
        float [][][] imagewidefield1 =new float [nbImage*magnification][width*magnification][height*magnification];//*2 for padding
        float [][][] imagewidefield2 =null;
        float [][][] imagewidefield3 =null;
        ImageStack ims2=null;
        if (im2!=null){
            imagewidefield2 =new float [nbImage*magnification][width*magnification][height*magnification];//*2 for padding
            ims2= im2.getStack();
        }
        
        ImageStack ims3=null;
        if (im3!=null){
            imagewidefield3 =new float [nbImage*magnification][width*magnification][height*magnification];//*2 for padding
            ims3= im3.getStack();
        }
        
        for (int i=0;i<nbImage;i++){
            ImageProcessor ip;
            ImageProcessor ip2=null;
            ImageProcessor ip3=null;
            if (stacked){
                ip = ims.getProcessor(i+1);
                if (im2!=null){
                    ip2 = ims2.getProcessor(i+1);
                }
                if (im3!=null){
                    ip3 = ims3.getProcessor(i+1);
                }
            }
            else{
                im1.setSlice(i+1);
                ip = im1.getProcessor();
                if (im2!=null){
                    im2.setSlice(i+1);
                    ip2 = im2.getProcessor();
                }
                if (im3!=null){
                    im3.setSlice(i+1);
                    ip3 = im3.getProcessor();
                }
            }
            for (int ii=0,jj=0;ii<width;ii++,jj++){
                for (int iii=0,jjj=0;iii<height;iii++,jjj++){
                    for (int u=0;u<magnificationConvolution;u++){
                        for (int uu=0;uu<magnificationConvolution;uu++){
                            for (int uuu=0;uuu<magnificationConvolution;uuu++){
                                image[i*magnificationConvolution+nbImage*magnificationConvolution/2+u][ii*magnificationConvolution+width*magnificationConvolution/2+uu][iii*magnificationConvolution+height*magnificationConvolution/2+uuu]=ip.getPixelValue(jj, jjj);
                            }
                        }
                    }
                    for (int u=0;u<magnification;u++){
                        for (int uu=0;uu<magnification;uu++){
                            for (int uuu=0;uuu<magnification;uuu++){
                                imagewidefield1[i*magnification+u][ii*magnification+uu][iii*magnification+uuu]=ip.getPixelValue(jj, jjj);
                                if (im2!=null){
                                    imagewidefield2[i*magnification+u][ii*magnification+uu][iii*magnification+uuu]=ip2.getPixelValue(jj, jjj);
                                }
                                if (im3!=null){
                                    imagewidefield3[i*magnification+u][ii*magnification+uu][iii*magnification+uuu]=ip3.getPixelValue(jj, jjj);
                                }
                            }
                        }
                    }
                    
                }
            }
        }
        
        IJ.showProgress(.1);
        
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
        for (int i=0;i<sl2.fl.size();i++){
            for (int j=0;j<sl2.fl.get(i).loc.size();j++){
                if (sl2.fl.get(i).loc.get(j).exists){
                    x=sl2.fl.get(i).loc.get(j).X;
                    y=sl2.fl.get(i).loc.get(j).Y;
                    z=sl2.fl.get(i).loc.get(j).Z;
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
        
        
        float [][][] image_ =new float [sizeFFTZ][sizeFFTX][sizeFFTY];
        
        ImagePlus imp=org.pasteur.imagej.postprocess.ZRendering.hist3D(sl2, pixsize/(double)magnificationConvolution, zstep/(double)magnificationConvolution, minX, maxX, minY, maxY, minZ, maxZ, 0);
        
        for (int i=0;i<sizeFFTZ/2;i++){
            int idimage=i+1;
            if (idimage<=imp.getNSlices()){
                imp.setSlice(idimage);
                ImageProcessor ip = imp.getProcessor();
                for (int ii=0,jj=0;ii<sizeFFTX/2;ii++,jj++){
                    for (int iii=0,jjj=0;iii<sizeFFTY/2;iii++,jjj++){
                        if (jj<imp.getWidth()&&jjj<imp.getHeight()){
                            image_[i+nbImage*magnificationConvolution/2][ii+width*magnificationConvolution/2][iii+height*magnificationConvolution/2]=ip.getPixelValue(jj, jjj);
                        }
                    }
                }
            }
                
        }
        
        
        
        float [][][] tmp1 =new float [sizeFFTZ][sizeFFTX][sizeFFTY];
        float [][][] tmp2 =new float [sizeFFTZ][sizeFFTX][sizeFFTY];
        float [][][] tmp3 =new float [sizeFFTZ][sizeFFTX][sizeFFTY];
        float [][][] tmp4 =new float [sizeFFTZ][sizeFFTX][sizeFFTY];
        float [][][] tmp5 =new float [sizeFFTZ][sizeFFTX][sizeFFTY];
        
        IJ.showProgress(.2);
        IJ.log("images created");
        
        double sigma=sigmaGaussian/pixsize;
        double sigmaFourierPowx=((1./(sigma))*(double)sizeFFTX/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTX/(Math.PI*2.));
        double sigmaFourierPowy=((1./(sigma))*(double)sizeFFTY/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTY/(Math.PI*2.));
        double sigmaFourierPowz=((1./(sigma))*(double)sizeFFTZ/(Math.PI*2.))*((1./(sigma))*(double)sizeFFTZ/(Math.PI*2.));
        
        for (int zzz=0;zzz<sizeFFTZ;zzz++){
            for (int xxx=0;xxx<sizeFFTX;xxx++){
                for (int yyy=0;yyy<sizeFFTY;yyy++){
                    
                    tmp2[zzz][xxx][yyy]=(float)Math.exp(-.5* ((xxx-sizeFFTX/2)*(xxx-sizeFFTX/2)/sigmaFourierPowx + (yyy-sizeFFTY/2)*(yyy-sizeFFTY/2)/sigmaFourierPowy +(zzz-sizeFFTZ/2)*(zzz-sizeFFTZ/2)/sigmaFourierPowz));

                    
                }
            }
        }
        FourierTransform ft = new FourierTransform();
        
        
        
        ft.shift3D(tmp2, tmp3);
        
        IJ.showProgress(.3);
        
        //org.pasteur.imagej.utils.ImageShow.imshow(image,"image");//////////////////////////////////:*****************************
        //org.pasteur.imagej.utils.ImageShow.imshow(image_,"image_");//////////////////////////////////:*****************************
        
        ft.fft(image_, tmp1,tmp4, tmp5);
        IJ.showProgress(.4);
        IJ.log("image SR FFT ok");
        
        for (int ii=0;ii<sizeFFTZ;ii++){
            for (int iii=0;iii<sizeFFTX;iii++){
                for (int iiii=0;iiii<sizeFFTY;iiii++){
                    tmp4[ii][iii][iiii]=tmp4[ii][iii][iiii]*tmp3[ii][iii][iiii];
                    tmp5[ii][iii][iiii]=tmp5[ii][iii][iiii]*tmp3[ii][iii][iiii];
                }
            }
        }
        IJ.showProgress(.5);
        
        IJ.log("gaussian blur ok");
        
        
        
        ft.fft(image, tmp1,tmp2, tmp3);
        IJ.log("image LR FFT ok");
        
        
        IJ.showProgress(.6);
        
        for (int ii=0;ii<sizeFFTZ;ii++){
            for (int iii=0;iii<sizeFFTX;iii++){
                for (int iiii=0;iiii<sizeFFTY;iiii++){
                    image[ii][iii][iiii]=tmp2[ii][iii][iiii]*tmp4[ii][iii][iiii]+tmp3[ii][iii][iiii]*tmp5[ii][iii][iiii];
                    image_[ii][iii][iiii]=tmp3[ii][iii][iiii]*tmp4[ii][iii][iiii]-tmp2[ii][iii][iiii]*tmp5[ii][iii][iiii];
                    
                    
                    
                }
            }
        }
        IJ.showProgress(.7);
        
        ft.ifft(image, image_,tmp2, tmp1);
        IJ.showProgress(.8);
        ft.shift3D(tmp2,tmp1);
        IJ.log("image 1 IFFT");
        
        //org.pasteur.imagej.utils.ImageShow.imshow(tmp1,"conv");//////////////////////////////////:*****************************
        IJ.showProgress(.9);
        LocalMaxima ml = new LocalMaxima(tmp1);
        ml.run(.95);
        
        double xs=ml.getShiftXinPix()/magnificationConvolution*pixsize+minX;
        double ys=ml.getShiftYinPix()/magnificationConvolution*pixsize+minY;
        double zs=ml.getShiftZinPix()/magnificationConvolution*zstep+minZ;
        
        IJ.log("shift: "+xs+"  "+ys+"  "+zs+"  XS:"+ml.getShiftXinPix()/magnificationConvolution+"  minX:"+minX);
        
        for (int i=0;i<sl2.fl.size();i++){
            for (int ii=0;ii<sl2.fl.get(i).loc.size();ii++){
                PLocalization p = sl2.fl.get(i).loc.get(ii);
                p.X-=xs;
                p.Y-=ys;
                p.Z-=zs;
            }
        }
        
        imp.close();
        
        
        sl2.saveCSV(pathRes);
        IJ.log("localization table registered saved:"+pathRes);
        IJ.log("To render the image with the same size as widefield image, please set the parameters as follow:");
        
        IJ.log("Pixel size: "+ pixsize);
        IJ.log("Axial pixel size: "+ zstep);
        IJ.log("shift histogram: "+ 0);
        IJ.log("3d rendering: ");
        IJ.log("min X: "+0);
        IJ.log("max X: "+width*pixsize);
        IJ.log("min Y: "+0);
        IJ.log("max Y: "+height*pixsize);
        IJ.log("min Z: "+0);
        IJ.log("max Z: "+nbImage*zstep);
        
        ImageShow.imshow(imagewidefield1,"wf");
        ImagePlus impfinal=org.pasteur.imagej.postprocess.ZRendering.hist3D(sl2, pixsize/magnification, zstep/magnification, 0, width*pixsize, 0, height*pixsize, 0, nbImage*zstep, 0);
        impfinal.setTitle("hr");
        try{Thread.sleep(500);}catch(Exception e){}
        IJ.run("Merge Channels...", "c1=wf c2=hr create");
        
        if (im2!=null){
            ImageShow.imshow(imagewidefield2,"wf");
            ImagePlus impfinal2=org.pasteur.imagej.postprocess.ZRendering.hist3D(sl2, pixsize/magnification, zstep/magnification, 0, width*pixsize, 0, height*pixsize, 0, nbImage*zstep, 0);
            impfinal2.setTitle("hr");
            if (im3!=null){
                ImageShow.imshow(imagewidefield3,"wfb");
                try{Thread.sleep(500);}catch(Exception e){}
                IJ.run("Merge Channels...", "c2=wf c1=hr c3=wfb create");
            }
            else{
                try{Thread.sleep(500);}catch(Exception e){}
                IJ.run("Merge Channels...", "c2=wf c1=hr create");
            }
            
        }
        IJ.showProgress(1);
        
    }
    
    
    //registration to table
    public void run1(){
        
        IJ.showProgress(0);
        
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
            IJ.log("wow, there is a surprisingly high number of negative coordinates: "+(nbNegative/2));
        }
        
        PolynomialFit pf=register();
        
        pf.save(pathRes);
        
        IJ.showProgress(0);
    }
    
    
    
    
    
    public PolynomialFit register(){
        
        
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
                if (Math.max(data1[0][t],data2[0][t])>maxX){
                    maxX=Math.max(data1[0][t],data2[0][t]);
                }
                if (Math.max(data1[1][t],data2[1][t])>maxY){
                    maxY=Math.max(data1[1][t],data2[1][t]);
                }
                if (Math.min(data1[2][t],data2[2][t])<minZ){
                    minZ=Math.min(data1[2][t],data2[2][t]);
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
        PolynomialFit pf=null;
        if (notApparied){
            IJ.log("oops : the two files are not apparied -> registration not possible");
        }
        else{
            pf = new PolynomialFit(order,data1,data2);
            boolean bool=pf.run();
            if (!bool){
                pf=null;
            }
            else{
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
            }
            
            
        }
        
        
        
        ZRendering.hist2D(sl1fuse, 20, minX, maxX, minY, maxY, minZ, maxZ, 1,"camera1");
        ZRendering.hist2D(sl2fuse, 20, minX, maxX, minY, maxY, minZ, maxZ, 1,"camera2");
        
        
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
        plotHist(axisAfterXY,histAfterXY,"Histogram of X/Y distances (cam1 vs. cam2) after registration", "distance (nm)", "occurrence #");
        plotHist(axisAfterZ,histAfterZ,"Histogram of Z distances (cam1 vs. cam2) after registration", "distance (nm)", "occurrence #");
        
        return pf;
    }
    
    
    
    
    public void plotHist(double [] x, double [] y,String title,String xlabel,String ylabel){
        Plot p = new Plot(""+title,xlabel,ylabel,x,y,3);
        p.setFont(0, 18);
        p.setLineWidth(2);
        p.show();
        
    }
    
    
    
    
    
    
    
    class LocalMaxima{
        float [][][] image;
        double x,y,z;//position to estimate -> drift
        int sizePatchX=11;//odd value;
        int sizePatchY=11;//odd value;
        int sizePatchZ=11;//odd value;
        int width;
        int height;
        int depth;
        
        LocalMaxima(float [][][] image){
            this.image=image;
            width=image.length;
            height=image[0].length;
            depth=image[0][0].length;
            
            sizePatchX=Math.min(sizePatchX, width-1);
            sizePatchY=Math.min(sizePatchY, height-1);
            sizePatchZ=Math.min(sizePatchZ, depth-1);
            if ((width%2==1)||(height%2==1)||(depth%2==1)){
                IJ.log("WARNING, drift image computation should have even size");
            }
            
            if (sizePatchX%2==0){
                sizePatchX=sizePatchX-1;
            }
            if (sizePatchY%2==0){
                sizePatchY=sizePatchY-1;
            }
            if (sizePatchZ%2==0){
                sizePatchZ=sizePatchZ-1;
            }
            
        }
        
        
        public void run(double threshold){//threshold=.75 -> 75% of data used to compute center
            
            
            double max=Double.MIN_VALUE;
            
            int xx=0,yy=0,zz=0;
            
            for (int i=0,t=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    for (int iii=0;iii<depth;iii++){
                        if (max<image[i][ii][iii]){
                            max=image[i][ii][iii];
                            xx=i;
                            yy=ii;
                            zz=iii;
                        }
                    }
                }
            }
            
            
            float [] vect=new float[sizePatchX*sizePatchY*sizePatchZ];
            int number=0;
            for (int i=0;i<sizePatchX;i++){
                for (int ii=0;ii<sizePatchY;ii++){
                    for (int iii=0;iii<sizePatchZ;iii++){
                        if ((i+xx-sizePatchX/2>=0)&&(ii+yy-sizePatchY/2>=0)&&(iii+zz-sizePatchZ/2>=0)&&(i+xx-sizePatchX/2<width)&&(ii+yy-sizePatchY/2<height)&&(iii+zz-sizePatchZ/2<depth)){
                            vect[number++]=image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                        }
                    }
                }
            }
            //IJ.log("number "+number+"  "+xx+"  "+yy+"  "+zz);
            Arrays.sort(vect);
            
            //IJ.log("threshold "+vect[(int)((number-1)*threshold)]+"   "+max);
            
            x=0;
            y=0;
            z=0;
            double sum=0;
            for (int i=0;i<sizePatchX;i++){
                for (int ii=0;ii<sizePatchY;ii++){
                    for (int iii=0;iii<sizePatchZ;iii++){
                        if ((i+xx-sizePatchX/2>=0)&&(ii+yy-sizePatchY/2>=0)&&(iii+zz-sizePatchZ/2>=0)&&(i+xx-sizePatchX/2<width)&&(ii+yy-sizePatchY/2<height)&&(iii+zz-sizePatchZ/2<depth)){
                            if (image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2]>vect[(int)((number-1)*threshold)]){
                                sum+=image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                                x+=(i+xx-sizePatchX/2)*image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                                y+=(ii+yy-sizePatchY/2)*image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                                z+=(iii+zz-sizePatchZ/2)*image[i+xx-sizePatchX/2][ii+yy-sizePatchY/2][iii+zz-sizePatchZ/2];
                            }
                        }
                    }
                }
            }
            //IJ.log("posit "+x+"  "+y+"  "+z);
            x/=sum;
            y/=sum;
            z/=sum;
            //IJ.log("posit "+x+"  "+y+"  "+z);
        }
        
        
        //WARNING reversed because matrix order...
        public double getShiftZinPix(){
            return -(x-width/2);
        }
        
        public double getShiftXinPix(){
            return -(y-depth/2);
        }
        
        public double getShiftYinPix(){
            return -(z-height/2);
        }
        
        
        
    }
    
    
}
