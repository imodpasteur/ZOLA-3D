/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;

/**
 *
 * @author benoit
 */

public class ImageShow {
    public static void imshow(double [][][] im){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus("imshow",ims);
        imp.show();
    }
    
    
    
    
    
    public static void imshow(int [][][] im){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus("imshow",ims);
        imp.show();
    }
    
    
    
    
    public static void imshow(short [][][] im,String s){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new ShortProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, im[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus(""+s,ims);
        imp.show();
    }
    
    
    
    public static void imshow(int [][][] im,String s){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus(""+s,ims);
        imp.show();
    }
    
    
    public static void  imshow(ArrayList<double [][]> image,String s){
        if (image.size()!=0){
            ImageStack ims = new ImageStack(image.get(0).length,image.get(0)[0].length);
            for (int k=0;k<image.size();k++){
                ImageProcessor ip = new FloatProcessor(image.get(0).length,image.get(0)[0].length);
                for (int i=0;i<image.get(0).length;i++){
                    for (int ii=0;ii<image.get(0)[0].length;ii++){
                        ip.putPixelValue(i, ii, image.get(k)[i][ii]);
                    }
                }
                ims.addSlice("array",ip);
            }
            ImagePlus imp=new ImagePlus(""+s,ims);
            imp.show();
        }
    }
    
    public static void imshow(double [][][] im,String s){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus(""+s,ims);
        imp.show();
    }
    
    
    
    
    public static void imshow(float [][][] im,String s){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[k][i][ii]);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus(""+s,ims);
        imp.show();
    }
    
    
    public static void imshow(boolean [][][] im,String s){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    if (im[k][i][ii]){
                        ip.putPixelValue(i, ii, 1);
                    }
                    else{
                        ip.putPixelValue(i, ii, 0);
                    }
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus(""+s,ims);
        imp.show();
    }
    
    public static void imshow(double [][][] im,double [][][] imprime,String s){
        if ((im[0][0].length==imprime[0][0].length)&&(im.length==imprime.length)){
            ImageStack ims = new ImageStack(im[0].length+imprime[0].length,im[0][0].length);
            for (int k=0;k<im.length;k++){
                ImageProcessor ip = new FloatProcessor(im[0].length+imprime[0].length,im[0][0].length);
                for (int i=0;i<im[0].length;i++){
                    for (int ii=0;ii<im[0][0].length;ii++){
                        ip.putPixelValue(i, ii, im[k][i][ii]);
                        
                    }
                }
                for (int i=0;i<imprime[0].length;i++){
                    for (int ii=0;ii<imprime[0][0].length;ii++){
                        ip.putPixelValue(i+imprime[0].length, ii, imprime[k][i][ii]);
                        
                    }
                }
                ims.addSlice(ip);
            }
            ImagePlus imp=new ImagePlus(""+s,ims);
            imp.show();
        }
        else{
            IJ.log("imshow of 2 stacks at the same time impossible due to different dimensions");
        }
    }
    
    public static void imshow(double [][] im){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[i][ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus("imshow",ip);
        imp.show();
    }
    
    
    
    public static void imshow(boolean [][][] im){
        ImageStack ims = new ImageStack(im[0].length,im[0][0].length);
        for (int k=0;k<im.length;k++){
            ImageProcessor ip = new FloatProcessor(im[0].length,im[0][0].length);
            for (int i=0;i<im[0].length;i++){
                for (int ii=0;ii<im[0][0].length;ii++){
                    if (im[k][i][ii])
                        ip.putPixelValue(i, ii, 1);
                    else
                        ip.putPixelValue(i, ii, 0);
                }
            }
            ims.addSlice(ip);
        }
        ImagePlus imp=new ImagePlus("imshow",ims);
        imp.show();
    }
    
    
    public static void imshow(boolean [][] im){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    if (im[i][ii])
                        ip.putPixelValue(i, ii, 1);
                    else
                        ip.putPixelValue(i, ii, 0);
                }
            }
        
        ImagePlus imp=new ImagePlus("imshow",ip);
        imp.show();
    }
    
    public static void imshow(boolean [][] im,String title){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    if (im[i][ii])
                        ip.putPixelValue(i, ii, 1);
                    else
                        ip.putPixelValue(i, ii, 0);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    public static void imshow(int [][] im){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    ip.putPixelValue(i, ii, im[i][ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus("imshow",ip);
        imp.show();
    }
    
    public static void imshow(double [][] im,String title){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    ip.putPixelValue(i, ii, (float)im[i][ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    
    public static void imshow(double [] im,String title){
        
            ImageProcessor ip = new FloatProcessor((int)Math.sqrt(im.length),(int)Math.sqrt(im.length));
            for (int i=0;i<(int)Math.sqrt(im.length);i++){
                for (int ii=0;ii<(int)Math.sqrt(im.length);ii++){
                    ip.putPixelValue(i, ii, (float)im[i*(int)Math.sqrt(im.length)+ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    
    
    
    public static void imshow(int [][] im,String title){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    ip.putPixelValue(i, ii, (int)im[i][ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    public static void imshow(float [][] im,String title){
        
            ImageProcessor ip = new FloatProcessor(im.length,im[0].length);
            for (int i=0;i<im.length;i++){
                for (int ii=0;ii<im[0].length;ii++){
                    ip.putPixelValue(i, ii, im[i][ii]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    
    
    
    
    
    
    
    public static void imshow(double [] im,int width,int height){
        imshow(im, width, height, "");
    }
    public static void imshow(double [] im,int width,int height,String title){
        
            ImageProcessor ip = new FloatProcessor(width,height);
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    ip.putPixelValue(i, ii, (float)im[ii*width+i]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    public static void imshow(float [] im,int width,int height){
        imshow(im, width, height, "");
    }
    public static void imshow(float [] im,int width,int height,String title){
        
            ImageProcessor ip = new FloatProcessor(width,height);
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    ip.putPixelValue(i, ii, (float)im[ii*width+i]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    public static void imshow(int [] im,int width,int height){
        imshow(im, width, height, "");
    }
    public static void imshow(int [] im,int width,int height,String title){
        
            ImageProcessor ip = new FloatProcessor(width,height);
            for (int i=0;i<width;i++){
                for (int ii=0;ii<height;ii++){
                    ip.putPixelValue(i, ii, (float)im[ii*width+i]);
                }
            }
        
        ImagePlus imp=new ImagePlus(""+title,ip);
        imp.show();
    }
    
    
    
    
    
    
    
    
}
