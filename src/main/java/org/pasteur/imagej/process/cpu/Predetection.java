/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.process.cpu;
import org.pasteur.imagej.utils.FastFourierTransform;
import ij.IJ;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
/**
 *
 * @author benoit
 */
public class Predetection {
    int width;
    int height;
    int sizeFullImageX;
    int sizeFullImageY;
    
    
    int sizeFullImagePadX;
    int sizeFullImagePadY;
    
    int sizeFullImageFFTX;
    int sizeFullImageFFTY;
    
    DataPhase dparam;
    double threshold;
    double zstep;
    double xystep;
    
    
    double [] varPSF;
    
    double [][] fullImage;
    double [][] meanImage;
    double [][] sumImage;
    double [][] sumSquareImage;
    double [][] sumImageSecond;
    double [][] normImage;
    double [][] sumSquareImageSecond;
    double [][] varImage;
    
    double [][][] dictionnaryR;
    double [][][] dictionnaryI;
    
    double [][][] resultCorrelation;
    
    int sizePSF;
    
    public double [] range;
    
    double [][][] psf;
    
    FastFourierTransform fft;
    FastFourierTransform fft2;
    
    double thresholdCrossCorrelation;
    
    
    public Predetection(int sizeFullImageX,int sizeFullImageY, DataPhase dparam,double mini, double maxi, double step,double thresholdCrossCorrelation){
        
        this.thresholdCrossCorrelation=thresholdCrossCorrelation;
        this.sizeFullImageX=sizeFullImageX;
        this.sizeFullImageY=sizeFullImageY;
        sizeFullImagePadX=sizeFullImageX*2;
        sizeFullImagePadY=sizeFullImageY*2;
        this.dparam=dparam;
        threshold=thresholdCrossCorrelation;
        
        
        
        zstep=step;
        xystep=dparam.param.xystep;
        
        
        
        meanImage=new double [sizeFullImageX][sizeFullImageY];
        sumImage=new double [sizeFullImageX][sizeFullImageY];
        sumSquareImage=new double [sizeFullImageX][sizeFullImageY];
        sumImageSecond=new double [sizeFullImageX][sizeFullImageY];
        normImage=new double [sizeFullImageX][sizeFullImageY];
        sumSquareImageSecond=new double [sizeFullImageX][sizeFullImageY];
        varImage=new double [sizeFullImageX][sizeFullImageY];
        fullImage=new double [sizeFullImageX][sizeFullImageY];
        
        
        
        sizePSF=dparam.param.sizeoutput;
        
        
        for (int i=0;i<sizeFullImageX;i++){
            for (int ii=0;ii<sizeFullImageY;ii++){
                int valx1=sizePSF/2;
                int valx2=sizePSF/2;
                if (i<sizePSF/2-1){
                    valx1=i+1;
                }
                if (i>sizeFullImageX-sizePSF/2-2){
                    valx2=sizeFullImageY-i-1;
                }
                int valy1=sizePSF/2;
                int valy2=sizePSF/2;
                if (ii<sizePSF/2-1){
                    valy1=ii+1;
                }
                if (ii>sizeFullImageX-sizePSF/2-2){
                    valy2=sizeFullImageY-ii-2;
                }
                normImage[i][ii]=(valx1+valx2)*(valy1+valy2);
            }
        }
        
        
        
        int number=(int)Math.ceil((maxi-mini)/step);
        range = new double [number];
        double z=mini;
        
        for (int i=0;i<range.length;i++){
            range[i]=z;
            z+=step;
        }
        
        
        fft = new FastFourierTransform(sizeFullImagePadX,sizeFullImagePadY);
        fft2 = new FastFourierTransform(sizeFullImagePadX,sizeFullImagePadY);
        
        sizeFullImageFFTX=fft.getWidth();
        sizeFullImageFFTY=fft.getHeight();
                
        dictionnaryR=new double [range.length][sizeFullImageFFTX][sizeFullImageFFTY];
        dictionnaryI=new double [range.length][sizeFullImageFFTX][sizeFullImageFFTY];
        resultCorrelation=null;
        
        
        
        makeDictionnary();
        
         
         
         
    }
    
    
    
    
    
    
    
    
    private void makeDictionnary(){
        
        
//        int [] sparseIndexPaddingShift2DEvenDic = new int[sizeFullImage*sizeFullImage];
//        for (int i=0;i<dparam.param.sizeoutput;i++){
//            for (int ii=0;ii<dparam.param.sizeoutput;ii++){
//                
//                int x=((i+sizeFullImageFFT/2-dparam.param.sizeoutput/2));
//                int y=((ii+sizeFullImageFFT/2-dparam.param.sizeoutput/2));
//                int p=(((x+sizeFullImageFFT/2)%sizeFullImageFFT)*sizeFullImageFFT+((y+sizeFullImageFFT/2)%sizeFullImageFFT));
//                sparseIndexPaddingShift2DEvenDic[i*dparam.param.sizeoutput+ii]=p*2;
//                
//            }
//        }
        
//        Pointer device_sparseIndexPaddingShift2DEvenDic;
//        device_sparseIndexPaddingShift2DEvenDic=new Pointer();
//        cudaResult =cudaMalloc(device_sparseIndexPaddingShift2DEvenDic, dparam.param.sizeoutput*dparam.param.sizeoutput * Sizeof.INT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 6 "+cudaResult);}
//        cudaResult =cudaMemcpyAsync(device_sparseIndexPaddingShift2DEvenDic, Pointer.to(sparseIndexPaddingShift2DEvenDic), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.INT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 7 "+cudaResult);}
//        
//        
//        Pointer device_dic;
//        device_dic = new Pointer();
//        cudaResult =cudaMalloc(device_dic, dparam.param.sizeoutput*dparam.param.sizeoutput*range.length * Sizeof.FLOAT);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 182 "+cudaResult);}
        
        
        varPSF=new double [range.length];
//        cudaResult =JCuda.cudaMemsetAsync(device_dicFFT, 0, this.totalSize*2*range.length * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));
        
        double [][] psf_d = new double [dparam.param.sizeoutput][dparam.param.sizeoutput];
        psf= new double [range.length][dparam.param.sizeoutput][dparam.param.sizeoutput];
        
        for (int i=0;i<range.length;i++){
           
            
            
            dparam.psf.computePSF(0,0,dparam.param.Zfocus,range[i]);
            
            //conversion float
            double [][] pp=dparam.psf.getPSFPointer();
            for (int u=0;u<sizePSF;u++){
                for (int uu=0;uu<sizePSF;uu++){
                    psf_d[u][uu]=pp[u][uu];
                }
            }
            double mean=0;
            for (int u=0;u<sizePSF;u++){
                for (int uu=0;uu<sizePSF;uu++){
                    psf[i][u][uu]=psf_d[u][uu];
                    mean+=psf_d[u][uu];
                }
            }
            mean/=(double)sizePSF*sizePSF;
            varPSF[i]=0;
            for (int u=0;u<sizePSF;u++){
                for (int uu=0;uu<sizePSF;uu++){
                    psf_d[u][uu]-=mean;
                    varPSF[i]+=psf_d[u][uu]*psf_d[u][uu];
                }
            }
            varPSF[i]/=(double)(sizePSF*sizePSF);
            
            for (int u=0;u<sizePSF;u++){
                for (int uu=0;uu<sizePSF;uu++){
                    int x=(u-sizePSF/2+sizeFullImageFFTX/2+sizeFullImageFFTX/2)%sizeFullImageFFTX;
                    int y=(uu-sizePSF/2+sizeFullImageFFTY/2+sizeFullImageFFTY/2)%sizeFullImageFFTY;
                    dictionnaryR[i][x][y]=psf_d[u][uu];
                    dictionnaryI[i][x][y]=0;
                }
            }
            
            
            
            fft.setReal(dictionnaryR[i]);
            fft.setImag(dictionnaryI[i]);
            fft.fft2D();
            double [][] a=fft.getPointerRealOut2D();
            double [][] b=fft.getPointerImagOut2D();
            for (int u=0;u<sizeFullImageFFTX;u++){
                for (int uu=0;uu<sizeFullImageFFTY;uu++){
                    dictionnaryR[i][u][uu]=a[u][uu];
                    dictionnaryI[i][u][uu]=b[u][uu];
                }
            }
            //ImageShow.imshow(dictionnaryR[i], "dictionnaryR_"+range[i]);
            //ImageShow.imshow(dictionnaryI[i], "dictionnaryI_"+range[i]);
            //ImageShow.imshow(dictionnaryR, "dictionnaryR");
            //ImageShow.imshow(dictionnaryI, "dictionnaryI");
            
//            cudaResult =cudaMemcpyAsync(device_dic, Pointer.to(psf_f), dparam.param.sizeoutput*dparam.param.sizeoutput*Sizeof.FLOAT, cudaMemcpyHostToDevice,MyCudaStream.getCudaStream_t(streamId));
//            
//            
//            cudaResult =cusparseSsctr(MyCudaStream.getHandleCuSparse(this.streamId), dparam.param.sizeoutput*dparam.param.sizeoutput, device_dic, device_sparseIndexPaddingShift2DEvenDic,device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT), CUSPARSE_INDEX_BASE_ZERO);
//            
//            
//            cudaResult =JCufft.cufftExecC2C(plan, device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),device_dicFFT.withByteOffset(i*this.totalSize*2*Sizeof.FLOAT),JCufft.CUFFT_FORWARD);if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR cuda 28 "+cudaResult);}

        
        }
        
//        JCuda.cudaFree(device_dic) ;
//        JCuda.cudaFree(device_sparseIndexPaddingShift2DEvenDic) ;
//        MyVecDouble.divScalarFloat(MyCudaStream.getCUstream(streamId), this.totalSize*2*this.range.length, device_dicFFT, device_dicFFT, (float)Math.sqrt(totalSize));
        
        
        
    }
    
    
    public double [][][] getPSFNonNormalized(){
        double [][][] psf = new double [this.psf.length][this.psf[0].length][this.psf[0][0].length];
        for (int i=0;i<psf.length;i++){
            for (int ii=0;ii<this.psf[0].length;ii++){
                for (int iii=0;iii<this.psf[0][0].length;iii++){
                    psf[i][ii][iii]=this.psf[i][ii][iii];
                }
            }
        }
        return psf;
    }
    
    
    public double [] getRange(){
        double [] range = new double [this.range.length];
        for (int i=0;i<this.range.length;i++){
            
                range[i]=this.range[i];
            
        }
        return range;
    }
        
    
    
    
    
    
    
    public boolean setImage(double [][] image){
        
        this.width=image.length;
        this.height=image[0].length;
        
        //cudaResult=JCuda.cudaMemsetAsync(device_fullImageFFT, 0, this.totalSize*2 * Sizeof.FLOAT,MyCudaStream.getCudaStream_t(streamId));if (cudaResult != cudaError.cudaSuccess){IJ.log("ERROR memset cuda 0");}
        
        if ((width!=sizeFullImageX)||(height!=sizeFullImageY)){
            IJ.log("OOPS: problems are coming: image size should be set in PreDetection "+width+"  "+height+"  "+sizeFullImageX+"  "+sizeFullImageY);
            return false;
        }
        
        double [][] realPointer=fft.getPointerReal2D();
        double [][] imagPointer=fft.getPointerImag2D();
        for (int i=0,j=0;i<sizeFullImageFFTX;i++){
            for (int ii=0;ii<sizeFullImageFFTY;ii++){
                realPointer[i][ii]=0;
                imagPointer[i][ii]=0;
            }
        }
        int x,y;
        for (int i=0;i<sizeFullImageX;i++){
            for (int ii=0;ii<sizeFullImageY;ii++){
//                if ((i<width)&&(ii<height)){
                
                fullImage[i][ii]=(float)image[i][ii];
                x=(i-sizeFullImageX/2+sizeFullImageFFTX/2+sizeFullImageFFTX/2)%sizeFullImageFFTX;
                y=(ii-sizeFullImageY/2+sizeFullImageFFTY/2+sizeFullImageFFTY/2)%sizeFullImageFFTY;
                realPointer[x][y]=image[i][ii];
                
//                }
//                else{
//                    fullImage[i][ii]=0;
//                }
            }
        }
        
        fft.fft2D();
        
        
        
       //double [][] tmp0 = new double [sizeFullImage][sizeFullImage];
       //double [][] tmp1 = new double [sizeFullImage][sizeFullImage];
       //double [][] tmp2 = new double [sizeFullImage][sizeFullImage];
       
        for (int i=0;i<sizeFullImageX;i++){
            sumImage[i][0]=0;
            sumSquareImage[i][0]=0;
            for (int t=0;t<sizePSF/2+1;t++){//before 0 it is outside
                if (t<sizeFullImageY){
                    sumImage[i][0]+=fullImage[i][t];
                    sumSquareImage[i][0]+=fullImage[i][t]*fullImage[i][t];
                }
            }
        
            for (int ii=1;ii<sizeFullImageY;ii++){
                sumImage[i][ii]=sumImage[i][ii-1];
                if (ii-sizePSF/2>=0){
                    sumImage[i][ii]-=fullImage[i][ii-sizePSF/2];
                }
                else{
                    //outside
                }
                if (ii+sizePSF/2<sizeFullImageY){
                    sumImage[i][ii]+=fullImage[i][ii+sizePSF/2];
                }
                else{
                    //outside
                }
                sumSquareImage[i][ii]=sumSquareImage[i][ii-1];
                if (ii-sizePSF/2>=0){
                    sumSquareImage[i][ii]-=fullImage[i][ii-sizePSF/2]*fullImage[i][ii-sizePSF/2];
                }
                if (ii+sizePSF/2<sizeFullImageY){
                    sumSquareImage[i][ii]+=fullImage[i][ii+sizePSF/2]*fullImage[i][ii+sizePSF/2];
                }
                
            }
        }
        
        for (int ii=0;ii<sizeFullImageY;ii++){
            sumImageSecond[0][ii]=0;
            sumSquareImageSecond[0][ii]=0;
            for (int t=0;t<sizePSF/2+1;t++){
                if (t<sizeFullImageX){
                    sumImageSecond[0][ii]+=sumImage[t][ii];
                    sumSquareImageSecond[0][ii]+=sumSquareImage[t][ii];
                }
            }
            for (int i=1;i<sizeFullImageX;i++){
                sumImageSecond[i][ii]=sumImageSecond[(i-1)][ii];
                if (i-sizePSF/2>=0){
                    sumImageSecond[i][ii]-=sumImage[(i-sizePSF/2)][ii];
                }
                if (i+sizePSF/2<sizeFullImageX){
                    sumImageSecond[i][ii]+=sumImage[(i+sizePSF/2)][ii];
                }
                
                sumSquareImageSecond[i][ii]=sumSquareImageSecond[(i-1)][ii];
                if (i-sizePSF/2>=0){
                    sumSquareImageSecond[i][ii]-=sumSquareImage[(i-sizePSF/2)][ii];
                }
                if (i+sizePSF/2<sizeFullImageX){
                    sumSquareImageSecond[i][ii]+=sumSquareImage[(i+sizePSF/2)][ii];
                }
                
            }
            for (int i=0;i<sizeFullImageX;i++){
                meanImage[i][ii]=sumImageSecond[i][ii]/normImage[i][ii];
            }
        }
        
        
        
        for (int i=0;i<sizeFullImageX;i++){
            for (int ii=0;ii<sizeFullImageY;ii++){
                varImage[i][ii]=(sumSquareImageSecond[i][ii]/normImage[i][ii])-meanImage[i][ii]*meanImage[i][ii];
            }
        }
        //ImageShow.imshow(varImage);
        //for (int ii=0;ii<sizeFullImage;ii++){
        //    for (int i=0;i<sizeFullImage;i++){
        //        tmp0[i][ii]=meanImage[i][ii];
        //        tmp1[i][ii]=varImage[i][ii];
        //        tmp2[i][ii]=normImage[i][ii];
        //    }
        //}
        //classes.ImageShow.imshow(tmp0,"mean_image");
        //classes.ImageShow.imshow(tmp1,"var_image");
        //classes.ImageShow.imshow(tmp2,"norm_image");
        //dparam.psf.imshow(sizeFullImage, sumSquareImageSecond, "sumSquareImageSecond");
        
        
       
       
        return true;
    }
    
    
    
    
    
    
    
    public double [][] convolveNormalizedInFourierDomain(){
        if (resultCorrelation==null){
            resultCorrelation=new double [range.length][sizeFullImageX][sizeFullImageY];
        }
        return convolveNormalizedInFourierDomain(this.resultCorrelation);
    }
    
    public double [][] convolveNormalizedInFourierDomain(double [][][] resultCorrelation){
        double [][] inputR=fft.getPointerRealOut2D();
        double [][] inputI=fft.getPointerImagOut2D();
        double [][] inputR2=fft2.getPointerReal2D();
        double [][] inputI2=fft2.getPointerImag2D();
        double [][] outputR2=fft2.getPointerRealOut2D();
        double [][] outputI2=fft2.getPointerImagOut2D();
        int x,y;
        
        for (int i=0;i<range.length;i++){
            for (int ii=0;ii<sizeFullImageFFTX;ii++){
                for (int iii=0;iii<sizeFullImageFFTY;iii++){
                    //complexeConjugate
                    inputR2[ii][iii]=inputI[ii][iii]*dictionnaryI[i][ii][iii]+inputR[ii][iii]*dictionnaryR[i][ii][iii];
                    inputI2[ii][iii]=inputI[ii][iii]*dictionnaryR[i][ii][iii]-inputR[ii][iii]*dictionnaryI[i][ii][iii];
                }
            }
            fft2.ifft2D();
            
            for (int ii=0;ii<sizeFullImageX;ii++){
                for (int iii=0;iii<sizeFullImageY;iii++){
                    x=(ii-sizeFullImageX/2+sizeFullImageFFTX/2+sizeFullImageFFTX/2)%sizeFullImageFFTX;
                    y=(iii-sizeFullImageY/2+sizeFullImageFFTY/2+sizeFullImageFFTY/2)%sizeFullImageFFTY;
                    
                    if (Math.sqrt(varPSF[i]*varImage[ii][iii])>0){
                        resultCorrelation[i][ii][iii]=outputR2[x][y]/Math.sqrt(varPSF[i]*varImage[ii][iii]);
                    }
                    else if (varPSF[i]<=0){
                        resultCorrelation[i][ii][iii]=-1;
                    }
                    else if (varImage[ii][iii]<=0){
                        resultCorrelation[i][ii][iii]=-2;
                    }
                }
            }
        }
        
        
        
        
        
        //ImageShow.imshow(resultCorrelation, "cross");
        
        
        
        
        double [][] vect= getVectMaxPosition(resultCorrelation);
        
        Arrays.sort(vect, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return ((Double) o2[0]).compareTo(o1[0]);
            }
        });
        return vect;
    }
    
    
    
    
    
    
        
        double [][] getVectMaxPosition(double [][][] resultConvolution){


            
            int decalFilter=3;
            ArrayList<double []> al = new ArrayList<double []>();
            double [][] tmp = new double [width][height];
            for (int r=0;r<resultConvolution.length;r++){
                for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                    for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                        if ((resultConvolution[r][i][j]>0)){
                            tmp[i][j]=Double.NEGATIVE_INFINITY;
                            for (int a=-decalFilter;a<=decalFilter;a++){
                                tmp[i][j]=Math.max(resultConvolution[r][i][j+a],tmp[i][j]);
                            }
                        }
                    }
                }
                double max=0;
                //double min=0;
                for (int j=Math.max(decalFilter,sizePSF/2);j<height-Math.max(decalFilter,sizePSF/2);j++){
                    for (int i=Math.max(decalFilter,sizePSF/2);i<width-Math.max(decalFilter,sizePSF/2);i++){
                        if ((resultConvolution[r][i][j]>thresholdCrossCorrelation)){
                            
                            if (resultConvolution[r][i][j]==tmp[i][j]){
                                max=Double.NEGATIVE_INFINITY;
                                for (int a=-decalFilter;a<=decalFilter;a++){
                                    max=Math.max(tmp[i+a][j],max);
                                }
                                if (max==resultConvolution[r][i][j]){
                                    
//                                    min=Double.POSITIVE_INFINITY;
//                                    for (int a=-decalFilter;a<=decalFilter;a++){
//                                        for (int aa=-decalFilter;aa<=decalFilter;aa++){
//                                            min=Math.min(resultConvolution[r][i+a][j+aa],min);
//                                        }
//                                    }
//                                    if((min>0)){
                                        double [] machin = new double[7];
                                        
                                        double topZ=0;
                                        if ((r>1)&&(r<resultConvolution.length-1)){
                                            topZ=topOfParabola(-1.*zstep,resultConvolution[r-1][i][j],0,resultConvolution[r][i][j],1.*zstep,resultConvolution[r+1][i][j]);
                                        }
                                        
                                        double topX=0;
                                        if ((i>1)&&(i<resultConvolution[0].length-1)){
                                            topX=topOfParabola(-1.*xystep,resultConvolution[r][i-1][j],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i+1][j]);
                                        }
                                        
                                        double topY=0;
                                        if ((j>1)&&(r<resultConvolution[0][0].length-1)){
                                            topY=topOfParabola(-1.*xystep,resultConvolution[r][i][j-1],0,resultConvolution[r][i][j],1.*xystep,resultConvolution[r][i][j+1]);
                                        }
                                        
                                        machin[0]=resultConvolution[r][i][j];
                                        machin[1]=r;
                                        machin[2]=i;
                                        machin[3]=j;
                                        machin[4]=topZ;
                                        machin[5]=topX;
                                        machin[6]=topY;
                                        al.add(machin);
//                                    }
                                }
                                
                                
                            }
                        }
                    }
                }
            }
            
            double [][] vect = new double [al.size()][7];
            for (int i=0;i<al.size();i++){
                vect[i]=al.get(i);
            }
            
            al.clear();
            return vect;
        }
    
    
        private double topOfParabola(double xa,double ya,double xb,double yb,double xc,double yc){
            double a=(yc-ya)/((xc-xa)*(xc-xb))-(yb-ya)/((xb-xa)*(xc-xb));
            double b=((yb-ya)/(xb-xa))-a*(xb+xa);
            return(-b/(2*a));
        }
    
    
    
    
}
