/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;

/**
 *
 * @author benoit
 */
public class FastFourierTransform {
    
    int width;
    int height;
    int depth;
    int widthpow;
    int heightpow;
    int depthpow;
    
    double [] x ;
    double [] y ;
    double [] z ;
    
    int side=0;
    
    double [][][] real ;
    double [][][] imag;
    double [][][] realOut ;
    double [][][] imagOut;
    double [][][] realTmp ;
    double [][][] imagTmp;
    double [][][] realTmp2 ;
    double [][][] imagTmp2;
    
    
    
    double[] vrH; 
    double[] viH; 
    double[] vfrH; 
    double[] vfiH; 
    double[] vrW; 
    double[] viW; 
    double[] vfrW; 
    double[] vfiW; 
    double[] vrD; 
    double[] viD; 
    double[] vfrD; 
    double[] vfiD; 
    
    boolean _2D;
    
    FastFT1D fft1d = new FastFT1D();
    
    public FastFourierTransform(int width,int height){
        
        _2D=true;
        
        this.width=1;
        this.height=width;
        this.depth=height;
        this.widthpow=1;
        this.heightpow=(int)Math.pow(2, (int)Math.ceil(Math.log(width)/Math.log(2)));
        this.depthpow=(int)Math.pow(2, (int)Math.ceil(Math.log(height)/Math.log(2)));
        
        real= new double [widthpow][heightpow][depthpow];
        realOut= new double [widthpow][heightpow][depthpow];
        imag= new double [widthpow][heightpow][depthpow];
        imagOut= new double [widthpow][heightpow][depthpow];
        
        realTmp= new double [widthpow][heightpow][depthpow];
        imagTmp= new double [widthpow][heightpow][depthpow];
        
        
        x = new double[widthpow*2+1];
        y = new double[heightpow*2+1];
        z = new double[depthpow*2+1];
    
        vrW = null; 
        viW = null; 
        vfrW = null; 
        vfiW = null; 
        
        vrH = new double[heightpow]; 
        viH = new double[heightpow]; 
        vfrH = new double[heightpow]; 
        vfiH = new double[heightpow]; 
        
        vrD = new double[depthpow]; 
        viD = new double[depthpow]; 
        vfrD = new double[depthpow]; 
        vfiD = new double[depthpow]; 
        
        
    }
    
    
    public int getWidth(){
        if (_2D){
            return heightpow;
        }
        else{
            return widthpow;
        }
    }
    
    public int getHeight(){
        if (_2D){
            return depthpow;
        }
        else{
            return heightpow;
        }
    }
    
    public int getDepth(){
        if (_2D){
            //IJ.log("sorry, depth does not exists for 2D fft");
            return 1;
        }
        else{
            return depthpow;
        }
    }
    
    public FastFourierTransform(int width,int height,int depth){
        
        _2D=false;
        
        this.width=width;
        this.height=height;
        this.depth=depth;
        
        this.widthpow=(int)Math.pow(2, (int)Math.ceil(Math.log(width)/Math.log(2)));
        this.heightpow=(int)Math.pow(2, (int)Math.ceil(Math.log(height)/Math.log(2)));
        this.depthpow=(int)Math.pow(2, (int)Math.ceil(Math.log(depth)/Math.log(2)));
        
        real= new double [widthpow][heightpow][depthpow];
        realOut= new double [widthpow][heightpow][depthpow];
        imag= new double [widthpow][heightpow][depthpow];
        imagOut= new double [widthpow][heightpow][depthpow];
        
        realTmp= new double [widthpow][heightpow][depthpow];
        imagTmp= new double [widthpow][heightpow][depthpow];
        
        realTmp2= new double [widthpow][heightpow][depthpow];
        imagTmp2= new double [widthpow][heightpow][depthpow];
        
        
        x = new double[widthpow*2+1];
        y = new double[heightpow*2+1];
        z = new double[depthpow*2+1];
        
        vrH = new double[heightpow]; 
        viH = new double[heightpow]; 
        vfrH = new double[heightpow]; 
        vfiH = new double[heightpow]; 
        
        vrW = new double[widthpow]; 
        viW = new double[widthpow]; 
        vfrW = new double[widthpow]; 
        vfiW = new double[widthpow];
        
        vrD = new double[depthpow]; 
        viD = new double[depthpow]; 
        vfrD = new double[depthpow]; 
        vfiD = new double[depthpow];
        
        
    }
    
    
    
    
    public void setReal(double [][] rIn){
        
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                real[0][i][ii]=rIn[i][ii];
            }
        }
        for (int i=rIn.length;i<heightpow;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                real[0][i][ii]=0;
            }
        }
        for (int i=0;i<heightpow;i++){
            for (int ii=rIn[0].length;ii<depthpow;ii++){
                real[0][i][ii]=0;
            }
        }
        
        
        
    }
    public void setImag(double [][] iIn){
        
        for (int i=0;i<iIn.length;i++){
            for (int ii=0;ii<iIn[0].length;ii++){
                imag[0][i][ii]=iIn[i][ii];
            }
        }
        for (int i=iIn.length;i<heightpow;i++){
            for (int ii=0;ii<iIn[0].length;ii++){
                imag[0][i][ii]=0;
            }
        }
        for (int i=0;i<heightpow;i++){
            for (int ii=iIn[0].length;ii<depthpow;ii++){
                imag[0][i][ii]=0;
            }
        }
        
    }
    
    
    
    
    
    
    
    
    public double [][] getPointerImag2D(){
        return imag[0]; 
    }
    
    public double [][] getPointerReal2D(){
        return real[0]; 
    }
    
    public double [][] getPointerImagOut2D(){
        return imagOut[0]; 
    }
    
    public double [][] getPointerRealOut2D(){
        return realOut[0]; 
    }
    
    
    
    
    
    
    
    
    
    
    
    public void setReal(double [][][] rIn){
        
        
        for (int i=0;i<rIn.length;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                for (int iii=0;iii<rIn[0][0].length;iii++){
                    real[i][ii][iii]=rIn[i][ii][iii];
                }
            }
        }
        for (int i=rIn.length;i<widthpow;i++){
            for (int ii=0;ii<rIn[0].length;ii++){
                for (int iii=0;iii<rIn[0][0].length;iii++){
                    real[i][ii][iii]=0;
                }
            }
        }
        for (int i=0;i<widthpow;i++){
            for (int ii=rIn[0].length;ii<heightpow;ii++){
                for (int iii=0;iii<rIn[0][0].length;iii++){
                    real[i][ii][iii]=0;
                }
            }
        }
        for (int i=0;i<widthpow;i++){
            for (int ii=0;ii<heightpow;ii++){
                for (int iii=rIn[0][0].length;iii<depthpow;iii++){
                    real[i][ii][iii]=0;
                }
            }
        }
        
    }
    public void setImag(double [][][] iIn){
        
        for (int i=0;i<iIn.length;i++){
            for (int ii=0;ii<iIn[0].length;ii++){
                for (int iii=0;iii<iIn[0][0].length;iii++){
                    imag[i][ii][iii]=iIn[i][ii][iii];
                }
            }
        }
        for (int i=iIn.length;i<widthpow;i++){
            for (int ii=0;ii<iIn[0].length;ii++){
                for (int iii=0;iii<iIn[0][0].length;iii++){
                    imag[i][ii][iii]=0;
                }
            }
        }
        for (int i=0;i<widthpow;i++){
            for (int ii=iIn[0].length;ii<heightpow;ii++){
                for (int iii=0;iii<iIn[0][0].length;iii++){
                    imag[i][ii][iii]=0;
                }
            }
        }
        for (int i=0;i<widthpow;i++){
            for (int ii=0;ii<heightpow;ii++){
                for (int iii=iIn[0][0].length;iii<depthpow;iii++){
                    imag[i][ii][iii]=0;
                }
            }
        }
        
    }
    
    
    
    
    
    
    
    
    public double [][][] getPointerImag3D(){
        return imag; 
    }
    
    public double [][][] getPointerReal3D(){
        return real; 
    }
    
    public double [][][] getPointerImagOut3D(){
        return imagOut; 
    }
    
    public double [][][] getPointerRealOut3D(){
        return realOut; 
    }
    
    
    
    
    public void fft2D(){
        run1(true);
        
        realTmp2=realTmp;
        imagTmp2=imagTmp;
        run3(true);
        
        
    }
    
    
    public void ifft2D(){
        
        run1(false);
        
        realTmp2=realTmp;
        imagTmp2=imagTmp;
        run3(false);
        
    }
    
    
    
    public void fft3D(){
        
        
        run1(true);
        
        run2(true);
        
        run3(true);
        
    }
    
    
    public void ifft3D(){
        
        
        run1(false);
        
        run2(false);
        
        run3(false);
        
    }
    
    
    







    
    //first part
    private void run1(boolean direct){


        for (int ii = 0; ii < heightpow; ii++){
            vfrH[ii]=0;
            vfiH[ii]=0;
        }
        
        if (direct){
            for (int i=0;i<widthpow;i++){
                for (int iii = 0; iii < depthpow; iii++){
                    for (int ii = 0; ii < heightpow; ii++){
                        vrH[ii] = (real[i][ii][iii]);
                        //IJ.log("realIn "+realIn[i][ii]);
                        viH[ii] = (imag[i][ii][iii]);
                    }
                    fft1d.fft(vrH,viH,vfrH,vfiH,1); 
                    for (int n = 0; n < heightpow; n++){
                            realTmp[i][n][iii] = (vfrH[n]);
                            imagTmp[i][n][iii] = (vfiH[n]);
                    }
                }
            }


        }
        else{
            for (int i=0;i<widthpow;i++){
                for (int iii = 0; iii < depthpow; iii++){
                    for (int ii = 0; ii < heightpow; ii++){
                        vrH[ii] = (real[i][ii][iii]);
                        //IJ.log("realIn "+realIn[i][ii]);
                        viH[ii] = (imag[i][ii][iii]);
                    }
                    fft1d.ifft(vrH,viH,vfrH,vfiH,1); 
                    for (int n = 0; n < heightpow; n++){
                            realTmp[i][n][iii] = (vfrH[n]);
                            imagTmp[i][n][iii] = (vfiH[n]);
                    }
                }

            }


        }
    }


    //second part
    private void run2(boolean direct){

        for (int ii = 0; ii < widthpow; ii++){
            vfrW[ii]=0;
            vfiW[ii]=0;
        }
        
        if (direct){
            
            for (int i=0;i<heightpow;i++){
                for (int iii = 0; iii < depthpow; iii++){
                    for (int ii = 0; ii < widthpow; ii++){
                        vrW[ii] = (realTmp[ii][i][iii]);
                        viW[ii] = (imagTmp[ii][i][iii]);
                    }
                    fft1d.fft(vrW,viW,vfrW,vfiW,0); 
                    for (int n = 0; n < widthpow; n++){
                            realTmp2[n][i][iii] = (vfrW[n]);
                            imagTmp2[n][i][iii] = (vfiW[n]);
                    }
                }
            }
            
        }
        else{
            for (int i=0;i<heightpow;i++){
                for (int iii = 0; iii < depthpow; iii++){
                    for (int ii = 0; ii < widthpow; ii++){
                        vrW[ii] = (realTmp[ii][i][iii]);
                        viW[ii] = (imagTmp[ii][i][iii]);
                    }
                    fft1d.ifft(vrW,viW,vfrW,vfiW,0); 
                    for (int n = 0; n < widthpow; n++){
                            realTmp2[n][i][iii] = (vfrW[n]);
                            imagTmp2[n][i][iii] = (vfiW[n]);
                    }
                }
            }
        }
    }
    
    
    
    
    
    //third part
    private void run3(boolean direct){

        for (int ii = 0; ii < depthpow; ii++){
            vfrD[ii]=0;
            vfiD[ii]=0;
        }
        
        if (direct){
            for (int i=0;i<heightpow;i++){
                for (int ii = 0; ii < widthpow; ii++){
                    for (int iii = 0; iii < depthpow; iii++){
                        vrD[iii] = (realTmp2[ii][i][iii]);
                        viD[iii] = (imagTmp2[ii][i][iii]);
                    }
                    fft1d.fft(vrD,viD,vfrD,vfiD,2); 
                    for (int n = 0; n < depthpow; n++){
                            realOut[ii][i][n] = (vfrD[n]);
                            imagOut[ii][i][n] = (vfiD[n]);
                    }
                }
            }

        }
        else{
            for (int i=0;i<heightpow;i++){
                for (int ii = 0; ii < widthpow; ii++){
                    for (int iii = 0; iii < depthpow; iii++){
                        vrD[iii] = (realTmp2[ii][i][iii]);
                        viD[iii] = (imagTmp2[ii][i][iii]);
                    }
                    fft1d.ifft(vrD,viD,vfrD,vfiD,2); 
                    for (int n = 0; n < depthpow; n++){
                            realOut[ii][i][n] = (vfrD[n]);
                            imagOut[ii][i][n] = (vfiD[n]);
                    }
                }
            }
        }
    }
    
    
    
    
    
    






    
    
    class FastFT1D {
        double TWOPI=2*Math.PI;
	/** constructor for the use of the 1D transformation routines */ 
	public FastFT1D() {
	}
        
        
        
        /** constructor for the use of the 1D transformation routines */ 
	public void fft(double [] realIn,double [] imagIn,double [] realOut,double [] imagOut,int side) {
            
            boolean allZero=true;
            loop:for (int i=0;i<realIn.length;i++){
                if ((realIn[i]!=0)||(imagIn[i]!=0)){
                    allZero=false;
                    break loop;
                }
            }
            if (!allZero){
                this.transform(realIn,imagIn,realOut,imagOut,1,side);
            }
            else{
                for (int i=0;i<realIn.length;i++){
                    realOut[i]=0;
                    imagOut[i]=0;
                }
            }
	}
        
        /** constructor for the use of the 1D transformation routines */ 
	public void ifft(double [] realIn,double [] imagIn,double [] realOut,double [] imagOut,int side) {
            
            boolean allZero=true;
            loop:for (int i=0;i<realIn.length;i++){
                if ((realIn[i]!=0)||(imagIn[i]!=0)){
                    allZero=false;
                    break loop;
                }
            }
            if (!allZero){
                this.transform(realIn,imagIn,realOut,imagOut,-1,side);
            }
            else{
                for (int i=0;i<realIn.length;i++){
                    realOut[i]=0;
                    imagOut[i]=0;
                }
            }
	}
        
        
        /** iterative version of fft */  
	private void transform(double [] realIn,double [] imagIn,double [] realOut,double [] imagOut,int inverseFFT,int side){
            
                int N = realIn.length*2+1;
                
                double [] p;
                if (side==0){
                    for (int t=0;t<realIn.length;t++){
                        x[2*t+1]=realIn[t];
                        x[2*t+2]=imagIn[t];
                    }
                    fft(x,realIn.length,inverseFFT);
                    for (int t=0;t<realIn.length;t++){
                        realOut[t]=x[2*t+1]/Math.sqrt(realIn.length);
                        imagOut[t]=x[2*t+2]/Math.sqrt(realIn.length);
                    }
                }
                else if (side==1){
                    for (int t=0;t<realIn.length;t++){
                        y[2*t+1]=realIn[t];
                        y[2*t+2]=imagIn[t];
                    }
                    fft(y,realIn.length,inverseFFT);
                    for (int t=0;t<realIn.length;t++){
                        realOut[t]=y[2*t+1]/Math.sqrt(realIn.length);
                        imagOut[t]=y[2*t+2]/Math.sqrt(realIn.length);
                    }
                }
                else if (side==2){
                    for (int t=0;t<realIn.length;t++){
                        z[2*t+1]=realIn[t];
                        z[2*t+2]=imagIn[t];
                    }
                    fft(z,realIn.length,inverseFFT);
                    for (int t=0;t<realIn.length;t++){
                        realOut[t]=z[2*t+1]/Math.sqrt(realIn.length);
                        imagOut[t]=z[2*t+2]/Math.sqrt(realIn.length);
                    }
                }
                
                
                
                
                
                
                
	}
        
        
        // compute the FFT of x[], assuming its length is a power of 2
        public void fft(double [] data, int nn, int isign){
            int n, mmax, m, j, istep, i;
            double wtemp, wr, wpr, wpi, wi, theta;
            double tempr, tempi;

            n = nn << 1;
            j = 1;
            for (i = 1; i < n; i += 2) {
            if (j > i) {
                tempr = data[j];     
                data[j] = data[i];     
                data[i] = tempr;

                tempr = data[j+1]; 
                data[j+1] = data[i+1]; 
                data[i+1] = tempr;
            }
            m = n >> 1;
            while (m >= 2 && j > m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        }
        mmax = 2;
        while (n > mmax) {
            istep = 2*mmax;
            theta = -TWOPI/(isign*mmax);

            wtemp = Math.sin(0.5*theta);
            wpi = Math.sin(theta);

            wpr = -2.0*wtemp*wtemp;

            wr = 1.0;
            wi = 0.0;
            for (m = 1; m < mmax; m += 2) {
                for (i = m; i <= n; i += istep) {
                    j =i + mmax;
                    tempr = wr*data[j]   - wi*data[j+1];
                    tempi = wr*data[j+1] + wi*data[j];
                    data[j]   = data[i]   - tempr;
                    data[j+1] = data[i+1] - tempi;
                    data[i] += tempr;
                    data[i+1] += tempi;
                }
                wr = (wtemp = wr)*wpr - wi*wpi + wr;
                wi = wi*wpr + wtemp*wpi + wi;
            }
            mmax = istep;
        }





        }






    }
    
    
}

