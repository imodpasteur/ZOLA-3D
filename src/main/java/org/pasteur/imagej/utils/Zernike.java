/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;



import ij.IJ;

/**
 *
 * @author benoit
 */
public class Zernike {
    
    double radiusRing;
    int sizeImagePx;
    Factorial fact = new Factorial();
    public int [][] poly;
    public int [] complexity;//sum of rang and index (abs(m)+abs(n))
    public double [][][] Z;
    double [] W;
    public Zernike(int sizeImagePx, double radiusRing,int polyNumber){
        this.radiusRing=radiusRing;
        this.sizeImagePx=sizeImagePx;
        poly = new int[polyNumber][2];
        complexity = new int[polyNumber];
        
        
        int index=0;
        loop:for (int i=0,rang=0;i<polyNumber;rang++){
            //String s="";
            for (int u=index;u<=-index;u+=2){
                //s+=" "+u+"     ";
                poly[i][0]=rang;
                poly[i][1]=u;
                
                i++;
                if (i>=polyNumber){
                    break loop;
                }
            }
            index--;
            //IJ.log("u  "+s);
        }
        
        Z = new double [polyNumber][][];
        for (int i=0;i<polyNumber;i++){
            //IJ.log("  "+poly[0][i]+"   "+poly[1][i]);
            Z[i]=getZ(poly[i][0],poly[i][1]);
            complexity[i]=Math.abs(poly[i][0])+Math.abs(poly[i][1]);
            
            //if (i==polyNumber-1){
            //    Z[i]=getZRot();
            //}
            //IJ.log("myComp= "+i+"  "+complexity[i]+"  "+poly[i][0]+"  "+poly[i][1]);
        }
        
        //ImageShow.imshow(Z,"Z");
        
        
        
        
        
        /*IJ.log("WARNING: REMOVE FOLLOWING LINES IN ZERNIKE CLASS ************************************");
        for (int i=polyNumber-20,j=0;i<polyNumber;i+=2,j++){
            Z[i]=getART(0,j,false);
            Z[i+1]=getART(0,j,true);

        }
        ImageShow.imshow(Z,"Z");
        IJ.log("END WARNING: REMOVE PREVIOUS LINES IN ZERNIKE CLASS ************************************");*/
        
        
        
    }
    
    
    public Zernike(int sizeImagePx, double radiusRing,boolean iscomplex){
        
            int maxM=10;
            int maxN=4;
            int polyNumber=maxM*maxN;
            this.radiusRing=radiusRing;
            this.sizeImagePx=sizeImagePx;
            poly = new int[polyNumber][2];
            complexity = new int[polyNumber];


            int index=0;
            loop:for (int i=0,rang=0;rang<maxN;rang++){
                //String s="";
                for (int u=0;u<maxM;u++){
                    //s+=" "+u+"     ";
                    poly[i][0]=rang;
                    poly[i][1]=u;

                    i++;
                    if (i>=polyNumber){
                        break loop;
                    }
                }
                index--;
                //IJ.log("u  "+s);
            }

            Z = new double [polyNumber][][];
            for (int i=0;i<polyNumber;i++){
                //IJ.log("  "+poly[0][i]+"   "+poly[1][i]);
                Z[i]=getART(poly[i][0],poly[i][1],iscomplex);
                complexity[i]=Math.abs(poly[i][0])+Math.abs(poly[i][1]);

                //if (i==polyNumber-1){
                //    Z[i]=getZRot();
                //}
                //IJ.log("myComp= "+i+"  "+complexity[i]+"  "+poly[i][0]+"  "+poly[i][1]);
            }
            //ImageShow.imshow(Z,"Z");
        
        
        
        
    }
    
    
    
    public Zernike(int sizeImagePx, double radiusRing,int polyNumber,int polyNumberReversed){
        
        this.radiusRing=radiusRing;
        this.sizeImagePx=sizeImagePx;
        poly = new int[polyNumber+polyNumberReversed][2];
        complexity = new int[polyNumber];
        
        
        int index=0;
        loop:for (int i=0,rang=0;i<polyNumber;rang++){
            //String s="";
            for (int u=index;u<=-index;u+=2){
                //s+=" "+u+"     ";
                poly[i][0]=rang;
                poly[i][1]=u;
                complexity[i]=Math.abs(rang)+Math.abs(u);
                i++;
                if (i>=polyNumber){
                    break loop;
                }
            }
            index--;
            //IJ.log("u  "+s);
        }
        
        index=0;
        loop:for (int i=polyNumber,rang=0;i<polyNumber+polyNumberReversed;rang++){
            //String s="";
            for (int u=index;u<=-index;u+=2){
                //s+=" "+u+"     ";
                poly[i][0]=rang;
                poly[i][1]=u;
                i++;
                if (i>=polyNumber+polyNumberReversed){
                    break loop;
                }
            }
            index--;
            //IJ.log("u  "+s);
        }
        
        Z = new double [polyNumber+polyNumberReversed][][];
        for (int i=0;i<polyNumber;i++){
            //IJ.log("  "+poly[0][i]+"   "+poly[1][i]);
            Z[i]=getZ(poly[i][0],poly[i][1]);
            complexity[i]=Math.abs(poly[i][0])+Math.abs(poly[i][1]);
        }
        for (int i=polyNumber;i<polyNumberReversed+polyNumber;i++){
            //IJ.log("  "+poly[0][i]+"   "+poly[1][i]);
            Z[i]=getZReversed(poly[i][0],poly[i][1]);
        }
        ImageShow.imshow(Z,"Z");
        
        
        
        
    }
    
    
    public Zernike(int sizeImagePx, double radiusRing,int polyNumber,boolean added){
        
        this.radiusRing=radiusRing;
        this.sizeImagePx=sizeImagePx;
        poly = new int[polyNumber][2];
        complexity = new int[polyNumber];
        
        
        int index=0;
        loop:for (int i=0,rang=0;i<polyNumber;rang++){
            //String s="";
            for (int u=index;u<=-index;u+=2){
                //s+=" "+u+"     ";
                poly[i][0]=rang;
                poly[i][1]=u;
                complexity[i]=Math.abs(rang)+Math.abs(u);
                i++;
                if (i>=polyNumber){
                    break loop;
                }
            }
            index--;
            //IJ.log("u  "+s);
        }
        
        Z = new double [polyNumber+1][][];
        for (int i=0;i<polyNumber;i++){
            //IJ.log("  "+poly[0][i]+"   "+poly[1][i]);
            Z[i]=getZ(poly[i][0],poly[i][1]);
            complexity[i]=Math.abs(poly[i][0])+Math.abs(poly[i][1]);
        }
        for (int i=polyNumber;i<1+polyNumber;i++){
            //IJ.log("  "+poly[0][i]+"   "+poly[1][i]);
            Z[i]=getZRot();
        }
        ImageShow.imshow(Z,"Z");
        
        
        
        
    }
    
    public void fit(double [][] f){
        int number=0;
        for (int i=0;i<f.length;i++){
            for (int ii=0;ii<f[0].length;ii++){
                if (Z[0][i][ii]>.5){
                    number++;
                }
            }
        }
        
        double [][] A=new double [number][poly.length];
        double [] B=new double [number];
        W=new double [poly.length];
        int t=0;
        for (int i=0;i<f.length;i++){
            for (int ii=0;ii<f[0].length;ii++){
                if (Z[0][i][ii]>.5){
                    for (int p=0;p<poly.length;p++){
                        A[t][p]=Z[p][i][ii];
                    }
                    B[t]=f[i][ii];
                    t++;
                }
            }
        }
        try{
            Matrixe AA = new Matrixe(A);
            Matrixe BB = new Matrixe(B);
            Matrixe AAt = Matrixe.transpose(AA);
            Matrixe AA1 = (Matrixe.times(AAt, AA));
            Matrixe BB1 = (Matrixe.times(AAt, BB));
            double [][] myAA1 = AA1.getMatrixe();
            double [][] myBB1 = BB1.getMatrixe();
            double [] myBB2 = new double[myBB1.length];
            for (int i=0;i<myBB1.length;i++){
                myBB2[i]=myBB1[i][0];
            }
            W=GaussianElimination.lsolve(myAA1, myBB2);
            //Matrixe AA1 = Matrixe.times(AAinv, AAt);
            //Matrixe WW=Matrixe.times(AA1, BB);
            //myW = WW.getMatrixe();
        }catch(Exception ee){IJ.log("some problems in matrix inversion... "+ee);}
        
        for (int i=0;i<W.length;i++){
            IJ.log("W "+W[i]);
        }
        IJ.log("");
        
    }
    
    
    double [][] unwrap(double [][] phase){
        return phase;
    }
    
    
    public double [] fitCosSin(double [][] f, double [][] weight){
        W=new double [poly.length];
        int iter=20;
        for (int t=0;t<iter;t++){
            double lik=getLik(f,weight);
            IJ.log("iter fit "+(iter-t)+"  "+lik);
            
            for (int i=0;i<W.length;i++){
                lik=getLik(f,weight);
                double save=W[i];
                
                double grad=getGrad(f,weight,i);
                loop:for (double a=.1;a>.00001;a/=10){
                    W[i]-=a*grad;
                    double likNew=getLik(f,weight);
                    if (true){//(likNew<lik){
                        lik=likNew;
                        break loop;
                    }
                    else{
                        W[i]=save;
                    }
                }
                //IJ.log("iter fit "+(iter-t)+"   "+(W.length-i)+"  "+lik);
            }
        }
        double [][] res = this.sumZ();
        double [][] resc = this.sumZ();
        double [][] ress = this.sumZ();
        double [][] fc = this.sumZ();
        double [][] fs = this.sumZ();
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                resc[u][uu]=Math.cos(res[u][uu]);
                ress[u][uu]=Math.sin(res[u][uu]);
                fc[u][uu]=Math.cos(f[u][uu]);
                fs[u][uu]=Math.sin(f[u][uu]);
            }
        }
        ImageShow.imshow(f,"f");
        ImageShow.imshow(res,"res");
        ImageShow.imshow(fs,"fs");
        ImageShow.imshow(ress,"ress");
        ImageShow.imshow(fc,"fc");
        ImageShow.imshow(resc,"resc");
        
        for (int i=0;i<W.length;i++){
            IJ.log("weight ("+i+") = "+W[i]);
        }
        
        
        return W;
    }
    
    public double  getLik(double [][] f,double [][] weight){
        double [][] mod=new double [sizeImagePx][sizeImagePx];
        
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                mod[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    mod[u][uu]+=W[i]*Z[i][u][uu];
                }
            }
        }
        
        double som=0;
        double x,y;
        //double [][] im1= new double [mod.length][mod[0].length];
        //double [][] im2= new double [mod.length][mod[0].length];
        //double [][] imgrad= new double [mod.length][mod[0].length];
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                //im1[u][uu]=Math.cos(f[u][uu]);
                //im1[u][uu]=Math.cos(mod[u][uu]);
                x=(Math.cos(f[u][uu])-Math.cos(mod[u][uu]));
                y=(Math.sin(f[u][uu])-Math.sin(mod[u][uu]));
                
                som+=(x*x+y*y);
                //imgrad[u][uu]=weight[u][uu]*(x*x);
            }
        }
        
        //ImageShow.imshow(im1,"im1");
        //ImageShow.imshow(im2,"im2");
        //ImageShow.imshow(imgrad,"grad");
        
        return som;
        //ImageShow.imshow(im,"weighted sum of Zernike coef");
    }
     
    
    
    public double  getGrad(double [][] f,double [][] weight,int idW){
        double [][] mod=new double [sizeImagePx][sizeImagePx];
        double h=.0001;
        double save=W[idW];
        W[idW]-=h;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                mod[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    mod[u][uu]+=W[i]*Z[i][u][uu];
                }
            }
        }
        
        double som0=0;
        double x,y;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                
                x=(Math.cos(f[u][uu])-Math.cos(mod[u][uu]));
                y=(Math.sin(f[u][uu])-Math.sin(mod[u][uu]));
                som0+=(x*x+y*y);
            }
        }
        W[idW]+=h;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                mod[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    mod[u][uu]+=W[i]*Z[i][u][uu];
                }
            }
        }
        
        double som1=0;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                
                x=(Math.cos(f[u][uu])-Math.cos(mod[u][uu]));
                y=(Math.sin(f[u][uu])-Math.sin(mod[u][uu]));
                som1+=(x*x+y*y);
            }
        }
        
        
        W[idW]+=h;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                mod[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    mod[u][uu]+=W[i]*Z[i][u][uu];
                }
            }
        }
        
        double som2=0;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                
                x=(Math.cos(f[u][uu])-Math.cos(mod[u][uu]));
                y=(Math.sin(f[u][uu])-Math.sin(mod[u][uu]));
                som2+=(x*x+y*y);
            }
        }
        
        
        W[idW]=save;
        double grad=0;
        if (Math.abs((som2+som0-2*som1)/(h*h))==0){
            grad=((som2-som0)/(2*h));
        }
        else{
            grad=((som2-som0)/(2*h))/Math.abs((som2+som0-2*som1)/(h*h));
        }
        return grad;
        //ImageShow.imshow(im,"weighted sum of Zernike coef");
    }
    
    
    
    //this function does not work...error ?
    public double  getDerivative(double [][] f,int idW){
        double [][] mod=new double [sizeImagePx][sizeImagePx];
        
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                mod[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    mod[u][uu]+=W[i]*Z[i][u][uu];
                }
            }
        }
        
        double som=0;
        double x,y;
        double c,s;
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                c=Math.cos(mod[u][uu]);
                s=Math.sin(mod[u][uu]);
                x=2*(Math.cos(f[u][uu])-c)*(s*Z[idW][u][uu]);
                y=2*(Math.sin(f[u][uu])-s)*(-c*Z[idW][u][uu]);
                
                som+=x+y;
            }
        }
        
        
        return som;
        //ImageShow.imshow(im,"weighted sum of Zernike coef");
    }
    
    
    public double [][] sumZ(){
        double [][] im=new double [sizeImagePx][sizeImagePx];
        
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                im[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    im[u][uu]+=W[i]*Z[i][u][uu];
                }
            }
        }
        return im;
        //ImageShow.imshow(im,"weighted sum of Zernike coef");
    }
    
    
    
    public void sumZcos(){
        double [][] im=new double [sizeImagePx][sizeImagePx];
        
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                im[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    im[u][uu]+=W[i]*Math.cos(Z[i][u][uu]);
                }
            }
        }
        
        ImageShow.imshow(im,"cos weighted sum of Zernike coef");
    }
    
    public void sumZsin(){
        double [][] im=new double [sizeImagePx][sizeImagePx];
        
        for (int u=0;u<sizeImagePx;u++){
            for (int uu=0;uu<sizeImagePx;uu++){
                im[u][uu]=0;
                for (int i=0;i<W.length;i++){
                    im[u][uu]+=W[i]*Math.sin(Z[i][u][uu]);
                }
            }
        }
        
        ImageShow.imshow(im,"sin weighted sum of Zernike coef");
    }
    
    
    
    
    
    
    public double [][] getZ(int n, int m){
        
        double [][] image = new double [sizeImagePx][sizeImagePx];
        
        double centerX=((double)sizeImagePx)/2.-.5;
        double centerY=((double)sizeImagePx)/2.-.5;
        for (int i=0;i<sizeImagePx;i++){
            for (int ii=0;ii<sizeImagePx;ii++){
                double x=((double)(i-centerX)/radiusRing);
                double y=((double)(ii-centerY)/radiusRing);
                double theta=Math.atan2(y,x);
                double rho=Math.sqrt(x*x+y*y);
                if (rho<1){
                    image[i][ii]=this.getZ(n, m, theta, rho);
                    //IJ.log("im "+i+"  "+ii+"  "+image[i][ii]);
                }
            }
        }
        
        
        
        
        
        
        //ImageShow.imshow(image,"zernike "+n+"  "+m);
        
        
        return image;
        
    }
    
    
    
    
    public double getZ(int n, int m,double theta, double rho){
        if (n<Math.abs(m)){
            IJ.log("error in zernike coefficient : n<m");
            return -1;
        }
        else{
            if (m>=0){
                double som=0;
                for (int k=0;k<=(n-m)/2;k++){
                    double truc=(Math.pow(-1, k)*fact.get(n-k));
                    double trac=fact.get(k)*fact.get(((n+m)/2)-k)*fact.get(((n-m)/2)-k);
                    som+=(truc/trac)*Math.pow(rho,n-2*k);
                }
                return som*(Math.cos(theta*m));
            }
            else{
                m=-m;
                double som=0;
                for (int k=0;k<=(n-m)/2;k++){
                    double truc=(Math.pow(-1, k)*fact.get(n-k));
                    double trac=fact.get(k)*fact.get(((n+m)/2)-k)*fact.get(((n-m)/2)-k);
                    som+=(truc/trac)*Math.pow(rho,n-2*k);
                }
                return som*Math.sin(theta*m);
            }
        }
    }
    
    
    
    
    
    public double [][] getZRot(){
        
        double [][] image = new double [sizeImagePx][sizeImagePx];
        
        double centerX=((double)sizeImagePx)/2.-.5;
        double centerY=((double)sizeImagePx)/2.-.5;
        for (int i=0;i<sizeImagePx;i++){
            for (int ii=0;ii<sizeImagePx;ii++){
                double x=((double)(i-centerX)/radiusRing);
                double y=((double)(ii-centerY)/radiusRing);
                double theta=Math.atan2(y,x);
                double rho=Math.sqrt(x*x+y*y);
                if (rho<1){
                    image[i][ii]=theta/Math.PI;
                    //IJ.log("im "+i+"  "+ii+"  "+image[i][ii]);
                }
            }
        }
        
        //ImageShow.imshow(image,"zernike "+n+"  "+m);
        
        
        return image;
        
    }
    
    
    
    
        public double [][] getZReversed(int n, int m){
        
        double [][] image = new double [sizeImagePx][sizeImagePx];
        
        double centerX=((double)sizeImagePx)/2.-.5;
        double centerY=((double)sizeImagePx)/2.-.5;
        for (int i=0;i<sizeImagePx;i++){
            for (int ii=0;ii<sizeImagePx;ii++){
                double x=((double)(i-centerX)/radiusRing);
                double y=((double)(ii-centerY)/radiusRing);
                double theta=Math.atan2(y,x);
                double rho=Math.sqrt(x*x+y*y);
                if (rho<1){
                    image[i][ii]=this.getZReversed(n, m, theta, rho);
                    //IJ.log("im "+i+"  "+ii+"  "+image[i][ii]);
                }
            }
        }
        
        //ImageShow.imshow(image,"zernike "+n+"  "+m);
        
        
        return image;
        
    }
    
        
        
    public double getZReversed(int n, int m,double theta, double rho){
        if (n<Math.abs(m)){
            IJ.log("error in zernike coefficient : n<m");
            return -1;
        }
        else{
            if (m>=0){
                double som=0;
                for (int k=0;k<=(n-m)/2;k++){
                    double truc=(Math.pow(-1, k)*fact.get(n-k));
                    double trac=fact.get(k)*fact.get(((n+m)/2)-k)*fact.get(((n-m)/2)-k);
                    som+=(truc/trac)*Math.pow(rho,n-2*k);
                }
                return som*(Math.cos(theta*m))*Math.signum(Math.sin(theta*m));
            }
            else{
                m=-m;
                double som=0;
                for (int k=0;k<=(n-m)/2;k++){
                    double truc=(Math.pow(-1, k)*fact.get(n-k));
                    double trac=fact.get(k)*fact.get(((n+m)/2)-k)*fact.get(((n-m)/2)-k);
                    som+=(truc/trac)*Math.pow(rho,n-2*k);
                }
                return som*Math.sin(theta*m)*Math.signum(Math.cos(theta*m));
            }
        }
    }
    
    
    
    
    
    
    
    //angular radial transform
    public double [][] getART(int n, int m,boolean iscomplex){
        
        double [][] image = new double [sizeImagePx][sizeImagePx];
        
        double centerX=((double)sizeImagePx)/2.-.5;
        double centerY=((double)sizeImagePx)/2.-.5;
        for (int i=0;i<sizeImagePx;i++){
            for (int ii=0;ii<sizeImagePx;ii++){
                double x=((double)(i-centerX)/radiusRing);
                double y=((double)(ii-centerY)/radiusRing);
                double theta=Math.atan2(y,x);
                double rho=Math.sqrt(x*x+y*y);
                
                if (rho<1){
                    image[i][ii]=this.getART(n, m, theta, rho,iscomplex);
                    //IJ.log("im "+i+"  "+ii+"  "+image[i][ii]);
                }
            }
        }
        
        
        
        
        
        
        //ImageShow.imshow(image,"zernike "+n+"  "+m);
        
        
        return image;
        
    }
    
    
    
    
    //art angular radial transform
    public double getART(int n, int m,double theta, double rho,boolean iscomplex){
        
        if (n==0){
            if (iscomplex){
                return (Math.sin(theta*m));
            }
            else{
                return (Math.cos(theta*m));
            }
        }
        else{
            if (iscomplex){
                return (Math.sin(theta*m))*2*Math.cos(Math.PI*n*rho);
            }
            else{
                return (Math.cos(theta*m))*2*Math.cos(Math.PI*n*rho);
                
            }

        }
        
    }
    
    
    
    
    
    
    
    class Factorial{
        int number=20;
        long [] f;
        Factorial(){
            f = new long[number];
            f[0]=1;
            for (int i=1;i<number;i++){
                f[i]=f[i-1]*i;
            }
        }
        double get(int num){
            if (num<number){
                return f[num];
            }
            else{
                double fact=f[number-1];
                for (int i = number; i <= num; i++) {
                    fact *= i;
                }
                return fact;
            }
        }
    }
}
