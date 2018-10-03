/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;


import ij.IJ;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 *
 * @author benoit
 */

public class PolynomialFit {
    
    public double [][] a;//parameters
    
    public int ordre;
    public int dimensionRef;
    public int dimensionFit;
    int [][] indices;//[size of a][dimension] : index ef exposant
    int size;
    int base;
    int nbdata;
    double[][] X;
    double[][] Y;
    //dataRef : premiere dimension : la dimension X, Y, Z, ...
    //dataRef : deuxieme dimension : les données
    
    public PolynomialFit(int ordre,double [][] dataRef,double [][] dataToFit){
        dimensionRef=dataRef.length;
        dimensionFit=dataToFit.length;
        
        //this.dimension=dataRef.length;
        this.ordre=ordre;
        
        this.Y=dataRef;
        this.X=dataToFit;
        nbdata=dataRef[0].length;
        a=new double [dimensionRef][];
        base=this.ordre+1;
        constructIndices();
        size=this.indices.length;
        /*for (int i=0;i<size;i++){
            String s="";
            for (int ii=0;ii<dimension;ii++){
                s+=indices[i][ii]+" , ";
            }
            IJ.log("u:"+i+"    "+s+"      val:"+this.getIndice(indices[i]));
        }*/
    }
    
    public PolynomialFit(int ordre,int dimensionRef,int dimensionFit){
        this.dimensionRef=dimensionRef;
        this.dimensionFit=dimensionFit;
        
        this.ordre=ordre;
        base=ordre+1;
        constructIndices();
        a=new double [dimensionRef][indices.length];
        size=indices.length;
    }
    
    public PolynomialFit(String path){
        this.load(path);
    }
    
    
    public PolynomialFit(double [][] a,int dimensionFit,int ordre){
        this.dimensionFit=dimensionFit;
        dimensionRef=a.length;
        this.ordre=ordre;
        base=ordre+1;
        constructIndices();
        size=this.indices.length;
        this.a=a;
    }
    
    
    public boolean run(){
        
        
        
        for (int d=0;d<dimensionRef;d++){
            //long dt = System.currentTimeMillis();
            double [][] XX=new double[size][size];

            double [] YY=new double[size];
            
            
            
            int nbThread=8;
            Computation [] cp = new Computation[nbThread];
            for (int t=0;t<cp.length;t++){
                cp[t]=null;
            }
            int k=0;
            for (int i=0;i<size;i++){
                boolean unfree=true;
                int positFree=-1;
                while (unfree){
                    bucl:for (int t=0;t<cp.length;t++){
                        if ((cp[t]==null)||(!cp[t].isAlive())){
                            unfree=false;
                            positFree=t;
                        }
                    }
                    try{
                        Thread.sleep(10);
                    }catch(Exception er){}
                }
                if (!unfree){
                    cp[positFree]=new Computation(i,d,XX,YY);
                    cp[positFree].start();
                }
                //IJ.log("process "+(size-i));
            }
            for (int t=0;t<cp.length;t++){
                try{
                    cp[t].join();
                }catch(Exception er){}
            }
            
            //IJ.log("time construct    " +(System.currentTimeMillis()-dt));
            //dt = System.currentTimeMillis();
            //if (XX.length>10){
            //    IJ.log("inversion matrix size "+XX.length+" ...");
            //}
            try{
                GaussianElimination ge= new GaussianElimination();
                a[d]=ge.lsolve(XX,YY);
            }catch(Exception e){IJ.log("OOPS, matrix non inversible. Try using less or more data points");a=null;return false;}
            //IJ.log("time inv    " +(System.currentTimeMillis()-dt));
        }
        return true;
        
        
    }
    
    
    
    
    public void removeParameter(int id){
        int nbMajore=(int)Math.pow(base,dimensionFit);
        if (id<this.indices.length){
            int [][] indicesNew= new int[this.indices.length-1][dimensionFit];
            for (int i=0,k=0;i<this.indices.length;i++){
                if (i==id){
                    continue;
                }
                for (int j=0;j<dimensionFit;j++){  
                    indicesNew[k][j]=indices[i][j];
                }
                k++;
            }
            indices=indicesNew;
            size=indicesNew.length;
            for (int i=0;i<a.length;i++){
                if (a[i]!=null){
                    double [] newA=new double[a[i].length-1];
                    for (int j=0,k=0;j<a[i].length;j++){
                        if (j==id){
                            continue;
                        }
                        newA[k++]=a[i][j];
                    }
                    a[i]=newA;
                }
            }
        }
        
    }
    
    
    
    
    
    
    public double [] transform(double [] X){//transform one position of N dimensions
        if (X.length!=dimensionFit){
            System.out.println("problem dimension registration transform");
        }
        double [] res = new double [dimensionRef];
        for (int d=0;d<a.length;d++){
            res[d]=0;
            for (int i=0;i<a[d].length;i++){
                double prod=a[d][i];
                //IJ.log("a "+a[d][i]);
                int [] indiceLigne = this.indices[i];
                for (int ii=0;ii<X.length;ii++){
//                    if (indiceLigne[ii]==0){
//                        
//                    }
//                    else if (indiceLigne[ii]==1){
//                        prod*=X[ii];
//                    }
//                    else{
                        prod*=Math.pow(X[ii], indiceLigne[ii]);
//                    }
                    //IJ.log("* "+X[ii]+" ^ "+indiceLigne[ii]);
                }
                res[d]+=prod;
                //IJ.log("res "+res[d]+"   "+prod);
            }
        }
        return res;
    }
    
    
    
    //indexDerivative corresponds to the variable : 0:x ; 1:y ; 2:z ...
    public double [] transformPartial1(double [] X,int indexDerivative){//transform one position of N dimensions
        if (X.length!=dimensionFit){
            System.out.println("problem dimension registration transform");
        }
        if ((indexDerivative>=X.length)||(indexDerivative<0)){
            System.out.println("problem : partial derivative impossible to compute according to index "+indexDerivative);
        }
        double [] res = new double [dimensionRef];
        for (int d=0;d<a.length;d++){
            res[d]=0;
            for (int i=0;i<a[d].length;i++){
                double prod=a[d][i];
                //IJ.log("a "+a[d][i]);
                int [] indiceLigne = this.indices[i];
                for (int ii=0;ii<X.length;ii++){
                    if (ii==indexDerivative){
                        if (indiceLigne[ii]>0){
                            prod*=indiceLigne[ii]*Math.pow(X[ii], indiceLigne[ii]-1);
                        }
                        else{
                            prod=0;
                        }
                    }
                    else{
                        prod*=Math.pow(X[ii], indiceLigne[ii]);
                    }
                    //IJ.log("* "+X[ii]+" ^ "+indiceLigne[ii]);
                }
                res[d]+=prod;
                //IJ.log("res "+res[d]+"   "+prod);
            }
        }
        return res;
    }
    
    
    
    //indexDerivative corresponds to the variable : 0:x ; 1:y ; 2:z ...
    public double [] transformPartial2(double [] X,int indexDerivative){//transform one position of N dimensions
        if (X.length!=dimensionFit){
            System.out.println("problem dimension registration transform");
        }
        if ((indexDerivative>=X.length)||(indexDerivative<0)){
            System.out.println("problem : partial derivative impossible to compute according to index "+indexDerivative);
        }
        double [] res = new double [dimensionRef];
        for (int d=0;d<a.length;d++){
            res[d]=0;
            for (int i=0;i<a[d].length;i++){
                double prod=a[d][i];
                //IJ.log("a "+a[d][i]);
                int [] indiceLigne = this.indices[i];
                for (int ii=0;ii<X.length;ii++){
                    if (ii==indexDerivative){
                        if (indiceLigne[ii]>1){
                            prod*=indiceLigne[ii]*(indiceLigne[ii]-1)*Math.pow(X[ii], indiceLigne[ii]-2);
                        }
                        else{
                            prod=0;
                        }
                    }
                    else{
                        prod*=Math.pow(X[ii], indiceLigne[ii]);
                    }
                    //IJ.log("* "+X[ii]+" ^ "+indiceLigne[ii]);
                }
                res[d]+=prod;
                //IJ.log("res "+res[d]+"   "+prod);
            }
        }
        return res;
    }
    
    
    
    
    
    
    //return position from i, j, k, ...
    int getIndice(int [] i){
        int posit=0;
        for (int t=0;t<i.length;t++){
            posit+=i[t]*((int)Math.pow(base,t));
        }
        return posit;
    }
    
    
    
    
    //inverse of getIndice : return i j k... from position
    int [] getIndices(int val){
        int [] i = new int[dimensionFit];
        for (int t=dimensionFit-1;t>=0;t--){
            i[t]=val/((int)Math.pow(base,t));
            val=val%((int)Math.pow(base,t));
        }
        return i;
    }
    
    
    
    void constructIndices(){
        int nbMajore=(int)Math.pow(base,dimensionFit);
        boolean [] ok = new boolean[nbMajore];
        int count=0;
        for (int u=0;u<nbMajore;u++){
            int [] p=this.getIndices(u);
            int valid=0;
            for (int uu=0;uu<p.length;uu++){
                valid+=p[uu];
            }
            if (valid<base){
                count++;
                ok[u]=true;
            }
            else{
                ok[u]=false;
            }
        }
        this.indices=new int[count][dimensionFit];
        int j=0;
        for (int u=0;u<nbMajore;u++){
            if (ok[u]){
                int [] p=this.getIndices(u);
                this.indices[j++]=p;
            }
        }
    }
    
    
    public void log(){
        String [] vect = new String [7];
        vect[0]="x^";
        vect[1]="y^";
        vect[2]="z^";
        vect[3]="t^";
        vect[4]="u^";
        vect[5]="v^";
        vect[6]="w^";
        
        if (dimensionFit<7){
            for (int d=0;d<dimensionRef;d++){
                String s="";
                for (int u=0;u<size;u++){
                    s+=a[d][u];
                    for (int uu=0;uu<dimensionFit;uu++){
                        if (indices[u][uu]!=0){
                            if (indices[u][uu]!=1){
                                s+="*("+vect[uu]+indices[u][uu]+")";
                            }
                            else{
                                s+="*"+vect[uu].substring(0, 1);
                            }
                        }
                    }
                    if ((u<size-1))
                        s+=" + ";
                }
                IJ.log("p("+d+")="+s);
                
            }
        }
        else{
            IJ.log("print yet not implemented with dimension > 7");
        }
        
        
    }
    
    
    
    
    public void save(String path){
        try {
            if (path==null){
                path=FileVectorLoader.chooseSavingPath("Select the file for saving registration parameters");
            }
            if (path!=null){
                PrintWriter sortie;


                    sortie = new PrintWriter(new FileWriter(path, false));
                    sortie.println(""+ordre);
                    sortie.println(""+dimensionRef);
                    sortie.println(""+dimensionFit);
                    for (int i=0;i<a.length;i++){
                        String s="";
                        for (int ii=0;ii<a[i].length&&ii<1;ii++){
                            s+=a[i][ii];
                        }
                        for (int ii=1;ii<a[i].length;ii++){
                            s+=","+a[i][ii];
                        }
                        sortie.println(""+s);
                    }
                    


                    IJ.log("Save well done !");
                    sortie.close();
            }
        
            

        } catch (Exception e) {
                e.printStackTrace();
        }
    }
    
    
    
    public String toString(){
        
        String out="";
        out+=(""+ordre)+"\n";
        out+=(""+dimensionRef)+"\n";
        out+=(""+dimensionFit)+"\n";
        for (int i=0;i<a.length;i++){
            String s="";
            for (int ii=0;ii<a[i].length&&ii<1;ii++){
                s+=a[i][ii];
            }
            for (int ii=1;ii<a[i].length;ii++){
                s+=","+a[i][ii];
            }
            out+=(""+s)+"\n";
        }
        
        return out;
    }
    
    
    
    void load(String path){
        
        if (path==null){
            path=FileVectorLoader.chooseLoadingPath("Select the file containing registration parameters");
        }
        
        if (path!=null){
            ArrayList<Double> al = new ArrayList<Double>();
            try{
                    InputStream ips=new FileInputStream(path); 
                    InputStreamReader ipsr=new InputStreamReader(ips);
                    BufferedReader br=new BufferedReader(ipsr);
                    String ligne;
                    ligne=br.readLine();
                    this.ordre=Integer.parseInt(ligne);
                    ligne=br.readLine();
                    this.dimensionRef=Integer.parseInt(ligne);
                    ligne=br.readLine();
                    this.dimensionFit=Integer.parseInt(ligne);
                    this.base=ordre+1;
                    this.constructIndices();
                    size=indices.length;
                    a=new double[dimensionFit][size];
                    int t=0;
                    try{
                        theloopa:while ((ligne=br.readLine())!=null){
                            String [] lin=ligne.split(",");
                            for (int u=0;u<lin.length;u++){
                                a[t][u]=Double.parseDouble(lin[u]);
                            }
                            t++;
                            if (t==dimensionRef){
                                break theloopa;
                            }
                        }
                        IJ.log("registration parameters loaded");
                    }catch(Exception er){IJ.log("error loading polynomial function... wrong size ? "+er);}
                    
                    br.close();
            }		
            catch (Exception e){
                IJ.log("error loading file of polynomial function...");
                    System.out.println(e.toString());
            }
        }

        
    }
    
    
    class Computation extends Thread{
        //int [] indices;
        int i;
        double [][] XX; 
        double [] YY;
        int d;
        Computation(int index,int dimension,double [][] XX, double [] YY){
            this.XX=XX;
            this.YY=YY;
            this.i=index;
            this.d=dimension;
        }
        
        public void run(){
            //compute XX
                long dtt1 = System.currentTimeMillis();
                int [] indiceLigne = indices[i];
                
                
                
                
                double [] prodLigne=new double [nbdata];//compute product for each line
                for (int xii=0;xii<nbdata;xii++){//for each data
                    prodLigne[xii]=1;
                    for (int xi=0;xi<dimensionFit;xi++){//for each dimension
                        prodLigne[xii]*=Math.pow(X[xi][xii],indiceLigne[xi]);
                    }
                }
                long dtt2 = System.currentTimeMillis();
                for (int ai=0;ai<size;ai++){//compute product for each line/column
                    int [] indiceColumn = indices[ai];
                    XX[i][ai]=0;
                    for (int xii=0;xii<nbdata;xii++){//for each data
                        double prod=prodLigne[xii];
                        for (int xi=0;xi<dimensionFit;xi++){//for each dimension
                            prod*=Math.pow(X[xi][xii],indiceColumn[xi]);
                        }
                        XX[i][ai]+=prod;
                    }
                }
                long dtt3 = System.currentTimeMillis();
                
                
                
                
                
                
                
                
                //compute YY
                YY[i]=0;
                for (int xii=0;xii<nbdata;xii++){//for each data
                    double prod=1;
                    double prod2=1;
                    for (int xi=0;xi<dimensionFit;xi++){//for each dimension
                        prod*=Math.pow(X[xi][xii],indiceLigne[xi]);
                    }
                    prod*=Y[d][xii];
                    YY[i]+=prod;
                }
                //IJ.log("YY[i] "+YY[i]);
                
                
        }
    }
    
            
    
}
