/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;


import ij.IJ;
/**
 *
 * code with Licence CPOL (free)
 */
public class Matrixe {
 
    private int nrows;
    private int ncols;
    private double[][] data;
 
    public Matrixe(double[][] dat) {
        this.data = dat;
        this.nrows = dat.length;
        this.ncols = dat[0].length;
    }
    
    public Matrixe(double[] vect) {
        double [][] dat = new double [vect.length][1];
        for (int i=0;i<vect.length;i++){
            dat[i][0]=vect[i];
        }
        this.data = dat;
        this.nrows = dat.length;
        this.ncols = dat[0].length;
    }
 
    public Matrixe(int nrow, int ncol) {
        this.nrows = nrow;
        this.ncols = ncol;
        data = new double[nrow][ncol];
    }

    public double [][] getMatrixe(){
        return data;
    }
    
    
    public void getMatrixe(double [][] mat){
        for (int i=0;i<getNrows();i++) {
            for (int j=0;j<getNcols();j++) {
                mat[i][j]=data[i][j];
            }
        }
    }
    
    
    
    public void free(){
        this.data=null;
    }
    
    public int getNrows(){
        return nrows;
    }
    public int getNcols(){
        return ncols;
    }
    
    
    public double getValueAt(int i,int j){
        return data[i][j];
    }
    
    public void setValueAt(int i, int j,double v){
        data[i][j]=v;
    }
    
    private int size(){
        return getNrows();
    }
    
    private boolean isSquare(){
        if (getNrows()==getNcols()){
            return true;
        }
        else{
            return false;
        }
    }
    
    
    double [][] getData(){
        return data;
    }
    
    // return C = A * B
    public Matrixe times(Matrixe B) {
        Matrixe A = this;
        if (A.getNcols() != B.getNrows()) throw new RuntimeException("Illegal matrix dimensions.");
        Matrixe C = new Matrixe(A.getNrows(), B.getNcols());
        for (int i = 0; i < C.getNrows(); i++)
            for (int j = 0; j < C.getNcols(); j++)
                for (int k = 0; k < A.getNcols(); k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }
    
    
    public static Matrixe times(Matrixe A,Matrixe B) {
        if (A.getNcols() != B.getNrows()) throw new RuntimeException("Illegal matrix dimensions.");
        Matrixe C = new Matrixe(A.getNrows(), B.getNcols());
        for (int i = 0; i < C.getNrows(); i++)
            for (int j = 0; j < C.getNcols(); j++)
                for (int k = 0; k < A.getNcols(); k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }
    
    
    public static Matrixe transpose(Matrixe matrix) {
        Matrixe transposedMatrix = new Matrixe(matrix.getNcols(), matrix.getNrows());
        for (int i=0;i<matrix.getNrows();i++) {
            for (int j=0;j<matrix.getNcols();j++) {
                transposedMatrix.setValueAt(j, i, matrix.getValueAt(i, j));
            }
        }
        return transposedMatrix;
    }

    
    
    public static void transpose(Matrixe matrix,Matrixe transposedMatrix) {
        for (int i=0;i<matrix.getNrows();i++) {
            for (int j=0;j<matrix.getNcols();j++) {
                transposedMatrix.setValueAt(j, i, matrix.getValueAt(i, j));
            }
        }
    }


    
    
    
    public static double determinant(Matrixe matrix) throws Exception {
    if (!matrix.isSquare())
        throw new Exception("matrix need to be square.");
    if (matrix.size()==1) {
        return (matrix.getValueAt(0, 0));
    }
    if (matrix.size()==2) {
        return (matrix.getValueAt(0, 0) * matrix.getValueAt(1, 1)) - ( matrix.getValueAt(0, 1) * matrix.getValueAt(1, 0));
    }
    double sum = 0.0;
    for (int i=0; i<matrix.getNcols(); i++) {
        sum += Math.pow(-1,i+1+1) * matrix.getValueAt(0, i) * determinant(createSubMatrix(matrix, 0, i));
    }
    return sum;
}


    
    
    
    public static Matrixe createSubMatrix(Matrixe matrix, int excluding_row, int excluding_col) {
    Matrixe mat = new Matrixe(matrix.getNrows()-1, matrix.getNcols()-1);
    int r = -1;
    for (int i=0;i<matrix.getNrows();i++) {
        if (i==excluding_row)
            continue;
            r++;
            int c = -1;
        for (int j=0;j<matrix.getNcols();j++) {
            if (j==excluding_col)
                continue;
            mat.setValueAt(r, ++c, matrix.getValueAt(i, j));
        }
    }
    return mat;
}


    
    
    
    
    public static Matrixe cofactor(Matrixe matrix) throws Exception {
        Matrixe mat = new Matrixe(matrix.getNrows(), matrix.getNcols());
        for (int i=0;i<matrix.getNrows();i++) {
            for (int j=0; j<matrix.getNcols();j++) {
                mat.setValueAt(i, j, Math.pow(-1,i+j+2) * determinant(createSubMatrix(matrix, i, j)));
            }
        }

        return mat;
    }
    
    
    public static void cofactor(Matrixe matrix,Matrixe matrixResult) throws Exception {
        
        for (int i=0;i<matrix.getNrows();i++) {
            for (int j=0; j<matrix.getNcols();j++) {
                matrixResult.setValueAt(i, j, Math.pow(-1,i+j+2) * determinant(createSubMatrix(matrix, i, j)));
            }
        }
    }

    private Matrixe multiplyByConstant(double v) throws Exception {
        Matrixe mat = new Matrixe(getNrows(), getNcols());
        for (int i=0;i<getNrows();i++) {
            for (int j=0; j<getNcols();j++) {
                mat.setValueAt(i, j, this.data[i][j]*v);
            }
        }

        return mat;
    }
    
    
    public static Matrixe inverse(Matrixe matrix) throws Exception {
        double det=1.0/determinant(matrix);
        Matrixe mmm = Matrixe.cofactor(matrix);
        Matrixe mat = matrix.transpose(mmm);
        for (int i=0;i<mat.getNrows();i++) {
            for (int j=0; j<mat.getNcols();j++) {

                mat.setValueAt(i, j, mat.getValueAt(i, j)*det);
            }
        }

        return mat;
    }
    
    
    public static void inverse(Matrixe matrix,Matrixe matrixresult) throws Exception {

        double det=1.0/determinant(matrix);
        Matrixe.cofactor(matrix,matrixresult);
        Matrixe mtmp=Matrixe.transpose(matrixresult);

        for (int i=0;i<matrix.getNrows();i++) {
            for (int j=0; j<matrix.getNcols();j++) {

                matrixresult.setValueAt(i, j, mtmp.getValueAt(i, j)*det);
            }
        }
        mtmp.free();
    }
    
    
    
    
    

    
}





