/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;


import ij.IJ;

public class GaussianElimination {
    private static final double EPSILON = 1e-20;
    private static final float EPSILONF = 0;
    
    
    
    //çif non square matrix... multiply A and b by transpose of A --> At.A.X=At.B
    
    // Gaussian elimination with partial pivoting
    public static double[] lsolve(double[][] A, double[] b) {
        int N  = b.length;
        
        for (int p = 0; p < N; p++) {
            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (Math.abs(A[p][p]) <= EPSILON) {
                throw new RuntimeException("Matrix is singular or nearly singular ");
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }
        // back substitution
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
    
    
    public static float[][] lsolve(float[][] A, float[][] b) {
        int N  = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            float[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            float   t    = b[p][0]; b[p][0] = b[max][0]; b[max][0] = t;

            // singular or nearly singular
            if (Math.abs(A[p][p]) <= EPSILONF) {
                throw new RuntimeException("Matrix is singular or nearly singular "+A[p][p]);
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                float alpha = A[i][p] / A[p][p];
                b[i][0] -= alpha * b[p][0];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        float[][] x = new float[N][1];
        for (int i = N - 1; i >= 0; i--) {
            float sum = (float)0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j][0];
            }
            x[i][0] = (b[i][0] - sum) / A[i][i];
        }
        return x;
    }
    
    
    
    
    
    public static double[][] lsolve(double[][] A, double[][] b) {
        int N  = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p][0]; b[p][0] = b[max][0]; b[max][0] = t;

            // singular or nearly singular
            if (Math.abs(A[p][p]) <= EPSILONF) {
                throw new RuntimeException("Matrix is singular or nearly singular "+A[p][p]);
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i][0] -= alpha * b[p][0];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[][] x = new double[N][1];
        for (int i = N - 1; i >= 0; i--) {
            double sum = (double)0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j][0];
            }
            x[i][0] = (b[i][0] - sum) / A[i][i];
        }
        return x;
    }
    
    
    
    
}
