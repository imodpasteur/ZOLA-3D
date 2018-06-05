/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.utils;


import ij.IJ;
import ij.io.OpenDialog;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import javax.swing.JFileChooser;

/**
 *
 * @author benoit
 */
public class FileVectorLoader {
    
    public static double [] getColumnFromFile(String path,int numColumn,String regex){
        
        if (path==null){
                path=chooseSavingPath("Select the file containing the table to load");
            }
        
        
        ArrayList<Double> al = new ArrayList<Double>();
        int nbLine=0;
        try{
                InputStream ips=new FileInputStream(path); 
                InputStreamReader ipsr=new InputStreamReader(ips);
                BufferedReader br=new BufferedReader(ipsr);
                String ligne;
                
                while ((ligne=br.readLine())!=null){
                    if (ligne.length()>2){
                        String [] lin=ligne.split(regex);
                        try{
                            if (lin.length>numColumn)
                                al.add(new Double(lin[numColumn]));
                        }catch(Exception syubcyfz){
                            //IJ.log("impossible to convert ... in double");
                        }
                    }

                }
                nbLine++;
                br.close();
        }		
        catch (Exception e){
                System.out.println(e.toString());
        }
        
        //IJ.log(""+al.size()+" element readable for "+nbLine+" lines");
        double [] r = new double[al.size()];
        for (int i=0;i<al.size();i++){
            r[i]=(double)al.get(i);
        }
        return r;
    }
    
    
    
    public static double [][] getTableFromFile(String path,String regex){
        
        if (path==null){
                path=chooseSavingPath("Select the file containing the table to load");
            }
        
        
        ArrayList<double []> al = new ArrayList<double []>();
        int nbLine=0;
        try{
                InputStream ips=new FileInputStream(path); 
                InputStreamReader ipsr=new InputStreamReader(ips);
                BufferedReader br=new BufferedReader(ipsr);
                String ligne;
                
                while ((ligne=br.readLine())!=null){
                    if (ligne.length()>2){
                        String [] lin=ligne.split(regex);
                        double [] v= new double [lin.length];
                        try{
                            for (int i=0;i<v.length;i++){
                                v[i]=Double.parseDouble(lin[i]);
                            }
                            al.add(v);
                        }catch(Exception syubcyfz){
                            //IJ.log("impossible to convert ... in double");
                        }
                    }

                }
                nbLine++;
                br.close();
        }		
        catch (Exception e){
                System.out.println(e.toString());
        }
        
        //IJ.log(""+al.size()+" element readable for "+nbLine+" lines");
        double [][] r = new double[al.size()][];
        for (int i=0;i<al.size();i++){
            r[i]=new double [al.get(i).length];
            for (int ii=0;ii<al.get(i).length;ii++){
                r[i][ii]=al.get(i)[ii];
            }
        }
        return r;
    }
    
    
    
    
    public static void saveTableInFile(String path,double [][] table,String regex){
        saveTableInFile(path,table,regex,false);
    }
    public static void saveTableInFile(String path,double [][] table,String regex,boolean transpose){

        try {
            if (path==null){
                path=chooseSavingPath("Select the file where you want to save the table");
            }

            PrintWriter sortie;
            

                sortie = new PrintWriter(new FileWriter(path, false));
                
                if (!transpose){
                    for (int i=0;i<table.length;i++){

                        if (i%1000==0){
                            //IJ.log("left : "+(i/1000));
                        }

                        String s="";
                        for (int ii=0;ii<table[i].length&&ii<1;ii++){
                            s+=table[i][ii];
                        }
                        for (int ii=1;ii<table[i].length;ii++){
                            s+=","+table[i][ii];
                        }
                        sortie.println(""+s);
                    }
                }
                else{
                    if (table[0]!=null){
                        for (int i=0;i<table[0].length;i++){

                            if (i%1000==0){
                                //IJ.log("left : "+(table[0].length/1000-i/1000));
                            }

                            String s="";
                            for (int ii=0;ii<table.length&&ii<1;ii++){
                                s+=table[ii][i];
                            }
                            for (int ii=1;ii<table.length;ii++){
                                s+=","+table[ii][i];
                            }
                            sortie.println(""+s);
                        }
                    }
                }
                
                
                
                
                IJ.log("Save well done !");
                sortie.close();
            
        
            

        } catch (Exception e) {
                e.printStackTrace();
        }
    }
    
    
    
    
    
    public static String chooseSavingPath(String message){
        IJ.log("Saving: "+message);
        try {
            JFileChooser dialogue ;

            dialogue = new JFileChooser(new File(OpenDialog.getDefaultDirectory()));

            File fichier;
            int x=dialogue.showSaveDialog(null);
            if (x==JFileChooser.APPROVE_OPTION) {
                fichier = dialogue.getSelectedFile();
                OpenDialog.setDefaultDirectory(fichier.getParent());
                return fichier.getPath();
                
            }
            else if (x==JFileChooser.CANCEL_OPTION){
                IJ.log("Save canceled !");
                return null;
            }
        
        
            

        } catch (Exception e) {
                e.printStackTrace();
        }
        return null;
                
    }
    
    public static String chooseLoadingPath(String message){
        IJ.log("Loading: "+message);
        try {
            JFileChooser dialogue ;

            dialogue = new JFileChooser(new File(OpenDialog.getDefaultDirectory()));

            File fichier;
            int x=dialogue.showOpenDialog(null);
            if (x==JFileChooser.APPROVE_OPTION) {
                fichier = dialogue.getSelectedFile();
                
                OpenDialog.setDefaultDirectory(fichier.getParent());
                return fichier.getPath();
                
            }
            else if (x==JFileChooser.CANCEL_OPTION){
                IJ.log("Load canceled !");
                return null;
            }
        
        
            

        } catch (Exception e) {
                e.printStackTrace();
        }
        return null;
                
    }
    
    
}
