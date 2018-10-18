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
import java.text.ParseException;
import java.util.ArrayList;

/**
 *
 * @author benoit
 */
public class JSONformat {
    
    String file="{}";
    
    
    public JSONformat(){
        
    }
    
    
    public void load(String path){
        
        file=readFile(path);
        
    }
    
    public void save(String path){
        
        saveFile(path);
        
    }
    
    
    public void put(String s,int value){
        if (file.endsWith("{}")){
            file="{\""+s+"\":"+value+"}";
        }
        else{
            file=file.substring(0, file.length()-1)+",\""+s+"\":"+value+"}";
        }
    }
    
    public void put(String s,double value){
        if (file.endsWith("{}")){
            file="{\""+s+"\":"+value+"}";
        }
        else{
            file=file.substring(0, file.length()-1)+",\""+s+"\":"+value+"}";
        }
    }
    
    public void put(String s,boolean value){
        if (file.endsWith("{}")){
            file="{\""+s+"\":"+value+"}";
        }
        else{
            file=file.substring(0, file.length()-1)+",\""+s+"\":"+value+"}";
        }
    }
    
    public void put(String s,int [] value){
        if (file.endsWith("{}")){
            file="{\""+s+"\":[";
            for (int i=0;i<value.length;i++){
                if (i==0){
                    file+=value[i];
                }
                else{
                    file+=","+value[i];
                }
            }
            file+="]}";
        }
        else{
            file=file.substring(0, file.length()-1)+",\""+s+"\":[";//+value+"}";
            for (int i=0;i<value.length;i++){
                if (i==0){
                    file+=value[i];
                }
                else{
                    file+=","+value[i];
                }
            }
            file+="]}";
        }
    }
    
    public void put(String s,double [] value){
        if (file.endsWith("{}")){
            file="{\""+s+"\":[";
            for (int i=0;i<value.length;i++){
                if (i==0){
                    file+=value[i];
                }
                else{
                    file+=","+value[i];
                }
            }
            file+="]}";
        }
        else{
            file=file.substring(0, file.length()-1)+",\""+s+"\":[";//+value+"}";
            for (int i=0;i<value.length;i++){
                if (i==0){
                    file+=value[i];
                }
                else{
                    file+=","+value[i];
                }
            }
            file+="]}";
        }
    }
    
    
    public void put(String s,boolean [] value){
        if (file.endsWith("{}")){
            file="{\""+s+"\":[";
            for (int i=0;i<value.length;i++){
                if (i==0){
                    file+=value[i];
                }
                else{
                    file+=","+value[i];
                }
            }
            file+="]}";
        }
        else{
            file=file.substring(0, file.length()-1)+",\""+s+"\":[";//+value+"}";
            for (int i=0;i<value.length;i++){
                if (i==0){
                    file+=value[i];
                }
                else{
                    file+=","+value[i];
                }
            }
            file+="]}";
        }
    }
    
    public String get(String s){
        int id=file.indexOf(s);
        if (id>=0){
            int start=1+file.indexOf(":", id);
            int end=file.indexOf(",", start);
            if (end<0){
                end=file.indexOf("}");
            }
            //IJ.log(""+s+"  "+id+"  "+start+"  "+end+"   "+file.substring(start));
            return file.substring(start, end);
        }
        else{
            return null;
        }
    }
    
    
    public String [] getVect(String s){
        int id=file.indexOf(s);
        int start=1+file.indexOf("[", id);
        int end=file.indexOf("]", id);
        String vect=file.substring(start, end);
        return vect.split(",");
    }
    
    
    
    
    String readFile(String path){
        
        
        String s="";
        try{
            InputStream ips=new FileInputStream(path); 
            
            InputStreamReader ipsr=new InputStreamReader(ips);
            
            BufferedReader br=new BufferedReader(ipsr);
            String ligne;

            while ((ligne=br.readLine())!=null){
                s+=ligne;

            }
            br.close();
        }		
        catch (Exception e){
                System.out.println(e.toString());
                return null;
        }
        return s;
    }
    
    
    void saveFile(String path){
        try{
            PrintWriter sortie;
            
            sortie = new PrintWriter(new FileWriter(path, false));
            
            sortie.println(file);
            sortie.close();
            
        }		
        catch (Exception e){
            e.printStackTrace();
            System.out.println(e.toString());
        }
    }
    
}
