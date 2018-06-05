/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.data;


import ij.IJ;
import java.util.ArrayList;

/**
 *
 * @author benoit
 */


public class PLocalization{
    
    public int frame= -1;
    //public int frameEnd=-1;//because a particle can be present on multiple frames
    public int id=-1;
    public double X=-1;
    public double Y=-1;
    public double Z=-1;
    
    public double drift_X=0;
    public double drift_Y=0;
    public double drift_Z=0;
    
    public double crlb_X=-1;
    public double crlb_Y=-1;
    public double crlb_Z=-1;
    
    public double score=-1;
    
    public double I=-1;
    public double BG=-1;
    public int occurrence=1;
    
    public ArrayList<String> otherVariableName=new ArrayList<String>();
    public double [] otherVariable=null;
    
    public boolean exists=false;
    
            
    public PLocalization(){
    }
    
    
    public PLocalization(int id,int frame,double x,double y, double z, double I,double BG,double score,double crlbX,double crlbY,double crlbZ){
        this.frame= frame;
        this.X=Math.floor(100* (x))/100;
        this.Y=Math.floor(100* (y))/100;
        this.Z=Math.floor(100* (z))/100;
        this.id=id;
        this.I=Math.floor(100* (I))/100;
        this.BG=Math.floor(100* (BG))/100;
        this.score=Math.floor(100* (score))/100;
        this.crlb_X=Math.floor(100* (crlbX))/100;
        this.crlb_Y=Math.floor(100* (crlbY))/100;
        this.crlb_Z=Math.floor(100* (crlbZ))/100;
        exists=true;
    }
    
    
    
    
    public PLocalization(int frame,double x,double y, double z){
        this.frame= frame;
        this.X=Math.floor(100* (x))/100;
        this.Y=Math.floor(100* (y))/100;
        this.Z=Math.floor(100* (z))/100;
        exists=true;
    }
    
    
    public PLocalization(int id,int frame,double x,double y, double z, double I,double BG,double score,double crlbX,double crlbY,double crlbZ,double driftX,double driftY,double driftZ){
        this.frame= frame;
        this.X=Math.floor(100* (x))/100;
        this.Y=Math.floor(100* (y))/100;
        this.Z=Math.floor(100* (z))/100;
        this.id=id;
        this.I=Math.floor(100* (I))/100;
        this.BG=Math.floor(100* (BG))/100;
        this.score=Math.floor(100* (score))/100;
        this.crlb_X=Math.floor(100* (crlbX))/100;
        this.crlb_Y=Math.floor(100* (crlbY))/100;
        this.crlb_Z=Math.floor(100* (crlbZ))/100;
        
        this.drift_X=Math.floor(100* (driftX))/100;
        this.drift_Y=Math.floor(100* (driftY))/100;
        this.drift_Z=Math.floor(100* (driftZ))/100;
        exists=true;
    }
    
    
    
    
    
    public PLocalization(int id,int frame,double x,double y, double z, double I,double BG,double score,double crlbX,double crlbY,double crlbZ,double driftX,double driftY,double driftZ,int occ){
        this.frame= frame;
        this.X=Math.floor(100* (x))/100;
        this.Y=Math.floor(100* (y))/100;
        this.Z=Math.floor(100* (z))/100;
        this.id=id;
        this.I=Math.floor(100* (I))/100;
        this.BG=Math.floor(100* (BG))/100;
        this.score=Math.floor(100* (score))/100;
        this.crlb_X=Math.floor(100* (crlbX))/100;
        this.crlb_Y=Math.floor(100* (crlbY))/100;
        this.crlb_Z=Math.floor(100* (crlbZ))/100;
        
        this.drift_X=Math.floor(100* (driftX))/100;
        this.drift_Y=Math.floor(100* (driftY))/100;
        this.drift_Z=Math.floor(100* (driftZ))/100;
        this.occurrence=occ;
        exists=true;
    }
    
    
    public void addListOfVariables(ArrayList<String> al){
        for (int i=0;i<al.size();i++){
            this.otherVariableName.add(al.get(i));
        }
        this.otherVariable=new double[al.size()];
    }
    
    
    
    /*public PLocalization(String s){
        
        String [] r = s.split(",");
        
        if (r.length==15){
            try{
                frame=Integer.parseInt(r[0]);
                id=Integer.parseInt(r[1]);
                X=Double.parseDouble(r[2]);
                Y=Double.parseDouble(r[3]);
                Z=Double.parseDouble(r[4]);
                I=Double.parseDouble(r[5]);
                BG=Double.parseDouble(r[6]);
                score=Double.parseDouble(r[7]);
                crlb_X=Double.parseDouble(r[8]);
                crlb_Y=Double.parseDouble(r[9]);
                crlb_Z=Double.parseDouble(r[10]);
                drift_X=Double.parseDouble(r[11]);
                drift_Y=Double.parseDouble(r[12]);
                drift_Z=Double.parseDouble(r[13]);
                occurrence=Double.parseDouble(r[14]);
                exists=true;
            }
            catch(Exception ee){exists=false;}
        }
        
        
    }*/
    
    
    
    public void clear(){
        
        
        frame= -1;
        id=-1;
        X=-1;
        Y=-1;
        Z=-1;
        I=-1;
        BG=-1;
        score=-1;
        
        crlb_X=-1;
        crlb_Y=-1;
        crlb_Z=-1;
        
        drift_X=0;
        drift_Y=0;
        drift_Z=0;
        occurrence=1;
        exists=false;
        
        
    }
    
    
    
    
    public String toString(){
        
        String s=""+id+getLabel_regex()+frame+getLabel_regex()+X+getLabel_regex()+Y+getLabel_regex()+Z+getLabel_regex()+I+getLabel_regex()+BG+getLabel_regex()+(score)+getLabel_regex()+(crlb_X)+getLabel_regex()+crlb_Y+getLabel_regex()+crlb_Z+getLabel_regex()+(drift_X)+getLabel_regex()+drift_Y+getLabel_regex()+drift_Z+getLabel_regex()+occurrence;
        for (int i=0;i<this.otherVariableName.size();i++){
            s+=getLabel_regex()+otherVariable[i];
        }
        return s;
        
        
    }
    
    
    
    
    
    public void setDrift_X(double drift_X){
        this.drift_X=drift_X;
        this.X-=drift_X;
        
    }
    
    
    public void setDrift_Y(double drift_Y){
        this.drift_Y=drift_Y;
        this.Y-=drift_Y;
    }
    
    
    public void setDrift_Z(double drift_Z){
        this.drift_Z=drift_Z;
        this.Z-=drift_Z;
        
    }
    
    
    public void setDrift(double drift_X,double drift_Y,double drift_Z){
        this.drift_X=drift_X;
        this.X-=drift_X;
        this.drift_Y=drift_Y;
        this.Y-=drift_Y;
        this.drift_Z=drift_Z;
        this.Z-=drift_Z;
    }
    
    
    public void removeDrift(){
        
        this.X+=drift_X;
        this.Y+=drift_Y;
        this.Z+=drift_Z;
        
        this.drift_X=0;
        this.drift_Y=0;
        this.drift_Z=0;
    }
    
    
    
    public String toStringName(){
        //return "\"id\",\"frame\",\"x [nm]\",\"y [nm]\",\"z [nm]\",\"intensity [photon]\",\"background\",\"score\",\"crlbX\",\"crlbY\",\"crlbZ\",\"driftX\",\"driftY\",\"driftZ\"";
        
        String s=getLabel_id()+getLabel_regex()+getLabel_frame()+getLabel_regex()+getLabel_x()+" [nm]"+getLabel_regex()+getLabel_y()+" [nm]"+getLabel_regex()+getLabel_z()+" [nm]"+getLabel_regex()+getLabel_A()+getLabel_regex()+getLabel_B()+getLabel_regex()+getLabel_score()+getLabel_regex()+getLabel_crlbX()+getLabel_regex()+getLabel_crlbY()+getLabel_regex()+getLabel_crlbZ()+getLabel_regex()+getLabel_driftX()+getLabel_regex()+getLabel_driftY()+getLabel_regex()+getLabel_driftZ()+getLabel_regex()+getLabel_occurrence();
        for (int i=0;i<this.otherVariableName.size();i++){
            s+=getLabel_regex()+otherVariableName.get(i);
        }
        return s;
        
    }
    
    
    
    public String getLabel(int number){
        switch (number) {
            case 0:  return "id";
            case 1:  return "frame";
            case 2:  return "x";
            case 3:  return "y";
            case 4:  return "z";
            case 5:  return "intensity";
            case 6:  return "background";
            case 7:  return "chi2";
            case 8:  return "crlbX";
            case 9:  return "crlbY";
            case 10: return "crlbZ";
            case 11: return "driftX";
            case 12: return "driftY";
            case 13: return "driftZ";
            case 14: return "occurrenceMerging";
        }
        for (int i=0;i<this.otherVariableName.size();i++){
            if ((number-15)==i){
                return otherVariableName.get(i);
            }
        }
        return null;
    }
    
    
    public int getNumberVariable(){
        return 15+otherVariableName.size();
    }
    
    public void setValueVariable(int variable,double value){
        switch (variable) {
            case 0:  this.id=(int)value;
            case 1:  this.frame=(int)value;
            case 2:  this.X=value;
            case 3:  this.Y=value;
            case 4:  this.Z=value;
            case 5:  this.I=value;
            case 6:  this.BG=value;
            case 7:  this.score=value;
            case 8:  this.crlb_X=value;
            case 9:  this.crlb_Y=value;
            case 10: this.crlb_Z=value;
            case 11: this.drift_X=value;
            case 12: this.drift_Y=value;
            case 13: this.drift_Z=value;
            case 14: this.occurrence=(int)value;
        }
        for (int i=0;i<this.otherVariableName.size();i++){
            if ((variable-15)==i){
                otherVariable[i]=value;
            }
        }
    }
    
    public double getValueVariable(int variable){
        switch (variable) {
            case 0:  return this.id;
            case 1:  return this.frame;
            case 2:  return this.X;
            case 3:  return this.Y;
            case 4:  return this.Z;
            case 5:  return this.I;
            case 6:  return this.BG;
            case 7:  return this.score;
            case 8:  return this.crlb_X;
            case 9:  return this.crlb_Y;
            case 10: return this.crlb_Z;
            case 11: return this.drift_X;
            case 12: return this.drift_Y;
            case 13: return this.drift_Z;
            case 14: return this.occurrence;
        }
        for (int i=0;i<this.otherVariableName.size();i++){
            if ((variable-15)==i){
                return otherVariable[i];
            }
        }
        IJ.log("variable not found in PLocalization class");
        return -1;
    }
    
    public int getNumberOtherVariable(){
        return this.otherVariableName.size();
    }
    
    public void setValueOtherVariable(int numberOtherVariable,double value){
        otherVariable[numberOtherVariable]=value;
    }
    
    public double getValueOtherVariable(int numberOtherVariable){
        if (numberOtherVariable<otherVariableName.size())
            return otherVariable[numberOtherVariable];
        else{
            IJ.log("variable not found in PLocalization class");
            return -1;
        }
            
    }
    
    
    public String getLabel_id(){
        return "id";
    }
    
    public String getLabel_frame(){
        return "frame";
    }
    
    public String getLabel_x(){
        return "x";
    }
    
    public String getLabel_y(){
        return "y";
    }
    
    public String getLabel_z(){
        return "z";
    }
    
    public String getLabel_A(){
        return "intensity";
    }
    
    public String getLabel_B(){
        return "background";
    }
    
    public String getLabel_score(){
        return "chi2";
    }
    
    public String getLabel_crlbX(){
        return "crlbX";
    }
    
    public String getLabel_crlbY(){
        return "crlbY";
    }
    
    public String getLabel_crlbZ(){
        return "crlbZ";
    }
    
    public String getLabel_driftX(){
        return "driftX";
    }
    
    public String getLabel_driftY(){
        return "driftY";
    }
    
    public String getLabel_driftZ(){
        return "driftZ";
    }
    
    
    public String getLabel_Xcam2(){
        return "X_cam2";
    }
    
    public String getLabel_Ycam2(){
        return "Y_cam2";
    }
    
    public String getLabel_Zcam2(){
        return "Z_cam2";
    }
    
    public String getLabel_occurrence(){
        return "occurrenceMerging";
    }
    
    public String getLabel_regex(){
        return ",";
    }
    
    
    
    
    public PLocalization copy(){
        PLocalization p = new PLocalization(id,frame,X,Y, Z, I,BG,score,crlb_X,crlb_Y,crlb_Z,drift_X,drift_Y,drift_Z,occurrence);
        if (otherVariableName.size()>0){
            p.addListOfVariables(otherVariableName);
            for (int i=0;i<otherVariableName.size();i++){
                p.setValueOtherVariable(i, otherVariable[i]);
            }
        }
        return p;
    }
    
    
    
    
}
