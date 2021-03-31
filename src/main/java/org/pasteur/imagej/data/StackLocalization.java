/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.data;

import ij.IJ;
import java.util.ArrayList;

import ij.io.OpenDialog;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.concurrent.locks.ReentrantLock;
import javax.swing.JFileChooser;


/*import org.pasteur.smlm.io.SMLMio;*/

/**
 *
 * @author benoit
 */
public class StackLocalization {
    
    boolean merged=false;
    boolean driftCorrected=false;
    
    static Save save;
    
    public ArrayList<FrameLocalization> fl = new ArrayList<FrameLocalization>();
    
    public StackLocalization(){
        
    }
    
    public StackLocalization(StackLocalization sl){
        for (int i=0;i<sl.fl.size();i++){
            FrameLocalization fl = new FrameLocalization(sl.fl.get(i).numFrame);
            for (int ii=0;ii<sl.fl.get(i).loc.size();ii++){
                fl.loc.add(sl.fl.get(i).loc.get(ii).copy());
            }
            this.fl.add(fl);
        }
    }
    
    
    
    //load stack localization
    public StackLocalization(String path){
        
        if (path.toLowerCase().endsWith(".csv")){
            this.loadCSV(path);
        }
        else if (path.toLowerCase().endsWith(".smlm")){
            this.loadSMLM(path);
        }
        else{
            IJ.log("file extension not recognized");
        }
        
    }
    
    
    
    
    private void loadSMLM(String path){
        
        /*
        
        SMLMio smlmio = new SMLMio();
        
        SMLMio.Data d=smlmio.import_smlm(path);
        
        String [] key = d.getTableKeys();
        
        if (key.length>1){
            IJ.log("WARNING: ZOLA can import only one file");
        }
        else if (key.length==0){
            IJ.log("no table found");
        }
        else{
            SMLMio.Table t=d.getTable(key[0]);
            ArrayList<PLocalization> al = new ArrayList<PLocalization>();
            String [] head=t.getHeaders();
            
            int id_id=-1;
            int id_frame=-1;
            int id_x=-1;
            int id_y=-1;
            int id_z=-1;
            int id_a=-1;
            int id_b=-1;
            int id_score=-1;
            int id_crlbx=-1;
            int id_crlby=-1;
            int id_crlbz=-1;
            int id_driftx=-1;
            int id_drifty=-1;
            int id_driftz=-1;
            int id_occ=-1;

            ArrayList<String> otherVariable=new ArrayList<String>();
            ArrayList<Integer> otherVariableID=new ArrayList<Integer>();

            int nbParam=head.length;
            for (int i=0;i<nbParam;i++){

                int xpt;
                if (((xpt=head[i].indexOf(PLocalization.getLabel_id()))!=-1)&&(xpt<5)){
                    id_id=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_frame()))!=-1)&&(xpt<5)){
                    id_frame=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_x()))!=-1)&&(xpt<5)){
                    id_x=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_y()))!=-1)&&(xpt<5)){
                    id_y=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_z()))!=-1)&&(xpt<5)){
                    id_z=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_A()))!=-1)&&(xpt<5)){
                    id_a=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_B()))!=-1)&&(xpt<5)){
                    id_b=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_score()))!=-1)&&(xpt<5)){
                    id_score=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_crlbX()))!=-1)&&(xpt<5)){
                    id_crlbx=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_crlbY()))!=-1)&&(xpt<5)){
                    id_crlby=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_crlbZ()))!=-1)&&(xpt<5)){
                    id_crlbz=i;

                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_driftX()))!=-1)&&(xpt<5)){
                    id_driftx=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_driftY()))!=-1)&&(xpt<5)){
                    id_drifty=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_driftZ()))!=-1)&&(xpt<5)){
                    id_driftz=i;
                }
                else if (((xpt=head[i].indexOf(PLocalization.getLabel_occurrence()))!=-1)&&(xpt<5)){
                    id_occ=i;
                }
                else if (((xpt=head[i].indexOf("score"))!=-1)&&(xpt<5)){
                    id_score=i;
                }
                else{

                    otherVariable.add(head[i]);
                    otherVariableID.add(i);
                }
            }

            if (id_occ==-1){
                this.merged=false;
            }
            else{
                this.merged=true;
            }
            if ((id_driftx==-1)&&(id_drifty==-1)&&(id_driftz==-1)){
                this.driftCorrected=false;
            }
            else{
                this.driftCorrected=true;
            }
            
            
            
            
            double [][] data = t.getTable();
            
            
            
            
            Hashtable<Integer, Integer> frame = new Hashtable<Integer, Integer>();
            ArrayList<Integer> frameId = new ArrayList<Integer>();
            Hashtable<Integer, Integer> frameIdReverse = new Hashtable<Integer, Integer>();
            if ((id_frame>=0)){
                
                
                for (int i=0;i<data.length;i++){
                    if (frame.get(data[i][id_frame])==null){
                        frame.put((int)data[i][id_frame],1);
                        frameId.add((int)data[i][id_frame]);
                    }
                    else{
                        frame.put((int)data[i][id_frame], frame.get(data[i][id_frame])+1);
                    }
                }
            }
            else{
                IJ.log("frame data not recognized");
                this.fl.add(new FrameLocalization(0));
                frameIdReverse.put(0, 0);
            }
            
            int length=data.length;
            
            
            
            for (int f=0;f<frame.size();f++){
                this.fl.add(new FrameLocalization(frameId.get(f)));
                frameIdReverse.put(frameId.get(f), f);
            }
            
            
            
            for (int i=0;i<length;i++){
                
                try{
                    PLocalization p= new PLocalization();
                    p.addListOfVariables(otherVariable);
                    
                    
                    
                    
                    p.exists=true;
                    if (id_id>=0){
                        //IJ.log("id ok");
                        p.id=(int)data[i][id_id];
                        
                    }
                    else{
                        p.id=i;
                    }
                    
                    
                    if (id_frame>=0){
                        //IJ.log("fr ok "+split[id_frame]);
                        p.frame=(int)data[i][id_frame];
                        
                    }
                    else{
                        p.frame=0;
                    }
                    
                    
                    if (id_x>=0){
                        //IJ.log("x ok");
                        p.X=(double)data[i][id_x];
                        
                    }
                    else{
                        p.X=0;
                    }
                    if (id_y>=0){
                        //IJ.log("y ok");
                        p.Y=(double)data[i][id_y];
                        
                    }
                    else{
                        p.Y=0;
                    }
                    if (id_z>=0){
                        //IJ.log("z ok");
                        p.Z=(double)data[i][id_z];
                        
                    }
                    else{
                        p.Z=0;
                    }
                    if (id_a>=0){
                        //IJ.log("a ok");
                        p.I=(double)data[i][id_a];
                        
                    }
                    else{
                        p.I=-1;
                    }
                    if (id_b>=0){
                        //IJ.log("b ok");
                        p.BG=(double)data[i][id_b];
                        
                    }
                    else{
                        p.BG=-1;
                    }
                    if (id_score>=0){
                        //IJ.log("score ok");
                        p.score=(double)data[i][id_score];
                        
                    }
                    else{
                        p.score=-1;
                    }
                    if (id_crlbx>=0){
                        p.crlb_X=(double)data[i][id_crlbx];
                        
                    }
                    else{
                        p.crlb_X=-1;
                    }
                    if (id_crlby>=0){
                        p.crlb_Y=(double)data[i][id_crlby];
                        
                    }
                    else{
                        p.crlb_Y=-1;
                    }
                    if (id_crlbz>=0){
                        p.crlb_Z=(double)data[i][id_crlbz];
                        
                    }
                    else{
                        p.crlb_Z=-1;
                    }
                    if (id_driftx>=0){
                        p.drift_X=(double)data[i][id_driftx];
                        
                    }
                    else{
                        p.drift_X=0;
                    }
                    if (id_drifty>=0){
                        p.drift_Y=(double)data[i][id_drifty];
                        
                    }
                    else{
                        p.drift_Y=0;
                    }
                    if (id_driftz>=0){
                        p.drift_Z=(double)data[i][id_driftz];
                        
                    }
                    else{
                        p.drift_Z=0;
                    }
                    if (id_occ>=0){
                        p.occurrence=(int)data[i][id_occ];
                        
                    }
                    else{
                        p.occurrence=1;
                    }
                    
                    
                    long time2=System.currentTimeMillis();
                    for (int ii=0;ii<otherVariable.size();ii++){
                        p.setValueOtherVariable(ii,(double)otherVariableID.get(ii));
                        
                    }
                    
                    
                    this.fl.get(frameIdReverse.get(p.frame)).loc.add(p);
                }catch(Exception e){}
                //IJ.log("value found "+p.id+"  "+p.X+"  "+p.Y+"   "+p.frame);
                //IJ.log("value found "+p.id+"  "+p.frame+"  "+p.drift_X+"  "+p.drift_Y);
                
                
                
                
                
            }
            
            
            IJ.log("Frame number: "+frame.size());
            IJ.log("Localization number: "+length);
            
            t=null;
            frame=null;
        }
        */
        
        
    }
    
    private void loadCSV(String path){
        
        ArrayList<PLocalization> al = new ArrayList<PLocalization>();
        
        try{
                InputStream ips=new FileInputStream(path); 
                InputStreamReader ipsr=new InputStreamReader(ips);
                BufferedReader br=new BufferedReader(ipsr);
                String ligne;
                ArrayList<String> otherVariable=new ArrayList<String>();
                ArrayList<Integer> otherVariableID=new ArrayList<Integer>();
                
                int id_id=-1;
                int id_frame=-1;
                int id_x=-1;
                int id_y=-1;
                int id_z=-1;
                int id_a=-1;
                int id_b=-1;
                int id_score=-1;
                int id_crlbx=-1;
                int id_crlby=-1;
                int id_crlbz=-1;
                int id_driftx=-1;
                int id_drifty=-1;
                int id_driftz=-1;
                int id_occ=-1;
                //read first line
                int nbParam=0;
                String regex=",";
                if ((ligne=br.readLine())!=null){
                    PLocalization p= new PLocalization();
                    regex=p.getLabel_regex();
                    String [] split=ligne.split(regex);
                    nbParam=split.length;
                    for (int i=0;i<nbParam;i++){
                    //for (int i=nbParam-1;i>=0;i--){
                        //IJ.log("split-"+split[i]);
                        int xpt;
                        if (((xpt=split[i].indexOf(p.getLabel_id()))!=-1)&&(xpt<5)&&(id_id==-1)){
                            id_id=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_frame()))!=-1)&&(xpt<5)&&(id_frame==-1)){
                            id_frame=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_x()))!=-1)&&(xpt<5)&&(id_x==-1)){
                            id_x=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_y()))!=-1)&&(xpt<5)&&(id_y==-1)){
                            id_y=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_z()))!=-1)&&(xpt<5)&&(id_z==-1)){
                            id_z=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_A()))!=-1)&&(xpt<5)&&(id_a==-1)){
                            id_a=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_B()))!=-1)&&(xpt<5)&&(id_b==-1)){
                            id_b=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_score()))!=-1)&&(xpt<5)&&(id_score==-1)){
                            id_score=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbX()))!=-1)&&(xpt<5)&&(id_crlbx==-1)){
                            id_crlbx=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbY()))!=-1)&&(xpt<5)&&(id_crlby==-1)){
                            id_crlby=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbZ()))!=-1)&&(xpt<5)&&(id_crlbz==-1)){
                            id_crlbz=i;
                            
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftX()))!=-1)&&(xpt<5)&&(id_driftx==-1)){
                            id_driftx=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftY()))!=-1)&&(xpt<5)&&(id_drifty==-1)){
                            id_drifty=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftZ()))!=-1)&&(xpt<5)&&(id_driftz==-1)){
                            id_driftz=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_occurrence()))!=-1)&&(xpt<5)&&(id_occ==-1)){
                            id_occ=i;
                        }
                        else if (((xpt=split[i].indexOf("score"))!=-1)&&(xpt<5)&&(id_score==-1)){
                            id_score=i;
                        }
                        else{
                            
                            otherVariable.add(split[i]);
                            otherVariableID.add(i);
                            
                        }
                        
                    }
                }
                
                
                
                if (id_occ==-1){
                    this.merged=false;
                }
                else{
                    this.merged=true;
                }
                if ((id_driftx==-1)&&(id_drifty==-1)&&(id_driftz==-1)){
                    this.driftCorrected=false;
                }
                else{
                    this.driftCorrected=true;
                }
                
                //IJ.log("id frame "+id_frame+"  "+id_x+"  "+id_y+"  "+id_z);
                
                
                int nbmol=0;
                
                int maxFrame=0;
                int num=0;
                boolean firstLineValue=true;
                toto:while ((ligne=br.readLine())!=null){
                    
                    if (ligne.length()>3){
                        
                        if (firstLineValue){
                            String [] spliter=ligne.split(regex);
                            
                            if (spliter.length!=nbParam){
                                regex=" ";
                                spliter=ligne.split(regex);
                                
                                if (spliter.length!=nbParam){
                                    regex=";";
                                    spliter=ligne.split(regex);
                                    
                                    if (spliter.length!=nbParam){
                                        IJ.log("regex not recognized: it should be coma");
                                        regex=",";
                                    }
                                }
                            }
                            firstLineValue=false;
                        }
                        String [] split=ligne.split(regex);
                        //PLocalization p= new PLocalization(ligne);
                        PLocalization p= new PLocalization();
                        p.addListOfVariables(otherVariable);
                        p.exists=true;
                        try{
                            if (id_id>=0){
                                //IJ.log("id ok");
                                p.id=(int)Double.parseDouble(split[id_id]);
                                
                            }
                            else{
                                p.id=num++;
                            }
                            if (id_frame>=0){
                                //IJ.log("fr ok "+split[id_frame]);
                                p.frame=(int)Double.parseDouble(split[id_frame]);
                            }
                            else{
                                p.frame=0;
                            }
                            if (id_x>=0){
                                //IJ.log("x ok");
                                p.X=Double.parseDouble(split[id_x]);
                                nbmol++;
                            }
                            else{
                                p.X=0;
                            }
                            if (id_y>=0){
                                //IJ.log("y ok");
                                p.Y=Double.parseDouble(split[id_y]);
                            }
                            else{
                                p.Y=0;
                            }
                            if (id_z>=0){
                                //IJ.log("z ok");
                                p.Z=Double.parseDouble(split[id_z]);
                            }
                            else{
                                p.Z=0;
                            }
                            if (id_a>=0){
                                //IJ.log("a ok");
                                p.I=Double.parseDouble(split[id_a]);
                            }
                            else{
                                p.I=-1;
                            }
                            if (id_b>=0){
                                //IJ.log("b ok");
                                p.BG=Double.parseDouble(split[id_b]);
                            }
                            else{
                                p.BG=-1;
                            }
                            if (id_score>=0){
                                //IJ.log("score ok");
                                p.score=Double.parseDouble(split[id_score]);
                            }
                            else{
                                p.score=-1;
                            }
                            if (id_crlbx>=0){
                                p.crlb_X=Double.parseDouble(split[id_crlbx]);
                            }
                            else{
                                p.crlb_X=-1;
                            }
                            if (id_crlby>=0){
                                p.crlb_Y=Double.parseDouble(split[id_crlby]);
                            }
                            else{
                                p.crlb_Y=-1;
                            }
                            if (id_crlbz>=0){
                                p.crlb_Z=Double.parseDouble(split[id_crlbz]);
                            }
                            else{
                                p.crlb_Z=-1;
                            }
                            if (id_driftx>=0){
                                p.drift_X=Double.parseDouble(split[id_driftx]);
                            }
                            else{
                                p.drift_X=0;
                            }
                            if (id_drifty>=0){
                                p.drift_Y=Double.parseDouble(split[id_drifty]);
                            }
                            else{
                                p.drift_Y=0;
                            }
                            if (id_driftz>=0){
                                p.drift_Z=Double.parseDouble(split[id_driftz]);
                            }
                            else{
                                p.drift_Z=0;
                            }
                            if (id_occ>=0){
                                p.occurrence=(int)Double.parseDouble(split[id_occ]);
                            }
                            else{
                                p.occurrence=1;
                            }
                            for (int i=0;i<otherVariable.size();i++){
                                try{
                                    p.setValueOtherVariable(i,Double.parseDouble(split[otherVariableID.get(i)]));
                                }catch(Exception e){
                                    
                                }
                            }
                            
                            //IJ.log("value found "+p.id+"  "+p.frame+"  "+p.drift_X+"  "+p.drift_Y);
                        }
                        catch(Exception ee){p.exists=false;}
                        if (p.exists){
                            al.add(p);
                            if (p.frame>maxFrame){
                                maxFrame=p.frame;
                                
                                
                            }
                        }
                    }
                    
                }
                br.close();
                int [] frame = new int[maxFrame+1];
                for (int u=0;u<frame.length;u++){
                    frame[u]=-1;
                }
                for (int i=0;i<al.size();i++){
                    if (frame[al.get(i).frame]==-1){
                        this.fl.add(new FrameLocalization(al.get(i).frame));
                        frame[al.get(i).frame]=this.fl.size()-1;
                        this.fl.get(frame[al.get(i).frame]).loc.add(al.get(i));
                    }
                    else{
                        this.fl.get(frame[al.get(i).frame]).loc.add(al.get(i));
                    }
                    //IJ.write(""+this.fl.get(frame[al.get(i).frame]).loc.get(this.fl.get(frame[al.get(i).frame]).loc.size()-1).toString());
                    
                }
                
            IJ.log("Frame number: "+fl.size());
            IJ.log("Localization number: "+nbmol);
        }		
        catch (Exception e){
                System.out.println(e.toString());
        }
        
        
    }
    
    
    
    //load stack localization
    public void append(String path){
        
        int maxFrameAlreadyPresent=-1;
        for (int i=0;i<this.fl.size();i++){
            if (this.fl.get(i).numFrame>maxFrameAlreadyPresent){
                maxFrameAlreadyPresent=this.fl.get(i).numFrame;
            }
        }
        
        
        
        
        ArrayList<PLocalization> al = new ArrayList<PLocalization>();
        
        try{
                InputStream ips=new FileInputStream(path); 
                InputStreamReader ipsr=new InputStreamReader(ips);
                BufferedReader br=new BufferedReader(ipsr);
                String ligne;
                ArrayList<String> otherVariable=new ArrayList<String>();
                ArrayList<Integer> otherVariableID=new ArrayList<Integer>();
                
                int id_id=-1;
                int id_frame=-1;
                int id_x=-1;
                int id_y=-1;
                int id_z=-1;
                int id_a=-1;
                int id_b=-1;
                int id_score=-1;
                int id_crlbx=-1;
                int id_crlby=-1;
                int id_crlbz=-1;
                int id_driftx=-1;
                int id_drifty=-1;
                int id_driftz=-1;
                int id_occ=-1;
                //read first line
                int nbParam=0;
                String regex=",";
                if ((ligne=br.readLine())!=null){
                    PLocalization p= new PLocalization();
                    regex=p.getLabel_regex();
                    String [] split=ligne.split(regex);
                    nbParam=split.length;
                    for (int i=0;i<nbParam;i++){
                    //for (int i=nbParam-1;i>=0;i--){
                        //IJ.log("split-"+split[i]);
                        int xpt;
                        if (((xpt=split[i].indexOf(p.getLabel_id()))!=-1)&&(xpt<5)){
                            id_id=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_frame()))!=-1)&&(xpt<5)){
                            id_frame=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_x()))!=-1)&&(xpt<5)){
                            id_x=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_y()))!=-1)&&(xpt<5)){
                            id_y=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_z()))!=-1)&&(xpt<5)){
                            id_z=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_A()))!=-1)&&(xpt<5)){
                            id_a=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_B()))!=-1)&&(xpt<5)){
                            id_b=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_score()))!=-1)&&(xpt<5)){
                            id_score=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbX()))!=-1)&&(xpt<5)){
                            id_crlbx=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbY()))!=-1)&&(xpt<5)){
                            id_crlby=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbZ()))!=-1)&&(xpt<5)){
                            id_crlbz=i;
                            
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftX()))!=-1)&&(xpt<5)){
                            id_driftx=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftY()))!=-1)&&(xpt<5)){
                            id_drifty=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftZ()))!=-1)&&(xpt<5)){
                            id_driftz=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_occurrence()))!=-1)&&(xpt<5)){
                            id_occ=i;
                        }
                        else{
                            otherVariable.add(split[i]);
                            otherVariableID.add(i);
                        }
                    }
                }
                
                
                if (id_occ==-1){
                    this.merged=false;
                }
                else{
                    this.merged=true;
                }
                if ((id_driftx==-1)&&(id_drifty==-1)&&(id_driftz==-1)){
                    this.driftCorrected=false;
                }
                else{
                    this.driftCorrected=true;
                }
                //IJ.log("id frame "+id_frame+"  "+id_x+"  "+id_y+"  "+id_z);
                
                int nbmol=0;
                int nbframe=0;
                
                int maxFrame=0;
                int num=0;
                boolean firstLineValue=true;
                toto:while ((ligne=br.readLine())!=null){
                    
                    if (ligne.length()>3){
                        
                        if (firstLineValue){
                            String [] spliter=ligne.split(regex);
                            
                            if (spliter.length!=nbParam){
                                regex=" ";
                                spliter=ligne.split(regex);
                                
                                if (spliter.length!=nbParam){
                                    regex=";";
                                    spliter=ligne.split(regex);
                                    
                                    if (spliter.length!=nbParam){
                                        IJ.log("regex not recognized: it should be coma");
                                        regex=",";
                                    }
                                }
                            }
                            firstLineValue=false;
                        }
                        String [] split=ligne.split(regex);
                        //PLocalization p= new PLocalization(ligne);
                        PLocalization p= new PLocalization();
                        p.addListOfVariables(otherVariable);
                        p.exists=true;
                        try{
                            if (id_id>=0){
                                //IJ.log("id ok");
                                p.id=(int)Double.parseDouble(split[id_id]);
                                
                            }
                            else{
                                p.id=num++;
                            }
                            if (id_frame>=0){
                                //IJ.log("fr ok "+split[id_frame]);
                                p.frame=(int)Double.parseDouble(split[id_frame])+maxFrameAlreadyPresent+1;
                            }
                            else{
                                p.frame=maxFrameAlreadyPresent+1;
                            }
                            if (id_x>=0){
                                //IJ.log("x ok");
                                p.X=Double.parseDouble(split[id_x]);
                                nbmol++;
                            }
                            else{
                                p.X=0;
                            }
                            if (id_y>=0){
                                //IJ.log("y ok");
                                p.Y=Double.parseDouble(split[id_y]);
                            }
                            else{
                                p.Y=0;
                            }
                            if (id_z>=0){
                                //IJ.log("z ok");
                                p.Z=Double.parseDouble(split[id_z]);
                            }
                            else{
                                p.Z=0;
                            }
                            if (id_a>=0){
                                //IJ.log("a ok");
                                p.I=Double.parseDouble(split[id_a]);
                            }
                            else{
                                p.I=0;
                            }
                            if (id_b>=0){
                                //IJ.log("b ok");
                                p.BG=Double.parseDouble(split[id_b]);
                            }
                            else{
                                p.BG=0;
                            }
                            if (id_score>=0){
                                //IJ.log("score ok");
                                p.score=Double.parseDouble(split[id_score]);
                            }
                            else{
                                p.score=0;
                            }
                            if (id_crlbx>=0){
                                p.crlb_X=Double.parseDouble(split[id_crlbx]);
                            }
                            else{
                                p.crlb_X=0;
                            }
                            if (id_crlby>=0){
                                p.crlb_Y=Double.parseDouble(split[id_crlby]);
                            }
                            else{
                                p.crlb_Y=0;
                            }
                            if (id_crlbz>=0){
                                p.crlb_Z=Double.parseDouble(split[id_crlbz]);
                            }
                            else{
                                p.crlb_Z=0;
                            }
                            if (id_driftx>=0){
                                p.drift_X=Double.parseDouble(split[id_driftx]);
                            }
                            else{
                                p.drift_X=0;
                            }
                            if (id_drifty>=0){
                                p.drift_Y=Double.parseDouble(split[id_drifty]);
                            }
                            else{
                                p.drift_Y=0;
                            }
                            if (id_driftz>=0){
                                p.drift_Z=Double.parseDouble(split[id_driftz]);
                            }
                            else{
                                p.drift_Z=0;
                            }
                            if (id_occ>=0){
                                p.occurrence=(int)Double.parseDouble(split[id_occ]);
                            }
                            else{
                                p.occurrence=1;
                            }
                            for (int i=0;i<otherVariable.size();i++){
                                p.setValueOtherVariable(i,Double.parseDouble(split[otherVariableID.get(i)]));
                            }
                            //IJ.log("value found "+p.id+"  "+p.X+"  "+p.Y+"   "+p.frame);
                            //IJ.log("value found "+p.id+"  "+p.frame+"  "+p.drift_X+"  "+p.drift_Y);
                        }
                        catch(Exception ee){p.exists=false;}
                        if (p.exists){
                            al.add(p);
                            if (p.frame>maxFrame){
                                maxFrame=p.frame;
                            }
                        }
                    }

                }
                br.close();
                int [] frame = new int[maxFrame+1];
                for (int u=0;u<frame.length;u++){
                    frame[u]=-1;
                }
                for (int i=0;i<al.size();i++){
                    if (frame[al.get(i).frame]==-1){
                        nbframe++;
                        this.fl.add(new FrameLocalization(al.get(i).frame));
                        frame[al.get(i).frame]=this.fl.size()-1;
                        this.fl.get(frame[al.get(i).frame]).loc.add(al.get(i));
                    }
                    else{
                        this.fl.get(frame[al.get(i).frame]).loc.add(al.get(i));
                    }
                    //IJ.write(""+this.fl.get(frame[al.get(i).frame]).loc.get(this.fl.get(frame[al.get(i).frame]).loc.size()-1).toString());
                    
                }
                
                
            //IJ.log("Frame number: "+nbframe);
            IJ.log("Localization number: "+nbmol);
                
        }		
        catch (Exception e){
                System.out.println(e.toString());
        }
        
        
        
    }
    
    
    
    
    
    
    //load stack localization
    public void fuse(String path){
        
        
        
        int maxFrameAlreadyPresent=-1;
        for (int i=0;i<this.fl.size();i++){
            if (this.fl.get(i).numFrame>maxFrameAlreadyPresent){
                maxFrameAlreadyPresent=this.fl.get(i).numFrame;
            }
        }
        
        
        
        ArrayList<PLocalization> al = new ArrayList<PLocalization>();
        
        try{
                InputStream ips=new FileInputStream(path); 
                InputStreamReader ipsr=new InputStreamReader(ips);
                BufferedReader br=new BufferedReader(ipsr);
                String ligne;
                ArrayList<String> otherVariable=new ArrayList<String>();
                ArrayList<Integer> otherVariableID=new ArrayList<Integer>();
                
                int id_id=-1;
                int id_frame=-1;
                int id_x=-1;
                int id_y=-1;
                int id_z=-1;
                int id_a=-1;
                int id_b=-1;
                int id_score=-1;
                int id_crlbx=-1;
                int id_crlby=-1;
                int id_crlbz=-1;
                int id_driftx=-1;
                int id_drifty=-1;
                int id_driftz=-1;
                int id_occ=-1;
                //read first line
                int nbParam=0;
                String regex=",";
                if ((ligne=br.readLine())!=null){
                    PLocalization p= new PLocalization();
                    regex=p.getLabel_regex();
                    String [] split=ligne.split(regex);
                    nbParam=split.length;
                    for (int i=0;i<nbParam;i++){
                    //for (int i=nbParam-1;i>=0;i--){
                        //IJ.log("split-"+split[i]);
                        int xpt;
                        if (((xpt=split[i].indexOf(p.getLabel_id()))!=-1)&&(xpt<5)){
                            id_id=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_frame()))!=-1)&&(xpt<5)){
                            id_frame=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_x()))!=-1)&&(xpt<5)){
                            id_x=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_y()))!=-1)&&(xpt<5)){
                            id_y=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_z()))!=-1)&&(xpt<5)){
                            id_z=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_A()))!=-1)&&(xpt<5)){
                            id_a=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_B()))!=-1)&&(xpt<5)){
                            id_b=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_score()))!=-1)&&(xpt<5)){
                            id_score=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbX()))!=-1)&&(xpt<5)){
                            id_crlbx=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbY()))!=-1)&&(xpt<5)){
                            id_crlby=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_crlbZ()))!=-1)&&(xpt<5)){
                            id_crlbz=i;
                            
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftX()))!=-1)&&(xpt<5)){
                            id_driftx=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftY()))!=-1)&&(xpt<5)){
                            id_drifty=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_driftZ()))!=-1)&&(xpt<5)){
                            id_driftz=i;
                        }
                        else if (((xpt=split[i].indexOf(p.getLabel_occurrence()))!=-1)&&(xpt<5)){
                            id_occ=i;
                        }
                        else{
                            otherVariable.add(split[i]);
                            otherVariableID.add(i);
                        }
                        
                    }
                }
                
                
                if (id_occ==-1){
                    this.merged=false;
                }
                else{
                    this.merged=true;
                }
                if ((id_driftx==-1)&&(id_drifty==-1)&&(id_driftz==-1)){
                    this.driftCorrected=false;
                }
                else{
                    this.driftCorrected=true;
                }
                //IJ.log("id frame "+id_frame+"  "+id_x+"  "+id_y+"  "+id_z);
                
                int nbmol=0;
                int nbframe=0;
                
                int maxFrameSecond=0;
                int num=0;
                boolean firstLineValue=true;
                toto:while ((ligne=br.readLine())!=null){
                    
                    if (ligne.length()>3){
                        
                        if (firstLineValue){
                            String [] spliter=ligne.split(regex);
                            
                            if (spliter.length!=nbParam){
                                regex=" ";
                                spliter=ligne.split(regex);
                                
                                if (spliter.length!=nbParam){
                                    regex=";";
                                    spliter=ligne.split(regex);
                                    
                                    if (spliter.length!=nbParam){
                                        IJ.log("regex not recognized: it should be coma");
                                        regex=",";
                                    }
                                }
                            }
                            firstLineValue=false;
                        }
                        String [] split=ligne.split(regex);
                        //PLocalization p= new PLocalization(ligne);
                        PLocalization p= new PLocalization();
                        p.addListOfVariables(otherVariable);
                        p.exists=true;
                        try{
                            if (id_id>=0){
                                //IJ.log("id ok");
                                p.id=(int)Double.parseDouble(split[id_id]);
                                
                            }
                            else{
                                p.id=num++;
                            }
                            if (id_frame>=0){
                                //IJ.log("fr ok "+split[id_frame]);
                                p.frame=(int)Double.parseDouble(split[id_frame]);
                                if (maxFrameSecond<p.frame){
                                    maxFrameSecond=p.frame;
                                }
                            }
                            else{
                                p.frame=1;
                                if (maxFrameSecond<p.frame){
                                    maxFrameSecond=p.frame;
                                }
                            }
                            if (id_x>=0){
                                //IJ.log("x ok");
                                p.X=Double.parseDouble(split[id_x]);
                                nbmol++;
                            }
                            else{
                                p.X=0;
                            }
                            if (id_y>=0){
                                //IJ.log("y ok");
                                p.Y=Double.parseDouble(split[id_y]);
                            }
                            else{
                                p.Y=0;
                            }
                            if (id_z>=0){
                                //IJ.log("z ok");
                                p.Z=Double.parseDouble(split[id_z]);
                            }
                            else{
                                p.Z=0;
                            }
                            if (id_a>=0){
                                //IJ.log("a ok");
                                p.I=Double.parseDouble(split[id_a]);
                            }
                            else{
                                p.I=0;
                            }
                            if (id_b>=0){
                                //IJ.log("b ok");
                                p.BG=Double.parseDouble(split[id_b]);
                            }
                            else{
                                p.BG=0;
                            }
                            if (id_score>=0){
                                //IJ.log("score ok");
                                p.score=Double.parseDouble(split[id_score]);
                            }
                            else{
                                p.score=0;
                            }
                            if (id_crlbx>=0){
                                p.crlb_X=Double.parseDouble(split[id_crlbx]);
                            }
                            else{
                                p.crlb_X=0;
                            }
                            if (id_crlby>=0){
                                p.crlb_Y=Double.parseDouble(split[id_crlby]);
                            }
                            else{
                                p.crlb_Y=0;
                            }
                            if (id_crlbz>=0){
                                p.crlb_Z=Double.parseDouble(split[id_crlbz]);
                            }
                            else{
                                p.crlb_Z=0;
                            }
                            if (id_driftx>=0){
                                p.drift_X=Double.parseDouble(split[id_driftx]);
                            }
                            else{
                                p.drift_X=0;
                            }
                            if (id_drifty>=0){
                                p.drift_Y=Double.parseDouble(split[id_drifty]);
                            }
                            else{
                                p.drift_Y=0;
                            }
                            if (id_driftz>=0){
                                p.drift_Z=Double.parseDouble(split[id_driftz]);
                            }
                            else{
                                p.drift_Z=0;
                            }
                            if (id_occ>=0){
                                p.occurrence=(int)Double.parseDouble(split[id_occ]);
                            }
                            else{
                                p.occurrence=1;
                            }
                            for (int i=0;i<otherVariable.size();i++){
                                p.setValueOtherVariable(i,Double.parseDouble(split[otherVariableID.get(i)]));
                            }
                            //IJ.log("value found "+p.id+"  "+p.X+"  "+p.Y+"   "+p.frame);
                            //IJ.log("value found "+p.id+"  "+p.frame+"  "+p.drift_X+"  "+p.drift_Y);
                        }
                        catch(Exception ee){p.exists=false;}
                        if (p.exists){
                            al.add(p);
                        }
                    }

                }
                br.close();
                //IJ.log("frame numb "+maxFrameAlreadyPresent+"  "+maxFrameSecond);
                int [] frameNumberPosition = new int[1+Math.max(maxFrameAlreadyPresent,maxFrameSecond)];
                for (int i=0;i<frameNumberPosition.length;i++){
                    frameNumberPosition[i]=-1;
                }

                for (int i=0;i<this.fl.size();i++){
                    frameNumberPosition[this.fl.get(i).numFrame]=i;
                }
                
                for (int i=0;i<al.size();i++){
                    PLocalization p=al.get(i);
                    if (frameNumberPosition[p.frame]==-1){
                        this.fl.add(new FrameLocalization(p.frame));
                        frameNumberPosition[p.frame]=this.fl.size()-1;
                    }
                    this.fl.get(frameNumberPosition[p.frame]).loc.add(p);
                    //IJ.log("add "+p.frame);
                }
                
                for (int i=0;i<fl.size();i++){
                    IJ.log("frame "+fl.get(i).numFrame+"  "+fl.get(i).loc.size());
                }
                
            
                
        }		
        catch (Exception e){
                System.out.println(e.toString());
        }
        
        
        
    }
    
    
    public String save(String path){
        if (path.length()>4){
            if (path.toLowerCase().endsWith(".csv")){
                this.saveCSV(path);
            }
            else if (path.toLowerCase().endsWith(".smlm")){
                this.saveSMLM(path);
            }
            else{
                path=path+".csv";
                this.saveCSV(path);
            }
            return path;
        }
        else{
            IJ.log("file not saved. Please, select a path with at least 4 characters");
            return null;
        }
        
        
    }
    
    
    public void saveCSV(String path){
        
        
        try {


            PrintWriter sortie;


                sortie = new PrintWriter(new FileWriter(path, false));

                //firstline
                loop1:for (int i=0;i<fl.size();i++){
                    for (int j=0;j<fl.get(i).loc.size();j++){
                        
                        sortie.println(""+fl.get(i).loc.get(j).toStringName());
                        break loop1;
                    }
                }


                for (int i=0;i<fl.size();i++){
                    for (int j=0;j<fl.get(i).loc.size();j++){
                        if (fl.get(i).loc.get(j).exists){
                            sortie.println(""+fl.get(i).loc.get(j).toString());
                        }
                    }
                }


                sortie.close();




        } catch (Exception e) {
                e.printStackTrace();
                IJ.log("error saving process");
        }
        
        
    }
    
    
    public void saveSMLM(String path){
        /*
        try {
        


            SMLMio smlmio = new SMLMio();
            SMLMio.Format format = smlmio.new Format();

            String name="ZOLA_localization_table";
            String nameBin="ZOLA_localization_table.bin";
            String nameFormat="smlm-table(binary)";

            ArrayList<String> head = new ArrayList();
            int variableNumber=0;
            loop1:for (int i=0;i<fl.size();i++){
                for (int j=0;j<fl.get(i).loc.size();j++){
                    variableNumber=fl.get(i).loc.get(j).getNumberVariable();
                    for (int v=0;v<variableNumber;v++){
                        String header=fl.get(i).loc.get(j).getLabel(v);

                        if (header.equals(PLocalization.getLabel_x())){
                            header="x";
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_y())){
                            header="y";
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_z())){
                            header="z";
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_crlbX())){
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_crlbY())){
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_crlbZ())){
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_driftX())){
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_driftY())){
                            format.addUnit(header, "nm");
                        }
                        else if (header.equals(PLocalization.getLabel_driftZ())){
                            format.addUnit(header, "nm");
                        }
                        else{
                            format.addUnit(header, "a.u.");
                        }
                        format.addShape(header,1);
                        if (header.equals(PLocalization.getLabel_frame())||header.equals(PLocalization.getLabel_id())||header.equals(PLocalization.getLabel_occurrence())){
                            format.addType(header, (SMLMio.type_uint32));
                        }
                        else{
                            format.addType(header, (SMLMio.type_float32));
                        }
                        head.add(header);
                    }

                    break loop1;
                }
            }


            format.addMeta_data("name", nameFormat);
            format.addMeta_data("mode", "binary");
            format.addMeta_data("type", "table");
            format.addMeta_data("extension", ".bin");
            format.addMeta_data("columns", variableNumber);

            SMLMio.Table table = smlmio.new Table(format);

            int elNumber=0;

            for (int i=0;i<fl.size();i++){
                elNumber+=fl.get(i).loc.size();
            }
            
            double [][] tableData=new double[elNumber][variableNumber];
            
            for (int c=0;c<variableNumber;c++) {
                int [] index = new int [1];
                index[0]=c;
                table.setIndex(head.get(c),index);
            }
            
            for (int v=0;v<variableNumber;v++){
                String header=head.get(v);
                
                
                
                double min=Double.POSITIVE_INFINITY;
                double max=Double.NEGATIVE_INFINITY;
                double avg=0;
                double value;
                for (int i=0,t=0;i<fl.size();i++){
                    for (int j=0;j<fl.get(i).loc.size();j++,t++){
                        value=(double)fl.get(i).loc.get(j).getValueVariable(v);
                        tableData[t][v]=value;
                        if (value<min){
                            min=value;
                        }
                        if (value>max){
                            max=value;
                        }
                        avg+=value;
                    }
                }
                avg/=(float)elNumber;
                table.setMin(header,min);
                table.setMax(header,max);
                table.setAverage(header,avg);
                


            }
            table.addTable(tableData);

            table.addMeta_data("name", nameBin);
            table.addMeta_data("format", nameFormat);
            table.addMeta_data("rows", elNumber);
            table.addMeta_data("type", "table");





            SMLMio.Data data = smlmio.new Data();
            data.format_version="0.2";
            data.addTable(nameBin, table);

            smlmio.export_smlm(path,data);
            
        } catch (Exception e) {
                e.printStackTrace();
                IJ.log("error saving process");
        }
        
        */
    }
    
    
    public void saveOnTheFly(String path){
        
        if (path!=null){
            save=new Save(path);
        }
        
    }
    
    
    public void stopSaveOnTheFly(){
        save.stopSave();
        
    }
    
    
    
    
    
    
    
    class Save extends Thread{
        
        ReentrantLock lock = new ReentrantLock();
        
        String path;
        boolean stopped=false;
        
        PrintWriter sortie;
        
        int treated=0;
                    
        int flSize;
        
        boolean firstlineOK=false;
        
        Save(){}
        
        Save(String path){
            
            this.path=path;
            
            if (path!=null){
                
                try {
                    int nn=(path.lastIndexOf("."));
                    if (nn>path.length()-5){
                        path=path.substring(0,nn);
                    }
                    path=path+".csv";




                        sortie = new PrintWriter(new FileWriter(path, false));


                        start();

                } catch (Exception e) {
                        e.printStackTrace();
                }
            }
            
            
        }
        
        
        public void stopSave(){
            
            
            try{
                lock.lock();
                {
                    lock.lock();
                    

                    flSize=fl.size();
                    //firstline
                    if (!firstlineOK){
                        loop1:for (int i=treated;i<flSize;i++){
                            for (int j=0;j<fl.get(i).loc.size();j++){
                                sortie.println(""+fl.get(i).loc.get(j).toStringName());
                                firstlineOK=true;
                                break loop1;
                            }
                        }
                    }


                    for (int i=treated;i<flSize;i++){
                        for (int j=0;j<fl.get(i).loc.size();j++){
                            sortie.println(""+fl.get(i).loc.get(j).toString());
                        }
                    }

                    treated=flSize;
                }
                sortie.close();
                stopped=true;
            }
            finally{
                lock.unlock();
            }
        }
        
        
        public void run(){
            
            
            firstlineOK=false;

            while (!stopped){
                
                try{
                    Thread.sleep(1000);
                } catch (Exception ey) {
                }
                
                try{
                    lock.lock();
                    
                    if (!stopped){
                        flSize=fl.size();
                        //firstline
                        if (!firstlineOK){
                            loop1:for (int i=treated;i<flSize;i++){
                                for (int j=0;j<fl.get(i).loc.size();j++){
                                    sortie.println(""+fl.get(i).loc.get(j).toStringName());
                                    firstlineOK=true;
                                    break loop1;
                                }
                            }
                        }


                        for (int i=treated;i<flSize;i++){
                            for (int j=0;j<fl.get(i).loc.size();j++){
                                sortie.println(""+fl.get(i).loc.get(j).toString());
                            }
                        }

                        treated=flSize;
                    }
                }
                finally{
                    lock.unlock();
                }
                
                
                
                
                
                

            }
                    
                    
        }
        
    }
    
    
}
