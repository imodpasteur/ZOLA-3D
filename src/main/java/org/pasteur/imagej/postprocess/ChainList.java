/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.postprocess;

import org.pasteur.imagej.data.*;
import ij.IJ;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/**
 *
 * @author benoit
 */
public class ChainList {
    
        public int framePrevious;
        public int frameNext;
        public int locPrevious;
        public int locNext;
        public double distPrevious;
        public double distNext;
        public boolean foundNext=false;
        public boolean foundPrevious=false;
        
        ChainList(){}
        
        
        
    
    
    
    
    public static ChainList [][] getChainList(StackLocalization sl1,double maxDistanceMergingXY,double maxDistanceMergingZ,int offFrame){
        
        
        
        
        int [][] idFrame=new int[sl1.fl.size()][2];
        for (int i=0;i<sl1.fl.size();i++){
            idFrame[i][0]=i;
            idFrame[i][1]=sl1.fl.get(i).numFrame;
        }
        Arrays.sort(idFrame, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                return ((Integer) o1[1]).compareTo(o2[1]);
            }
        });
        
        ChainList [][] condens=new ChainList[sl1.fl.size()][];
        for (int i=0;i<sl1.fl.size();i++){
            condens[i]=new ChainList[sl1.fl.get(i).loc.size()];
            for (int ii=0;ii<condens[i].length;ii++){
                condens[i][ii]=new ChainList();
            }
        }
        
        int totNumber=0;
        int condenseNumber=0;
        
        
        double distXY,dXY;
        double distZ,dZ;
        double dist,d;
        int [] idclosest = new int[2];
        double xx,yy,zz,x,y,z;
        ArrayList<int []> positionToRecompute = new ArrayList<int[]>();
        for (int i=0;i<idFrame.length-1;i++){
            //IJ.showProgress(.5*((float)i/(float)idFrame.length));
            int id=idFrame[i][0];
            for (int ii=0;ii<sl1.fl.get(id).loc.size();ii++){
                totNumber++;
                PLocalization p=sl1.fl.get(id).loc.get(ii);
                x=p.X;
                y=p.Y;
                z=p.Z;
                distXY=Double.MAX_VALUE;
                distZ=Double.MAX_VALUE;
                dist=Double.MAX_VALUE;
                boolean found=false;
                for (int dec=0;dec<=offFrame&&!found;dec++){
                    if (i+1+dec<idFrame.length){
                        if (idFrame[i+1+dec][1]-idFrame[i][1]<=1+offFrame){//added
                            int idNext=idFrame[i+1+dec][0];
                            for (int jj=0;jj<sl1.fl.get(idNext).loc.size();jj++){
                                PLocalization pp=sl1.fl.get(idNext).loc.get(jj);
                                xx=pp.X;
                                yy=pp.Y;
                                zz=pp.Z;
                                dXY=Math.sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy));
                                if (dXY<=maxDistanceMergingXY){
                                    dZ=Math.sqrt((z-zz)*(z-zz));
                                    if ((dZ<=maxDistanceMergingZ)||(maxDistanceMergingZ<0)){
                                        d=Math.sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz));
                                        if ((dist>d)){
                                            if (!condens[idNext][jj].foundPrevious){
                                                dist=d;
                                                distZ=dZ;
                                                distXY=dXY;
                                                idclosest[0]=i+1;
                                                idclosest[1]=jj;
                                                condens[idNext][jj].foundPrevious=true;
                                                condens[idNext][jj].framePrevious=id;
                                                condens[idNext][jj].locPrevious=ii;
                                                condens[idNext][jj].distPrevious=dist;

                                                condens[id][ii].foundNext=true;
                                                condens[id][ii].frameNext=idNext;
                                                condens[id][ii].locNext=jj;
                                                condens[id][ii].distNext=dist;
                                                found=true;
                                                condenseNumber++;
                                            }
                                            else if (condens[idNext][jj].distPrevious>dist){
                                                int [] pToRecompute = new int[2];
                                                pToRecompute[0]=condens[idNext][jj].framePrevious;
                                                pToRecompute[1]=condens[idNext][jj].locPrevious;
                                                positionToRecompute.add(pToRecompute);
                                                condens[pToRecompute[0]][pToRecompute[1]].foundNext=false;
                                                dist=d;
                                                distZ=dZ;
                                                distXY=dXY;
                                                idclosest[0]=i+1;
                                                idclosest[1]=jj;
                                                condens[idNext][jj].foundPrevious=true;
                                                condens[idNext][jj].framePrevious=id;
                                                condens[idNext][jj].locPrevious=ii;
                                                condens[idNext][jj].distPrevious=dist;

                                                condens[id][ii].foundNext=true;
                                                condens[id][ii].frameNext=idNext;
                                                condens[id][ii].locNext=jj;
                                                condens[id][ii].distNext=dist;
                                                found=true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
            }
        }
        
        
        //IJ.log("Condensation step 1 : "+condenseNumber+"/"+totNumber);
        
        //IJ.showProgress(.5);
        
        for (int i=0;i<positionToRecompute.size();i++){
            //IJ.showProgress(.5+.1*(float)i/(float)positionToRecompute.size());
            int id=positionToRecompute.get(i)[0];
            int ii=positionToRecompute.get(i)[1];
            
            PLocalization p=sl1.fl.get(id).loc.get(ii);
            x=p.X;
            y=p.Y;
            z=p.Z;
            distXY=Double.MAX_VALUE;
            distZ=Double.MAX_VALUE;
            dist=Double.MAX_VALUE;
            boolean found=false;
            for (int dec=0;dec<=offFrame&&!found;dec++){
                if (i+1+dec<idFrame.length){
                    if (idFrame[i+1+dec][1]-idFrame[i][1]<=1+offFrame){//added
                        int idNext=idFrame[i+1+dec][0];
                        for (int jj=0;jj<sl1.fl.get(idNext).loc.size();jj++){
                            PLocalization pp=sl1.fl.get(idNext).loc.get(jj);
                            xx=pp.X;
                            yy=pp.Y;
                            zz=pp.Z;
                            dXY=Math.sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy));
                            if (dXY<=maxDistanceMergingXY){
                                dZ=Math.sqrt((z-zz)*(z-zz));
                                if ((dZ<=maxDistanceMergingZ)||(maxDistanceMergingZ<0)){
                                    d=Math.sqrt((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz));
                                    if ((dist>d)){
                                        if (!condens[idNext][jj].foundPrevious){
                                            dist=d;
                                            distZ=dZ;
                                            distXY=dXY;
                                            idclosest[0]=i+1;
                                            idclosest[1]=jj;
                                            condens[idNext][jj].foundPrevious=true;
                                            condens[idNext][jj].framePrevious=id;
                                            condens[idNext][jj].locPrevious=ii;
                                            condens[idNext][jj].distPrevious=dist;

                                            condens[id][ii].foundNext=true;
                                            condens[id][ii].frameNext=idNext;
                                            condens[id][ii].locNext=jj;
                                            condens[id][ii].distNext=dist;
                                            found=true;
                                        }
                                        else if (condens[idNext][jj].distPrevious>dist){
                                            int [] pToRecompute = new int[2];
                                            pToRecompute[0]=condens[idNext][jj].framePrevious;
                                            pToRecompute[1]=condens[idNext][jj].locPrevious;
                                            positionToRecompute.add(pToRecompute);
                                            condens[pToRecompute[0]][pToRecompute[1]].foundNext=false;
                                            dist=d;
                                            distZ=dZ;
                                            distXY=dXY;
                                            idclosest[0]=i+1;
                                            idclosest[1]=jj;
                                            condens[idNext][jj].foundPrevious=true;
                                            condens[idNext][jj].framePrevious=id;
                                            condens[idNext][jj].locPrevious=ii;
                                            condens[idNext][jj].distPrevious=dist;

                                            condens[id][ii].foundNext=true;
                                            condens[id][ii].frameNext=idNext;
                                            condens[id][ii].locNext=jj;
                                            condens[id][ii].distNext=dist;
                                            found=true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            
            
        }
        
        //IJ.log("conflicts managed "+positionToRecompute.size());
            
        //IJ.showProgress(.6);
                
        
        
        return condens;
    }
        
        
    }

