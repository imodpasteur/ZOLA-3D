/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.pasteur.imagej.data;


import org.pasteur.imagej.data.PLocalization;
import ij.IJ;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;



public class FrameLocalization{ 
    
    public ArrayList<PLocalization> loc = new ArrayList<PLocalization>();
    public int numFrame;
    public FrameLocalization(int numFrame){
        this.numFrame=numFrame;
    }
    
    
    
}

