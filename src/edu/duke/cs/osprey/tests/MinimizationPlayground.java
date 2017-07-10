/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;

/**
 *
 * Playing with minimization
 * 
 * @author mhall44
 */
public class MinimizationPlayground {
    
    public static void main(String args[]){
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);//args are configuration files
        cfp.loadData();
        SearchProblem sp = cfp.getSearchProblem();
        
        int conf1[] = new int[15];//start with all wt rots
        conf1[12] = 5;//TRP
        System.out.println("Energy: "+sp.minimizedEnergy(conf1));
        sp.outputMinimizedStruct(conf1, "VOXMIN.pdb");
    }
    
}
