/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.control;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SearchProblem;

/**
 *
 * This class outputs info (structures, energies) on specified conformations 
 * in a conformation file
 * 
 * @author mhall44
 */
public class ConfInfo {
    
    ConfigFileParser cfp;
    int numRes;
    boolean outputPDBs;
    int numConfsToOutput;
    
    boolean doFullOutput;//If true, describe each conf in full, not just energy
    //(written to stdout and a conf file as in findGMEC)
    String outputConfFileName;//If doFullOutput, will write to this file
    
    
    public ConfInfo(ConfigFileParser cfp){
        this.cfp = cfp;
        numRes = cfp.getFlexRes().size();
        outputPDBs = cfp.params.getBool("outputPDBs");
        numConfsToOutput = cfp.params.getInt("NumConfsToOutput");
        
        doFullOutput = cfp.params.getBool("doFullOutput");
        outputConfFileName = cfp.params.getRunSpecificFileName("OutputConfFileName", ".outputconfs.txt");
    }
    
    
    public void outputConfInfo(){
        
        String confFileName = cfp.params.getRunSpecificFileName("CONFFILENAME", ".confs.txt");
        try (BufferedReader br = new BufferedReader(new FileReader(confFileName))) {
            
            SearchProblem searchProb = cfp.getSearchProblem();
            
            if(searchProb.useERef || doFullOutput)
                searchProb.loadEnergyMatrix();
            
            System.out.println("OUTPUTTING CONF INFO...");
            
            ConfPrinter confPrinter = null;
            if(outputPDBs)
                setupPDBDir();
            if(doFullOutput)
                confPrinter = new ConfPrinter(searchProb, outputConfFileName, false);
            
            double bestESoFar = Double.POSITIVE_INFINITY;

            for(int confCount=0; confCount<numConfsToOutput; confCount++){
                //Each line in conf file describes a conf
                String confLine = br.readLine();
                int[] conf = parseConfLine(confLine);
                double confE = searchProb.minimizedEnergy(conf);
                
                bestESoFar = Math.min(confE,bestESoFar);

                if(doFullOutput){
                	
                    //let's just show the standard lower bound here
                    EnergiedConf econf = new EnergiedConf(conf, searchProb.lowerBound(conf), confE);
                    
                    confPrinter.printConf(econf);
                    System.out.println("ENUMERATING CONFORMATION:");
                    System.out.print(confPrinter.getConfReport(econf));
                }
                else
                    System.out.println("Conf "+confCount+" energy: "+confE);

                if(outputPDBs){
                    searchProb.outputMinimizedStruct(conf, "pdbs/"+searchProb.name+".STRUCT"+confCount+".pdb");
                }
            }
            
            if(doFullOutput)
                confPrinter.closeConfFile();
            
            System.out.println("DONE OUTPUTTING CONF INFO.");
        }
        catch(IOException e){
            throw new RuntimeException(e);
        }
    }
    
    
    void setupPDBDir(){
        //Make sure there is a pdbs directory, creating it if needed
        File pdbsHandle = new File("pdbs");
        if(!pdbsHandle.isDirectory()){
            if(pdbsHandle.exists())
                throw new RuntimeException("ERROR: Can't create pdbs directory (there is already a file called pdbs)");
            else
                pdbsHandle.mkdir();
        }
        //else write to the existing directory (note: this may overwrite files there,
        //which we assume is desired behavior)
    }
    
    int[] parseConfLine(String confLine){
        //Parse out the conformation (list of RCs) from a line in the conf file
        //format: "0 CONF: 1 2 3 4 RESTYPES..."
        //where RCs for the four flexible residues are 1,2,3,4
        
        if(confLine==null)
            throw new RuntimeException("ERROR: Conformation file ends prematurely");
        
        StringTokenizer st = new StringTokenizer(confLine," ");
        if(st.countTokens()<numRes+2)
            throw confLineException(confLine);
        
        st.nextToken();
        if(!st.nextToken().equalsIgnoreCase("CONF:"))//check for CONF: token
            throw confLineException(confLine);
        
        int[] ans = new int[numRes];
        try {
            for(int res=0; res<numRes; res++)
                ans[res] = Integer.valueOf(st.nextToken());
        }
        catch(NumberFormatException e){
            throw confLineException(confLine);
        }
        
        return ans;
    }
    
    
    RuntimeException confLineException(String confLine){
        //an exception in parsing a line in the conformation file
        throw new RuntimeException( "ERROR: Expected conformation file line in format '"
                + "conf_num CONF: rc_number1 ... rc_number"+numRes+
                " ...' but found this: " + confLine);
    }
    
    
}
