/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.SearchProblem;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

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
    
    public ConfInfo(ConfigFileParser cfp){
        this.cfp = cfp;
        numRes = cfp.getFlexRes().size();
        outputPDBs = cfp.params.getBool("outputPDBs");
        numConfsToOutput = cfp.params.getInt("NumConfsToOutput");
    }
    
    
    public void outputConfInfo(){
        
        try {
            String confFileName = cfp.params.getRunSpecificFileName("CONFFILENAME", ".confs.txt");
            BufferedReader br = new BufferedReader(new FileReader(confFileName));
            SearchProblem searchProb = cfp.getSearchProblem();

            System.out.println("OUTPUTTING CONF INFO...");

            if(outputPDBs)
                setupPDBDir();

            for(int confCount=0; confCount<numConfsToOutput; confCount++){
                //Each line in conf file describes a conf
                String confLine = br.readLine();
                int[] conf = parseConfLine(confLine);
                double confE = searchProb.minimizedEnergy(conf);
                System.out.println("Conf "+confCount+" energy: "+confE);

                if(outputPDBs){
                    searchProb.outputMinimizedStruct(conf, "pdbs/"+searchProb.name+".STRUCT"+confCount+".pdb");
                }
            }

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
