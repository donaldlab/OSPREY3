/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.gmec;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.control.ConfigFileParser;

/**
 *
 * This is for when you have a bunch of sequences and you want the GMEC for each
 * Commands: -c KStar.cfg findSeqGMECs System.cfg DEE.cfg sequences.mut
 * 
 * @author mhall44
 */
public class SeqGMECFinder {
    
    ConfigFileParser cfp;
    File mutFile;
    
    public SeqGMECFinder(ConfigFileParser cfp){
        this.cfp = cfp;
        this.mutFile = cfp.params.getFile("MutFile");
    }
    
    
    public void calcAllSeqGMECs(){
        try{
            FileInputStream is = new FileInputStream(mutFile);
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
            
            int seqNum=0;
            
            for(String curLine=bufread.readLine(); curLine!=null; curLine=bufread.readLine()){
                StringTokenizer st = new StringTokenizer(curLine);
                int numPos = st.countTokens();
                String seq[] = new String[numPos];
                for(int pos=0; pos<numPos; pos++)//take only the 3-letter AA type
                    seq[pos] = st.nextToken().substring(0,3);
                
                calcSeqGMEC(seq, seqNum);
                seqNum++;
            }
            
            bufread.close();
        }
        catch(IOException ex){
            throw new RuntimeException(ex);
        }
    }
    
    
    private double calcSeqGMEC(String[] AATypes, int seqNum){
        //Calculate the GMEC for sequence AATypes
        //This is based on calcStateGMEC from COMETSDoer
        
        System.out.println();
        System.out.print("CALCULATING GMEC FOR SEQUENCE "+seqNum+":");
        for(String aaType : AATypes)
            System.out.print(" "+aaType);
        System.out.println();
        System.out.println();
        
        ConfigFileParser cfp = new ConfigFileParser(this.cfp);
        cfp.loadData();
        
        if(cfp.getFlexRes().size()!=AATypes.length){
            throw new RuntimeException("ERROR: Wrong number of residues in sequence, "+AATypes.length+
                    " instead of "+cfp.getFlexRes().size());
        }
        
        //set up for particular AA types
        String panSeqRunName = cfp.params.getValue("RUNNAME");
        cfp.params.setValue("RUNNAME", panSeqRunName+"_SEQ"+seqNum);
        //also want to restrict to only this seq: addWT = false
        cfp.params.setValue("ADDWT", "false");        
        
        //matching format to ConfigFileParser.getAllowedAAs
        int posCount = 0;
        
        for(int str=0; str<10; str++){
            ArrayList<String> resAllowedRecords = cfp.params.searchParams("RESALLOWED"+str);
            int numRecordsInStrand = resAllowedRecords.size();
            
            //must go through residues in numerical order
            for(int recNum=0; recNum<numRecordsInStrand; recNum++){
                String param = "RESALLOWED" + str + "_" + recNum;
                cfp.params.setValue(param, AATypes[posCount]);
                posCount++;
            }
        }
        
        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        double E = gf.calcGMEC().get(0).getEnergy();
        System.out.print("DONE CALCULATING GMEC FOR SEQUENCE "+seqNum);
        return E;
    }
    
}
