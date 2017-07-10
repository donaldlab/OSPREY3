/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

/**
 *
 * A space of sequences that we want to limit a GMEC search to
 * Defined by a mut file
 * 
 * @author mhall44
 */
public class GMECMutSpace {
    
    //We'll want to look up what residue types can follow given partial sequences
    HashMap<String,HashSet<String>> acceptableNextResType = new HashMap<>();
    
    ArrayList<ArrayList<String>> rcAATypes = new ArrayList<>();//what res type each RC is (indexed by pos, rc)
    
    
    public GMECMutSpace(String mutFileName, ConfSpace confSpace){
        
        //First read the mut file, prepare acceptableNextResType
        readFile(mutFileName);
        
        //Now record the aa types based on confSpace
        for(int pos=0; pos<confSpace.numPos; pos++){
            ArrayList<String> aaTypesAtPos = new ArrayList<>();
            ArrayList<RC> RCs = confSpace.posFlex.get(pos).RCs;
            for(int rc=0; rc<RCs.size(); rc++){
                aaTypesAtPos.add( RCs.get(rc).AAType );
            }
            rcAATypes.add(aaTypesAtPos);
        }
    }
    
    
    public GMECMutSpace(String mutFileName, SimpleConfSpace confSpace){
        //First read the mut file, prepare acceptableNextResType
        readFile(mutFileName);
        
        //Now record the aa types based on confSpace
        for(int pos=0; pos<confSpace.getNumPos(); pos++){
            ArrayList<String> aaTypesAtPos = new ArrayList<>();
            List<ResidueConf> RCs = confSpace.positions.get(pos).resConfs;//posFlex.get(pos).RCs;
            for(int rc=0; rc<RCs.size(); rc++){
                aaTypesAtPos.add( RCs.get(rc).template.name );
            }
            rcAATypes.add(aaTypesAtPos);
        }
    }
    
    
    private void readFile(String mutFileName){
        try{
            FileInputStream is = new FileInputStream(mutFileName);
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
                        
            for(String curLine=bufread.readLine(); curLine!=null; curLine=bufread.readLine()){
                StringTokenizer st = new StringTokenizer(curLine);
                int numPos = st.countTokens();
                String seq = "";
                for(int pos=0; pos<numPos; pos++)//take only the 3-letter AA type
                    seq = seq + st.nextToken().substring(0,3).toUpperCase();
                
                for(int pos=0; pos<numPos; pos++){
                    String partialSeq = seq.substring(0, 3*pos);
                    if(!acceptableNextResType.containsKey(partialSeq))
                        acceptableNextResType.put(partialSeq, new HashSet<>());
                    
                    acceptableNextResType.get(partialSeq).add(seq.substring(3*pos, 3*(pos+1)));
                }
            }
            
            bufread.close();
        }
        catch(FileNotFoundException e){
            throw new RuntimeException("ERROR: Couldn't find mut file "+mutFileName);
        }
        catch(Exception e){
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }
    }
    
    public boolean isNewRCAllowed(int[] confSoFar, int numDefined, int newRC){
        //we've defined the first numDefined RC's in confSoFar;
        //if we add newRC then will this stay within our sequence space?
        String partialSeq = writePartialSeq(confSoFar,numDefined);
        
        if(!acceptableNextResType.containsKey(partialSeq))
            throw new RuntimeException("ERROR: partial conf has forbidden sequence");
        
        String newResType = rcAATypes.get(numDefined).get(newRC);
        return acceptableNextResType.get(partialSeq).contains(newResType);
    }
    
    String writePartialSeq(int[] confSoFar, int numDefined){
        String ans = "";
        for(int pos=0; pos<numDefined; pos++)
            ans += rcAATypes.get(pos).get(confSoFar[pos]);
        return ans;
    }
    
    
}
