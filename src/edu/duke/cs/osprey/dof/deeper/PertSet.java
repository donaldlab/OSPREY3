/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper;

import edu.duke.cs.osprey.dof.deeper.perts.PartialStructureSwitch;
import edu.duke.cs.osprey.dof.deeper.perts.LoopClosureAdjustment;
import edu.duke.cs.osprey.dof.deeper.perts.Shear;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.dof.deeper.perts.Backrub;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.dof.deeper.perts.PerturbationBlock;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.StringParsing;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 * 
 * Description of the perturbations in our system,
 * not set up in a molecule (describes the content of a .pert file)
 * 
 * @author mhall44
 */
public class PertSet implements Serializable {
    
    ArrayList<String> pertTypes = new ArrayList<>();//List of types ("BACKRUB", etc.) for the perturbations
    ArrayList<ArrayList<String>> resNums = new ArrayList<>();//corresponding perturbations'
    //lists of directly-affected residue numbers (PDB numbering)
    
    ArrayList<ArrayList<double[]>> pertIntervals = new ArrayList<>();
    //for each of the perturbations, a list of intervals for its parameter to be in
    //these will be used to define RCs
    
    ArrayList<ArrayList<ArrayList<int[]>>> pertStates = new ArrayList<>();
    //for each flexible residue position, for each "perturbation state" of the residue,
    //a list of pairs (pert#, interval#) defining what perturbations (indexed in pertNames)
    //have what intervals (defined in pertIntervals)
    
    
    //Some perturbations (now, just partial structure switches) 
    //need additional info to define them
    //For each perturbation, the lines of additional info are listed here (null if none)
    ArrayList<ArrayList<String>> additionalInfo = new ArrayList<>();
    
    
    public boolean loadPertFile(String pertFileName, boolean loadStates, ResidueTermini termini){
        //load perturbations from the pert file
        //Return whether we found the file or not
        //Load perturbations and their intervals; if loadStates then residue pert states too
        try{

            BufferedReader br=new BufferedReader(new FileReader(pertFileName));
            StringTokenizer st;
            br.readLine();//Title
            readPerts(br, termini);
            
            if(loadStates){
            
                pertStates = new ArrayList<>();


                while(br.readLine() != null){//Read residue perturbation states.  Skipping the "RES" line

                    //Removing residue number line (from OSPREY 2.x)
                    
                    int numStates = Integer.valueOf( StringParsing.getToken(br.readLine(), 1) );

                    ArrayList<Integer> resPerts = new ArrayList<>();
                    //indices of perturbations that can change the conf of this residue
                    //(to be listed presently)

                    //Read the perturbations affecting this residue
                    //This section is included in the file because it specifies if some indirect effects should be neglected
                    st = new StringTokenizer(br.readLine()," ");
                    int resNumPerts = st.countTokens() - 1;
                    st.nextToken();//"PERTURBATIONS"
                    for(int a=0; a<resNumPerts; a++)
                        resPerts.add( Integer.valueOf(st.nextToken()) );

                    ArrayList<ArrayList<int[]>> resPertStates = new ArrayList<>();

                    for(int a=0; a<numStates; a++){//Read perturbation states
                        ArrayList<int[]> pertState = new ArrayList<>();

                        st = new StringTokenizer(br.readLine()," ");
                        for(int b=0; b<resNumPerts; b++){
                            int[] p = new int[] {resPerts.get(b),Integer.valueOf(st.nextToken())};//pert #, interval #
                            pertState.add(p);
                        }
                        resPertStates.add(pertState);
                    }

                    pertStates.add(resPertStates);
                    //The RCs section in OSPREY 2.x is now omitted
                }
            }
            

            br.close();
        }
        catch(FileNotFoundException e){//basically means we should select perturbations
            return false;
        }
        catch(Exception e){//this is actually an error
            e.printStackTrace();
            throw new RuntimeException("ERROR READING PERTURBATION FILE: "+e.getMessage());
        }
        
        return true;
    }
    
    
    public void readPerts(BufferedReader br) throws Exception {
        //read the actual perturbations, including the residues they affect
        //and the parameter intervals we're using for them

        int numPerts = Integer.valueOf(br.readLine().trim());

        StringTokenizer st;
        
        pertTypes = new ArrayList<>();
        resNums = new ArrayList<>();
        pertIntervals = new ArrayList<>();
        additionalInfo = new ArrayList<>();
        

        for(int a=0;a<numPerts;a++){//Read perturbations
            String pertType = br.readLine();
            pertTypes.add(pertType);
            
            st = new StringTokenizer(br.readLine()," ");
            int numAffectedRes = st.countTokens();
            ArrayList<String> pertResNums = new ArrayList<>();

            for(int b=0;b<numAffectedRes;b++){
                String inputNumber = st.nextToken();
                pertResNums.add(inputNumber);
            }
            
            resNums.add(pertResNums);

            recordAdditionalInfo(br, pertType);
            
            st = new StringTokenizer(br.readLine()," ");
            int numStates = Integer.valueOf(st.nextToken());
            ArrayList<double[]> curPertIntervals = new ArrayList<>();

            for(int b=0;b<numStates;b++){
                st = new StringTokenizer(br.readLine()," ");
                if(st.countTokens() != 2)
                    throw new java.lang.Exception("Bad formatting of perturbation "+a);
                //read the interval
                double lo = Double.valueOf(st.nextToken());
                double hi = Double.valueOf(st.nextToken());
                curPertIntervals.add( new double[] {lo,hi} );
            }
            
            pertIntervals.add(curPertIntervals);
        }
    }
    
    
    public void readPerts(BufferedReader br, ResidueTermini termini) throws Exception {
        //read the actual perturbations, including the residues they affect
        //and the parameter intervals we're using for them

        int numPerts = Integer.valueOf(br.readLine().trim());

        StringTokenizer st;
        
        pertTypes = new ArrayList<>();
        resNums = new ArrayList<>();
        pertIntervals = new ArrayList<>();
        additionalInfo = new ArrayList<>();
        

        for(int a=0;a<numPerts;a++){//Read perturbations
            String pertType = br.readLine();
            pertTypes.add(pertType);
            
            st = new StringTokenizer(br.readLine()," ");
            int numAffectedRes = st.countTokens();
            ArrayList<String> pertResNums = new ArrayList<>();

            for(int b=0;b<numAffectedRes;b++){
                String inputNumber = st.nextToken();
                if(termini == null || termini.contains(inputNumber)) {
                	pertResNums.add(inputNumber);
                }
            }
            
            if(pertResNums.isEmpty()) {
            	pertTypes.remove(pertTypes.size()-1);
            	continue;
            }
            
            resNums.add(pertResNums);

            recordAdditionalInfo(br, pertType);
            
            st = new StringTokenizer(br.readLine()," ");
            int numStates = Integer.valueOf(st.nextToken());
            ArrayList<double[]> curPertIntervals = new ArrayList<>();

            for(int b=0;b<numStates;b++){
                st = new StringTokenizer(br.readLine()," ");
                if(st.countTokens() != 2)
                    throw new java.lang.Exception("Bad formatting of perturbation "+a);
                //read the interval
                double lo = Double.valueOf(st.nextToken());
                double hi = Double.valueOf(st.nextToken());
                curPertIntervals.add( new double[] {lo,hi} );
            }
            
            pertIntervals.add(curPertIntervals);
        }
        
        // pertfile does not apply to this strand. advance the pertfile to end
        if(pertTypes.isEmpty()) while(br.readLine() != null);
    }
    
    
    void recordAdditionalInfo(BufferedReader br, String pertType) throws IOException {
        
        if(pertType.equalsIgnoreCase("PARTIAL STRUCTURE SWITCH")){
            //record the PDB file names
            //Additional info should be like this:
            //2 structures
            //1ABC.pdb
            //(JUST ONE PDB LISTED, BECAUSE THE FIRST STRUCTURE IS ALWAYS THE ORIGINAL ONE)
            int numStructs = Integer.valueOf( StringParsing.getToken( br.readLine(), 1) );
            
            ArrayList<String> altPDBs = new ArrayList<>();
            
            for(int structNum=1; structNum<numStructs; structNum++)
                altPDBs.add(br.readLine().trim());
            
            additionalInfo.add(altPDBs);
        }
        else//no additional info needed
            additionalInfo.add(null);
    }
    
    
    public void writePertFile(String pertFileName){
        try{
            BufferedWriter bw=new BufferedWriter(new FileWriter(pertFileName));
            bw.append("PERTURBATIONS");
            bw.newLine();
            
            int numPerts = pertTypes.size();
            bw.append(String.valueOf(numPerts));
            bw.newLine();

            for(int pertNum=0; pertNum<numPerts; pertNum++){//Perturbation info
                bw.append(pertTypes.get(pertNum));
                //Perturbation pert = m.perts[pertNum];
                //bw.append(pert.type);

                bw.newLine();
                for(String resNum : resNums.get(pertNum))
                    bw.append(resNum+" ");
                

                bw.newLine();

                writeAdditionalInfo(bw, additionalInfo.get(pertNum));
                
                int numIntervals = pertIntervals.get(pertNum).size();
                bw.append(numIntervals+" states");
                bw.newLine();
                for(int state=0; state<numIntervals; state++){
                    double interv[] = pertIntervals.get(pertNum).get(state);
                    bw.append(interv[0]+" "+interv[1]);
                    bw.newLine();
                }
            }
            
            
            //OK now go through the residues to write their perturbation states
            for(int pos=0;pos<pertStates.size(); pos++){//Residue perturbation state, RC info
                //Residue res=m.residue[pos];
                //int posInStrand = res.strandResidueNumber;
                
                ArrayList<ArrayList<int[]>> resPertStates = pertStates.get(pos);
                //if( res.perts.length > 0 ){
                //UNLIKE OSPREY 2.X WILL HAVE RES PERT STATES RECORD FOR EVERY RESIDUE
                    //curStrandRCs = (StrandRCs)strandRot[res.strandNumber];

                    bw.append("RES");
                    bw.newLine();

                    bw.append(resPertStates.size() + " states ");
                    bw.newLine();

                    bw.append("PERTURBATIONS ");
                    if(resPertStates.size()>0){
                        ArrayList<int[]> firstState = resPertStates.get(0);
                        for(int a=0;a<firstState.size();a++)
                            bw.append(firstState.get(a)[0] + " ");//firstState consists of pairs (pert #, interval #)
                    }
                    bw.newLine();

                    for(int state=0; state<resPertStates.size(); state++){
                        ArrayList<int[]> pertState = resPertStates.get(state);
                        
                        for(int pertInd=0; pertInd<pertState.size(); pertInd++)
                            bw.append(pertState.get(pertInd)[1] + " ");

                        bw.newLine();
                    }

                //}
            }

            bw.close();
        }

        catch(IOException e){
            e.printStackTrace();
            throw new RuntimeException("ERROR WRITING PERTURBATION FILE: "+e.getMessage());
        }
    }
    
    
    private void writeAdditionalInfo(BufferedWriter bw, ArrayList<String> altPDBs) throws IOException {
        //write additional info needed to reconstruct a perturbation
        if(altPDBs!=null){//there is some info needed (only for partial structure switch currently)
            int numStructs = altPDBs.size()+1;//include original structure
            bw.append(numStructs+" structures");
            bw.newLine();
            
            for(String fileName : altPDBs){
                bw.append(fileName);
                bw.newLine();
            }
        }
    }
    
    
    ArrayList<Perturbation> makePerturbations(Molecule m){
        //Generate these perturbations in the molecule of interest
        //and build a block out of them
        
        ArrayList<Perturbation> ans = new ArrayList<>();
        
        for(int pertNum=0; pertNum<pertTypes.size(); pertNum++){
            
            String type = pertTypes.get(pertNum);
            ArrayList<Residue> directlyAffectedResidues = new ArrayList<>();
            for(String resNum : resNums.get(pertNum))
                directlyAffectedResidues.add( m.getResByPDBResNumber(resNum) );
            
            Perturbation pert;
            
            if(type.equalsIgnoreCase("BACKRUB"))
                pert = new Backrub(directlyAffectedResidues);
            else if(type.equalsIgnoreCase("SHEAR"))
                pert = new Shear(directlyAffectedResidues);
            else if(type.equalsIgnoreCase("LOOP CLOSURE ADJUSTMENT"))
                pert = new LoopClosureAdjustment(directlyAffectedResidues);
            else if(type.equalsIgnoreCase("PARTIAL STRUCTURE SWITCH"))
                pert = new PartialStructureSwitch(directlyAffectedResidues, additionalInfo.get(pertNum));
            else
                throw new RuntimeException("ERROR: Unrecognized perturbation type: "+type);
            
            ans.add(pert);
        }
        
        PerturbationBlock pblock = new PerturbationBlock(ans);
        //pblock will automatically be stored in each of the perturbations
        
        return ans;
    }
    
    
    
    
    PertSet makeDiscreteVersion(){
        //Make a discrete version of this PertSet.  Shallow copying when no change
        PertSet discrSet = new PertSet();
        discrSet.pertTypes = pertTypes;
        discrSet.resNums = resNums;
        discrSet.pertStates = pertStates;
        discrSet.additionalInfo = additionalInfo;
        
        for(ArrayList<double[]> curPertIntervals : pertIntervals){
            ArrayList<double[]> curDiscr = new ArrayList<>();
            for(double[] bounds : curPertIntervals){
                if(bounds[0]==bounds[1])
                    curDiscr.add(bounds);
                else{
                    double discrValue = (bounds[0]+bounds[1])/2;
                    curDiscr.add(new double[] {discrValue,discrValue});
                }
            }
            discrSet.pertIntervals.add(curDiscr);
        }
        
        return discrSet;
    }
    
    
}
