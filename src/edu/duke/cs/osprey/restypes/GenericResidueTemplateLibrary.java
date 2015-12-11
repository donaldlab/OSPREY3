/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.StringParsing;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 *
 * @author mhall44
 */
public class GenericResidueTemplateLibrary extends ResidueTemplateLibrary {
    //This library of residue templates defines what types of residues we can model
    //and what flexibility and energy parameters come with each type
    //NAMING: We assume each distinct residue (AA or otherwise) has its own name
    //however, many residues will have multiple slightly different forms (N-terminal, etc.) 
    //and these will share a name and a rotamer library entry
    
    

    
    public int totalNumRotamers = 0;//total number of rotamers read in from rotamer library file(s)
    //starts at 0
    
    HashMap<String,Double> resEntropy = new HashMap<>();//We will look up residue entropy
    //by the name of the residue
    
    //the templates contain forcefield information like atom types, so they need to go with 
    //a set of forcefield parameters
    public ForcefieldParams ffParams;
    
    public GenericResidueTemplateLibrary(String[] templateFiles, ForcefieldParams fp){
        //create a library based on template files
        //we can then load coordinates and rotamer libraries for these templates separately, if we have these

        ffParams = fp;
        
        for(String fileName : templateFiles){
            loadTemplates(fileName);
        }
    }
        
        
    public void loadTemplates(String templateFile){//load templates from the indicated files,
        //and add them to our list of templates
        try {

            FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir().concat(templateFile) );
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));

            // Skip over first 2 lines of header info
            bufread.readLine();
            bufread.readLine();

            while(true) {//read all the templates
                //at the beginning of this loop, curLine is the long amino acid name,
                //which we discard.  readTemplate will thus start reading the next line
                ResidueTemplate newTemplate = readTemplate(bufread);
                if(newTemplate==null)//null newTemplate means lines needed to be skipped or something
                    break;
                else {
                    templates.add(newTemplate);
                }
            }
            bufread.close();
        }
        
        catch (FileNotFoundException e) {
                System.out.println("ERROR: Template File Not Found: "+e);
                e.printStackTrace();
                System.exit(0);
        }
        catch (IOException e) {
                System.out.println("ERROR reading template file: "+e);
                e.printStackTrace();
                System.exit(0);
        }
    }
    

    
    private ResidueTemplate readTemplate (BufferedReader bufread) throws IOException {
        //read a template from the BufferedReader provided (it reads a template file)
        //null means all templates already read from file
        
        String curLine = bufread.readLine();
        if(curLine==null)//file ended!
            return null;
        else if (curLine.length() >= 4){
                if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
                        curLine = bufread.readLine();
                        return null;//finished reading file!
                }
        }
        // Skip blank line
        curLine = bufread.readLine();
        // The next line contains the 3 letter amino acid name
        curLine = bufread.readLine();
        String templateName = StringParsing.getToken(curLine,1);
        
        // Skip next 2 lines
        curLine = bufread.readLine();
        curLine = bufread.readLine();
        // Now we're into the section with atoms
        curLine = bufread.readLine();
        // Skip the dummy atoms
        int dumPresent = 0;
        while (StringParsing.getToken(curLine,2).equalsIgnoreCase("DUMM")) {
                dumPresent++;
                curLine = bufread.readLine();
        }
        dumPresent++; // to adjust for 0-based
        
        ArrayList<Atom> atomList = new ArrayList<>();
        
        while (!StringParsing.getToken(curLine,2).equals("")) {//read info on atoms
                String atomName = StringParsing.getToken(curLine,2);
                Atom at = new Atom (atomName);
                
                at.forceFieldType = StringParsing.getToken(curLine,3);
                at.type = ffParams.atomTypeToInt(at.forceFieldType);
                at.charge = (double) (new Double(StringParsing.getToken(curLine,11)).doubleValue());

                //the template only records bonds within the residue
                
                //KER: The first atom is bonded to a dummy atom so we can't include that
                //KER: in the bond list, so check atom is >= 0
                int atomBondedTo = ((new Integer(StringParsing.getToken(curLine,5))).intValue())-dumPresent;
                if(atomBondedTo >=0){
                        at.addBond(atomList.get(atomBondedTo));
                }
                
                atomList.add(at);
                curLine = bufread.readLine();  // read next line
        }

        
        Residue templateRes = new Residue(atomList,null,templateName,null);//no molecule or coordinates yets


        do {//we expect one or more blank lines before the LOOP and IMPROPER records
            curLine = bufread.readLine();
        }
        while(curLine.trim().isEmpty());
        
        //KER: Read LOOP data if any
        if (curLine.length() >= 4){
                if(StringParsing.getToken(curLine, 1).equalsIgnoreCase("LOOP")){
                        curLine = bufread.readLine();
                        while(!StringParsing.getToken(curLine,2).equals("")){
                                //find atom1
                                for(Atom a : atomList){
                                        if(a.name.equalsIgnoreCase(StringParsing.getToken(curLine,1))){
                                                //find atom2
                                                for(Atom b : atomList){
                                                        if(b.name.equalsIgnoreCase(StringParsing.getToken(curLine,2))){
                                                                a.addBond(b);
                                                        }
                                                }
                                        }
                                }
                                curLine = bufread.readLine();
                        }
                }
        }
        
        //at this point templateRes has all its intra-res bonds all marked
        templateRes.intraResBondsMarked = true;

        // Eventually we might want to be able to handle the improper
        //  torsions listed here

        
        // Read until the end of the residue
        boolean atDone = false;
        if (curLine.length() >= 4)
                atDone = curLine.substring(0,4).equalsIgnoreCase("done");
        else
                atDone = false;
        while (!atDone) {
                curLine = bufread.readLine();
                if (curLine.length() >= 4)
                        atDone = curLine.substring(0,4).equalsIgnoreCase("done");
        }
        
        
        ResidueTemplate newTemplate = new ResidueTemplate(templateRes,templateName);
        return newTemplate;
    }
    
    
    
    public void loadTemplateCoords(String fileName){
        //load coordinates for templates, given in the following file
        //they should match templates that are already loaded (same residue and atom names)
        //we will need coordinates for any residue type that we mutate to
        
        try {
                
            FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir().concat(fileName) );
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
            String curLine = null, tmpName = null;
            int tmpCtr = 0;

            curLine = bufread.readLine();
            while ( curLine.startsWith("#") ){
                    curLine = bufread.readLine();
            }
            boolean foundRes = false;
            boolean foundAtom = false;
            while (curLine != null ){
                    String resName = StringParsing.getToken(curLine,1);
                    int numAtoms = new Integer(StringParsing.getToken(curLine,2)); 
                    foundRes = false;
                    for(ResidueTemplate template : templates){//find the template to assign these coordinates to
            
                        if(template.name.equalsIgnoreCase(resName)){//found it
                            
                                Residue r = template.templateRes;
                                if(r.atoms.size()!=numAtoms)
                                    throw new RuntimeException("ERROR: Coords file has wrong number of atoms for "+r.fullName);

                                r.coords = new double[3*numAtoms];//allocate coordinates

                                foundRes = true;
                                for(int i=0;i<numAtoms;i++){
                                        curLine = bufread.readLine();
                                        //Find the current atom in the residue
                                        foundAtom = false;
                                        for(int atNum=0; atNum<numAtoms; atNum++){
                                                Atom at = r.atoms.get(atNum);
                                                if(at.name.equalsIgnoreCase(StringParsing.getToken(curLine,1))){
                                                        foundAtom = true;
                                                        //read coordinates for this atom
                                                        double x = new Double(StringParsing.getToken(curLine,2));
                                                        double y = new Double(StringParsing.getToken(curLine,3)); 
                                                        double z = new Double(StringParsing.getToken(curLine,4));
                                                        r.coords[3*atNum] = x;
                                                        r.coords[3*atNum+1] = y;
                                                        r.coords[3*atNum+2] = z;
                                                        break;
                                                }
                                        }
                                        if(!foundAtom){
                                                System.out.println("Residue coord and forcefield templates did not match up.");
                                                System.exit(0);
                                                //Possible issue: first occurrence of this name in the template residues is the wrong form for these coordinates?
                                        }		
                                }
                                break;
                        }
                    }
                    //If we didn't find a match we need to get rid of those 
                    //lines from the file
                    if(!foundRes){
                            for(int i=0; i<numAtoms;i++){
                                    curLine=bufread.readLine();
                            }
                    }
                    //Read to catch the ENDRES line and then
                    //get the start of the next AA
                    curLine = bufread.readLine();
                    curLine = bufread.readLine();
            }
            bufread.close();
        }
        catch (FileNotFoundException e) {
                System.out.println("ERROR: Template File Not Found: "+e);
                e.printStackTrace();
                System.exit(0);
        }
        catch (IOException e) {
                System.out.println("ERROR reading template file: "+e);
                e.printStackTrace();
                System.exit(0);
        }
        
    }
    
    /** 
     * PGC 2015: Supporting two rotamer libraries right now, Lovell and Dunbrack backbone dependent rotamers.
     * load rotamer information into templates.
     * @param filename for Lovell-style or Dunbrack Rotamer Library.
     * @param backbone_dependent_rotamers Use Dunbrack Rotamer Library? 
     */
    public void loadRotamerLibrary(String fileName, boolean dunbrack_backbone_dependent_rotamers){

        if(dunbrack_backbone_dependent_rotamers){
        	RotamerLibraryReader.readDunbrackRotamerLibraryForResiduePosition(fileName, this);
        }
        else{
        	RotamerLibraryReader.readRotLibrary(fileName, this);
        }
        //read volume here too?
        
    }
    
    
    
    public ResidueTemplate getTemplateForMutation(String resTypeName, Residue res, boolean errorIfNone){
        //We want to mutate res to type resTypeName.  Get the appropriate template.
        //Currently only one template capable of being mutated to (i.e., having coordinates)
        //is available for each residue type.  If this changes update here!
        for(ResidueTemplate template : templates){
            if(template.name.equalsIgnoreCase(resTypeName)){
                if(template.templateRes.coords!=null){
                    //we have coordinates for templateRes, so can mutate to it
                    return template;
                }
            }
        }
        
        if(errorIfNone){//actually trying to mutate...throw an error if can't get a mutation
            throw new RuntimeException("ERROR: Couldn't find a template for mutating "+res.fullName
                    +" to "+resTypeName);
        }
        else//just checking if template available for mutation...return null to indicate not possible
            return null;
    }
    
    
    public void makeDAminoAcidTemplates(){
        //Make a D-amino acid template for every standard L-amino acid template in the library
        //This lets us do D-amino acids without having to have separate templates
        //A D-amino acid is just an inverted L-amino acid (including Ile and Thr,
        //which have a second chiral center)
        ArrayList<ResidueTemplate> DTemplates = new ArrayList<>();
        
        for(ResidueTemplate template : templates){
            
            if( DAminoAcidHandler.isStandardLAminoAcid(template.name) ){
                //template recognized as standard L-amino acid
                
                DTemplates.add( DAminoAcidHandler.makeDTemplate(template) );
            }
        }
        
        templates.addAll(DTemplates);
    }


    @Override
    public int numRotForResType(int pos, String resType, double phi, double psi) {
        return firstTemplate(resType).getNumRotamers(phi, psi);
    }
    
    
    public void loadResEntropy(String entropyFile){
        //It is convenient to load residue entropies into a hash map, rather than
        //into template objects, because they correspond to template names
        try {
            FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir().concat(entropyFile) );
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
            
            String curLine = bufread.readLine();
            
            while (curLine != null ){
                String resType = StringParsing.getToken(curLine,1);
                double entropy = new Double(StringParsing.getToken(curLine,2)); 
                resEntropy.put(resType.toUpperCase(), entropy);
                curLine = bufread.readLine();
            }
            bufread.close();
        }
        catch (FileNotFoundException e) {
                System.out.println("ERROR: Residue entropy file not found: "+e);
                e.printStackTrace();
                System.exit(0);
        }
        catch (IOException e) {
                System.out.println("ERROR reading residue entropy file: "+e);
                e.printStackTrace();
                System.exit(0);
        }
    }
    
    
    public double getResEntropy(String resType){
        if(resEntropy.containsKey(resType.toUpperCase()))
            return resEntropy.get(resType.toUpperCase());
        else//default
            return 0;
    }


    
}

