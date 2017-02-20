/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import edu.duke.cs.osprey.control.Defaults;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.Forcefield;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 * This library of residue templates defines what types of residues Osprey can model
 * and what flexibility and energy parameters come with each type.
 *
 * @author mhall44
 */
public class GenericResidueTemplateLibrary extends ResidueTemplateLibrary implements Serializable {
	
    //NAMING: We assume each distinct residue (AA or otherwise) has its own name
    //however, many residues will have multiple slightly different forms (N-terminal, etc.) 
    //and these will share a name and a rotamer library entry
	
	private static final long serialVersionUID = 2643396106081986558L;

	public static class Builder {
		
		private static final String LovellRotamersPath = "/config/LovellRotamer.dat";
		
		/** used to match molecule residues to templates */
		private Forcefield forcefield;
		
		/** text of file with template coordinates */
		private String templateCoordsText;
		
		/** text of file with rotamer dihedral angles */
		private String rotamersText;
		
		/** text of file with backbone dependent rotamer dihedral angles */
		private String backboneDependentRotamersText;
		private String entropyText;
		private boolean makeDAminoAcidTemplates;
		
		public Builder() {
			forcefield = Defaults.forcefield;
			templateCoordsText = FileTools.readResource("/config/all_amino_coords.in");
			rotamersText = FileTools.readResource(LovellRotamersPath);
			backboneDependentRotamersText = null;
			entropyText = FileTools.readResource("/config/ResEntropy.dat");
			makeDAminoAcidTemplates = true;
		}
		
		public Builder setForcefield(Forcefield val) {
			forcefield = val;
			return this;
		}
		
		public Builder setTemplateCoords(String text) {
			templateCoordsText = text;
			return this;
		}
		
		public Builder setRotamers(String text) {
			rotamersText = text;
			return this;
		}
		
		public Builder setLovellRotamers() {
			return setRotamers(FileTools.readResource(LovellRotamersPath));
		}
		
		public Builder setBackboneDependentRotamers(String text) {
			backboneDependentRotamersText = text;
			return this;
		}
		
		public Builder setEntropy(String text) {
			entropyText = text;
			return this;
		}
		
		public Builder setMakeDAminoAcidTemplates(boolean val) {
			makeDAminoAcidTemplates = val;
			return this;
		}
		
		public GenericResidueTemplateLibrary build() {
			GenericResidueTemplateLibrary lib = new GenericResidueTemplateLibrary(forcefield);
			lib.loadTemplateCoords(templateCoordsText);
			if (rotamersText != null) {
				lib.loadRotamerLibrary(rotamersText);
			}
			if (backboneDependentRotamersText != null) {
				lib.loadBackboneDependentRotamerLibrary(backboneDependentRotamersText);
			}
			if (makeDAminoAcidTemplates) {
				lib.makeDAminoAcidTemplates();
			}
			if (entropyText != null) {
				lib.loadResEntropy(entropyText);
			}
			return lib;
		}
	}
	
	public static Builder builder() {
		return new Builder();
	}
    
    public int totalNumRotamers = 0;//total number of rotamers read in from rotamer library file(s)
    //starts at 0
    
    HashMap<String,Double> resEntropy = new HashMap<>();//We will look up residue entropy
    //by the name of the residue
    
    //the templates contain forcefield information like atom types, so they need to go with 
    //a set of forcefield parameters
    public final ForcefieldParams ffParams;
    
    public GenericResidueTemplateLibrary(Forcefield forcefield) {
        //create a library based on template files
        //we can then load coordinates and rotamer libraries for these templates separately, if we have these

        ffParams = new ForcefieldParams(forcefield);
        
        loadTemplates(FileTools.readResource(ffParams.forcefld.aaPath));
        loadTemplates(FileTools.readResource(ffParams.forcefld.aaNTPath));
        loadTemplates(FileTools.readResource(ffParams.forcefld.aaCTPath));
        loadTemplates(FileTools.readResource(ffParams.forcefld.grPath));
    }
    
    public void loadTemplates(String text) {
        
        Iterator<String> lines = FileTools.parseLines(text).iterator();
        
        // Skip over first 2 lines of header info
        lines.next();
        lines.next();

        while (true) {//read all the templates
            //at the beginning of this loop, curLine is the long amino acid name,
            //which we discard.  readTemplate will thus start reading the next line
            ResidueTemplate newTemplate = readTemplate(lines);
            if (newTemplate==null) {//null newTemplate means lines needed to be skipped or something
                break;
            } else {
                templates.add(newTemplate);
            }
        }
    }
    

    
    private ResidueTemplate readTemplate(Iterator<String> lines) {
        //read a template from the BufferedReader provided (it reads a template file)
        //null means all templates already read from file
        
        String curLine = lines.next();
        if(curLine==null)//file ended!
            return null;
        else if (curLine.length() >= 4){
                if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
                        curLine = lines.next();
                        return null;//finished reading file!
                }
        }
        // Skip blank line
        curLine = lines.next();
        // The next line contains the 3 letter amino acid name
        curLine = lines.next();
        String templateName = StringParsing.getToken(curLine,1);
        
        // Skip next 2 lines
        curLine = lines.next();
        curLine = lines.next();
        // Now we're into the section with atoms
        curLine = lines.next();
        // Skip the dummy atoms
        int dumPresent = 0;
        while (StringParsing.getToken(curLine,2).equalsIgnoreCase("DUMM")) {
                dumPresent++;
                curLine = lines.next();
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
                curLine = lines.next();
        }

        
        Residue templateRes = new Residue(atomList,(double[])null,templateName,null);//no molecule or coordinates yets


        do {//we expect one or more blank lines before the LOOP and IMPROPER records
            curLine = lines.next();
        }
        while(curLine.trim().isEmpty());
        
        //KER: Read LOOP data if any
        if (curLine.length() >= 4){
                if(StringParsing.getToken(curLine, 1).equalsIgnoreCase("LOOP")){
                        curLine = lines.next();
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
                                curLine = lines.next();
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
                curLine = lines.next();
                if (curLine.length() >= 4)
                        atDone = curLine.substring(0,4).equalsIgnoreCase("done");
        }
        
        
        ResidueTemplate newTemplate = new ResidueTemplate(templateRes,templateName);
        return newTemplate;
    }
    
    
    
    public void loadTemplateCoords(String text) {
             
        Iterator<String> lines = FileTools.parseLines(text).iterator();
        
        String curLine = null;

            curLine = lines.next();
            while ( curLine.startsWith("#") ){
                    curLine = lines.next();
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
                                        curLine = lines.next();
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
                                                throw new Error("Residue coord and forcefield templates did not match up.");
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
                                    curLine = lines.next();
                            }
                    }
                    //Read to catch the ENDRES line and then
                    //get the start of the next AA
                    curLine = lines.next();
                    curLine = lines.next();
            }
    }
    
    public void loadRotamerLibrary(String text) {
        RotamerLibraryReader.readRotLibrary(text, this);
    }
    
    public void loadBackboneDependentRotamerLibrary(String text) {
        RotamerLibraryReader.readDunbrackRotamerLibraryForResiduePosition(text, this);
    }
    
    
    
    public ResidueTemplate getTemplate(String resTypeName) {
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
        return null;
    }
    
    public ResidueTemplate getTemplateForMutation(String resTypeName, Residue res) {
        //We want to mutate res to type resTypeName.  Get the appropriate template.
        ResidueTemplate template = getTemplate(resTypeName);
        if (template != null) {
            return template;
        }
        
        //actually trying to mutate...throw an error if can't get a mutation
        throw new RuntimeException("ERROR: Couldn't find a template for mutating " + res.fullName + " to " + resTypeName);
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
    
    
    public void loadResEntropy(String text){
        //It is convenient to load residue entropies into a hash map, rather than
        //into template objects, because they correspond to template names
        for (String line : FileTools.parseLines(text)) {
            
            // skip comments
            if (line.startsWith("%")) {
                continue;
            }
            
            String resType = StringParsing.getToken(line,1);
            double entropy = new Double(StringParsing.getToken(line,2)); 
            resEntropy.put(resType.toUpperCase(), entropy);
        }
    }
    
    
    public double getResEntropy(String resType){
        if(resEntropy.containsKey(resType.toUpperCase()))
            return resEntropy.get(resType.toUpperCase());
        else//default
            return 0;
    }


    
}

