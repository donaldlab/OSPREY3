/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

import java.io.*;
import java.util.*;

/**
 * This class implements a rotamer library reader. It reads from an input file that contains rotamer information
 * for amino acid types or other generic residues.
 */
public class RotamerLibraryReader implements Serializable {
	
	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = false;
	
	
	//Read in all of the rotamers for all amino acids from the rotFilename file
        //into the ResidueTemplateLibrary
	public static int readRotLibrary(String text, List<ResidueTemplate> templates) {
		
		//String volFilename = rotFile + ".vol";

		int numRotamersRead = 0;

		// HANDLE THE NORMAL AAs
		Iterator<String> lines = FileTools.parseLines(text).iterator();
		int numLines = 0;
		while (lines.hasNext()) {
			String line = lines.next();
			
			// Skip over comments (lines starting with !)
			if (line.startsWith("!")) {
				continue;
			}

			// skip over blank lines
			if (line.isEmpty()) {
				continue;
			}
			
			// skip over the first line (number of residues in library)
			if (numLines++ == 0) {
				continue;
			}
			
			String aaName = StringParsing.getToken(line,1);
			int numDihedrals = Integer.parseInt(StringParsing.getToken(line,2));
			int numRotamers = Integer.parseInt(StringParsing.getToken(line,3));
			if (numRotamers<0) {
				numRotamers = 0;
			}
		
			String dihedralAtomNames[][] = new String[numDihedrals][4];//atoms involved in each dihedral
			double rotamerValues[][] = new double[numRotamers][numDihedrals];//dihedral values at each rotamer

			// Read in the actual dihedrals
			for(int q=0;q<numDihedrals;q++) {
				line = lines.next();
				dihedralAtomNames[q][0] = StringParsing.getToken(line,1);
				dihedralAtomNames[q][1] = StringParsing.getToken(line,2);
				dihedralAtomNames[q][2] = StringParsing.getToken(line,3);
				dihedralAtomNames[q][3] = StringParsing.getToken(line,4);
			}
			// Read in the actual rotamers
			for(int q=0;q<numRotamers;q++) {
				line = lines.next();
				for(int w=0;w<numDihedrals;w++) {
					rotamerValues[q][w] = Integer.parseInt(StringParsing.getToken(line,(w+1)));
				}
			}
			
            //we've now got all the rotamer information for the residue type aaName
            //record it in all templates for aaName in the template library
            boolean foundTemplate = false;
            
            for(ResidueTemplate template : templates){
                if(template.name.equalsIgnoreCase(aaName)){
                    
                    //record information in template
                    template.numDihedrals = numDihedrals;
                    // Backbone independent rotamer libraries have a single "bin" for backbone dihedrals.
                    template.setNumberOfPhiPsiBins(1);
                    // And the resolution for a backbone independent rotamer library's bin is 360 (i.e. there is a single bin)
                    template.setRLphiPsiResolution(360);
                    template.initializeRotamerArrays();
                    
                    // Use 0 and 0 for the phi and psi angles, to indicate a backbone independent rotamer library.
                    template.setNumRotamers(numRotamers, 0, 0);
                    template.setRotamericDihedrals(rotamerValues, 0, 0);
                    
                    //convert atom names to indices for the atoms in templateRes
                    template.dihedral4Atoms = new int[dihedralAtomNames.length][4];
                    for(int dihedNum=0; dihedNum<numDihedrals; dihedNum++){
                        for(int a=0; a<4; a++){
                            String atomName = dihedralAtomNames[dihedNum][a];
                            template.dihedral4Atoms[dihedNum][a] = template.templateRes.getAtomIndexByName(atomName);
                        }
                    }
                    
                    //compute what atoms move when the dihedral is adjusted
                    template.computeDihedralMovingAtoms();
                    foundTemplate = true;
                }
            }
            
            if(!foundTemplate){
                throw new RuntimeException("ERROR: Have rotamer information for residue type "
                        +aaName+" but can't find a template for it");
            }
            
			//update total number of rotamers read into templateLib
			numRotamersRead = Math.max(1, numRotamers);
		}

		return numRotamersRead;
	}
	/**
	 * 
	 * @author pablo
	 * BBDepRotamersForTypePhiAndPsi: Backbone Dependent Rotamers For Amino Acid Type given a Phi and Psi angle value.
	 * This class is a data structure that stores all the rotamers for a specific phi and psi bin, given an angle.
	 */
	private static class BBDepRotamersForTypePhiAndPsi{
		ArrayList <double[]>rotamersForPhiAndPsi;
		public double phi;
		public double psi;
		public BBDepRotamersForTypePhiAndPsi(double phi, double psi){
			this.rotamersForPhiAndPsi = new ArrayList<double[]>();
			this.phi = phi;
			this.psi = psi;
		}
		/**
		 * Adds a new rotamer to list of rotamers for this type, phi and psi values.
		 * @param dihedralValues
		 */
		public void addRotamerForTypePhiAndPsi(double dihedralValues[]){
			rotamersForPhiAndPsi.add(dihedralValues);
		}
		/**
		 * Converts the array list rotamersForPhiAndPsi to a 2D array of type double
		 * @return a 2D array of rotamers and dihedrals
		 */
		public double [][] convertRotamersToArray(){
			double arrayOfRotamers[][] = new double[this.rotamersForPhiAndPsi.size()][];
			for (int rot = 0; rot < this.rotamersForPhiAndPsi.size(); rot++){
				arrayOfRotamers[rot] =  this.rotamersForPhiAndPsi.get(rot);
			}
			return arrayOfRotamers;
		}
	}
	/**
	 * @author PGC 2015
	 * 
	 * BBDepRotamersForType: Backbone dependent rotamers for type.
	 * This class is a data structure that stores all the backbone rotamers for a specific amino acid type.
	 * It is used by the Dunbrack parser.
	 */
	private static class BBDepRotamersForType{
		// aaType: amino acid type for which we are storing the rotamers.
		public String aaType;
		// resolution: The increments on which the rotamer library classifies rotamers. For Dunbrack this is 10.
		public double resolution;
		// numBins: The number of "bins" for each phi and psi. For Dunbrack this is 37 (e.g. at -180, -170, ..., 0, ..., 170, 180)
		public int numBins;
		// The array structure that stores the rotamers for each combination of phi and psi angles. The size of this 2D array is numBins*numBins.
		public BBDepRotamersForTypePhiAndPsi rotamersForType[][];
		
		/**
		 * 
		 * @param aaType amino acid type for which we are storing the rotamers.
		 * @param numBins The increments on which the rotamer library classifies rotamers. For Dunbrack this is 10.
		 * @param resolution The number of "bins" for each phi and psi. For Dunbrack this is 37 (e.g. at -180, -170, ..., 0, ..., 170, 180)
		 */
		public BBDepRotamersForType(String aaType, int numBins, double resolution){
			this.aaType = aaType;
			this.numBins = numBins;
			this.resolution = resolution;

			rotamersForType = new BBDepRotamersForTypePhiAndPsi[numBins][];
			for(int phiBin = 0; phiBin < numBins; phiBin++){
				rotamersForType[phiBin] = new BBDepRotamersForTypePhiAndPsi[numBins];
				// Each bin corresponds to a phi and psi value.  We can compute each using the resolution and numBins info.
				double phi = (phiBin - numBins/2)*resolution ;
				for(int psiBin = 0; psiBin < numBins; psiBin++){
					double psi = (psiBin- numBins/2)*resolution ;
					// Each field in rotamersForType stores a data structure with all the rotamers for that phi and psi bin.
					rotamersForType[phiBin][psiBin] = new BBDepRotamersForTypePhiAndPsi(phi, psi);				
				}				
			}					
		}
		/** 
		 * Add a new rotamer to its respective bin for this aaType.
		 * @param phi The phi angle for this rotamer
		 * @param psi The psi angle for this rotamer.
		 * @param dihedralValues
		 */
		public void addNewRotamer(double phi, double psi, double dihedralValues[]){
	    	int phiBin = (int)((Math.round(phi/resolution))) + this.numBins/2;    	
			int psiBin = (int)((Math.round(psi/resolution))) + this.numBins/2;
			rotamersForType[phiBin][psiBin].addRotamerForTypePhiAndPsi(dihedralValues);
		}

	} 
	/** PGC 2015:  Read all the rotamer entries for the dunbrack rotamer library and add the corresponding ones to the residue templates in 
	 *    the template library templateLib. 
	 */
	public static void readDunbrackRotamerLibraryForResiduePosition(String text, List<ResidueTemplate> templates){

		// Dunbrack has rotamers in increments of 10 for phi and psi.
		double dunbrackResolution = 10;
		// Number of bins in the dunbrack rotamer library.
		int binsDunbrack = 37;
		// Create a hash map to map every amino acid in the template library to its rotamers
		HashMap<String, BBDepRotamersForType> allDunbrackRotamersMap = new HashMap<String, BBDepRotamersForType>();
		for(ResidueTemplate template : templates){
			BBDepRotamersForType bbDepRotamerForTypeOfTemplate = new RotamerLibraryReader.BBDepRotamersForType(template.name, binsDunbrack, dunbrackResolution);
			allDunbrackRotamersMap.put(template.name, bbDepRotamerForTypeOfTemplate);
		
		}
		
		// Read in the entire dunbrack library
		// Parse the whole file first and add each rotamer entry to the array list
		for (String line : FileTools.parseLines(text)) {
			
			// Skip comment lines: those that start with #
			if (line.startsWith("#")) {
				continue;
			}
			
			// Convert the line into a char array.
			char dbLine [] = line.toCharArray(); 
			// First three letters are the amino acid.
			String aaNameNewRot = ""+dbLine[0]+dbLine[1]+dbLine[2];
			// We compute the number of dihedrals for this AA by counting the number of none zeroes in columns r1-r4.
			int r1 = Integer.parseInt((""+dbLine[24]+dbLine[25]).trim());
			int r2 = Integer.parseInt((""+dbLine[27]+dbLine[28]).trim());
			int r3 = Integer.parseInt((""+dbLine[30]+dbLine[31]).trim());
			int r4 = Integer.parseInt((""+dbLine[33]+dbLine[34]).trim());
			int numDih = 1;
			if(r2 > 0){
				numDih++;
			} 
			if(r3 > 0){
				numDih++;
			}
			if(r4 > 0){
				numDih++;
			}
			// Parse the backbone phi and psi for this rotamer.
			// Phi is in chars 5, 6, 7, and 8 
			String phiForThisEntryAsStr = ""+dbLine[5]+dbLine[6]+dbLine[7]+dbLine[8];
			double phiForThisEntry = Double.valueOf(phiForThisEntryAsStr.trim());
			// Psi is in chars 10, 11, 12, 13
			String psiForThisEntryAsStr = ""+dbLine[10]+dbLine[11]+dbLine[12]+dbLine[13];
			double psiForThisEntry = Double.valueOf(psiForThisEntryAsStr.trim());
			
			// Parse the rotamer probability which is located at position 37 to 44.
			String probabilityStr = ""+dbLine[37]+dbLine[38]+dbLine[39]+dbLine[40]+dbLine[41]+dbLine[42]+dbLine[43]+dbLine[44];
			double myProbability = Double.parseDouble(probabilityStr.trim());
			
			// The first char for each of the dihedrals is, respectively at 47, 55, 63, 71
			int dihedralIndices [] = {47, 55, 63, 71}; 
			// Parse all the dihedrals for each entry
			double chiAngles[] = new double [numDih];
			for(int chiIx = 0; chiIx < numDih; chiIx++ ){
				String dihedralStr = ""+dbLine[dihedralIndices[chiIx]]+dbLine[dihedralIndices[chiIx]+1]+dbLine[dihedralIndices[chiIx]+2]+
						dbLine[dihedralIndices[chiIx]+3]+dbLine[dihedralIndices[chiIx]+4]+dbLine[dihedralIndices[chiIx]+5];
				chiAngles[chiIx] = Double.parseDouble(dihedralStr.trim());
			}
			// Add a new Rotamer.	
			if(myProbability >= EnvironmentVars.DUNBRACK_PROBABILTY_CUTOFF){
				BBDepRotamersForType bbRotForType = allDunbrackRotamersMap.get(aaNameNewRot);
				// The dunbrack rotamer library contains some rotamers not defined by our template.
				// Those rotamers are ignored.
				if(bbRotForType != null){						
					bbRotForType.addNewRotamer(phiForThisEntry, psiForThisEntry, chiAngles);
				}
			}
		}
	
		// For each amino acid in the template				
        for(ResidueTemplate myTemplate : templates){
        	myTemplate.setRLphiPsiResolution(dunbrackResolution); 
        	myTemplate.setNumberOfPhiPsiBins(binsDunbrack);
        	myTemplate.initializeRotamerArrays();
        	for(int phiBin = 0; phiBin < binsDunbrack; phiBin++){
        		for(int psiBin = 0; psiBin < binsDunbrack; psiBin++){
        			double rotamers [][] = (allDunbrackRotamersMap.get(myTemplate.name)).rotamersForType[phiBin][psiBin].convertRotamersToArray();
        			myTemplate.setRotamericDihedrals(rotamers, phiBin, psiBin); 
        			myTemplate.setNumRotamers(rotamers.length, phiBin, psiBin);
        			if(rotamers.length > 0 ){
        				myTemplate.setNumDihedrals(rotamers[0].length);
        			}
            	}
        	}
        	

        }
	}        
        
        //Will need volume stuff for K*...
        /*
	//reads in the rotamer volume data from volFilename
	public void loadVolFile (String volFilename){
		
		try {
			readRotVol(volFilename);
		}
		catch (Exception e){
			System.out.println("Rotamer volumes file not found. Computing it..");
			computeAAVolumes(volFilename);
			try {
				readRotVol(volFilename);
			}
			catch (Exception ex){
				throw new Error(ex);
			}
		}
	}
	
	//Read in all of the rotamer volumes for all amino acids from the volFilename file
	private void readRotVol(String volFilename) throws Exception {
		
		FileInputStream is = new FileInputStream( volFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;
		
		rotamerVolumes = new double[numAAallowed][];
		
		is = new FileInputStream( volFilename );
		bufread = new BufferedReader(new InputStreamReader(is));
		
		// Skip over comments (lines starting with !)
		curLine = bufread.readLine();
		while( curLine.charAt(0) == '!' )
			curLine = bufread.readLine();
  		
		int curResult = 0;
		while( curLine != null ) {
							
			int aaIndex = getAARotamerIndex(StringParsing.getToken(curLine,1));
			
			int numRotVol = numRotamers[aaIndex];
			if (numRotamers[aaIndex]==0)
				numRotVol++;
			
			rotamerVolumes[aaIndex] = new double[numRotVol];
			for (int j=0; j<numRotVol; j++)
				rotamerVolumes[aaIndex][j] = new Double(StringParsing.getToken(curLine,j+2)).doubleValue();
		
			curLine = bufread.readLine();
			curResult++;
		}
		
		bufread.close();
		
		if (curResult!=numAAallowed){
			throw new Error("not all amino acid types read from rotamer volumes file");
		}
	}	
	
	// Uses the VolModule class to calculate the volume of each rotamer of each amino acid
	public void computeAAVolumes(String volFileName) {

		Molecule m = new Molecule();

		Amber96PolyPeptideResidue ppr = new Amber96PolyPeptideResidue();
		StrandRotamers LR = null;

		Residue res = ppr.getResidue("ala");
		//res.fullName = "ALA  "; 
		m.addResidue(0,res);
		
		VolModule sm = new VolModule(m);
		sm.setResidueTreatment(0,1);
		
		LR = new StrandRotamers(rotFile,m.strand[0]);		

		PrintStream printStream = setupOutputFile(volFileName);

		String aanames[] = getAAtypesAllowed();
		int numAAs = getNumAAallowed();
		
		for(int i=0;i<numAAs;i++){
			if(canMutate)
				LR.changeResidueType(m,0,aanames[i],true);
			printStream.print(aanames[i] + " ");
			System.out.println(aanames[i] + " ");
			if(getNumRotamers(aanames[i])==0){		// ALA or GLY
				double vol = sm.getMoleculeVolume(0.25f,0.0f);
				printStream.print(vol + " ");
				System.out.println(vol + " ");
			}			
			for(int j=0;j<getNumRotamers(aanames[i]);j++){
				LR.applyRotamer(m,0,j);
				double vol = sm.getMoleculeVolume(0.25f,0.0f);
				printStream.print(vol + " ");
				System.out.println(vol + " ");
			}
			printStream.println();
		}
		printStream.close();		
	}
        */
	
}
