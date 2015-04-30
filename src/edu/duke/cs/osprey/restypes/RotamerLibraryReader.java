/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.control.EnvironmentVars;
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
	public static void readRotLibrary(String rotFilename, ResidueTemplateLibrary templateLib) {
            
            try {
                
		//String volFilename = rotFile + ".vol";
		
		// HANDLE THE NORMAL AAs	
		FileInputStream is = new FileInputStream( EnvironmentVars.getDataDir() + rotFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;

		// Skip over comments (lines starting with !)
		curLine = bufread.readLine();
		while( curLine.charAt(0) == '!' )
			curLine = bufread.readLine();
                
                curLine = bufread.readLine();//skip over number of residues in library
                
		
	  	while( curLine != null ) {
			if(curLine.charAt(0) == '!'){
				curLine = bufread.readLine();
				continue;
			}
			
	  		String aaName = StringParsing.getToken(curLine,1);
			int numDihedrals = (new Integer(StringParsing.getToken(curLine,2))).intValue();
			int numRotamers = (new Integer(StringParsing.getToken(curLine,3))).intValue();
			if (numRotamers<0)
                            numRotamers = 0;
			
			String dihedralAtomNames[][] = new String[numDihedrals][4];//atoms involved in each dihedral
			double rotamerValues[][] = new double[numRotamers][numDihedrals];//dihedral values at each rotamer
                        

			// Read in the actual dihedrals
			for(int q=0;q<numDihedrals;q++) {
				curLine = bufread.readLine();
				dihedralAtomNames[q][0] = StringParsing.getToken(curLine,1);
				dihedralAtomNames[q][1] = StringParsing.getToken(curLine,2);
				dihedralAtomNames[q][2] = StringParsing.getToken(curLine,3);
				dihedralAtomNames[q][3] = StringParsing.getToken(curLine,4);
			}
			// Read in the actual rotamers
			for(int q=0;q<numRotamers;q++) {
				curLine = bufread.readLine();
				for(int w=0;w<numDihedrals;w++) {
					rotamerValues[q][w] = (new Integer(StringParsing.getToken(curLine,(w+1)))).intValue();
				}
			}
			
                        //we've now got all the rotamer information for the residue type aaName
                        //record it in all templates for aaName in the template library
                        boolean foundTemplate = false;
                        
                        for(ResidueTemplate template : templateLib.templates){
                            if(template.name.equalsIgnoreCase(aaName)){
                                
                                //record information in template
                                template.numDihedrals = numDihedrals;
                                template.numRotamers = numRotamers;
                                template.rotamericDihedrals = rotamerValues;
                                
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
			templateLib.totalNumRotamers += numRotamers;
			if (numRotamers<=0) //ALA or GLY
				templateLib.totalNumRotamers += 1;
                        
			curLine = bufread.readLine();
		}
	  	
		bufread.close();
            }
            catch(Exception e){
                e.printStackTrace();
                throw new RuntimeException("ERROR reading rotamer library: "+e.getMessage());
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
			catch (Exception e1){
				System.out.println("ERROR: "+e1);
				System.exit(1);
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
			System.out.println("ERROR: not all amino acid types read from rotamer volumes file");
			System.exit(1);
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
