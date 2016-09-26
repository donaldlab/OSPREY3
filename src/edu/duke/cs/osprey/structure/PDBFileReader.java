/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.kstar.KSTermini;
import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class PDBFileReader {

	public static Molecule readPDBFile(String PDBFile) {
		return readPDBFile(PDBFile, null);
	}
	
	public static Molecule readPDBFile( String PDBFile, KSTermini termini ){
		//Take pretty much verbatim from PDBChemModel
		//if templates not null, four things we may decide to do (should give options):
		//1. Assign templates to residues 2. Rename atoms in matching residues to match templates
		//(may even want to rename residues, for example to handle his better)
		//3. Delete residues not matching any templates 4. Reconstruct residues partially matching templates
		//need to thing hard about BB type genericness
		//uses of atom renaming: 1. template matching for forcefield (can be done w/o renaming, just record matching/atom type)
		//2. DOF applications may use atom names.  Two ways this happens: 
		//-backbone atoms are used for mutations and perturbations.  So need standardized backbone naming scheme
		//Can create templates for these.  
		//-Sidechain atoms must match the rotamer library--for use when applying sidechain dihedrals.  
		//(and bond angle changes if we end up implementing those)
		//-SO LETS HAVE A SEPARATE TEMPLATE SCHEME FOR FLEXIBILITY, PROBABLY HARD-CODED FOR BB
		//AND THEN WE CAN HAVE FORCE FIELD TEMPLATES, WHICH MUST MATCH THIS SCHEME


		//read from PDB file...
		//make each new residue as list of atoms
		//assign a template to it (have flag to either delete or leave as non-interacting MISC res if can't assign)
		//rename any templated atoms (though might store PDBName too...)
		//idealize sidechains if indicated

		ArrayList<Integer> filter = new ArrayList<>();
		
		Molecule m = new Molecule();

		try {

			FileInputStream is = new FileInputStream(PDBFile);
			BufferedReader bufread = new BufferedReader(new InputStreamReader(is));

			String curLine = bufread.readLine();

			ArrayList<String> helixStarts = new ArrayList<>();//Residues where helices start
			ArrayList<String> helixEnds = new ArrayList<>();//Residues where they end
			ArrayList<Character> helixChains = new ArrayList<>();
			ArrayList<String> sheetStarts = new ArrayList<>();
			ArrayList<String> sheetEnds = new ArrayList<>();
			ArrayList<Character> sheetChains = new ArrayList<>();
			//So for each helix/sheet, 
			//we record its starting and ending residue numbers and its chain


			ArrayList<Atom> curResAtoms = new ArrayList<>();
			ArrayList<double[]> curResCoords = new ArrayList<>();//coordinates for these atoms
			String curResFullName = "NONE";

			while(curLine!=null){

				// First pad line to 80 characters
				int lineLen = curLine.length();
				for (int i=0; i < (80-lineLen); i++)
					curLine += " ";
                                
                                if(curLine.startsWith("MODEL"))
                                    System.out.println("Warning: OSPREY doesn't understand PDB files with multiple models. ");

				if ( (curLine.regionMatches(true,0,"ATOM  ",0,6)) || (curLine.regionMatches(true,0,"HETATM",0,6)) ){

					if( EnvironmentVars.deleteNonTemplateResidues ){//Ignore alternates other than A; treat A alternates as the real structure
						char alt = curLine.charAt(16);//This specifies which alternate the atom is (space if not an alternate)
						if( ( alt != ' ' ) && ( alt != 'A' ) ){
							curLine = bufread.readLine();
							continue;//skip the line and  go to the next one
						}
					}

					String fullResName = fullResidueName(curLine);

					if( (!fullResName.equalsIgnoreCase(curResFullName)) && !curResAtoms.isEmpty() ){
						
						Residue newRes = new Residue( curResAtoms, curResCoords, curResFullName, m );

						if(termini != null && !termini.contains(newRes)) filter.add(m.residues.size());
						
						m.appendResidue(newRes);
						
						curResAtoms = new ArrayList<>();
						curResCoords = new ArrayList<>();
					}

					curResFullName = fullResName;

					readAtom(curLine,curResAtoms,curResCoords);
				}

				else if(curLine.regionMatches(true,0,"HELIX  ",0,7)){//Read helix records
					helixStarts.add( curLine.substring(21,25).trim() );
					helixEnds.add( curLine.substring(33,37).trim() );
					helixChains.add(curLine.charAt(19));
				}
				else if(curLine.regionMatches(true,0,"SHEET  ",0,7)){
					sheetStarts.add( curLine.substring(22,26).trim() );
					sheetEnds.add( curLine.substring(33,37).trim() );
					sheetChains.add(curLine.charAt(21));
				}

				curLine = bufread.readLine(); 
			}

			//make last residue
			if( ! curResAtoms.isEmpty() ){
				Residue newRes = new Residue( curResAtoms, curResCoords, curResFullName, m );
				if(termini != null && !termini.contains(newRes)) filter.add(m.residues.size());
				m.appendResidue(newRes);
			}


			bufread.close();  // close the buffer

			deleteFilteredResidues(m, filter);
			
			//Assign the secondary structure we have read
			assignSecStruct(m, helixStarts, helixEnds, helixChains, sheetStarts, sheetEnds, sheetChains);
		}
		catch(FileNotFoundException e){
			System.err.println("ERROR: PDB file not found: "+PDBFile);
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
		catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		//assign proline puckers?  if treated as DOF might handle along with regular dihedrals
		if(EnvironmentVars.assignTemplatesToStruct)
			assignTemplates(m);

		return m;
	}
	
	
	private static void deleteFilteredResidues(Molecule m, ArrayList<Integer> filter) {
		for( int index = 0; index < filter.size(); ++index ) {
			int resIndex = filter.get(index);
			
			m.deleteResidue(resIndex);
			
			for( int index2 = index; index2 < filter.size(); ++index2 )
				filter.set(index2, filter.get(index2)-1);
		}
	}



	static void assignTemplates(Molecule m){

		for(int resNum=m.residues.size()-1; resNum>=0; resNum--){
			//go through residues backwards so we can delete some if needed

			Residue res = m.residues.get(resNum);

			DAminoAcidHandler.tryRenamingAsD(res);//We accept D-amino acid named using the usual L names, but must change them here
			//so the right template name is used

			boolean templateAssigned = res.assignTemplate();

			if(EnvironmentVars.deleteNonTemplateResidues && !templateAssigned){
				//residue unrecognized or incomplete...delete it
				m.deleteResidue(resNum);
			}
		}

		HardCodedResidueInfo.markInterResBonds(m);//assigning templates marks intra-res bonds; we can now mark inter-res too
	}




	//parsing ATOM lines from a PDB file
	static String fullResidueName(String line){
		return line.substring(17,27);
	}


	static void readAtom(String curLine, ArrayList<Atom> atomList, ArrayList<double[]> coordList) {
		//read an ATOM line and store it in the list of atoms and of atom coordinates

		String tmpStg = curLine.substring(6,11).trim();  // Snag atom serial number
		int modelAtomNumber = (new Integer(tmpStg)).intValue();
		String atomName = curLine.substring(12,16).trim();  // Snag atom name

		tmpStg = curLine.substring(30,38);  // Snag x coord
		double x = (double) new Double(tmpStg).doubleValue();
		tmpStg = curLine.substring(38,46);  // Snag y coord
		double y = (double) new Double(tmpStg).doubleValue();
		tmpStg = curLine.substring(46,54);  // Snag z coord
		double z = (double) new Double(tmpStg).doubleValue();

		tmpStg = curLine.substring(60,66).trim();  // Snag B-factor
		double BFactor=0;
		if(!tmpStg.isEmpty())
			BFactor = (double) new Double(tmpStg).doubleValue();

		String elementType = curLine.substring(76,78).trim();  // Snag atom elementType
		// If we can't get element type from substring(76,78) snag
		//  the first character of the atom name
		if (elementType.equalsIgnoreCase(""))
			elementType = getEleType(curLine.substring(12,15));

		Atom newAtom = new Atom(atomName, elementType, BFactor, modelAtomNumber);
		double coords[] = new double[] {x,y,z};
		atomList.add(newAtom);
		coordList.add(coords);
	}



	// This function pulls the element type from
	//  the atom name
	private static String getEleType(String str){

		int start=0, end=-1;
		int i=0;
		while( (str.charAt(i)==' ') || ((str.charAt(i)>='0') && (str.charAt(i)<='9')) ) {
			i++;
		}
		start = i;
		end = i++;
		if (i<str.length())
			if((str.charAt(i)>='a') && (str.charAt(i)<='z'))
				end = i;
		return(str.substring(start,end+1));
	}


	private static void assignSecStruct(Molecule m, ArrayList<String> helixStarts, ArrayList<String> helixEnds,
			ArrayList<Character> helixChains, ArrayList<String> sheetStarts, ArrayList<String>sheetEnds,
			ArrayList<Character> sheetChains){

		//this isn't the most efficient, but I doubt it'll be a bottleneck
		int curSecStruct = Residue.LOOP;
		int curHelixNum = -1;//if we're in a helix or sheet, keep track of which one
		int curSheetNum = -1;

		for(Residue res : m.residues){
			//go through residues in molecule in order
			String curResNum = res.getPDBResNumber();
			char curChain = res.fullName.charAt(4);

			if(curSecStruct == Residue.LOOP){
				//look for the start of a helix or sheet
				for(int helixNum=0; helixNum<helixStarts.size(); helixNum++){
					if(helixChains.get(helixNum) == curChain){
						if(helixStarts.get(helixNum).equalsIgnoreCase(curResNum)){
							curHelixNum = helixNum;
							curSecStruct = Residue.HELIX;
						}
					}
				}
				for(int sheetNum=0; sheetNum<sheetStarts.size(); sheetNum++){
					if(sheetChains.get(sheetNum) == curChain){
						if(sheetStarts.get(sheetNum).equalsIgnoreCase(curResNum)){
							curSheetNum = sheetNum;
							curSecStruct = Residue.SHEET;
						}
					}
				}
			}

			res.secondaryStruct = curSecStruct;
			//if we're at the end residue of a helix or loop, stay in it for this residue
			//now check if we're at the end so the next residue can be something else

			if(curSecStruct == Residue.HELIX){
				//look for end of helix
				if(helixEnds.get(curHelixNum).equalsIgnoreCase(curResNum)){
					curHelixNum = -1;
					curSecStruct = Residue.LOOP;
				}
			}
			else if (curSecStruct == Residue.SHEET){ //sheet
				if(sheetEnds.get(curSheetNum).equalsIgnoreCase(curResNum)){
					curSheetNum = -1;
					curSecStruct = Residue.LOOP;
				}
			}
		}
	}


}
