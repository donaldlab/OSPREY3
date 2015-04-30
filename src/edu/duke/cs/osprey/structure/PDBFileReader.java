/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class PDBFileReader {
    
    public static Molecule readPDBFile( String PDBFile ){
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
        
        
        Molecule m = new Molecule();
        
        try {
            
            FileInputStream is = new FileInputStream(PDBFile);
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
            
            String curLine = bufread.readLine();

            ArrayList<Integer> helixStarts = new ArrayList<Integer>();//Residues where helices start
            ArrayList<Integer> helixEnds = new ArrayList<Integer>();//Residues where they end
            ArrayList<Character> helixStrands = new ArrayList<Character>();
            ArrayList<Integer> sheetStarts = new ArrayList<Integer>();
            ArrayList<Integer> sheetEnds = new ArrayList<Integer>();
            ArrayList<Character> sheetStrands = new ArrayList<Character>();
            //COMPUTE RES 2NDARY STRUCT FROM THESE?


            ArrayList<Atom> curResAtoms = new ArrayList<>();
            ArrayList<double[]> curResCoords = new ArrayList<>();//coordinates for these atoms
            String curResFullName = "NONE";

            while(curLine!=null){
                
                // First pad line to 80 characters
                int lineLen = curLine.length();
                for (int i=0; i < (80-lineLen); i++)
                        curLine += " ";
                

                if ( (curLine.regionMatches(true,0,"ATOM  ",0,6)) || (curLine.regionMatches(true,0,"HETATM",0,6)) ){

                    if( EnvironmentVars.fixNonTemplateResidues ){//Ignore alternates other than A; treat A alternates as the real structure
                        char alt = curLine.charAt(16);//This specifies which alternate the atom is (space if not an alternate)
                        if( ( alt != ' ' ) && ( alt != 'A' ) ){
                            curLine = bufread.readLine();
                            continue;//skip the line and  go to the next one
                        }
                    }
                    
                    String fullResName = fullResidueName(curLine);
                    
                    if( (!fullResName.equalsIgnoreCase(curResFullName)) && !curResAtoms.isEmpty() ){
                        
                        Residue newRes = new Residue( curResAtoms, curResCoords, curResFullName, m );
                        newRes.indexInMolecule = m.residues.size();
                        m.residues.add(newRes);
                        
                        curResAtoms = new ArrayList<>();
                        curResCoords = new ArrayList<>();
                    }
                    
                    curResFullName = fullResName;

                    readAtom(curLine,curResAtoms,curResCoords);
                }

                else if(curLine.regionMatches(true,0,"HELIX  ",0,7)){//Read helix records
                    helixStarts.add(new Integer(curLine.substring(21,25).trim()));
                    helixEnds.add(new Integer(curLine.substring(33,37).trim()));
                    helixStrands.add(curLine.charAt(19));
                }
                else if(curLine.regionMatches(true,0,"SHEET  ",0,7)){
                    sheetStarts.add(new Integer(curLine.substring(22,26).trim()));
                    sheetEnds.add(new Integer(curLine.substring(33,37).trim()));
                    sheetStrands.add(curLine.charAt(21));
                }

                curLine = bufread.readLine(); 
            }


            bufread.close();  // close the buffer
        }
        catch(IOException e){
            System.err.println(e.getMessage());
            e.printStackTrace();
        }

        
                //assign proline puckers?  if treated as DOF might handle along with regular dihedrals

        
        if(EnvironmentVars.assignTemplatesToStruct){
            for(Residue res : m.residues)
                res.assignTemplate();//using ResidueTemplateLibrary
            
            HardCodedResidueInfo.markInterResBonds(m);//assigning templates marks intra-res bonds; we can now mark inter-res too
        }
        
        return m;
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
}
