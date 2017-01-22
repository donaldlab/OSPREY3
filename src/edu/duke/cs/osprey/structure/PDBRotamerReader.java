package edu.duke.cs.osprey.structure;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.restypes.PositionSpecificRotamerLibrary;
import edu.duke.cs.osprey.restypes.ResidueTemplate;

/***
 * Simple class to read in a PDB file and spit out rotamers for 
 * both the main conformation and alternates, if any.
 * @author JJ
 *
 */
public class PDBRotamerReader {
    

    public static void createTemplates(Molecule m, PositionSpecificRotamerLibrary library, String PDBFileName)
    {
        GenericResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        
        Map<Integer, Map<String, ArrayList<Residue>>> positionSpecificRotamers = new HashMap<>();
        try {

            FileInputStream is = new FileInputStream(PDBFileName);
            BufferedReader bufread = new BufferedReader(new InputStreamReader(is));

            String curLine = bufread.readLine();

            ArrayList<Atom> curResAtoms = new ArrayList<>();
            Map<Character, ArrayList<Atom>> alternateAtoms = new HashMap<>();
            ArrayList<double[]> curResCoords = new ArrayList<>();//coordinates for these atoms
            Map<Character, ArrayList<double[]>> alternateResidueCoords = new HashMap<>();

            String curResFullName = "NONE";

            while(curLine!=null){

                // First pad line to 80 characters
                int lineLen = curLine.length();
                for (int i=0; i < (80-lineLen); i++)
                    curLine += " ";
                System.out.println("Line: "+curLine);

                if ( (curLine.regionMatches(true,0,"ATOM  ",0,6)) || (curLine.regionMatches(true,0,"HETATM",0,6)) ){

                    char alt = curLine.charAt(16);//This specifies which alternate the atom is (space if not an alternate)
                    

                    if(alt != ' ')
                    {
                        String fullResName = fullResidueName(curLine);
                        int residueIndex = getResidueIndex(curLine);
                        System.out.println("Residue "+residueIndex+", alternate "+alt);
                        
                        if( (!fullResName.equalsIgnoreCase(curResFullName)) && !curResAtoms.isEmpty() ){

                            

                            for(char c: alternateAtoms.keySet())
                            {
                                if(alternateAtoms.get(c).size() > 0)
                                {
                                    Residue alternateConformation = new Residue(alternateAtoms.get(c), 
                                            alternateResidueCoords.get(c), curResFullName, m);
                                    try
                                    {
                                        boolean success = alternateConformation.assignTemplate(templateLib);
                                        System.out.println("Template assignment succeeded?:"+success);
                                        if(!success)
                                        {
                                            System.out.println("Assignment failed.");
                                            alternateConformation.assignTemplate(templateLib);
                                            continue;
                                        }
                                        m.addAlternate(residueIndex, alternateConformation);
                                        //library.addRotamer(residueIndex, alternateConformation.template.name, alternateConformation);
                                        addRotamer(residueIndex, positionSpecificRotamers,alternateConformation); 

                                    }
                                    catch (Exception e)
                                    {
                                        e.printStackTrace();
                                    }


                                }
                                alternateAtoms.put(c, new ArrayList<>());
                                alternateResidueCoords.put(c,new ArrayList<>());
                            }

                            curResAtoms = new ArrayList<>();
                            curResCoords = new ArrayList<>();
                        }

                        curResFullName = fullResName;

                        readAtom(curLine,curResAtoms,curResCoords);
                        if(!alternateAtoms.containsKey(alt))
                        {
                            alternateAtoms.put(alt, new ArrayList<>());
                            alternateResidueCoords.put(alt, new ArrayList<>());
                        }
                        readAtom(curLine,alternateAtoms.get(alt),alternateResidueCoords.get(alt));
                        System.out.println("Added atoms to "+alt+":");
                        System.out.println(alternateAtoms.get(alt));
                        System.out.println(alternateResidueCoords.get(alt));

                    }
                }


                curLine = bufread.readLine(); 
            }

            //make last residue
            if( ! curResAtoms.isEmpty() ){
                Residue newRes = new Residue( curResAtoms, curResCoords, curResFullName, m );
                m.appendResidue(newRes);
            }



            bufread.close();  // close the buffer
            


        }
        catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
        
        // At this point, we should conver the stored conformations into a set of position-specific
        // ResidueTemplates, assigning the number of rotamers, number of phipsibins, number of dihedrals,
        // and then storing each template at its appropriate position.
        for(Integer residueIndex : positionSpecificRotamers.keySet())
        {
            for(String resType : positionSpecificRotamers.get(residueIndex).keySet())
            {
            	List<Residue> residues = positionSpecificRotamers.get(residueIndex).get(resType);
                library.addResidueTemplate(residueIndex, resType, ResidueTemplate.makeFromResidueConfs(residues));
            }

        }
    }
    
    private static void addRotamer (int residueIndex,
            Map<Integer, Map<String, ArrayList<Residue>>> positionSpecificRotamers,
            Residue alternateConformation) {
        if(!positionSpecificRotamers.containsKey(residueIndex))
            positionSpecificRotamers.put(residueIndex, new HashMap<>());
        Map<String, ArrayList<Residue>> residuesAtPosition = positionSpecificRotamers.get(residueIndex);
        if(!residuesAtPosition.containsKey(alternateConformation.template.name))
            residuesAtPosition.put(alternateConformation.template.name, new ArrayList<Residue>());
        List<Residue> residuesAtPositionForAA = residuesAtPosition.get(alternateConformation.template.name);
        residuesAtPositionForAA.add(alternateConformation);
        System.out.println(residuesAtPosition.get(alternateConformation.template.name).size());
        
    }

    private static int getResidueIndex (String curLine) {
        return Integer.valueOf(curLine.substring(23,26).trim());
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

        Atom newAtom = new Atom(atomName, elementType);
        newAtom.BFactor = BFactor;
        newAtom.modelAtomNumber = modelAtomNumber;
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



    //parsing ATOM lines from a PDB file
    static String fullResidueName(String line){
        return line.substring(17,27);
    }

}
