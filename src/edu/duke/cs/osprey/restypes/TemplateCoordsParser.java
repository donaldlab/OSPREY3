package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

import java.util.Iterator;
import java.util.List;

public class TemplateCoordsParser {

	public final List<ResidueTemplate> templates;

	public TemplateCoordsParser(List<ResidueTemplate> templates) {
		this.templates = templates;
	}

	public void parse(String text) {

		Iterator<String> lines = FileTools.parseLines(text).iterator();

		String curLine = null;

		curLine = lines.next();
		while (curLine.startsWith("#") || curLine.isEmpty()) {
			curLine = lines.next();
		}
		boolean foundRes = false;
		boolean foundAtom = false;
		while (curLine != null) {
			String resName = StringParsing.getToken(curLine,1);
			int numAtoms = Integer.parseInt(StringParsing.getToken(curLine,2));
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


	public static String writeTemplateCoords(Residue res, String templateName){
		//Write out a residue's coords in template format
		StringBuilder sb = new StringBuilder();
		sb.append(templateName+" "+res.atoms.size()+"\n");
		for(Atom at : res.atoms){
			double coords[] = at.getCoords();
			sb.append(at.name);
			for(double c : coords)
				sb.append(" "+c);
			sb.append("\n");
		}
		sb.append("ENDRES\n");
		return sb.toString();
	}

}
