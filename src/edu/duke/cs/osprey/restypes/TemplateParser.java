package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class TemplateParser {

	public final ForcefieldParams ffparams;

	public TemplateParser(ForcefieldParams ffparams) {
		this.ffparams = ffparams;
	}

	public List<ResidueTemplate> parse(String text) {

		List<ResidueTemplate> templates = new ArrayList<>();

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

		return templates;
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
			at.type = ffparams.atomTypeToInt(at.forceFieldType);
			at.charge = Double.parseDouble(StringParsing.getToken(curLine,11));

			//the template only records bonds within the residue

			//KER: The first atom is bonded to a dummy atom so we can't include that
			//KER: in the bond list, so check atom is >= 0
			int atomBondedTo = Integer.parseInt(StringParsing.getToken(curLine,5)) - dumPresent;
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


		return new ResidueTemplate(templateRes, templateName);
	}
}