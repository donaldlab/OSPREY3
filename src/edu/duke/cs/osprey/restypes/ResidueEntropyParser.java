package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

public class ResidueEntropyParser {

	private ResidueEntropies resent;

	public ResidueEntropyParser(ResidueEntropies resent) {
		this.resent = resent;
	}

	public void parse(String text) {
		//It is convenient to load residue entropies into a hash map, rather than
		//into template objects, because they correspond to template names
		for (String line : FileTools.parseLines(text)) {

			// skip comments
			if (line.startsWith("%")) {
				continue;
			}

			String resType = StringParsing.getToken(line,1);
			double entropy = Double.parseDouble(StringParsing.getToken(line,2));
			resent.set(resType, entropy);
		}
	}
}
