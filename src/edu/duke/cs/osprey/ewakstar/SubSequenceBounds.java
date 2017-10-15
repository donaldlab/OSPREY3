package edu.duke.cs.osprey.ewakstar;

import java.util.StringTokenizer;

import edu.duke.cs.osprey.control.ConfigFileParser;

public class SubSequenceBounds {

	public int start;
	public int end;
	public int strand;
	private int offset;
	
	public SubSequenceBounds(ConfigFileParser cfp, int strand) {
		this.strand = strand;
		start = end = offset = 0;
		
		for(int i = 0; i <= strand; ++i) {
			start = end;
			
			String flexRes = cfp.getParams().getValue("STRANDMUT"+i);
			StringTokenizer st = new StringTokenizer(flexRes);
			while(st.hasMoreTokens()) {
				offset += st.nextToken().length()+4;//residue+-+3 letter aa code
				offset += 1;//add space
			}
			
			end = offset;
		}
	}
	
}
