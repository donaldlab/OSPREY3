package edu.duke.cs.osprey.kstar.pfunc;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFNew03;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTrad;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFNew02;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFNew01;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFnew00;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFNew04;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Chooses the type of partition function implementation, depending upon user
 * supplied value
 */
public class PFFactory {

	public static PFAbstract getPartitionFunction( String implementation, int strand, 
			ArrayList<String> sequence, ArrayList<Integer> absolutePos, String checkPointPath, 
			String searchProblemName, ConfigFileParser cfp, SearchProblem sp ) {

		switch( implementation.toLowerCase() ) {

		case "trad":
			return new PFTrad( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
		
		case "new00":
			return new PFnew00( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
			
		case "new01":
			return new PFNew01( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "new02":
			return new PFNew02( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "new03":
			return new PFNew03( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "new04":
			return new PFNew04( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		default:
			throw new RuntimeException("ERROR: specified value of parameter pFuncMethod is invalid");

		}
	}

}
