package edu.duke.cs.osprey.kstar.pfunc;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.pfunc.impl.PF1NMTPCPMCache;
import edu.duke.cs.osprey.kstar.pfunc.impl.PF1NNoCache;
import edu.duke.cs.osprey.kstar.pfunc.impl.PF1NPCPMCache;
import edu.duke.cs.osprey.kstar.pfunc.impl.PF1NPMCache;
import edu.duke.cs.osprey.kstar.pfunc.impl.PF1NUBNM;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFMNPCPMCache;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Chooses the type of partition function implementation, depending upon user
 * supplied value
 */
public class PFFactory {

	public static PFAbstract getPartitionFunction( String implementation,
			ArrayList<String> sequence, ConfigFileParser cfp, SearchProblem sp,
			double EW_I0 ) {

		switch( implementation.toLowerCase() ) {

		case "1nnocache":
			return new PF1NNoCache( sequence, cfp, sp, EW_I0 );
		
		case "1nubnm":
			return new PF1NUBNM( sequence, cfp, sp, EW_I0 );
			
		case "1npmcache":
			return new PF1NPMCache( sequence, cfp, sp, EW_I0 );

		case "1npcpmcache":
			return new PF1NPCPMCache( sequence, cfp, sp, EW_I0 );

		case "1nmtpcpmcache":
			return new PF1NMTPCPMCache( sequence, cfp, sp, EW_I0 );

		case "mnpcpmcache":
			return new PFMNPCPMCache( sequence, cfp, sp, EW_I0 );

		default:
			throw new RuntimeException("ERROR: specified value of parameter pFuncMethod is invalid");

		}
	}

}
