package edu.duke.cs.osprey.tools.branchdecomposition;
import java.util.StringTokenizer;


public class BranchDecomposition {

	public static void main(String[] args) {
		
		if (args.length<=0) {
			
			byte bytebuff[] = new byte[150];
			System.out.print("> ");
			try {
				System.in.read(bytebuff);
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(0);
			}
			String s = new String(bytebuff).trim();  // create a string from bytebuff
			args = new String[2];
			args[0] = getToken(s,1);
			args[1] = getToken(s,2);
		}
		
		else if ( (args.length==1) || (args.length>=3) ){
			System.out.println("ERROR: arguments must be: 1)input filename, and 2)output filename");
			System.exit(1);
		}
		
		try{
		    BranchDecompositionH bd = new BranchDecompositionH(args);
		}
		catch( Exception e)
		{
		    e.printStackTrace();
		}
	}

	// This function returns the xth token in string s
	private static String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				System.out.println("ERROR: Unable to access argument " + x + " from input string");
				System.exit(1);
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		
		else {
			System.out.println("ERROR: Unable to access argument " + x + " from input string");
			System.exit(1);
			return null;
		}

	} // end getToken
}
