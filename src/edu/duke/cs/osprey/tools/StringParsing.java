/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.tools;

import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *
 * @author mhall44
 */

//tools for string parsing

public class StringParsing {
    
    	// This function returns the xth token in string s
	public static String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken
        
        
        public static boolean containsIgnoreCase(ArrayList<String> list, String s){
            //Does list contain s, ignoring case?
            for(String a : list){
                if(a.equalsIgnoreCase(s))
                    return true;
            }
            
            return false;
        }
        
     
    	public static int ordinalIndexOf(String str, String c, int n) {
    	    int pos = str.indexOf(c, 0);
    	    while (n-- > 0 && pos != -1)
    	        pos = str.indexOf(c, pos+1);
    	    return pos;
    	}
    
}
