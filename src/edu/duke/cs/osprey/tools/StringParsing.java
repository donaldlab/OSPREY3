/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

