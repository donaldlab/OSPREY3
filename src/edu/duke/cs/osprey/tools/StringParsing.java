/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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
