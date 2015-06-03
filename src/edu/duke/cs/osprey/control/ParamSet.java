/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

//This class stores parameters from input files
//Control package classes will probably each want a param set
//which will control its choice of algorithms, and will be used to initialize other classes 
//like EnergyMatrix and A* according to user settings


import edu.duke.cs.osprey.handlempi.MPIMaster;
import edu.duke.cs.osprey.tools.StringParsing;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * Handles reading in and managing parameter/value pairs from the input configuration files
 */
public class ParamSet implements Serializable {
	
	private TreeMap<String,String> params = new TreeMap<>();//map parameter/value pairs
        //parameter names will be stored as all upper-case, to avoid confusion
	
	//constructor
	ParamSet(){
            
	}
	
	//Reads in all parameter pairs from the file fName and updates the params
	public void addParamsFromFile(String fName){
		
		BufferedReader bufread = null;
		String curLine = null;
		boolean done = false;
		
		// First attempt to open and read the config file
		try{
			FileInputStream is = new FileInputStream(fName);
			bufread = new BufferedReader(new InputStreamReader(is));

			curLine = bufread.readLine();

			while (curLine != null){
				done = false;
				while (!done) {
                                        if(curLine.isEmpty()){//skip blank lines
                                            curLine = bufread.readLine();
                                        }
                                        else if (curLine.charAt(0) == '%'){
						curLine = bufread.readLine();
					}
					else {
						done = true;
					}
					if (curLine == null)
						done = true;
				}
				if (curLine != null){
					String paramName = StringParsing.getToken(curLine,1).trim();
                                        paramName = paramName.toUpperCase();
                                        
                                        if(params.containsKey(paramName))
                                            throw new RuntimeException("ERROR: parameter "+paramName+" already read");
					else { // new parameter
                                            String paramVal = curLine.substring(paramName.length()+1);
                                            params.put(paramName, paramVal);
                                            curLine = bufread.readLine();
					}
				}
			}
			bufread.close();
		}
		catch(Exception ex)
		{
			System.out.println("ERROR: An error occurred reading configuration file "+fName);
                        System.out.println(ex.getMessage());
                        ex.printStackTrace();
			System.exit(1);
		}
	}

	
	//Sets the value of parameter paramName to newValue
	public void setValue(String paramName, String newValue){
            params.put(paramName, newValue);
	}

        
        
        //Methods to get parameter values
        //Default methods can be used if not set; if null default then return an error
        
	public String getValue(String paramName, String defaultVal){
            
            paramName = paramName.toUpperCase();
            String val = params.get(paramName);
            
            if(val==null){
                if(defaultVal==null)//no default...must be set
                    throw new RuntimeException("ERROR: Parameter "+paramName+" not found");
                
                val = defaultVal;
                MPIMaster.printIfMaster("Parameter "+paramName+" not set. Using default value "+defaultVal);
            }
            else
                MPIMaster.printIfMaster("Parameter "+paramName+" set to "+val);
            
            return val.trim();
	}
        
        //getting values with no defaults
        public String getValue(String paramName){
            return getValue(paramName,null);
        }
        
        //The following methods return a parameter expected to be a certain non-string type
        //It's an error if they're not that type
        public int getInt(String paramName, int defaultVal){
            String val = getValue(paramName,String.valueOf(defaultVal));
            try {
                return Integer.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as an integer");
            }
        }
        
        public boolean getBool(String paramName, boolean defaultVal){
            String val = getValue(paramName,String.valueOf(defaultVal));
            try {
                return Boolean.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a boolean");
            }
        }
        
        public double getDouble(String paramName, double defaultVal){
            String val = getValue(paramName,String.valueOf(defaultVal));
            try {
                return Double.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a double");
            }
        }
        
        
        
        
        //searching for a parameter
        public ArrayList<String> searchParams(String searchTerm){
            //return any parameter names that include searchTerm
            searchTerm = searchTerm.toUpperCase();//all parameters are upper-case
            
            ArrayList<String> matches = new ArrayList<>();//parameters matching the searchTerm
            
            for(String paramName : params.keySet()){
                if(paramName.contains(searchTerm)){
                    matches.add(paramName);
                }
            }
            
            return matches;
        }
        
        

}
