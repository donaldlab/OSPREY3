/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

//This class stores parameters from input files
//Control package classes will probably each want a param set
//which will control its choice of algorithms, and will be used to initialize other classes 
//like EnergyMatrix and A* according to user settings


import edu.duke.cs.osprey.handlempi.MPIMaster;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.FileTools.FilePathRoot;
import edu.duke.cs.osprey.tools.FileTools.PathRoot;
import edu.duke.cs.osprey.tools.FileTools.ResourcePathRoot;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 * Handles reading in and managing parameter/value pairs from the input configuration files
 */
public class ParamSet implements Serializable {
	
	private static final long serialVersionUID = 4364963601242324780L;
	
	private static class Entry implements Serializable {
		
		private static final long serialVersionUID = 521766449139228629L;
		
		public final String value;
		public final PathRoot root;
		
		public Entry(String value, PathRoot root) {
			this.value = value;
			this.root = root;
		}
	}
	
	private Map<String,Entry> params;//map parameter/value pairs
	//parameter names will be stored as all upper-case, to avoid confusion

	private Map<String,Entry> defaultParams;
	
	private Map<String,Entry> wildcardDefaults;
	//Handles default params of the form "BLABLA* 1" (would map BLABLA -> 1)
	//then for example if we needed a default for parameter BLABLABLA, it would return 1

	private boolean isVerbose;
	
	public ParamSet() {
		
		params = new HashMap<>();
		defaultParams = new HashMap<>();
		wildcardDefaults = new HashMap<>();
		isVerbose = true;
		
		// load the defaults
		ResourcePathRoot configRoot = new ResourcePathRoot("/config");
		loadParams(configRoot, "defaults.cfg", defaultParams);
		
		// handle wildcard defaults
		for (String param : defaultParams.keySet()) {
			if (param.endsWith("*")) {
				String wildcard = param.substring(0, param.length() - 1);
				wildcardDefaults.put(wildcard, defaultParams.get(param));
			}
		}
	}
	
	public ParamSet(ParamSet other) {
		params = new HashMap<>(other.params);
		defaultParams = new HashMap<>(other.defaultParams);
		wildcardDefaults = new HashMap<>(other.wildcardDefaults);
	}

	public void setVerbosity(boolean val) {
		isVerbose = val;
	}
	
	public boolean isVerbose() {
		return isVerbose;
	}

	
	public void addParamsFromFile(String path) {
		addParamsFromFile(new File(path));
	}
	
	public void addParamsFromFile(File file) {
		FilePathRoot root = new FilePathRoot(file.getAbsoluteFile().getParentFile());
		file = root.makeRelative(file);
		addParams(root, file.getPath());
	}
	
	public void addParamsFromResource(String path) {
		ResourcePathRoot root = ResourcePathRoot.parentOf(path);
		path = root.makeRelative(path);
		addParams(root, path);
	}
	
	public void addParams(PathRoot root, String path) {
		loadParams(root, path, params);
	}
	
	private static void loadParams(PathRoot root, String path, Map<String,Entry> paramMap) {
		// load all parameters for cfg file fName into the map paramMap
		
		// read the file
		String text = root.read(path);
		
		// First attempt to open and read the config file
		// read line-by-line
		for (String line : FileTools.parseLines(text)) {
			
			// strip comments
			int commentStartPos = line.indexOf('%');
			if (commentStartPos >= 0) {
				line = line.substring(0, commentStartPos);
			}
			
			// skip blank lines
			if (line.isEmpty()) {
				continue;
			}
			
			// parse the param
			String paramName = StringParsing.getToken(line, 1).trim();
			paramName = paramName.toUpperCase();
			if (paramMap.containsKey(paramName)) {
			
				// duplicate param
				throw new RuntimeException("parameter " + paramName + " already read");
				
			} else {
				
				// new param
				String paramVal = line.substring(paramName.length() + 1);
				paramMap.put(paramName, new Entry(paramVal, root));
			}
		}
	}

	
	//Sets the value of parameter paramName to newValue
	public void setValue(String paramName, String newValue) {
		setValue(paramName, newValue, null);
	}
	
	public void setValue(String paramName, String newValue, PathRoot root) {
		paramName = paramName.toUpperCase();
		params.put(paramName, new Entry(newValue, root));
	}

        
        
        //Methods to get parameter values
        //We first try to get a value from params (loaded from user's cfg files)
        //If there is none, we try to get a value from defaultParams (loaded from default cfg)
        //If there is none there either, this is an error
        
    public String getValue(String paramName){
            
            paramName = paramName.toUpperCase();
            Entry entry = params.get(paramName);
            
            if(entry==null){
                entry = defaultParams.get(paramName);
                
                if(entry==null){//See if this is a wildcard param
                   entry = findWildcard(paramName);
                }
                
                if(entry==null)//still null
                    throw new RuntimeException("ERROR: Parameter "+paramName+" not found");

                if (isVerbose) {
                    MPIMaster.printIfMaster("Parameter "+paramName+" not set. Using default value "+entry.value);
                }
            }
            else {
                if (isVerbose) {
                    MPIMaster.printIfMaster("Parameter "+paramName+" set to "+entry.value);
                }
            }
            
            return entry.value.trim();
    }
    
    private Entry findWildcard(String paramName) {
        for(String wildcard : wildcardDefaults.keySet()){
            if(paramName.startsWith(wildcard)){
                return wildcardDefaults.get(wildcard);
            }
        }
        return null;
    }
    
    public PathRoot getRoot(String paramName) {
        
        paramName = paramName.toUpperCase();
        Entry entry = params.get(paramName);
        
        if (entry == null) {
            entry = defaultParams.get(paramName);
        }
        
        if (entry == null) {
            entry = findWildcard(paramName);
        }
        
        if (entry != null) {
            return entry.root;
        }
        
        return null;
    }
        
        
        //The following methods return a parameter expected to be a certain non-string type
        //It's an error if they're not that type
        public int getInt(String paramName){
            String val = getValue(paramName);
            try {
                return Integer.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as an integer");
            }
        }
        
        public boolean getBool(String paramName){
            String val = getValue(paramName);
            try {
                return Boolean.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a boolean");
            }
        }
        
        public double getDouble(String paramName){
            String val = getValue(paramName);
            try {
                return Double.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a double");
            }
        }
        
        
        public String getRunSpecificFileName(String paramName, String suffix){
            //This is for a special kind of parameter whose value is 
            //the name of a run-specific file.  So by default the value is runName + suffix
            paramName = paramName.toUpperCase();
            
            if(params.containsKey(paramName))
                return params.get(paramName).value;
            else
                return getValue("RUNNAME") + suffix;
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
        
        

      //Methods to get parameter values
        //Default methods can be used if not set; if null default then return an error
        
	public String getValue(String paramName, String defaultVal){
            
            paramName = paramName.toUpperCase();
            Entry entry = params.get(paramName);
            String val;
            
            if(entry==null){
                if(defaultVal==null)//no default...must be set
                    throw new RuntimeException("ERROR: Parameter "+paramName+" not found");
                
                val = defaultVal;
                if(isVerbose) {
                    MPIMaster.printIfMaster("Parameter "+paramName+" not set. Using default value "+val);
                }
            }
            else {
                val = entry.value;
                if(isVerbose) {
                    MPIMaster.printIfMaster("Parameter "+paramName+" set to "+val);
                }
            }
            
            return val.trim();
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
        
        public File getFile(String paramName) {
            String path = getValue(paramName);
            PathRoot root = getRoot(paramName);
            if (root instanceof FilePathRoot) {
                return ((FilePathRoot)root).resolve(new File(path));
            }
            throw new Error("param " + paramName + " was loaded from a resource, so it can't reference files");
        }

		public String readPath(String paramName) {
			String path = getValue(paramName);
			PathRoot root = getRoot(paramName);
			return root.read(path);
		}
        
}
