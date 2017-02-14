/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tools;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InvalidClassException;
import java.io.NotSerializableException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OptionalDataException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.commons.io.FileUtils;

/**
 *
 * @author mhall44
 */

//tools for storing objects, like an energy matrix
//We can also deep-copy objects by passing them through a stream, as if to store to/load from a file

public class ObjectIO {

    public static class BadFileException extends Exception {
        
        private static final long serialVersionUID = -9191066141220686516L;
        
        private static final String MsgFormat = "Can't load file data from %s: %s";
        
        public BadFileException(File file, String msg) {
            super(String.format(MsgFormat, file.getAbsolutePath(), msg));
        }
        
        public BadFileException(File file, String msg, Throwable cause) {
            super(String.format(MsgFormat, file.getAbsolutePath(), msg), cause);
        }
    }
    
    public static class CantWriteException extends Exception {
        
        private static final long serialVersionUID = 8152920955998224621L;
        
        private static final String MsgFormat = "Can't write file data to %s";
        
        public CantWriteException(File file, Throwable cause) {
            super(String.format(MsgFormat, file.getAbsolutePath()), cause);
        }
    }
    
    public static interface Validator<T> {
    	boolean isValid(T thing);
    }
    
    public static Object readObject( String fileName, boolean allowNull ){
        //Read the object from the file
        //If cannot (e.g., because object doesn't exist), then return null if allowNull, else 
        //raise an error
        Object inObj = null;
        
        try{
                ObjectInputStream in = new ObjectInputStream(new FileInputStream(fileName));
                inObj = in.readObject();
                in.close();
        }
        catch (Exception e){//couldn't read the object
                if (allowNull)
                    return null;
                else
                	throw new RuntimeException("ERROR: Failed to read object from file " + fileName + "\n" + e.getMessage());
        }

        return inObj;
    }
    
    
    // these read methods have some more error-checking, so callers better decide how to recover from errors
    
    public static <T> T read(String path, Class<T> type)
    throws BadFileException {
        return read(new File(path), type);
    }
    
    public static <T> T read(File file, Class<T> type)
    throws BadFileException {
        
        if (!file.exists()) {
            return null;
        }
        
        // reading objects using Java's built-in serialization is Fraught With Peril!
        // recovery from failures is possible, but the caller should decide how to do it.
        // for things that can't be saved (like stack overflows), we'll just have to hard exit the process
        // but punt problems with files to the caller using a simpler friendlier exception format
        
        try (ObjectInputStream in = new ObjectInputStream(new FileInputStream(file))) {
        
            // try to read the object and see what happens
            Object obj = in.readObject();
            
            // did we read the object we expected?
            if (obj == null) {
                throw new BadFileException(file, "file contained no data");
            }
            
            if (!type.isAssignableFrom(obj.getClass())) {
                throw new BadFileException(file, "file did not contain a " + type.getName());
            }
            
            @SuppressWarnings("unchecked")
            // no really, we just checked this above
            T castObj = (T)obj;
            
            return castObj;
            
        // can't recover from this one, need to bail hard
        } catch(StackOverflowError ex) {
            throw new Error("stack overflow, consider increasing -Xss", ex);
        
        // the caller can decide how to recover from this, so throw a friendlier exception
        } catch (ClassNotFoundException | InvalidClassException | OptionalDataException ex) {
            throw new BadFileException(file, "Osprey version doesn't match file", ex);
        } catch (IOException ex) {
            throw new BadFileException(file, "file is unreadable or corrupt", ex);
        }
    }
    
    
    
    public static void writeObject( Object outObj, String outFile ){
        try{
                FileOutputStream fout = new FileOutputStream(outFile);
                ObjectOutputStream out = new ObjectOutputStream(fout);
                out.writeObject(outObj);
                out.close();
                
        } catch (StackOverflowError ex) {
            //For objects with lots of links, might need to raise stack size, 
            //since writing is recursive
              throw new Error("Stack overflow in deepCopy.  Consider increasing -Xss");
        } catch (Exception ex) {
            throw new Error("can't write object", ex);
        }
    }
    
    
    public static void write(Object obj, String path)
    throws CantWriteException {
        write(obj, new File(path));
    }
    
    public static void write(Object obj, File file)
    throws CantWriteException {
        
        // writing is much less fraught with peril, and errors are much less likely to happen
        // so just pass along the errors
        
        try (ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(file))) {
            
            // try to write the object and see what happens
            out.writeObject(obj);
        
        // if this happens, bail hard and tell a programmer (hopefully it won't happen in production)
        } catch (InvalidClassException | NotSerializableException ex) {
            throw new Error("Can't write file data, classes not configured for serialization. This is a bug", ex);
        
        // sadly, unrecoverable, so bail hard
        } catch (StackOverflowError ex) {
            throw new Error("stack overflow, consider increasing -Xss", ex);
        
        // nothing we can do about the rest of the IOExceptions, just pass them up
        } catch (IOException ex) {
            throw new CantWriteException(file, ex);
        }
    }
    
    
    public static <T> T readOrMake(File file, Class<T> type, String name, Factory<T,Void> factory) {
    	return readOrMake(file, type, name, (thing) -> true, factory);
    }
    
    public static <T> T readOrMake(File file, Class<T> type, String name, Validator<T> validator, Factory<T,Void> factory) {
        
        // try to read from the cache
        try {
            
            T thing = read(file, type);
            if (thing != null) {
                System.out.println("read " + name + " from file: " + file.getAbsolutePath());
                
                // make sure it's valid
                if (validator.isValid(thing)) {
                	return thing;
                }
                
                System.out.println("WARNING: " + name + " from file is invalid, will create new one");
            }
            
        } catch (BadFileException ex) {
            ex.printStackTrace(System.out);
            System.out.println("WARNING: can't read " + name + ", will create new one");
        }
        
        // make the thing
        T thing = factory.make(null);
        
        // try to write to the cache
        try {
            
            write(thing, file);
            System.out.println("wrote " + name + " to file: " + file.getAbsolutePath());
            
        } catch (CantWriteException ex) {
            ex.printStackTrace(System.out);
            System.out.println("WARNING: can't write " + name + ", will have to be created again next time");
        }
        
        return thing;
    }
    
    
    
    	//Function adapted from: http://www.javaworld.com/javaworld/javatips/jw-javatip76.html?page=2
	//Java Tip 76: An alternative to the deep copy technique
	//Author: Dave Miller
	static public Object deepCopy(Object oldObj) {
          ObjectOutputStream oos = null;
          ObjectInputStream ois = null;
          try
          {
             ByteArrayOutputStream bos = 
                   new ByteArrayOutputStream(); // A
             oos = new ObjectOutputStream(bos); // B
             // serialize and pass the object
             oos.writeObject(oldObj);   // C
             oos.flush();               // D
             ByteArrayInputStream bin = 
                   new ByteArrayInputStream(bos.toByteArray()); // E
             ois = new ObjectInputStream(bin);                  // F
             
             // return the new object
             Object ans = ois.readObject(); // G
             
             oos.close();
             ois.close();
             
             return ans;
             
          } catch(StackOverflowError ex){
              //For objects with lots of links, might need to raise stack size, 
              //since writing is recursive
              throw new Error("Stack overflow in deepCopy.  Consider increasing -Xss");
          } catch(Exception ex) {
              throw new Error("can't deep-copy object", ex);
          }
       }
    
	

	public static File[] getFilesInDir( String path ) {

		File[] listOfFiles = null;

		File loc = new File(path);
		if( loc.isDirectory() )
			listOfFiles = loc.listFiles();
			
		return listOfFiles;
	}
	
	
	public static void delete( String path ) {
		try {
			File loc = new File(path);
			if( loc.exists() ) {
				if( loc.isDirectory() ) FileUtils.deleteDirectory(loc);
				else loc.delete();
			}
		} 
		catch (IOException ex) {
			throw new Error("can't delete file: " + path, ex);
		}
	}


	public static void makeDir( String path, boolean deleteExisting ) {
		File dir = new File(path);
		if( deleteExisting && dir.exists() ) delete(path);
		if( !dir.exists() ) dir.mkdirs();
	}
	
	
	public static String formatBigDecimal(BigDecimal x, int decimalPoints) {
		NumberFormat formatter = new DecimalFormat("0.0E0");
		formatter.setRoundingMode(RoundingMode.CEILING);
		formatter.setMinimumFractionDigits(decimalPoints);
		return formatter.format(x);
	}
}
