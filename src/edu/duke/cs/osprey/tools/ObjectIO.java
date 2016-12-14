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
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
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
