/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tools;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

/**
 *
 * @author mhall44
 */

//tools for storing objects, like an energy matrix

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
                    throw new RuntimeException("ERROR: Failed to read object from file "+fileName);
        }

        return inObj;
    }
    
    
    
    public static void writeObject( Object outObj, String outFile ){
        try{
                FileOutputStream fout = new FileOutputStream(outFile);
                ObjectOutputStream out = new ObjectOutputStream(fout);
                out.writeObject(outObj);
                out.close();
        }
        catch (Exception e){
                System.out.println(e.toString());
                System.out.println("ERROR: An exception occurred while writing object file");
                System.exit(0);
        }
    }
    
}
