package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.util.ArrayList;

import edu.duke.cs.osprey.structure.Residue;

@SuppressWarnings("serial")
public class Strand implements Serializable {
	
	public static final int PROTEIN = 0;
	public static final int LIGAND = 1;
	public static final int COMPLEX = 2;
	
	private int strand = -1;
	private int numFlexible = -1;
	private int begin = -1;
	private int end = -1;
	
	public static String getStrandString( int strand ) {
		switch ( strand ) {
        
		case Strand.COMPLEX:
        	return "complex";
        	
        case Strand.PROTEIN:
        	return "protein";
        
        case Strand.LIGAND:
        	return "ligand";
        
        default:
        	throw new RuntimeException("Error: invalid strand ID " + strand);
        }
		
	}
	
	public Strand( int strand, int numFlexible, ArrayList<String> limits ) {
		this.strand = strand;
		this.numFlexible = numFlexible;
		this.begin = Integer.parseInt(limits.get(0));
		this.end = Integer.parseInt(limits.get(1));
	}
	
	public int getStrandID() {
		return strand;
	}
	
	public int getNumFlexible() {
		return numFlexible;
	}
	
	public int getStrandBegin() {
		return begin;
	}
	
	public int getStrandEnd() {
		return end;
	}
	
	public boolean contains( Residue res ) {
		int pdbResNumber  = Integer.parseInt( res.getPDBResNumber() );
		return contains( pdbResNumber );
	}
	
	public boolean contains( int pdbResNumber ) {
		return pdbResNumber >= begin && pdbResNumber <= end ? true : false;
	}
}
