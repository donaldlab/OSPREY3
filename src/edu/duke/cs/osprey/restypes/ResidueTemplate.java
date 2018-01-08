/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class ResidueTemplate implements Serializable {
    //This class defines a residue type, like ARG or H2O or whatever
    //structures are composed entirely of residues. 
    //We need to identify the type of each residue whenever we want to move it it a residue-specific way
    //(e.g. to define sidechain dihedrals and rotamers),
    //or whenever we want to evaluate its energy (except by an ab initio method that takes only element 
    //types and coordinates)
    //We will also need these templates to mutate from one type to another

    //any of the information below can be null (i.e., not provided), but if so 
    //then attempts to call it will produce an error

    private static final long serialVersionUID = 4568917747972477569L;

    public String name;//e.g. ARG

    //"standard" residue
    //info from here will mostly be copied into any residue in the structure to which this template is assigned
    //(charges, standardized atom names, bonds, etc.)
    public Residue templateRes;

    public InterResBondingTemplate interResBonding;

    //dihedral information
    public int numDihedrals = 0;
    public int dihedral4Atoms[][];//for each dihedral, list of 4 atoms defining it
    //these are indices in all our atom-wise arrays
    public ArrayList<ArrayList<Integer>> dihedralMovingAtoms;//list of atoms that move for each dihedral


    // Rotameric information for this residue template. ResidueTemplate supports both backbone dependent and backbone independent rotamer libraries.
    // PGC 2015: If backbone dependent rotamer libraries are used, the number of rotamers is dependent on the backbone dihedrals.
    // Dimensions: (# phi increments)x(# psi increments)
    public int numRotamers[][] = null;
    //For each rotamer, gives the ideal values of all dihedrals
    //note: this is just the standard set of ideal rotamers, we can modify this as desired
    double rotamericDihedrals[][][][];// (# phi increments)x(# psi increments) x numRotamers x numDihedrals array.
    // Number of "bins" for Phi and Psi angles for backbone dependent rotamer libraries. For example, in the Dunbrack library there are 37 bins for each dihedral.
    // We assume that the number of bins will always be the same for phi or psi.
    public int numberOfPhiPsiBins = -1;
    // Resolution of phi and psi backbone dihedrals in the backbone dependent rotamer library. For example, in the Dunbrack rotamer library the resolution is 10. 
    public double phiPsiResolution = -1;

    public String CAEquivalent = null;//The name of the atom that the sidechain branches off of (e.g. CA in an amino acid)
    //if null we can't mutate
    
    public static ResidueTemplate makeFromResidueConfs(Residue ... residues) {
    	return makeFromResidueConfs(Arrays.asList(residues));
    }
    
    public static ResidueTemplate makeFromResidueConfs(List<Residue> residues) {
    	
    	// copy the first residue conformation to embed into the template
    	Residue firstRes = residues.get(0);
    	ResidueTemplate oldTemplate = firstRes.template;
    	Residue templateRes = new Residue(firstRes);
    	templateRes.copyIntraBondsFrom(firstRes);
        ResidueTemplate newTemplate = new ResidueTemplate(templateRes, oldTemplate.name, oldTemplate.interResBonding, oldTemplate.CAEquivalent);
    	
    	// copy template info
		newTemplate.setNumberOfPhiPsiBins(oldTemplate.numberOfPhiPsiBins);
		newTemplate.initializeRotamerArrays();
		newTemplate.dihedral4Atoms = oldTemplate.dihedral4Atoms;
		newTemplate.numDihedrals = oldTemplate.numDihedrals;
		newTemplate.setNumRotamers(residues.size());
		
		// measure the dihedrals for all conformations
		double[][] dihedrals = new double[residues.size()][newTemplate.numDihedrals];
		for (int i=0; i<residues.size(); i++) {
			Residue res = residues.get(i);
			assert (res.template.numDihedrals == oldTemplate.numDihedrals);
			for (int j=0; j<res.getNumDihedrals(); j++) {
				dihedrals[i][j] = res.getDihedralAngle(j);
			}
		}
		newTemplate.setRotamericDihedrals(dihedrals);
		newTemplate.computeDihedralMovingAtoms();
		
		// TODO: the second conformation doesn't necessarily differ from the first only by dihedrals
		// if these are eg crystallographic alternates, they could differ in eg bond angles too
		// ideally, the template should allow arbitrary conformations for each rotamer
		// but that's a much bigger refactor than we have to do today =)
		// for now, we'll get the i=0 rotamer conformation right, and the i>0 conformations will be approximations
		
		return newTemplate;
    }

    public ResidueTemplate(Residue res, String name, InterResBondingTemplate templ, String CAEquivalent){
        //initializes only with info from a template file, +interResBonding (and CAEquivalent if want to mutate)
        //res won't even have coordinates yet
        templateRes = res;
        this.name = name;
        this.CAEquivalent = CAEquivalent;
        interResBonding = templ;
    }

    //information on dihedrals
    public int[] getDihedralDefiningAtoms(int dihedralNum){
        //numbers (indices among the atoms in this template) of the 4 atoms defining the dihedral
        return dihedral4Atoms[dihedralNum];
    }

    public ArrayList<Integer> getDihedralRotatedAtoms(int dihedralNum){
        //the atoms that actually move (including the 4th of the dihedral-defining atoms)
        return dihedralMovingAtoms.get(dihedralNum);
    }

    public void computeDihedralMovingAtoms(){
        //compute what atoms move when a dihedral is changed
        //Compute from dihedral4Atoms, and write answer in dihedralMovingAtoms
        //we assume atom 4 (from the dihedral4Atoms) moves and the others are fixed; 
        //thus the moving atoms consist of the atoms bonded to atom 3 (excluding atom 2)
        //and anything beyond them (compute recursively).  If any of these moving atoms bond back to atom 2 
        //or atom 1 or to another residue,
        //then the dihedral can't move freely, and we return an error.  
        dihedralMovingAtoms = new ArrayList<>();

        for(int dihedNum=0; dihedNum<numDihedrals; dihedNum++){

            ArrayList<Integer> movingAtoms = new ArrayList<>();//indices of atoms that move (in templateRes.atoms)
            Atom atom2 = templateRes.atoms.get(dihedral4Atoms[dihedNum][1]);
            Atom atom3 = templateRes.atoms.get(dihedral4Atoms[dihedNum][2]);

            for(Atom atomBeyond : atom3.bonds){
                if(atomBeyond!=atom2){//atom 2 of dihedral doesn't move

                    movingAtomsDFS( atomBeyond, movingAtoms, dihedral4Atoms[dihedNum]);
                    //search recursively for other moving atoms
                }
            }

            dihedralMovingAtoms.add(movingAtoms);
        }
    }


    private void movingAtomsDFS(Atom atom, ArrayList<Integer> movingAtoms, int[] dihedAtoms){
        //Given an atom that moves when a dihedral is changed,
        //add it and any atoms distal to it to movingAtoms
        //dihedAtoms are the indices of the 4 atoms defining the dihedral
        if(atom==null)//atom is in another residue...this is an error
            throw new RuntimeException("ERROR: Trying to define non-free dihedral in rotamer library (bonds to another residue)");

        int atomIndex = templateRes.getAtomIndexByName(atom.name);

        if(atomIndex==dihedAtoms[0] || atomIndex==dihedAtoms[1])//there is a cycle within the residue
            throw new RuntimeException("ERROR: Trying to define non-free dihedral in rotamer library (dihedral is in a ring)");

        if(atomIndex == dihedAtoms[2])//dihedAtoms[2] doesn't move, though it is bonded to moving atoms
            return;

        //check that we haven't already counted atom as moving (since cycles are allowed within the set of moving atoms)
        if(movingAtoms.contains(atomIndex))
            return;

        //if we get here, atom is a new moving atom
        //add it to the list of movingAtoms, and search for atoms distal to it
        movingAtoms.add(atomIndex);
        for(Atom distalAtom : atom.bonds)
            movingAtomsDFS(distalAtom,movingAtoms,dihedAtoms);
    }
    /**
     * Returns the dihedral angle dihedralNum, for rotamer rotNum, for the backbone closest to the values phi and psi.
     * If the rotamer library is backbone independent, then phi and psi are ignored.
     * @param phi The backbone phi angle 
     * @param psi The backbone psi angle
     * @param rotNum The rotamer's index number
     * @param dihedralNum The dihedral's index number.
     * @return the angle value for the desired dihedral.
     */
    public double getRotamericDihedrals(double phi, double psi, int rotNum, int dihedralNum) {
    	return getRotamericDihedrals(phi, psi, rotNum)[dihedralNum];
    }
    
    public double getRotamericDihedrals(int rotNum, int dihedralNum) {
        return getRotamericDihedrals(0,0,rotNum,dihedralNum);
    }
    
    public double[] getRotamericDihedrals(int rotNum) {
    	return getRotamericDihedrals(0, 0, rotNum);
    }
    
    public double[] getRotamericDihedrals(double phi, double psi, int rotNum) {
    	
        if(Double.isNaN(phi) || Double.isNaN(psi)){//dihedrals not defined for this residue
            if(numberOfPhiPsiBins > 1)
                throw new RuntimeException("ERROR: Can't use Dunbrack library on residues w/o phi/psi defined");

            return rotamericDihedrals[0][0][rotNum];
        }

        // Under the dunbrack rotamer library, backbone dependent rotamers have a resolution of 10 degrees, 
        //    while in backbone-independent rotamer libraries they have a resolution of 360.    		
        //  Therefore, we round to the closest "bin" and add numberOfPhiPsiBins/2 (to make all numbers positive)
        int phiBin = (int)(Math.round(phi/this.phiPsiResolution)) + this.numberOfPhiPsiBins/2;
        int psiBin = (int)(Math.round(psi/this.phiPsiResolution)) + this.numberOfPhiPsiBins/2;
        return rotamericDihedrals[phiBin][psiBin][rotNum];
    }
    
    /**
     * Returns the number of rotamers for a backbone-dependent or backbone independent rotamer library.
     * @param phi The backbone phi angle, ignored if backbone independent rotamer library.
     * @param psi The backbone psi angle, ignored if backbone independent rotamer library.
     * @return The number of rotamers.
     */

    public int getNumRotamers(double phi, double psi){

        if(numberOfPhiPsiBins==-1)
            throw new RuntimeException("ERROR: Rotamers not set up for residue "+name);
        
        if(Double.isNaN(phi) || Double.isNaN(psi)){//dihedrals not defined for this residue
            if(numberOfPhiPsiBins > 1)
                throw new RuntimeException("ERROR: Can't use Dunbrack library on residues w/o phi/psi defined");

            return numRotamers[0][0];
        }


        // Under the dunbrack rotamer library, backbone dependent rotamers have a resolution of 10 degrees, 
        //    while in backbone-independent rotamer libraries they have a resolution of 360.    		
        //  Therefore, we round to the closest "bin" and add numberOfPhiPsiBins/2 (to make all numbers positive)
        int phiBin = (int)((Math.round(phi/this.phiPsiResolution))) + this.numberOfPhiPsiBins/2;
        int psiBin = (int)((Math.round(psi/this.phiPsiResolution))) + this.numberOfPhiPsiBins/2;
        return this.numRotamers[phiBin][psiBin];
    }
    
    public int getNumRotamers() {
        return getNumRotamers(0,0);
    }
    
    /**
     * Number of "bins" for Phi/Psi angles for backbone dependent rotamer libraries. For example, in the Dunbrack library there are 37 bins.
     * @param numberOfPsiBins
     */
    public void setNumberOfPhiPsiBins(int numberOfPhiPsiBins) {
        this.numberOfPhiPsiBins = numberOfPhiPsiBins;
    }
    /**
     * Resolution of backbone "bins" in the rotamer library. For example in a backbone independent rotamer library this should be "360". 
     * For the Dunbrack rotamer library this is 10. Set in the corresponding parser of the rotamer library. 
     * @param phiPsiResolution.
     */
    public void setRLphiPsiResolution(double phiPsiResolution) {
        this.phiPsiResolution = phiPsiResolution;
    }

    /**
     * Initialize the numRotamers and rotamericDihedrals arrays to the number of bins according to each rotamer library. 
     */
    public void initializeRotamerArrays(){
        assert (numberOfPhiPsiBins > 0);
        this.numRotamers= new int [numberOfPhiPsiBins][];
        this.rotamericDihedrals = new double[numberOfPhiPsiBins][][][];
        for (int i = 0; i < numberOfPhiPsiBins; i++){
            this.numRotamers[i] = new int [numberOfPhiPsiBins];
            this.rotamericDihedrals[i] = new double[numberOfPhiPsiBins][][];
        }   	

    }

    /**
     * PGC 2015: 
     * Sets the number of rotamers for a residue template; phi and psi are only 
     * used when using a backbone dependent rotamer library.
     * @param numRotamers Number of rotamers
     * @param phiBin Bin in the rotamer array where the rotameric dihedrals will be stored. 
     * @param psiBin Bin in the rotamer array where the rotameric dihedrals will be stored. 
     */

    public void setNumRotamers(int numRotamers, int phiBin, int psiBin){

        this.numRotamers[phiBin][psiBin] = numRotamers;

    }
    
    public void setNumRotamers(int numRotamers){
        setNumRotamers(numRotamers,0,0);
    }
    
    /**
     * PGC 2015:
     * Set the number of dihedrals for this Residue type 
     */
    public void setNumDihedrals(int numDihedrals){
        this.numDihedrals = numDihedrals;
    }
    /**
     * PGC 2015: 
     * Sets the number of rotamers for a residue template; phi and psi are only 
     * used when using a backbone dependent rotamer library.
     * @param numRotamers Number of rotamers
     * @param phiBin Bin in the rotamer array where the rotameric dihedrals will be stored. 
     * @param psiBin Bin in the rotamer array where the rotameric dihedrals will be stored. 
     */    
    public void setRotamericDihedrals(double newRotamericDihedrals[][], int phiBin, int psiBin){
        this.rotamericDihedrals[phiBin][psiBin] = newRotamericDihedrals;
    }
    
    public void setRotamericDihedrals(double newRotamericDihedrals[][]){
        setRotamericDihedrals(newRotamericDihedrals,0,0);
     }

    @Override
    public String toString() {
        // there are tons of templates with the same name, so add the hash id too
        return name + ":" + System.identityHashCode(this);
    }
}
