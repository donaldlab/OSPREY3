/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy.forcefield;

import java.io.Serializable;
import java.util.Iterator;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SimpleConfSpace.Builder;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 * Options for configuring forcefields for energy calculation.
 * @author mhall44
 */
public class ForcefieldParams implements Serializable {
    
    private static final long serialVersionUID = 3124964506851762586L;
    
    public static final double coulombConstant = 332.0;
    public static final double solvCutoff = 9.0;
    public static final double solvCutoff2 = solvCutoff*solvCutoff;

    public static boolean printWarnings = true;

    final int atomTypeX = -2; //the atom type number for the X wildcard atom type
    private final int noMatchInt = 9999;

    String[] atomTypeNames = null;
    double[] atomAtomicMasses = null;
    int[]   bondAtomType1 = null;
    int[]   bondAtomType2 = null;
    double[] bondHFC = null; // Harmonic Force Constant
    double[] bondEBL = null; // Equilibrium Bond Length
    int[]		angleAtomType1 = null;
    int[]		angleAtomType2 = null;
    int[]		angleAtomType3 = null;
    double[] angleHFC = null; // Harmonic Force Constant
    double[] angleEBA = null; // Equilibrium Bond Angle
    int numGeneralDihedParams = 0; //the number of generic dihedral parameters (that have wildcard X atoms)
    int[]		dihedAtomType1 = null;
    int[]		dihedAtomType2 = null;
    int[]		dihedAtomType3 = null;
    int[]		dihedAtomType4 = null;
    // Dihedral: (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
    double[] dihedTerm1 = null; // (PK/IDIVF)
    int[] dihedPN = null; // Periodicity
    double[] dihedPhase = null; // Phase
    int[]		impDihedAtomType1 = null;
    int[]		impDihedAtomType2 = null;
    int[]		impDihedAtomType3 = null;
    int[]		impDihedAtomType4 = null;
    // Imprper Dihedral: PK * (1 + cos(PN*phi - PHASE))
    double[] impDihedTerm1 = null; // PK
    int[] impDihedPN = null; // Periodicity
    double[] impDihedPhase = null; // Phase
    int[]		vdwAtomType1 = null;
    double[] vdwR = null;
    double[] vdwE = null;
    int[][] equivAtoms = null;
    // For two atoms i & j
    //  B(ij) = 2 * (Ri + Rj)^6 * ei * ej
    //  A(ij) = (Ri + Rj)^12 * ei * ej
    // Those that are 1-4 separated are scaled 


    
    //what forcefield are these parameters for?
    public Forcefield forcefld = Forcefield.AMBER;
    
    //The solvation parameters object
    public EEF1 eef1parms = null;
    
    // default values copied from resources/config/defalts.cfg
    public double vdwMultiplier = 0.95;
    public double solvScale = 0.5; //the scale factor for the solvation energies
    public double dielectric = 6;	
    public boolean distDepDielect = true;
    public boolean hElect = true;
    public boolean hVDW = true;
    /** for new code, use {@link Builder#shellDist} instead **/
    @Deprecated
    public double shellDistCutoff = Double.POSITIVE_INFINITY; //distance cutoff for interactions (angstroms)
    public SolvationForcefield solvationForcefield = SolvationForcefield.EEF1;
    
    public enum Forcefield {
        
        AMBER(
            "/config/parm96a.dat",
            "/config/all_amino94.in",
            "/config/all_aminont94.in",
            "/config/all_aminoct94.in",
            "/config/all_nuc94_and_gr.in",
            0.5,
            1.0,
            1.0/1.2
        ),
        CHARMM22(
            "/config/parmcharm22.dat",
            "/config/all_amino_charmm22.txt",
            "/config/all_amino_charmm22_nt.txt",
            "/config/all_amino_charmm22_ct.txt",
            "/config/all_nuc_and_gr_charmm.in",
            Double.NaN,
            Double.NaN,
            Double.NaN
        ),
        CHARMM19NEUTRAL(
            "/config/parmcharm19.dat",
            "/config/all_amino_charmm19_neutral.in",
            "/config/all_amino_charmm19_neutral_nt.in",
            "/config/all_amino_charmm19_neutral_ct.in",
            "/config/all_nuc_and_gr_charmm.in",
            1.0,
            2.0,
            0.4
        ) {
        	@Override
        	public void modifyNBParams(Atom atom, AtomNeighbors.Type neighborType, NBParams nbparams) {
        		CHARMM19.modifyNBParams(atom, neighborType, nbparams);
        	}
        },
        CHARMM19(
            "/config/parmcharm19.dat",
            "/config/all_amino_charmm19.in",
            "/config/all_amino_charmm19_nt.in",
            "/config/all_amino_charmm19_ct.in",
            "/config/all_nuc_and_gr_charmm.in",
            1.0,
            2.0,
            0.4
        ) {
        	@Override
        	public void modifyNBParams(Atom atom, AtomNeighbors.Type neighborType, NBParams nbparams) {
        		
        		if (atom.isCarbon() && neighborType == AtomNeighbors.Type.BONDED14) {
					nbparams.epsilon = 0.1;
					nbparams.r = 1.9;
        		}
        	}
        };
        
        public final String paramsPath;
        public final String aaPath;
        public final String aaNTPath;
        public final String aaCTPath;
        public final String grPath;
        public final double Aij14Factor;
        public final double Bij14Factor;
        public final double coulombScaling;
        
        private Forcefield(String paramsPath, String aaPath, String aaNTPath, String aaCTPath, String grPath, double Aij14Factor, double Bij14Factor, double coulombScaling) {
            this.paramsPath = paramsPath;
            this.aaPath = aaPath;
            this.aaNTPath = aaNTPath;
            this.aaCTPath = aaCTPath;
            this.grPath = grPath;
            this.Aij14Factor = Aij14Factor;
            this.Bij14Factor = Bij14Factor;
            this.coulombScaling = coulombScaling;
        }
        
        public static Forcefield get(String name) {
        	return valueOf(name.toUpperCase());
        }
        
        public void modifyNBParams(Atom atom, AtomNeighbors.Type neighborType, NBParams nbparams) {
        	// by default, don't modify anything
        }
    }
    
    public static enum SolvationForcefield {
    	
        EEF1 {
        	@Override
        	public ResiduesInfo makeInfo(ForcefieldParams ffparams, Residues residues) {
        		return ffparams.eef1parms.new ResiduesInfo(residues);
        	}
        },
        PoissonBoltzmann {
        	@Override
        	public ResiduesInfo makeInfo(ForcefieldParams ffparams, Residues residues) {
        		// TODO: implmement PB solvation
        		return null;
        	}
        };
    	
    	public interface ResiduesInfo {
			int getNumPrecomputedPerAtomPair();
			double getResPairEnergy(Residue res1, Residue res2);
			void putPrecomputed(double[] out, int i, Residue res1, int atomIndex1, Residue res2, int atomIndex2, double scale);
    	}
    	
    	public abstract ResiduesInfo makeInfo(ForcefieldParams ffparams, Residues residues);
    }
    
    public ForcefieldParams() {
        init();
    }
    
    public ForcefieldParams(String frcefld) {
        this.forcefld = Forcefield.valueOf(frcefld.toUpperCase());
        init();
    }
    
    public ForcefieldParams(Forcefield frcefld) {
        this.forcefld = frcefld;
        init();
    }
    
    private void init() {
        // Read in AMBER forcefield parameters
        // parm96a.dat
        try {
            readParm96(forcefld);
        }
        catch (Exception ex) {
            throw new Error("can't read forcefield params", ex);
        }
                
        // Read in the EEF1 solvation parameters
        // TODO: lazy loading of EEF1 params?
        try {
            eef1parms = new EEF1();
            eef1parms.readEEF1parm();
        } catch (Exception ex) {
            throw new Error("can't read solvation params", ex);
        }
    }
    
    public ForcefieldParams(ForcefieldParams other) {
        this(other.forcefld);
        
        vdwMultiplier = other.vdwMultiplier;
        solvScale = other.solvScale;
        dielectric = other.dielectric;	
        distDepDielect = other.distDepDielect;
        hElect = other.hElect;
        hVDW = other.hVDW;
        shellDistCutoff = other.shellDistCutoff;
        solvationForcefield = other.solvationForcefield;
    }
    
    
    
    //************************************
	// This function reads the AMBER forcefield parameter file
	//  parm96a.dat
	// Although this function is meant to be relatively generic
	//  it is optimized for reading parm96a.dat. Slight changes
	//  will most likely be required to read other parameter
	//  files. Reading of other files should be done in other
	//  functions
	private void readParm96(Forcefield ff) throws Exception {
	
		Iterator<String> lines = FileTools.parseLines(FileTools.readResource(ff.paramsPath)).iterator();
			
		String curLine = null;
		int tmpInt = 0;
		
		final int initSize = 10; //the initial size of the arrays to store the data that is read

		// Skip over the first line of header info
		curLine = lines.next();
		
		// 1. Read atom names and atomic masses
		atomTypeNames = new String[initSize];
		atomAtomicMasses = new double[initSize];
		curLine = lines.next();
		tmpInt = 0; // temporary integer
		// Until we're at a blank line (or until we've read numAtomTypes)
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=atomTypeNames.length){ //double the array sizes
				atomTypeNames = doubleArraySize(atomTypeNames);
				atomAtomicMasses = doubleArraySize(atomAtomicMasses);
			}
			
			atomTypeNames[tmpInt] = StringParsing.getToken(curLine,1);  // snag atom name
			atomAtomicMasses[tmpInt] = (new Double(StringParsing.getToken(curLine,2))).doubleValue();
			tmpInt++;
			curLine = lines.next();
		}
		atomTypeNames = reduceArraySize(atomTypeNames,tmpInt);
		atomAtomicMasses = reduceArraySize(atomAtomicMasses,tmpInt);
		

		// Skip unknown line
		curLine = lines.next();

		// 2. Read Bonds
		bondAtomType1 = new int[initSize];
		bondAtomType2 = new int[initSize];
		bondHFC = new double[initSize];
		bondEBL = new double[initSize];
		curLine = lines.next();
		tmpInt = 0;
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=bondAtomType1.length){
				bondAtomType1 = doubleArraySize(bondAtomType1);
				bondAtomType2 = doubleArraySize(bondAtomType2);
				bondHFC = doubleArraySize(bondHFC);
				bondEBL = doubleArraySize(bondEBL);
			}
			
			//tmpStr = curLine.substring(0,5);
			bondAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine, 1));
			bondAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine, 2));
			bondHFC[tmpInt] = (new Double(getDashedToken(curLine,3))).doubleValue();
			bondEBL[tmpInt] = (new Double(getDashedToken(curLine,4))).doubleValue();
			tmpInt++;
			curLine = lines.next();
		}
		bondAtomType1 = reduceArraySize(bondAtomType1,tmpInt);
		bondAtomType2 = reduceArraySize(bondAtomType2,tmpInt);
		bondHFC = reduceArraySize(bondHFC,tmpInt);
		bondEBL = reduceArraySize(bondEBL,tmpInt);
		

		// 3. Read Angles
		angleAtomType1 = new int[initSize];
		angleAtomType2 = new int[initSize];
		angleAtomType3 = new int[initSize];
		angleHFC = new double[initSize];
		angleEBA = new double[initSize];
		curLine = lines.next();
		tmpInt = 0;
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=angleAtomType1.length){
				angleAtomType1 = doubleArraySize(angleAtomType1);
				angleAtomType2 = doubleArraySize(angleAtomType2);
				angleAtomType3 = doubleArraySize(angleAtomType3);
				angleHFC = doubleArraySize(angleHFC);
				angleEBA = doubleArraySize(angleEBA);
			}
			
			//tmpStr = curLine.substring(0,8);
			angleAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine,1));
			angleAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine,2));
			angleAtomType3[tmpInt] = atomTypeToInt(getDashedToken(curLine,3));
			angleHFC[tmpInt] = (new Double(getDashedToken(curLine,4))).doubleValue();
			angleEBA[tmpInt] = (new Double(getDashedToken(curLine,5))).doubleValue();
			tmpInt++;
			curLine = lines.next();
		}
		angleAtomType1 = reduceArraySize(angleAtomType1,tmpInt);
		angleAtomType2 = reduceArraySize(angleAtomType2,tmpInt);
		angleAtomType3 = reduceArraySize(angleAtomType3,tmpInt);
		angleHFC = reduceArraySize(angleHFC,tmpInt);
		angleEBA = reduceArraySize(angleEBA,tmpInt);
		
		
		// 4. Read Dihedrals
		numGeneralDihedParams = 0;
		dihedAtomType1 = new int[initSize];
		dihedAtomType2 = new int[initSize];
		dihedAtomType3 = new int[initSize];
		dihedAtomType4 = new int[initSize];
		dihedTerm1 = new double[initSize];
		dihedPhase = new double[initSize];
		dihedPN = new int[initSize];
		curLine = lines.next();
		tmpInt = 0;
		double tmpFlt = 0.0f;
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=dihedAtomType1.length){
				dihedAtomType1 = doubleArraySize(dihedAtomType1);
				dihedAtomType2 = doubleArraySize(dihedAtomType2);
				dihedAtomType3 = doubleArraySize(dihedAtomType3);
				dihedAtomType4 = doubleArraySize(dihedAtomType4);
				dihedTerm1 = doubleArraySize(dihedTerm1);
				dihedPhase = doubleArraySize(dihedPhase);
				dihedPN = doubleArraySize(dihedPN);
			}
			
			//tmpStr = curLine.substring(0,11);
			dihedAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine,1));
			dihedAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine,2));
			dihedAtomType3[tmpInt] = atomTypeToInt(getDashedToken(curLine,3));
			dihedAtomType4[tmpInt] = atomTypeToInt(getDashedToken(curLine,4));
			
			if ( dihedAtomType1[tmpInt]==atomTypeX || dihedAtomType2[tmpInt]==atomTypeX || dihedAtomType3[tmpInt]==atomTypeX || dihedAtomType4[tmpInt]==atomTypeX ) //at least one of the atoms is a wildcard
				numGeneralDihedParams++;
			
			tmpFlt = (new Double(getDashedToken(curLine,5))).doubleValue();
			dihedTerm1[tmpInt] = (new Double(getDashedToken(curLine,6))).doubleValue() / tmpFlt;
			dihedPhase[tmpInt] = (new Double(getDashedToken(curLine,7))).doubleValue();
			dihedPN[tmpInt] = (new Double(getDashedToken(curLine,8))).intValue();
			// If dihedPN is negative then there are one or more additional terms
			//  nothing fancy needs to be done because they will all be read in anyway
			//  but we do need to correct the sign
			if (dihedPN[tmpInt] < 0)
				dihedPN[tmpInt] = -dihedPN[tmpInt];
			tmpInt++;
			curLine = lines.next();
		}
		dihedAtomType1 = reduceArraySize(dihedAtomType1,tmpInt);
		dihedAtomType2 = reduceArraySize(dihedAtomType2,tmpInt);
		dihedAtomType3 = reduceArraySize(dihedAtomType3,tmpInt);
		dihedAtomType4 = reduceArraySize(dihedAtomType4,tmpInt);
		dihedTerm1 = reduceArraySize(dihedTerm1,tmpInt);
		dihedPhase = reduceArraySize(dihedPhase,tmpInt);
		dihedPN = reduceArraySize(dihedPN,tmpInt);
		

		// 5. Read Improper Dihedrals
		impDihedAtomType1 = new int[initSize];
		impDihedAtomType2 = new int[initSize];
		impDihedAtomType3 = new int[initSize];
		impDihedAtomType4 = new int[initSize];
		impDihedTerm1 = new double[initSize];
		impDihedPhase = new double[initSize];
		impDihedPN = new int[initSize];
		curLine = lines.next();
		tmpInt = 0;
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=impDihedAtomType1.length){
				impDihedAtomType1 = doubleArraySize(impDihedAtomType1);
				impDihedAtomType2 = doubleArraySize(impDihedAtomType2);
				impDihedAtomType3 = doubleArraySize(impDihedAtomType3);
				impDihedAtomType4 = doubleArraySize(impDihedAtomType4);
				impDihedTerm1 = doubleArraySize(impDihedTerm1);
				impDihedPhase = doubleArraySize(impDihedPhase);
				impDihedPN = doubleArraySize(impDihedPN);
			}
			
			//tmpStr = curLine.substring(0,11);
			impDihedAtomType1[tmpInt] = atomTypeToInt(getDashedToken(curLine,1));
			impDihedAtomType2[tmpInt] = atomTypeToInt(getDashedToken(curLine,2));
			impDihedAtomType3[tmpInt] = atomTypeToInt(getDashedToken(curLine,3));
			impDihedAtomType4[tmpInt] = atomTypeToInt(getDashedToken(curLine,4));
			impDihedTerm1[tmpInt] = (new Double(getDashedToken(curLine,5))).doubleValue();
			impDihedPhase[tmpInt] = (new Double(getDashedToken(curLine,6))).doubleValue();
			impDihedPN[tmpInt] = (new Double(getDashedToken(curLine,7))).intValue();
			tmpInt++;
			curLine = lines.next();
		}
		impDihedAtomType1 = reduceArraySize(impDihedAtomType1,tmpInt);
		impDihedAtomType2 = reduceArraySize(impDihedAtomType2,tmpInt);
		impDihedAtomType3 = reduceArraySize(impDihedAtomType3,tmpInt);
		impDihedAtomType4 = reduceArraySize(impDihedAtomType4,tmpInt);
		impDihedTerm1 = reduceArraySize(impDihedTerm1,tmpInt);
		impDihedPhase = reduceArraySize(impDihedPhase,tmpInt);
		impDihedPN = reduceArraySize(impDihedPN,tmpInt);
		

		// Skip 2 lines (we might also be able to go until the keyword MOD4
		curLine = lines.next();
		curLine = lines.next();

		// Read the equivalence lines
		// The first atomnum in equivAtoms is the main atom and numbers with index 1..n are
		//  equivalent to the atom in index 0.
		equivAtoms = new int[initSize][];
		curLine = lines.next();
		tmpInt = 0;
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=equivAtoms.length){
				equivAtoms = doubleArraySize(equivAtoms);
			}
			
			int numEquivAtoms = (new StringTokenizer(curLine," ,;\t\n\r\f")).countTokens();
			equivAtoms[tmpInt] = new int[numEquivAtoms];
			
			for(int q=0;q<equivAtoms[tmpInt].length;q++)
				equivAtoms[tmpInt][q] = -noMatchInt;
			int tmpInt2=1;
			while (!StringParsing.getToken(curLine,tmpInt2).equalsIgnoreCase("")) { 
				equivAtoms[tmpInt][tmpInt2-1] = atomTypeToInt( StringParsing.getToken(curLine,tmpInt2) );
				tmpInt2++;
			}
			tmpInt++;
			curLine = lines.next();
		}
		equivAtoms = reduceArraySize(equivAtoms,tmpInt);
		
				
		// Skip a line (we might also be able to go until the keyword MOD4
		curLine = lines.next();

		// 6. Read vdw
		vdwAtomType1 = new int[initSize];
		vdwR = new double[initSize];
		vdwE = new double[initSize];
		curLine = lines.next();
		tmpInt = 0;
		while (!(StringParsing.getToken(curLine,1).equals(""))) {
			
			if (tmpInt>=vdwAtomType1.length){
				vdwAtomType1 = doubleArraySize(vdwAtomType1);
				vdwR = doubleArraySize(vdwR);
				vdwE = doubleArraySize(vdwE);
			}
			
			vdwAtomType1[tmpInt] = atomTypeToInt(StringParsing.getToken(curLine,1));
			vdwR[tmpInt] = (new Double(StringParsing.getToken(curLine,2))).doubleValue();
			vdwE[tmpInt] = (new Double(StringParsing.getToken(curLine,3))).doubleValue();
			tmpInt++;
			curLine = lines.next();
		}
		vdwAtomType1 = reduceArraySize(vdwAtomType1,tmpInt);
		vdwR = reduceArraySize(vdwR,tmpInt);
		vdwE = reduceArraySize(vdwE,tmpInt);		
			
	// DEBUG START * Good to keep for when the parameter file changes to
	//  make sure you're reading it correctly
/*	System.out.println("ATOM TYPES");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + atomTypeNames[q] + " mass:" + atomAtomicMasses[q]);
	}
	System.out.println("BONDS");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + bondAtomType1[q] + " " + bondAtomType2[q] + " HFC:" + 
			bondHFC[q] + " EBL:" + bondEBL[q]);
	}
	System.out.println("ANGLES");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + angleAtomType1[q] + " " + angleAtomType2[q] + " " + angleAtomType3[q]
		 + " HFC:" + angleHFC[q] + " EBA:" + angleEBA[q]);
	}
	System.out.println("DIHED");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + dihedAtomType1[q] + " " + dihedAtomType2[q] + " " + dihedAtomType3[q]
		 + " " + dihedAtomType4[q] + " term1:" + dihedTerm1[q] + " PN:" + dihedPN[q] + " Phase:" + dihedPhase[q]);
	}
	System.out.println("IMPROP DIHED");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + impDihedAtomType1[q] + " " + impDihedAtomType2[q] + " " + impDihedAtomType3[q]
		 + " " + impDihedAtomType4[q] + " term1:" + impDihedTerm1[q] + " PN:" + impDihedPN[q] + " Phase:" + impDihedPhase[q]);
	}
	System.out.println("VDW");
	for(int q=0;q<5;q++) {
		System.out.println(q + " " + vdwAtomType1[q] + " " + " R:" + vdwR[q] + " E:" + vdwE[q]);
	}
	// DEBUG END
*/
	
	}


	private String getDashedToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f-");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end StringParsing.getToken
	

	// This function returns the numeric atom type based on the string atom type
	// If atom type is 'x' then return atomTypeX which means it's a wildcard
	public int atomTypeToInt(String s) {
		s = s.trim();
		if (s.equalsIgnoreCase("x"))
			return atomTypeX;
		for(int q=0;q<atomTypeNames.length;q++) {
			if (atomTypeNames[q].equalsIgnoreCase(s))
				return q;
		}
		return -1;
	}
	
	
	// This function searches the bond constants for the atoms specified
	//  and returns the approprate constants 
	public boolean getStretchParameters(int atomType1, int atomType2,
		double forceConstant[],	double equilibriumDistance[]){

		int numGeneric = 5, tmpInt = 0;
		double tmpFC = 0.0f, tmpED = 0.0f;
		boolean matched = false;
		for(int q=0;q<bondAtomType1.length;q++) {
			if (match2Atoms(atomType1, atomType2, bondAtomType1[q], bondAtomType2[q])) {
				tmpInt = 0;
				matched = true;
				if (bondAtomType1[q] == atomTypeX)
					tmpInt++;
				if (bondAtomType2[q] == atomTypeX)
					tmpInt++;
				if (tmpInt < numGeneric) {
					numGeneric = tmpInt;
					tmpFC = bondHFC[q];
					tmpED = bondEBL[q];
				}
			}
		}
		
		if (matched) {
			forceConstant[0] = tmpFC;
			equilibriumDistance[0] = tmpED;
			return(true);
		}
		else {
			forceConstant[0] = 317;
      equilibriumDistance[0] = 1.522;
			System.out.println("Ambstretch DEFAULTING TO C-CT");
			return(false);
		}
	}
	
	// This function searches the angle constants for the atoms specified
	//  and returns the approprate constants
	public boolean getBendParameters(int atomType1, int atomType2, int atomType3, 
		double forceConstant[], double equilibriumAngle[]){

		if (atomType3 < atomType1){
			int temp = atomType3;
			atomType3 = atomType1;
			atomType1 = temp;
		}

		int numGeneric = 5, tmpInt = 0;
		double tmpFC = 0.0f, tmpEA = 0.0f;
		boolean matched = false;
		for(int q=0;q<angleAtomType1.length;q++) {
			if (match3Atoms(atomType1, atomType2, atomType3, angleAtomType1[q], 
				angleAtomType2[q], angleAtomType3[q])) {
				tmpInt = 0;
				matched = true;
				if (angleAtomType1[q] == atomTypeX)
					tmpInt++;
				if (angleAtomType2[q] == atomTypeX)
					tmpInt++;
				if (angleAtomType3[q] == atomTypeX)
					tmpInt++;
				if (tmpInt < numGeneric) {
					numGeneric = tmpInt;
					tmpFC = angleHFC[q];
					tmpEA = angleEBA[q];
				}
			}
		}
		
		if (matched) {
			forceConstant[0] = tmpFC;
			equilibriumAngle[0] = tmpEA;
			return(true);
		}
		else {
			forceConstant[0] = 63.000000;
			equilibriumAngle[0] = 113.100000;
			System.out.println( "AmberBend: Could not find correct angle, defaulting to CC-CT-CT");
			return( false );
		}
	}

	// This function searches the dihedral constants for the atoms specified
	//  and returns the approprate constants
	public boolean getTorsionParameters(int atomType1, int atomType2, int atomType3, 
		int atomType4, double forceConstant[], double equilibriumAngle[], int terms[], 
		int multiplicity[]){

		if (atomType3 < atomType2){
			int temp = atomType3;
			atomType3 = atomType2;
			atomType2 = temp;
			temp = atomType4;
			atomType4 = atomType1;
			atomType1 = temp;
		}

		// First search through generic torsions (those terms containing wildcards)
		boolean matched = false;
		for(int q=0;q<numGeneralDihedParams;q++) {
			if (match2Atoms(atomType2, atomType3, dihedAtomType2[q], dihedAtomType3[q])) {
				matched = true;
				forceConstant[0] = dihedTerm1[q];
				equilibriumAngle[0] = dihedPhase[q];
				terms[0] = dihedPN[q];
				multiplicity[0] = 0;
			}
		}
		
		// According to the original paper "any specific parameter, such as OS-CH-CH-OS,
		//  overrides any general parameter."
		int forceCounter = 0, eqCounter = 0, multCounter = 0, termCounter = 0;
		for(int q=numGeneralDihedParams;q<dihedAtomType1.length;q++) {
			if (match4Atoms(atomType1, atomType2, atomType3, atomType4, dihedAtomType1[q],
						dihedAtomType2[q], dihedAtomType3[q], dihedAtomType4[q])) {
				matched = true;
				forceConstant[forceCounter++] = dihedTerm1[q];
				equilibriumAngle[eqCounter++] = dihedPhase[q];
				terms[termCounter++] = dihedPN[q];
				multiplicity[0] = multCounter++;
			}
		}
		
		if (matched) {
			return(true);
		}
		else {
			System.out.println("AmberDihed: Could not find correct torsion");
			return( false );
		}
	}
	
	// This function attempts to match two atoms for a bond. An atom type
	//  of atomTypeX is a generic term and can match anything
	private boolean match2Atoms(int atType1, int atType2, int known1, int known2) {
	
		if ((atType1 == known1) && (atType2 == known2))
			return(true);
		if ((atType1 == known2) && (atType2 == known1))
			return(true);
		if ((atType1 == known1) && (known2 == atomTypeX))
			return(true);
		if ((atType1 == known2) && (known1 == atomTypeX))
			return(true);
		if ((atType2 == known1) && (known2 == atomTypeX))
			return(true);
		if ((atType2 == known2) && (known1 == atomTypeX))
			return(true);
		return(false);
	}
	
	// This function attempts to match three atoms for an angle. An atom
	//  type of atomTypeX is a generic term and can match anything
	private boolean match3Atoms(int atType1, int atType2, int atType3,
		int known1, int known2, int known3) {
		
		if ((atType1 == known1) && (atType2 == known2) && (atType3 == known3))
			return(true);
		if ((atType1 == known3) && (atType2 == known2) && (atType3 == known1))
			return(true);
		if ((known1 == atomTypeX) && (atType2 == known2) && (atType3 == known3))
			return(true);
		if ((known1 == atomTypeX) && (atType2 == known2) && (atType1 == known3))
			return(true);
		if ((atType1 == known1) && (atType2 == known2) && (known3 == atomTypeX))
			return(true);
		if ((atType3 == known1) && (atType2 == known2) && (known3 == atomTypeX))
			return(true);
		if ((atType1 == known1) && (known2 == atomTypeX) && (atType3 == known3))
			return(true);
		if ((atType3 == known1) && (known2 == atomTypeX) && (atType1 == known3))
			return(true);
		return(false);
	}

	// This function attempts to match four atoms for a dihedral (no generic atoms)
	private boolean match4Atoms(int atType1, int atType2, int atType3, int atType4,
		int known1, int known2, int known3, int known4) {
	
		if ((atType1 == known1) && (atType2 == known2) &&
			(atType3 == known3) && (atType4 == known4))
			return(true);
		else if ((atType1 == known4) && (atType2 == known3) &&
			(atType3 == known2) && (atType4 == known1))
					return(true);

		return(false);
	}

	
	// This function returns the equivalent class for the given atomtype
	private int getEquivalentType(int atomType) {

		for(int i=0;i<equivAtoms.length;i++) {
			for(int j=1;j<equivAtoms[i].length;j++) {
				if (atomType == equivAtoms[i][j])
					return(equivAtoms[i][0]);
			}
		}
		return(-1);
	}


	public static class NBParams {
		public double r;
		public double epsilon;
	}
	
	// This function returns the r and epsilon paramters for a given atom type
	public boolean getNonBondedParameters(int atomType, NBParams out) {
		
		for(int q=0;q<vdwAtomType1.length;q++) {
			if (vdwAtomType1[q]==atomType) {
				out.r=vdwR[q];
				out.epsilon=vdwE[q];
				return (true);
			}
		}

		// Check for equivalent atoms
		int equivType = getEquivalentType(atomType);
		for(int q=0;q<vdwAtomType1.length;q++) {
			if (vdwAtomType1[q]==equivType) {
				out.r=vdwR[q];
				out.epsilon=vdwE[q];
				return (true);
			}
		}
		
		return(false);
	}
	
	public boolean getNonBondedParameters(Atom atom, AtomNeighbors.Type neighborType, NBParams out) {
		boolean success = getNonBondedParameters(atom.type, out);
		if (!success) {
			return false;
		}
		forcefld.modifyNBParams(atom, neighborType, out);
		return true;
	}
	
	public void getNonBondedParametersOrThrow(Atom atom, AtomNeighbors.Type neighborType, NBParams out) {
		boolean success = getNonBondedParameters(atom, neighborType, out);
		if (!success) {
			throw new Error("couldn't find non-bonded parameters for atom type: " + atom.forceFieldType);
		}
	}
	
	
	public static class VdwParams {
		public double Aij = 0.0;
		public double Bij = 0.0;
	}
	
	public void getVdwParams(Atom atom1, Atom atom2, AtomNeighbors.Type neighborType, VdwParams out) {
		
		// calc vdW params
		// Aij = (ri+rj)^12 * sqrt(ei*ej)
		// Bij = (ri+rj)^6 * sqrt(ei*ej)
		
		// TODO: optimize out this allocation? or will escape analysis during JIT use stack allocation?
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		getNonBondedParametersOrThrow(atom1, neighborType, nbparams1);
		getNonBondedParametersOrThrow(atom2, neighborType, nbparams2);
		double epsilon = Math.sqrt(nbparams1.epsilon*nbparams2.epsilon);
		double radiusSum = nbparams1.r + nbparams2.r;
		double Bij = radiusSum*radiusSum*vdwMultiplier*vdwMultiplier;
		Bij = Bij*Bij*Bij;
		double Aij = Bij*Bij;
		
		Aij *= epsilon;
		Bij *= epsilon;
		
		// vdW scaling by connectivity
		switch (neighborType) {
			case BONDED14:
				Aij *= forcefld.Aij14Factor;
				Bij *= forcefld.Bij14Factor;
			break;
			case NONBONDED:
				Bij *= 2;
			break;
			default:
				throw new IllegalArgumentException("no van der Waals params for closely bonded atoms");
		}
		
		out.Aij = Aij;
		out.Bij = Bij;
	}
        
        
        public double getBondEBL(int atomType1, int atomType2){
            //get the equilibrium bond length for the specified atom types
            for(int bondTypeNum=0; bondTypeNum<bondEBL.length; bondTypeNum++){
                if( (bondAtomType1[bondTypeNum]==atomType1) && (bondAtomType2[bondTypeNum]==atomType2) ){
                    return bondEBL[bondTypeNum];
                }
                if( (bondAtomType2[bondTypeNum]==atomType1) && (bondAtomType1[bondTypeNum]==atomType2) ){
                    return bondEBL[bondTypeNum];//could be listed in either order...
                }
            }
            
            if(printWarnings) {
            System.out.println("Warning: No equilibrium bond length listed for atom types "
                    +atomType1+" and "+atomType2);
            }
            //this is used to get an estimated bond distance matrix, in which case
            //we can estimate using other atom types for the same elements
            
            return Double.NaN;
        }
        
        
        public double estBondEBL(int atomType1, int atomType2){
            //if we don't know the true equilibrium bond length, we can try to estimate it
            //based on another pair of atom types from the same element types
            
            double mass1 = atomAtomicMasses[atomType1];
            double mass2 = atomAtomicMasses[atomType2];
            
            for(int altAtomType1=0; altAtomType1<atomTypeNames.length; altAtomType1++){
                if(atomAtomicMasses[altAtomType1] == mass1){//same element
                    
                    for(int altAtomType2=0; altAtomType2<atomTypeNames.length; altAtomType2++){
                        if(atomAtomicMasses[altAtomType2] == mass2){
                            //Let's try using the alt types as a substitute
                            double bondEBL = getBondEBL(altAtomType1, altAtomType2);
                            if(!Double.isNaN(bondEBL))
                                return bondEBL;
                        }
                    }
                }
            }
            
            throw new RuntimeException("ERROR: Couldn't find any equilibrium bond length"
                    + " for atoms with masses " + mass1 + " and " + mass2);
        }
        
        
        
        //Doubles the size of the a[] String array
	private String [] doubleArraySize(String a[]){
		String tmp[] = new String[a.length*2];
		System.arraycopy(a, 0, tmp, 0, a.length);
		return tmp;
	}
	
	//Doubles the size of the a[] double array
	private double [] doubleArraySize(double a[]){
		double tmp[] = new double[a.length*2];
		System.arraycopy(a, 0, tmp, 0, a.length);
		return tmp;
	}
	
	//Doubles the size of the a[] int array
	private int [] doubleArraySize(int a[]){
		int tmp[] = new int[a.length*2];
		System.arraycopy(a, 0, tmp, 0, a.length);
		return tmp;
	}
	
	//Doubles the size (first dimension only) of the a[] int array
	private int [][] doubleArraySize(int a[][]){
		int tmp[][] = new int[a.length*2][];
		for (int i=0; i<a.length; i++){
			tmp[i] = new int[a[i].length];
			System.arraycopy(a[i], 0, tmp[i], 0, a[i].length);
		}
		return tmp;
	}
	
	//Reduce the a[] String array to keep only the first newSize elements
	private String [] reduceArraySize(String a[], int newSize){
		String tmp[] = new String[newSize];
		System.arraycopy(a, 0, tmp, 0, tmp.length);
		return tmp;
	}
	
	//Reduce the a[] double array to keep only the first newSize elements
	private double [] reduceArraySize(double a[], int newSize){
		double tmp[] = new double[newSize];
		System.arraycopy(a, 0, tmp, 0, tmp.length);
		return tmp;
	}
	
	//Reduce the a[] int array to keep only the first newSize elements
	private int [] reduceArraySize(int a[], int newSize){
		int tmp[] = new int[newSize];
		System.arraycopy(a, 0, tmp, 0, tmp.length);
		return tmp;
	}
	
	private int [][] reduceArraySize(int a[][], int newSize){
		int tmp[][] = new int[newSize][];
		for (int i=0; i<newSize; i++){
			tmp[i] = new int[a[i].length];
			System.arraycopy(a[i], 0, tmp[i], 0, a[i].length);
		}
		return tmp;
	}
}
