/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.energy.forcefield;

import edu.duke.cs.osprey.confspace.SimpleConfSpace.Builder;
import edu.duke.cs.osprey.energy.forcefield.amber.*;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import one.util.streamex.EntryStream;
import one.util.streamex.IntStreamEx;
import one.util.streamex.StreamEx;

import java.io.IOException;
import java.io.Serializable;
import java.util.List;
import java.util.Map;

import static edu.duke.cs.osprey.tools.Log.log;

/**
 * Contains forcefield parameters for a variety of different published forcefields.  By default, {@link ForcefieldParams.Forcefield#AMBER} is used.
 */
public class ForcefieldParams implements Serializable {
    
    private static final long serialVersionUID = 3124964506851762586L;
    
    public static final double coulombConstant = 332.0;
    public static final double solvCutoff = 9.0;
    public static final double solvCutoff2 = solvCutoff*solvCutoff;

    public static boolean printWarnings = true;

    final int atomTypeX = -2; //the atom type number for the X wildcard atom type
    private final int noMatchInt = 9999;
	private final ForcefieldFileParser parameters;

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
    double[] dihedPN = null; // Periodicity
    double[] dihedPhase = null; // Phase
    int[]		impDihedAtomType1 = null;
    int[]		impDihedAtomType2 = null;
    int[]		impDihedAtomType3 = null;
    int[]		impDihedAtomType4 = null;
    // Imprper Dihedral: PK * (1 + cos(PN*phi - PHASE))
    double[] impDihedTerm1 = null; // PK
    double[] impDihedPN = null; // Periodicity
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
    public Forcefield forcefld;
    
    //The solvation parameters object
    public EEF1 eef1parms = null;
    
    // default values copied from resources/config/defalts.cfg
    public double vdwMultiplier = 0.95;
    public double solvScale = 0.5; //the scale factor for the solvation energies
    public double dielectric = 6;	
    public boolean distDepDielect = true;
    public boolean hElect = true;
    public boolean hVDW = true;
    /** for new code, use {@link Builder#setShellDistance} instead **/
    @Deprecated
    public double shellDistCutoff = Double.POSITIVE_INFINITY; //distance cutoff for interactions (angstroms)
    public SolvationForcefield solvationForcefield = SolvationForcefield.EEF1;
	private Map<String, VanDerWaalsRadius> vanDerWaalsMap;
	private Map<String, BondLengthParameter> bondLengthByName;
	private Map<String, AtomSymbolAndMass> atomNamesMap;

	public enum Forcefield {
        AMBER(
            "/config/parm96.dat",
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

    /** creates an Amber 96 forcefield */
    public ForcefieldParams() {
    	this(Forcefield.AMBER);
    }

    public ForcefieldParams(String frcefld) {
        this(Forcefield.valueOf(frcefld.toUpperCase()));
    }

    public ForcefieldParams(Forcefield frcefld) {
		this(frcefld, new ForcefieldFileParser(ForcefieldParams.class.getResourceAsStream(frcefld.paramsPath)));
	}

	public ForcefieldParams(Forcefield ffChoice, ForcefieldFileParser parameterFile) {
    	this.forcefld = ffChoice;
    	this.parameters = parameterFile;
		try {
			parameterFile.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
		assignFromForcefieldParamFile();
        readEEF();
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

	private void readEEF() {
		// Read in the EEF1 solvation parameters
		// TODO: lazy loading of EEF1 params?
		try {
			eef1parms = new EEF1();
			eef1parms.readEEF1parm();
		} catch (Exception ex) {
			throw new Error("can't read solvation params", ex);
		}
	}

	private static String bondLengthKey(AtomSymbolAndMass atom1, AtomSymbolAndMass atom2) {
    	return String.format("%s-%s", atom1.KNDSYM(), atom2.KNDSYM());
	}

	private static String bondLengthKey(Atom first, Atom second) {
    	return String.format("%s-%s", first.forceFieldType, second.forceFieldType);
	}

	//************************************
	// This function reads the AMBER forcefield parameter file
	//  parm96a.dat
	// Although this function is meant to be relatively generic
	//  it is optimized for reading parm96a.dat. Slight changes
	//  will most likely be required to read other parameter
	//  files. Reading of other files should be done in other
	//  functions
	private void assignFromForcefieldParamFile() {

		atomTypeNames = StreamEx.of(parameters.atomSymbolsAndMasses())
				.map(AtomSymbolAndMass::KNDSYM)
                .toArray(String[]::new);

		atomNamesMap = StreamEx.of(parameters.atomSymbolsAndMasses())
                .toMap(AtomSymbolAndMass::KNDSYM, a -> a);

		atomAtomicMasses = StreamEx.of(parameters.atomSymbolsAndMasses())
                .mapToDouble(AtomSymbolAndMass::AMASS)
				.toArray();

		bondAtomType1 = StreamEx.of(parameters.bondLengthParameters())
				.mapToInt(x -> atomTypeToInt(x.IBT().KNDSYM()))
				.toArray();

		bondAtomType2 = StreamEx.of(parameters.bondLengthParameters())
				.mapToInt(x -> atomTypeToInt(x.JBT().KNDSYM()))
				.toArray();

		bondHFC = StreamEx.of(parameters.bondLengthParameters())
				.mapToDouble(BondLengthParameter::RK)
				.toArray();

		bondEBL = StreamEx.of(parameters.bondLengthParameters())
				.mapToDouble(BondLengthParameter::REQ)
				.toArray();

		bondLengthByName = StreamEx.of(parameters.bondLengthParameters())
				.toMap(a -> bondLengthKey(a.IBT(), a.JBT()), a -> a);

		angleAtomType1 = StreamEx.of(parameters.bondAngleParameters())
				.mapToInt(x -> atomTypeToInt(x.ITT().KNDSYM()))
				.toArray();

		angleAtomType2 = StreamEx.of(parameters.bondAngleParameters())
				.mapToInt(x -> atomTypeToInt(x.JTT().KNDSYM()))
				.toArray();

		angleAtomType3 = StreamEx.of(parameters.bondAngleParameters())
				.mapToInt(x -> atomTypeToInt(x.KTT().KNDSYM()))
				.toArray();

		angleHFC = StreamEx.of(parameters.bondAngleParameters())
				.mapToDouble(BondAngleParameter::TK)
				.toArray();

		angleEBA = StreamEx.of(parameters.bondAngleParameters())
				.mapToDouble(BondAngleParameter::TEQ)
				.toArray();

		dihedAtomType1 = StreamEx.of(parameters.dihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.IPT().KNDSYM()))
				.toArray();

		dihedAtomType2 = StreamEx.of(parameters.dihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.JPT().KNDSYM()))
				.toArray();

		dihedAtomType3 = StreamEx.of(parameters.dihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.KPT().KNDSYM()))
				.toArray();

		dihedAtomType4 = StreamEx.of(parameters.dihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.LPT().KNDSYM()))
				.toArray();

		dihedTerm1 = StreamEx.of(parameters.dihederalParameters())
				.mapToDouble(x -> x.PK() / x.IDIVF())
				.toArray();

		dihedPhase = StreamEx.of(parameters.dihederalParameters())
				.mapToDouble(DihederalParameter::PHASE)
				.toArray();

		// If dihedPN is negative then there are one or more additional terms
		//  nothing fancy needs to be done because they will all be read in anyway
		//  but we do need to correct the sign
		dihedPN = StreamEx.of(parameters.dihederalParameters())
				.mapToDouble(x -> Math.abs(x.PN()))
                .toArray();

		numGeneralDihedParams = StreamEx.of(parameters.dihederalParameters())
				.mapToInt(dihedral -> List.of(dihedral.IPT(), dihedral.JPT(), dihedral.KPT(), dihedral.LPT()).contains(ForcefieldFileParser.WildcardAtom) ? 1 : 0)
                .sum();

		impDihedAtomType1 = StreamEx.of(parameters.improperDihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.IPT().KNDSYM()))
				.toArray();

		impDihedAtomType2 = StreamEx.of(parameters.improperDihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.JPT().KNDSYM()))
				.toArray();

		impDihedAtomType3 = StreamEx.of(parameters.improperDihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.KPT().KNDSYM()))
				.toArray();

		impDihedAtomType4 = StreamEx.of(parameters.improperDihederalParameters())
				.mapToInt(x -> atomTypeToInt(x.LPT().KNDSYM()))
				.toArray();

		impDihedTerm1 = StreamEx.of(parameters.improperDihederalParameters())
				.mapToDouble(ImproperDihederalParameter::PK)
				.toArray();

		impDihedPhase = StreamEx.of(parameters.improperDihederalParameters())
				.mapToDouble(ImproperDihederalParameter::PHASE)
				.toArray();

		impDihedPN = StreamEx.of(parameters.improperDihederalParameters())
				.mapToDouble(ImproperDihederalParameter::PN)
				.toArray();

		var map =
				StreamEx.of(parameters.equivalencingAtomsForNonBonded6_12PotentialParameters())
						.groupingBy(EquivalencingAtom::IORG);
		equivAtoms = StreamEx
				.of(EntryStream.of(map)
						.map(entry ->
								IntStreamEx.of(atomTypeToInt(entry.getKey().KNDSYM()))
										.append(entry.getValue().stream().mapToInt(x -> atomTypeToInt(x.IEQV().KNDSYM())))
										.toArray()
						)).toArray(int[][]::new);

		// TODO: handle when these aren't vdW parameters.
		vdwAtomType1 = StreamEx.of(parameters.vanDerWaalsRadii())
				.mapToInt(x -> atomTypeToInt(x.LTYNB().KNDSYM()))
				.toArray();

		vdwR = StreamEx.of(parameters.vanDerWaalsRadii())
				.mapToDouble(VanDerWaalsRadius::R)
				.toArray();

		vdwE = StreamEx.of(parameters.vanDerWaalsRadii())
				.mapToDouble(VanDerWaalsRadius::EDEP)
				.toArray();

		var vdwRadii = StreamEx.of(parameters.vanDerWaalsRadii())
				.toMap(a -> a.LTYNB().KNDSYM(), a -> a);
		vanDerWaalsMap = StreamEx.of(parameters.equivalencingAtomsForNonBonded6_12PotentialParameters())
				.mapToEntry(a -> a.IEQV().KNDSYM(), a -> vdwRadii.get(a.IORG().KNDSYM()))
				.filterKeys(a -> !vdwRadii.containsKey(a)) // there are equivalencing atoms that are also listed explicitly
                .append(vdwRadii)
				.toMap();

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

	// This function returns the r and epsilon parameters for a given atom type
	public boolean getNonBondedParameters(Atom atom, NBParams out) {

        if (vanDerWaalsMap.containsKey(atom.forceFieldType)) {
			var vdwRadii = vanDerWaalsMap.get(atom.forceFieldType);
			out.epsilon = vdwRadii.EDEP();
			out.r = vdwRadii.R();
			return true;
		}

        return false;
	}
	
	public boolean getNonBondedParameters(Atom atom, AtomNeighbors.Type neighborType, NBParams out) {
		boolean success = getNonBondedParameters(atom, out);
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
		public final double Aij;
		public final double Bij;

		VdwParams(double aij, double bij) {
		    Aij = aij;
		    Bij = bij;
		}
	}
	
	public VdwParams getVdwParams(Atom atom1, Atom atom2, AtomNeighbors.Type neighborType) {
		
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

		return new VdwParams(Aij, Bij);
	}

	private double getBondEquilibriumLength(AtomSymbolAndMass atom1, AtomSymbolAndMass atom2) {
		var forwardKey = bondLengthKey(atom1, atom2);
		var reverseKey = bondLengthKey(atom2, atom1);

		if (bondLengthByName.containsKey(forwardKey)) {
			return bondLengthByName.get(forwardKey).REQ();
		} else if (bondLengthByName.containsKey(reverseKey)) {
			return bondLengthByName.get(reverseKey).REQ();
		}

		if (printWarnings) {
			System.out.println("Warning: No equilibrium bond length listed for atom types " + atom1 + " and " + atom2);
		}

		//this is used to get an estimated bond distance matrix, in which case
		//we can estimate using other atom types for the same elements
		return Double.NaN;
	}

	//get the equilibrium bond length for the specified atom types
	public double getBondEquilibriumLength(Atom atom1, Atom atom2) {
		var atomSymbolAndMass1 = atomNamesMap.get(atom1.forceFieldType);
		var atomSymbolAndMass2 = atomNamesMap.get(atom2.forceFieldType);
		return getBondEquilibriumLength(atomSymbolAndMass1, atomSymbolAndMass2);
	}


	public double estBondEBL(Atom atom1, Atom atom2) {
		//if we don't know the true equilibrium bond length, we can try to estimate it
		//based on another pair of atom types from the same element types

		var mass1 = atomNamesMap.get(atom1.forceFieldType).AMASS();
		var mass2 = atomNamesMap.get(atom2.forceFieldType).AMASS();

		for (var firstAtom : parameters.atomSymbolsAndMasses()) {
			if (firstAtom.AMASS() == mass1) {
				for (var secondAtom : parameters.atomSymbolsAndMasses()) {
					if (secondAtom.AMASS() == mass2) {
						double bondEBL = getBondEquilibriumLength(firstAtom, secondAtom);
						if (!Double.isNaN(bondEBL))
							return bondEBL;
					}
				}
			}
		}

		throw new RuntimeException("ERROR: Couldn't find any equilibrium bond length"
				+ " for atoms with masses " + mass1 + " and " + mass2);
	}
}
