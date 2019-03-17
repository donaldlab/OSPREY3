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

package edu.duke.cs.osprey;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.DihedralRotation;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.*;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.analysis.*;
import edu.duke.cs.osprey.tools.*;

import java.io.*;
import java.util.*;
import java.util.function.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


public class FlexLab {

	public static void main(String[] args)
	throws Exception {
		checkRotamerClashes();
		//top8000Dihedrals("arg");
		//top8000Angles("leu");
		//top8000Tetrahedrals("leu");
		//top8000Methyls();
		//top8000Clashes();
		//lovellRotamers("asp");
		//energyLandscape("leu");
		//top8000RotamerStats("asp");
		//checkTemplates();
	}

	private static PDBScanner scanner = new PDBScanner(
		new File("/home/jeff/dlab/top8000"),

		// these all have duplicated residues for some reason
		"3o2rFH_D.pdb",
		"1xpmFH_D.pdb",
		"2hmqFH_D.pdb",
		"3gjuFH_A.pdb",
		"2zs0FH_D.pdb",
		"3lf3FH_A.pdb",
		"3blnFH_A.pdb",
		"2wtaFH_A.pdb",
		"1eisFH_A.pdb",
		"1yslFH_B.pdb",
		"3piuFH_A.pdb",
		"2j6vFH_B.pdb",
		"2vqrFH_A.pdb"
	);

	private static ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();
	private static ForcefieldParams ffparams = new ForcefieldParams();

	// define all the dihedrals
	private static MeasurementLibrary chiLib = new MeasurementLibrary();
	static {

		// arg
		chiLib.add("ARG", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("ARG", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD"));
		chiLib.add("ARG", new MeasurementLibrary.DihedralAngle("chi3", "CB", "CG", "CD", "NE"));
		chiLib.add("ARG", new MeasurementLibrary.DihedralAngle("chi4", "CG", "CD", "NE", "CZ"));
		// guanidine group is pretty planar, not much variation in chi5
		//chiLib.add("ARG", new MeasurementLibrary.DihedralAnglesMinDist("chi5", 0, "CD", "NE", "CZ", "NH1", "NH2"));

		// his
		for (String type : Arrays.asList("HIP", "HID", "HIE")) {
			chiLib.add(type, new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
			chiLib.add(type, new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "ND1"));
		}

		// lys
		chiLib.add("LYS", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("LYS", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD"));
		chiLib.add("LYS", new MeasurementLibrary.DihedralAngle("chi3", "CB", "CG", "CD", "CE"));
		chiLib.add("LYS", new MeasurementLibrary.DihedralAngle("chi4", "CG", "CD", "CE", "NZ"));

		// asp
		chiLib.add("ASP", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("ASP", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "OD1"));

		// glu
		chiLib.add("GLU", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("GLU", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD"));
		chiLib.add("GLU", new MeasurementLibrary.DihedralAngle("chi3", "CB", "CG", "CD", "OE1"));

		// cys
		chiLib.add("CYS", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "SG"));

		// val
		chiLib.add("VAL", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG1"));

		// ile
		chiLib.add("ILE", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG1"));
		chiLib.add("ILE", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG1", "CD1"));

		// leu
		chiLib.add("LEU", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("LEU", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD1"));

		// met
		chiLib.add("MET", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("MET", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "SD"));
		chiLib.add("MET", new MeasurementLibrary.DihedralAngle("chi3", "CB", "CG", "SD", "CE"));

		// phe
		chiLib.add("PHE", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("PHE", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD1"));

		// tyr
		chiLib.add("TYR", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("TYR", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD1"));
		chiLib.add("TYR", new MeasurementLibrary.DihedralAngle("OH", "CE1", "CZ", "OH", "HH"));

		// trp
		chiLib.add("TRP", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("TRP", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD1"));

		// ser
		chiLib.add("SER", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "OG"));
		chiLib.add("SER", new MeasurementLibrary.DihedralAngle("OH", "CA", "CB", "OG", "HG"));

		// thr
		chiLib.add("THR", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "OG1"));
		chiLib.add("THR", new MeasurementLibrary.DihedralAngle("OH", "CA", "CB", "OG1", "HG1"));

		// asn
		chiLib.add("ASN", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("ASN", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "OD1"));

		// gln
		chiLib.add("GLN", new MeasurementLibrary.DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("GLN", new MeasurementLibrary.DihedralAngle("chi2", "CA", "CB", "CG", "CD"));
		chiLib.add("GLN", new MeasurementLibrary.DihedralAngle("chi3", "CB", "CG", "CD", "OE1"));
	}

	// define all the methyl groups (and NH3 groups too)
	private static MeasurementLibrary methylLib = new MeasurementLibrary();
	static {

		// lys
		methylLib.add("LYS", new MeasurementLibrary.DihedralAnglesMinDist("NH3", 180, "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3"));

		// ala
		methylLib.add("ALA", new MeasurementLibrary.DihedralAnglesMinDist("CH3", 180, "N", "CA", "CB", "HB1", "HB2", "HB3"));

		// val
		methylLib.add("VAL", new MeasurementLibrary.DihedralAnglesMinDist("CH3-1", 180, "CA", "CB", "CG1", "HG11", "HG12", "HG13"));
		methylLib.add("VAL", new MeasurementLibrary.DihedralAnglesMinDist("CH3-2", 180, "CA", "CB", "CG2", "HG21", "HG22", "HG23"));

		// ile
		methylLib.add("ILE", new MeasurementLibrary.DihedralAnglesMinDist("CH3-1", 180, "CB", "CG1", "CD1", "HD11", "HD12", "HD13"));
		methylLib.add("ILE", new MeasurementLibrary.DihedralAnglesMinDist("CH3-2", 180, "CA", "CB", "CG2", "HG21", "HG22", "HG23"));

		// leu
		methylLib.add("LEU", new MeasurementLibrary.DihedralAnglesMinDist("CH3-1", 180, "CB", "CG", "CD1", "HD11", "HD12", "HD13"));
		methylLib.add("LEU", new MeasurementLibrary.DihedralAnglesMinDist("CH3-2", 180, "CB", "CG", "CD2", "HD21", "HD22", "HD23"));

		// met
		methylLib.add("MET", new MeasurementLibrary.DihedralAnglesMinDist("CH3", 180, "CG", "SD", "CE", "HE1", "HE2", "HE3"));

		// thr
		methylLib.add("THR", new MeasurementLibrary.DihedralAnglesMinDist("CH3", 180, "CA", "CB", "CG2", "HG21", "HG22", "HG23"));
	}

	// define all the bond angles
	private static MeasurementLibrary angleLib = new MeasurementLibrary();
	static {
		angleLib.add("LEU", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("LEU", new MeasurementLibrary.BondAngle("C", "CA", "CB"));
		angleLib.add("LEU", new MeasurementLibrary.BondAngle("N", "CA", "C"));

		angleLib.add("LEU", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));
		angleLib.add("LEU", new MeasurementLibrary.BondAngle("CB", "CG", "CD1"));
		angleLib.add("LEU", new MeasurementLibrary.BondAngle("CB", "CG", "CD2"));

		// aromatics
		angleLib.add("HIS", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("HIS", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("PHE", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("PHE", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("TYR", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("TYR", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("TRP", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("TRP", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		// charged
		angleLib.add("ARG", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("ARG", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("LYS", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("LYS", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("ASP", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("ASP", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("GLU", new MeasurementLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("GLU", new MeasurementLibrary.BondAngle("CA", "CB", "CG"));
	}

	// define all the tetrahedral angles
	private static MeasurementLibrary tetraLib = new MeasurementLibrary();
	static {
		tetraLib.add("LEU", new MeasurementLibrary.TetrahedralInPlaneAngle("N", "CA", "C", "CB"));
		tetraLib.add("LEU", new MeasurementLibrary.TetrahedralOutOfPlaneAngle("N", "CA", "C", "CB"));
	}

	static class Rotamer {

		public final String name;
		public final double percent;
		public final SmallAngleVoxel voxel;

		public Rotamer(String name, double percent, SmallAngleVoxel voxel) {
			this.name = name;
			this.percent = percent;
			this.voxel = voxel;
		}

		public static Rotamer ofModes(String name, double percent, double ... modalAngles) {
			return new Rotamer(name, percent, new SmallAngleVoxel(
				Arrays.stream(modalAngles)
					.mapToObj(angle -> new SmallAngleVoxel.Interval(-9, angle, 9))
					.toArray(size -> new SmallAngleVoxel.Interval[size])
			));
		}

		public boolean matches(double[] dihedrals) {
			if (voxel == null) {
				return true;
			}
			return voxel.contains(dihedrals);
		}
	}

	static class RotamerLibrary {

		public final Map<String,List<Rotamer>> rotamers = new HashMap<>();

		public List<Rotamer> get(String type) {
			return rotamers.computeIfAbsent(type.toUpperCase(), t -> new ArrayList<>());
		}

		public void add(String type, Rotamer rot) {
			get(type).add(rot);
		}

		public Rotamer find(String type, double[] dihedrals) {
			for (Rotamer rot : get(type)) {
				if (rot.matches(dihedrals)) {
					return rot;
				}
			}
			return null;
		}
	}

	// define a rotamer library based on voxels
	private static RotamerLibrary rotamers = new RotamerLibrary();
	static {

		// re-define the rotamer library based on the Lovell data, but with names and populations
		rotamers.add("VAL", Rotamer.ofModes("p",  6.20,   64));
		rotamers.add("VAL", Rotamer.ofModes("t", 75.56,  175));
		rotamers.add("VAL", Rotamer.ofModes("m", 17.94,  -60));

		// leu
		rotamers.add("LEU", Rotamer.ofModes("pp",  0.45,   62,   80));
		rotamers.add("LEU", Rotamer.ofModes("tp", 30.12, -177,   65));
		rotamers.add("LEU", Rotamer.ofModes("tt",  1.37, -172,  145));
		rotamers.add("LEU", Rotamer.ofModes("mp",  2.36,  -85,   65));
		rotamers.add("LEU", Rotamer.ofModes("mt", 64.54,  -65,  175));

		// ile
		rotamers.add("ILE", Rotamer.ofModes("pp",  0.35,   62,  100));
		rotamers.add("ILE", Rotamer.ofModes("pt", 12.33,   62,  170));
		rotamers.add("ILE", Rotamer.ofModes("tp",  2.61, -177,  66));
		rotamers.add("ILE", Rotamer.ofModes("tt",  5.81, -177,  165));
		rotamers.add("ILE", Rotamer.ofModes("mp",  0.87,  -65,  100));
		rotamers.add("ILE", Rotamer.ofModes("mt", 62.03,  -65,  170));
		rotamers.add("ILE", Rotamer.ofModes("mm", 15.70,  -57,  -60));

		// phe
		rotamers.add("PHE", Rotamer.ofModes("p90",  11.17,   62,   90));
		rotamers.add("PHE", Rotamer.ofModes("t80",  34.27, -177,   80));
		rotamers.add("PHE", Rotamer.ofModes("m-80", 47.38,  -65,  -85));
		rotamers.add("PHE", Rotamer.ofModes("m-10",  6.89,  -65,  -30));

		// tyr
		rotamers.add("TYR", Rotamer.ofModes("p90:0",    11.57,   62,   90,   0));
		rotamers.add("TYR", Rotamer.ofModes("p90:180",  11.57,   62,   90, 180));
		rotamers.add("TYR", Rotamer.ofModes("t80:0",    34.53, -177,   80,   0));
		rotamers.add("TYR", Rotamer.ofModes("t80:180",  34.53, -177,   80, 180));
		rotamers.add("TYR", Rotamer.ofModes("m-80:0",   48.01,  -65,  -85,   0));
		rotamers.add("TYR", Rotamer.ofModes("m-80:180", 48.01,  -65,  -85, 180));
		rotamers.add("TYR", Rotamer.ofModes("m-10:0",    5.55,  -65,  -30,   0));
		rotamers.add("TYR", Rotamer.ofModes("m-10:180",  5.55,  -65,  -30, 180));

		// trp
		rotamers.add("TRP", Rotamer.ofModes("p-90",  10.35,   62,  -90));
		rotamers.add("TRP", Rotamer.ofModes("p90",    5.19,   62,   90));
		rotamers.add("TRP", Rotamer.ofModes("t-100", 15.46, -177, -105));
		rotamers.add("TRP", Rotamer.ofModes("t60",   18.09, -177,   90));
		rotamers.add("TRP", Rotamer.ofModes("m-90",   5.13,  -65,  -90));
		rotamers.add("TRP", Rotamer.ofModes("m-10",  11.73,  -65,   -5));
		rotamers.add("TRP", Rotamer.ofModes("m100",  33.76,  -65,   95));

		// cys
		rotamers.add("CYS", Rotamer.ofModes("p", 17.73,   62));
		rotamers.add("CYS", Rotamer.ofModes("t", 26.33, -177));
		rotamers.add("CYS", Rotamer.ofModes("m", 55.67,  -65));

		// met
		rotamers.add("MET", Rotamer.ofModes("ptp",  2.41,   62,  180,   75));
		rotamers.add("MET", Rotamer.ofModes("ptm",  2.23,   62,  180,  -75));
		rotamers.add("MET", Rotamer.ofModes("tpp",  6.78, -177,   65,   75));
		rotamers.add("MET", Rotamer.ofModes("tpt",  2.31, -177,   65,  180));
		rotamers.add("MET", Rotamer.ofModes("ttp",  7.42, -177,  180,   75));
		rotamers.add("MET", Rotamer.ofModes("ttt",  3.39, -177,  180,  180));
		rotamers.add("MET", Rotamer.ofModes("ttm",  6.69, -177,  180,  -75));
		rotamers.add("MET", Rotamer.ofModes("mtp", 16.76,  -67,  180,   75));
		rotamers.add("MET", Rotamer.ofModes("mtt",  9.18,  -67,  180,  180));
		rotamers.add("MET", Rotamer.ofModes("mtm", 11.02,  -67,  180,  -75));
		rotamers.add("MET", Rotamer.ofModes("mmp",  3.10,  -65,  -65,  103));
		rotamers.add("MET", Rotamer.ofModes("mmt",  3.55,  -65,  -65,  180));
		rotamers.add("MET", Rotamer.ofModes("mmm", 19.97,  -65,  -65,  -70));

		// ser
		rotamers.add("SER", Rotamer.ofModes("p:-60",  48.44,   62,  -60));
		rotamers.add("SER", Rotamer.ofModes("p:60",   48.44,   62,   60));
		rotamers.add("SER", Rotamer.ofModes("p:180",  48.44,   62,  180));
		rotamers.add("SER", Rotamer.ofModes("p:0",    48.44,   62,    0));
		rotamers.add("SER", Rotamer.ofModes("p:120",  48.44,   62,  120));
		rotamers.add("SER", Rotamer.ofModes("p:-120", 48.44,   62, -120));
		rotamers.add("SER", Rotamer.ofModes("t:-60",  22.97, -177,  -60));
		rotamers.add("SER", Rotamer.ofModes("t:60",   22.97, -177,   60));
		rotamers.add("SER", Rotamer.ofModes("t:180",  22.97, -177,  180));
		rotamers.add("SER", Rotamer.ofModes("t:0",    22.97, -177,    0));
		rotamers.add("SER", Rotamer.ofModes("t:120",  22.97, -177,  120));
		rotamers.add("SER", Rotamer.ofModes("t:-120", 22.97, -177, -120));
		rotamers.add("SER", Rotamer.ofModes("m:-60",  28.30,  -65,  -60));
		rotamers.add("SER", Rotamer.ofModes("m:60",   28.30,  -65,   60));
		rotamers.add("SER", Rotamer.ofModes("m:180",  28.30,  -65,  180));
		rotamers.add("SER", Rotamer.ofModes("m:0",    28.30,  -65,    0));
		rotamers.add("SER", Rotamer.ofModes("m:120",  28.30,  -65,  120));
		rotamers.add("SER", Rotamer.ofModes("m:-120", 28.30,  -65, -120));

		// thr
		rotamers.add("THR", Rotamer.ofModes("p:-60",  48.14,   62,  -60));
		rotamers.add("THR", Rotamer.ofModes("p:60",   48.14,   62,   60));
		rotamers.add("THR", Rotamer.ofModes("p:180",  48.14,   62,  180));
		rotamers.add("THR", Rotamer.ofModes("p:0",    48.14,   62,    0));
		rotamers.add("THR", Rotamer.ofModes("p:120",  48.14,   62,  120));
		rotamers.add("THR", Rotamer.ofModes("p:-120", 48.14,   62, -120));
		rotamers.add("THR", Rotamer.ofModes("p:-60",   6.91, -175,  -60));
		rotamers.add("THR", Rotamer.ofModes("p:60",    6.91, -175,   60));
		rotamers.add("THR", Rotamer.ofModes("p:180",   6.91, -175,  180));
		rotamers.add("THR", Rotamer.ofModes("p:0",     6.91, -175,    0));
		rotamers.add("THR", Rotamer.ofModes("p:120",   6.91, -175,  120));
		rotamers.add("THR", Rotamer.ofModes("p:-120",  6.91, -175, -120));
		rotamers.add("THR", Rotamer.ofModes("p:-60",  44.64,  -65,  -60));
		rotamers.add("THR", Rotamer.ofModes("p:60",   44.64,  -65,   60));
		rotamers.add("THR", Rotamer.ofModes("p:180",  44.64,  -65,  180));
		rotamers.add("THR", Rotamer.ofModes("p:0",    44.64,  -65,    0));
		rotamers.add("THR", Rotamer.ofModes("p:120",  44.64,  -65,  120));
		rotamers.add("THR", Rotamer.ofModes("p:-120", 44.64,  -65, -120));

		// lys
		rotamers.add("LYS", Rotamer.ofModes("ptpt",  0.42,  62,  180,   68,  180));
		rotamers.add("LYS", Rotamer.ofModes("pttp",  0.69,  62,  180,  180,   65));
		rotamers.add("LYS", Rotamer.ofModes("pttt",  3.98,  62,  180,  180,  180));
		rotamers.add("LYS", Rotamer.ofModes("pttt",  0.77,  62,  180,  180,  -65));
		rotamers.add("LYS", Rotamer.ofModes("ptmt",  0.54,  62,  180,  -68,  180));
		rotamers.add("LYS", Rotamer.ofModes("tptp",  1.17,-177,   68,  180,   65));
		rotamers.add("LYS", Rotamer.ofModes("tptt",  3.53,-177,   68,  180,  180));
		rotamers.add("LYS", Rotamer.ofModes("tptm",  0.57,-177,   68,  180,  -65));
		rotamers.add("LYS", Rotamer.ofModes("ttpp",  0.66,-177,  180,   68,   65));
		rotamers.add("LYS", Rotamer.ofModes("ttpt",  2.54,-177,  180,   68,  180));
		rotamers.add("LYS", Rotamer.ofModes("tttp",  3.54,-177,  180,  180,   65));
		rotamers.add("LYS", Rotamer.ofModes("tttt", 14.48,-177,  180,  180,  180));
		rotamers.add("LYS", Rotamer.ofModes("tttm",  3.38,-177,  180,  180,  -65));
		rotamers.add("LYS", Rotamer.ofModes("ttmt",  1.94,-177,  180,  -68,  180));
		rotamers.add("LYS", Rotamer.ofModes("ttmm",  0.57,-177,  180,  -68,  -65));
		rotamers.add("LYS", Rotamer.ofModes("mptt",  0.36, -90,   68,  180,  180));
		rotamers.add("LYS", Rotamer.ofModes("mtpp",  1.13, -67,  180,   68,   65));
		rotamers.add("LYS", Rotamer.ofModes("mtpt",  3.90, -67,  180,   68,  180));
		rotamers.add("LYS", Rotamer.ofModes("mttp",  4.06, -67,  180,  180,   65));
		rotamers.add("LYS", Rotamer.ofModes("mttt", 24.68, -67,  180,  180,  180));
		rotamers.add("LYS", Rotamer.ofModes("mttm",  5.25, -67,  180,  180,  -65));
		rotamers.add("LYS", Rotamer.ofModes("mtmt",  3.77, -67,  180,  -68,  180));
		rotamers.add("LYS", Rotamer.ofModes("mtmm",  1.22, -67,  180,  -68,  -65));
		rotamers.add("LYS", Rotamer.ofModes("mmtp",  1.33, -62,  -68,  180,   65));
		rotamers.add("LYS", Rotamer.ofModes("mmtt",  9.01, -62,  -68,  180,  180));
		rotamers.add("LYS", Rotamer.ofModes("mmtm",  2.09, -62,  -68,  180,  -65));
		rotamers.add("LYS", Rotamer.ofModes("mmmt",  1.56, -62,  -68,  -68,  180));

		// arg
		rotamers.add("ARG", Rotamer.ofModes("ptp90",    0.48,   62,  180,   65,   85));
		rotamers.add("ARG", Rotamer.ofModes("ptp-170",  0.84,  62,  180,   65, -175));
		rotamers.add("ARG", Rotamer.ofModes("ptt90",    1.76,  62,  180,  180,   85));
		rotamers.add("ARG", Rotamer.ofModes("ptt180",   1.77,  62,  180,  180,  180));
		rotamers.add("ARG", Rotamer.ofModes("ptt-90",   1.57,  62,  180,  180,  -85));
		rotamers.add("ARG", Rotamer.ofModes("ptm160",   1.08,  62,  180,  -65,  175));
		rotamers.add("ARG", Rotamer.ofModes("ptm-80",   0.46,  62,  180,  -65,  -85));
		rotamers.add("ARG", Rotamer.ofModes("tpp80",    0.78,-177,   65,   65,   85));
		rotamers.add("ARG", Rotamer.ofModes("tpp-160",  1.07,-177,   65,   65, -175));
		rotamers.add("ARG", Rotamer.ofModes("tpt90",    1.41,-177,   65,  180,   85));
		rotamers.add("ARG", Rotamer.ofModes("tpt170",   1.78,-177,   65,  180,  180));
		rotamers.add("ARG", Rotamer.ofModes("ttp80",    4.09,-177,  180,   65,   85));
		rotamers.add("ARG", Rotamer.ofModes("ttp-170",  3.31,-177,  180,   65, -175));
		rotamers.add("ARG", Rotamer.ofModes("ttp-110",  1.34,-177,  180,   65, -105));
		rotamers.add("ARG", Rotamer.ofModes("ttt90",    2.28,-177,  180,  180,   85));
		rotamers.add("ARG", Rotamer.ofModes("ttt180",   5.04,-177,  180,  180,  180));
		rotamers.add("ARG", Rotamer.ofModes("ttt-90",   2.98,-177,  180,  180,  -85));
		rotamers.add("ARG", Rotamer.ofModes("ttm110",   1.56,-177,  180,  -65,  105));
		rotamers.add("ARG", Rotamer.ofModes("ttm170",   2.84,-177,  180,  -65,  175));
		rotamers.add("ARG", Rotamer.ofModes("ttm-80",   3.24,-177,  180,  -65,  -85));
		rotamers.add("ARG", Rotamer.ofModes("mtp85",    4.00, -67,  180,   65,   85));
		rotamers.add("ARG", Rotamer.ofModes("mtp180",   5.40, -67,  180,   65, -175));
		rotamers.add("ARG", Rotamer.ofModes("mtp-110",  1.01, -67,  180,   65, -105));
		rotamers.add("ARG", Rotamer.ofModes("mtt90",    5.30, -67,  180,  180,   85));
		rotamers.add("ARG", Rotamer.ofModes("mtt180",   9.90, -67,  180,  180,  180));
		rotamers.add("ARG", Rotamer.ofModes("mtt-85",   6.13, -67,  180,  180,  -85));
		rotamers.add("ARG", Rotamer.ofModes("mtm110",   1.68, -67,  180,  -65,  105));
		rotamers.add("ARG", Rotamer.ofModes("mtm180",   5.19, -67,  180,  -65,  175));
		rotamers.add("ARG", Rotamer.ofModes("mtm-85",   6.14, -67, -167,  -65,  -85));
		rotamers.add("ARG", Rotamer.ofModes("mmt90",    1.22, -62,  -68,  180,   85));
		rotamers.add("ARG", Rotamer.ofModes("mmt180",   2.59, -62,  -68,  180,  180));
		rotamers.add("ARG", Rotamer.ofModes("mmt-90",   3.08, -62,  -68,  180,  -85));
		rotamers.add("ARG", Rotamer.ofModes("mmm160",   2.05, -62,  -68,  -65,  175));
		rotamers.add("ARG", Rotamer.ofModes("mmm-85",   2.20, -62,  -68,  -65,  -85));

		// his
		for (String type : Arrays.asList("HIP", "HID", "HIE")) {
			rotamers.add(type, Rotamer.ofModes("p-80",   7.39,   62,  -75));
			rotamers.add(type, Rotamer.ofModes("p90",    5.01,   62,   80));
			rotamers.add(type, Rotamer.ofModes("t-170",  4.47, -177, -165));
			rotamers.add(type, Rotamer.ofModes("t-90",  11.93, -177,  -80));
			rotamers.add(type, Rotamer.ofModes("t70",   17.01, -177,   60));
			rotamers.add(type, Rotamer.ofModes("m-70",  31.73,  -65,  -70));
			rotamers.add(type, Rotamer.ofModes("m170",   9.05,  -65,  165));
			rotamers.add(type, Rotamer.ofModes("m90",   13.14,  -65,   80));
		}

		// asp
		rotamers.add("ASP", Rotamer.ofModes("p0:1", 16.24/2,   62,  -10));
		rotamers.add("ASP", Rotamer.ofModes("p0:2", 16.24/2,   62,   30));
		rotamers.add("ASP", Rotamer.ofModes("t0",     23.65, -177,    0));
		rotamers.add("ASP", Rotamer.ofModes("t70",     8.33, -177,   65));
		rotamers.add("ASP", Rotamer.ofModes("m-30",   51.48,  -70,  -15));

		// glu
		rotamers.add("GLU", Rotamer.ofModes("pt0",    4.87,   62,  180,  -20));
		rotamers.add("GLU", Rotamer.ofModes("pm20",   2.58,   70,  -80,    0));
		rotamers.add("GLU", Rotamer.ofModes("tp30",   8.03, -177,   65,   10));
		rotamers.add("GLU", Rotamer.ofModes("tt0",   23.69, -177,  180,    0));
		rotamers.add("GLU", Rotamer.ofModes("tm-30",  1.50, -177,  -80,  -25));
		rotamers.add("GLU", Rotamer.ofModes("mp0",    6.39,  -65,   85,    0));
		rotamers.add("GLU", Rotamer.ofModes("mt-10", 36.58,  -67,  180,  -10));
		rotamers.add("GLU", Rotamer.ofModes("mm-30", 15.80,  -65,  -65,  -40));

		// asn
		rotamers.add("ASN", Rotamer.ofModes("p0:1",   14.00/2,   62,  -10));
		rotamers.add("ASN", Rotamer.ofModes("p0:2",   14.00/2,   62,   30));
		rotamers.add("ASN", Rotamer.ofModes("t0:1",   29.10/2, -174,  -20));
		rotamers.add("ASN", Rotamer.ofModes("t0:2",   29.10/2, -177,   30));
		rotamers.add("ASN", Rotamer.ofModes("m-40:1", 49.01/2,  -65,  -20));
		rotamers.add("ASN", Rotamer.ofModes("m-40:2", 49.01/2,  -65,  -75));
		rotamers.add("ASN", Rotamer.ofModes("m110",      7.46,  -65,  120));

		// gln
		rotamers.add("GLN", Rotamer.ofModes("pt0",     5.08,   62,  180,   20));
		rotamers.add("GLN", Rotamer.ofModes("pm20",    1.31,   70,  -75,    0));
		rotamers.add("GLN", Rotamer.ofModes("tp-100",  1.44, -177,   65, -100));
		rotamers.add("GLN", Rotamer.ofModes("tp40",    9.75, -177,   65,   60));
		rotamers.add("GLN", Rotamer.ofModes("tt0",    18.69, -177,  180,    0));
		rotamers.add("GLN", Rotamer.ofModes("mp10",    3.25,  -65,   85,    0));
		rotamers.add("GLN", Rotamer.ofModes("mt0",    38.71,  -67,  180,  -25));
		rotamers.add("GLN", Rotamer.ofModes("mm-40",  16.05,  -65,  -65,  -40));
		rotamers.add("GLN", Rotamer.ofModes("mm110",   3.09,  -65,  -65,  100));

		/* TODO: build a better rotamer library
		// LEU: +-9 voxels <= 50% (ish) coverage (Lovell voxels get 51.5%)
		rotamers.add("LEU", new Rotamer("pp", 0.6, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -15.0,  62.0,  14.0),
			new SmallAngleVoxel.Interval( -19.0,  81.0,  13.0)
		)));
		rotamers.add("LEU", new Rotamer("pt", 0.2, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -10.0,  75.0,   7.0),
			new SmallAngleVoxel.Interval( -10.0, 169.0,  10.0)
		)));
		rotamers.add("LEU", new Rotamer("mm", 0.5, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -17.0, -82.0,  18.0),
			new SmallAngleVoxel.Interval( -13.0, -64.0,  20.0)
		)));
		rotamers.add("LEU", new Rotamer("tt", 2.0, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -19.0,-173.0,  20.0),
			new SmallAngleVoxel.Interval( -21.0, 155.0,  21.0)
		)));
		rotamers.add("LEU", new Rotamer("tp", 28.9, new SmallAngleVoxel( // big!
			new SmallAngleVoxel.Interval( -32.0,-172.0,  43.0),
			new SmallAngleVoxel.Interval( -26.0,  62.0,  30.0)
		)));
		rotamers.add("LEU", new Rotamer("mt", 61.7, new SmallAngleVoxel( // biggest!
			new SmallAngleVoxel.Interval( -42.0, -70.0,  33.0),
			new SmallAngleVoxel.Interval( -36.0, 175.0,  32.0)
		)));
		rotamers.add("LEU", new Rotamer("tm", 0.2, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -10.0,-174.0,   9.0),
			new SmallAngleVoxel.Interval(  -9.0, -76.0,   9.0)
		)));
		rotamers.add("LEU", new Rotamer("mp", 3.5, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -28.0, -82.0,  32.0),
			new SmallAngleVoxel.Interval( -52.0,  60.0,  43.0)
		)));
		rotamers.add("LEU", new Rotamer("xx", 0.2, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval(  -8.0,-146.0,  13.0),
			new SmallAngleVoxel.Interval(  -8.0,-152.0,   8.0)
		)));
		rotamers.add("LEU", new Rotamer("NA", 0.0, null));

		// TRP: +-9 voxels <= 36.2% coverage (Lovell voxels get 29.3%)
		rotamers.add("TRP", new Rotamer("t60", 15.9, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -19.0, 179.0,  19.0),
			new SmallAngleVoxel.Interval( -70.0,  76.0,  26.0)
		)));
		rotamers.add("TRP", new Rotamer("m-10", 11.0, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -16.0, -68.0,  18.0),
			new SmallAngleVoxel.Interval( -32.0,  -8.0,  56.0)
		)));
		rotamers.add("TRP", new Rotamer("m100", 31.2, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -26.0, -68.0,  26.0),
			new SmallAngleVoxel.Interval( -46.0,  96.0,  33.0)
		)));
		rotamers.add("TRP", new Rotamer("t-100", 13.5, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -20.0,-179.0,  24.0),
			new SmallAngleVoxel.Interval( -25.0,-104.0,  33.0)
		)));
		rotamers.add("TRP", new Rotamer("p-90", 8.9, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -20.0,  61.0,  18.0),
			new SmallAngleVoxel.Interval( -17.0, -91.0,  19.0)
		)));
		rotamers.add("TRP", new Rotamer("m-90", 3.6, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -16.0, -67.0,  15.0),
			new SmallAngleVoxel.Interval( -14.0, -92.0,  19.0)
		)));
		rotamers.add("TRP", new Rotamer("p90", 3.7, new SmallAngleVoxel(
			new SmallAngleVoxel.Interval( -16.0,  59.0,  14.0),
			new SmallAngleVoxel.Interval( -14.0,  90.0,  14.0)
		)));
		rotamers.add("TRP", new Rotamer("NA", 0.0, null));

		rotamers.add("PHE", new Rotamer("NA", 0.0, null));
		rotamers.add("ASP", new Rotamer("NA", 0.0, null));
		rotamers.add("TYR", new Rotamer("NA", 0.0, null));
		*/
	}

	static enum TemplateType {

		Current,
		Idealized,
		IdealizedV2;

		public static class Map<T> extends java.util.EnumMap<TemplateType,T> {

			public Map() {
				super(TemplateType.class);
			}

			public <R> Map<R> map(Function<T,R> f) {
				return mapOf(
					f.apply(get(Current)),
					f.apply(get(Idealized)),
					f.apply(get(IdealizedV2))
				);
			}
		}

		public static <T> TemplateType.Map<T> mapOf(T current, T idealized, T idealizedV2) {
			TemplateType.Map<T> map = new TemplateType.Map<>();
			map.put(Current, current);
			map.put(Idealized, idealized);
			map.put(IdealizedV2, idealizedV2);
			return map;
		}
	}

	public static void checkRotamerClashes() {

		List<String> resTypes = Arrays.asList(
			"ARG", "HIP", "HID", "HIE", "LYS", "ASP", "GLU",
			"CYS", "GLY", "PRO",
			"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP",
			"SER", "THR", "ASN", "GLN"
		);

		// make the various template libs
		TemplateType.Map<Map<String,ResidueTemplate>> templateLibs = TemplateType.mapOf(
			pickTemplates(
				resTypes,
				new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
					.build()
			),
			pickTemplates(
				resTypes,
				new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
					.clearTemplateCoords()
					.addTemplateCoords(FileTools.readFile("template coords.txt"))
					.build()
			),
			pickTemplates(
				resTypes,
				new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
					.clearTemplateCoords()
					.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
					.build()
			)
		);

		class TypeStats {

			final ClusterR1 energies = new ClusterR1();
			final ClusterR1 clashes = new ClusterR1();

			MathTools.DoubleBounds clashBounds() {
				if (clashes.isEmpty()) {
					return new MathTools.DoubleBounds(0, 0);
				}
				return clashes.bounds();
			}
		}
		Map<String,TemplateType.Map<TypeStats>> statsByType = new HashMap<>();
		Map<String,TemplateType.Map<TypeStats>> statsOptByType = new HashMap<>();

		for (String resType : resTypes) {
			for (TemplateType templateType : TemplateType.values()) {

				Map<String,ResidueTemplate> templates = templateLibs.get(templateType);
				if (templates == null) {
					continue;
				}

				ResidueTemplate templ = templates.get(resType);

				Probe probe = new Probe();
				probe.matchTemplate(templ);
				AtomConnectivity connectivity = new AtomConnectivity.Builder()
					.addTemplates(Arrays.asList(templ))
					.set15HasNonBonded(false)
					.build();

				// make a molecule from just the template
				Molecule mol = new Molecule();
				Residue res = new Residue(templ.templateRes);
				res.fullName = resType + " A   1";
				res.template = templ;
				res.markIntraResBondsByTemplate();
				mol.residues.add(res);
				res.indexInMolecule = 0;

				Runnable analyzeClashes = () -> {
					probe.getInteractions(res, connectivity).stream()
						.filter(interaction -> interaction.getViolation(0.0) > 0)
						.sorted(Comparator.comparing(interaction -> interaction.getViolation(0.0)))
						.forEach(interaction -> {
							double vdw = calcVDW(interaction.atomPair.a, interaction.atomPair.b, ffparams);
							// TEMP
							//log("\t\t%s   %s     vdw=%.4f", interaction.atomPair, interaction, vdw);
							log("        %4s <-> %-4s   clash=%8.3f A   %12s   %10s   vdw=%.4f kcal/mol",
								interaction.atomPair.a.name,
								interaction.atomPair.b.name,
								interaction.getViolation(0.0),
								interaction.contact,
								interaction.atomPair.attraction == Probe.Attraction.None ? "" : interaction.atomPair.attraction,
								vdw
							);
						});
				};

				Supplier<Double> maxClash = () ->
					probe.getInteractions(res, connectivity).stream()
						.filter(interaction -> interaction.getViolation(0.0) > 0)
						.mapToDouble(interaction -> interaction.getViolation(0.0))
						.max()
						.orElse(0);

				try (EnergyCalculator ecalc = new EnergyCalculator.Builder(mol.residues, ffparams).build()) {

					// make residue interactions for the ecalc
					ResidueInteractions inters = new ResidueInteractions();
					inters.addSingle(res.getPDBResNumber());

					// make the dihedral dofs
					VoxelShape.Rect voxel = new VoxelShape.Rect();
					List<DegreeOfFreedom> dofs = voxel.makeDihedralDOFs(res);

					// get all the dihedral dof bounds
					List<ObjectiveFunction.DofBounds> boundsByRot = new ArrayList<>();
					if (templ.getNumRotamers() > 0) {
						for (int i=0; i<templ.getNumRotamers(); i++) {
							boundsByRot.add(voxel.makeDihedralBounds(templ, i));
						}
					} else {
						boundsByRot.add(new ObjectiveFunction.DofBounds(0));
					}

					// for each rotamer
					for (ObjectiveFunction.DofBounds dofBounds : boundsByRot) {

						Runnable setRotamericPose = () -> {
							for (int d=0; d<dofs.size(); d++) {
								dofs.get(d).apply(dofBounds.getCenter(d));
							}
						};

						setRotamericPose.run();

						// get the rotamer name
						String rotName = null;
						Double rotPercent = null;
						Rotamer rot = rotamers.find(resType, chiLib.measure(res, resType));
						if (rot != null) {
							rotName = rot.name;
							rotPercent = rot.percent;
						}

						TypeStats stats = statsByType
							.computeIfAbsent(resType, t -> new TemplateType.Map<>())
							.computeIfAbsent(templateType, t -> new TypeStats());

						// calc energy and clashes with backbone minimization

						/* TEMP
						log("\t%10s %12s: amber=%7.3f avgClash=%5.3f",
							templateType,
							"rigid",
							rigidEcalc.calcEnergy(pmol, inters).energy,
							avgClash.apply(pmol.mol)
						);
						*/

						setRotamericPose.run();
						double amberEnergy = ecalc.calcEnergy(new ParametricMolecule(mol), inters).energy;
						stats.energies.add(amberEnergy);
						stats.clashes.add(maxClash.get());
						/* TEMP
						log("%3s\t%s\t%.1f\t%s\t%.3f\t%.3f",
							type,
							rot != null ? rot.name : "(none)",
							rot != null ? rot.percent : 100,
							templateType.name(),
							epmol.energy,
							maxClash.apply(epmol.pmol.mol)
						);
						*/
						log("%3s %6s %5.1f%% population   coords=%-10s   AMBER=%7.3f   worstClash=%5.3f",
							resType,
							rot != null ? rot.name : "(none)",
							rot != null ? rot.percent : 100,
							templateType.name(),
							amberEnergy,
							maxClash.get()
						);
						analyzeClashes.run();

						TypeStats statsOpt = statsOptByType
							.computeIfAbsent(resType, t -> new TemplateType.Map<>())
							.computeIfAbsent(templateType, t -> new TypeStats());

						// minimize methyls if possible
						List<MeasurementLibrary.Measurement> methylMeasurements = methylLib.get(resType);
						if (methylMeasurements != null) {
							setRotamericPose.run();
							ParametricMolecule pmol = new ParametricMolecule(mol);
							pmol = addDofs(
								pmol,
								methylMeasurements.stream()
									.map(measurement -> {
										MeasurementLibrary.DihedralAnglesMinDist m = (MeasurementLibrary.DihedralAnglesMinDist)measurement;
										return new BoundedDof(
											new MethylRotation(res, m.a, m.b, m.c, m.d[0], m.d[1], m.d[2]),
											180 - 40,
											180 + 40
										);
									})
									.collect(Collectors.toList())
							);
							EnergyCalculator.EnergiedParametricMolecule epmol = ecalc.calcEnergy(pmol, inters);
							/* TEMP
							double[] methyls = methylLib.measure(res, type);
							log("\t%10s %12s: amber=%7.3f avgClash=%5.3f methyls=%s",
								templateType,
								"min methyls",
								epmol.energy,
								avgClash.apply(pmol.mol),
								degreesToString(methyls)
							);
							*/
							statsOpt.energies.add(epmol.energy);
							statsOpt.clashes.add(maxClash.get());
						}
					}
				}
			}
		}

		log("\nres type stats for");
		for (String type : resTypes) {

			TypeStats statsCurrent = statsByType.get(type).get(TemplateType.Current);
			TypeStats statsIdealized = statsByType.get(type).get(TemplateType.Idealized);
			TypeStats statsCurrentOpt = statsOptByType.get(type).get(TemplateType.Current);
			TypeStats statsIdealizedOpt = statsOptByType.get(type).get(TemplateType.Idealized);

			MathTools.DoubleBounds energiesCurrent = statsCurrent.energies.bounds();
			MathTools.DoubleBounds energiesIdealized = statsIdealized.energies.bounds();
			MathTools.DoubleBounds energiesCurrentOpt = statsCurrentOpt.energies.bounds();
			MathTools.DoubleBounds energiesIdealizedOpt = statsIdealizedOpt.energies.bounds();

			MathTools.DoubleBounds clashesCurrent = statsCurrent.clashBounds();
			MathTools.DoubleBounds clashesIdealized = statsIdealized.clashBounds();
			MathTools.DoubleBounds clashesCurrentOpt = statsCurrentOpt.clashBounds();
			MathTools.DoubleBounds clashesIdealizedOpt = statsIdealizedOpt.clashBounds();

			log("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
				type,
				energiesCurrent.lower, energiesCurrent.upper,
				clashesCurrent.lower, clashesCurrent.upper,
				energiesIdealized.lower, energiesIdealized.upper,
				clashesIdealized.lower, clashesIdealized.upper,
				energiesCurrentOpt.lower, energiesCurrentOpt.upper,
				clashesCurrentOpt.lower, clashesCurrentOpt.upper,
				energiesIdealizedOpt.lower, energiesIdealizedOpt.upper,
				clashesIdealizedOpt.lower, clashesIdealizedOpt.upper
			);
		}
	}

	private static ParametricMolecule makeMolecule(SimpleConfSpace confSpace, SimpleConfSpace.ResidueConf rc) {

		// easy mode
		//return confSpace.makeMolecule(new RCTuple(0, rc.index));

		// get the name of the original res
		SimpleConfSpace.Position pos = confSpace.positions.get(0);

		// make the molecule from the template res
		Molecule mol = new Molecule();
		Residue res = new Residue(rc.template.templateRes);
		res.fullName = pos.strand.mol.getResByPDBResNumber(pos.resNum).fullName;
		res.template = rc.template;
		res.markIntraResBondsByTemplate();
		res.molec = mol;
		res.indexInMolecule = mol.residues.size();
		mol.residues.add(res);
		mol.markInterResBonds();

		// make the residue DOFs and bounds
		VoxelShape.Rect voxel = new VoxelShape.Rect();
		List<DegreeOfFreedom> dofs = voxel.makeDihedralDOFs(res);
		ObjectiveFunction.DofBounds dofBounds;
		if (rc.rotamerIndex != null) {
			dofBounds = voxel.makeDihedralBounds(rc.template, rc.rotamerIndex);
		} else {
			dofBounds = new ObjectiveFunction.DofBounds(0);
		}

		// center all the dofs
		for (int d=0; d<dofs.size(); d++) {
			dofs.get(d).apply(dofBounds.getCenter(d));
		}

		return new ParametricMolecule(mol, dofs, dofBounds);
	}

	private static class BoundedDof {

		final DegreeOfFreedom dof;
		final double min;
		final double max;

		public BoundedDof(DegreeOfFreedom dof, double min, double max) {
			this.dof = dof;
			this.min = min;
			this.max = max;
		}
	}

	private static ParametricMolecule wipeDofs(ParametricMolecule pmol) {
		return new ParametricMolecule(
			pmol.mol,
			Collections.emptyList(),
			new ObjectiveFunction.DofBounds(new DoubleMatrix1D[] {
				DoubleFactory1D.dense.make(0),
				DoubleFactory1D.dense.make(0)
			})
		);
	}

	private static ParametricMolecule addDofs(ParametricMolecule pmol, List<BoundedDof> newDofs) {

		// expand the dofs
		List<DegreeOfFreedom> combinedDofs = new ArrayList<>(pmol.dofs);
		for (BoundedDof bdof : newDofs) {
			combinedDofs.add(bdof.dof);
		}

		// expand the dof bounds
		DoubleMatrix1D mins = DoubleFactory1D.dense.make(pmol.dofs.size() + newDofs.size());
		DoubleMatrix1D maxs = mins.copy();
		for (int i=0; i<pmol.dofs.size(); i++) {
			mins.set(i, pmol.dofBounds.getMin(i));
			maxs.set(i, pmol.dofBounds.getMax(i));
		}
		for (int i=0; i<newDofs.size(); i++) {
			mins.set(i + pmol.dofs.size(), newDofs.get(i).min);
			maxs.set(i + pmol.dofs.size(), newDofs.get(i).max);
		}

		// build a new pmol
		return new ParametricMolecule(
			pmol.mol,
			combinedDofs,
			new ObjectiveFunction.DofBounds(new DoubleMatrix1D[] { mins, maxs })
		);
	}

	private static String degreesToString(double[] degrees) {

		if (degrees == null) {
			return "null";
		}

		StringBuilder buf = new StringBuilder();
		buf.append("[");
		for (int i=0; i<degrees.length; i++) {
			if (i > 0) {
				buf.append(",");
			}
			buf.append(String.format("%6.1f", degrees[i]));
		}
		buf.append("]");
		return buf.toString();
	}

	private static class DihedralDof extends DegreeOfFreedom {

		public final Residue res;

		private final int a;
		private final int b;
		private final int c;
		private final int d;
		private final int[] r;

		public DihedralDof(Residue res, String a, String b, String c, String d, String ... r) {
			this.res = res;
			this.a = res.getAtomByName(a).indexInRes;
			this.b = res.getAtomByName(b).indexInRes;
			this.c = res.getAtomByName(c).indexInRes;
			this.d = res.getAtomByName(d).indexInRes;
			this.r = Arrays.stream(r)
				.mapToInt(atomName -> res.getAtomByName(atomName).indexInRes)
				.toArray();
		}

		public double measure() {
			return Protractor.measureDihedral(res.coords, a, b, c, d);
		}

		@Override
		public void apply(double angleDegrees) {

			// compute the target dihedral
			double angleRadians = Math.toRadians(angleDegrees);
			double sin = Math.sin(angleRadians);
			double cos = Math.cos(angleRadians);

			// measure the current dihedral
			double measuredSinCos[] = Protractor.measureDihedralSinCos(res.coords, a, b, c, d);

			// calc the dihedral rotation as a rigid body transformation relative to the current pose
			double dsin = sin*measuredSinCos[1] - cos*measuredSinCos[0];
			double dcos = cos*measuredSinCos[1] + sin*measuredSinCos[0];

			double[] bcoords = new double[] {
				res.coords[b*3    ],
				res.coords[b*3 + 1],
				res.coords[b*3 + 2],
			};
			double[] ccoords = new double[] {
				res.coords[c*3    ],
				res.coords[c*3 + 1],
				res.coords[c*3 + 2],
			};
			RigidBodyMotion dihRotation = new DihedralRotation(bcoords, ccoords, dsin, dcos);

			// rotate all the downstream atoms
			for (int r : this.r) {
				dihRotation.transform(res.coords, r);
			}
		}

		@Override
		public DOFBlock getBlock() {
			return null;
		}

		@Override
		public String getName() {
			return "dihedral";
		}
	}

	private static class MethylRotation extends DegreeOfFreedom {

		public final Residue res;

		private final int a;
		private final int b;
		private final int C;
		private final int H1;
		private final int H2;
		private final int H3;

		public MethylRotation(Residue res, String a, String b, String C, String H1, String H2, String H3) {
			this.res = res;
			this.a = res.getAtomByName(a).indexInRes;
			this.b = res.getAtomByName(b).indexInRes;
			this.C = res.getAtomByName(C).indexInRes;
			this.H1 = res.getAtomByName(H1).indexInRes;
			this.H2 = res.getAtomByName(H2).indexInRes;
			this.H3 = res.getAtomByName(H3).indexInRes;
		}

		public double measure() {
			return Protractor.measureDihedral(res.coords, a, b, C, H1);
		}

		@Override
		public void apply(double angleDegrees) {

			// compute the target dihedral
			double angleRadians = Math.toRadians(angleDegrees);
			double sin = Math.sin(angleRadians);
			double cos = Math.cos(angleRadians);

			// measure the current dihedral
			double measuredSinCos[] = Protractor.measureDihedralSinCos(res.coords, a, b, C, H1);

			// calc the dihedral rotation as a rigid body transformation relative to the current pose
			double dsin = sin*measuredSinCos[1] - cos*measuredSinCos[0];
			double dcos = cos*measuredSinCos[1] + sin*measuredSinCos[0];

			double[] bcoords = new double[] {
				res.coords[b*3    ],
				res.coords[b*3 + 1],
				res.coords[b*3 + 2],
			};
			double[] Ccoords = new double[] {
				res.coords[C*3    ],
				res.coords[C*3 + 1],
				res.coords[C*3 + 2],
			};
			RigidBodyMotion dihRotation = new DihedralRotation(bcoords, Ccoords, dsin, dcos);

			// rotate all the H atoms
			dihRotation.transform(res.coords, H1);
			dihRotation.transform(res.coords, H2);
			dihRotation.transform(res.coords, H3);
		}

		@Override
		public DOFBlock getBlock() {
			return null;
		}

		@Override
		public String getName() {
			return "methyl";
		}
	}

	private static double calcVDW(Atom a1, Atom a2, ForcefieldParams ffparams) {

		// get the radius
		// read atom coords
		double x1 = a1.res.coords[a1.indexInRes*3];
		double y1 = a1.res.coords[a1.indexInRes*3 + 1];
		double z1 = a1.res.coords[a1.indexInRes*3 + 2];
		double x2 = a2.res.coords[a2.indexInRes*3];
		double y2 = a2.res.coords[a2.indexInRes*3 + 1];
		double z2 = a2.res.coords[a2.indexInRes*3 + 2];

		// compute r2
		double dx = x1 - x2;
		double dy = y1 - y2;
		double dz = z1 - z2;
		double r2 = dx*dx + dy*dy + dz*dz;

		// get vdw params for this atom pair
		AtomNeighbors.Type type = new AtomNeighbors(a1).classifyAtom(a2);
		ForcefieldParams.VdwParams vdw = new ForcefieldParams.VdwParams();
		ffparams.getVdwParams(a1, a2, type, vdw);

		// compute vdw
		double r6 = r2*r2*r2;
		double r12 = r6*r6;
		return vdw.Aij/r12 - vdw.Bij/r6;
	}

	public static void top8000Dihedrals(String type)
	throws Exception {

		// clustering settings
		int densityWindowRadius = 2;
		int densityWindowCountThreshold = 20;
		double clusterDistThreshold = 50.0; // TODO: is this really a distance? it seems too high

		int numAngles = chiLib.get(type).size();
		assert (numAngles > 0);

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", true);
		log("%s dihedrals: %d", type, dihedralsByKey.size());

		VisIt.writeAngles2D(dihedralsByKey.values(), 0, 1, new File(type + ".dihedrals.vtk"));

		// make a histogram
		log("building histogram");
		DegreesHistogram hist = new DegreesHistogram(numAngles);
		for (double[] dihedrals : dihedralsByKey.values()) {
			hist.add(dihedrals);
		}
		log("histogram done");

		// filter the histogram
		{
			int count = hist.count();
			log("before filtering: %d", count);
			hist.filterDensityWindow(densityWindowRadius, densityWindowCountThreshold);
			int keptCount = hist.count();
			log("after filtering:  %d  %.1f%%", keptCount, 100f*keptCount/count);
		}

		// convert the histogram back to dihedrals
		List<double[]> keptDihedrals = new ArrayList<>();
		for (long key : hist.buckets.keySet()) {
			keptDihedrals.add(hist.makeDihedrals(key));
		}

		VisIt.writeAngles2D(keptDihedrals, 0, 1, new File(type + ".keptDihedrals.vtk"));

		// cluster the points
		List<List<double[]>> clusters = AngleClustering.cluster(keptDihedrals, numAngles, clusterDistThreshold);

		//VisIt.writeAngles2D(clusters.get(5), 0, 1, new File(type + ".cluster.vtk"));

		// calc bounding voxels for the clusters
		List<SmallAngleVoxel> voxels = clusters.stream()
			.map(cluster -> AngleClustering.calcVoxel(cluster))
			.collect(Collectors.toList());

		VisIt.writeVoxels(voxels, 0, 1, new File(type + ".voxels.vtk"));

		// fix fixed-size voxels to the voxels to maximize conf coverage
		List<SmallAngleVoxel> fixedVoxels = new ArrayList<>();
		for (int i=0; i<voxels.size(); i++) {
			SmallAngleVoxel voxel = voxels.get(i);

			// get the points in the voxel
			List<double[]> points = voxel.filter(dihedralsByKey.values());
			log("\nvoxel %2d   %6d side chains  %.1f%% of total:\n%s",
				i,
				points.size(),
				100f*points.size()/dihedralsByKey.size(),
				voxel
			);

			// fit the fixed-size voxel
			SmallAngleVoxel fixedVoxel = AngleClustering.fitFixedVoxel(points, voxel, 9);
			fixedVoxels.add(fixedVoxel);

			int count = fixedVoxel.filter(points).size();
			log("%6d selected   %.1f%% of cluster   %.1f%% of total:\n%s",
				count,
				100f*count/points.size(),
				100f*count/dihedralsByKey.size(),
				fixedVoxel
			);
		}

		VisIt.writeVoxels(fixedVoxels, 0, 1, new File(type + ".fixedVoxels.vtk"));
	}

	public static void lovellRotamers(String type)
	throws Exception {

		File dir = new File("/home/jeff/dlab/top8000");
		List<String> types = Arrays.asList(type.toUpperCase());

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", false);

		// analyze the Lovell rotamers for leucine
		ArrayList<SmallAngleVoxel> voxels = new ArrayList<>();
		double totalPercent = 0.0;
		for (Rotamer rot : rotamers.get(type)) {

			voxels.add(rot.voxel);

			// get dihedral counts
			long count = dihedralsByKey.values().stream()
				.filter(p -> rot.voxel.contains(p))
				.count();
			double percent = 100f*count/dihedralsByKey.size();
			totalPercent += percent;

			log("%s %8s: %6d   %4.1f%% of total", type, rot.name, count, percent);
		}

		log("total percent: %.1f%%", totalPercent);

		VisIt.writeVoxels(voxels, 0, 1, new File(type + ".lovellVoxels.vtk"));

		// analyze the voxel stats
		for (Rotamer rot : rotamers.get(type)) {

			Predicate<Residue> filter = res -> {
				double[] dihedrals = chiLib.measure(res, type);
				return dihedrals != null && rot.voxel.contains(dihedrals);

			};

			log("\n\n%s %8s", type, rot.name);

			Map<String,List<TemplateChooser.MeasuredRes>> measurementsByType = TemplateChooser.measureResidues(dir, types, filter);
			Map<String,double[]> modes = TemplateChooser.calcModes(measurementsByType, types);
		}
	}

	public static void top8000Angles(String type)
	throws Exception {

		// what does our template have?
		for (ResidueTemplate template : templateLib.templates) {
			if (template.name.equalsIgnoreCase(type) && template.templateRes.coords != null) {

				double[] angles = angleLib.measure(template.templateRes, type);
				log("template %s angles:", template.name);
				for (int d=0; d<angles.length; d++) {
					log("%20s = %.1f", angleLib.get(type).get(d).name, angles[d]);
				}
			}
		}

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", false);
		Map<ResKey,double[]> anglesByKey = readAngles(type, angleLib, type + ".angles.dat", false);

		// get the angle distributions for each rotamer
		Map<Rotamer,List<double[]>> anglesByRot = new HashMap<>();
		for (ResKey key : anglesByKey.keySet()) {

			// get the dihedrals for this res, if any
			double[] dihedrals = dihedralsByKey.get(key);
			if (dihedrals == null) {
				continue;
			}

			// what rotamer is this?
			Rotamer rot = rotamers.find(type, dihedrals);

			double[] angles = anglesByKey.get(key);

			anglesByRot.computeIfAbsent(rot, (r) -> new ArrayList<>()).add(angles);
		}

		List<double[]> allAngles = anglesByRot.values().stream().flatMap(it -> it.stream()).collect(Collectors.toList());
		VisIt.writeAngles2D(allAngles, 0, 1, new File(type + ".angles.vtk"));

		log("angles histograms:\n%s", new DegreesHistogram(allAngles).dump());

		// get detailed stats on each angle
		List<MeasurementLibrary.Measurement> libAngles = angleLib.get(type);
		for (int a=0; a<libAngles.size(); a++) {

			ClusterS1 cluster = new ClusterS1();
			for (double[] angles : allAngles) {
				cluster.add(angles[a]);
			}

			// TEMP
			//log("angle %s: %s", libAngles.get(a).name, cluster.new Stats());
		}

		/*
		for (Rotamer rot : anglesByRot.keySet()) {
			VisIt.writeAngles2D(anglesByRot.get(rot), 0, 1, new File(type + ".angles." + rot.name + ".vtk"));
		}
		*/
	}

	public static void top8000Tetrahedrals(String type)
	throws Exception {

		// what does our template have?
		for (ResidueTemplate template : templateLib.templates) {
			if (template.name.equalsIgnoreCase(type) && template.templateRes.coords != null) {

				double[] angles = tetraLib.measure(template.templateRes, type);
				log("template %s angles:", template.name);
				for (int d=0; d<angles.length; d++) {
					log("%20s = %.1f", tetraLib.get(type).get(d).name, angles[d]);
				}
			}
		}

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", false);
		Map<ResKey,double[]> tetrasByKey = readAngles(type, tetraLib, type + ".tetras.dat", true);

		// get the angle distributions for each rotamer
		Map<Rotamer,List<double[]>> anglesByRot = new HashMap<>();
		for (ResKey key : tetrasByKey.keySet()) {

			// get the dihedrals for this res, if any
			double[] dihedrals = dihedralsByKey.get(key);
			if (dihedrals == null) {
				continue;
			}

			// what rotamer is this?
			Rotamer rot = rotamers.find(type, dihedrals);

			double[] angles = tetrasByKey.get(key);

			anglesByRot.computeIfAbsent(rot, (r) -> new ArrayList<>()).add(angles);
		}

		List<double[]> allTetras = anglesByRot.values().stream().flatMap(it -> it.stream()).collect(Collectors.toList());
		VisIt.writeAngles2D(allTetras, 0, 1, new File(type + ".tetras.vtk"));

		log("angles histograms:\n%s", new DegreesHistogram(allTetras).dump());

		// get detailed stats on each angle
		List<MeasurementLibrary.Measurement> libTetras = tetraLib.get(type);
		for (int a=0; a<libTetras.size(); a++) {

			ClusterS1 cluster = new ClusterS1();
			for (double[] tetras : allTetras) {
				cluster.add(tetras[a]);
			}

			// TEMP
			//log("tetra %s: %s", libTetras.get(a).name, cluster.new Stats());
		}

		/*
		for (Rotamer rot : anglesByRot.keySet()) {
			VisIt.writeAngles2D(anglesByRot.get(rot), 0, 1, new File(type + ".angles." + rot.name + ".vtk"));
		}
		*/
	}

	public static void top8000Methyls(String type)
	throws Exception {

		// NOTE: reduce apparently puts all methyls so H1 is at 180 degrees! (ie, mod 3 is 60 degrees)
		// so there's really no information about methyl orientation in the structure at all!

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", false);
		Map<ResKey,double[]> anglesByKey = readAngles(type, angleLib, type + ".angles.dat", false);
		Map<ResKey,double[]> methylsByKey = readAngles(type, methylLib, type + ".methyls.dat", true);

		// what does our template have?
		for (ResidueTemplate template : templateLib.templates) {
			if (template.name.equalsIgnoreCase(type) && template.templateRes.coords != null) {

				double[] dihedrals = chiLib.measure(template.templateRes, type);
				Rotamer rot = rotamers.find(type, dihedrals);

				double[] methyls = methylLib.measure(template.templateRes, type);
				log("template rot: %s methyls: %.1f, %.1f", rot.name, methyls[0], methyls[1]);
			}
		}

		// get the methyl dihedral distributions for each rotamer
		Map<Rotamer,List<double[]>> methylsByRot = new HashMap<>();
		for (ResKey key : anglesByKey.keySet()) {

			// get the dihedrals for this res, if any
			double[] dihedrals = dihedralsByKey.get(key);
			if (dihedrals == null) {
				continue;
			}

			// what rotamer is this?
			Rotamer rot = rotamers.find(type, dihedrals);

			double[] methyls = methylsByKey.get(key);

			methylsByRot.computeIfAbsent(rot, (r) -> new ArrayList<>()).add(methyls);
		}

		List<double[]> allMethyls = methylsByRot.values().stream().flatMap(it -> it.stream()).collect(Collectors.toList());
		VisIt.writeAngles2D(allMethyls, 0, 1, new File(type + ".methyls.vtk"));

		for (Rotamer rot : methylsByRot.keySet()) {
			VisIt.writeAngles2D(methylsByRot.get(rot), 0, 1, new File(String.format(type + ".methyls.%s.vtk", rot.name)));
		}
	}

	public static void top8000Clashes() {

		Probe probe = new Probe();
		probe.matchTemplates(templateLib);

		logf("building atom connectivity...");
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(templateLib)
			.set15HasNonBonded(false)
			.build();
		log(" done!");

		List<Double> clashes = new ArrayList<>();

		scanner.scan((file, mol) -> {

			// analyze each residue
			for (Residue res : mol.residues) {

				// TEMP: just leucines for now
				String type = res.getType();
				if (!type.equals("LEU")) {
					continue;
				}

				// use a custom template assigner for now
				res.assignTemplateSimple(templateLib);
				if (res.template == null) {
					//log("WARN: no template assigned for %s %s", file.getName(), res.fullName);
					continue;
				}

				clashes.add(getProbeViolationsAllRotamers(probe, connectivity, res, type));
			}
		});

		log("results: %s", clashes.stream().mapToDouble(d -> d).summaryStatistics());
		log("histogram:\n%s", new Histogram(clashes, 20).toString(5, 3, 20));

		//VisIt.writeAngles2D(clashingDihedrals, 0, 1, new File("leu.clashingDihedrals.vtk"));
	}

	private static double getProbeViolationsAllRotamers(Probe probe, AtomConnectivity connectivity, Residue res, String type) {

		// get the dihedrals for this res, if any
		double[] dihedrals = chiLib.measure(res, type);
		assert (dihedrals != null);

		List<Double> violations = new ArrayList<>();
		Runnable collectViolations = () -> {
			probe.getInteractions(res, connectivity).stream()
				.mapToDouble(inter -> inter.getViolation(0.0))
				.filter(v -> v > 0.0)
				.forEach(v -> violations.add(v));
		};

		// for each rotamer
		List<Rotamer> rots = rotamers.get(type);
		if (rots.isEmpty()) {

			// no rotamer, eg gly, pro, ala
			collectViolations.run();

		} else {
			for (Rotamer rot : rotamers.get(type)) {

				// center this conformation on the rotamer voxel
				// NOTE: assumes rotamer voxel dihedral angles have same order as template dihedral angles
				for (int d=0; d<rot.voxel.intervals.length; d++) {
					new FreeDihedral(res, d).apply(rot.voxel.intervals[d].center);
				}

				collectViolations.run();
			}
		}

		// return the RMS violation
		if (violations.isEmpty()) {
			return 0.0;
		}
		double sum = 0.0;
		for (double v : violations) {
			sum += v*v;
		}
		return Math.sqrt(sum/violations.size());
	}

	private static Map<String,ResidueTemplate> pickTemplates(List<String> resTypes, ResidueTemplateLibrary templateLib) {
		Map<String,ResidueTemplate> out = new HashMap<>();
		for (String resType : resTypes) {
			out.put(resType, templateLib.getTemplate(resType));
		}
		return out;
	}

	private static void checkTemplates() {

		List<String> resTypes = Arrays.asList(
			"ARG", "HIP", "HID", "HIE", "LYS", "ASP", "GLU",
			"CYS", "GLY", "PRO",
			"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP",
			"SER", "THR", "ASN", "GLN"
		);

		// make the various template libs
		TemplateType.Map<Map<String,ResidueTemplate>> templateLibs = TemplateType.mapOf(
			pickTemplates(
				resTypes,
				new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
					.build()
			),
			pickTemplates(
				resTypes,
				new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
					.clearTemplateCoords()
					.addTemplateCoords(FileTools.readFile("template coords.txt"))
					.build()
			),
			pickTemplates(
				resTypes,
				new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
					.clearTemplateCoords()
					.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
					.build()
			)
		);

		for (TemplateType templateType : TemplateType.values()) {

			Map<String,ResidueTemplate> templates = templateLibs.get(templateType);
			if (templates == null) {
				continue;
			}
			log("%s", templateType);

			Probe probe = new Probe();
			probe.matchTemplates(templates.values());
			AtomConnectivity connectivity = new AtomConnectivity.Builder()
				.addTemplates(templates.values())
				.set15HasNonBonded(false)
				.build();

			for (String resType : resTypes) {
				ResidueTemplate templ = templates.get(resType);

				log("\t%s: %s", resType, templ);

				// make a molecule from just the template
				Molecule mol = new Molecule();
				Residue res = new Residue(templ.templateRes);
				res.fullName = resType + " A   1";
				res.template = templ;
				res.markIntraResBondsByTemplate();
				mol.residues.add(res);
				res.indexInMolecule = 0;

				// make the dihedral dofs
				VoxelShape.Rect voxel = new VoxelShape.Rect();
				List<DegreeOfFreedom> dofs = voxel.makeDihedralDOFs(res);

				// get all the dof bounds
				List<ObjectiveFunction.DofBounds> boundsByRot = new ArrayList<>();
				if (templ.getNumRotamers() > 0) {
					for (int i=0; i<templ.getNumRotamers(); i++) {
						boundsByRot.add(voxel.makeDihedralBounds(templ, i));
					}
				} else {
					boundsByRot.add(new ObjectiveFunction.DofBounds(0));
				}

				// for each rotamer
				for (ObjectiveFunction.DofBounds dofBounds : boundsByRot) {

					// set the rotameric pose
					for (int d=0; d<dofs.size(); d++) {
						dofs.get(d).apply(dofBounds.getCenter(d));
					}

					// get the rotamer name
					String rotName = null;
					Double rotPercent = null;
					Rotamer rot = rotamers.find(resType, chiLib.measure(res, resType));
					if (rot != null) {
						rotName = rot.name;
						rotPercent = rot.percent;
					}

					log("\t\trotamer %8s: %4.1f%%  %s", rotName, rotPercent, dofBounds.toString(4, 0));

					// check for probe clashes
					probe.getInteractions(res, connectivity).stream()
						.filter(inter -> inter.getViolation(0.0) > 0.0)
						.forEach(inter ->
							log("\t\t\t%s %s", inter.atomPair, inter)
						);
				}
			}
		}
	}

	private static Map<ResKey,double[]> readAngles(String type, MeasurementLibrary lib, String filename, boolean recalc)
	throws Exception {

		if (!recalc) {
			logf("load");

			//noinspection unchecked
			HashMap<ResKey,double[]> anglesByKey = ObjectIO.read(filename, HashMap.class);

			log("ed %d angles from %s", anglesByKey.size(), filename);

			return anglesByKey;
		}

		Map<ResKey,double[]> anglesByKey = new HashMap<>();

		scanner.scan((file, mol) -> {

			// analyze each residue
			for (Residue res : mol.residues) {

				// filter by type
				if (!res.getType().equalsIgnoreCase(type)) {
					continue;
				}

				ResKey key = new ResKey(type, file.getName(), res.getPDBResNumber());

				// measure the angles
				double[] angles = lib.measure(res, type);
				if (angles == null) {
					continue;
				}

				// make sure we got good angles
				boolean isGood = true;
				for (double a : angles) {
					if (!Double.isFinite(a)) {
						isGood = false;
						break;
					}
				}
				if (!isGood) {
					log("warning: bad angles for %s, ignoring", key);
					continue;
				}

				/* TEMP: make sure template assignment doesn't change any angles
				res.assignTemplate(templateLib);

				// the two measurements should be exactly identical
				assert (Arrays.equals(angles, lib.measure(res))) :
					String.format("%s\n\texpected: %s\n\tobserved: %s",
						key, Arrays.toString(angles), Arrays.toString(lib.measure(res))
					);
				*/

				anglesByKey.put(key, angles);
			}
		});

		// write all the dihedrals to a file
		ObjectIO.write(anglesByKey, filename);
		log("wrote %d angles to %s", anglesByKey.size(), filename);

		return anglesByKey;
	}

	static class ResKey implements Serializable {

		private static final long serialVersionUID = 7858382121141521600L;

		public final String type;
		public final String filename;
		public final String resNum;

		public ResKey(String type, String filename, String resNum) {
			this.type = type;
			this.filename = filename;
			this.resNum = resNum;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				type.hashCode(),
				filename.hashCode(),
				resNum.hashCode()
			);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof ResKey && equals((ResKey)other);
		}

		public boolean equals(ResKey other) {
			return this.type.equals(other.type)
				&& this.filename.equals(other.filename)
				&& this.resNum.equals(other.resNum);
		}

		@Override
		public String toString() {
			return String.format("%s:%s:%s", type, filename, resNum);
		}
	}

	public static void energyLandscape(String type) {

		TemplateType.Map<ResidueTemplateLibrary> templateLibs = TemplateType.mapOf(
			new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
				.build(),
			new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
				.clearTemplateCoords()
				.addTemplateCoords(FileTools.readFile("template coords.txt"))
				.build(),
			new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
				.clearTemplateCoords()
				.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
				.build()
		);

		// make a mechanism to change template coords on-the-fly
		BiConsumer<SimpleConfSpace.ResidueConf,TemplateType> setTemplateCoords = (rc, templateType) -> {
			System.arraycopy(
				templateLibs.get(templateType).getTemplate(rc.template.name, true).templateRes.coords, 0,
				rc.template.templateRes.coords, 0,
				rc.template.templateRes.coords.length
			);
		};

		// make a conspace with all the rotamers for this type
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		String resNum = "A37";
		strand.flexibility.get(resNum).setLibraryRotamers(type).setContinuous(); // asn
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// pick just intra interactions
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntra(0)
			.make();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build()) {

			// pick the first rc, doesn't matter
			SimpleConfSpace.ResidueConf rc = confSpace.positions.get(0).resConfs.get(0);

			// set the template coords
			TemplateType templateType = TemplateType.Idealized;
			setTemplateCoords.accept(rc, templateType);

			// grab the first two chi angles from the dofs
			ParametricMolecule pmol = makeMolecule(confSpace, rc);
			assert (pmol.dofs.size() >= 2);
			FreeDihedral[] chi = new FreeDihedral[] {
				(FreeDihedral)pmol.dofs.get(0),
				(FreeDihedral)pmol.dofs.get(1)
			};

			// wipe all the dofs
			pmol = wipeDofs(pmol);

			// set the backbone groups to be out of the way
			Residue res = pmol.mol.getResByPDBResNumber(resNum);
			new DihedralDof(res, "C", "CA", "N", "H", "H").apply(0);
			new DihedralDof(res, "N", "CA", "C", "O", "O").apply(0);

			// add methyl rotations
			assert (type.equalsIgnoreCase("LEU"));
			pmol = addDofs(pmol, Arrays.asList(
				new BoundedDof(new MethylRotation(res, "CB", "CG", "CD1", "HD11", "HD12", "HD13"), 0, 120),
				new BoundedDof(new MethylRotation(res, "CB", "CG", "CD2", "HD21", "HD22", "HD23"), 0, 120)
			));

			// allocate space for the energies
			double[][] energies = new double[360][360];
			for (int i=0; i<360; i++) {
				Arrays.fill(energies[i], Double.NaN);
			}

			// sweep over chi space and sample the energy function
			Progress progress = new Progress(360*360);
			for (int[] indices : new MathTools.GridIterable(new int[] { 360, 360 })) {

				// set dihedral angles
				chi[0].apply(indices[0]);
				chi[1].apply(indices[1]);

				// calc the energy here
				double energy = ecalc.calcEnergy(pmol, inters).energy;

				// store results in y-major order
				int ix = indices[0];
				int iy = indices[1];
				energies[iy][ix] = energy;

				progress.incrementProgress();
			}

			// sample the axes
			double[] xaxis = new double[360];
			double[] yaxis = new double[360];
			for (int i=0; i<360; i++) {
				xaxis[i] = i;
				yaxis[i] = i;
			}

			VisIt.writeGrid2D(xaxis, yaxis, energies, new File(type + ".dihedralEnergy." + templateType.name().toLowerCase() + ".minbb.methyls.vtk"));
		}
	}

	private static class BondAngleRotation extends DegreeOfFreedom {

		public final Residue res;

		private final int a;
		private final int b;
		private final int c;
		private final List<Integer> rotatedIndices;

		public BondAngleRotation(Residue res, String a, String b, String c, List<String> rotatedAtoms) {
			this.res = res;
			this.a = res.getAtomByName(a).indexInRes;
			this.b = res.getAtomByName(b).indexInRes;
			this.c = res.getAtomByName(c).indexInRes;
			this.rotatedIndices = rotatedAtoms.stream()
				.map(atomName -> res.getAtomByName(atomName).indexInRes)
				.collect(Collectors.toList());
		}

		@Override
		public void apply(double angleDegrees) {

			// measure the current angle
			double measuredAngleDegrees = Protractor.measureBondAngle(res.coords, a, b, c);

			int a3 = a*3;
			int b3 = b*3;
			int c3 = c*3;
			double[] coords = res.coords;

			// the center is the b atom
			double[] center = new double[] {
				coords[b3    ],
				coords[b3 + 1],
				coords[b3 + 2]
			};

			// the rotation axis is ba x bc
			// ba = a - b
			double[] ba = new double[] {
				coords[a3    ] - coords[b3    ],
				coords[a3 + 1] - coords[b3 + 1],
				coords[a3 + 2] - coords[b3 + 2]
			};
			// bc = c - b
			double[] bc = new double[] {
				coords[c3    ] - coords[b3    ],
				coords[c3 + 1] - coords[b3 + 1],
				coords[c3 + 2] - coords[b3 + 2]
			};
			double[] axis = VectorAlgebra.cross(ba, bc);

			// the angle is easy
			double theta = Protractor.getDeltaDegrees(measuredAngleDegrees, angleDegrees);

			// rotate all the H atoms
			RigidBodyMotion rotation = new RigidBodyMotion(center, axis, theta, false);
			for (int i : rotatedIndices) {
				rotation.transform(res.coords, i);
			}
		}

		@Override
		public DOFBlock getBlock() {
			return null;
		}

		@Override
		public String getName() {
			return "methyl";
		}
	}

	private static class TetrahedralInPlaneRotation extends DegreeOfFreedom {

		public final Residue res;

		private final int a;
		private final int b;
		private final int c;
		private final int d;
		private final List<Integer> rotatedIndices;

		private Protractor.TetrahedralGeometry tetra = new Protractor.TetrahedralGeometry();

		public TetrahedralInPlaneRotation(Residue res, String a, String b, String c, String d, List<String> rotatedAtoms) {
			this.res = res;
			this.a = res.getAtomByName(a).indexInRes;
			this.b = res.getAtomByName(b).indexInRes;
			this.c = res.getAtomByName(c).indexInRes;
			this.d = res.getAtomByName(d).indexInRes;
			this.rotatedIndices = rotatedAtoms.stream()
				.map(atomName -> res.getAtomByName(atomName).indexInRes)
				.collect(Collectors.toList());
		}

		public double measure() {
			tetra.update(res.coords, a, b, c, d);
			return tetra.inPlaneDegrees;
		}

		@Override
		public void apply(double angleDegrees) {

			tetra.update(res.coords, a, b, c, d);

			// the angle is easy
			double theta = Protractor.getDeltaDegrees(tetra.inPlaneDegrees, angleDegrees);

			// rotate all the downstream atoms
			RigidBodyMotion rotation = new RigidBodyMotion(tetra.center, tetra.outOfPlaneAxis, theta, false);
			for (int i : rotatedIndices) {
				rotation.transform(res.coords, i);
			}
		}

		@Override
		public DOFBlock getBlock() {
			return null;
		}

		@Override
		public String getName() {
			return "Tetra-IP";
		}
	}

	private static class TetrahedralOutOfPlaneRotation extends DegreeOfFreedom {

		public final Residue res;

		private final int a;
		private final int b;
		private final int c;
		private final int d;
		private final List<Integer> rotatedIndices;

		private Protractor.TetrahedralGeometry tetra = new Protractor.TetrahedralGeometry();

		public TetrahedralOutOfPlaneRotation(Residue res, String a, String b, String c, String d, List<String> rotatedAtoms) {
			this.res = res;
			this.a = res.getAtomByName(a).indexInRes;
			this.b = res.getAtomByName(b).indexInRes;
			this.c = res.getAtomByName(c).indexInRes;
			this.d = res.getAtomByName(d).indexInRes;
			this.rotatedIndices = rotatedAtoms.stream()
				.map(atomName -> res.getAtomByName(atomName).indexInRes)
				.collect(Collectors.toList());
		}

		public double measure() {
			tetra.update(res.coords, a, b, c, d);
			return tetra.outOfPlaneDegrees;
		}

		@Override
		public void apply(double angleDegrees) {

			tetra.update(res.coords, a, b, c, d);

			// the angle is easy
			double theta = Protractor.getDeltaDegrees(tetra.outOfPlaneDegrees, angleDegrees);

			// get the rotation axis
			double[] axis = VectorAlgebra.cross(tetra.bd, tetra.outOfPlaneAxis);

			// rotate all the downstream atoms
			RigidBodyMotion rotation = new RigidBodyMotion(tetra.center, axis, theta, false);
			for (int i : rotatedIndices) {
				rotation.transform(res.coords, i);
			}
		}

		@Override
		public DOFBlock getBlock() {
			return null;
		}

		@Override
		public String getName() {
			return "Tetra-OOP";
		}
	}

	private static void top8000RotamerStats(String type)
	throws Exception {

		// clustering settings
		int densityWindowRadius = 2;
		int densityWindowCountThreshold = 20;
		double clusterDistThreshold = 50.0; // TODO: is this really a distance? it seems too high

		int numAngles = chiLib.get(type).size();
		assert (numAngles > 0);

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", true);
		log("%s dihedrals: %d", type, dihedralsByKey.size());

		VisIt.writeAngles2D(dihedralsByKey.values(), 0, 1, new File(type + ".dihedrals.vtk"));

		// make a histogram
		log("building histogram");
		DegreesHistogram hist = new DegreesHistogram(numAngles);
		for (double[] dihedrals : dihedralsByKey.values()) {
			hist.add(dihedrals);
		}
		log("histogram done");

		// filter the histogram
		{
			int count = hist.count();
			log("before filtering: %d", count);
			hist.filterDensityWindow(densityWindowRadius, densityWindowCountThreshold);
			int keptCount = hist.count();
			log("after filtering:  %d  %.1f%%", keptCount, 100f*keptCount/count);
		}

		// convert the histogram back to dihedrals
		List<double[]> keptDihedrals = new ArrayList<>();
		for (long key : hist.buckets.keySet()) {
			keptDihedrals.add(hist.makeDihedrals(key));
		}

		VisIt.writeAngles2D(keptDihedrals, 0, 1, new File(type + ".keptDihedrals.vtk"));

		// cluster the points
		List<List<double[]>> clusters = AngleClustering.cluster(keptDihedrals, numAngles, clusterDistThreshold);

		// calc bounding voxels for the clusters
		List<SmallAngleVoxel> voxels = clusters.stream()
			.map(cluster -> AngleClustering.calcVoxel(cluster))
			.collect(Collectors.toList());

		VisIt.writeVoxels(voxels, 0, 1, new File(type + ".voxels.vtk"));
	}
}
