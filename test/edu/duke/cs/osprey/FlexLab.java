package edu.duke.cs.osprey;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
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
import java.util.function.Consumer;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


public class FlexLab {

	public static void main(String[] args)
	throws Exception {
		//checkRotamerClashes();
		//top8000Dihedrals("leu");
		//top8000Angles("leu");
		//top8000Tetrahedrals("leu");
		//top8000Methyls();
		//top8000Clashes();
		//lovellRotamers("leu");
		energyLandscape("leu");
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
	private static AngleLibrary chiLib = new AngleLibrary();
	static {
		chiLib.add("LEU", new AngleLibrary.Dihedral("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("LEU", new AngleLibrary.Dihedral("chi2", "CA", "CB", "CG", "CD1"));

		chiLib.add("TRP", new AngleLibrary.Dihedral("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("TRP", new AngleLibrary.Dihedral("chi2", "CA", "CB", "CG", "CD1"));
	}

	// define all the methyl groups
	private static AngleLibrary methylLib = new AngleLibrary();
	static {
		methylLib.add("LEU", new AngleLibrary.DihedralMod3("methyl1", "CB", "CG", "CD1", "HD11"));
		methylLib.add("LEU", new AngleLibrary.DihedralMod3("methyl2", "CB", "CG", "CD2", "HD21"));
	}

	// define all the bond angles
	private static AngleLibrary angleLib = new AngleLibrary();
	static {
		angleLib.add("LEU", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("LEU", new AngleLibrary.BondAngle("C", "CA", "CB"));
		angleLib.add("LEU", new AngleLibrary.BondAngle("N", "CA", "C"));

		angleLib.add("LEU", new AngleLibrary.BondAngle("CA", "CB", "CG"));
		angleLib.add("LEU", new AngleLibrary.BondAngle("CB", "CG", "CD1"));
		angleLib.add("LEU", new AngleLibrary.BondAngle("CB", "CG", "CD2"));

		// aromatics
		angleLib.add("HIS", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("HIS", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("PHE", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("PHE", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("TYR", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("TYR", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("TRP", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("TRP", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		// charged
		angleLib.add("ARG", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("ARG", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("LYS", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("LYS", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("ASP", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("ASP", new AngleLibrary.BondAngle("CA", "CB", "CG"));

		angleLib.add("GLU", new AngleLibrary.BondAngle("N", "CA", "CB"));
		angleLib.add("GLU", new AngleLibrary.BondAngle("CA", "CB", "CG"));
	}

	// define all the tetrahedral angles
	private static AngleLibrary tetraLib = new AngleLibrary();
	static {
		tetraLib.add("LEU", new AngleLibrary.TetrahedralInPlaneAngle("N", "CA", "C", "CB"));
		tetraLib.add("LEU", new AngleLibrary.TetrahedralOutOfPlaneAngle("N", "CA", "C", "CB"));
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
			throw new NoSuchElementException("no matching " + type + " rotamer for " + Arrays.toString(dihedrals));
		}
	}

	// define a rotamer library based on voxels
	private static RotamerLibrary rotamers = new RotamerLibrary();
	static {

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
	}

	public static void checkRotamerClashes() {

		// make a conspace with all the rotamers
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A23").setLibraryRotamers("ARG", "HIS", "LYS", "ASP", "GLU", "CYS", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER", "THR", "ASN", "GLN").setContinuous(); // asn
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		Probe probe = new Probe();
		probe.matchTemplates(confSpace);

		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.set15HasNonBonded(false)
			.build();

		// pick just intra interactions
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntra(0)
			.make();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build()) {

			EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(ecalc)
				.setIsMinimizing(false)
				.build();

			// for each rotamer
			for (SimpleConfSpace.ResidueConf rc : confSpace.positions.get(0).resConfs) {

				// TEMP: only leucines
				if (!rc.template.name.equals("LEU")) {
					continue;
				}
				String type = "LEU";

				RCTuple tuple = new RCTuple(0, rc.index);
				log("%s %d - %3d", rc.template.name, rc.rotamerIndex, rc.index);

				Consumer<Molecule> analyzeClashes = (mol) -> {
					probe.getInteractions(mol.residues, inters, connectivity).stream()
						.filter(interaction -> interaction.contact.isClash)
						.sorted(Comparator.comparing(interaction -> interaction.getViolation(0.0)))
						.forEach(interaction -> {
							double vdw = calcVDW(interaction.atomPair.a, interaction.atomPair.b, ffparams);
							log("\t\t%s   %s     vdw=%.4f", interaction.atomPair, interaction, vdw);
						});
				};

				// measure dihedrals and bond angles
				{
					ParametricMolecule pmol = confSpace.makeMolecule(tuple);
					Residue res = pmol.mol.getResByPDBResNumber("A23");
					double[] dihedrals = chiLib.measure(res);
					Rotamer rot = rotamers.find(type, dihedrals);
					log("\trotamer:     %s  %.1f%%", rot.name, rot.percent);
					log("\tdihedrals:   %s", Arrays.toString(dihedrals));
					log("\tbond angles: %s", Arrays.toString(angleLib.measure(res)));
					log("\tmethyls:     %s", Arrays.toString(methylLib.measure(res)));
					log("\trigid: %.3f kcal/mol", rigidEcalc.calcEnergy(pmol, inters).energy);
					analyzeClashes.accept(pmol.mol);
				}

				// minimize dihedrals
				{
					EnergyCalculator.EnergiedParametricMolecule epmol = ecalc.calcEnergy(confSpace.makeMolecule(tuple), inters);
					log("\tminimzed dihedrals: %.3f kcal/mol", epmol.energy);
					Residue res = epmol.pmol.mol.getResByPDBResNumber("A23");
					log("\t\tdihedrals: %s", Arrays.toString(chiLib.measure(res)));
					analyzeClashes.accept(epmol.pmol.mol);
				}

				// minimize methyls
				{
					Molecule mol = confSpace.makeMolecule(tuple).mol;
					Residue res = mol.getResByPDBResNumber("A23");
					ParametricMolecule pmol = new ParametricMolecule(
						mol,
						Arrays.asList(
							new MethylRotation(res, "CB", "CG", "CD1", "HD11", "HD12", "HD13"),
							new MethylRotation(res, "CB", "CG", "CD2", "HD21", "HD22", "HD23")
						),
						new ObjectiveFunction.DofBounds(new DoubleMatrix1D[] {
							DoubleFactory1D.dense.make(new double[] { -60.0, -60.0 }),
							DoubleFactory1D.dense.make(new double[] { 60.0, 60.0 })
						})
					);
					EnergyCalculator.EnergiedParametricMolecule epmol = ecalc.calcEnergy(pmol, inters);
					log("\tminimzed methyls: %.3f kcal/mol", epmol.energy);
					log("\t\tmethyls: %s", Arrays.toString(methylLib.measure(res)));
					analyzeClashes.accept(epmol.pmol.mol);
				}

				// minimize dihedrals and methyls
				{
					ParametricMolecule pmol = confSpace.makeMolecule(tuple);
					assert (pmol.dofs.size() == 2);
					Residue res = pmol.mol.getResByPDBResNumber("A23");
					ParametricMolecule pmol2 = new ParametricMolecule(
						pmol.mol,
						Arrays.asList(
							pmol.dofs.get(0), // chi1
							pmol.dofs.get(1), // chi2
							new MethylRotation(res, "CB", "CG", "CD1", "HD11", "HD12", "HD13"),
							new MethylRotation(res, "CB", "CG", "CD2", "HD21", "HD22", "HD23")
						),
						new ObjectiveFunction.DofBounds(new DoubleMatrix1D[] {
							DoubleFactory1D.dense.make(new double[] { pmol.dofBounds.getMin(0), pmol.dofBounds.getMin(1), -60.0, -60.0 }),
							DoubleFactory1D.dense.make(new double[] { pmol.dofBounds.getMax(0), pmol.dofBounds.getMax(1), 60.0, 60.0 })
						})
					);
					EnergyCalculator.EnergiedParametricMolecule epmol = ecalc.calcEnergy(pmol2, inters);
					log("\tminimzed dihedrals and methyls: %.3f kcal/mol", epmol.energy);
					log("\t\tdihedrals: %s", Arrays.toString(chiLib.measure(res)));
					log("\t\tmethyls: %s", Arrays.toString(methylLib.measure(res)));
					analyzeClashes.accept(epmol.pmol.mol);
				}
			}

			/* TEMP
			// look at individual structures
			RCTuple tuple = new RCTuple(0, 99);
			ParametricMolecule pmol = confSpace.makeMolecule(tuple);
			PDBIO.writeFile(pmol.mol, new File("rotamer.pdb"));
			*/
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

		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", false);
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

		int numAngles = chiLib.get(type).size();
		Map<ResKey,double[]> dihedralsByKey = readAngles(type, chiLib, type + ".dihedrals.dat", false);

		// analyze the Lovell rotamers for leucine
		List<SmallAngleVoxel> lovellVoxels;
		if (type.equalsIgnoreCase("leu")) {
			lovellVoxels = Arrays.asList(
				new SmallAngleVoxel(new double[] { 62, 80 }),
				new SmallAngleVoxel(new double[] { -177, 65 }),
				new SmallAngleVoxel(new double[] { -172, 145 }),
				new SmallAngleVoxel(new double[] { -85, 65 }),
				new SmallAngleVoxel(new double[] { -65, 175 })
			);
		} else if (type.equalsIgnoreCase("trp")) {
			lovellVoxels = Arrays.asList(
				new SmallAngleVoxel(new double[] { 62, -90 }),
				new SmallAngleVoxel(new double[] { 62, 90 }),
				new SmallAngleVoxel(new double[] { -177, -105 }),
				new SmallAngleVoxel(new double[] { -177, 90 }),
				new SmallAngleVoxel(new double[] { -65, -90 }),
				new SmallAngleVoxel(new double[] { -65, -5 }),
				new SmallAngleVoxel(new double[] { -65, 95 })
			);
		} else {
			throw new IllegalArgumentException("unknown lovell voxels for " + type);
		}

		for (SmallAngleVoxel voxel : lovellVoxels) {
			for (int d=0; d<numAngles; d++) {
				voxel.intervals[d].less = -9;
				voxel.intervals[d].more = 9;
			}

			// get dihedral counts
			long count = dihedralsByKey.values().stream()
				.filter(p -> voxel.contains(p))
				.count();

			log("\nLovell voxel: %6d   %.1f%% of total\n%s",
				count,
				100f*count/dihedralsByKey.size(),
				voxel
			);
		}

		VisIt.writeVoxels(lovellVoxels, 0, 1, new File(type + ".lovellVoxels.vtk"));
	}

	public static void top8000Angles(String type)
	throws Exception {

		// what does our template have?
		for (ResidueTemplate template : templateLib.templates) {
			if (template.name.equalsIgnoreCase(type) && template.templateRes.coords != null) {

				double[] angles = angleLib.measure(template.templateRes);
				log("template %s angles:", template.name);
				for (int d=0; d<angles.length; d++) {
					log("%20s = %.1f", new ArrayList<>(angleLib.get(type).values()).get(d).name, angles[d]);
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
		List<AngleLibrary.Angle> libAngles = new ArrayList<>(angleLib.get(type).values());
		for (int a=0; a<libAngles.size(); a++) {

			SmallAngleCluster cluster = new SmallAngleCluster();
			for (double[] angles : allAngles) {
				cluster.add(angles[a]);
			}

			log("angle %s: %s", libAngles.get(a).name, cluster.new Stats());
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

				double[] angles = tetraLib.measure(template.templateRes);
				log("template %s angles:", template.name);
				for (int d=0; d<angles.length; d++) {
					log("%20s = %.1f", new ArrayList<>(tetraLib.get(type).values()).get(d).name, angles[d]);
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
		List<AngleLibrary.Angle> libTetras = new ArrayList<>(tetraLib.get(type).values());
		for (int a=0; a<libTetras.size(); a++) {

			SmallAngleCluster cluster = new SmallAngleCluster();
			for (double[] tetras : allTetras) {
				cluster.add(tetras[a]);
			}

			log("tetra %s: %s", libTetras.get(a).name, cluster.new Stats());
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

				double[] dihedrals = chiLib.measure(template.templateRes);
				Rotamer rot = rotamers.find(type, dihedrals);

				double[] methyls = methylLib.measure(template.templateRes);
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

		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(templateLib)
			.set15HasNonBonded(false)
			.build();

		List<double[]> clashingDihedrals = new ArrayList<>();

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

				//ResKey key = new ResKey(type, file.getName(), res.getPDBResNumber());

				// get the dihedrals for this res, if any
				double[] dihedrals = chiLib.measure(res);
				if (dihedrals == null) {
					continue;
				}

				/*
				// what rotamer is this?
				Rotamer rot = rotamers.find(dihedrals);

				// just look at the popular rotamers for now
				if (rot.percent < 50) {
					continue;
				}
				*/

				/*
				// is this atom pair clashing?
				Probe.AtomPair atomPair = probe.new AtomPair(
					res.getAtomByName("C"),
					res.getAtomByName("HD12")
				);
				Probe.AtomPair.Interaction interaction = atomPair.getInteraction();
				if (interaction.contact == Probe.Contact.BadClash) {
					double vdw = calcVDW(interaction.atomPair.a, interaction.atomPair.b, ffparams);
					log("\t%24s   %s   %s     vdw=%.4f", key, rot.name, interaction, vdw);
				}
				*/

				// is this example clashing?
				for (int[] atomPair : connectivity.getAtomPairs(res, res).getPairs(AtomNeighbors.Type.NONBONDED)) {

					Probe.AtomPair probePair = probe.new AtomPair(
						res.atoms.get(atomPair[0]),
						res.atoms.get(atomPair[1])
					);

					if (probePair.getInteraction().contact == Probe.Contact.BadClash) {
						clashingDihedrals.add(dihedrals);
						break;
					}
				}
			}
		});

		VisIt.writeAngles2D(clashingDihedrals, 0, 1, new File("leu.clashingDihedrals.vtk"));
	}

	private static Map<ResKey,double[]> readAngles(String type, AngleLibrary lib, String filename, boolean recalc)
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
				double[] angles = lib.measure(res);
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

		// make a conspace with all the rotamers for this type
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A23").setLibraryRotamers(type).setContinuous(); // asn
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// pick just intra interactions
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntra(0)
			.make();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build()) {

			EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(ecalc)
				.setIsMinimizing(false)
				.build();

			// sweep space and sample the energy function at every point
			RCTuple tuple = new RCTuple(0, 0);
			ParametricMolecule pmol = confSpace.makeMolecule(tuple);

			// look at the bond angles of the template
			Residue res = pmol.mol.residues.getOrThrow("A23");
			log("%s template angles: %s", type, Arrays.toString(angleLib.measure(res)));

			// make tetrahedral geometry dofs
			List<String> sideChain = Arrays.asList("CB", "CG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23");
			TetrahedralInPlaneRotation tetraIPCB = new TetrahedralInPlaneRotation(
				res,
				"N", "CA", "C", "CB",
				sideChain
			);
			TetrahedralOutOfPlaneRotation tetraOOPCB = new TetrahedralOutOfPlaneRotation(
				res,
				"N", "CA", "C", "CB",
				sideChain
			);

			// make bond angle dofs
			BondAngleRotation rotCACBCG = new BondAngleRotation(
				res,
				"CA", "CB", "CG",
				Arrays.asList("CG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23")
			);

			// make methyl dofs
			MethylRotation[] methyl = new MethylRotation[] {
				new MethylRotation(res, "CB", "CG", "CD1", "HD11", "HD12", "HD13"),
				new MethylRotation(res, "CB", "CG", "CD2", "HD21", "HD22", "HD23")
			};

			Consumer<String> dumpAngles = (label) ->
				log("%s:"
					+ "\n\ttetra IP:  %6.1f"
					+ "\n\ttetra OOP: %6.1f"
					+ "\n\tCA-CB-CG:  %6.1f"
					+ "\n\tmethyl1:   %6.1f"
					+ "\n\tmethyl2:   %6.1f",
					label,
					tetraLib.get(type).get("TetraIP-CB").measure(res),
					tetraLib.get(type).get("TetraOOP-CB").measure(res),
					angleLib.get(type).get("CA-CB-CG").measure(res),
					methylLib.get(type).get("methyl1").measure(res),
					methylLib.get(type).get("methyl2").measure(res)
				);
			dumpAngles.accept("stock angles");

			// default angles are -3.2, -53.9, 109.4, 60.0, 60.0
			//tetraIPCB.apply(-3.2);
			//tetraOOPCB.apply(-53.9);
			//rotCACBCG.apply(109.4);
			//methyl[0].apply(60.0);
			//methyl[1].apply(60.0);

			// change the bond angles to leu modals: -0.6, -52.4, 116.4
			tetraIPCB.apply(-0.6);
			tetraOOPCB.apply(-52.4);
			rotCACBCG.apply(116.4);

			dumpAngles.accept("adjusted angles");

			// sweep over chi angles
			assert (pmol.dofs.size() == 2);
			FreeDihedral[] chi = new FreeDihedral[] {
				(FreeDihedral)pmol.dofs.get(0),
				(FreeDihedral)pmol.dofs.get(1)
			};

			// allow minimizing methyls
			ParametricMolecule pmolMethyls = new ParametricMolecule(
				pmol.mol,
				Arrays.asList(methyl[0], methyl[1]),
				new ObjectiveFunction.DofBounds(new DoubleMatrix1D[] {
					DoubleFactory1D.dense.make(new double[] { 0, 0 }),
					DoubleFactory1D.dense.make(new double[] { 120, 120 })
				})
			);

			double[][] energies = new double[360][360];
			double[][] methylAngles1 = new double[360][360];
			double[][] methylAngles2 = new double[360][360];
			for (int i=0; i<360; i++) {
				Arrays.fill(energies[i], Double.NaN);
				Arrays.fill(methylAngles1[i], Double.NaN);
				Arrays.fill(methylAngles2[i], Double.NaN);
			}

			Progress progress = new Progress(360*360);
			for (int[] indices : new MathTools.GridIterable(new int[] { 360, 360 })) {

				// set dihedral angles
				chi[0].apply(indices[0]);
				chi[1].apply(indices[1]);

				// calc the rigid energy
				//double energy = rigidEcalc.calcEnergy(pmol, inters).energy;

				// calc the minimized energy
				double energy = ecalc.calcEnergy(pmolMethyls, inters).energy;

				// store results in y-major order
				int ix = indices[0];
				int iy = indices[1];
				energies[iy][ix] = energy;

				if (energy > 5.0) {
					methylAngles1[iy][ix] = 0.0;
					methylAngles2[iy][ix] = 0.0;
				} else {
					methylAngles1[iy][ix] = Protractor.getDistDegrees(methyl[0].measure(), 60);
					methylAngles2[iy][ix] = Protractor.getDistDegrees(methyl[1].measure(), 60);
				}

				progress.incrementProgress();
			}

			// sample the axes
			double[] xaxis = new double[360];
			double[] yaxis = new double[360];
			for (int i=0; i<360; i++) {
				xaxis[i] = i;
				yaxis[i] = i;
			}

			VisIt.writeGrid2D(xaxis, yaxis, energies, new File(type + ".dihedralEnergy.fixedTetraCACBCG.minMethyl.vtk"));
			VisIt.writeGrid2D(xaxis, yaxis, methylAngles1, new File(type + ".methylAngles1.vtk"));
			VisIt.writeGrid2D(xaxis, yaxis, methylAngles2, new File(type + ".methylAngles2.vtk"));
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
}
