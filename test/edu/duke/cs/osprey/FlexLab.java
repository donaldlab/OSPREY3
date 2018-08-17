package edu.duke.cs.osprey;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.*;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.analysis.*;
import edu.duke.cs.osprey.tools.*;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class FlexLab {

	public static void main(String[] args)
	throws Exception {
		//checkRotamerClashes();
		top8000();
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

		ForcefieldParams ffparams = new ForcefieldParams();
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build()) {

			/* TEMP
			// for each rotamer
			for (SimpleConfSpace.ResidueConf rc : confSpace.positions.get(0).resConfs) {

				RCTuple tuple = new RCTuple(0, rc.index);
				ParametricMolecule pmol = confSpace.makeMolecule(tuple);

				//EnergyCalculator.EnergiedParametricMolecule epmol = ecalc.calcEnergy(confSpace.makeMolecule(tuple), inters);

				// look at probe contacts
				log("%s %d - %3d", rc.template.name, rc.rotamerIndex, rc.index);
				List<Probe.AtomPair.Interaction> interactions = probe.getInteractions(pmol.mol.residues, inters, connectivity).stream()
					.filter(interaction -> interaction.contact.isContact)
					.sorted(Comparator.comparing(interaction -> interaction.getViolation(0.0)))
					.collect(Collectors.toList());

				for (Probe.AtomPair.Interaction interaction : interactions) {
					if (interaction.contact.isClash) {
						double vdw = calcVDW(interaction.atomPair.a, interaction.atomPair.b, ffparams);
						log("\t%s   %s     vdw=%.4f", interaction.atomPair, interaction, vdw);
					}
				}
			}
			*/

			// look at individual structures
			RCTuple tuple = new RCTuple(0, 99);
			ParametricMolecule pmol = confSpace.makeMolecule(tuple);
			PDBIO.writeFile(pmol.mol, new File("rotamer.pdb"));
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

	public static void top8000()
	throws Exception {

		// define all the dihedrals
		AngleLibrary angleLib = new AngleLibrary();
		angleLib.add("LEU", new AngleLibrary.Dihedral("chi1", "N", "CA", "CB", "CG"));
		angleLib.add("LEU", new AngleLibrary.Dihedral("chi2", "CA", "CB", "CG", "CD1"));

		// TODO: fix template matching system
		//ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();

		final Map<DihedralKey,double[]> dihedralsByKey;
		if (false) {

			dihedralsByKey = new HashMap<>();

			// these all have duplicated residues for some reason
			new PDBScanner(new File("/home/jeff/dlab/top8000"),
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
			).scan((file, mol) -> {

				// analyze each residue
				for (Residue res : mol.residues) {

					// TEMP: just leucines for now
					String type = res.getType();
					if (!type.equals("LEU")) {
						continue;
					}

					// TODO: filter residues by probe clashes?

					/* TODO: fix template assignment!!
					// DANGER WILL ROBINSON!!! assigning templates causes dihedral angle measurement to be incorrect!!
					// assign a template
					res.assignTemplate(templateLib);
					if (res.template == null) {
						log("WARN: no template assigned for %s %s", file.getName(), res.fullName);
						continue;
					}
					*/

					DihedralKey key = new DihedralKey(type, file.getName(), res.getPDBResNumber());

					// measure the dihedrals
					double[] dihedrals = angleLib.measure(res);
					if (dihedrals == null) {
						continue;
					}

					dihedralsByKey.put(key, dihedrals);
				}
			});

			// write all the dihedrals to a file
			ObjectIO.write(dihedralsByKey, "dihedrals.leu.dat");

		} else {

			log("loading dihedrals");
			//noinspection unchecked
			dihedralsByKey = ObjectIO.read("dihedrals.leu.dat", HashMap.class);
		}

		log("LEU dihedrals: %d", dihedralsByKey.size());
		int numAngles = 2;

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
			hist.filterDensityWindow(2, 20);
			int keptCount = hist.count();
			log("after filtering:  %d  %.1f%%", keptCount, 100f*keptCount/count);
		}

		// convert the histogram back to dihedrals
		List<double[]> keptDihedrals = new ArrayList<>();
		for (long key : hist.buckets.keySet()) {
			keptDihedrals.add(hist.makeDihedrals(key));
		}

		// cluster the points
		List<List<double[]>> clusters = AngleClustering.cluster(keptDihedrals, numAngles, 50.0);

		// calc bounding voxels for the clusters
		List<SmallAngleVoxel> voxels = clusters.stream()
			.map(cluster -> AngleClustering.calcVoxel(cluster))
			.collect(Collectors.toList());

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

		// analyze the Lovell rotamers
		List<SmallAngleVoxel> lovellVoxels = Arrays.asList(
			new SmallAngleVoxel(new double[] { 62, 80 }),
			new SmallAngleVoxel(new double[] { -177, 65 }),
			new SmallAngleVoxel(new double[] { -172, 145 }),
			new SmallAngleVoxel(new double[] { -85, 65 }),
			new SmallAngleVoxel(new double[] { -65, 175 })
		);
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

		// vis it
		VisIt.writeDihedrals(dihedralsByKey.values(), new File("leu.dihedrals.vtk"));
		VisIt.writeDihedrals(keptDihedrals, new File("leu.keptDihedrals.vtk"));
		VisIt.writeDihedrals(clusters.get(5), new File("leu.cluster.vtk"));
		VisIt.writeVoxels(voxels, new File("leu.voxels.vtk"));
		VisIt.writeVoxels(fixedVoxels, new File("leu.fixedVoxels.vtk"));
		VisIt.writeVoxels(lovellVoxels, new File("leu.lovellVoxels.vtk"));
	}

	static class DihedralKey implements Serializable {

		private static final long serialVersionUID = 7858382121141521600L;

		public final String type;
		public final String filename;
		public final String resNum;

		public DihedralKey(String type, String filename, String resNum) {
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
			return other instanceof DihedralKey && equals((DihedralKey)other);
		}

		public boolean equals(DihedralKey other) {
			return this.type.equals(other.type)
				&& this.filename.equals(other.filename)
				&& this.resNum.equals(other.resNum);
		}

		@Override
		public String toString() {
			return String.format("%s:%s:%s", type, filename, resNum);
		}
	}
}
