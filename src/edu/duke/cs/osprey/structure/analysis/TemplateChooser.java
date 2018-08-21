package edu.duke.cs.osprey.structure.analysis;


import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.analysis.MeasurementLibrary.*;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Protractor;

import java.io.File;
import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


/** chooses the most ideal template for each rotamer from a folder full of PDB files */
public class TemplateChooser {
	
	public static void main(String[] args) {

		Map<String,List<MeasuredRes>> measurementsByType = new HashMap<>();

		// scan all the PDB files to collect the measurements
		PDBScanner scanner = new PDBScanner(new File(args[0]));
		scanResidues(scanner, (filename, res, type) -> {

			// gather all the measurements
			double[] measurements = mlib.measure(res);
			if (measurements == null) {
				return;
			}

			// make sure we don't have any bad measurements
			for (double measurement : measurements) {
				if (!Double.isFinite(measurement)) {
					List<Measurement> libMeasurements = mlib.get(type);
					StringBuilder buf = new StringBuilder();
					buf.append(String.format("invalid measurements for %s %s %s, this is a bug!", filename, res.getPDBResNumber(), type));
					for (int m=0; m<measurements.length; m++) {
						buf.append(String.format("\n\t%28s: %8.3f", libMeasurements.get(m).name, measurements[m]));
					}
					throw new Error(buf.toString());
				}
			}

			// group all the measurements by type,voxel
			measurementsByType
				.computeIfAbsent(type, t -> new ArrayList<>())
				.add(new MeasuredRes(filename, res.getPDBResNumber(), measurements));
		});

		// analyze all the measurements to find the modes
		log("calculating modal values...");
		Map<String,double[]> modesByType = new HashMap<>();
		for (Map.Entry<String,List<MeasuredRes>> entry : measurementsByType.entrySet()) {

			String type = entry.getKey();
			List<MeasuredRes> measuredResidues = entry.getValue();
			log("%s: %d samples", type, measuredResidues.size());

			List<Measurement> measurements = mlib.get(type);

			double[] modes = new double[measurements.size()];
			Arrays.fill(modes, Double.NaN);

			// analyze the distributions of each measurement
			for (int m=0; m<measurements.size(); m++) {

				Measurement measurement = measurements.get(m);
				logf("%28s  ", measurement.name);

				switch (measurement.space) {

					case R: {
						ClusterR1 cluster = new ClusterR1();
						for (MeasuredRes mres : measuredResidues) {
							cluster.add(mres.measurements[m]);
						}
						ClusterR1.Stats stats = cluster.new Stats(90, 10);
						MathTools.DoubleBounds bounds90 = stats.getInterval(stats.mode, 90);
						logf(" mean=%6.3f  mode=%6.3f  bounds=[%6.3f,%6.3f]:%6.3f  90%%bounds=[%6.3f,%6.3f]:%6.3f",
							stats.mean, stats.mode,
							stats.bounds.lower, stats.bounds.upper, stats.bounds.size(),
							bounds90.lower, bounds90.upper, bounds90.size()
						);
						modes[m] = stats.mode;
					} break;

					case S: {
						ClusterS1 cluster = new ClusterS1();
						for (MeasuredRes mres : measuredResidues) {
							cluster.add(mres.measurements[m]);
						}
						ClusterS1.Stats stats = cluster.new Stats(90, 40);
						SmallAngleVoxel.Interval bounds90 = stats.getInterval(stats.mode, 90);
						logf(" mean=%6.1f  mode=%6.1f  bounds=[%6.1f,%6.1f]:%6.1f  90%%bounds=[%6.1f,%6.1f]:%6.1f",
							stats.mean, stats.mode,
							stats.bounds.min(), stats.bounds.max(), stats.bounds.size(),
							bounds90.min(), bounds90.max(), bounds90.size()
						);
						modes[m] = stats.mode;
					} break;
				}

				log("");
			}

			modesByType.put(type, modes);
		}

		// define a distance function for residues
		DistanceFunction measurementDist = (measurements, a, b) -> {

			// let's do manhattan distance
			// but try to normalize for the size of each space
			double dist = 0;
			for (int m=0; m<measurements.size(); m++) {
				switch (measurements.get(m).space) {

					case R: {
						dist += Math.abs(a[m] - b[m]);
					} break;

					case S: {
						dist += Protractor.getDistDegrees(a[m], b[m])/10;
					} break;
				}
			}
			return dist;
		};

		// find the residue that most matches the modes
		log("finding ideal residues...");
		Map<String,MeasuredRes> idealResidues = new HashMap<>();
		scanResidues(scanner, (filename, res, type) -> {

			double[] measurements = mlib.measure(res);
			if (measurements == null) {
				return;
			}

			// do we have an ideal-est residue yet?
			MeasuredRes oldMres = idealResidues.get(type);
			if (oldMres == null) {

				// nope, assume this one is most ideal for now
				idealResidues.put(type, new MeasuredRes(filename, res.getPDBResNumber(), measurements));

			} else {

				// yup, is the new one more ideal?
				List<Measurement> libMeasurements = mlib.get(type);
				double[] modes = modesByType.get(type);
				double oldDist = measurementDist.calc(libMeasurements, modes, oldMres.measurements);
				double newDist = measurementDist.calc(libMeasurements, modes, measurements);
				if (newDist < oldDist) {

					// yup, make the new residue the current ideal residue
					idealResidues.put(type, new MeasuredRes(filename, res.getPDBResNumber(), measurements));
				}
			}
		});

		// show the ideal residues
		for (Map.Entry<String,MeasuredRes> entry : idealResidues.entrySet()) {
			String type = entry.getKey();
			MeasuredRes mres = entry.getValue();

			List<Measurement> measurements = mlib.get(type);
			double[] modes = modesByType.get(type);
			double dist = measurementDist.calc(measurements, modes, mres.measurements);

			log("%s: %s %s   d=%8.4f", type, mres.filename, mres.resNum, dist);
			for (int m=0; m<measurements.size(); m++) {

				Measurement measurement = measurements.get(m);
				logf("%28s  ", measurement.name);

				switch (measurement.space) {

					case R: {
						logf("%6.3f  mode=%6.3f  dist=%6.3f",
							mres.measurements[m],
							modes[m],
							Math.abs(modes[m] - mres.measurements[m])
						);
					} break;

					case S: {
						logf("%6.1f  mode=%6.1f  dist=%6.1f (%6.1f)",
							mres.measurements[m],
							modes[m],
							Protractor.getDistDegrees(modes[m], mres.measurements[m]),
							Protractor.getDistDegrees(modes[m], mres.measurements[m])/10
						);
					} break;
				}

				log("");
			}
		}
	}

	private static void scanResidues(PDBScanner pdbScanner, ResidueScanner resScanner) {

		// scan all the PDB files to collect the measurements
		pdbScanner.scan((file, mol) -> {

			// analyze each residue
			for (Residue res : mol.residues) {
				resScanner.onResidue(file.getName(), res, res.getType().toUpperCase());
			}
		});
	}

	/* TODO: try this if osprey ever supports per-voxel template coords
	private static void scanResidues(PDBScanner pdbScanner, ResidueVoxelScanner resScanner) {

		// scan all the PDB files to collect the measurements
		pdbScanner.scan((file, mol) -> {

			// analyze each residue
			for (Residue res : mol.residues) {
				String type = res.getType().toUpperCase();

				// what voxel is this?
				double[] chiAngles = chiLib.measure(res);
				if (chiAngles == null) {
					continue;
				}
				Voxel voxel = voxelLib.find(type, chiAngles);
				if (voxel == null) {
					continue;
				}

				resScanner.onResidue(file.getName(), res, type, voxel);
			}
		});
	}


	// define the amino acid libraries
	private static final MeasurementLibrary chiLib = new MeasurementLibrary();
	private static final VoxelLibrary voxelLib = new VoxelLibrary();
	static {

		// leucine
		chiLib.add("LEU", new DihedralAngle("chi1", "N", "CA", "CB", "CG"));
		chiLib.add("LEU", new DihedralAngle("chi2", "CA", "CB", "CG", "CD1"));
		voxelLib.add("LEU", "pp", 0.6,   -15, 62, 14,    -19, 81, 13);
		voxelLib.add("LEU", "pt", 0.2,   -10, 75, 7,     -10, 169, 10);
		voxelLib.add("LEU", "mm", 0.5,   -17, -82, 18,   -13, -64, 20);
		voxelLib.add("LEU", "tt", 2.0,   -19, -173, 20,  -21, 155, 21);
		voxelLib.add("LEU", "tp", 28.9,  -32, -172, 43,  -26, 62, 30);
		voxelLib.add("LEU", "mt", 61.7,  -42, -70, 33,   -36, 175, 32);
		voxelLib.add("LEU", "tm", 0.2,   -10, -174, 9,   -9, -76, 9);
		voxelLib.add("LEU", "mp", 3.5,   -28, -82, 32,   -52, 60, 43);
		voxelLib.add("LEU", "xx", 0.2,   -8, -146, 13,   -8, -152, 8); // unnamed by rlab
	}
	*/

	private static final MeasurementLibrary mlib = new MeasurementLibrary();
	static {

		// backbone measurements common to all amino acids
		for (String type : Arrays.asList("LEU")) {
			mlib.add(type, new BondAngle("N", "CA", "C"));
			mlib.add(type, new BondLength("N", "CA"));
			mlib.add(type, new BondLength("CA", "C"));
		}

		// backbone measurements common to all amino acids with a beta carbon
		for (String type : Arrays.asList("LEU")) {
			mlib.add(type, new TetrahedralInPlaneAngle("N", "CA", "C", "CB"));
			mlib.add(type, new TetrahedralOutOfPlaneAngle("N", "CA", "C", "CB"));
			mlib.add(type, new BondLength("CA", "CB"));
		}

		// leucine
		mlib.add("LEU", new BondAngle("CA", "CB", "CG"));
		mlib.add("LEU", new BondAngle("CB", "CG", "CD1"));
		mlib.add("LEU", new BondAngle("CB", "CG", "CD2"));
		mlib.add("LEU", new BondLength("CB", "CG"));
		mlib.add("LEU", new BondLength("CG", "CD1"));
		mlib.add("LEU", new BondLength("CG", "CD2"));
		addTetrahedralH2Measurements("LEU", "CA", "CB", "CG");
		addH3Measurements("LEU", "CB", "CG", "CD1");
		addH3Measurements("LEU", "CB", "CG", "CD2");
	}

	private static void addTetrahedralH2Measurements(String type, String b, String c, String d) {
		String cpos = c.substring(1);
		String[] h = new String[] {
			"H" + cpos + "2",
			"H" + cpos + "3"
		};
		for (int i=0; i<2; i++) {
			mlib.add(type, new TetrahedralInPlaneAngles2(b, c, d, h[0], h[1], i));
			mlib.add(type, new TetrahedralOutOfPlaneAngles2(b, c, d, h[0], h[1], i));
			mlib.add(type, new BondLength(c, h[i]));
		}
	}

	private static void addH3Measurements(String type, String a, String b, String c) {
		String cpos = c.substring(1);
		String[] h = new String[] {
			"H" + cpos + "1",
			"H" + cpos + "2",
			"H" + cpos + "3"
		};
		for (int i=0; i<3; i++) {
			mlib.add(type, new DihedralAngles3(a, b, c, h[0], h[1], h[2], i));
			mlib.add(type, new BondAngle(b, c, h[i]));
			mlib.add(type, new BondLength(c, h[i]));
		}
	}


	private static interface DistanceFunction {
		double calc(List<Measurement> measurements, double[] a, double[] b);
	}

	private static interface ResidueScanner {
		void onResidue(String filename, Residue res, String type);
	}

	private static interface ResidueVoxelScanner {
		void onResidue(String filename, Residue res, String type, Voxel voxel);
	}

	private static class MeasuredRes {

		public final String filename;
		public final String resNum;
		public final double[] measurements;

		public MeasuredRes(String filename, String resNum, double[] measurements) {
			this.filename = filename;
			this.resNum = resNum;
			this.measurements = measurements;
		}
	}

	private static class Voxel {

		public final String name;
		public final double percent;
		public final SmallAngleVoxel voxel;

		public Voxel(String name, double percent, SmallAngleVoxel voxel) {
			this.name = name;
			this.percent = percent;
			this.voxel = voxel;
		}

		public Voxel(String name, double percent, double ... bounds) {
			this(name, percent, SmallAngleVoxel.makeFromBounds(bounds));
		}

		public boolean contains(double[] dihedrals) {
			if (voxel == null) {
				return true;
			}
			return voxel.contains(dihedrals);
		}
	}

	private static class VoxelLibrary {

		private final Map<String,List<Voxel>> voxels = new HashMap<>();

		public List<Voxel> get(String type) {
			return voxels.computeIfAbsent(type.toUpperCase(), t -> new ArrayList<>());
		}

		public void add(String type, String name, double percent, double ... bounds) {
			add(type, new Voxel(name, percent, bounds));
		}

		public void add(String type, Voxel voxel) {
			get(type).add(voxel);
		}

		public Voxel find(String type, double[] p) {
			for (Voxel voxel : get(type)) {
				if (voxel.contains(p)) {
					return voxel;
				}
			}
			return null;
		}
	}
}
