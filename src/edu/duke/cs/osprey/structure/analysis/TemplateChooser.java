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

package edu.duke.cs.osprey.structure.analysis;


import edu.duke.cs.osprey.confspace.VoxelShape;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Probe;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.analysis.MeasurementLibrary.*;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Protractor;

import java.io.File;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


/** chooses the most ideal template for each residue type from a folder full of PDB files */
public class TemplateChooser {
	
	public static void main(String[] args) {

		File dir = new File(args[0]);

		// what residues types do we expect to see?
		List<String> expectedTypes = Arrays.asList(
			"ARG", "HIP", "HID", "HIE", "LYS",
			"ASP", "GLU",
			"CYS", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP",
			"SER", "THR", "ASN", "GLN"
		);

		// measure all the residues in the PDB dir
		Map<String,List<MeasuredRes>> measurementsByType = measureResidues(dir, expectedTypes);

		// analyze all the measurements to find the modes
		log("calculating modal values...");
		Map<String,double[]> modesByType = calcModes(measurementsByType, expectedTypes);

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
			return 100.0*dist/measurements.size();
		};

		// define a distance function for probe violations over all the rotamers
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();
		Probe probe = new Probe();
		probe.matchTemplates(templateLib);
		logf("building atom connectivity...");
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(templateLib)
			.set15HasNonBonded(false)
			.build();
		log(" done!");
		BiFunction<Residue,String,Double> rmsViolations = (res, resType) -> {

			// try to assign a template to the residue
			res.assignTemplateSimple(templateLib, resType);
			if (res.template == null) {
				return Double.POSITIVE_INFINITY;
			}

			List<Double> violations = new ArrayList<>();
			Runnable collectViolations = () -> {
				probe.getInteractions(res, connectivity).stream()
					.mapToDouble(inter -> inter.getViolation(0.0))
					.filter(v -> v > 0.0)
					.forEach(v -> violations.add(v));
			};

			// for each rotamer
			if (res.template.getNumRotamers() <= 0) {

				// no rotamer, eg gly, pro, ala
				collectViolations.run();

			} else {
				VoxelShape voxel = new VoxelShape.Rect();
				List<DegreeOfFreedom> dofs = voxel.makeDihedralDOFs(res);
				for (int i=0; i<res.template.getNumRotamers(); i++) {

					// center this conformation on the rotamer voxel
					ObjectiveFunction.DofBounds dofBounds = voxel.makeDihedralBounds(res.template, i);
					for (int d=0; d<dofs.size(); d++) {
						dofs.get(d).apply(dofBounds.getCenter(d));
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
		};

		// find the residue that most matches the modes
		log("finding ideal residues...");
		Map<String,MeasuredRes> idealResidues = new HashMap<>();
		scanResidues(dir, (filename, res, type) -> {

			// skip terminal residues
			// osprey can't design with terminal residues yet
			if (res.getAtomByName("H3") != null || res.getAtomByName("OXT") != null) {
				return;
			}

			double[] measurements = mlib.measure(res, type);
			if (measurements == null) {
				return;
			}

			double rmsViolation = rmsViolations.apply(res, type);

			// do we have an ideal-est residue yet?
			MeasuredRes oldMres = idealResidues.get(type);
			if (oldMres == null) {

				// nope, assume this one is most ideal for now
				idealResidues.put(type, new MeasuredRes(filename, res.getPDBResNumber(), measurements, rmsViolation));

			} else {

				final double idealWeight = 1.0;
				final double probeWeight = 12.0;

				// yup, is the new one more ideal?
				List<Measurement> libMeasurements = mlib.get(type);
				double[] modes = modesByType.get(type);
				double oldDist = measurementDist.calc(libMeasurements, modes, oldMres.measurements)*idealWeight + oldMres.rmsViolation*probeWeight;
				double newDist = measurementDist.calc(libMeasurements, modes, measurements)*idealWeight + rmsViolation*probeWeight;
				if (newDist < oldDist) {

					// yup, make the new residue the current ideal residue
					idealResidues.put(type, new MeasuredRes(filename, res.getPDBResNumber(), measurements, rmsViolation));
				}
			}
		});

		// show the ideal residues and the current templates
		for (String type : expectedTypes) {
			MeasuredRes mres = idealResidues.get(type);
			int count = measurementsByType.get(type).size();

			List<Measurement> measurements = mlib.get(type);
			double[] modes = modesByType.get(type);
			double dist = measurementDist.calc(measurements, modes, mres.measurements);

			ResidueTemplate template = templateLib.getTemplate(type, true);
			double[] templateMeasurements = mlib.measure(template.templateRes, type);
			double templateDist = measurementDist.calc(measurements, modes, templateMeasurements);
			double templateRmsViolation = rmsViolations.apply(template.templateRes, type);

			log("%s: %s %s   n=%7d   thisd=%8.4f   currentd=%8.4f   thisRmsV=%5.3f   currentRmsV=%5.3f",
				type, mres.filename, mres.resNum, count, dist, templateDist, mres.rmsViolation, templateRmsViolation
			);
			for (int m=0; m<measurements.size(); m++) {

				Measurement measurement = measurements.get(m);
				logf("%30s  ", measurement.name);

				switch (measurement.space) {

					case R: {
						logf("mode=%6.3f  this=%6.3f d=%6.3f  current=%6.3f d=%6.3f  delta=%6.3f",
							modes[m],
							mres.measurements[m],
							Math.abs(modes[m] - mres.measurements[m]),
							templateMeasurements[m],
							Math.abs(modes[m] - templateMeasurements[m]),
							Math.abs(mres.measurements[m] - templateMeasurements[m])
						);
					} break;

					case S: {
						logf("mode=%6.1f  this=%6.1f d=%6.1f  current=%6.1f d=%6.1f  delta=%6.1f",
							modes[m],
							mres.measurements[m],
							Protractor.getDistDegrees(modes[m], mres.measurements[m]),
							templateMeasurements[m],
							Protractor.getDistDegrees(modes[m], templateMeasurements[m]),
							Protractor.getDistDegrees(mres.measurements[m], templateMeasurements[m])
						);
					} break;
				}

				log("");
			}
		}

		/* write out the ideal coords in osprey format, eg:
			GLY 7
			N  2.027f  1.358f  0.000f
			H  1.697f  1.839f  0.824f
			CA  1.522f  0.000f  0.000f
			HA2  1.886f  -0.523f  -0.885f
			HA3  1.874f  -0.506f  0.899f
			C  0.000f  0.000f  0.000f
			O  -0.624f  1.058f  0.000f
			ENDRES
		*/
		log("template coords:");
		for (String type : expectedTypes) {
			MeasuredRes mres = idealResidues.get(type);

			new PDBScanner(dir).scan(mres.filename, (file, mol) -> {

				Residue res = mol.getResByPDBResNumber(mres.resNum);

				log("%s %d", getResType(res), res.atoms.size());
				for (Atom atom : res.atoms) {
					int i3 = atom.indexInRes*3;
					log("%s  %.3ff  %.3ff  %.3ff",
						atom.name,
						res.coords[i3    ],
						res.coords[i3 + 1],
						res.coords[i3 + 2]
					);
				}
				log("ENDRES");
			});
		}
	}

	public static Map<String,List<MeasuredRes>> measureResidues(File dir, Collection<String> types) {
		return measureResidues(dir, types, res -> true);
	}

	public static Map<String,List<MeasuredRes>> measureResidues(File dir, Collection<String> types, Predicate<Residue> filter) {

		Map<String,List<MeasuredRes>> measurementsByType = new HashMap<>();

		// scan all the PDB files to collect the measurements
		scanResidues(dir, (filename, res, type) -> {

			// only check types we care about
			if (!types.contains(type)) {
				return;
			}

			// gather all the measurements
			double[] measurements = mlib.measure(res, type);
			if (measurements == null) {
				return;
			}

			// apply the filter
			if (!filter.test(res)) {
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

			// group all the measurements by type
			measurementsByType
				.computeIfAbsent(type, t -> new ArrayList<>())
				.add(new MeasuredRes(filename, res.getPDBResNumber(), measurements, Double.NaN));
		});

		return measurementsByType;
	}

	public static Map<String,double[]> calcModes(Map<String,List<MeasuredRes>> measurementsByType, List<String> types) {

		boolean foundAllExpected = true;
		Map<String,double[]> modesByType = new HashMap<>();
		for (String type : types) {

			List<MeasuredRes> measuredResidues = measurementsByType.get(type);
			if (measuredResidues == null) {
				log("%s: no samples! are the measurements wrong?", type);
				foundAllExpected = false;
				continue;
			}

			log("%s: %d samples", type, measuredResidues.size());

			List<Measurement> measurements = mlib.get(type);

			double[] modes = new double[measurements.size()];
			Arrays.fill(modes, Double.NaN);

			// analyze the distributions of each measurement
			for (int m=0; m<measurements.size(); m++) {

				Measurement measurement = measurements.get(m);
				logf("%30s  ", measurement.name);

				switch (measurement.space) {

					case R: {
						ClusterR1 cluster = new ClusterR1();
						for (MeasuredRes mres : measuredResidues) {
							cluster.add(mres.measurements[m]);
						}
						ClusterR1.Stats stats = cluster.new Stats(90, 10);
						MathTools.DoubleBounds bounds90 = stats.getInterval(stats.mode, 90);
						logf(" mean=%6.3f  mode=%6.3f  90%%bounds=[%6.3f,%6.3f]:%6.3f  100%%bounds=[%6.3f,%6.3f]:%6.3f",
							stats.mean, stats.mode,
							bounds90.lower, bounds90.upper, bounds90.size(),
							stats.bounds.lower, stats.bounds.upper, stats.bounds.size()
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
						logf(" mean=%6.1f  mode=%6.1f  90%%bounds=[%6.1f,%6.1f]:%6.1f  100%%bounds=[%6.1f,%6.1f]:%6.1f",
							stats.mean, stats.mode,
							bounds90.min(), bounds90.max(), bounds90.size(),
							stats.bounds.min(), stats.bounds.max(), stats.bounds.size()
						);
						modes[m] = stats.mode;
					} break;
				}

				log("");
			}

			modesByType.put(type, modes);
		}

		if (!foundAllExpected) {

			// something's missing... stop here and find out what
			throw new NoSuchElementException("missing data for some residue types, aborting: "
				+ types.stream()
					.filter(type -> !measurementsByType.containsKey(type))
					.collect(Collectors.toList())
			);
		}

		return modesByType;
	}

	private static void scanResidues(File dir, ResidueScanner resScanner) {

		// scan all the PDB files to collect the measurements
		new PDBScanner(dir).scan((file, mol) -> {

			// analyze each residue
			for (Residue res : mol.residues) {

				// classify the histidines based on protonation state
				String type = getResType(res);
				if (type == null) {
					System.err.println("error reading residue at " + file.getName() + " " + res.fullName + ", skipping it: unprotonated HIS");
					continue;
				}

				resScanner.onResidue(file.getName(), res, type);
			}
		});
	}

	private static String getResType(Residue res) {

		// classify the histidines based on protonation state
		String type = res.getType().toUpperCase();
		if (type.equals("HIS")) {

			boolean hasD = res.getAtomByName("HD1") != null;
			boolean hasE = res.getAtomByName("HE2") != null;

			if (hasD && hasE) {
				return "HIP";
			} else if (hasD) {
				return "HID";
			} else if (hasE) {
				return "HIE";
			} else {
				return null;
			}
		}

		// everything else can just pass through
		return type;
	}

	/* TODO: try this if osprey ever supports per-voxel template coords
	private static void scanResidues(PDBScanner pdbScanner, ResidueVoxelScanner resScanner) {

		// scan all the PDB files to collect the measurements
		pdbScanner.scan((file, mol) -> {

			// analyze each residue
			for (Residue res : mol.residues) {
				String type = res.getType().toUpperCase();

				// what voxel is this?
				double[] chiAngles = chiLib.measure(res, type);
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
		for (String type : Arrays.asList("ARG", "HIP", "HID", "HIE", "LYS", "ASP", "GLU", "CYS", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER", "THR", "ASN", "GLN")) {
			mlib.add(type, new BondAngle("N", "CA", "C"));
			mlib.add(type, new BondLength("N", "CA"));
			mlib.add(type, new BondLength("CA", "C"));
		}

		// measurements common to all amino acids with a beta carbon
		for (String type : Arrays.asList("ARG", "HIP", "HID", "HIE", "LYS", "ASP", "GLU", "CYS", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "SER", "THR", "ASN", "GLN")) {
			addTetrahedralMeasurements(type, "N", "CA", "C", "CB");
			addTetrahedralMeasurements(type, "N", "CA", "C", "HA");
		}

		// check means vs mode: big differences could be problems?

		// arginine
		addTetrahedralH2Measurements("ARG", "CA", "CB", "CG");
		addLinkMeasurements(         "ARG", "CA", "CB", "CG");
		addTetrahedralH2Measurements("ARG", "CB", "CG", "CD");
		addLinkMeasurements(         "ARG", "CB", "CG", "CD");
		addTetrahedralH2Measurements("ARG", "CG", "CD", "NE");
		addLinkMeasurements(         "ARG", "CG", "CD", "NE");
		addTrigonalMeasurements(     "ARG", "CD", "NE", "CZ", "HE");
		addLinkMeasurements(         "ARG", "CD", "NE", "CZ");
		addBranchMeasurements(       "ARG", "NE", "CZ", "NH1", "NH2");
		addH2Measurements(           "ARG", "NE", "CZ", "NH1");
		addH2Measurements(           "ARG", "NE", "CZ", "NH2");

		// histidine w/ HD1 and HE2
		addTetrahedralH2Measurements("HIP", "CA", "CB", "CG");
		addLinkMeasurements(         "HIP", "CA", "CB", "CG");
		addBranchMeasurements(       "HIP", "CB", "CG", "ND1", "CD2");
		addTrigonalMeasurements(     "HIP", "CG", "ND1", "CE1", "HD1"); // HD1
		addTrigonalMeasurements(     "HIP", "CG", "CD2", "NE2", "HD2");
		addRingLinkMeasurements(     "HIP", "CB", "CG", "ND1", "CE1");
		addTrigonalMeasurements(     "HIP", "ND1", "CE1", "NE2", "HE1");
		addRingLinkMeasurements(     "HIP", "CB", "CG", "CD2", "NE2");
		addTrigonalMeasurements(     "HIP", "CD2", "NE2", "CE1", "HE2"); // HE2
		addRingLinkMeasurements(     "HIP", "CG", "ND1", "CE1", "NE2");

		// histidine w/ only HD1
		addTetrahedralH2Measurements("HID", "CA", "CB", "CG");
		addLinkMeasurements(         "HID", "CA", "CB", "CG");
		addBranchMeasurements(       "HID", "CB", "CG", "ND1", "CD2");
		addTrigonalMeasurements(     "HID", "CG", "ND1", "CE1", "HD1"); // HD1
		addTrigonalMeasurements(     "HID", "CG", "CD2", "NE2", "HD2");
		addRingLinkMeasurements(     "HID", "CB", "CG", "ND1", "CE1");
		addTrigonalMeasurements(     "HIP", "ND1", "CE1", "NE2", "HE1");
		addRingLinkMeasurements(     "HID", "CB", "CG", "CD2", "NE2");
		addRingLinkMeasurements(     "HID", "CG", "ND1", "CE1", "NE2");

		// histidine w/ only HE2
		addTetrahedralH2Measurements("HIE", "CA", "CB", "CG");
		addLinkMeasurements(         "HIE", "CA", "CB", "CG");
		addBranchMeasurements(       "HIE", "CB", "CG", "ND1", "CD2");
		addTrigonalMeasurements(     "HIE", "CG", "CD2", "NE2", "HD2");
		addRingLinkMeasurements(     "HIE", "CB", "CG", "ND1", "CE1");
		addTrigonalMeasurements(     "HIE", "ND1", "CE1", "NE2", "HE1");
		addRingLinkMeasurements(     "HIE", "CB", "CG", "CD2", "NE2");
		addTrigonalMeasurements(     "HIE", "CD2", "NE2", "CE1", "HE2"); // HE2
		addRingLinkMeasurements(     "HIE", "CG", "ND1", "CE1", "NE2");

		// lysine
		addTetrahedralH2Measurements("LYS", "CA", "CB", "CG");
		addLinkMeasurements(         "LYS", "CA", "CB", "CG");
		addTetrahedralH2Measurements("LYS", "CB", "CG", "CD");
		addLinkMeasurements(         "LYS", "CB", "CG", "CD");
		addTetrahedralH2Measurements("LYS", "CG", "CD", "CE");
		addLinkMeasurements(         "LYS", "CG", "CD", "CE");
		addTetrahedralH2Measurements("LYS", "CD", "CE", "NZ");
		addLinkMeasurements(         "LYS", "CD", "CE", "NZ");
		addH3Measurements(           "LYS", "CD", "CE", "NZ");

		// apartic acid
		addTetrahedralH2Measurements("ASP", "CA", "CB", "CG");
		addLinkMeasurements(         "ASP", "CA", "CB", "CG");
		addBranchMeasurements(       "ASP", "CB", "CG", "OD1", "OD2");

		// glutamic acid
		addTetrahedralH2Measurements("GLU", "CA", "CB", "CG");
		addLinkMeasurements(         "GLU", "CA", "CB", "CG");
		addTetrahedralH2Measurements("GLU", "CB", "CG", "CD");
		addLinkMeasurements(         "GLU", "CB", "CG", "CD");
		addBranchMeasurements(       "GLU", "CG", "CD", "OE1", "OE2");

		// cysteine
		addTetrahedralH2Measurements("CYS", "CA", "CB", "SG");
		addLinkMeasurements(         "CYS", "CA", "CB", "SG");
		addLinkMeasurements(         "CYS", "CB", "SG", "HG");

		// glycine
		addTetrahedralH2Measurements("GLY", "N", "CA", "C");

		// proline
		addTetrahedralH2Measurements("PRO", "CA", "CB", "CG");
		addLinkMeasurements(         "PRO", "CA", "CB", "CG");
		addTetrahedralH2Measurements("PRO", "CB", "CG", "CD");
		addLinkMeasurements(         "PRO", "CB", "CG", "CD");
		addTetrahedralH2Measurements("PRO", "CG", "CD", "N");
		addLinkMeasurements(         "PRO", "CG", "CD", "N");

		// alanine
		addH3Measurements(           "ALA", "N", "CA", "CB");

		// valine
		addTetrahedralMeasurements(  "VAL", "CA", "CB", "CG1", "HB");
		addBranchMeasurements(       "VAL", "CA", "CB", "CG1", "CG2");
		addH3Measurements(           "VAL", "CA", "CB", "CG1");
		addH3Measurements(           "VAL", "CA", "CB", "CG2");

		// isoleucine
		addTetrahedralMeasurements(  "ILE", "CA", "CB", "CG1", "HB");
		addBranchMeasurements(       "ILE", "CA", "CB", "CG1", "CG2");
		addH3Measurements(           "ILE", "CA", "CB", "CG2");
		addTetrahedralH2Measurements("ILE", "CB", "CG1", "CD1");
		addLinkMeasurements(         "ILE", "CB", "CG1", "CD1");
		addH3Measurements(           "ILE", "CB", "CG1", "CD1");

		// leucine
		addTetrahedralH2Measurements("LEU", "CA", "CB", "CG");
		addLinkMeasurements(         "LEU", "CA", "CB", "CG");
		addTetrahedralMeasurements(  "LEU", "CB", "CG", "CD1", "HG");
		addBranchMeasurements(       "LEU", "CB", "CG", "CD1", "CD2");
		addH3Measurements(           "LEU", "CB", "CG", "CD1");
		addH3Measurements(           "LEU", "CB", "CG", "CD2");

		// methionine
		addTetrahedralH2Measurements("MET", "CA", "CB", "CG");
		addLinkMeasurements(         "MET", "CA", "CB", "CG");
		addTetrahedralH2Measurements("MET", "CB", "CG", "SD");
		addLinkMeasurements(         "MET", "CB", "CG", "SD");
		addLinkMeasurements(         "MET", "CG", "SD", "CE");
		addH3Measurements(           "MET", "CG", "SD", "CE");

		// phenylalanine
		addTetrahedralH2Measurements("PHE", "CA", "CB", "CG");
		addLinkMeasurements(         "PHE", "CA", "CB", "CG");
		addBranchMeasurements(       "PHE", "CB", "CG", "CD1", "CD2");
		addTrigonalMeasurements(     "PHE", "CG", "CD1", "CE1", "HD1");
		addTrigonalMeasurements(     "PHE", "CG", "CD2", "CE2", "HD2");
		addRingLinkMeasurements(     "PHE", "CB", "CG", "CD1", "CE1");
		addRingLinkMeasurements(     "PHE", "CB", "CG", "CD2", "CE2");
		addTrigonalMeasurements(     "PHE", "CD1", "CE1", "CZ", "HE1");
		addTrigonalMeasurements(     "PHE", "CD2", "CE2", "CZ", "HE2");
		addRingLinkMeasurements(     "PHE", "CG", "CD1", "CE1", "CZ");
		addRingLinkMeasurements(     "PHE", "CG", "CD2", "CE2", "CZ");
		addTrigonalMeasurements(     "PHE", "CE1", "CZ", "CE2", "HZ");

		// tyrosine
		addTetrahedralH2Measurements("TYR", "CA", "CB", "CG");
		addLinkMeasurements(         "TYR", "CA", "CB", "CG");
		addBranchMeasurements(       "TYR", "CB", "CG", "CD1", "CD2");
		addTrigonalMeasurements(     "TYR", "CG", "CD1", "CE1", "HD1");
		addTrigonalMeasurements(     "TYR", "CG", "CD2", "CE2", "HD2");
		addRingLinkMeasurements(     "TYR", "CB", "CG", "CD1", "CE1");
		addRingLinkMeasurements(     "TYR", "CB", "CG", "CD2", "CE2");
		addTrigonalMeasurements(     "TYR", "CD1", "CE1", "CZ", "HE1");
		addTrigonalMeasurements(     "TYR", "CD2", "CE2", "CZ", "HE2");
		addRingLinkMeasurements(     "TYR", "CG", "CD1", "CE1", "CZ");
		addRingLinkMeasurements(     "TYR", "CG", "CD2", "CE2", "CZ");
		addTrigonalMeasurements(     "TYR", "CE1", "CZ", "CE2", "OH");
		addLinkMeasurements(         "TYR", "CZ", "OH", "HH");

		// tryptophan
		addTetrahedralH2Measurements("TRP", "CA", "CB", "CG");
		addLinkMeasurements(         "TRP", "CA", "CB", "CG");
		addBranchMeasurements(       "TRP", "CB", "CG", "CD1", "CD2");
		addTrigonalMeasurements(     "TRP", "CG", "CD1", "NE1", "HD1");
		addRingLinkMeasurements(     "TRP", "CB", "CG", "CD1", "NE1");
		addTrigonalMeasurements(     "TRP", "CD1", "NE1", "CE2", "HE1");
		addRingLinkMeasurements(     "TRP", "CB", "CG", "CD2", "CE2");
		addRingLinkMeasurements(     "TRP", "CB", "CG", "CD2", "CE3");
		addTrigonalMeasurements(     "TRP", "CD2", "CE3", "CZ3", "HE3");
		addRingLinkMeasurements(     "TRP", "CG", "CD1", "NE1", "CE2");
		addRingLinkMeasurements(     "TRP", "CD1", "NE1", "CE2", "CZ2");
		addTrigonalMeasurements(     "TRP", "CE2", "CZ2", "CH2", "HZ2");
		addRingLinkMeasurements(     "TRP", "NE1", "CE2", "CZ2", "CH2");
		addTrigonalMeasurements(     "TRP", "CZ2", "CH2", "CZ3", "HH2");
		addRingLinkMeasurements(     "TRP", "CG", "CD2", "CE3", "CZ3");
		addTrigonalMeasurements(     "TRP", "CE3", "CZ3", "CH2", "HZ3");
		addRingLinkMeasurements(     "TRP", "CD2", "CE3", "CZ3", "CH2");

		// serine
		addTetrahedralH2Measurements("SER", "CA", "CB", "OG");
		addLinkMeasurements(         "SER", "CA", "CB", "OG");
		addLinkMeasurements(         "SER", "CB", "OG", "HG");

		// thrionine
		addTetrahedralMeasurements(  "THR", "CA", "CB", "OG1", "HB");
		addBranchMeasurements(       "THR", "CA", "CB", "OG1", "CG2");
		addH3Measurements(           "THR", "CA", "CB", "CG2");
		addLinkMeasurements(         "THR", "CB", "OG1", "HG1");

		// aparagine
		addTetrahedralH2Measurements("ASN", "CA", "CB", "CG");
		addLinkMeasurements(         "ASN", "CA", "CB", "CG");
		addBranchMeasurements(       "ASN", "CB", "CG", "OD1", "ND2");
		addH2Measurements(           "ASN", "CB", "CG", "ND2");

		// glutamine
		addTetrahedralH2Measurements("GLN", "CA", "CB", "CG");
		addLinkMeasurements(         "GLN", "CA", "CB", "CG");
		addTetrahedralH2Measurements("GLN", "CB", "CG", "CD");
		addLinkMeasurements(         "GLN", "CB", "CG", "CD");
		addBranchMeasurements(       "GLN", "CG", "CD", "OE1", "NE2");
		addH2Measurements(           "GLN", "CG", "CD", "NE2");
	}

	private static void addLinkMeasurements(String type, String a, String b, String c) {
		mlib.add(type, new BondAngle(a, b, c));
		mlib.add(type, new BondLength(b, c));
	}

	private static void addRingLinkMeasurements(String type, String a, String b, String c, String d) {
		mlib.add(type, new DihedralAngle(a, b, c, d));
		mlib.add(type, new BondAngle(b, c, d));
		mlib.add(type, new BondLength(c, d));
	}

	private static void addTetrahedralMeasurements(String type, String a, String b, String c, String d) {
		mlib.add(type, new TetrahedralInPlaneAngle(a, b, c, d));
		mlib.add(type, new TetrahedralOutOfPlaneAngle(a, b, c, d));
		mlib.add(type, new BondLength(b, d));
	}

	private static void addTrigonalMeasurements(String type, String a, String b, String c, String d) {
		// can use tetrahedral model on trigonal geometry too
		addTetrahedralMeasurements(type, a, b, c, d);
	}

	private static void addTetrahedralH2Measurements(String type, String a, String b, String c) {
		String bpos = b.substring(1);
		String[] h = new String[] {
			"H" + bpos + "2",
			"H" + bpos + "3"
		};
		for (int i=0; i<2; i++) {
			mlib.add(type, new TetrahedralInPlaneAngles2(a, b, c, h[0], h[1], i));
			mlib.add(type, new TetrahedralOutOfPlaneAngles2(a, b, c, h[0], h[1], i));
			mlib.add(type, new BondLength(b, h[i]));
		}
	}

	private static void addH2Measurements(String type, String a, String b, String c) {
		String cpos = c.substring(1);
		String[] h = new String[] {
			"H" + cpos + "1",
			"H" + cpos + "2"
		};
		mlib.add(type, new DeltaDihedralAngle(a, b, c, h[0], h[1]));
		for (int i=0; i<2; i++) {
			mlib.add(type, new BondAngle(b, c, h[i]));
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
		mlib.add(type, new DeltaDihedralAngle(a, b, c, h[0], h[1]));
		mlib.add(type, new DeltaDihedralAngle(a, b, c, h[1], h[2]));
		mlib.add(type, new DeltaDihedralAngle(a, b, c, h[2], h[0]));
		for (int i=0; i<3; i++) {
			mlib.add(type, new BondAngle(b, c, h[i]));
			mlib.add(type, new BondLength(c, h[i]));
		}
	}

	private static void addBranchMeasurements(String type, String a, String b, String c1, String c2) {
		mlib.add(type, new BondAngle(a, b, c1));
		mlib.add(type, new BondAngle(a, b, c2));
		mlib.add(type, new BondAngle(c1, b, c2));
		mlib.add(type, new BondLength(b, c1));
		mlib.add(type, new BondLength(b, c2));
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

	public static class MeasuredRes {

		public final String filename;
		public final String resNum;
		public final double[] measurements;
		public final double rmsViolation;

		public MeasuredRes(String filename, String resNum, double[] measurements, double rmsViolation) {
			this.filename = filename;
			this.resNum = resNum;
			this.measurements = measurements;
			this.rmsViolation = rmsViolation;
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
