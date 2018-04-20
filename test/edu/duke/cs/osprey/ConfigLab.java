package edu.duke.cs.osprey;

import com.google.common.collect.Iterators;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;


public class ConfigLab {

	private static final File inDir = new File("/home/jeff/Desktop/temp");
	private static final File outDir = new File("/home/jeff/Desktop/tempout");

	private static StringBuilder logbuf = null;

	public static void main(String[] args) {

		// these look tricky, skip them for now
		Set<String> skips = new HashSet<>();
		skips.add("3U7Y_KS");
		skips.add("3U7Y_KS_tr");

		for (String code : inDir.list()) {
			if (skips.contains(code)) {
				continue;
			}
			System.out.println(code);
			go(code);
		}
	}

	public static File whicheverExists(File ... files) {
		for (File file : files) {
			if (file.exists()) {
				return file;
			}
		}
		throw new Error("no files exist: " + Arrays.toString(files));
	}

	public static String[] tokenize(String line) {
		return Arrays.stream(line.split("\\s"))
			.filter((part) -> !part.isEmpty())
			.toArray((size) -> new String[size]);
	}

	public static String[] getMutations(ConfigFileParser cfp, int strandIndex, int i) {

		String nope = "__NOPE__";

		String mutsString = cfp.params.getValue("resAllowed" + strandIndex + "_" + i, nope);
		if (mutsString == nope && strandIndex == 0) {
			mutsString = cfp.params.getValue("resAllowed" + i, nope);
		}
		if (mutsString == nope) {
			throw new Error("can't find mutations for strand " + strandIndex + " at pos " + i);
		}

		return tokenize(mutsString);
	}

	public static String findResNum(Molecule mol, String resNum) {

		for (Residue res : mol.residues) {
			String num = res.fullName.substring(5).trim();
			if (resNum.equalsIgnoreCase(num)) {
				return res.getPDBResNumber();
			}
		}

		throw new Error("can't find res num: " + resNum);
	}

	public static void go(String code) {

		logbuf = new StringBuilder();

		File inDir = new File(ConfigLab.inDir, code);

		// concat the config files
		ConfigFileParser cfp = ConfigFileParser.makeFromFiles(Arrays.asList(
			new File(inDir, "KStar.cfg"),
			new File(inDir, "System.cfg"),
			whicheverExists(new File(inDir, "DEE.cfg"), new File(inDir, "MutSearch.cfg"))
		));

		// read some general settings
		boolean continuousFlex = cfp.params.getBool("doMinimize");
		boolean addWtAA = cfp.params.getBool("addWT");
		boolean addWtRots = cfp.params.getBool("addWTRots");
		boolean useEref = cfp.params.getBool("useEref");
		double distCutoff = cfp.params.getDouble("distCutoff", cfp.params.getDouble("shellDistCutoff"));

		logln("");
		logln("import osprey");
		logln("osprey.start()");
		logln("");
		logln("stopwatch = osprey.c.tools.Stopwatch().start()");

		String pdbFilename = cfp.params.getValue("pdbName");
		Molecule mol = PDBIO.readFile(new File(inDir, pdbFilename));

		logln("mol = osprey.readPdb('%s')", pdbFilename);

		int numStrands = cfp.params.getInt("numOfStrands");
		for (int strandIndex=0; strandIndex<numStrands; strandIndex++) {

			// read the strand residues
			String[] parts = tokenize(cfp.params.getValue("strand" + strandIndex));
			String startRes = findResNum(mol, parts[0]);
			String stopRes = findResNum(mol, parts[1]);

			// make sure they actually exist
			mol.residues.getOrThrow(startRes);
			mol.residues.getOrThrow(stopRes);

			logln("strand%d = osprey.Strand(mol, residues=['%s', '%s'])", strandIndex, startRes, stopRes);

			// read the flexibility
			String[] resNums = tokenize(cfp.params.getValue("strandMut" + strandIndex));
			for (int i=0; i<resNums.length; i++) {

				String[] mutations = getMutations(cfp, strandIndex, i);

				String resNum = findResNum(mol, resNums[i]);
				log("strand%d.flexibility['%s'].setLibraryRotamers(", strandIndex, resNum);
				if (addWtAA) {
					log("osprey.WILD_TYPE");
				}

				for (int m=0; m<mutations.length; m++) {
					if (m > 0 || addWtAA) {
						log(", ");
					}
					log("'%s'", mutations[m]);
				}
				log(")");

				if (addWtRots) {
					log(".addWildTypeRotamers()");
				}

				if (continuousFlex) {
					log(".setContinuous()");
				}

				logln("");
			}
		}

		// make the conf space
		log("confSpace = osprey.ConfSpace([");
		for (int strandIndex=0; strandIndex<numStrands; strandIndex++) {
			if (strandIndex > 0) {
				log(", ");
			}
			log("strand%d", strandIndex);
		}
		log("]");
		if (Double.isFinite(distCutoff)) {
			log(", shellDist=%.3f", distCutoff);
		}
		logln(")");

		logln("");
		logln("ffparams = osprey.ForcefieldParams()");
		logln("ecalc = osprey.EnergyCalculator(confSpace, ffparams)");
		if (useEref) {
			logln("eref = osprey.ReferenceEnergies(confSpace, ecalc)");
			logln("confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)");
		} else {
			logln("confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)");
		}
		logln("emat = osprey.EnergyMatrix(confEcalc)");
		logln("astar = osprey.AStarTraditional(emat, confSpace)");
		logln("");

		// read the run script to see if this is a GMEC or a pfunc run
		File execFile = whicheverExists(
			new File(inDir, code + ".q"),
			new File(inDir, code.substring(0, Math.max(0, code.indexOf('_'))) + ".q"),
			new File(inDir, code.substring(0, Math.max(0, code.indexOf('+'))) + ".q")
		);
		String execLine = Iterators.getLast(FileTools.parseLines(FileTools.readFile(execFile)).iterator());
		if (execLine.contains("doDEE")) {

			double Ew = cfp.params.getDouble("initEw");

			logln("# GMEC");
			logln("osprey.GMECFinder(astar, confEcalc).find(%f)", Ew);

		} else if (execLine.contains("KSMaster")) {

			logln("# PFUNC");

			// I think these are all single-sequence partition functions, right?
			for (int strandIndex=0; strandIndex<numStrands; strandIndex++) {
				int numPos = tokenize(cfp.params.getValue("strandMut" + strandIndex)).length;
				for (int i=0; i<numPos; i++) {
					String[] mutations = getMutations(cfp, strandIndex, i);
					if (mutations.length > 0) {
						throw new Error("multiple mutations possible here, need to check");
					}
				}
			}

			double epsilon = cfp.params.getDouble("epsilon");

			logln("pfunc = osprey.c.kstar.pfunc.GradientDescentPfunc(astar, confEcalc)");
			logln("pfunc.setReportProgress(True)");
			logln("pfunc.init(%.4f)", epsilon);
			logln("pfunc.compute(1000000000)");

		} else {
			throw new Error("unknown exec: " + execLine);
		}

		logln("");
		logln("print('\\n\\nfinished %s in %%s ms\\n' %% stopwatch.stop().getTimeMs())", code);


		System.out.println("\n\n");
		System.out.println(logbuf.toString());

		File outDir = new File(ConfigLab.outDir, code);
		outDir.mkdirs();
		FileTools.writeFile(logbuf.toString(), new File(outDir, "run.py"));
		PDBIO.writeFile(mol, new File(outDir, pdbFilename));
	}

	private static void log(String msg, Object ... args) {
		logbuf.append(String.format(msg, args));
	}

	private static void logln(String msg, Object ... args) {
		log(msg, args);
		logbuf.append("\n");
	}
}
