package edu.duke.cs.osprey;


import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.AtomPairVanDerWaalsContribution;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Streams;

import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;

public class Test {

	/* NOTES

		Differences between Osprey 2 and Osprey 3 so far:
			- harder to customize the Amber forcefield params
			- slightly different minimizers
			- slightly different implementation for WT rots

		protein energies are very similar

		custom amber params for small molecule are hugely inconsistent!!
		  results in almost all small molecule atoms having very wrong vdW params
		  no visible errors to the user!!

		Osprey 2 does not attempt to compute solvation energies for small molecules!!
	*/

	public static void main(String[] args) {

		var dir = Paths.get("/home/jeff/dlab/osprey test cases/osprey 2 thanatin/EGFR/endogenous_ligand/osprey 3");

		// pick a forcefield
		//var paramsName = "parm96a.original.dat";
		var paramsName = "parm96a.fixed.dat";
		var amberParams = FileTools.readFile(dir.resolve(paramsName).toFile());
		var ffparams = new ForcefieldParams(ForcefieldParams.Forcefield.AMBER, amberParams);

		var mol = PDBIO.read(FileTools.readFile(dir.resolve("2ITX-prep-12A.pdb").toFile()));

		// make a customized template library
		var ampTemplate =
			"\n" +
				"Ignored comment 1\n" +
				"Ignored comment 2\n" +
				"ANP\n" +
				"ANP   INT  0\n" +
				"CORRECT     OMIT DU   BEG\n" +
				"  0.0000\n" +
				"   1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000\n" +
				"   2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000\n" +
				"   3  DUMM  DU    M    2   1   0     1.523   111.21       .0      .00000\n" +
				"   4  O1G   o     M    3   2   1     1.540   111.208  -180.000 -0.642728\n" +
				"   5  PG    p5    M    4   3   2     1.515   132.522    -5.756 -0.002544\n" +
				"   6  O2G   o     E    5   4   3     1.514   112.329    59.122 -0.642728\n" +
				"   7  O3G   o     E    5   4   3     1.517   112.384   -68.646 -0.642728\n" +
				"   8  N3B   n3    M    5   4   3     1.782   105.140   173.031 -0.233631\n" +
				"   9  HN3B  hn    E    8   5   4     1.016   112.476   -52.191  0.136211\n" +
				"  10  PB    p5    M    8   5   4     1.768   119.919    84.086  0.204543\n" +
				"  11  O1B   o     E   10   8   5     1.481   110.071    74.360 -0.518931\n" +
				"  12  O2B   o     E   10   8   5     1.483   107.791  -154.598 -0.518931\n" +
				"  13  O3A   os    M   10   8   5     1.691   100.496   -40.009 -0.265932\n" +
				"  14  PA    p5    M   13  10   8     1.689   138.637   -92.719  0.267985\n" +
				"  15  O1A   o     E   14  13  10     1.478   108.453   149.678 -0.509812\n" +
				"  16  O2A   o     E   14  13  10     1.482   110.097    17.324 -0.509812\n" +
				"  17  O5'   os    M   14  13  10     1.707   101.399   -96.022 -0.318532\n" +
				"  18  C5'   c3    M   17  14  13     1.420   123.165   -42.964  0.083774\n" +
				"  19  H5'1  h1    E   18  17  14     1.091   105.940   138.663  0.060057\n" +
				"  20  H5'2  h1    E   18  17  14     1.090   111.250    24.006  0.060057\n" +
				"  21  C4'   c3    M   18  17  14     1.539   114.470  -101.842  0.112054\n" +
				"  22  C3'   c3    3   21  18  17     1.530   117.827    80.849  0.112250\n" +
				"  23  O3'   oh    S   22  21  18     1.417   113.939    79.644 -0.387451\n" +
				"  24  HO3'  ho    E   23  22  21     0.946   108.562    67.581  0.210727\n" +
				"  25  C2'   c3    B   22  21  18     1.523   101.761  -159.268  0.122247\n" +
				"  26  O2'   oh    S   25  22  21     1.418   113.174   -87.328 -0.386146\n" +
				"  27  HO2'  ho    E   26  25  22     0.945   108.556    59.363  0.210797\n" +
				"  28  H2'   h1    E   25  22  21     1.092   111.738   150.982  0.067196\n" +
				"  29  H3'   h1    E   22  21  18     1.084   110.929   -43.700  0.065853\n" +
				"  30  H4'   h1    E   21  18  17     1.092   107.935  -157.277  0.065790\n" +
				"  31  O4'   os    M   21  18  17     1.425   111.084   -40.779 -0.351695\n" +
				"  32  C1'   c3    M   31  21  18     1.398   110.105   151.830  0.142621\n" +
				"  33  H1'   h2    E   32  31  21     1.091   111.552   125.420  0.081023\n" +
				"  34  N9    n3    M   32  31  21     1.495   108.064  -119.424 -0.231091\n" +
				"  35  C8    c3    M   34  32  31     1.377   129.589    30.276  0.111551\n" +
				"  36  H8    h3    E   35  34  32     1.073   122.328    -5.298  0.068712\n" +
				"  37  N7    n3    M   35  34  32     1.305   115.773   174.873 -0.203494\n" +
				"  38  C5    c3    M   37  35  34     1.385   103.507     2.477  0.136185\n" +
				"  39  C6    c3    M   38  37  35     1.416   134.339  -179.498  0.144994\n" +
				"  40  N6    n3    B   39  38  37     1.339   123.160     3.576 -0.307466\n" +
				"  41  H61   hn    E   40  39  38     1.012   120.157  -176.009  0.120433\n" +
				"  42  H62   hn    E   40  39  38     1.011   119.721    -1.048  0.120433\n" +
				"  43  N1    n3    M   39  38  37     1.352   118.251  -175.884 -0.190967\n" +
				"  44  C2    c3    M   43  39  38     1.336   119.658    -2.300  0.125477\n" +
				"  45  H2    h3    E   44  43  39     1.078   116.094   176.793  0.070068\n" +
				"  46  N3    n3    M   44  43  39     1.333   128.911    -1.902 -0.190343\n" +
				"  47  C4    c3    M   46  44  43     1.363   111.118     3.490  0.153923\n" +
				"\n" +
				"\n" +
				"LOOP\n" +
				"  C1'  C2'\n" +
				"   C4   N9\n" +
				"   C4   C5\n" +
				"   C4   C5\n" +
				"\n" +
				"IMPROPER\n" +
				"                                                                \n" +
				"DONE                                                            \n";
		var templateLib = new ResidueTemplateLibrary.Builder(ffparams)
			.addTemplates(ampTemplate)
			.addRotamers(FileTools.readFile("GenericRotamers.dat"))
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readResource("/config/all_amino_coords.in"))
			//.addTemplateCoords(FileTools.readResource("/config/template_coords_v2.txt"))
			.build();

		// define the protein strand
		var proteinStrand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("707", "1014")
			.build();
		var allAA = Arrays.asList(
			Strand.WildType,
			"ALA", "GLY", "LEU", "ILE", "MET", "VAL", "SER", "THR", "TYR",
			"HIE", "HID", "HIP",
			"GLN", "ASN", "LYS", "ASP", "GLU", "CYS", "TRP", "ARG", "PHE"
		);
		proteinStrand.flexibility.get("718").setLibraryRotamers(allAA).addWildTypeRotamers().setContinuous();
		proteinStrand.flexibility.get("719").setLibraryRotamers(allAA).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		var ligandStrand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("2020", "2020")
			// not enough bond params defined in the forcefield to match templates with the usual method
			.setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames)
			.build();
		ligandStrand.flexibility.get("2020").addWildTypeRotamers().setContinuous();
		// no flexibility

		var protein = new SimpleConfSpace.Builder()
			.addStrand(proteinStrand)
			.build();
		var ligand = new SimpleConfSpace.Builder()
			.addStrand(ligandStrand)
			.build();
		var complex = new SimpleConfSpace.Builder()
			.addStrand(proteinStrand)
			//.addStrand(ligandStrand)
			.addStrand(ligandStrand, new StrandFlex.TranslateRotate())
			.build();

		/* TEMP
		//var confSpace = protein;
		var confSpace = ligand;

		IntFunction<Integer> findWildType = (posi) ->
			confSpace.positions.get(posi).resConfs.stream()
				.filter(rc -> rc.type == SimpleConfSpace.ResidueConf.Type.WildType)
				.findFirst()
				.orElseThrow()
				.index;

		// get the wild type conf
		int[] conf = confSpace.positions.stream()
			.mapToInt(pos -> findWildType.apply(pos.index))
			.toArray();
		*/

		// compute an energy for one conf
		try (var ecalc = new EnergyCalculator.Builder(Arrays.asList(protein, ligand, complex), ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			//.setIsMinimizing(false)
			.build()
		) {

			var epart = EnergyPartition.Traditional;

			/* TEMP
			var confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(epart)
				.build();

			// calculate the intra-shell energies
			var intraShellInters = new ResidueInteractions();
			var shell = new ArrayList<>(confSpace.shellResNumbers);
			for (int i=0; i<shell.size(); i++) {
				for (int j=0; j<=i; j++) {
					intraShellInters.addPair(shell.get(i), shell.get(j));
				}
			}
			var pmol = confSpace.makeMolecule(new int [] {});
			double intraShellEnergy = ecalc.calcEnergy(pmol, intraShellInters).energy;
			log("intra-shell energy: %f", intraShellEnergy);

			// calc the energy
			var tuple = new RCTuple(conf);
			double energy = confEcalc.calcEnergy(tuple).energy;
			log("conf = %s   energy = %f", Conf.toString(conf), energy);

			var modConfEcalc = new IntraShellConfEcalc(confSpace, ecalc, epart);
			double modEnergy = modConfEcalc.calcEnergy(tuple).energy;
			log("conf = %s   energy w/ intra-shell = %f", Conf.toString(conf), modEnergy);
			*/

			/* TEMP

			// show the RCs
			for (int posi=0; posi<confSpace.numPos(); posi++) {
				var rc = confSpace.positions.get(posi).resConfs.get(conf[posi]);
				log("\tpos %d: %s",
					posi,
					IntStream.range(0, rc.template.numDihedrals)
						.mapToObj(dihedrali -> {

							String atomNames = Arrays.stream(rc.template.getDihedralDefiningAtoms(dihedrali))
								.mapToObj(atomi -> rc.template.templateRes.atoms.get(atomi).name)
								.collect(Collectors.joining(","));

							double degrees = rc.template.getRotamericDihedrals(rc.rotamerIndex, dihedrali);

							return String.format("%s:%s=%.2f", rc.template.name, atomNames, degrees);
						})
						.collect(Collectors.joining(", "))
				);
			}

			// break down the energies by residue
			var analysis = new ConfAnalyzer(modConfEcalc).analyze(conf);
			var analysisByRes = analysis.breakdownEnergyByResidue();
			//log("conf energies: %f\n%s",
			//	analysis.epmol.energy,
			//	analysisByRes.breakdownForcefield(ResidueForcefieldBreakdown.Type.VanDerWaals).toStringScientific()
			//);
			log("\t%8s = %f\n\t%8s = %f\n\t%8s = %f",
				"es", analysisByRes.efunc.getElectrostaticsEnergy(),
				"vdW", analysisByRes.efunc.getVanDerWaalsEnergy(),
				"solv", analysisByRes.efunc.getSolvationEnergy()
			);
			*/

			/* TEMP
			// break down the energies by atom
			ResPairCache.ResPair[] resPairs = Arrays.stream(analysisByRes.efunc.resPairs)
				.filter(resPair ->
					//resPair.res1.getPDBResNumber().equals("722")
					//&& resPair.res2.getPDBResNumber().equals("722")
					true
				)
				.toArray(ResPairCache.ResPair[]::new);
			var pairs = analysisByRes.efunc.getVanDerWaalsEnergyContributions(resPairs).stream()
				.flatMap(resPair -> resPair.getAtomPairEnergyContributions().stream())
				.filter(pair ->
					true
					//pair.getEnergy() > 0.0
				)
				.collect(Collectors.toList());
			Function<Atom,String> atomNamer = atom -> String.format("%s:%s", atom.res.fullName, atom.name);
			Collections.sort(pairs, Comparator.comparing(pair -> {
				var a = atomNamer.apply(pair.getAtom1());
				var b = atomNamer.apply(pair.getAtom2());
				if (a.compareTo(b) > 0) {
					var swap = a;
					a = b;
					b = swap;
				}
				return String.format("%s %s", a, b);
			}));
			//final int numPairs = 1000;
			final int numPairs = Integer.MAX_VALUE;
			int n = 0;
			for (var p : pairs) {
				var pair = (AtomPairVanDerWaalsContribution)p;
				var label1 = atomNamer.apply(pair.getAtom1());
				var label2 = atomNamer.apply(pair.getAtom2());
				if (label1.compareTo(label2) > 0) {
					var l = label1;
					label1 = label2;
					label2 = l;
				}
				log("%14s - %-14s   t %s(%d=%s)-%s(%d=%s)   r %.6f   e %.6f   0 %.6f   1 %.6f",
					label1,
					label2,
					pair.getAtom1().forceFieldType, pair.getAtom1().type, ffparams.atomType(pair.getAtom1().type),
					pair.getAtom2().forceFieldType, pair.getAtom2().type, ffparams.atomType(pair.getAtom2().type),
					Math.sqrt(pair.getR2()),
					pair.getEnergy(),
					pair.getAij(), pair.getBij()
				);
				if (++n >= numPairs) {
					break;
				}
			}
			*/

			/* TEMP
			// print out the residue cheat sheet
			for (int i=0; i<analysisByRes.efunc.residues.size(); i++) {
				var res = analysisByRes.efunc.residues.get(i);
				log("res %2d = %s %s", i, res.getPDBResNumber(), res.template.name);
			}
			*/

			// run K*
			/*
			var kstarSettings = new KStar.Settings.Builder()
				.setEpsilon(0.3)
				.setMaxSimultaneousMutations(0)
				.setShowPfuncProgress(true)
				.addScoreConsoleWriter(info -> {

					// focus just on the protein state
					return String.format("sequence [%s]:\n\tprotein: [%e,%e]\n\tligand: [%e,%e]\n\tcomplex: [%e,%e]",
						info.sequence,
						info.kstarScore.protein.values.calcLowerBound(), info.kstarScore.protein.values.calcUpperBound(),
						info.kstarScore.ligand.values.calcLowerBound(), info.kstarScore.ligand.values.calcUpperBound(),
						info.kstarScore.complex.values.calcLowerBound(), info.kstarScore.complex.values.calcUpperBound()
					);
				})
				.build();
			*/

			var kstarSettings = new KStar.Settings.Builder()
				.setEpsilon(0.0)
				.setMaxSimultaneousMutations(1)
				.setStabilityThreshold(null)
				.addScoreWriter(new KStarScoreWriter.ToConsole(info ->
					String.format("sequence %4d/%4d   %s   : %s in [%s,%s]   protein %d:%s   ligand %d:%s   complex %d:%s",
						info.sequenceNumber + 1,
						info.numSequences,
						info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1),
						edu.duke.cs.osprey.tools.Log.formatBigEngineering(info.kstarScore.score),
						edu.duke.cs.osprey.tools.Log.formatBigEngineering(info.kstarScore.lowerBound),
						edu.duke.cs.osprey.tools.Log.formatBigEngineering(info.kstarScore.upperBound),
						info.kstarScore.protein.numConfs,
						edu.duke.cs.osprey.tools.Log.formatBigEngineering(info.kstarScore.protein.values.calcLowerBound()),
						info.kstarScore.ligand.numConfs,
						edu.duke.cs.osprey.tools.Log.formatBigEngineering(info.kstarScore.ligand.values.calcLowerBound()),
						info.kstarScore.complex.numConfs,
						edu.duke.cs.osprey.tools.Log.formatBigEngineering(info.kstarScore.complex.values.calcLowerBound())
					)
				))
				.build();

			var kstar = new KStar(protein, ligand, complex, kstarSettings);
			for (var info : kstar.confSpaceInfos()) {

				var cs = (SimpleConfSpace)info.confSpace;
				//info.confEcalc = new ConfEnergyCalculator.Builder(cs, ecalc).build();
				info.confEcalc = new IntraShellConfEcalc(cs, ecalc, epart);

				// compute the energy matrix
				var emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

				// make a pfunc calculator
				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalc,
					new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
					new ConfAStarTree.Builder(emat, rcs).setTraditional().build(),
					rcs.getNumConformations()
				);

				// cleanup the confDB
				info.confDBFile.delete();
			}
			kstar.run(ecalc.tasks);
		}

		log("Done! =)");
	}

	// TEMP: try to mod the confEcalc to add the intra-shell energies
	// TODO: add an option to the real conf energy calculator to add intra-shell inters?
	static class IntraShellConfEcalc extends ConfEnergyCalculator {

		public IntraShellConfEcalc(SimpleConfSpace confSpace, EnergyCalculator ecalc, EnergyPartition epart) {
			super(confSpace, ecalc, epart, null, false, null, 0.0);
		}

		@Override
		public ResidueInteractions makeFragInters(RCTuple frag) {

			var inters = EnergyPartition.makeFragment(confSpace, eref, addResEntropy, frag);

			// add the intra-shell inters too
			var shell = new ArrayList<>(confSpace.shellResNumbers);
			for (int i=0; i<shell.size(); i++) {
				for (int j=0; j<=i; j++) {
					inters.addPair(shell.get(i), shell.get(j));
				}
			}

			return inters;
		}
	}
}
