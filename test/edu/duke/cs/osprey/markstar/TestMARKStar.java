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

package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestBBKStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.TestKStar.ConfSpaces;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.visualizer.KStarTreeManipulator;
import edu.duke.cs.osprey.markstar.visualizer.KStarTreeNode;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.math.BigInteger;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.kstar.TestBBKStar.runBBKStar;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.assertThat;

//import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;

public class TestMARKStar {

	public static final int NUM_CPUs = 4;
	public static boolean REUDCE_MINIMIZATIONS = true;
	public static final EnergyPartition ENERGY_PARTITION = EnergyPartition.Traditional;

	public static class Result {
		public MARKStar markstar;
		public List<MARKStar.ScoredSequence> scores;
	}


	protected  static void TestComparison(TestKStar.ConfSpaces confSpaces, double epsilon, boolean runkstar) {
		Stopwatch runtime = new Stopwatch().start();
		String kstartime = "(Not run)";
		List<KStar.ScoredSequence> seqs = null;
		if(runkstar) {
			seqs = runKStar(confSpaces, epsilon);
			runtime.stop();
			kstartime = runtime.getTime(2);
			runtime.reset();
			runtime.start();
		}
		Result result = runMARKStar(confSpaces, epsilon);
		runtime.stop();
		String markstartime = runtime.getTime(2);
		for(MARKStar.ScoredSequence seq: result.scores)
			printMARKStarComputationStats(seq);
		if(seqs != null)
		for(KStar.ScoredSequence seq: seqs)
			printKStarComputationStats(seq);
		System.out.println("MARK* time: "+markstartime+", K* time: "+kstartime);
	}


	protected  static ConfSpaces makeConfSpaces(String pdb, String[] proteinDef, String[] proteinFlex, String[] ligandDef,
								String[] ligandFlex) {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile(pdb);

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(proteinDef[0], proteinDef[1])
			.build();
		for(String resName: proteinFlex)
            protein.flexibility.get(resName).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(ligandDef[0], ligandDef[1])
			.build();
		for(String resName: ligandFlex)
			ligand.flexibility.get(resName).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		return confSpaces;
	}



    protected  static TestKStar.ConfSpaces loadSSFromCFS(String cfsFileName) throws FileNotFoundException {
	    String fileContents = FileTools.readFile(cfsFileName);
	    TestKStar.ConfSpaces confSpaces = new ConfSpaces();
		confSpaces.ffparams = new ForcefieldParams();
		Map<String, Strand.Builder> strandBuilderMap = new HashMap<>();
		Map<String, Strand> strandMap = new HashMap<>();
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		for (String line : FileTools.parseLines(fileContents)) {
			String[] parts = line.split(" = ");
		    String lineType = parts[0];
		    switch(lineType) {
				case "mol":
					String pdbName = parts[1].substring(parts[1].lastIndexOf('/')+1, parts[1].length()-1);
					Molecule mol = PDBIO.readFile("examples/python.KStar/"+pdbName);
					strandBuilderMap.put("strand0", new Strand.Builder(mol).setTemplateLibrary(templateLib));
					strandBuilderMap.put("strand1", new Strand.Builder(mol).setTemplateLibrary(templateLib));
					break;
				case "strand_defs":
					Pattern definitionPattern = Pattern.compile("\\'(strand\\d)\\': \\[(\\'\\w+\\', ?\\'\\w+\\')\\]");
					Matcher matcher = definitionPattern.matcher(parts[1]);
					while(matcher.find()) {
						String strandNum = matcher.group(1);
						String strandRange = matcher.group(2);
						strandRange = strandRange.replaceAll("'","");
						String[] startAndEnd = strandRange.split(", ?");
						Strand.Builder strandBuilder = strandBuilderMap.get(strandNum);
						strandBuilder.setResidues(startAndEnd[0],startAndEnd[1]);
						strandMap.put(strandNum, strandBuilder.build());
						//Strand 0: protein
						//Strant 1: ligand
					}
					break;
				case "strand_flex":
					String content = parts[1];
					boolean match0 = content.matches(".*\\'strand\\d\\'.*");
					boolean match1 = content.matches(".*'strand\\d': ?\\{'\\w+': ?\\['\\w+',.*");
					boolean match2 = content.matches(".*'strand\\d': ?\\{'\\w+': ?\\[(('\\w+', ?)*)'\\w+'\\].*");
					definitionPattern = Pattern.compile("'(strand\\d)': ?\\{(('\\w+': ?\\[(('\\w+', ?)*)'\\w+'], ?)*'\\w+': ?\\[(('\\w+', ?)*)'\\w+']+\\s*)}");
					matcher = definitionPattern.matcher(parts[1]);
					while(matcher.find()) {
					    String strandNum = matcher.group(1);
					    String strandResDefs = matcher.group(2);
					    Pattern resDefPattern = Pattern.compile("'(\\w+)': ?\\[(('\\w+', ?)*'\\w+')]");
					    Matcher resDefMatcher = resDefPattern.matcher(strandResDefs);
					    while(resDefMatcher.find()) {
							String resNum = resDefMatcher.group(1);
							strandMap.get(strandNum).flexibility.get(resNum)
									.setLibraryRotamers(Strand.WildType)
									.setContinuous()
									.addWildTypeRotamers();
						}
						//Strand 0: protein
						//Strant 1: ligand
					}
					break;
			}
		}
		confSpaces.protein = new SimpleConfSpace.Builder().addStrand(strandMap.get("strand0")).build();
		confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(strandMap.get("strand1")).build();
		confSpaces.complex= new SimpleConfSpace.Builder().addStrands(strandMap.get("strand0"),
				strandMap.get("strand1")).build();


	    return confSpaces;
	}

	@Test
	public void testConfSpaceParse() {
		try {
			loadFromCFS("examples/python.KStar/2rl0_A_13res_6.837E+28.cfs");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}


	protected static TestKStar.ConfSpaces loadFromCFS(String cfsFileName) throws FileNotFoundException {
	    String fileContents = FileTools.readFile(cfsFileName);
	    TestKStar.ConfSpaces confSpaces = new ConfSpaces();
		confSpaces.ffparams = new ForcefieldParams();
		Map<String, Strand.Builder> strandBuilderMap = new HashMap<>();
		Map<String, Strand> strandMap = new HashMap<>();
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		for (String line : FileTools.parseLines(fileContents)) {
			String[] parts = line.split(" = ");
		    String lineType = parts[0];
		    switch(lineType) {
				case "mol":
					String pdbName = parts[1].substring(parts[1].lastIndexOf('/')+1, parts[1].length()-1);
					Molecule mol = PDBIO.readFile("examples/python.KStar/"+pdbName);
					strandBuilderMap.put("strand0", new Strand.Builder(mol).setTemplateLibrary(templateLib));
					strandBuilderMap.put("strand1", new Strand.Builder(mol).setTemplateLibrary(templateLib));
					break;
				case "strand_defs":
					Pattern definitionPattern = Pattern.compile("\\'(strand\\d)\\': \\[(\\'\\w+\\', ?\\'\\w+\\')\\]");
					Matcher matcher = definitionPattern.matcher(parts[1]);
					while(matcher.find()) {
						String strandNum = matcher.group(1);
						String strandRange = matcher.group(2);
						strandRange = strandRange.replaceAll("'","");
						String[] startAndEnd = strandRange.split(", ?");
						Strand.Builder strandBuilder = strandBuilderMap.get(strandNum);
						strandBuilder.setResidues(startAndEnd[0],startAndEnd[1]);
						strandMap.put(strandNum, strandBuilder.build());
						//Strand 0: protein
						//Strant 1: ligand
					}
					break;
				case "strand_flex":
					String content = parts[1];
					boolean match0 = content.matches(".*\\'strand\\d\\'.*");
					boolean match1 = content.matches(".*'strand\\d': ?\\{'\\w+': ?\\['\\w+',.*");
					boolean match2 = content.matches(".*'strand\\d': ?\\{'\\w+': ?\\[(('\\w+', ?)*)'\\w+'\\].*");
					definitionPattern = Pattern.compile("'(strand\\d)': ?\\{(('\\w+': ?\\[(('\\w+', ?)*)'\\w+'], ?)*'\\w+': ?\\[(('\\w+', ?)*)'\\w+']+\\s*)}");
					matcher = definitionPattern.matcher(parts[1]);
					while(matcher.find()) {
					    String strandNum = matcher.group(1);
					    String strandResDefs = matcher.group(2);
					    Pattern resDefPattern = Pattern.compile("'(\\w+)': ?\\[(('\\w+', ?)*'\\w+')]");
					    Matcher resDefMatcher = resDefPattern.matcher(strandResDefs);
					    while(resDefMatcher.find()) {
							String resNum = resDefMatcher.group(1);
							String allowedAAs = resDefMatcher.group(2);
							allowedAAs = allowedAAs.replaceAll("'","");
							strandMap.get(strandNum).flexibility.get(resNum)
									.setLibraryRotamers(allowedAAs.split(", ?"))
									.addWildTypeRotamers()
									.setContinuous();
						}
						//Strand 0: protein
						//Strant 1: ligand
					}
					break;
			}
		}
		confSpaces.protein = new SimpleConfSpace.Builder().addStrand(strandMap.get("strand0")).build();
		confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(strandMap.get("strand1")).build();
		confSpaces.complex= new SimpleConfSpace.Builder().addStrands(strandMap.get("strand0"),
				strandMap.get("strand1")).build();


	    return confSpaces;
	}


	@Test
    public void test2RL0DEEper() {
		ConfSpaces confSpaces = make2RL0DEEPer();
		runMARKStar(confSpaces, 0.99);

    }

    private ConfSpaces make2RL0DEEPer(){

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/2RL0.min.reduce.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();


		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("G649", "G654")
				.build();
//        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
//		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
//		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
//		protein.flexibility.get("G652").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
//		protein.flexibility.get("G653").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
//		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ArrayList<String> bbflexlist = new ArrayList<>();
		for(int startIndex = 649; startIndex < 654; startIndex++) {
            protein.flexibility.get("G"+startIndex).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
            bbflexlist.add("G"+startIndex);
		}



		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("155", "194")
				.build();
		ligand.flexibility.get("A189").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A190").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A191").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A194").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		String perturbationFileName = "examples/python.KStar/STR0.2rl0.pert.pert";
		DEEPerSettings deepersettings = new DEEPerSettings(true, perturbationFileName,
				false, perturbationFileName, true, 2.5,
				2.5, false, bbflexlist, "examples/python.KStar/shell.pdb", false, templateLib);
		DEEPerStrandFlex protein_bbflex = new DEEPerStrandFlex(protein, deepersettings);

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein, protein_bbflex)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrand(protein, protein_bbflex)
				.addStrand(ligand)
				.build();

		return confSpaces;
	}

	@Test
	public void testConsolidateTree() {
		KStarTreeNode root = KStarTreeNode.parseTree("ComplexConfTreeBounds.txt");
		List<String> proteinResidues = Arrays.stream(new String[]{
				"A177",
				"A178",
				"A179",
				"A180"
		}).collect(Collectors.toList());
		KStarTreeNode newRoot = KStarTreeManipulator.consolidateTree(root, proteinResidues);
		newRoot.printTree();
	}

	@Test
	public void test1UBQDEEPer() {
	    ConfSpaces confSpaces = make1UBQDEEPer();
	    runMARKStar(confSpaces, 0.9);

	}

	private ConfSpaces make1UBQDEEPer() {
		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readFile("examples/python.KStar/1ubqdefault.prepared.minimized.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

//		A6
//		A7
//		A8
//		A43
//		A44
//		A67
//		A68
//		A69


		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A1", "A65")
				.build();
		protein.flexibility.get("A3").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A4").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A5").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A6").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A7").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A8").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A9").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A10").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A11").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A12").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A13").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A14").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A15").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A43").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A44").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ArrayList<String> bbflexlist = new ArrayList<>();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A66", "A71")
				.build();
		ligand.flexibility.get("A67").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A68").setLibraryRotamers("ARG").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A69").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		String perturbationFileName = "examples/python.KStar/1ubqLoopStabilization.pert";
		DEEPerSettings deepersettings = new DEEPerSettings(true, perturbationFileName,
				false, perturbationFileName, true, 2.5,
				2.5, false, bbflexlist, "examples/python.KStar/1ubqdefault.prepared.minimized.pdb", false, templateLib);
		DEEPerStrandFlex protein_bbflex = new DEEPerStrandFlex(protein, deepersettings);
		//CATSStrandFlex protein_contflex = new CATSStrandFlex(protein, "A3", "A14");

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein, protein_bbflex)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrand(protein, protein_bbflex)
				.addStrand(ligand)
				.build();

		return confSpaces;
	}


	@Test
	public void testComputeTreeError() {
		KStarTreeNode proteinRoot = KStarTreeNode.parseTree("2XXM10ResComplexConfTreeBounds.txt");
		proteinRoot.preprocess();
		double initialDiff = 0.1;
		for(double i = initialDiff; i < 10; i+= 0.1) {
			double diffFromGMEC = i;
			double[] range = proteinRoot.computeEnergyErrorWithinEnergyRange(diffFromGMEC);
			int numConfs = proteinRoot.numConfsWithin(diffFromGMEC);
			System.out.println("For "+numConfs+" confs within " + diffFromGMEC + " kcal/mol of GMEC: Min: " + range[0] + ", max: " + range[1]);
		}

	}

	@Test
	public void testComputeCrossTreeEnthalpy() {
		KStarTreeNode proteinRoot = KStarTreeNode.parseTree("ProteinConfTreeBounds.txt");
		int proteinLevel = 6;
		int ligandLevel = 4;
		System.out.println("Protein enthalpy: "+proteinRoot.computeEnthalpy(proteinLevel));
		System.out.println("Protein entropy: "+proteinRoot.computeEntropy(proteinLevel));
		System.out.println("Protein states: "+proteinRoot.numStatesAtLevel(proteinLevel));

		KStarTreeNode ligandRoot = KStarTreeNode.parseTree("LigandConfTreeBounds.txt");
		System.out.println("Ligand enthalpy: "+ligandRoot.computeEnthalpy(ligandLevel));
		System.out.println("Ligand entropy: "+ligandRoot.computeEntropy(ligandLevel));
		System.out.println("Ligand states: "+ligandRoot.numStatesAtLevel(ligandLevel));

		KStarTreeNode complexRoot = KStarTreeNode.parseTree("2XXM10ResComplexConfTreeBounds.txt");
		System.out.println("Complex enthalpy: "+complexRoot.computeEnthalpy(proteinLevel));
		System.out.println("Complex (protein) entropy: "+complexRoot.computeEntropy(proteinLevel));
		System.out.println("Complex (protein) states: "+complexRoot.numStatesAtLevel(proteinLevel));
		KStarTreeNode complexLigandFirstRoot = KStarTreeNode.parseTree("2XXM10ResComplexLigandFirstBounds.txt");
		System.out.println("Complex (ligand) entropy: "+complexLigandFirstRoot.computeEntropy(ligandLevel));
		System.out.println("Complex (ligand) states: "+complexLigandFirstRoot.numStatesAtLevel(ligandLevel));
		double crossTreeEnthalpy = complexRoot.computeEnthalpyWithEnergiesFrom(proteinRoot.computeEnergyMap(proteinLevel), proteinLevel);
		System.out.println("Cross tree protein enthalpy: "+crossTreeEnthalpy);
		crossTreeEnthalpy = complexLigandFirstRoot.computeEnthalpyWithEnergiesFrom(ligandRoot.computeEnergyMap(ligandLevel), ligandLevel);
		System.out.println("Cross tree ligand enthalpy: "+crossTreeEnthalpy);
		System.out.println("Complex full entropy:"+complexRoot.computeEntropy());
		System.out.println("Complex full enthalpy:"+complexRoot.computeEnthalpy());
		System.out.println("Complex max conf error bound: "+complexRoot.maxConfErrorBound());
		System.out.println(String.format("Complex max weighted error bound: %12.4e, largest upper bound %12.4e",
				complexRoot.maxWeightedErrorBound(),
				new BoltzmannCalculator(PartitionFunction.decimalPrecision).calc(complexRoot.getConfLowerBound())));


	}

	@Test
	public void testGenerateEnsemble() {
		ConfSpaces confSpaces = make2XXM10Res();
        KStarTreeNode root = KStarTreeNode.parseTree("2XXM10ResComplexConfTreeBounds.txt");
        int numConfs = 150;
        int levelThreshold = 4;
        Map<KStarTreeNode, List<KStarTreeNode>> samples = root.getTopSamples(numConfs, levelThreshold);
        Set<KStarTreeNode> levelNodes = root.getLevelNodes(1);
        System.out.println("Tried for "+numConfs+" confs, got "+samples.size()+" lists");
        for(KStarTreeNode subtreeRoot:samples.keySet()) {
			System.out.println("Under " + subtreeRoot + ":");
			for (KStarTreeNode conf : samples.get(subtreeRoot)) {
				System.out.println(conf.toString());
			}
		}
		Map<KStarTreeNode, List<ConfSearch.ScoredConf>> organizedConfs = new HashMap<>();
        Map<ConfSearch.ScoredConf, String> confNames = new HashMap<>();
        for(KStarTreeNode levelNode: levelNodes)
        	organizedConfs.put(levelNode, new ArrayList<>());

		boolean printPDBs = true;
		if(printPDBs) {
			EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
					.setParallelism(new Parallelism(4, 0, 0))
					.build();
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpaces.complex, minimizingEcalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaces.complex, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();

			// calc energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.setCacheFile(new File("2XXM.ensemble.emat"))
					.build()
					.calcEnergyMatrix();
			ConfAnalyzer analyzer = new ConfAnalyzer(confEcalc);
			Map<KStarTreeNode, List<ConfSearch.ScoredConf>> confLists = new HashMap<>();
			for (KStarTreeNode subtreeRoot : samples.keySet()) {
				System.out.println("Under " + subtreeRoot + ":");
				confLists.put(subtreeRoot, new ArrayList<>());
				for (KStarTreeNode conf : samples.get(subtreeRoot)) {
					ConfSearch.ScoredConf scoredConf = new ConfSearch.ScoredConf(conf.getConfAssignments(), conf.getConfLowerBound());
					confLists.get(subtreeRoot).add(scoredConf);
					confNames.put(scoredConf, subtreeRoot.getEnsemblePDBName());
					for(KStarTreeNode levelNode:levelNodes) {
						if(conf.isChildOf(levelNode)) {
							organizedConfs.get(levelNode).add(scoredConf);
						}
					}
				}
			}

			File pdbDir = new File("pdb");
			if(!pdbDir.exists())
                pdbDir.mkdir();
			for (KStarTreeNode subtreeRoot : organizedConfs.keySet()) {
				System.out.println("Under " + subtreeRoot + ":");
				File subDir = new File("pdb/"+subtreeRoot.getEnsemblePDBName());
				if(!subDir.exists())
					subDir.mkdir();
				ConfAnalyzer.ConfAnalysis analysis2 = analyzer.analyze(organizedConfs.get(subtreeRoot).get(0));
				System.out.println(analysis2);
				System.out.println(analysis2.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All));
				ConfAnalyzer.EnsembleAnalysis analysis = analyzer.analyzeEnsemble(organizedConfs.get(subtreeRoot).iterator(), Integer.MAX_VALUE);
				analysis.writePdbs("pdb/"+subtreeRoot.getEnsemblePDBName()+"/"
						+"*.pdb");

			}
		}
	}

	@Test
	public void test3u7y() {
		try {
			ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/3u7y_L_15res_1.326E+48.cfs");
			runMARKStar(confSpaces, 0.9999);
			//runBBKStar(confSpaces, 5, 0.99, null, 5, true);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}


	@Test
	public void test2XXMSmaller() {
		ConfSpaces confSpaces = make2XXMSmaller();
		final double epsilon = 0.01;
		String kstartime = "(not run)";
		boolean runkstar = false;
		Stopwatch runtime = new Stopwatch().start();
		List<KStar.ScoredSequence> kStarSeqs = null;
		if(runkstar) {
			kStarSeqs = runKStar(confSpaces, epsilon);
			runtime.stop();
			kstartime = runtime.getTime(2);
			runtime.reset();
			runtime.start();
		}
		Result result = runMARKStar(confSpaces, epsilon);
		runtime.stop();
		String markstartime = runtime.getTime(2);
		System.out.println("MARK* time: "+markstartime+", K* time: "+kstartime);
		for(MARKStar.ScoredSequence seq: result.scores)
			printMARKStarComputationStats(seq);
		if(runkstar)
			for(KStar.ScoredSequence seq: kStarSeqs)
				printKStarComputationStats(seq);
	}

	@Test
	public void testApoPfuncs() {
		ConfSpaces proteinConfspaces = make3DS210Res();
		ConfSpaces ligandConfspaces = make2XXC10Res();
		ConfSpaces mergedConfspaces = new ConfSpaces();
		mergedConfspaces.ffparams = new ForcefieldParams();
		mergedConfspaces.protein = proteinConfspaces.protein;
		mergedConfspaces.ligand = ligandConfspaces.protein;
		mergedConfspaces.complex = new SimpleConfSpace.Builder()
				.addStrands(proteinConfspaces.protein.strands.get(0),
						ligandConfspaces.protein.strands.get(0))
				.build();
		runMARKStar(mergedConfspaces, 0.01);

	}

	private ConfSpaces make3DS210Res() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/3ds2.prepared.minimized.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A152", "A216")
				.build();
		protein.flexibility.get("A177").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A180").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A181").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A184").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A185").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A219", "A220")
				.build();
		ligand.flexibility.get("A219").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();


		return confSpaces;

	}

	private ConfSpaces make1A8O10Res() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/1a8o_mutated.prepared.minimized.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A152", "A216")
				.build();
		protein.flexibility.get("A177").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A180").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A181").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A184").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A185").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A219", "A220")
				.build();
		ligand.flexibility.get("A219").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();


		return confSpaces;
	}


	private ConfSpaces make2XXC10Res() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/2xxc.prepared.minimized.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B30", "B111")
				.build();
		protein.flexibility.get("B47").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("B58").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("B61").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("B64").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B4", "B29")
				.build();
		ligand.flexibility.get("B4").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();


		return confSpaces;
	}

	@Test
	public void generate2XXM10Res() {
	}

	private ConfSpaces make2XXM10Res() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/2xxm_prepped.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A146", "A218")
				.build();
		protein.flexibility.get("A177").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A180").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A181").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A184").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A185").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B4", "B113")
				.build();
		ligand.flexibility.get("B47").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B58").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B61").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B64").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();


		return confSpaces;
	}

	private ConfSpaces make2XXMSmaller() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/2xxm_prepped.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A146", "A218")
				.build();
		protein.flexibility.get("A177").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A179").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A180").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B4", "B113")
				.build();
		ligand.flexibility.get("B58").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B60").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B61").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B64").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();

		return confSpaces;
	}



	@Test
	public void test1GUA11MARKVsTraditional() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make1GUA11();
		final double epsilon = 0.999;
		final int numSequences = 6;


		String traditionalTime = "(Not run)";
		boolean runkstar = true;
		Stopwatch timer = new Stopwatch().start();
		if(runkstar) {
			TestBBKStar.Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
			timer.stop();
			traditionalTime = timer.getTime(2);
			timer.reset();
			timer.start();
		}
		runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
		String MARKStarTime = timer.getTime(2);
		timer.stop();

		//assert2RL0(results, numSequences);
		System.out.println("Traditional time: "+traditionalTime);
		System.out.println("MARK* time: "+MARKStarTime);
	}


	@Test
	public void timeMARKStarVsTraditional() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.68;
		final int numSequences = 10;
		Stopwatch timer = new Stopwatch().start();
		TestBBKStar.Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
		timer.stop();
		String traditionalTime = timer.getTime(2);
		timer.reset();
		timer.start();
		results = runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
		String MARKStarTime = timer.getTime(2);
		timer.stop();

		//assert2RL0(results, numSequences);
		System.out.println("Traditional time: "+traditionalTime);
		System.out.println("MARK* time: "+MARKStarTime);
	}

	@Test
	public void testMARKStarVsKStar() {
		int numFlex = 8;
		double epsilon = 0.68;
		compareMARKStarAndKStar(numFlex, epsilon);
	}

	@Test
	public void testMARKStarTinyEpsilon() {
		printMARKStarComputationStats(runMARKStar(5, 0.68).get(0));

	}

	@Test
	public void compareReducedMinimizationsVsNormal() {
		int numFlex = 14;
		double epsilon = 0.68;
		REUDCE_MINIMIZATIONS = false;
		System.out.println("Trying without reduced minimizations...");
		Stopwatch runTime = new Stopwatch().start();
		runMARKStar(numFlex, epsilon);
		String withoutRMTime= runTime.getTime(2);
		runTime.stop();
		runTime.reset();
		runTime.start();
		REUDCE_MINIMIZATIONS = true;
		System.out.println("Retrying with reduced minimizations...");
		runMARKStar(numFlex, epsilon);
		runTime.stop();
		System.out.println("Without Reduced Minimization time: "+withoutRMTime);
		System.out.println("Reduced minimization time: "+runTime.getTime(2));
	}

	@Test
    public void test1GUASmall() {
	    runMARKStar(3,0.99);
    }

	@Test
	public void test1GUASmallUpTo()
	{
		int maxNumFlex = 8;
		double epsilon = 0.68;
		for(int i = 1; i < maxNumFlex; i++)
			compareMARKStarAndKStar(i,epsilon);
	}

	@Test
	public void test4hem7resNanBug() {
		try {
			ConfSpaces confSpaces = loadFromCFS("examples/python.KStar/4hem_B_7res_1.131E+41.cfs");
			TestBBKStar.Results results = runBBKStar(confSpaces, 1, 0.999999999999, null, 2, true);
			for(KStar.ScoredSequence seq: results.sequences)
				System.out.println(seq);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private void compareMARKStarAndKStar(int numFlex, double epsilon) {
		Stopwatch runTime = new Stopwatch().start();
		String kstartime = "(not run)";
		List<KStar.ScoredSequence> kStarSeqs = null;
		boolean runkstar = true;
		if(runkstar) {
			kStarSeqs = runKStarComparison(numFlex, epsilon);
			runTime.stop();
			kstartime = runTime.getTime(2);
			runTime.reset();
			runTime.start();
		}
		List<MARKStar.ScoredSequence> markStarSeqs = runMARKStar(numFlex, epsilon);
		runTime.stop();
		for(MARKStar.ScoredSequence seq: markStarSeqs)
			printMARKStarComputationStats(seq);
		if(runkstar)
			for(KStar.ScoredSequence seq: kStarSeqs)
				printKStarComputationStats(seq);
		System.out.println("K* time: "+kstartime);
		System.out.println("MARK* time: "+runTime.getTime(2));
	}



	protected static void printKStarComputationStats(KStar.ScoredSequence result)
	{}

	protected static void printMARKStarComputationStats(MARKStar.ScoredSequence result)
	{}

	@Test
	public void KStarComparison() {
		List<KStar.ScoredSequence> results = runKStarComparison(5,0.68);
		for (int index = 0; index < results.size(); index++) {
			int totalConfsEnergied = results.get(index).score.complex.numConfs + results.get(index).score.protein.numConfs + results.get(index).score.ligand.numConfs;
			int totalConfsLooked = -1;
			System.out.println(String.format("score:%12e in [%12e,%12e], confs looked at:%4d, confs minimized:%4d",results.get(index).score.score, results.get(index).score.lowerBound,
					results.get(index).score.upperBound,totalConfsLooked,totalConfsEnergied));
		}
	}



	@Test
	public void testMARKStar() {
		List<MARKStar.ScoredSequence> results = runMARKStar(1, 0.01);
		for (int index = 0; index < results.size(); index++) {
			int totalConfsEnergied = results.get(index).score.complex.numConfs + results.get(index).score.protein.numConfs + results.get(index).score.ligand.numConfs;
			int totalConfsLooked = -1;
			System.out.println(String.format("score:%12e in [%12e,%12e], confs looked at:%4d, confs minimized:%4d",results.get(index).score.score, results.get(index).score.lowerBound,
					results.get(index).score.upperBound,totalConfsLooked,totalConfsEnergied));
		}
	}

	@Test
    public void test1GUASmallDEEP() {
	    ConfSpaces confSpaces = make1GUASmallDEEP(5);
	    runMARKStar(confSpaces, 0.99);
    }

	private static List<MARKStar.ScoredSequence> runMARKStar(int numFlex, double epsilon) {
		//ConfSpaces confSpaces = make1GUASmallCATS(numFlex);
		//ConfSpaces confSpaces = make1GUASmallDEEP(numFlex);
		ConfSpaces confSpaces = make1GUASmall(numFlex);
		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.setIsMinimizing(true)
				.build();
		// Define the rigid energy calculator
		EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
				.setIsMinimizing(false)
				.build();
		// how should we define energies of conformations?
		MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
			return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.setEnergyPartition(ENERGY_PARTITION)
					.build();
		};

		MARKStar.Settings settings = new MARKStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setEnergyMatrixCachePattern("*testmat.emat")
				.setShowPfuncProgress(true)
				.setParallelism(parallelism)
				.setReduceMinimizations(REUDCE_MINIMIZATIONS)
				.build();
		MARKStar run = new MARKStar(confSpaces.protein, confSpaces.ligand,
				confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
		return run.run();
	}

	public static class KstarResult {
		public KStar kstar;
		public List<KStar.ScoredSequence> scores;
	}

	public static List<KStar.ScoredSequence> runKStarComparison(int numFlex, double epsilon) {
		//ConfSpaces confSpaces = make1GUASmallCATS(numFlex);
		//ConfSpaces confSpaces = make1GUASmallDEEP(numFlex);
		ConfSpaces confSpaces = make1GUASmall(numFlex);
		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.build();
		// configure K*
		KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.setShowPfuncProgress(true)
				.build();
		KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
		for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

			// how should we define energies of conformations?
			info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, minimizingEcalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();

			// calc energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

			// how should confs be ordered and searched?
			info.confSearchFactory = (rcs) -> {
				ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
				return builder.build();
			};
		}

		// run K*
		KstarResult result = new KstarResult();
		result.kstar = kstar;
		result.scores = kstar.run();
		return result.scores;
	}

	public static List<KStar.ScoredSequence> runKStar(TestKStar.ConfSpaces confSpaces, double epsilon) {

		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.build();
		// configure K*
		KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.setShowPfuncProgress(true)
				.build();
		KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
		for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

			// how should we define energies of conformations?
			info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, minimizingEcalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();

			// calc energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.setCacheFile(new File(info.type+".kstar.emat"))
					.build()
					.calcEnergyMatrix();

			// how should confs be ordered and searched?
			info.confSearchFactory = (rcs) -> {
				ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
				return builder.build();
			};
		}

		// run K*
		KstarResult result = new KstarResult();
		result.kstar = kstar;
		result.scores = kstar.run();
		return result.scores;
	}

	private void runMARKStarWithLUTE(ConfSpaces confSpaces, double v) {
	}

	public static Result runMARKStar(TestKStar.ConfSpaces confSpaces, double epsilon){

		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.build();
		// Define the rigid energy calculator
		EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.setIsMinimizing(false)
				.build();
		// how should we define energies of conformations?
		MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
			return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();
		};

		Result result = new Result();

		MARKStar.Settings settings = new MARKStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setEnergyMatrixCachePattern("*.mark.emat")
				.setShowPfuncProgress(true)
				.setParallelism(parallelism)
				.setReduceMinimizations(REUDCE_MINIMIZATIONS)
				.build();

		result.markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
		result.markstar.precalcEmats();
		result.scores = result.markstar.run();
		return result;
	}


	public static ConfSpaces make1GUASmall(int numFlex) {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("1", "180")
				.build();
		int start = 21;
		for(int i = start; i < start+numFlex; i++) {
			protein.flexibility.get(i+"").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}



		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();

		return confSpaces;
	}
	public static ConfSpaces make1GUASmallDEEP(int numFlex) {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		ArrayList <String> bbflexlist = new ArrayList();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("1", "180")
				.build();
		int start = 21;
		for(int i = start; i < start+numFlex; i++) {
			protein.flexibility.get(i+"").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			bbflexlist.add(i+"");
		}



		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		bbflexlist.add("209");
		DEEPerSettings deepersettings = new DEEPerSettings(true, "test_deeper.pert", true, "None", false, 2.5,2.5, false, bbflexlist, "/1gua_adj.min.pdb", false, templateLib);
		DEEPerStrandFlex ligand_bbflex = new DEEPerStrandFlex(ligand, deepersettings);

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand, ligand_bbflex)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand, ligand_bbflex)
				.build();

		return confSpaces;
	}

	public static ConfSpaces make1GUA11() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("1", "180")
				.build();
		protein.flexibility.get("21").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("25").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("27").setLibraryRotamers(Strand.WildType, "HID").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("29").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("40").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("213").setLibraryRotamers(Strand.WildType, "HID", "HIE", "LYS", "ARG").addWildTypeRotamers().setContinuous();

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();

		return confSpaces;
	}

	public static void printSequence(Result result, int sequenceIndex){
		MARKStar.ScoredSequence scoredSequence =result.scores.get(sequenceIndex);
		String out = "Printing sequence "+sequenceIndex+": "+scoredSequence.sequence.toString(Sequence.Renderer.ResType)+"\n"+
				"Protein LB: "+String.format("%6.3e",scoredSequence.score.protein.values.qstar)+
				" Protein UB: "+String.format("%6.3e",scoredSequence.score.protein.values.pstar)+"\n"+
				"Ligand LB: "+String.format("%6.3e",scoredSequence.score.ligand.values.qstar)+
				" Ligand UB: "+String.format("%6.3e",scoredSequence.score.ligand.values.pstar)+"\n"+
				"Complex LB: "+String.format("%6.3e",scoredSequence.score.complex.values.qstar)+
				" Complex UB: "+String.format("%6.3e",scoredSequence.score.complex.values.pstar)+"\n"+
				"KStar Score: "+String.format("%6.3e",scoredSequence.score.complex.values.pstar.divide(scoredSequence.score.ligand.values.qstar.multiply(scoredSequence.score.protein.values.qstar), RoundingMode.HALF_UP));
		System.out.println(out);
	}


	public static void assertResult(PartitionFunction.Result result, Double qstar, double epsilon) {
		if (qstar != null) {
			double comparison = qstar*(1.0 - epsilon);
			if(comparison < 1e-5)
				comparison = 0;
			assertThat(result.status, is(PartitionFunction.Status.Estimated));
			assertThat(result.values.qstar.doubleValue(), greaterThanOrEqualTo(comparison));
			assertThat(result.values.getEffectiveEpsilon(), lessThanOrEqualTo(epsilon));
		} else {
			assertThat(result.status, is(not(PartitionFunction.Status.Estimated)));
		}
	}


}
