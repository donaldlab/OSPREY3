package edu.duke.cs.osprey.tests;

import static org.junit.Assert.*;
import java.math.BigInteger;
import org.junit.Test;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.sparse.BranchDecomposedProblem;
import edu.duke.cs.osprey.sparse.BranchTree;
import edu.duke.cs.osprey.sparse.ConformationProcessor;
import edu.duke.cs.osprey.sparse.ResidueInteractionGraph;
import edu.duke.cs.osprey.sparse.SparseKStarScoreEvaluator;
import edu.duke.cs.osprey.sparse.Subproblem;
import edu.duke.cs.osprey.sparse.TreeEdge;
import edu.duke.cs.osprey.sparse.TreeNode;
import edu.duke.cs.osprey.tools.ResidueIndexMap;
import edu.duke.cs.osprey.tools.branchdecomposition.BranchDecomposition;
import junit.framework.TestCase;

public class TestSparseAlgorithms  extends TestCase {

	ConfigFileParser cfp;
	String PDBFileLocation = "test/4NPD/4NPD.pdb";
	SearchProblem searchSpace;
	RCs fullRCSpace;
	ResidueIndexMap resMap;
	private static boolean bdGenerated = false;
	protected void setUp () throws Exception {
		super.setUp();


		String[] testArgs = new String[]
				{"-c", "test/1CC8Sparse/KStar.cfg", "Dummy command", 
						"test/1CC8Sparse/DEESparse.cfg", "test/1CC8Sparse/SystemSparse.cfg"};
		cfp = new ConfigFileParser(testArgs);//args 1, 3+ are configuration files
		cfp.loadData();
		searchSpace = cfp.getSearchProblem();
        double Ew = cfp.getParams().getDouble("Ew");
		searchSpace.loadEnergyMatrix();
		resMap = ResidueIndexMap.createResidueIndexMap(searchSpace.confSpace);
		double I0 = 0;
		
        boolean doIMinDEE = cfp.getParams().getBool("imindee");
        if(doIMinDEE){
            I0 = cfp.getParams().getDouble("Ival");
        }
		double pruningInterval = Ew + I0;
		
        double stericThresh = cfp.getParams().getDouble("StericThresh");
		
        PruningControl pruningControl = new PruningControl(
                searchSpace,
                0, // pruning interval, set by initPruning()
                cfp.getParams().getBool("TYPEDEP"), 
                cfp.getParams().getDouble("BOUNDSTHRESH"),
                cfp.getParams().getInt("ALGOPTION"), 
                cfp.getParams().getBool("USEFLAGS"),
                cfp.getParams().getBool("USETRIPLES"),
                false,
                false, // useEPIC, set by initPruning()
                false, // useTupExp, set by initPruning()
                stericThresh
            );

        System.out.println("Conformations before pruning: "+searchSpace.confSpace.getNumConformations());
        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        if(searchSpace.competitorPruneMat == null){
            System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
            initPruning(searchSpace, pruningControl, 0, false, false);
            pruningControl.setOnlyGoldstein(true);
            pruningControl.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
            System.out.println("COMPETITOR PRUNING DONE");
        }
        
        
        //Next, do DEE, which will fill in the pruning matrix
        initPruning(searchSpace, pruningControl, pruningInterval, false, false);
        pruningControl.prune();//pass in DEE options, and run the specified types of DEE    
        fullRCSpace = new RCs(searchSpace.pruneMat);
        System.out.println("Conformations after pruning: "+fullRCSpace.unprunedConfsFromRCs());
        fullRCSpace.breakDownRCSpace();
		
		
//		if(!bdGenerated)
//		{
//			String runName = cfp.getParams().getValue("runName");
//			String graphFileName = "test/1CC8Sparse/"+runName;
//			String bdFileName = "test/1CC8Sparse/"+runName+"_bd";
//			
//			EnergyFunction efunction = searchSpace.fullConfE;
//			ConfSpace conformationSpace = searchSpace.confSpace;
//			ResidueInteractionGraph graph = ResidueInteractionGraph.generateCompleteGraph(searchSpace);
//			graph.computeEdgeBounds(searchSpace, efunction);
//			graph.printStatistics();
//			graph.applyEnergyCutoff(0.2, searchSpace, efunction);
//			graph.writeGraph(graphFileName);
//
//			String[] args = new String[]{graphFileName, bdFileName};
//			long startBD = System.currentTimeMillis();
//			BranchDecomposition.main(args);
//			long endBD = System.currentTimeMillis();
//			long BDTime = endBD - startBD;
//
//			System.out.println("Branch Decomposition generation time: "+BDTime);
//			long start = System.currentTimeMillis();
//			System.out.println("Branch Decomposition generated. Calculating GMEC...");
//
//
//
//			long end = System.currentTimeMillis();
//			long time = end - start;
//			System.out.println("Total time BD generation time taken in ms: "+time);
//			bdGenerated = true;
//		}

	}

	protected void tearDown () throws Exception {
		super.tearDown();
	}


	/***
	 * This test determines if a constrained conformation space is working correctly. It will fail the moment
	 * confSpaces have the wrong size.
	 */
	@Test
	public void testConstrainConfSpace()
	{
		RCTuple initialConf = new RCTuple();
		for(int i = 0; i < searchSpace.confSpace.numPos/2; i++)
		{
			int RCAssigned = 0;
			initialConf = initialConf.addRC(i, RCAssigned);
		}
		RCs localConfSpace = fullRCSpace.returnSubspace(initialConf);
		if(localConfSpace.getNumPos() != fullRCSpace.getNumPos())
		{
			System.err.println("Test not designed to handle conformation spaces"
					+" of different sizes. (It should. Fix this.)");
			System.exit(-1);
		}

		// Make sure the ConfSpace is properly constrained.
		for(int tupleIndex = 0; tupleIndex < initialConf.size(); tupleIndex++)
		{
			int RCPosition = initialConf.pos.get(tupleIndex);
			int RCConf = initialConf.RCs.get(tupleIndex);

			int numRCs = localConfSpace.getNum(RCPosition);
			if(localConfSpace.getNum(RCPosition) != 1)
			{
				System.err.println("ERROR: Residue "+RCPosition
						+"hasn't been constrained. It still has "
						+localConfSpace.getNum(RCPosition)
						+" RCs allowed.");
			}
			assert(localConfSpace.getNum(RCPosition) == 1);
			int RCIndex = localConfSpace.get(RCPosition, 0);
			assert(RCIndex == RCConf);
			if(RCIndex != RCConf)
			{
				System.err.println("ERROR: Residue "+RCPosition
						+" doesn't match the constrained space: "
						+RCIndex
						+" isn't "+RCConf);
				System.exit(-1);
			}
			System.out.println("Residue "+RCPosition
					+" matches the constrained space: "
					+RCIndex
					+" == "+RCConf);
		}
		System.out.println("Test complete.");
	}
	
	@Test
	public void testSubproblemExhaustiveEnumeration()
	{

		String runName = cfp.getParams().getValue("runName");
		SearchProblem problem = cfp.getSearchProblem();
		EnergyFunction efunction = problem.fullConfE;
		ConfSpace conformationSpace = problem.confSpace;

		String[] args = new String[]{"test/1CC8Sparse/"+runName, "test/1CC8Sparse/"+runName+"_bd"};
		long startBD = System.currentTimeMillis();
		BranchDecomposition.main(args);
		long endBD = System.currentTimeMillis();
		long BDTime = endBD - startBD;

		System.out.println("Branch Decomposition generation time: "+BDTime);
		long start = System.currentTimeMillis();
		System.out.println("Branch Decomposition generated. Calculating GMEC...");

		String bdFile = "test/1CC8Sparse/"+runName+"_bd";

		BranchTree tree = new BranchTree(bdFile, problem);
		TreeNode root = tree.getRoot();
		TreeEdge rootEdge = root.getCofEdge();

		RCTuple initialConf = new RCTuple();
		for(int i = 0; i < problem.confSpace.numPos/2; i++)
		{
			int RCAssigned = (int)(Math.random()*problem.confSpace.posFlex.get(i).RCs.size());
			initialConf = initialConf.addRC(i, RCAssigned);
		}
		
		Subproblem sparseProblem = new TestSubproblem(new RCs(searchSpace.pruneMat), rootEdge, resMap, initialConf);
		ConformationCounter counter = new ConformationCounter();
		sparseProblem.addConformationProcessor(counter);
		sparseProblem.preprocess();
		BigInteger totalConfs = fullRCSpace.unprunedConfsFromRCs();
		BigInteger subproblemConfs = sparseProblem.getTotalLocalConformations();
		if(!counter.numConfs.equals(subproblemConfs))
		{
			System.err.println("Conformations not processed in subproblem: processed "+counter.numConfs+", expected "+subproblemConfs);
		}
		assert(counter.numConfs.equals(subproblemConfs));
		System.out.println("Num confs processed: "+counter.numConfs);
		System.out.println("Num subproblem confs possible: "+subproblemConfs);
		System.out.println("Num confs possible: "+totalConfs);
	}
	
	@Test
	public void testFullTreeExhaustiveEnumeration()
	{

		String runName = cfp.getParams().getValue("runName");
		SearchProblem problem = cfp.getSearchProblem();
		ConfSpace conformationSpace = problem.confSpace;


		String bdFile = "test/1CC8Sparse/"+runName+"_bd";

		BranchTree tree = new BranchTree(bdFile, problem);
		TreeEdge rootEdge = tree.getRootEdge();
		rootEdge.compactTree();
		rootEdge.printTreeMol("");
		

		RCTuple initialConf = new RCTuple();
		for(int i = 0; i < problem.confSpace.numPos/2; i++)
		{
			int RCAssigned = (int)(Math.random()*problem.confSpace.posFlex.get(i).RCs.size());
			initialConf = initialConf.addRC(i, RCAssigned);
		}
		Subproblem sparseProblem = new Subproblem(new RCs(searchSpace.pruneMat), rootEdge, resMap, initialConf);
		ConformationCounter counter = new ConformationCounter();
		sparseProblem.addConformationProcessor(counter);
		sparseProblem.preprocess();
		BigInteger totalConfs = conformationSpace.getNumConformations();
		BigInteger subproblemConfs = sparseProblem.getTotalLocalConformations();
		if(!counter.numConfs.equals(subproblemConfs))
		{
			System.err.println("Conformations not processed in subproblem: processed "+counter.numConfs+", expected "+subproblemConfs);
		}
		//assert(counter.numConfs.equals(subproblemConfs));
		System.out.println("Num confs processed: "+counter.numConfs);
		System.out.println("Num subproblem confs possible: "+subproblemConfs);
		System.out.println("Num confs possible: "+totalConfs);
	}

	@Test
	public void testComputeKStarScore()
	{

	}


	@Test
	public void testEnumerateConformations() throws Exception
	{

	}

	@Test
	public void testEnumerateConformationsByKStarScore() throws Exception
	{

	}
	
	private class ConformationCounter implements ConformationProcessor
	{
		BigInteger numConfs = BigInteger.ZERO;
		@Override
		public void processConformation (RCTuple conformation) {
			numConfs = numConfs.add(BigInteger.ONE);
		}
		
	}
	
	
    private void initPruning(SearchProblem searchSpace, PruningControl pruningControl,
    		double pruningInterval, boolean useEPIC, boolean useTupExp) {
        
        // init the pruning matrix if needed
        if(searchSpace.pruneMat == null || searchSpace.pruneMat.getPruningInterval() < pruningInterval) {
            searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        }
        
        // configure the pruner
        pruningControl.setOnlyGoldstein(false);
        pruningControl.setPruningInterval(pruningInterval);
        pruningControl.setUseEPIC(useEPIC);
        pruningControl.setUseTupExp(useTupExp);
    }
}
