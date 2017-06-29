package edu.duke.cs.osprey.sparse;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.TermECalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;

public class ResidueInteractionGraph {
	
	ArrayList<Integer> vertices = new ArrayList<>();
	Map<Integer,Set<Integer>> adjacencyMatrix = new HashMap<>();
	Map<Integer, Integer> internalMap = new HashMap<>();
	Map<Integer, Integer> PDBIndexMap = new HashMap<>();
	double distanceCutoff = Double.MAX_VALUE; // distance cutoff, in angstroms
	double energyCutoff = 0; // energy cutoff, in kcal/mol
	double[][] distanceBounds;
	double[][] energyBounds;
	
	public ResidueInteractionGraph()
	{
		
	}
	
	public Set<Integer> getVertices()
	{
		return new HashSet<Integer>(vertices);
	}
	
	private int translateToInternalIndex(int index)
	{
		return internalMap.get(index);
	}
	
	private int translateToPDBIndex(int index)
	{
		return PDBIndexMap.get(index);
	}
	public double getEnergyBound(int v1, int v2)
	{
		return energyBounds[translateToInternalIndex(v1)][translateToInternalIndex(v2)];
	}
	
	public double getDistanceBound(int v1, int v2)
	{
		return distanceBounds[translateToInternalIndex(v1)][translateToInternalIndex(v2)];
	}
	
	public void updateEnergyBound(int v1, int v2, double energy)
	{
		int internalV1 = translateToInternalIndex(v1);
		int internalV2 = translateToInternalIndex(v2);
		energyBounds[internalV1][internalV2] = Math.max(energy, energyBounds[internalV1][internalV2]);
	}
	
	public void updateDistanceBound(int v1, int v2, double distance)
	{
		int internalV1 = translateToInternalIndex(v1);
		int internalV2 = translateToInternalIndex(v2);
		distanceBounds[internalV1][internalV2] = Math.min(distance, distanceBounds[internalV1][internalV2]);
	}
	
	public static ResidueInteractionGraph generateCompleteGraph(int numResidues)
	{
		ResidueInteractionGraph outputGraph = new ResidueInteractionGraph();
		for(int i = 0; i < numResidues; i++)
		{
			outputGraph.addVertex(i);
			for(int j = 0; j < i; j++)
			{
				outputGraph.addEdge(i,j);
			}
		}
		return outputGraph;
	}
	
	public static ResidueInteractionGraph generateCompleteGraph(SearchProblem problem)
	{
		ResidueInteractionGraph graph = new ResidueInteractionGraph();
		Set<Integer> residueIndexSet = createResidueIndexSet(problem);
		graph.setMutableResidues(residueIndexSet);
		
		return graph;
	}
	
	public static ResidueInteractionGraph generateGraph(
			ArrayList<Residue> residues, Molecule m,
			SearchProblem problem, EnergyFunction termE,
			double distanceCutoff, double energyCutoff)
	{
		ResidueInteractionGraph graph = new ResidueInteractionGraph();
		Set<Integer> residueIndexSet = createResidueIndexSet(problem);
		graph.setMutableResidues(residueIndexSet);
		graph.applyDistanceAndEnergyCutoff(distanceCutoff, energyCutoff, problem, termE);
		
		return graph;
	}
	
	public static ResidueInteractionGraph generateGraph(
			ArrayList<Residue> residues, Molecule m,
			SearchProblem problem, EnergyFunction termE)
	{
		ResidueInteractionGraph graph = new ResidueInteractionGraph();
		Set<Integer> residueIndexSet = createResidueIndexSet(problem);
		graph.setMutableResidues(residueIndexSet);
		
		return graph;
	}
	
	private static Set<Integer> createResidueIndexSet (SearchProblem problem) {
		Set<Integer> residueIndexSet = new HashSet<>();
		for(PositionConfSpace space: problem.confSpace.posFlex)
		{
			residueIndexSet.add(space.res.getPDBIndex());
		}
		return residueIndexSet;
	}

	public void addVertex(int vertex)
	{
		internalMap.put(vertex, vertices.size());
		PDBIndexMap.put(vertices.size(), vertex);
		vertices.add(vertex);
	}
	
	public void addEdge(int v1, int v2)
	{
		int min = Math.min(v1,v2);
		int max = Math.max(v1, v2);
		if(!adjacencyMatrix.containsKey(min))
			adjacencyMatrix.put(min, new HashSet<>());
		adjacencyMatrix.get(min).add(max);
	}
	
	public void pruneEdge(int v1, int v2)
	{
		int min = Math.min(v1,v2);
		int max = Math.max(v1, v2);
		if(!adjacencyMatrix.containsKey(min))
			adjacencyMatrix.put(min, new HashSet<>());
		adjacencyMatrix.get(min).remove(max);
	}
	
	public boolean connected(int source, int target)
	{
		int min = Math.min(source,target);
		int max = Math.max(source,target);
		return adjacencyMatrix.containsKey(min) 
				&& adjacencyMatrix.get(min).contains(max);
	}
	

	public void applyDistanceCutoff(double cutoff, SearchProblem problem, EnergyFunction termE)
	{
		distanceCutoff = cutoff;
		computeGraph(problem, termE);
	}
	
	public void applyEnergyCutoff(double cutoff, SearchProblem problem, EnergyFunction termE)
	{
		energyCutoff = cutoff;
		computeGraph(problem, termE);
	}
	
	public void applyDistanceAndEnergyCutoff(double dCutoff, double eCutoff,
			SearchProblem problem, EnergyFunction termE)
	{
		energyCutoff = eCutoff;
		distanceCutoff = dCutoff;
		computeGraph(problem, termE);
	}
	
	public void computeEdgeBounds (SearchProblem problem, EnergyFunction termE)
	{
		ConfSpace conformations = problem.confSpace;
		RCs RCSpace = new RCs(problem.pruneMat);
		int numResidues = vertices.size();
		double[][] pairwiseEnergyMaxBounds = new double[numResidues][numResidues];
		double[][] pairwiseEnergyMinBounds = new double[numResidues][numResidues];
		distanceBounds = new double[numResidues][numResidues];
		energyBounds = new double[numResidues][numResidues];
		
		for(int i = 0; i < pairwiseEnergyMaxBounds.length; i++)
			for(int j = 0; j < pairwiseEnergyMaxBounds[i].length; j++)
			{
				pairwiseEnergyMaxBounds[i][j] = Double.NEGATIVE_INFINITY;
				pairwiseEnergyMinBounds[i][j] = Double.POSITIVE_INFINITY;
				distanceBounds[i][j] = Double.POSITIVE_INFINITY;
			}
		
		for(int i =0; i < RCSpace.getNumPos(); i++)
		{
			Residue resi = conformations.posFlex.get(i).res;
			for(int j = 0; j < RCSpace.getNum(i); j++)
			{
				int unprunedRCi = RCSpace.get(i, j);
				RCTuple conformationTupleTest = new RCTuple(i,unprunedRCi);
				MoleculeModifierAndScorer mofTest = new MoleculeModifierAndScorer(termE,conformations,conformationTupleTest);

	            DoubleMatrix1D bestDOFValsTest;

	            if(mofTest.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
	                CCDMinimizer ccdMin = new CCDMinimizer(mofTest,true);
	                bestDOFValsTest = ccdMin.minimize().dofValues;
	            }
	            else//molecule is already in the right, rigid conformation
	            	bestDOFValsTest = DoubleFactory1D.dense.make(0);


	            double oneBodyEnergy = mofTest.getEnergyAndReset(bestDOFValsTest);

				for(int k = i+1; k < RCSpace.getNumPos(); k++)
				{
					Residue resj = conformations.posFlex.get(k).res;
					for(int l = 0; l < RCSpace.getNum(k); l++)
					{
						int unprunedRCk = RCSpace.get(k, l);
						RCTuple conformationTuple = new RCTuple(i,unprunedRCi,k,unprunedRCk);
						MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(termE,conformations,conformationTuple);

			            DoubleMatrix1D bestDOFVals;

			            if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
			                CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
			                bestDOFVals = ccdMin.minimize().dofValues;
			            }
			            else//molecule is already in the right, rigid conformation
			                bestDOFVals = DoubleFactory1D.dense.make(0);


			            double pairwiseEnergy = mof.getEnergyAndReset(bestDOFVals);
						if(pairwiseEnergy > 100)
						{
							//System.out.println("Clash?");
						}
			            double distance = resi.distanceTo(resj);
			            distanceBounds[i][k] = Math.min(distance, distanceBounds[i][k]);
			            updateDistanceBound(vertices.get(i),vertices.get(k),distance);
			            pairwiseEnergyMaxBounds[i][k] = Math.max(pairwiseEnergy, pairwiseEnergyMaxBounds[i][k]);
			            pairwiseEnergyMinBounds[i][k] = Math.min(pairwiseEnergy, pairwiseEnergyMaxBounds[i][k]);
//			            System.out.println("Energy of ("+i+"-"+j+","+k+"-"+l+"):"+pairwiseEnergy);
//			            System.out.println("Distance between ("+i+"-"+j+","+k+"-"+l+"):"+distance);
			            	
					}

				}
			}

		}
		
		for(int i = 0; i < vertices.size(); i++)
			for(int j = i+1; j < vertices.size(); j++)
			{
				updateEnergyBound(vertices.get(i),vertices.get(j), pairwiseEnergyMaxBounds[i][j] - pairwiseEnergyMinBounds[i][j]);
				updateEnergyBound(vertices.get(j),vertices.get(i), pairwiseEnergyMaxBounds[i][j] - pairwiseEnergyMinBounds[i][j]);
			}
	}
	
	public void computeEdgeBoundsOld (SearchProblem problem, EnergyFunction termE)
	{
		ConfSpace conformations = problem.confSpace;
		int numResidues = vertices.size();
		double[][] pairwiseEnergyMaxBounds = new double[numResidues][numResidues];
		double[][] pairwiseEnergyMinBounds = new double[numResidues][numResidues];
		distanceBounds = new double[numResidues][numResidues];
		energyBounds = new double[numResidues][numResidues];
		
		for(int i = 0; i < pairwiseEnergyMaxBounds.length; i++)
			for(int j = 0; j < pairwiseEnergyMaxBounds[i].length; j++)
			{
				pairwiseEnergyMaxBounds[i][j] = Double.NEGATIVE_INFINITY;
				pairwiseEnergyMinBounds[i][j] = Double.POSITIVE_INFINITY;
				distanceBounds[i][j] = Double.POSITIVE_INFINITY;
			}
		
		for(int i =0; i < vertices.size(); i++)
		{
			Residue resi = conformations.posFlex.get(i).res;
			for(int j = 0; j < conformations.posFlex.get(i).RCs.size(); j++)
			{
				RCTuple conformationTupleTest = new RCTuple(i,j);
				MoleculeModifierAndScorer mofTest = new MoleculeModifierAndScorer(termE,conformations,conformationTupleTest);

	            DoubleMatrix1D bestDOFValsTest;

	            if(mofTest.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
	                CCDMinimizer ccdMin = new CCDMinimizer(mofTest,true);
	                bestDOFValsTest = ccdMin.minimize().dofValues;
	            }
	            else//molecule is already in the right, rigid conformation
	            	bestDOFValsTest = DoubleFactory1D.dense.make(0);


	            double oneBodyEnergy = mofTest.getEnergyAndReset(bestDOFValsTest);

				for(int k = i+1; k < vertices.size(); k++)
				{
					Residue resj = conformations.posFlex.get(k).res;
					for(int l = 0; l < conformations.posFlex.get(k).RCs.size(); l++)
					{
						RCTuple conformationTuple = new RCTuple(i,j,k,l);
						MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(termE,conformations,conformationTuple);

			            DoubleMatrix1D bestDOFVals;

			            if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
			                CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
			                bestDOFVals = ccdMin.minimize().dofValues;
			            }
			            else//molecule is already in the right, rigid conformation
			                bestDOFVals = DoubleFactory1D.dense.make(0);


			            double pairwiseEnergy = mof.getEnergyAndReset(bestDOFVals);
						if(pairwiseEnergy > 100)
						{
							//System.out.println("Clash?");
						}
			            double distance = resi.distanceTo(resj);
			            distanceBounds[i][k] = Math.min(distance, distanceBounds[i][k]);
			            updateDistanceBound(vertices.get(i),vertices.get(k),distance);
			            pairwiseEnergyMaxBounds[i][k] = Math.max(pairwiseEnergy, pairwiseEnergyMaxBounds[i][k]);
			            pairwiseEnergyMinBounds[i][k] = Math.min(pairwiseEnergy, pairwiseEnergyMaxBounds[i][k]);
//			            System.out.println("Energy of ("+i+"-"+j+","+k+"-"+l+"):"+pairwiseEnergy);
//			            System.out.println("Distance between ("+i+"-"+j+","+k+"-"+l+"):"+distance);
			            	
					}

				}
			}

		}
		
		for(int i = 0; i < vertices.size(); i++)
			for(int j = i+1; j < vertices.size(); j++)
			{
				updateEnergyBound(vertices.get(i),vertices.get(j), pairwiseEnergyMaxBounds[i][j] - pairwiseEnergyMinBounds[i][j]);
				updateEnergyBound(vertices.get(j),vertices.get(i), pairwiseEnergyMaxBounds[i][j] - pairwiseEnergyMinBounds[i][j]);
			}
	}
	
	public void computeGraph(SearchProblem problem, EnergyFunction termE)
	{
		
		computeEdgeBounds(problem, termE);
		
		int edgesPruned = 0;
		for(int i = 0; i < vertices.size(); i++)
			for(int j = i+1; j < vertices.size(); j++)
			{
				double minDistance = distanceBounds[i][j];
				if(energyBounds[i][j] < energyCutoff || minDistance > distanceCutoff)
				{
					int iPDBIndex = translateToPDBIndex(i);
					int jPDBIndex = translateToPDBIndex(j);
					System.out.println("Pruning edge ("+iPDBIndex+","+jPDBIndex+"), energy "+energyBounds[i][j]+", distance "+minDistance);
					pruneEdge(iPDBIndex,jPDBIndex);
					edgesPruned++;
				}
			}
		
		System.out.println("Edges pruned: "+edgesPruned);
				
	}
	
	
	public void setMutableResidues(Set<Integer> residues)
	{
		vertices.clear();
		for(Integer i : residues)
			addVertex(i);
		createCompleteGraph();
	}
	
	private void createCompleteGraph()
	{
		for(Integer i : vertices)
		{
			for(Integer j : vertices)
			{
				if(i!=j)
					addEdge(i,j);
			}
		}
	}

	public void writeGraph (String outputFileName) {
		// TODO Auto-generated method stub
		try
		{
			FileOutputStream fileOutputStream = new FileOutputStream(outputFileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream(
					fileOutputStream);
			PrintStream printStream = new PrintStream(bufferedOutputStream);


			for(Integer i : adjacencyMatrix.keySet())
			{
				for(Integer j : adjacencyMatrix.get(i))
				{
					printStream.println(i+","+j);
				}
				if(adjacencyMatrix.isEmpty())
					printStream.println(i+","+i);
			}
			printStream.close();
		}        
		catch (IOException e) {
			System.out.println("ERROR: An io exception occurred while writing file "+outputFileName);
			System.exit(0);
		}
		catch ( Exception e ){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing file");
			System.exit(0);
		}
	}
	
	public void printStatistics()
	{
		System.out.println("==================== Edge Statistics =========================");
		for(int i = 0; i < vertices.size(); i++)
			for(int j = i+1; j < vertices.size(); j++)
			{
				int iPDBIndex = translateToPDBIndex(i);
				int jPDBIndex = translateToPDBIndex(j);
					System.out.println("Edge ("+iPDBIndex+","+jPDBIndex+"): energy "
							+energyBounds[i][j]+", distance "+distanceBounds[i][j]);
				
			}
	}

}
