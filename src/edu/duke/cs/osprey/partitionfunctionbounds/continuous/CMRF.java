/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import Jama.Matrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.ToDoubleFunction;

/**
 * Continuous-label Markov Random Field
 * 
 * Includes an implementation of TRBP/SCMF upper/lower bounds on the log partition function 
 * 
 * @author aditya
 */
public class CMRF {

	public CMRFNode[] nodes;
	public int numNodes;
	public int[][] adjacencyMatrix;
	public double[] nodeWeights;

	public double[][] edgeProbs;
	public CMRFEdge[][] edges;

	double threshold = 0.0001;
	double constRT = PoissonBoltzmannEnergy.constRT;
	int maxIters = 1000000;
	double lambda = 0.7;

	boolean nodesAdded = false;
	boolean ranSCMF = false;

	/**
	 * Sets up an empty cMRF -- population is done later 
	 * @param numNodes 
	 */
	public CMRF(int numNodes) { 
		this.numNodes = numNodes;
		nodes = new CMRFNode[numNodes];
		edges = new CMRFEdge[numNodes][numNodes];
	}
	
	/**
	 * Runs a quick toy example using both SCMF and TRBP
	 * Actual logZ value for this example is roughly 171.526 
	 * @param size
	 * @param kernelMult
	 * @param iter
	 * @return
	 */
	public static double[] runToyCMRF (double size, double kernelMult, int iter) { 

		//Main.main(args);

		System.out.println("CMRF main");
				
		double[][] b1 = new double[1][2]; b1[0][0] = 0; b1[0][1] = size;
		double[] lb1 = new double[1]; lb1[0] = 0; 
		double[] ub1 = new double[1]; ub1[0] = size;
		Kernel k1 = new KernelGaussian( b1 , kernelMult*size );
		CMRFNodeDomain nd1 = new CMRFNodeDomain(lb1, ub1, k1, (point)->(-Math.pow(point[0]-5,2)));

		double[][] b2 = new double[2][2];
		b2[0][0] = 0; b2[0][1] = size;
		b2[1][0] = 0; b2[1][1] = size;
		double[] lb2 = {0, 0};
		double[] ub2 = {size, size};
		Kernel k2 = new KernelGaussian( b2, kernelMult*size);
		CMRFNodeDomain nd2 = new CMRFNodeDomain(lb2, ub2, k2, (point)->(-1*point[0]));

		HashMap<Integer, CMRFNodeDomain[]> h = new HashMap<>();
		h.put(0, new CMRFNodeDomain[]{nd1});
		h.put(1, new CMRFNodeDomain[]{nd2});

		ToDoubleFunction<double[]>f = (point)->(-1*point[2]);
		HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map1 = new HashMap<>();
		map1.put(nd2, f);
		HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map2 = new HashMap<>();
		map2.put(nd1, map1);
		HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map3 = new HashMap<>();
		map3.put(1, map2);
		HashMap<Integer, HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>>> map4 = new HashMap<>();
		map4.put(0, map3);
		
		HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map5 = new HashMap<>();
		map5.put(nd1, f);
		HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map6 = new HashMap<>();
		map6.put(nd2, map5);
		HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map7 = new HashMap<>();
		map7.put(0, map6);
		map4.put(1, map7);


		CMRF c = new CMRF(2);
		c.constRT = 1.0;
		c.addNodes(h, map4);
		
		SCMF s = new SCMF(c);
		TRBP t = new TRBP(c);
		
		double logZLB = s.runSCMF();
		double logZUB = t.runTRBP(0); // no iterations of LBP
		
		double[] ret = new double[2];
		ret[0] = logZLB; ret[1] = logZUB;
		return ret;
	}


	/**
	 * Adds a list of nodes to the cMRF -- each node is given a set of domains and edges are set up
	 * The eFuncMap is a map that gives you the pairwise energy functions. If N is the set of nodes and D_i is the set
	 * of domains for node i, then eFuncMap is a function from N x N -> D_i x D_j -> function
	 *
	 * @param domains
	 * @param eFuncMap
	 */
	public void addNodes(
			HashMap<Integer, CMRFNodeDomain[]> domains,
			// mapping is node 1 --> node 2 --> domain 1 --> domain 2 --> eFunc
			HashMap<Integer,
			HashMap<Integer,
			HashMap<CMRFNodeDomain,
			HashMap<CMRFNodeDomain,
			ToDoubleFunction<double[]>>>>> eFuncMap) {

		System.out.println("Adding nodes...");

		// make the nodes
		for (int i=0; i<numNodes; i++) {
			nodes[i] = new CMRFNode(domains.get(i));
		}

		// initialize outmessage maps
		for (int i=0; i<numNodes; i++) {
			System.out.print(i+"/"+numNodes+" ");
			CMRFNode node = nodes[i];
			node.outMessages = new HashMap<>();
			for (int j=0; j<numNodes; j++) {
				if (i==j) { continue; }
				node.outMessages.put(nodes[j], new HashMap<>());
			}
			System.out.println();
		}

		System.out.println("Adding edges...");
		for (int i=0; i<nodes.length; i++) {
			for (int j=0; j<nodes.length; j++) {
				if (i==j) { continue; }
				System.out.print(i+"-"+j+" ");
				edges[i][j] = new CMRFEdge(
						nodes[i],
						nodes[j],
						eFuncMap.get(i).get(j));
			}
			System.out.println();
		}

		nodesAdded = true;

		printNaivePDFs();
	}

	/** 
	 * Returns the sum of an array of RKHSFunctions evaluated at a given point
	 * @param point
	 * @param funcs
	 * @return 
	 */
	public double sumOverMessages(double[] point, RKHSFunction[] funcs) { 
		double res = 0.0;
		for (RKHSFunction func : funcs) {
			res += func.eval(point);
		}
		return res;
	}
	
	public double sumOverMessages(double[] point, ArrayList<ToDoubleFunction<double[]>> funcs) { 
		double res = 0.0;
		res = funcs.stream().map((func) -> func.applyAsDouble(point)).reduce(res, (accumulator, _item) -> accumulator + _item);
		return res;
	}

	// and introduce a variant that lets us evaluate LCs of RKHS functions
	public double sumOverMessages(double[] point, RKHSFunction[] funcs, double[] coeffs) { 
		if (funcs.length != coeffs.length) { 
			throw new RuntimeException("Coefficients don't match functions.");
		}
		double res = 0.0;
		for (int i=0; i<funcs.length; i++) { 
			res += coeffs[i] * funcs[i].eval(point);
		}
		return res;
	}

	/**
	 * Splits an array at the nth slot
	 * @param arr
	 * @param n
	 * @return 
	 */
	public static ArrayList<double[]> splitArray(double[] arr, int n) { 
		ArrayList<double[]> arrs = new ArrayList<>();
		arrs.add(Arrays.copyOfRange(arr, 0, n));
		arrs.add(Arrays.copyOfRange(arr, n, arr.length));
		return arrs;
	}

	/**
	 * Concatenates two arrays
	 * I feel really stupid having actually written this method 
	 * It's literally just a wrapper to save me some characters 
	 * @param arr1
	 * @param arr2
	 * @return 
	 */
	public static double[] concatArrays(double[] arr1, double[] arr2) { 
		return CMRFEdgeDomain.concatArrays(arr1, arr2);
	}

	/** 
	 * Gets the product of several RKHSFunctions, each raised to a possibly distinct power 
	 * @param funcs
	 * @param powers
	 * @param point
	 * @return 
	 */
	public double getProdOfFuncPowers(RKHSFunction[] funcs, double[] powers, double[] point) { 
		double result = 1.0;
		for (int i=0; i<funcs.length; i++) { 
			result *= Math.pow(funcs[i].eval(point), powers[i]);
		}
		return result;
	}

	/**
	 * Gets the index of an object in an array
	 * Honestly what the hell was I thinking I don't even know anymore 
	 * @param o
	 * @param arr
	 * @return 
	 */
	public int getIndexInArray(Object o, Object[] arr) { 
		for (int i=0; i<arr.length; i++) { 
			if (arr[i].equals(o)) { 
				return i;
			}
		}
		return -1;
	}


	/**
	 * Returns the PDF for a particular CMRF node's domain 
	 * @param v
	 * @param d
	 * @return 
	 */
	public RKHSFunction getPDF(CMRFNode v, CMRFNodeDomain d) { 
		// get the pdf for the node domain
		//   p(x_s) = \sum_{neibhors} {m_{t->s}{x_s) * \mu_{ts}}
		int recNodeIndex = this.getIndexInArray(v, nodes);
		RKHSFunction[] messages = new RKHSFunction[nodes.length-1];
		double[] vEdgeProbs = new double[nodes.length-1];

		for (CMRFNode n : nodes) {
			if (n.equals(v)) { continue; }
			int sendNodeIndex = this.getIndexInArray(n, nodes);
			messages[sendNodeIndex] = n.outMessages.get(v).get(d);
			vEdgeProbs[sendNodeIndex] = this.edgeProbs[sendNodeIndex][recNodeIndex];
		}
		RKHSFunction probabilityFunc = new RKHSFunction(
				d.k,
				d.domainLB,
				d.domainUB,
				(point) -> (
						this.sumOverMessages(point, messages, vEdgeProbs)));
		return probabilityFunc;
	}

	/**
	 * Returns a function that returns the product of the two PDFs -- note this is NOT the inter-rotamer probability
	 * density but the product of the intra-rotamer probability densities! 
	 * @param n1
	 * @param d1
	 * @param n2
	 * @param d2
	 * @return 
	 */
	public RKHSFunction getProductOverCrossPDF(CMRFNode n1, CMRFNodeDomain d1, CMRFNode n2, CMRFNodeDomain d2) { 
		RKHSFunction pdf1 = this.getPDF(n1, d1);
		RKHSFunction pdf2 = this.getPDF(n2, d2);

		int sendNodeInd = this.getIndexInArray(n1, nodes);
		int recNodeInd = this.getIndexInArray(n2, nodes);
		Kernel prodK = this.edges[sendNodeInd][recNodeInd].getEdgeDomain(d1, d2).resAllK;
		return RKHSFunction.getCartesianProductFunction(pdf1, pdf2, prodK);
	}

	public void printNaivePDFs() { 
		System.out.print("Printing naive pdfs... ");
		if (!nodesAdded) { return; } // do nothing if we don't have an MRF

		for (int i=0; i<nodes.length; i++) {
			CMRFNode node = nodes[i];
			for (int j=0; j<node.domains.length; j++) {

				// print unitary pdfs
				CMRFNodeDomain domain = node.domains[j];
				try {
					String filename = "cmrf_u_"+i+"-"+j+".dat";
					PrintWriter writer = new PrintWriter(filename, "UTF-8");
					Matrix m = domain.probabilityRKHS.dumpPoints();
					m.print(writer, 3, 5);
					writer.flush();
				} catch(FileNotFoundException | UnsupportedEncodingException e) {
					throw new RuntimeException("PrintWriting failed for "+
							"node " + i +", domain " + j +".\n" +e.getMessage());
				}

				// print pairwise pdfs
				for (int k=0; k<nodes.length; k++) {
					if (k==i) { continue; } 
					CMRFEdge edge = edges[i][k];
					for (int l=0; l<edge.domainLinks.length; l++) {
						CMRFEdgeDomain edgeDomain = edge.domainLinks[l];
						int x1 = getIndexInArray(edgeDomain.resOneDomain, node.domains);
						int x2 = getIndexInArray(edgeDomain.resTwoDomain, nodes[k].domains);

						try {
							String filename = "cmrf_p_"+i+k+"_"+x1+x2+".dat";
							PrintWriter writer = new PrintWriter(filename, "UTF-8");
							Matrix m = edgeDomain.pFuncRKHS.dumpPoints();
							m.print(writer, 3, 5);
							writer.flush();
						} catch (FileNotFoundException | UnsupportedEncodingException e) {
							throw new RuntimeException("PrintWriting failed for "+
									"edge " + i+"-"+k +", domain " + x1+"-"+x2 +".\n" +e.getMessage());
						}
					}
				}
			}
		}	
		System.out.println("done.");
	}



}
