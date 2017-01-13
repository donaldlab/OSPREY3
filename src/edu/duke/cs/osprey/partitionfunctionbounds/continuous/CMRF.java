/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import Jama.Matrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
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
    
    double logZ = Double.POSITIVE_INFINITY;
    double threshold = 1e-6;
    double constRT = PoissonBoltzmannEnergy.constRT;
    
    boolean nodesAdded = false;
    
    public CMRF(int numNodes) { 
	this.numNodes = numNodes;
	nodes = new CMRFNode[numNodes];
    }
    
    /**
     * Adds a list of nodes to the cMRF -- each node is given a set of domains and associated energy functions
     * @param domains
     * @param edgeEnergyFuncMap 
     */
    public void addNodes(
	    Map<Integer, CMRFNodeDomain[]> domains,
	    Map<Integer, Map<Integer, ToDoubleFunction<double[]>>> edgeEnergyFuncMap) { 
	for (int i=0; i<numNodes; i++) { 
	    nodes[i] = new CMRFNode(domains.get(i));
	}
	
	for (Integer i : edgeEnergyFuncMap.keySet()) { 
	    for (Integer j : edgeEnergyFuncMap.get(i).keySet()) { 
		edges[i][j] = new CMRFEdge(nodes[i], nodes[j], edgeEnergyFuncMap.get(i).get(j));
	    }
	}
	
	nodesAdded = true;
    }
    
    /**
     * Initializes edge probabilities -- this is shamelessly stolen from Hunter's code 
     */
    public void initializeEdgeProbsTRBP() { 
	if (! nodesAdded) { 
	    throw new RuntimeException("Can't initialize edge probabilities for the"
		    + "continuous-label MRF unless nodes have already been added. ");
	}
	
	Matrix adj = new Matrix(numNodes, numNodes, 1.0);
	Matrix deg = Jama.Matrix.identity(numNodes, numNodes).times(numNodes);
	Matrix laplacian = deg.minus(adj);
	Matrix lapMinOne = laplacian.minus(new Matrix(numNodes, numNodes, 1.0));
	Matrix invLap = lapMinOne.inverse();
	
	this.edgeProbs = new double[numNodes][numNodes];
	for (int i=0; i<edgeProbs.length; i++) { 
	    for (int j=0; j<edgeProbs[i].length; j++) { 
		edgeProbs[i][j] = 
			adj.get(i, j) * (invLap.get(i, i) + invLap.get(j, j) - 2*invLap.get(i, j));
	    }
	}
    }

    /**
     * Initializes TRBP update messages as uniform distributions 
     */
    public void initializeMessagesTRBP() { 
	for (int i=0; i<numNodes; i++) { 
	    CMRFNode sender = nodes[i];
	    for (int j=0; j<numNodes; j++) { 
		CMRFNode receiver = nodes[j];
		if (i == j) { continue; } // don't want to sender to be the receiver
		for (CMRFNodeDomain d : receiver.domains) {  // message is just a uniform distribution
		    RKHSFunction message = new RKHSFunction(
			    d.k,
			    d.domainLB,
			    d.domainUB,
			    (x)->(1.0/d.volume));
		    if (sender.outMessages.get(receiver) == null) {
			sender.outMessages.put(
				receiver, 
				new HashMap<CMRFNodeDomain, RKHSFunction>());
		    }
		    sender.outMessages.get(receiver).put(d, message);
		}
	    }
	}	
    }
    
    public void updateMessagesMeanField() { 
	
    }
    
    /**
     * TRBP message update procedure 
     */
    public void updateMessagesTRBP() {
	// first, construct the messages; then, update them all at once 

	// we'll store all the messages in a hashmap buffer
	// map from sender --> receiver --> receiver domain --> rkhsFunction
	HashMap<CMRFNode, HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>> messageMaps = 
		new HashMap<CMRFNode, HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>>();
	for (CMRFNode sender : nodes) { 
	    HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>> senderMap = 
		    new HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>();
	    for (CMRFNode receiver : nodes) {  
		if (sender.equals(receiver)) { continue; }
		HashMap<CMRFNodeDomain, RKHSFunction> domainMap = new HashMap<CMRFNodeDomain, RKHSFunction>();
		for (CMRFNodeDomain domain : receiver.domains) { 
		    domainMap.put(
			    domain,
			    new RKHSFunction(
				    domain.k,
				    domain.domainLB,
				    domain.domainUB,
				    (point)->(1.0/domain.volume)));
		}
		senderMap.put(receiver, domainMap);
	    }
	    messageMaps.put(sender, senderMap);
	}
	
	for (int i=0; i<numNodes; i++) {  // i is the sender
	    for (int j=0; j<numNodes; j++) { // j is the receiver
		if (i==j) { continue; }
		CMRFNode sender = nodes[i];
		CMRFNode receiver = nodes[j];
		final double edgeProb = edgeProbs[i][j];
		
		// each of the update functions is tied to a domain in the receiver
		// we'll collapse all of this when we're done constructing them
		RKHSFunction[] receiverDomainFuncs = new RKHSFunction[receiver.domains.length];
		
		for (CMRFNodeDomain recDom : receiver.domains) {
		    RKHSFunction[] senderFuncs = new RKHSFunction[sender.domains.length];
		    for (CMRFNodeDomain d : sender.domains) {
			// denominator of quotient
			RKHSFunction denomFunc = new RKHSFunction(
				d.k,
				d.domainLB,
				d.domainUB,
				(point) ->
					(Math.pow(
						sender.outMessages.get(receiver).get(d).eval(point),
						1-edgeProb)));
			
			// numerator of quotient
			ArrayList<CMRFNode> parents = new ArrayList<CMRFNode>();
			for (CMRFNode n : nodes) {
			    if (!n.equals(sender) && !n.equals(receiver)) {
				parents.add(n);
			    }
			}
			double[] parentProbs = new double[parents.size()];
			RKHSFunction[] parentMsgs = new RKHSFunction[parents.size()];
			for (int k=0; k<parentMsgs.length; k++) {
			    CMRFNode parent = parents.get(k);
			    parentMsgs[k] = parent.outMessages.get(sender).get(d);
			    parentProbs[k] = edgeProbs[i][this.getIndexInArray(parent, nodes)];
			}
			RKHSFunction numFunc = new RKHSFunction(
				d.k,
				d.domainLB,
				d.domainUB,
				(point) -> (getProdOfFuncPowers(parentMsgs, parentProbs, point)));
			
			// exponential component of update
			RKHSFunction updateFunc = new RKHSFunction(
				recDom.k,
				recDom.domainLB,
				recDom.domainUB,
				(point) -> (new RKHSFunction(
					d.k,
					d.domainLB,
					d.domainUB,
					(Xt)->(
						this.getModifiedExponentialFunction(sender, receiver, d, point).eval(Xt) 
							*  numFunc.eval(Xt)/denomFunc.eval(Xt))).computeIntegral())
				);
			senderFuncs[this.getIndexInArray(d, sender.domains)] = updateFunc;
		    }
		    receiverDomainFuncs[this.getIndexInArray(recDom, receiver.domains)] = 
			    new RKHSFunction(
				    recDom.k,
				    recDom.domainLB,
				    recDom.domainUB,
				    (point) -> (this.sumOverMessages(point, senderFuncs)));
		}
		
		// normalize the messages 
		double normalizingConstant = 0.0;
		for (RKHSFunction func : receiverDomainFuncs) { 
		    normalizingConstant += func.computeIntegral();
		}
		final double Q = normalizingConstant;
		for (int k=0; k< receiverDomainFuncs.length; k++) { 
		    RKHSFunction oldFunc = receiverDomainFuncs[k];
		    receiverDomainFuncs[k] = new RKHSFunction(
			    oldFunc.k,
			    oldFunc.domainLB,
			    oldFunc.domainUB,
			    (point) -> (oldFunc.eval(point)/Q));
		}
		
		// dump it all to the buffer hashmap
		HashMap<CMRFNodeDomain, RKHSFunction> funcMap = sender.outMessages.get(receiver);
		for (CMRFNodeDomain rD : receiver.domains) { 
		    RKHSFunction rDFunc = receiverDomainFuncs[this.getIndexInArray(rD, receiver.domains)];
		    funcMap.put(rD, rDFunc);
		}
		messageMaps.get(sender).put(receiver, funcMap);
	    }
	}
	
	// now update the messages from the buffer
	for (CMRFNode sender : nodes) { 
	    for (CMRFNode receiver : nodes) { 
		if (sender.equals(receiver)) { continue; }
		for (CMRFNodeDomain domain : receiver.domains) { 
		    sender.outMessages.get(receiver).put(domain, messageMaps.get(sender).get(receiver).get(domain));
		}
	    }
	}
	
	// let's make a poor stab at pretending we care about software engineering 	
	messageMaps.clear();
    }
    
    // returns the exponential function at a specific point for the TRBP message update
    // this is basically currying, but awkward as all hell because Java
    private RKHSFunction getModifiedExponentialFunction(
	    CMRFNode sender, 
	    CMRFNode receiver, 
	    CMRFNodeDomain domain,
	    double[] point) {
	int i = this.getIndexInArray(sender, nodes); 
	int j = this.getIndexInArray(receiver, nodes); 
	CMRFEdge edge = this.edges[i][j];
        // god i wish i could write all this crap in python
	return new RKHSFunction(
		domain.k,
		domain.domainLB,
		domain.domainUB,
		(Xt) -> (Math.exp(edge.getEnergyAtPoint(CMRFEdgeDomain.concatArrays(point, Xt))/edgeProbs[i][j] +
			sender.getDomainForPoint(Xt).getEnergyAtPoint(Xt))));
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
    
    /**
     * Splits an array at the nth slot
     * @param arr
     * @param n
     * @return 
     */
    public static ArrayList<double[]> splitArray(double[] arr, int n) { 
	ArrayList<double[]> arrs = new ArrayList<double[]>();
	arrs.add(Arrays.copyOfRange(arr, 0, n));
	arrs.add(Arrays.copyOfRange(arr, n, arr.length));
	return arrs;
    }
    
    /**
     * Concats two arrays
     * I feel really stupid having actually written this method 
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
    
    public double getLogZ() {
	return 0.0;
    }
    
    public double computeEnthalpy() { 
	return 0.0;
    }
    
    public double computeEntropy() { 
	return 0.0;
    }
}
