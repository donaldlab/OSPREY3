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
    
    public void updateMessagesTRBP() {
	// first, construct the messages; then, update them 
	// sender -> receiver -> domain -> function
	HashMap<CMRFNode, HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>> messageMap = 
		new HashMap<CMRFNode, HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>>();

	// initialize all the empty maps
	for (int i=0; i<numNodes; i++) { 
	    messageMap.put(nodes[i], new HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>());
	    HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>> senderMap = 
		    messageMap.get(nodes[i]);
	    for (int j=0; j<numNodes; j++) { 
		if (i==j) { continue; }
		senderMap.put(nodes[j], new HashMap<CMRFNodeDomain, RKHSFunction>());
	    }
	}
	
	for (int i=0; i<numNodes; i++) {  // i is the sender
	    for (int j=0; j<numNodes; j++) { // j is the receiver
		if (i==j) { continue; } 
		CMRFNode sender = nodes[i]; 
		CMRFNode receiver = nodes[j];
		final double edgeProb = edgeProbs[i][j];
		
		// eventually have to take a sum over domain functoins
		RKHSFunction[] domainFuncs = new RKHSFunction[nodes[j].domains.length];
		
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

		    RKHSFunction intEnergy = d.energyRKHS;
		    
		    // get function for the domain

		    domainFuncs[this.getIndexInArray(d, sender.domains)] = update;
		}
		
		RKHSFunction message = new RKHSFunction(
			
		)

	    } 
	}
	
    }
    
    public double getTRBPMessageUpdateValue(CMRFNode sender, CMRFNode receiver, double[] point) { 

    }
    
    public ToDoubleFunction<double[]> getTRBPMessageUpdateFunction(
	    CMRFNode sender, 
	    CMRFNode receiver, 
	    RKHSFunction numFunc,
	    RKHSFunction denomFunc,
	    double[] Xt) { 
	int sInd = this.getIndexInArray(sender, nodes);
	int rInd = this.getIndexInArray(receiver, nodes);
	CMRFEdge edge = edges[sInd][rInd];
	return (xs)->
		(Math.exp(edge.getEnergyAtPoint(concatArrays(xs, Xt))/edgeProbs[sInd][rInd] + 
			sender.getDomainForPoint(Xt).getEnergyAtPoint(Xt)) * numFunc/denomFunc);
    }
    
    public double sumOverMessages(double[] point, RKHSFunction[] funcs) { 
	double res = 0.0;
	for (RKHSFunction func : funcs) {
	    res += func.eval(point);
	}
	return res;
    }
    
    public static ArrayList<double[]> splitArray(double[] arr, int n) { 
	ArrayList<double[]> arrs = new ArrayList<double[]>();
	arrs.add(Arrays.copyOfRange(arr, 0, n));
	arrs.add(Arrays.copyOfRange(arr, n, arr.length));
	return arrs;
    }
    
    public static double[] concatArrays(double[] arr1, double[] arr2) { 
	return CMRFEdgeDomain.concatArrays(arr1, arr2);
    }
    
    public double getProdOfFuncPowers(RKHSFunction[] funcs, double[] powers, double[] point) { 
	double result = 1.0;
	for (int i=0; i<funcs.length; i++) { 
	    result *= Math.pow(funcs[i].eval(point), powers[i]);
	}
	return result;
    }
    
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
