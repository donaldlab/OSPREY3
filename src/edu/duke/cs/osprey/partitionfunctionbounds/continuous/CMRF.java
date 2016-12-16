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
    
    double threshold = 1e-6;
    double constRT = PoissonBoltzmannEnergy.constRT;
    
    boolean nodesAdded = false;
    
    /**
     * Sets up an empty cMRF -- population is done later 
     * @param numNodes 
     */
    public CMRF(int numNodes) { 
	this.numNodes = numNodes;
	nodes = new CMRFNode[numNodes];
    }
    


    /**
     * Skeleton stub for the running of TRBP 
     * @return 
     */
    public double runTRBP() { 
        this.initializeEdgeProbsTRBP();
        this.initializeMessagesTRBP();
        
        double oldEnth = Double.POSITIVE_INFINITY;
        double oldEntr = Double.NEGATIVE_INFINITY; 
        double oldLogZ = Double.POSITIVE_INFINITY; 
        
        while (true) { 
            this.updateMessagesTRBP();
            double enth = this.computeEnthalpyTRBP();
            double entr = this.computeEntropyTRBP();
            
            double freeEnergy = enth - this.constRT*entr;
            double logZ = Math.log(-freeEnergy/this.constRT);
            if (Math.abs(logZ-oldLogZ) <= this.threshold) { 
                return logZ;
            }
            
            oldEnth = enth;
            oldEntr = entr;
            oldLogZ = logZ;
        }
    }    

    /**
     * Initializes update messages as uniform distributions 
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
				new HashMap<>());
		    }
		    sender.outMessages.get(receiver).put(d, message);
		}
	    }
	}	
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
     * TRBP message update procedure 
     */
    public void updateMessagesTRBP() {
	// first, construct the messages; then, update them all at once 

	// we'll store all the messages in a hashmap buffer
	// map from sender --> receiver --> receiver domain --> rkhsFunction
	HashMap<CMRFNode, HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>>> messageMaps = 
		new HashMap<>();
	for (CMRFNode sender : nodes) { 
	    HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>> senderMap = 
		    new HashMap<>();
	    for (CMRFNode receiver : nodes) {  
		if (sender.equals(receiver)) { continue; }
		HashMap<CMRFNodeDomain, RKHSFunction> domainMap = new HashMap<>();
		for (CMRFNodeDomain domain : receiver.domains) { 
		    domainMap.put(
			    domain,
			    new RKHSFunction( // default to the uniform distribution, but this is just a placeholder
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
                        // note we're looking at nodes who send to the sender here 
			ArrayList<CMRFNode> parents = new ArrayList<>();
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
			
			// now let's put it all together
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
							*  numFunc.eval(Xt)/denomFunc.eval(Xt))).computeIntegral()));
                        // goddamn that is beautiful 
                        // we get a function for each sender domain, and sum over
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
		    sender.outMessages.get(receiver).put(
                            domain, 
                            messageMaps.get(sender).get(receiver).get(domain));
		}
	    }
	}
	
	// let's make a poor stab at pretending we care about software engineering 	
	messageMaps.clear();
    }
    
    /**
     * Computes the enthalpy of the cMRF in its current state 
     * @return 
     */
    public double computeEnthalpyTRBP() { 
        double totalEnthalpy = 0.0;
        // sum over nodes of p*E plus pariwise p*E
        for (CMRFNode v : nodes) { 
            int recNodeIndex = this.getIndexInArray(v, nodes);
            double nodeEnthalpy = 0.0; 
            
            for (CMRFNodeDomain d : v.domains) { 
                // compute single-node domain enthalpy 
                RKHSFunction probabilityFunc = this.getPDF(v, d);
                RKHSFunction enthalpyFunc = new RKHSFunction(
                        d.k,
                        d.domainLB,
                        d.domainUB,
                        (point) -> (
                                probabilityFunc.eval(point) * d.energyFunction.applyAsDouble(point)));
                nodeEnthalpy += enthalpyFunc.computeIntegral();
                
                for (CMRFNode neighbor : nodes) { 
                    if (neighbor.equals(v)) { continue; }
                    for (CMRFNodeDomain nd : neighbor.domains) { 
                        // get the pdf for the neighbor's domain 
                        int nRecNodeInd = this.getIndexInArray(neighbor, nodes);
                        CMRFEdgeDomain edgeDomain = this.edges[nRecNodeInd][recNodeIndex].getEdgeDomain(d, nd);
                        RKHSFunction pairwiseProbFunc = edgeDomain.pFuncRKHS;
                        RKHSFunction pairwiseEnergyFunc = edgeDomain.eFuncRKHS;
                        
                        // compute enthalpy, add it to single node enthalpy 
                        RKHSFunction pairwiseEnthalpyFunc = new RKHSFunction(
                                pairwiseProbFunc.k,
                                pairwiseProbFunc.domainLB,
                                pairwiseProbFunc.domainUB,
                                (point) -> (pairwiseProbFunc.eval(point) * pairwiseEnergyFunc.eval(point)));
                        nodeEnthalpy += pairwiseEnthalpyFunc.computeIntegral();
                    }
                }
            }
            totalEnthalpy += nodeEnthalpy;
        }
	return totalEnthalpy;
    }
    
    /**
     * Computes the entropy of the cMRF in its current state 
     * @return 
     */
    public double computeEntropyTRBP() { 
        double totalEntropy = 0.0;
        for (CMRFNode node : nodes) { 
            double nodeEntropy = 0.0;
            for (CMRFNodeDomain domain : node.domains) { 
                RKHSFunction domainPDF = this.getPDF(node, domain);
                RKHSFunction domainEntropyFunc = new RKHSFunction(
                        domainPDF.k,
                        domainPDF.domainLB,
                        domainPDF.domainUB,
                        (point)->(-1*domainPDF.eval(point)*Math.log(domainPDF.eval(point))));
                double domainEntropy = domainEntropyFunc.computeIntegral();
                nodeEntropy += domainEntropy;
                
                double edgeEntropy = 0.0;
                for (CMRFNode neighbor : nodes) { 
                    if (node.equals(neighbor)) { continue; } 
                    for (CMRFNodeDomain neighborDomain : neighbor.domains) { 
                        int nodeInd = this.getIndexInArray(node, nodes);
                        int neighborInd = this.getIndexInArray(neighbor, nodes);
                        CMRFEdgeDomain edgeDomain = 
                                this.edges[nodeInd][neighborInd].getEdgeDomain(domain, neighborDomain);
                        
                        RKHSFunction pairwiseProbFunc = edgeDomain.pFuncRKHS;
                        RKHSFunction pairwiseEntropy = new RKHSFunction(
                                pairwiseProbFunc.k,
                                pairwiseProbFunc.domainLB,
                                pairwiseProbFunc.domainUB,
                                (point) -> (-1*pairwiseProbFunc.eval(point)*Math.log(pairwiseProbFunc.eval(point))));
                        edgeEntropy += pairwiseEntropy.computeIntegral();
                    }
                }
            }
            totalEntropy += nodeEntropy;
        }
        
        return totalEntropy;
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
	
        edgeEnergyFuncMap.keySet().stream().forEach((i) -> { 
            edgeEnergyFuncMap.get(i).keySet().stream().forEach((j) -> {
                edges[i][j] = new CMRFEdge(nodes[i], nodes[j], edgeEnergyFuncMap.get(i).get(j));
            });
        });
	
	nodesAdded = true;
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
}
