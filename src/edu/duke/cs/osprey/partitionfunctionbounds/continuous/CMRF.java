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
    
    double threshold = 0.1;
    double constRT = PoissonBoltzmannEnergy.constRT;
    
    boolean nodesAdded = false;
    
    /**
     * Sets up an empty cMRF -- population is done later 
     * @param numNodes 
     */
    public CMRF(int numNodes) { 
	this.numNodes = numNodes;
	nodes = new CMRFNode[numNodes];
        edges = new CMRFEdge[numNodes][numNodes];
    }
    
    public static void main (String[] args) { 
        double[][] b1 = new double[1][2]; b1[0][0] = 0; b1[0][1] = 1;
        double[] lb1 = new double[1]; lb1[0] = 0; 
        double[] ub1 = new double[1]; ub1[0] = 1;
        Kernel k1 = new KernelGaussian( b1 , 0.25 );
        CMRFNodeDomain d1 = new CMRFNodeDomain(lb1, ub1, k1, (point)->(-10) );
                
        CMRFNodeDomain[] nd1 = new CMRFNodeDomain[1]; nd1[0] = d1;
        HashMap<Integer, CMRFNodeDomain[]> doms = new HashMap<>();
        doms.put(0, nd1); 
        
        HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map1 = new HashMap<>();
        HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map2 = new HashMap<>();
        HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map3 = new HashMap<>();
        HashMap<Integer, HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>>> map4 = new HashMap<>();
        
        CMRF c = new CMRF(2);
        c.addNodes(doms, map4);        
        c.runSCMF();
    }


    /**
     * Runs the TRBP algorithm and returns an upper bound on the log partition function 
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
                        edgeEntropy += edgeProbs[nodeInd][neighborInd]*pairwiseEntropy.computeIntegral();
                    }
                }
            }
            totalEntropy += nodeEntropy;
        }
        
        return totalEntropy;
    }

    /**
     * Runs a mean-field approximation to the partition function for a lower bound on the log partition function 
     * @return 
     */
    public double runSCMF() { 
        this.initializeMarginalsSCMF();
        
        double oldEnth = Double.NEGATIVE_INFINITY;
        double oldEntr = Double.POSITIVE_INFINITY; 
        double oldLogZ = Double.NEGATIVE_INFINITY; 
        
        while (true) { 
            this.updateMarginalsSCMF();
            double enth = this.computeEnthalpySCMF();
            double entr = this.computeEntropySCMF();
            
            double freeEnergy = enth + entr;
            double logZ = freeEnergy; // freeEnergy is a direct lower bound on logZ
            
            if (logZ > oldLogZ) { 
                
                System.out.println("enth: "+enth+", entr: "+entr+", logZLB: "+logZ);
                
                if (Math.abs(logZ-oldLogZ) <= this.threshold) {
                    System.out.println("DONE: logZLB: "+logZ);
                    return logZ;
                }
                
                oldEnth = enth;
                oldEntr = entr;
                oldLogZ = logZ;
            }
        }
        
    }
    
    /**
     * Initializes marginal beliefs for SCMF -- those beliefs are just the probabilities induced by the intra-rotamer
     * energy function
     */
    public void initializeMarginalsSCMF() { 
	if (! nodesAdded) { 
	    throw new RuntimeException("Can't initialize SCMF marginals for the"
		    + "continuous-label MRF unless nodes have already been added. ");
	}
        
        System.out.print("Initializing marginals... ");
        
        for (CMRFNode node : nodes) { 
            for (CMRFNodeDomain domain : node.domains) { 
                node.marginals.put(domain, domain.probabilityRKHS);
            }
        }
        System.out.println("Done.");
        
        double enth = this.computeEnthalpySCMF();
        double entr = this.computeEntropySCMF();
        
        double freeEnergy = enth + entr;
        double logZ = freeEnergy; // freeEnergy is a direct lower bound on logZ
        
        System.out.println("logZ LB: "+logZ);
    }
    
    /**
     * Updates SCMF marginals 
     */
    public void updateMarginalsSCMF() { 
        
	/**
	 * 
	 * oneBodyE
	 * meanFieldE
	 * logUnn = -(oneBodyE + meanFieldE) / RT
	 * meanFieldE = sum(averageE)
	 * 
	 * partFn \sum{logUnn}
	 * logNB = logUnn/partFn
	 * newB = oldB*lambda + (1-lambda)*logNb
	 */
        
        // System.out.println("Updating marginals... ");
        
        // first make a buffer to hold all the new beliefs 
        HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>> newBeliefs = new HashMap<>();
        for (CMRFNode node : nodes) { 
            System.out.print("    Starting node "+this.getIndexInArray(node, nodes)+": ");
            HashMap<CMRFNodeDomain, RKHSFunction> marginals = new HashMap<>();
            for (CMRFNodeDomain domain : node.domains) { 
                System.out.print("d"+this.getIndexInArray(domain, node.domains)+"/");
                HashMap<CMRFNode, RKHSFunction> neighborFuncs = new HashMap<>();
                
                for (CMRFNode neighbor : nodes) { 
                    if (neighbor == node) { continue; }
                    System.out.print("n"+this.getIndexInArray(neighbor, nodes)+"/");
                    
                    RKHSFunction[] domainFuncs = new RKHSFunction[neighbor.domains.length];
                    for (CMRFNodeDomain neighborDomain : neighbor.domains) { 
                        System.out.print("d"+this.getIndexInArray(neighborDomain, neighbor.domains)+"/");
                        
                        // expectation of the log energy 
                        
                        RKHSFunction domainFunc = new RKHSFunction(
                                domain.k,
                                domain.domainLB,
                                domain.domainUB,
                                (point)->
                                        (this.getSCMFUpdateValue(
                                                node, 
                                                neighbor, 
                                                neighborDomain, 
                                                point)));
                                                
                        domainFuncs[this.getIndexInArray(neighborDomain, neighbor.domains)] = domainFunc;
                    }
                    
                    // sum over domains
                    RKHSFunction neighborFunc = new RKHSFunction(
                            domain.k,
                            domain.domainLB,
                            domain.domainUB,
                            (point) -> this.sumOverMessages(point, domainFuncs));
                    
                    neighborFuncs.put(neighbor, neighborFunc);
                }
                
                // sum over neighbors
                ArrayList<CMRFNode> neighbors = new ArrayList(neighborFuncs.keySet());
                RKHSFunction[] neighborMessages = new RKHSFunction[neighbors.size()];
                for (int i=0; i<neighborMessages.length; i++) { 
                    neighborMessages[i] = neighborFuncs.get(neighbors.get(i));
                }
                RKHSFunction domainFunc = new RKHSFunction(
                        domain.k,
                        domain.domainLB,
                        domain.domainUB,
                        (point)->(Math.exp(this.sumOverMessages(point, neighborMessages))));
                final RKHSFunction df = domainFunc;
                
                // normalize
                final double norm = domainFunc.computeIntegral();
                domainFunc = new RKHSFunction(
                        domain.k,
                        domain.domainLB,
                        domain.domainUB,
                        (point)->(df.eval(point)/norm));
                
                marginals.put(domain, domainFunc);
                // System.out.println("~ ");
            }
            
            double normalizingConstant = 0.0;
            for (RKHSFunction f : marginals.values()) { 
                normalizingConstant += f.computeIntegral();
            }
            final double Q = normalizingConstant;
            
            marginals.keySet().stream().forEach((domain) -> { 
                marginals.put(
                        domain, 
                        new RKHSFunction(
                                domain.k,
                                domain.domainLB,
                                domain.domainUB,
                                (point)->(marginals.get(domain).eval(point)/Q)));
            });
                        
            newBeliefs.put(node, marginals);
        }
        
        // now we load the new beliefs into the nodes 
        for (CMRFNode node : nodes) { 
            node.marginals = newBeliefs.get(node);
        }
    }

    /**
     * Computes the enthalpy when running SCMF -- note that there are no pairwise terms 
     * @return 
     */
    public double computeEnthalpySCMF() { 
        double totalEnthalpy = 0.0;
        // sum over nodes of p*E plus pariwise p*E
        for (CMRFNode v : nodes) { 
            int recNodeIndex = this.getIndexInArray(v, nodes);
            double nodeEnthalpy = 0.0; 
            
            for (CMRFNodeDomain d : v.domains) { 
                // compute single-node domain enthalpy -- Ex~Q[\ln phi]
                RKHSFunction probabilityFunc = v.marginals.get(d);
                RKHSFunction enthalpyFunc = new RKHSFunction(
                        d.k,
                        d.domainLB,
                        d.domainUB,
                        (point) -> (
                                probabilityFunc.eval(point) * Math.log(-d.energyRKHS.eval(point))));
                double domEnth = enthalpyFunc.computeIntegral();
                if (!Double.isNaN(domEnth)) { 
                    nodeEnthalpy += domEnth;
                }
//                if (Double.isNaN(domEnth)) { 
//                    
//                    System.out.println("prob: ");
//                    Matrix m = probabilityFunc.dumpPoints();
//                    m.print(1, 3);
//                    
//                    System.out.println("energy: ");
//                    m = d.energyRKHS.dumpPoints();
//                    m.print(1,3);
//
//                    System.out.println("enthalpy: ");
//                    m = enthalpyFunc.dumpPoints();
//                    m.print(3, 3);
//                    
//                    
//                    throw new RuntimeException(
//                            "NaN Enthalpy: node "+recNodeIndex+" domain "+this.getIndexInArray(d, v.domains)); 
//                }
//                nodeEnthalpy += domEnth;
                
            }
            totalEnthalpy += nodeEnthalpy;
        }
	return totalEnthalpy;
    }
    
    /**
     * Computes the entropy when running SCMF -- again, there are no pairwise terms 
     * @return 
     */
    public double computeEntropySCMF() { 
        double totalEntropy = 0.0;
        for (CMRFNode node : nodes) { 
            double nodeEntropy = 0.0;
            for (CMRFNodeDomain domain : node.domains) { 
                RKHSFunction domainPDF = node.marginals.get(domain);
                RKHSFunction domainEntropyFunc = new RKHSFunction(
                        domainPDF.k,
                        domainPDF.domainLB,
                        domainPDF.domainUB,
                        (point)->(-domainPDF.eval(point)*
                                Math.log((domainPDF.eval(point)))));
                double domainEntropy = domainEntropyFunc.computeIntegral();
//                if (Double.isNaN(domainEntropy)) { 
//                    
//                    Matrix m = domainPDF.dumpPoints();
//                    m.print(1, 7);
//                    
//                    m = domainEntropyFunc.dumpPoints();
//                    m.print(1, 7);
//                    
//                    throw new RuntimeException("NaN Entropy");
//                }
                if (!Double.isNaN(domainEntropy)) { nodeEntropy += domainEntropy; }
            }
            totalEntropy += nodeEntropy;
        }
        return totalEntropy;
    }
    
    /**
     * Returns the modified log function for the SCMF belief updates 
     * @param node
     * @param neighbor
     * @param neighborDomain
     * @param point
     * @return 
     */    
    public double getSCMFUpdateValue(
            CMRFNode node,
            CMRFNode neighbor,
            CMRFNodeDomain neighborDomain, // domain should be from the neighbor
            double[] point) { 

        CMRFNodeDomain pointDomain = node.getDomainForPoint(point);
        
        RKHSFunction neighborMarginal = neighbor.marginals.get(neighborDomain);
        RKHSFunction nodeMarginal = node.marginals.get(pointDomain);
        
        int nodeIndex = this.getIndexInArray(node, nodes);
        int neighborIndex = this.getIndexInArray(neighbor, nodes);

        RKHSFunction pairwiseEFunc = 
                this.edges[nodeIndex][neighborIndex].getEdgeDomain(
                        pointDomain,
                        neighborDomain).eFuncRKHS;
        RKHSFunction pairwisePFunc = 
                this.edges[nodeIndex][neighborIndex].getEdgeDomain(
                        pointDomain,
                        neighborDomain).pFuncRKHS;
        
        // expectation w/r/t Q of ln(-E)
        return new RKHSFunction(
                neighborDomain.k,
                neighborDomain.domainLB,
                neighborDomain.domainUB,
                (nPoint) -> 
                        (Math.log(-1*(
				pairwiseEFunc.eval(concatArrays(point, nPoint)))))).computeExpectation();
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
	for (int i=0; i<numNodes; i++) { 
	    nodes[i] = new CMRFNode(domains.get(i));
	}
               
        for (int i=0; i<nodes.length; i++) { 
            for (int j=i; j<nodes.length; j++) { 
                if (i==j) { continue; } // no self loops
                edges[i][j] = new CMRFEdge(
                        nodes[i], 
                        nodes[j], 
                        eFuncMap.get(i).get(j));
                edges[j][i] = edges[i][j];
            }
        }
	
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
