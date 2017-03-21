/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import Jama.Matrix;
import edu.duke.cs.osprey.control.Main;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
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
    
    public static void main (String[] args) { 
        
        Main.main(args);
        
	System.exit(0);
	
        System.out.println("CMRF main");
        
        double size = 1;
        double[][] b1 = new double[1][2]; b1[0][0] = 0; b1[0][1] = size;
        double[] lb1 = new double[1]; lb1[0] = 0; 
        double[] ub1 = new double[1]; ub1[0] = size;
	Kernel k1 = new KernelGaussian( b1 , 0.25*size );
	CMRFNodeDomain nd1 = new CMRFNodeDomain(lb1, ub1, k1, (point)->(-10));

	double[][] b2 = new double[2][2];
	b2[0][0] = 0; b2[0][1] = size;
	b2[1][0] = 0; b2[1][1] = size;
	double[] lb2 = {0, 0};
	double[] ub2 = {size, size};
	Kernel k2 = new KernelGaussian( b2, 0.25*size);
        CMRFNodeDomain nd2 = new CMRFNodeDomain(lb2, ub2, k2, (point)->(-20));
        
        HashMap<Integer, CMRFNodeDomain[]> h = new HashMap<>();
	h.put(0, new CMRFNodeDomain[]{nd1});
        h.put(1, new CMRFNodeDomain[]{nd2});
        
        ToDoubleFunction<double[]>f = (point)->(-15);
        HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map1 = new HashMap<>();
        map1.put(nd2, f);
        HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map2 = new HashMap<>();
        map2.put(nd1, map1);
        HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map3 = new HashMap<>();
        map3.put(1, map2);
        HashMap<Integer, HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>>> map4 = new HashMap<>();
        map4.put(0, map3);
        
        CMRF c = new CMRF(2);
	c.constRT = 1.0;
        c.addNodes(h, map4);
	c.runSCMF();
        
        c.runTRBP();
    }


    /**
     * Runs the TRBP algorithm and returns an upper bound on the log partition function 
     * @return 
     */
    public double runTRBP() { 
	System.out.print("Initializing...");
        this.initializeEdgeProbsTRBP();
        this.initializeMessagesTRBP();
	System.out.println("done.");
        
        double oldEnth = Double.POSITIVE_INFINITY;
        double oldEntr = Double.NEGATIVE_INFINITY; 
	double oldEnrg = Double.POSITIVE_INFINITY;
        double oldLogZ = Double.POSITIVE_INFINITY; 
        
	int iter = 0; 
        boolean haveValidLogZ = false;
        
        while (true) { 
            this.updateMessagesTRBP();
	    
	    System.out.print("Updating pseudomarginals...");
	    this.computeSingletonPseudomarginalsTRBP();
	    this.computePairwiseMarginalsTRBP();
	    System.out.println("done.");
	    
	    
	    
            double enth = this.computeEnthalpyTRBP();
            double entr = this.computeEntropyTRBP();
            
            double enrg = enth - this.constRT*entr;
            double logZ = Math.log(-enrg/this.constRT);
                        
            if (Double.isNaN(logZ) && !haveValidLogZ) { 
                System.out.println("logZUB is NaN, restarting...");
                return this.runTRBP();
            } else { 
               haveValidLogZ = true; 
            }


            // break if the bound gets worse, i.e. we step over a local maximum
            if ((enrg > oldEnrg || logZ > oldLogZ) && !Double.isNaN(oldLogZ)) { 
                System.out.println("DONE: logZUB: "+oldLogZ);
		printMarginalsTRBP();
                return oldLogZ;
            }

            System.out.println("enth: "+enth+", entr: "+entr+", enrg: " + enrg + ", logZUB: "+logZ);

            // break if the other termination condition is reached
            if ((Math.abs(logZ-oldLogZ) <= this.threshold) || (iter >= maxIters)) { 
                System.out.println("DONE: logZUB: "+logZ);
		printMarginalsTRBP();
                return logZ;                
            }
	    
            oldEnth = enth;
            oldEntr = entr;
            oldLogZ = logZ;
	    oldEnrg = enrg;
	    iter++;
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
    
    public void computeSingletonPseudomarginalsTRBP() { 
        // calculate singleton pseudomarginals
        for (CMRFNode node : nodes) { 
            int nodeIndex = getIndexInArray(node, nodes);

            HashMap<CMRFNodeDomain, RKHSFunction> pseudomarginals = new HashMap<>();
            double partFn = 0.0;
            
            // calculate unnormalized psueodmarginals
            for (CMRFNodeDomain domain : node.domains) { 
                RKHSFunction pFunc = domain.probabilityRKHS; 

                ArrayList<Double> powers = new ArrayList<>();
                ArrayList<RKHSFunction> neighborFuncs = new ArrayList<>();         
                
                for (int i=0; i < nodes.length; i++) { 
                    if (nodes[i] == node) { continue; }
                    CMRFNode neighbor = nodes[i];
                    powers.add(edgeProbs[nodeIndex][i]);
                    neighborFuncs.add(neighbor.outMessages.get(node).get(domain));
                }
                
                double[] edgePowers = new double[powers.size()];
                for (int i=0; i<powers.size(); i++) { edgePowers[i] = powers.get(i); }
                
                RKHSFunction[] neighborFunctions = new RKHSFunction[neighborFuncs.size()];
                for (int i=0; i<neighborFuncs.size(); i++) { neighborFunctions[i] = neighborFuncs.get(i); }
                
                pseudomarginals.put(
                        domain,
                        new RKHSFunction(
                                domain.k,
                                domain.domainLB,
                                domain.domainUB,
                                (point)->(getProdOfFuncPowers(neighborFunctions, edgePowers, point))));
                partFn += pseudomarginals.get(domain).computeIntegral();
            }
            
            // normalize
            final double Z = partFn;
            for (CMRFNodeDomain domain : node.domains) { 
                node.pseudomarginals.put(
                        domain,
                        new RKHSFunction(
                                domain.k,
                                domain.domainLB,
                                domain.domainUB,
                                (point)->(pseudomarginals.get(domain).eval(point)/Z)));
            }
        }
    }
    
    public void computePairwiseMarginalsTRBP() { 
	for (int i=0; i<edges.length; i++) { 
	    for (int j=0; j<edges.length; j++) { 
		if (j==i) { continue; }
		
		CMRFEdge edge = edges[i][j];
		CMRFNode sender = nodes[i];
		CMRFNode receiver = nodes[j];
		
		HashMap<CMRFEdgeDomain, RKHSFunction> newPairwiseMarginals = new HashMap<>();
		
		for (CMRFEdgeDomain domain : edge.domainLinks) { 
		    CMRFNodeDomain domainOne = domain.resOneDomain;
		    CMRFNodeDomain domainTwo = domain.resTwoDomain;
		    
		    // get the modified probability function
		    RKHSFunction pairwiseEFunc = domain.eFuncRKHS;
		    RKHSFunction oneEFunc = domain.resOneDomain.energyRKHS;
		    RKHSFunction twoEFunc = domain.resTwoDomain.energyRKHS;
		    final double edgeProb = edgeProbs[i][j];				    

		    RKHSFunction phiFunc = new RKHSFunction(
			    pairwiseEFunc.k,
			    pairwiseEFunc.domainLB,
			    pairwiseEFunc.domainUB,
			    (point) -> (Math.exp(-1*pairwiseEFunc.eval(point)/edgeProb 
				                    - oneEFunc.eval(splitArray(point, domainOne.domainLB.length).get(0))
				                    - twoEFunc.eval(splitArray(point, domainOne.domainLB.length).get(1)))));
		    
		    // get neighbors of sender
		    RKHSFunction[] senderFuncs = new RKHSFunction[nodes.length-1];
		    double[] senderPows = new double[nodes.length-1];
		    int snIndex = 0;
		    for (CMRFNode node : nodes) { 
			if (node == receiver || node == sender) { 
			    continue; 
			}
			senderPows[snIndex] = edgeProbs[i][getIndexInArray(node, nodes)];
			senderFuncs[snIndex] = node.outMessages.get(sender).get(domainOne);
			snIndex++;
		    }
		    RKHSFunction senderDenom = new RKHSFunction(
			    domainOne.k,
			    domainOne.domainLB,
			    domainOne.domainUB,
			    (point)->(Math.pow(receiver.outMessages.get(sender).get(domainOne).eval(point), 1-edgeProb)));
		    RKHSFunction senderFunction;
		    if (nodes.length > 2) {
			senderFunction = new RKHSFunction(
				domainOne.k,
				domainOne.domainLB,
				domainOne.domainUB,
				(point)->(getProdOfFuncPowers(senderFuncs, senderPows, point)/senderDenom.eval(point)));
		    } else {
			senderFunction = new RKHSFunction(
				domainOne.k,
				domainOne.domainLB,
				domainOne.domainUB,
				(point)->(1.0/senderDenom.eval(point)));
		    }
		    
		    // get neighbors of receiver
		    RKHSFunction[] receiverFuncs = new RKHSFunction[nodes.length-1];
		    double[] receiverPows = new double[nodes.length-1];
		    int rnIndex = 0;
		    for (CMRFNode node : nodes) { 
			if (node == receiver || node==sender) { continue; }
			receiverPows[snIndex] = edgeProbs[getIndexInArray(node, nodes)][j];
			receiverFuncs[snIndex] = node.outMessages.get(receiver).get(domainTwo);
			rnIndex++;
		    }
		    RKHSFunction receiverDenom = new RKHSFunction(
			    domainTwo.k,
			    domainTwo.domainLB,
			    domainTwo.domainUB,
			    (point)->(Math.pow(sender.outMessages.get(receiver).get(domainTwo).eval(point), 1-edgeProb)));
		    RKHSFunction receiverFunction;
		    if (nodes.length > 2) {
			receiverFunction = new RKHSFunction(
				domainTwo.k,
				domainTwo.domainLB,
				domainTwo.domainUB,
				(point)->(getProdOfFuncPowers(receiverFuncs, receiverPows, point)/receiverDenom.eval(point)));
		    } else {
			receiverFunction = new RKHSFunction(
				domainTwo.k,
				domainTwo.domainLB,
				domainTwo.domainUB,
				(point)->(1.0/receiverDenom.eval(point)));
		    }
		    
		    newPairwiseMarginals.put(
			    domain, 
			    new RKHSFunction(
				    domain.resAllK,
				    domain.domainLB,
				    domain.domainUB,
				    (point)->(phiFunc.eval(point) + 
					    senderFunction.eval(splitArray(point, domain.resOneLB.length).get(0)) + 
					    receiverFunction.eval(splitArray(point, domain.resOneLB.length).get(1)))));
		    domain.pseudomarginal = newPairwiseMarginals.get(domain);
		}
				
		edge.pseudomarginals = newPairwiseMarginals;
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
	System.out.print("Updating messages...");
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
		System.out.print(i+"-"+j+" ");
		if (i==j) { continue; }
		CMRFNode sender = nodes[i];
		CMRFNode receiver = nodes[j];
		final double edgeProb = edgeProbs[i][j];
		
		// each of the update functions is tied to a domain in the receiver
		// we'll collapse all of this when we're done constructing them
		RKHSFunction[] receiverDomainFuncs = new RKHSFunction[receiver.domains.length];
		
		for (CMRFNodeDomain recDom : receiver.domains) {
		    RKHSFunction[] senderFuncs = new RKHSFunction[sender.domains.length];
		    
		    // Xt --> sender domain (senDom)
		    // xs --> receiver domain (recDom)
		    
		    for (CMRFNodeDomain senDom : sender.domains) {
			// denominator of quotient
			ToDoubleFunction<double[]> denomFunc = 
				(Xt)->Math.pow(receiver.outMessages.get(sender).get(senDom).eval(Xt), 1-edgeProb);
				
			
			// numerator of quotient
                        // note we're looking at nodes who send to the sender here 
			ArrayList<CMRFNode> parents = new ArrayList<>();
			for (CMRFNode n : nodes) { if (!n.equals(sender) && !n.equals(receiver)) { parents.add(n); } }
			
			double[] parentProbs = new double[parents.size()];
			RKHSFunction[] parentMsgs = new RKHSFunction[parents.size()];
			
			for (int k=0; k<parentMsgs.length; k++) {
			    CMRFNode parent = parents.get(k);
			    parentMsgs[k] = parent.outMessages.get(sender).get(senDom);
			    parentProbs[k] = edgeProbs[i][this.getIndexInArray(parent, nodes)];
			}
			ToDoubleFunction<double[]> numFunc = (point) -> (getProdOfFuncPowers(parentMsgs, parentProbs, point));
			
			// now let's put it all together
			RKHSFunction updateFunc = new RKHSFunction(
				recDom.k,
				recDom.domainLB,
				recDom.domainUB,
				(xs) -> (new RKHSFunction(
					senDom.k,
					senDom.domainLB,
					senDom.domainUB,
					(Xt)->(
						this.getModifiedExponentialFunction(sender, receiver, senDom, recDom, xs).eval(Xt) 
							*  numFunc.applyAsDouble(Xt)/denomFunc.applyAsDouble(Xt)
						)).computeIntegral()));
                        // goddamn that is beautiful 
                        // we get a function for each sender domain, and sum over
			senderFuncs[this.getIndexInArray(senDom, sender.domains)] = updateFunc;
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
	System.out.println("done.");
    }
    
    /**
     * Computes the enthalpy of the cMRF in its current state using pseudomarginals
     * NOTE: pseudomarginals must have been initialized (or be extant in some form)
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
                RKHSFunction probabilityFunc = v.pseudomarginals.get(d);
                RKHSFunction enthalpyFunc = new RKHSFunction(
                        d.k,
                        d.domainLB,
                        d.domainUB,
                        (point) -> (
                                probabilityFunc.eval(point) * d.energyFunction.applyAsDouble(point)));
		double domainEnthalpy = enthalpyFunc.computeIntegral();
		if (Double.isNaN(domainEnthalpy)) { throw new RuntimeException("NaN enthalpy"); }
                nodeEnthalpy += domainEnthalpy;
                
                for (CMRFNode neighbor : nodes) { 
                    if (neighbor.equals(v)) { continue; }
                    for (CMRFNodeDomain nd : neighbor.domains) { 
                        // get the pdf for the neighbor's domain 
                        int nRecNodeInd = this.getIndexInArray(neighbor, nodes);
			CMRFEdge edge = this.edges[nRecNodeInd][recNodeIndex];
                        CMRFEdgeDomain edgeDomain = edge.getEdgeDomain(d, nd);
			
                        RKHSFunction pairwiseProbFunc = edge.pseudomarginals.get(edgeDomain);
                        RKHSFunction pairwiseEnergyFunc = edgeDomain.eFuncRKHS;
			
                        // compute enthalpy, add it to single node enthalpy 
                        RKHSFunction pairwiseEnthalpyFunc = new RKHSFunction(
                                pairwiseProbFunc.k,
                                pairwiseProbFunc.domainLB,
                                pairwiseProbFunc.domainUB,
                                (point) -> (pairwiseProbFunc.eval(point) * pairwiseEnergyFunc.eval(point)));
			double pairwiseEnthalpy = pairwiseEnthalpyFunc.computeIntegral();
			if (Double.isNaN(pairwiseEnthalpy)) { throw new RuntimeException("NaN enthalpy"); }
                        nodeEnthalpy += pairwiseEnthalpy;
                    }
                }
            }
	    if (Double.isNaN(nodeEnthalpy)) { throw new RuntimeException("NaN enthalpy"); }
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
                RKHSFunction domainPDF = node.pseudomarginals.get(domain);
                RKHSFunction domainEntropyFunc = new RKHSFunction(
                        domainPDF.k,
                        domainPDF.domainLB,
                        domainPDF.domainUB,
                        (point)->(-1*domainPDF.eval(point)*Math.log(Math.max(domainPDF.eval(point), Double.MIN_VALUE))));
                double domainEntropy = domainEntropyFunc.computeIntegral();
		if (Double.isNaN(domainEntropy)) { 
		    Matrix m = domainPDF.dumpPoints();
		    m.print(3, 5);
		    m = domainEntropyFunc.dumpPoints();
		    m.print(3, 5);
		    throw new RuntimeException("NaN entropy"); 
		}
                nodeEntropy += domainEntropy;
		
		double edgeEntropy = 0.0;
		for (CMRFNode neighbor : nodes) {
		    if (node.equals(neighbor)) { continue; }
		    int nodeInd = getIndexInArray(node, nodes);
		    int neighborInd = getIndexInArray(neighbor, nodes);
		    CMRFEdge edge = edges[nodeInd][neighborInd];

                    for (CMRFNodeDomain neighborDomain : neighbor.domains) { 
			RKHSFunction neighborPDF = neighbor.pseudomarginals.get(neighborDomain);
			CMRFEdgeDomain edgeDomain = edge.getEdgeDomain(domain, neighborDomain);
			RKHSFunction pairwisePDF = edge.pseudomarginals.get(edgeDomain);
			
			RKHSFunction pairwiseEntropy = new RKHSFunction(
				edgeDomain.resAllK,
				edgeDomain.domainLB,
				edgeDomain.domainUB,
				(point) -> (pairwisePDF.eval(point) *
					Math.log(Math.max(
						pairwisePDF.eval(point)/
							(domainPDF.eval(splitArray(point, domainPDF.domainLB.length).get(0)) *
								neighborPDF.eval(splitArray(point, domainPDF.domainLB.length).get(1))),
						Double.MIN_VALUE))));
			
			double pairEntropy = edgeProbs[nodeInd][neighborInd]*pairwiseEntropy.computeIntegral();
			if (Double.isNaN(pairEntropy)) {
			    Matrix m = pairwiseEntropy.dumpPoints();
			    m.print(3, 5);
			    throw new RuntimeException("NaN entropy");
			}
			
			edgeEntropy += pairEntropy;
                    }
                }
		if (Double.isNaN(edgeEntropy)) { throw new RuntimeException("NaN entropy"); }
		nodeEntropy += edgeEntropy;
            }
	    if (Double.isNaN(nodeEntropy)) { throw new RuntimeException("NaN entropy"); }
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
        int iter = 0;
        
        double oldEnth = Double.NEGATIVE_INFINITY;
        double oldEntr = Double.POSITIVE_INFINITY;
        double oldLogZ = Double.NEGATIVE_INFINITY;
        
        while (true) {
            this.updateMarginalsSCMF();
            double enth = this.computeEnthalpySCMF();
            double entr = this.computeEntropySCMF();
            
            double freeEnergy = enth + entr;
            double logZ = freeEnergy; // freeEnergy is a direct lower bound on logZ
            
            // break if the bound gets worse, i.e. we step over a local maximum
            if (logZ < oldLogZ) { 
                ranSCMF = true;
                System.out.println("DONE: logZLB: "+oldLogZ);
		printMarginalsSCMF();
                return oldLogZ;
            }

            System.out.println("enth: "+enth+", entr: "+entr+", logZLB: "+logZ);

            // break if the other termination condition is reached
            if ((Math.abs(logZ-oldLogZ) <= this.threshold) || (iter >= maxIters)) { 
                ranSCMF = true;
                System.out.println("DONE: logZLB: "+logZ);
		printMarginalsSCMF();
                return logZ;                
            }

            oldEnth = enth;
            oldEntr = entr;
            oldLogZ = logZ;
            iter++;
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
        HashMap<CMRFNode, RKHSFunction[]> newBeliefs = new HashMap<>();
        for (CMRFNode node : nodes) {
            
            double partFn = 0.0;
            RKHSFunction[] domainMarginals = new RKHSFunction[node.domains.length];
            
            for (CMRFNodeDomain domain : node.domains) {
                int domainIndex = Arrays.asList(node.domains).indexOf(domain);
                RKHSFunction oneBodyEnergy = domain.energyRKHS;
                
                // get a list of average interaction energy functions
                ArrayList<RKHSFunction> avgInteractionEnergies = new ArrayList<>();
                for (CMRFNode neighbor : nodes) {
                    int nodeInd = Arrays.asList(nodes).indexOf(node);
                    int neighborInd = Arrays.asList(nodes).indexOf(node);
                    if (neighborInd!=nodeInd) {
                        for (CMRFNodeDomain nDomain : neighbor.domains) {
                            avgInteractionEnergies.add(
                                    new RKHSFunction(
                                            domain.k,
                                            domain.domainLB,
                                            domain.domainUB,
                                            (point)->(pairwiseExpectation(
                                                    point,
                                                    node,
                                                    domain,
                                                    neighbor,
                                                    nDomain))));
                        }
                    }
                }
                
                // mean field energy is the sum of these average interactions
		RKHSFunction[] avgInteractionE = 
			new RKHSFunction[avgInteractionEnergies.size()];
		for (int i=0; i<avgInteractionE.length; i++) { 
		    avgInteractionE[i] = avgInteractionEnergies.get(i); 
		}
		RKHSFunction meanFieldE = new RKHSFunction(
			domain.k,
			domain.domainLB,
			domain.domainUB,
			(point)->sumOverMessages(point, avgInteractionE));
		
		RKHSFunction logUnnormalizedBelief = new RKHSFunction(
			domain.k,
			domain.domainLB,
			domain.domainUB,
			(point) -> (-(oneBodyEnergy.eval(point) + meanFieldE.eval(point))/this.constRT));
		
		final double logBeliefNormalizer = logUnnormalizedBelief.computeIntegral();
		RKHSFunction logNormalizedBelief = new RKHSFunction(
			domain.k,
			domain.domainLB,
			domain.domainUB,
			(point)->(logUnnormalizedBelief.eval(point)/logBeliefNormalizer));
		
		RKHSFunction marginalBelief = new RKHSFunction(
			domain.k,
			domain.domainLB,
			domain.domainUB,
			(point)->(Math.exp(logNormalizedBelief.eval(point))));
		
		domainMarginals[domainIndex] = marginalBelief;
		partFn += marginalBelief.computeIntegral();
	    }
	    
	    // normalize marginals to get beliefs
	    for (int i=0; i<domainMarginals.length; i++) { 
		final double pF = partFn;
		final RKHSFunction[] dMs = domainMarginals;
		final RKHSFunction marg = domainMarginals[i];
		domainMarginals[i] = new RKHSFunction(
			dMs[i].k,
			dMs[i].domainLB,
			dMs[i].domainUB,
			(point)->(marg.eval(point)/pF));
	    }
	    newBeliefs.put(node, domainMarginals);
	}
	
	// now we update our beliefs
	for (CMRFNode node : nodes) { 
	    RKHSFunction[] domainMarginals = newBeliefs.get(node);
	    for (int i=0; i<node.domains.length; i++) { 
		final CMRFNodeDomain domain = node.domains[i];
		final RKHSFunction oldM = node.marginals.get(domain);
		final RKHSFunction newM = domainMarginals[i];
		final double alpha = lambda;
		node.marginals.put(domain, new RKHSFunction(
			domain.k,
			domain.domainLB,
			domain.domainUB,
			(point)->(alpha*newM.eval(point) + (1-alpha)*oldM.eval(point))));
	    }
	}
	
	lambda = lambda * (1 - 1.0/maxIters);
    }
    
    public double pairwiseExpectation(
	    double[] point, 
	    CMRFNode node,
	    CMRFNodeDomain domain,
	    CMRFNode neighbor,
	    CMRFNodeDomain nDomain) { 
	int nodeIndex = Arrays.asList(nodes).indexOf(node);
	int neighborIndex = Arrays.asList(nodes).indexOf(node);
	
	CMRFEdge e = edges[nodeIndex][neighborIndex];
	CMRFEdgeDomain ed = e.getEdgeDomain(domain, nDomain);
	
	RKHSFunction func = new RKHSFunction(
		nDomain.k,
		nDomain.domainLB,
		nDomain.domainUB,
		(nPoint)->(ed.eFunc.applyAsDouble(CMRFEdgeDomain.concatArrays(point, nPoint))));
	
	return func.computeExpectation();
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
                                probabilityFunc.eval(point) * Math.log(-d.energyFunction.applyAsDouble(point))));
                double domEnth = enthalpyFunc.computeIntegral();
                if (!Double.isNaN(domEnth)) { 
                    nodeEnthalpy += domEnth;
		}
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
                if (!Double.isNaN(domainEntropy)) { 
		    nodeEntropy += domainEntropy; 
		}
            }
            totalEntropy += nodeEntropy;
        }
        return totalEntropy;
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
    
    // returns the exponential function at a specific point for the TRBP message update
    // this is basically currying, but awkward as all hell because Java
    private RKHSFunction getModifiedExponentialFunction(
	    CMRFNode sender, 
	    CMRFNode receiver, 
	    CMRFNodeDomain senDomain,
	    CMRFNodeDomain recDomain,
	    double[] xs) {
	int i = this.getIndexInArray(sender, nodes); 
	int j = this.getIndexInArray(receiver, nodes); 
	CMRFEdge edge = this.edges[i][j];
        // god i wish i could write all this crap in python
	return new RKHSFunction(
		senDomain.k,
		senDomain.domainLB,
		senDomain.domainUB,
		(Xt) -> (Math.exp(
			edge.getEnergyAtPoint(
				CMRFEdgeDomain.concatArrays(Xt, xs))
				/edgeProbs[i][j] +
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
    
    public void printMarginalsSCMF() { 
	System.out.print("Printing SCMF marginals... ");
	if (!nodesAdded) { return; } // do nothing if we don't have an MRF 
	for (int i=0; i<nodes.length; i++) { 
	    CMRFNode node = nodes[i];
	    for (int j=0; j<node.domains.length; j++) { 
		 CMRFNodeDomain domain = node.domains[j];
		 String filename = "scmf_"+i+"-"+j+".dat";
		 try { 
		     PrintWriter writer = new PrintWriter(filename, "UTF-8");
		     Matrix m = node.marginals.get(domain).dumpPoints();
		     m.print(writer, 3, 5);
		     writer.flush();
		 } catch(FileNotFoundException | UnsupportedEncodingException e) { 
		     throw new RuntimeException("PrintWriting failed for "+ 
			     "node " + i +", domain " + j +".\n" +e.getMessage());
		 }
	    }
	}
	System.out.println("done.");
    }
    
    public void printMarginalsTRBP() {
	System.out.print("Printing TRBP pseudomarginals... ");
	if (!nodesAdded) { return; } // do nothing if we don't have an MRF
	
	for (int i=0; i<nodes.length; i++) {
	    CMRFNode node = nodes[i];
	    for (int j=0; j<node.domains.length; j++) {
		
		// print unitary pseudomarginals
		CMRFNodeDomain domain = node.domains[j];
		try {
		    String filename = "trbp_u_"+i+"-"+j+".dat";
		    PrintWriter writer = new PrintWriter(filename, "UTF-8");
		    Matrix m = node.pseudomarginals.get(domain).dumpPoints();
		    m.print(writer, 3, 5);
		    writer.flush();
		} catch(FileNotFoundException | UnsupportedEncodingException e) {
		    throw new RuntimeException("PrintWriting failed for "+
			    "node " + i +", domain " + j +".\n" +e.getMessage());
		}
		
		// print pairwise pseudomarginals
		for (int k=0; k<nodes.length; k++) {
		    CMRFEdge edge = edges[i][k];
		    for (int l=0; l<edge.domainLinks.length; k++) {
			CMRFEdgeDomain edgeDomain = edge.domainLinks[l];
			int x1 = getIndexInArray(edgeDomain.resOneDomain, node.domains);
			int x2 = getIndexInArray(edgeDomain.resTwoDomain, nodes[k].domains);
			
			try {
			    String filename = "trbp_p_"+i+k+"_"+x1+x2+".dat";
			    PrintWriter writer = new PrintWriter(filename, "UTF-8");
			    Matrix m = edge.pseudomarginals.get(edgeDomain).dumpPoints();
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
		    CMRFEdge edge = edges[i][k];
		    for (int l=0; l<edge.domainLinks.length; k++) {
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
