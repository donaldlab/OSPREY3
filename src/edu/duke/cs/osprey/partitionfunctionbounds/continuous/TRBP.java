package edu.duke.cs.osprey.partitionfunctionbounds.continuous;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.function.ToDoubleFunction;

import Jama.Matrix;

public class TRBP {

	CMRF cmrf;
	double entropyMult = 10E5; 
	
	public TRBP(CMRF cmrf) {
		this.cmrf = cmrf;
	}

	/**
	 * Runs the TRBP algorithm and returns an upper bound on the log partition function 
	 * @return 
	 */
	public double runTRBP(int it) { 
		System.out.print("Initializing TRBP...");
		this.initializeEdgeProbsTRBP();
		this.initializeEdgeWeightsTRBP();
		this.initializeMessagesTRBP(it);
		System.out.println("Done.");

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
			
			this.updateEdgeProbsTRBP(iter);

			double enth = this.computeEnthalpyTRBP();
			double entr = this.computeEntropyTRBP();

			double enrg = enth - cmrf.constRT*entr;
			double logZ = -enrg/cmrf.constRT;

			System.out.println("enth: "+enth+", entr: "+entr+", enrg: " + enrg + ", logZUB: "+logZ);			
			
			if (Double.isNaN(logZ) && !haveValidLogZ) { 
				System.out.println("logZUB is NaN, restarting...");
				return this.runTRBP(it);
			} else { 
				haveValidLogZ = true; 
			}


			// break if things are dying
			boolean energyWorse = (Math.abs(enrg) - Math.abs(oldEnrg) > 0);
			if ((energyWorse || logZ > oldLogZ) && Double.isNaN(oldLogZ)) {
				if (energyWorse) { System.out.println("energy got worse"); }
				if (Double.isNaN(logZ)) { System.out.println("Ended on a NaN"); }
				printMarginalsTRBP();
				System.out.println("DONE: logZUB: "+oldLogZ);
				System.out.println("Fenth: "+enth+", Fentr: "+entr+", Fenrg: " + enrg + ", FlogZUB: "+logZ);

				return oldLogZ;
			}


			// break if the other termination condition is reached
			if ((Math.abs(logZ-oldLogZ) <= cmrf.threshold) || (iter >= 50)) { 
				printMarginalsTRBP();
				System.out.println("DONE: logZUB: "+Math.min(logZ, oldLogZ));
				return Math.min(logZ, oldLogZ);                
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
	public void initializeMessagesTRBP(int lbpIters) { 
		for (int i=0; i<cmrf.numNodes; i++) { 
			CMRFNode sender = cmrf.nodes[i];
			for (int j=0; j<cmrf.numNodes; j++) { 
				CMRFNode receiver = cmrf.nodes[j];
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
		
		calculateLBPMessagesTRBP(lbpIters);
	}
	
	/**
	 * Initializes the TRBP messages to the messages obtained after 1000 runs of the 
	 * loopy BP message update procedure 
	 */
	public void calculateLBPMessagesTRBP(int iter) { 
		System.out.println();
		for (int it=0; it<iter; it++) { 
			System.out.print("Running LBP, iter " + (it+1) + " of " + iter + "... ");
			for (int i=0; i<cmrf.numNodes; i++) {
				System.out.print(i+"{");
				CMRFNode sender = cmrf.nodes[i];

				for (int j=0; j<cmrf.numNodes; j++) { 
					if (i==j) { continue; }
					
					System.out.print(j);
					
					CMRFNode receiver = cmrf.nodes[j];
					HashMap<CMRFNodeDomain, RKHSFunction> senderMessages = new HashMap<>();

					// compute unNormalized message updates
					for (CMRFNodeDomain Wt : receiver.domains) { 

						RKHSFunction[] senderFuncs = new RKHSFunction[sender.domains.length]; 
						int sendFuncInd = 0;
						for (CMRFNodeDomain Ws : sender.domains) {
							RKHSFunction pairwiseProbFunc = cmrf.edges[i][j].getEdgeDomain(Ws, Wt).pFuncRKHS;
							RKHSFunction senderProbFunc = Ws.probabilityRKHS;

							RKHSFunction[] parentFuncs = new RKHSFunction[cmrf.numNodes-2];
							double[] parentPows = new double[cmrf.numNodes-2];
							int pFuncInd = 0;
							for (CMRFNode parent : cmrf.nodes) { 
								if (parent==receiver || parent==sender) { continue; } 
								parentFuncs[pFuncInd] = parent.outMessages.get(sender).get(Ws);
								parentPows[pFuncInd] = 1.0;
								pFuncInd++;
							}

							RKHSFunction parentMsgFunc = new RKHSFunction(
									Ws.k,
									Ws.domainLB,
									Ws.domainUB,
									(point)->(cmrf.getProdOfFuncPowers(parentFuncs, parentPows, point)));

							RKHSFunction unNormalizedSenderFunc = new RKHSFunction(
									Wt.k,
									Wt.domainLB,
									Wt.domainUB,
									(Xt)->(new RKHSFunction(
											Ws.k,
											Ws.domainLB,
											Ws.domainUB,
											(Xs)->(pairwiseProbFunc.eval(cmrf.concatArrays(Xs, Xt))
													*senderProbFunc.eval(Xs)
													*parentMsgFunc.eval(Xs))).computeIntegral()));
							senderFuncs[sendFuncInd] = unNormalizedSenderFunc;
							sendFuncInd++;
						}
						RKHSFunction unNormalizedReceiverFunc = new RKHSFunction(
								Wt.k,
								Wt.domainLB,
								Wt.domainUB,
								(point)->(cmrf.sumOverMessages(point, senderFuncs)));
						senderMessages.put(Wt, unNormalizedReceiverFunc);
					}

					// normalize message updates
					double partFn = 0.0;
					for (CMRFNodeDomain Wt : receiver.domains) { 
						partFn += senderMessages.get(Wt).computeIntegral();
					}
					final double Q = partFn;
					for (CMRFNodeDomain Wt : receiver.domains) { 
						senderMessages.put(
								Wt,
								new RKHSFunction(
										Wt.k,
										Wt.domainLB,
										Wt.domainUB,
										(point)->(senderMessages.get(Wt).eval(point)/Q)));
					}
				}
				System.out.print("}");
			}
			System.out.println("done.");
		}
	}

	public void computeSingletonPseudomarginalsTRBP() { 
		// calculate singleton pseudomarginals
		for (CMRFNode node : cmrf.nodes) { 
			int nodeIndex = cmrf.getIndexInArray(node, cmrf.nodes);

			HashMap<CMRFNodeDomain, RKHSFunction> pseudomarginals = new HashMap<>();
			double partFn = 0.0;

			// calculate unnormalized psueodmarginals
			for (CMRFNodeDomain domain : node.domains) { 
				RKHSFunction pFunc = domain.probabilityRKHS; 

				ArrayList<Double> powers = new ArrayList<>();
				ArrayList<RKHSFunction> neighborFuncs = new ArrayList<>();         

				for (int i=0; i < cmrf.nodes.length; i++) { 
					if (cmrf.nodes[i] == node) { continue; }
					CMRFNode neighbor = cmrf.nodes[i];
					powers.add(cmrf.edgeProbs[nodeIndex][i]);
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
								(point)->(cmrf.getProdOfFuncPowers(neighborFunctions, edgePowers, point))));
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
		for (int i=0; i<cmrf.edges.length; i++) { 
			for (int j=0; j<cmrf.edges.length; j++) { 
				if (j==i) { continue; }

				CMRFEdge edge = cmrf.edges[i][j];
				CMRFNode sender = cmrf.nodes[i];
				CMRFNode receiver = cmrf.nodes[j];

				HashMap<CMRFEdgeDomain, RKHSFunction> newPairwiseMarginals = new HashMap<>();

				for (CMRFEdgeDomain domain : edge.domainLinks) { 
					CMRFNodeDomain domainOne = domain.resOneDomain;
					CMRFNodeDomain domainTwo = domain.resTwoDomain;

					// get the modified probability function
					RKHSFunction pairwiseEFunc = domain.eFuncRKHS;
					RKHSFunction oneEFunc = domain.resOneDomain.energyRKHS;
					RKHSFunction twoEFunc = domain.resTwoDomain.energyRKHS;
					final double edgeProb = cmrf.edgeProbs[i][j];				    

					RKHSFunction phiFunc = new RKHSFunction(
							pairwiseEFunc.k,
							pairwiseEFunc.domainLB,
							pairwiseEFunc.domainUB,
							(point) -> (Math.exp(-1*pairwiseEFunc.eval(point)/edgeProb
									- oneEFunc.eval(CMRF.splitArray(point, domainOne.domainLB.length).get(0))
									- twoEFunc.eval(CMRF.splitArray(point, domainOne.domainLB.length).get(1)))));

					// get neighbors of sender
					RKHSFunction[] senderFuncs = new RKHSFunction[cmrf.nodes.length-2];
					double[] senderPows = new double[cmrf.nodes.length-2];
					int snIndex = 0;
					for (CMRFNode node : cmrf.nodes) { 
						if (node == receiver || node == sender) { 
							continue; 
						}
						senderPows[snIndex] = cmrf.edgeProbs[i][cmrf.getIndexInArray(node, cmrf.nodes)];
						senderFuncs[snIndex] = node.outMessages.get(sender).get(domainOne);
						snIndex++;
					}
					RKHSFunction senderDenom = new RKHSFunction(
							domainOne.k,
							domainOne.domainLB,
							domainOne.domainUB,
							(point)->(Math.pow(receiver.outMessages.get(sender).get(domainOne).eval(point), 1-edgeProb)));
					RKHSFunction senderFunction;
					if (cmrf.nodes.length > 2) {
						senderFunction = new RKHSFunction(
								domainOne.k,
								domainOne.domainLB,
								domainOne.domainUB,
								(point)->(cmrf.getProdOfFuncPowers(senderFuncs, senderPows, point)/senderDenom.eval(point)));
					} else {
						senderFunction = new RKHSFunction(
								domainOne.k,
								domainOne.domainLB,
								domainOne.domainUB,
								(point)->(1.0/senderDenom.eval(point)));
					}

					// get neighbors of receiver
					RKHSFunction[] receiverFuncs = new RKHSFunction[cmrf.nodes.length-2];
					double[] receiverPows = new double[cmrf.nodes.length-2];
					int rnIndex = 0;
					for (CMRFNode node : cmrf.nodes) { 
						if (node == receiver || node==sender) { continue; }
						receiverPows[rnIndex] = cmrf.edgeProbs[cmrf.getIndexInArray(node, cmrf.nodes)][j];
						receiverFuncs[rnIndex] = node.outMessages.get(receiver).get(domainTwo);
						rnIndex++;
					}
					RKHSFunction receiverDenom = new RKHSFunction(
							domainTwo.k,
							domainTwo.domainLB,
							domainTwo.domainUB,
							(point)->(Math.pow(sender.outMessages.get(receiver).get(domainTwo).eval(point), 1-edgeProb)));
					RKHSFunction receiverFunction;
					if (cmrf.nodes.length > 2) {
						receiverFunction = new RKHSFunction(
								domainTwo.k,
								domainTwo.domainLB,
								domainTwo.domainUB,
								(point)->(cmrf.getProdOfFuncPowers(receiverFuncs, receiverPows, point)/receiverDenom.eval(point)));
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
									(point)->(phiFunc.eval(point) * 
											senderFunction.eval(CMRF.splitArray(point, domain.resOneLB.length).get(0)) * 
											receiverFunction.eval(CMRF.splitArray(point, domain.resOneLB.length).get(1)))));
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
		if (! cmrf.nodesAdded) { 
			throw new RuntimeException("Can't initialize edge probabilities for the"
					+ "continuous-label MRF unless nodes have already been added. ");
		}

		if (cmrf.numNodes==1) { 
			return;
		}

		double[][] adjacency = new double[cmrf.numNodes][cmrf.numNodes];
		for (int i=0; i<cmrf.numNodes; i++) { 
			for (int j=0; j<cmrf.numNodes; j++) { 
				if (i==j) { 
					adjacency[i][j] = 0; 
				} else { 
					adjacency[i][j] = 1;
				} 
			}
		}
		Matrix adj = Jama.Matrix.constructWithCopy(adjacency);
		
		Matrix deg = Jama.Matrix.identity(cmrf.numNodes, cmrf.numNodes).times(cmrf.numNodes);
		Matrix laplacian = deg.minus(adj);
		Matrix lapMinOne = laplacian.minus(new Matrix(cmrf.numNodes, cmrf.numNodes, 1.0));
		Matrix invLap = lapMinOne.inverse();

		cmrf.edgeProbs = new double[cmrf.numNodes][cmrf.numNodes];
		for (int i=0; i<cmrf.edgeProbs.length; i++) { 
			for (int j=0; j<cmrf.edgeProbs[i].length; j++) { 
				cmrf.edgeProbs[i][j] = 
						adj.get(i, j) * (invLap.get(i, i) + invLap.get(j, j) - 2*invLap.get(i, j));
				if (Double.isNaN(cmrf.edgeProbs[i][j])) { 
					throw new RuntimeException("Edge probability is NaN... how did you do this?");
				}
			}
		}
		
	}

	public void initializeEdgeWeightsTRBP() { 
		cmrf.edgeWeights = new double[cmrf.numNodes][cmrf.numNodes];
		for (int i=0; i<cmrf.edgeWeights.length; i++) { 
			for (int j=0; j<i; j++) { 
				CMRFEdge edge = cmrf.edges[i][j];
				double mutInf = 0.0;
				for (CMRFEdgeDomain edgeDomain : edge.domainLinks) { 
					RKHSFunction pairwisePDF = edgeDomain.pFuncRKHS;
					RKHSFunction domainPDF = edgeDomain.resOneDomain.probabilityRKHS;
					RKHSFunction neighborPDF = edgeDomain.resTwoDomain.probabilityRKHS;
					RKHSFunction pairwiseEntropy = new RKHSFunction(
							edgeDomain.resAllK,
							edgeDomain.domainLB,
							edgeDomain.domainUB,
							(point) -> (pairwisePDF.eval(point) *
									Math.log(Math.max(
											pairwisePDF.eval(point)/
											(domainPDF.eval(CMRF.splitArray(point, domainPDF.domainLB.length).get(0)) *
													neighborPDF.eval(CMRF.splitArray(point, domainPDF.domainLB.length).get(1))),
											Double.MIN_VALUE))));
					mutInf += pairwiseEntropy.computeAreaUnderCurve();
				}
				cmrf.edgeWeights[i][j] = mutInf;
				cmrf.edgeWeights[j][i] = mutInf;
			}
		}
	}
	
	public void updateEdgeProbsTRBP(int iterCount) { 
		
		System.out.print("Updating edge probabilities...");
		TRBPMinSpanningTree mst = new TRBPMinSpanningTree();
		
		double[][] vec = mst.getMinSpanningTree(cmrf.edgeWeights);
		double stepSize = 2.0/(iterCount+20); // from Hunter
		
		for (int i=0; i<cmrf.numNodes; i++) { 
			for (int j=0; j<cmrf.numNodes; j++) { 
				if (i==j) { continue; }
				cmrf.edgeProbs[i][j] = stepSize * vec[i][j] + (1-stepSize) * cmrf.edgeProbs[i][j];
			}
		}
		
		System.out.println("done.");
		
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
		for (CMRFNode sender : cmrf.nodes) { 
			HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>> senderMap = 
					new HashMap<>();
			for (CMRFNode receiver : cmrf.nodes) {  
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

		for (int i=0; i<cmrf.numNodes; i++) {  // i is the sender
			for (int j=0; j<cmrf.numNodes; j++) { // j is the receiver
				System.out.print(i+"-"+j+" ");
				if (i==j) { continue; }
				CMRFNode sender = cmrf.nodes[i];
				CMRFNode receiver = cmrf.nodes[j];
				final double edgeProb = cmrf.edgeProbs[i][j];

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
								for (CMRFNode n : cmrf.nodes) { if (!n.equals(sender) && !n.equals(receiver)) { parents.add(n); } }

								double[] parentProbs = new double[parents.size()];
								RKHSFunction[] parentMsgs = new RKHSFunction[parents.size()];

								for (int k=0; k<parentMsgs.length; k++) {
									CMRFNode parent = parents.get(k);
									parentMsgs[k] = parent.outMessages.get(sender).get(senDom);
									parentProbs[k] = cmrf.edgeProbs[i][cmrf.getIndexInArray(parent, cmrf.nodes)];
								}
								ToDoubleFunction<double[]> numFunc = (point) -> (cmrf.getProdOfFuncPowers(parentMsgs, parentProbs, point));

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
														getModifiedExponentialFunction(sender, receiver, senDom, recDom, xs).eval(Xt) 
														*  numFunc.applyAsDouble(Xt)/denomFunc.applyAsDouble(Xt)
														)).computeIntegral()));
								// goddamn that is beautiful 
								// we get a function for each sender domain, and sum over
								senderFuncs[cmrf.getIndexInArray(senDom, sender.domains)] = updateFunc;
					}

					receiverDomainFuncs[cmrf.getIndexInArray(recDom, receiver.domains)] = 
							new RKHSFunction(
									recDom.k,
									recDom.domainLB,
									recDom.domainUB,
									(point) -> (cmrf.sumOverMessages(point, senderFuncs)));
				}

				// normalize the messages 
				double normalizingConstant = 0.0;
				for (RKHSFunction func : receiverDomainFuncs) { 
					normalizingConstant += func.computeIntegral();
				}
				final double Q = normalizingConstant;
				for (int k=0; k< receiverDomainFuncs.length; k++) {
					if (Q != 0) { 
						RKHSFunction oldFunc = receiverDomainFuncs[k];
						receiverDomainFuncs[k] = new RKHSFunction(
								oldFunc.k,
								oldFunc.domainLB,
								oldFunc.domainUB,
								(point) -> (oldFunc.eval(point)/Q));
					} else {  // integral is zero, so 
						RKHSFunction oldFunc = receiverDomainFuncs[k];
						receiverDomainFuncs[k] = new RKHSFunction(
								oldFunc.k,
								oldFunc.domainLB,
								oldFunc.domainUB,
								(point) -> (1.0/oldFunc.computeDomainVolume()));
					}
				}

				// dump it all to the buffer hashmap
				HashMap<CMRFNodeDomain, RKHSFunction> funcMap = sender.outMessages.get(receiver);
				for (CMRFNodeDomain rD : receiver.domains) { 
					RKHSFunction rDFunc = receiverDomainFuncs[cmrf.getIndexInArray(rD, receiver.domains)];
					funcMap.put(rD, rDFunc);	
				}
				messageMaps.get(sender).put(receiver, funcMap);
			}
		}

		// now update the messages from the buffer
		for (CMRFNode sender : cmrf.nodes) { 
			for (CMRFNode receiver : cmrf.nodes) { 
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
		for (CMRFNode v : cmrf.nodes) { 
			int recNodeIndex = cmrf.getIndexInArray(v, cmrf.nodes);
			double nodeEnthalpy = 0.0; 

			for (CMRFNodeDomain d : v.domains) { 
				// compute single-node domain enthalpy 
				RKHSFunction probabilityFunc = v.pseudomarginals.get(d);
				RKHSFunction enthalpyFunc = new RKHSFunction(
						d.k,
						d.domainLB,
						d.domainUB,
						(point) -> (
								cmrf.functionFloor(
										probabilityFunc.eval(point) * d.energyRKHS.eval(point))));
				double domainEnthalpy = enthalpyFunc.computeIntegral();
				if (Double.isNaN(domainEnthalpy)) { throw new RuntimeException("NaN enthalpy"); }
				nodeEnthalpy += domainEnthalpy;

				for (CMRFNode neighbor : cmrf.nodes) { 
					if (neighbor.equals(v)) { continue; }
					for (CMRFNodeDomain nd : neighbor.domains) { 
						// get the pdf for the neighbor's domain 
						int nRecNodeInd = cmrf.getIndexInArray(neighbor, cmrf.nodes);
						CMRFEdge edge = cmrf.edges[nRecNodeInd][recNodeIndex];
						CMRFEdgeDomain edgeDomain = edge.getEdgeDomain(d, nd);

						RKHSFunction pairwiseProbFunc = edge.pseudomarginals.get(edgeDomain);
						RKHSFunction pairwiseEnergyFunc = edgeDomain.eFuncRKHS;

						// compute enthalpy, add it to single node enthalpy 
						RKHSFunction pairwiseEnthalpyFunc = new RKHSFunction(
								pairwiseProbFunc.k,
								pairwiseProbFunc.domainLB,
								pairwiseProbFunc.domainUB,
								(point) -> (
										cmrf.functionFloor(
												pairwiseProbFunc.eval(point) * pairwiseEnergyFunc.eval(point))));
						double pairwiseEnthalpy = pairwiseEnthalpyFunc.computeIntegral();
						if (Double.isNaN(pairwiseEnthalpy)) { throw new RuntimeException("NaN enthalpy"); }
						nodeEnthalpy += 0.5 * cmrf.edgeProbs[recNodeIndex][nRecNodeInd] * pairwiseEnthalpy;
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

		for (CMRFNode node: cmrf.nodes) { 
			double nodeEntropy = 0.0;

			for (CMRFNode neighbor : cmrf.nodes) {
				double edgeEntropy = 0.0;
				double mutualInf = 0.0;
				if (node.equals(neighbor)) { continue; }
				int nodeInd = cmrf.getIndexInArray(node, cmrf.nodes);
				int neighborInd = cmrf.getIndexInArray(neighbor, cmrf.nodes);
				CMRFEdge edge = cmrf.edges[nodeInd][neighborInd];


				for (CMRFNodeDomain domain : node.domains) { 
					// compute intra-node domain entropy
					RKHSFunction domainPDF = node.pseudomarginals.get(domain);
					RKHSFunction domainEntropyFunc = new RKHSFunction(
							domainPDF.k,
							domainPDF.domainLB,
							domainPDF.domainUB,
							(point)->(this.safeEntropyCalcSingleton(domainPDF, point)));
									
					double domainEntropy = domainEntropyFunc.computeAreaUnderCurve();
					if (Double.isNaN(domainEntropy)) { 
						Matrix m = domainPDF.dumpPoints();
						m.print(3, 5);
						m = domainEntropyFunc.dumpPoints();
						m.print(3, 5);
						throw new RuntimeException("NaN entropy"); 	
					}
					if (domainEntropy < 0) { 
						throw new RuntimeException("Negative entropy");
					}
					nodeEntropy += domainEntropy;

					// compute pairwise entropy
					for (CMRFNodeDomain neighborDomain : neighbor.domains) { 
						RKHSFunction neighborPDF = neighbor.pseudomarginals.get(neighborDomain);
						CMRFEdgeDomain edgeDomain = edge.getEdgeDomain(domain, neighborDomain);
						RKHSFunction pairwisePDF = edge.pseudomarginals.get(edgeDomain);

						RKHSFunction pairwiseEntropy = new RKHSFunction(
								edgeDomain.resAllK,
								edgeDomain.domainLB,
								edgeDomain.domainUB,
								(point) -> (this.safeEntropyCalcPairwise(
												pairwisePDF, 
												domainPDF, 
												neighborPDF, 
												point)));

						double pairEntropy = pairwiseEntropy.computeAreaUnderCurve();
						
						if (Double.isNaN(pairEntropy)) {
							Matrix m = pairwiseEntropy.dumpPoints();
							m.print(3, 5);
							throw new RuntimeException("NaN entropy");
						}

						edgeEntropy += pairEntropy * cmrf.edgeProbs[nodeInd][neighborInd];
						mutualInf += pairEntropy;
					}
					if (Double.isNaN(edgeEntropy)) { throw new RuntimeException("NaN entropy"); }

					cmrf.edgeWeights[nodeInd][neighborInd] = mutualInf;
					cmrf.edgeWeights[neighborInd][nodeInd] = mutualInf;

					nodeEntropy += 0.5 * edgeEntropy;
				}
			}
			if (Double.isNaN(nodeEntropy)) { throw new RuntimeException("NaN entropy"); }

			totalEntropy += nodeEntropy;

		}

		return totalEntropy;
	}
	
	/** returns the exponential function at a specific point for the TRBP message update 
	 * this is basically currying, but awkward as all hell because Java
	 * 
	 * @param sender
	 * @param receiver
	 * @param senDomain
	 * @param recDomain
	 * @param xs
	 * @return
	 */
	private RKHSFunction getModifiedExponentialFunction(
			CMRFNode sender, 
			CMRFNode receiver, 
			CMRFNodeDomain senDomain,
			CMRFNodeDomain recDomain,
			double[] xs) {
		int i = cmrf.getIndexInArray(sender, cmrf.nodes); 
		int j = cmrf.getIndexInArray(receiver, cmrf.nodes); 
		CMRFEdge edge = cmrf.edges[i][j];
		// god i wish i could write all this crap in python
		return new RKHSFunction(
				senDomain.k,
				senDomain.domainLB,
				senDomain.domainUB,
				(Xt) -> (Math.exp(
						edge.getEnergyAtPoint(
								CMRFEdgeDomain.concatArrays(Xt, xs))
						/cmrf.edgeProbs[i][j] +
						sender.getDomainForPoint(Xt).getEnergyAtPoint(Xt))));
	}

	/**
	 * Prints out the TRBP marginals to .dat files
	 */
	public void printMarginalsTRBP() {
		System.out.print("Printing TRBP pseudomarginals... ");
		if (!cmrf.nodesAdded) { return; } // do nothing if we don't have an MRF

		for (int i=0; i<cmrf.nodes.length; i++) {
			CMRFNode node = cmrf.nodes[i];
			for (int j=0; j<node.domains.length; j++) {

				// print unitary pseudomarginals
				CMRFNodeDomain domain = node.domains[j];
				try {
					String filename = "trbp_u_"+i+"-"+j+".dat";
					PrintWriter writer = new PrintWriter(filename, "UTF-8");
					Matrix m = node.pseudomarginals.get(domain).dumpPoints();
					m.print(writer, 10, 10);
					writer.flush();
				} catch(FileNotFoundException | UnsupportedEncodingException e) {
					throw new RuntimeException("PrintWriting failed for "+
							"node " + i +", domain " + j +".\n" +e.getMessage());
				}

				// print pairwise pseudomarginals
				for (int k=0; k<cmrf.nodes.length; k++) {
					if (i==k) { continue; }
					CMRFEdge edge = cmrf.edges[i][k];
					for (int l=0; l<edge.domainLinks.length; l++) {
						CMRFEdgeDomain edgeDomain = edge.domainLinks[l];
						int x1 = cmrf.getIndexInArray(edgeDomain.resOneDomain, node.domains);
						int x2 = cmrf.getIndexInArray(edgeDomain.resTwoDomain, cmrf.nodes[k].domains);

						try {
							String filename = "trbp_p_"+i+k+"_"+x1+x2+".dat";
							PrintWriter writer = new PrintWriter(filename, "UTF-8");
							Matrix m = edge.pseudomarginals.get(edgeDomain).dumpPoints();
							m.print(writer, 10, 10);
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
	
	/**
	 * Wraps entropy calculations in a "safe" way
	 * @param pdf
	 * @param point
	 * @return
	 */
	public double safeEntropyCalcSingleton(RKHSFunction pdf, double[] point) { 
		double prob = pdf.eval(point);
		if (prob < CMRF.logThreshold) { 
			return -1 * prob * (prob-1);
		} else { 
			return -1* prob * Math.log(prob);
		}
	}
	
	/**
	 * Wraps the pairwise TRBP entropy calculation in a "safe" way
	 * @param pairPDF
	 * @param nodePDF
	 * @param neighborPDF
	 * @param point
	 * @return
	 */
	public double safeEntropyCalcPairwise(
			RKHSFunction pairPDF,
			RKHSFunction nodePDF,
			RKHSFunction neighborPDF,
			double[] point) {
		ArrayList<double[]> indivPts = CMRF.splitArray(point, nodePDF.domainLB.length);
		
		double pairVal = pairPDF.eval(point);
		double nodeInt = new RKHSFunction(
				nodePDF.k,
				nodePDF.domainLB,
				nodePDF.domainUB,
				(nodePoint) -> (pairPDF.eval(CMRF.concatArrays(nodePoint, indivPts.get(1)))))
				.computeAreaUnderCurve();
		double neighborInt = new RKHSFunction(
				neighborPDF.k,
				neighborPDF.domainLB,
				neighborPDF.domainUB,
				(neighborPoint) -> (pairPDF.eval(CMRF.concatArrays(indivPts.get(0), neighborPoint))))
				.computeAreaUnderCurve();		
		double quotient = pairVal/(nodeInt * neighborInt);
		
		// approximate log(quotient) as quotient-1 for small quotient
		if (quotient < CMRF.logThreshold) { 
			return pairVal * (quotient - 1); 
		} else { 
			return pairVal * Math.log(quotient);
		}		
	}

}
