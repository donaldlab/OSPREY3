package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import Jama.Matrix;

public class SCMF {

	CMRF cmrf;
	
	public SCMF(CMRF cmrf) {
		this.cmrf = cmrf;
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

			double freeEnergy = enth - cmrf.constRT*entr;
			double logZ =  - freeEnergy/cmrf.constRT; // freeEnergy is a direct lower bound on logZ

			// break if the bound gets worse, i.e. we step over a local maximum
//			if (logZ < oldLogZ) { 
//				cmrf.ranSCMF = true;
//				System.out.println("Ran over a local maximum.");
//				System.out.println("DONE: logZLB: "+oldLogZ);
//				printMarginalsSCMF();
//				return oldLogZ;
//			}

			System.out.println("enth: "+enth+", entr: "+entr+", logZLB: "+logZ);

			// break if the other termination condition is reached
			if ((Math.abs(logZ-oldLogZ) <= cmrf.threshold) || (iter >= cmrf.maxIters)) {
				System.out.println("Converged.");
				cmrf.ranSCMF = true;
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
		if (! cmrf.nodesAdded) { 
			throw new RuntimeException("Can't initialize SCMF marginals for the"
					+ "continuous-label MRF unless nodes have already been added. ");
		}

		System.out.print("Initializing marginals... ");

		for (CMRFNode node : cmrf.nodes) { 
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
		for (CMRFNode node : cmrf.nodes) {

			double partFn = 0.0;
			RKHSFunction[] domainMarginals = new RKHSFunction[node.domains.length];

			for (CMRFNodeDomain domain : node.domains) {
				int domainIndex = Arrays.asList(node.domains).indexOf(domain);
				RKHSFunction oneBodyEnergy = domain.energyRKHS;

				// get a list of average interaction energy functions
				ArrayList<RKHSFunction> avgInteractionEnergies = new ArrayList<>();
				for (CMRFNode neighbor : cmrf.nodes) {
					int nodeInd = Arrays.asList(cmrf.nodes).indexOf(node);
					int neighborInd = Arrays.asList(cmrf.nodes).indexOf(node);
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
						(point)->cmrf.sumOverMessages(point, avgInteractionE));

				RKHSFunction logUnnormalizedBelief = new RKHSFunction(
						domain.k,
						domain.domainLB,
						domain.domainUB,
						(point) -> (-(oneBodyEnergy.eval(point) + meanFieldE.eval(point))/cmrf.constRT));

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
		for (CMRFNode node : cmrf.nodes) { 
			RKHSFunction[] domainMarginals = newBeliefs.get(node);
			for (int i=0; i<node.domains.length; i++) { 
				final CMRFNodeDomain domain = node.domains[i];
				final RKHSFunction oldM = node.marginals.get(domain);
				final RKHSFunction newM = domainMarginals[i];
				final double alpha = cmrf.lambda;
				node.marginals.put(domain, new RKHSFunction(
						domain.k,
						domain.domainLB,
						domain.domainUB,
						(point)->(alpha*newM.eval(point) + (1-alpha)*oldM.eval(point))));
			}
		}

		cmrf.lambda = cmrf.lambda * (1 - 1.0/cmrf.maxIters);
	}

	public double pairwiseExpectation(
			double[] point, 
			CMRFNode node,
			CMRFNodeDomain domain,
			CMRFNode neighbor,
			CMRFNodeDomain nDomain) { 
		int nodeIndex = Arrays.asList(cmrf.nodes).indexOf(node);
		int neighborIndex = Arrays.asList(cmrf.nodes).indexOf(node);

		CMRFEdge e = cmrf.edges[nodeIndex][neighborIndex];
		CMRFEdgeDomain ed = e.getEdgeDomain(domain, nDomain);

		RKHSFunction func = new RKHSFunction(
				nDomain.k,
				nDomain.domainLB,
				nDomain.domainUB,
				(nPoint)->(ed.eFuncRKHS.eval(CMRFEdgeDomain.concatArrays(point, nPoint))));

		return func.computeExpectation();
	}

	/**
	 * Computes the enthalpy when running SCMF -- note that there are no pairwise terms 
	 * @return 
	 */
	public double computeEnthalpySCMF() { 
		double totalEnthalpy = 0.0;
		// sum over nodes of p*E plus pariwise p*E
		for (CMRFNode v : cmrf.nodes) { 
			int recNodeIndex = cmrf.getIndexInArray(v, cmrf.nodes);
			double nodeEnthalpy = 0.0; 

			for (CMRFNodeDomain d : v.domains) { 
				// compute single-node domain enthalpy -- Ex~Q[\ln phi]
				RKHSFunction probabilityFunc = v.marginals.get(d);
				RKHSFunction enthalpyFunc = new RKHSFunction(
						d.k,
						d.domainLB,
						d.domainUB,
						(point) -> (
								cmrf.functionFloor(
										probabilityFunc.eval(point) * d.energyRKHS.eval(point))));
				double domEnth = enthalpyFunc.computeIntegral();
				//System.out.println("SCMF domEnth: " + domEnth);
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
		for (CMRFNode node : cmrf.nodes) { 
			double nodeEntropy = 0.0;
			for (CMRFNodeDomain domain : node.domains) { 
				RKHSFunction domainPDF = node.marginals.get(domain);
				RKHSFunction domainEntropyFunc = new RKHSFunction(
						domainPDF.k,
						domainPDF.domainLB,
						domainPDF.domainUB,
						(point)->(
								this.safeEntropyCalcSingleton(domainPDF, point)));
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
	 * Wraps entropy calculations in a "safe" way
	 * @param pdf
	 * @param point
	 * @return
	 */
	public double safeEntropyCalcSingleton(RKHSFunction pdf, double[] point) { 
		double prob = pdf.eval(point);
		if (prob < CMRF.logThreshold) { 
			return -1 * prob;
		} else { 
			return -1* prob * Math.log(prob);
		}
	}
	
	
	/**
	 * Prints out SCMF marginals to .dat files
	 */
	public void printMarginalsSCMF() { 
		System.out.print("Printing SCMF marginals... ");
		if (! cmrf.nodesAdded) { return; } // do nothing if we don't have an MRF 
		for (int i=0; i<cmrf.nodes.length; i++) { 
			CMRFNode node = cmrf.nodes[i];
			for (int j=0; j<node.domains.length; j++) { 
				CMRFNodeDomain domain = node.domains[j];
				String filename = "scmf_"+i+"-"+j+".dat";
				try { 
					PrintWriter writer = new PrintWriter(filename, "UTF-8");
					Matrix m = node.marginals.get(domain).dumpPoints();
					m.print(writer, 10, 10);
					writer.flush();
				} catch(FileNotFoundException | UnsupportedEncodingException e) { 
					throw new RuntimeException("PrintWriting failed for "+ 
							"node " + i +", domain " + j +".\n" +e.getMessage());
				}
			}
		}
		System.out.println("done.");
	}


	
	
}
