package edu.duke.cs.osprey.astar.conf.scoring;

import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class MPLPPairwiseHScorer implements AStarScorer {
	
	private static class BetaVars {
		
		private RCs rcs;
		
		// one beta var (decision variable, beta_ij(xi,xj))
		// for every node,rotamer,node,rotamer tuple, where (node,rotamer) pairs are ordered
		private double[][] pairs;
		
		public BetaVars(RCs rcs) {
			
			this.rcs = rcs;
			
			// allocate space for the vars
			int numPos = rcs.getNumPos();
			pairs = new double[numPos*numPos][];
			for (int i=0; i<numPos; i++) {
				for (int j=0; j<numPos; j++) {
					pairs[getNodePairIndex(i, j)] = new double[rcs.get(i).length*rcs.get(j).length];
				}
			}
		}
		
		public double get(int pos1, int rci1, int pos2, int rci2) {
			return pairs[getNodePairIndex(pos1, pos2)][getRCPairIndex(pos1, rci1, rci2)];
		}
		
		public void set(int pos1, int rci1, int pos2, int rci2, double val) {
			pairs[getNodePairIndex(pos1, pos2)][getRCPairIndex(pos1, rci1, rci2)] = val;
		}
		
		private int getNodePairIndex(int pos1, int pos2) {
			return pos1*rcs.getNumPos() + pos2;
		}
		
		private int getRCPairIndex(int pos1, int rci1, int rci2) {
			return rci2*rcs.get(pos1).length + rci1;
		}
	}
	
	private static class LambdaVars {
		
		private RCs rcs;
		
		// two edge messages (decision variables, lambda_ji(xi) and lambda_ij(xj))
		// for every node,node,rotamer tuple, where nodes are ordered
		private double[][] vars;
		
		public LambdaVars(RCs rcs) {
			
			this.rcs = rcs;
			
			// allocate space for the messages
			int numPos = rcs.getNumPos();
			vars = new double[numPos*numPos][];
			for (int i=0; i<numPos; i++) {
				for (int j=0; j<numPos; j++) {
					
					if (i == j) {
						continue;
					}
					
					vars[getNodePairIndex(i, j)] = new double[rcs.get(j).length];
				}
			}
		}
		
		public double get(int pos1, int pos2, int rci2) {
			return vars[getNodePairIndex(pos1, pos2)][rci2];
		}
		
		public void set(int pos1, int pos2, int rci2, double val) {
			vars[getNodePairIndex(pos1, pos2)][rci2] = val;
		}
		
		
		public double getSum(ConfIndex confIndex, int pos1, int rci1) {
			
			assert (confIndex.getNumUndefined() == rcs.getNumPos());
			
			double sum = 0;
			for (int k=0; k<confIndex.getNumUndefined(); k++) {
				int pos3 = confIndex.getUndefinedPos()[k];
				
				if (pos3 == pos1) {
					continue;
				}
				
				sum += get(pos3, pos1, rci1);
			}
			return sum;
		}
		
		public double getSumWithout(ConfIndex confIndex, int pos1, int rci1, int pos2) {
			return getSum(confIndex, pos1, rci1) - get(pos2, pos1, rci1);
		}

		private int getNodePairIndex(int pos1, int pos2) {
			assert (pos1 != pos2);
			return pos1*rcs.getNumPos() + pos2;
		}
	}
	
	private static class GammaVars {
		
		private RCs rcs;
		
		// one node message (decision variable, gamma_ij(xi,xj))
		// for every node,node,rotamer tuple, where nodes are ordered
		private double[][] vars;
		
		public GammaVars(RCs rcs) {
			
			this.rcs = rcs;
			
			// allocate space for the messages
			int numPos = rcs.getNumPos();
			vars = new double[numPos*numPos][];
			for (int i=0; i<numPos; i++) {
				for (int j=0; j<numPos; j++) {
					
					if (i == j) {
						continue;
					}
					
					vars[getNodePairIndex(i, j)] = new double[rcs.get(j).length];
				}
			}
		}
		
		public double get(int pos1, int pos2, int rci2) {
			return vars[getNodePairIndex(pos1, pos2)][rci2];
		}
		
		public void set(int pos1, int pos2, int rci2, double val) {
			vars[getNodePairIndex(pos1, pos2)][rci2] = val;
		}
		
		private int getNodePairIndex(int pos1, int pos2) {
			return pos1*rcs.getNumPos() + pos2;
		}
	}
	
	private EnergyMatrix emat;
	private BetaVars betas;
	private LambdaVars lambdas;
	private GammaVars gammas;
	
	public MPLPPairwiseHScorer(EnergyMatrix emat, RCs rcs) {
		this.emat = emat;
		this.betas = new BetaVars(rcs);
		this.lambdas = new LambdaVars(rcs);
		this.gammas = new GammaVars(rcs);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
		// run the traditional A* heuristic to get the initial conformation
		
		// init the undefined positions
		int[] undefinedRCis = new int[confIndex.getNumUndefined()];
		Arrays.fill(undefinedRCis, -1);
		
		double tradEnergy = 0;
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			double minPos1Energy = Double.POSITIVE_INFINITY;
			int minPos1RCi = -1;
			
			for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
				int rc1 = rcs.get(pos1)[rci1];
			
				// undefined single energy
				double pos1Energy = emat.getOneBody(pos1, rc1);
				
				// TEMP
				assert (confIndex.getNumDefined() == 0);
				/* undefined-defined pairwise energy
				for (int j=0; j<confIndex.getNumDefined(); j++) {
					int pos2 = confIndex.getDefinedPos()[j];
					int rc2 = confIndex.getDefinedRCs()[j];
					if (pos2 < pos1) {
						pos1Energy += emat.getPairwise(pos1, rc1, pos2, rc2);
					}
				}
				*/
				
				// undefined-undefined pairwise energy
				for (int j=0; j<confIndex.getNumUndefined(); j++) {
					int pos2 = confIndex.getUndefinedPos()[j];
					
					if (pos2 >= pos1) {
						continue;
					}
					
					double minPos2Energy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {
						minPos2Energy = Math.min(minPos2Energy, emat.getPairwise(pos1, rc1, pos2, rc2));
					}
					
					pos1Energy += minPos2Energy;
				}
				
				if (pos1Energy < minPos1Energy) {
					minPos1Energy = pos1Energy;
					minPos1RCi = rci1;
				}
			}
			
			// save the conf
			undefinedRCis[i] = minPos1RCi;
			
			tradEnergy += minPos1Energy;
		}
		
		// TEMP
		System.out.println(String.format("traditional H score check: %f", tradEnergy));
		System.out.println(String.format("traditional conf energy: %f", calcConfEnergy(rcs, confIndex, undefinedRCis)));
		
		// compute the beta vars
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos1 == pos2) {
					continue;
				}
					
				for (int rci1=0; rci1<rcs.get(i).length; rci1++) {
					for (int rci2=0; rci2<rcs.get(j).length; rci2++) {
						double theta = calcBeta(rcs, confIndex, pos1, rci1, pos2, rci2);
						betas.set(pos1, rci1, pos2, rci2, theta);
					}
				}
			}
		}
		
		// compute the dual LP objective function value
		double obj = 0;
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			double min1 = Double.POSITIVE_INFINITY;
			for (int rci1=0; rci1<rcs.get(i).length; rci1++) {
			
				double val1 = 0;
				
				for (int j=0; j<confIndex.getNumUndefined(); j++) {
					int pos2 = confIndex.getUndefinedPos()[j];
					
					double min2 = Double.POSITIVE_INFINITY;
					for (int rci2=0; rci2<rcs.get(j).length; rci2++) {
						min2 = Math.min(min2, betas.get(pos1, rci1, pos2, rci2));
					}
				
					val1 += min2;
				}
				
				if (val1 < min1) {
					min1 = val1;
				}
			}
			
			obj += min1;
		}
		
		System.out.println(String.format("dual obj: %f", obj));
		
		// NOTE: initializing lambda vars from the conf speeds convergence
		// on my DAGK test case by about 2 iterations of EMPLP
		// ALSO: we MUST initialize the vars for early stopping to be sound
		
		// compute the lambda vars
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos1 == pos2) {
					continue;
				}
				
				// lambda_ji(xi)
				for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
					lambdas.set(pos2, pos1, rci1, calcLambda(rcs, pos2, pos1, rci1));
				}
				
				// lambda_ij(xj)
				for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
					lambdas.set(pos1, pos2, rci2, calcLambda(rcs, pos1, pos2, rci2));
				}
			}
		}
		//
		
		// check the initial lambda energy
		double lambdaEnergy = calcLambdaEnergy(rcs, confIndex, undefinedRCis);
		
		// do edge-based MPLP
		for (int i=0; i<500; i++) {
			lambdas = stepEMPLP(rcs, confIndex);
			double newLambdaEnergy = calcLambdaEnergy(rcs, confIndex, undefinedRCis);
			double diff = Math.abs(newLambdaEnergy - lambdaEnergy);
			if (diff < 0.001) {
				System.out.println(String.format("converged in %d steps", i));
				break;
			}
			lambdaEnergy = newLambdaEnergy;
		}
		
		// compute the gamma vars
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos1 == pos2) {
					continue;
				}
				
				for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
					gammas.set(pos1, pos2, rci2, calcGamma(rcs, confIndex, pos1, pos2, rci2));
				}
			}
		}
		//
		
		// check initial gamma energy
		double gammaEnergy = calcGammaEnergy(rcs, confIndex, undefinedRCis);
		
		// do node-based MPLP
		for (int i=0; i<500; i++) {
			gammas = stepNMPLP(rcs, confIndex);
			double newGammaEnergy = calcGammaEnergy(rcs, confIndex, undefinedRCis);
			double diff = Math.abs(newGammaEnergy - gammaEnergy);
			if (diff < 0.001) {
				System.out.println(String.format("converged in %d steps", i));
				break;
			}
			gammaEnergy = newGammaEnergy;
		}
		
		return 0;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		throw new Error("not implemented");
	}
	
	private double calcTheta(RCs rcs, ConfIndex confIndex, int pos1, int rci1, int pos2, int rci2) {
		// theta_ij = single(i)/|N(i)| + single(j)/|N(j)| + pair(i,j)
		int rc1 = rcs.get(pos1)[rci1];
		int rc2 = rcs.get(pos2)[rci2];
		return emat.getOneBody(pos1, rc1)/(confIndex.getNumUndefined() - 1)
			+ emat.getOneBody(pos2, rc2)/(confIndex.getNumUndefined() - 1)
			+ emat.getPairwise(pos1, rc1, pos2, rc2);
	}
	
	private double calcBeta(RCs rcs, ConfIndex confIndex, int pos1, int rci1, int pos2, int rci2) {
		
		// calc beta based on splitting theta evenly across the two directions
		double theta = calcTheta(rcs, confIndex, pos1, rci1, pos2, rci2);
		return theta/2;
		
		/* match the traditional A* heuristic
		int rc1 = rcs.get(pos1)[rci1];
		int rc2 = rcs.get(pos2)[rci2];
		double e = emat.getOneBody(pos1, rc1)/(confIndex.getNumUndefined() - 1);
		if (pos2 < pos1) {
			e += emat.getPairwise(pos1, rc1, pos2, rc2);
		}
		return e;
		*/
	}
	
	private double calcLambda(RCs rcs, int pos3, int pos1, int rci1) {
		// lambda_ki(xi) = max_xk beta_ki(xk, xi)
		// lambda_{pos3,pos1}(rci1) = max_{rci3} beta_{pos3,pos1}(rci3, rci1)
		double minVal = Double.POSITIVE_INFINITY;
		for (int rci3=0; rci3<rcs.get(pos3).length; rci3++) {
			minVal = Math.min(minVal, betas.get(pos3, rci3, pos1, rci1));
		}
		return minVal;
	}
	
	private double calcLambdaMinus(RCs rcs, ConfIndex confIndex, int pos1, int rci1, int pos2) {
		// lambda_i-j(xi) = sum_{k in N(i) without j} lambda_ki(xi)
		// lambda_{pos1-pos2}(rci1) = sum_{pos3 !in pos1,pos2} lambda_{pos3,pos1}(rci1)
		double sum = 0;
		for (int k=0; k<confIndex.getNumUndefined(); k++) {
			int pos3 = confIndex.getUndefinedPos()[k];
			
			if (pos3 == pos1 || pos3 == pos2) {
				continue;
			}
			
			sum += calcLambda(rcs, pos3, pos1, rci1);
		}
		return sum;
	}
	
	private double calcGamma(RCs rcs, ConfIndex confIndex, int pos2, int pos1, int rci1) {
		// gamma_ji(xi) = max_xj [ theta_ij(xi,xj) + lambda_j-i(xj) ]
		// gamma_{pos2,pos1}(rci1) = max_rci2 [ theta_{pos1,pos2}(rci1, rci2) + lambda_{pos2-pos1}(rci2) ]
		double minVal = Double.POSITIVE_INFINITY;
		for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
			double theta = calcTheta(rcs, confIndex, pos1, rci1, pos2, rci2);
			double lambda = calcLambdaMinus(rcs, confIndex, pos2, rci2, pos1);
			minVal = Math.min(minVal, theta + lambda);
		}
		return minVal;
	}
	
	private double calcLambdaEnergy(RCs rcs, ConfIndex confIndex, int[] undefinedRCis) {
		
		int[] mplpRCis = new int[undefinedRCis.length];
		Arrays.fill(mplpRCis, -1);
		
		// compute energies from the lambda vars
		// b_i(xi) = sum_k lambda_ki(xi)
		// energy_pos1(rci1) = sum_{pos3 != pos1} lambda_{pos3,pos1}(rci1)
		double lambdaEnergy = 0;
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			double minPosEnergy = Double.POSITIVE_INFINITY;
			int minRCi1 = -1;
			
			for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {

				double posEnergy = 0;
				for (int k=0; k<confIndex.getNumUndefined(); k++) {
					int pos3 = confIndex.getUndefinedPos()[k];
					if (pos3 != pos1) {
						posEnergy += lambdas.get(pos3, pos1, rci1);
					}
				}
				
				if (posEnergy < minPosEnergy) {
					minPosEnergy = posEnergy;
					minRCi1 = rci1;
				}
			}
			
			lambdaEnergy += minPosEnergy;
			mplpRCis[i] = minRCi1;

			/* TEMP
			System.out.println(String.format("pos %3d: tradRC=%3d, mplpRC=%3d, lambda=%f",
				pos1,
				rcs.get(pos1)[undefinedRCis[pos1]],
				mplpRCis[i],
				minPosEnergy
			));
			*/
		}
		
		System.out.println(String.format("lambda energy: %f, confEnergy: %f",
			lambdaEnergy,
			calcConfEnergy(rcs, confIndex, mplpRCis)
		));
		
		return lambdaEnergy;
	}
	
	private double calcGammaEnergy(RCs rcs, ConfIndex confIndex, int[] undefinedRCis) {
		
		int[] mplpRCis = new int[undefinedRCis.length];
		Arrays.fill(mplpRCis, -1);
		
		// compute energies from the gamma vars
		// b_i(xi) = sum_{k in N(i)} gamma_ki(xi)
		// energy_pos1(rci1) = sum_{pos3 != pos1} gamma_{pos3,pos1}(rci1)
		double gammaEnergy = 0;
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			double minPosEnergy = Double.POSITIVE_INFINITY;
			int minRCi1 = -1;
			
			for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {

				double posEnergy = 0;
				for (int k=0; k<confIndex.getNumUndefined(); k++) {
					int pos3 = confIndex.getUndefinedPos()[k];
					if (pos3 != pos1) {
						posEnergy += gammas.get(pos3, pos1, rci1);
					}
				}
				
				if (posEnergy < minPosEnergy) {
					minPosEnergy = posEnergy;
					minRCi1 = rci1;
				}
			}
			
			gammaEnergy += minPosEnergy;
			mplpRCis[i] = minRCi1;

			/* TEMP
			System.out.println(String.format("pos %3d: tradRC=%3d, mplpRC=%3d, gamma=%f",
				pos1,
				rcs.get(pos1)[undefinedRCis[pos1]],
				mplpRCis[i],
				minPosEnergy
			));
			*/
		}
		
		gammaEnergy /= confIndex.getNumUndefined();
		
		System.out.println(String.format("gamma energy: %f, confEnergy: %f",
			gammaEnergy,
			calcConfEnergy(rcs, confIndex, mplpRCis)
		));
		
		return gammaEnergy;
	}
	
	private double calcConfEnergy(RCs rcs, ConfIndex confIndex, int[] rcis) {
		double confEnergy = 0;
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			int rc1 = rcs.get(pos1)[rcis[i]];
			
			confEnergy += emat.getOneBody(pos1, rc1);
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				int rc2 = rcs.get(pos2)[rcis[j]];
				
				if (pos2 < pos1) {
					confEnergy += emat.getPairwise(pos1, rc1, pos2, rc2);
				}
			}
		}
		return confEnergy;
	}
	
	private LambdaVars stepEMPLP(RCs rcs, ConfIndex confIndex) {
		
		// lambda_ji(xi) = -0.5*lambda_{i-j}(xi) + 0.5*max_xj [ lamda_{j-i}(xj) + theta_ij(xi,xj) ]
		// and with i,j reversed
		
		// lambda_{pos2,pos1}(rci1) = [
		//                              max_{rci2} [ lambda_{pos2-pos1}(rci2) + theta_{pos1,pos2)(rci1,rci2) ]
		//                              - lamba_{pos1-pos2}(rci1)
		//                            ]/2
		// and with pos1,pos2 reversed
		
		// time complexity
		// O(n*n*2(m*(m*n + n)))
		//  = O(2n^2*(nm^2 + nm))
		//  = O(2n^3*m^2 + 2n^3*m)
		
		// NOTE: don't copy-on-write to the lambda vars!
		// the algorithm is designed for immediate updates to the lambda vars
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos2 >= pos1) {
					continue;
				}
				
				for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
					
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
						double energy = lambdas.getSumWithout(confIndex, pos2, rci2, pos1)
							+ calcTheta(rcs, confIndex, pos1, rci1, pos2, rci2);
						minEnergy = Math.min(minEnergy, energy);
					}
					
					minEnergy = (minEnergy - lambdas.getSumWithout(confIndex, pos1, rci1, pos2))/2;
					
					lambdas.set(pos2, pos1, rci1, minEnergy);
				}
				
				for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
					
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
						double energy = lambdas.getSumWithout(confIndex, pos1, rci1, pos2)
							+ calcTheta(rcs, confIndex, pos2, rci2, pos1, rci1);
						minEnergy = Math.min(minEnergy, energy);
					}
					
					minEnergy = (minEnergy - lambdas.getSumWithout(confIndex, pos2, rci2, pos1))/2;
					
					lambdas.set(pos1, pos2, rci2, minEnergy);
				}
			}
		}
		
		return lambdas;
	}
	
	private GammaVars stepNMPLP(RCs rcs, ConfIndex confIndex) {
		
		// gamma_ij(xj) = max_xi [
		//                  theta_ij(xi,xj) - gamma_ji(xi)
		//                  + 2/(|N(i)|+1)*sum_{k in N(i)} gamma_ki(xi)
		//                ]
		// for all i and j in N(i)
		// gamma_{pos1,pos2}(rci2) = max_rci1 [
		//                             theta_{pos1,pos2}(rci1,rci2) - gamma_{pos2,pos1}(rci1)
		//                             + 2/numUndefined*sum_{pos3 != pos1} gamma_{pos3,pos1}(rci1)
		//                           ]
		// for all pos1 and pos2 != pos1
		
		// time complexity
		// O(n*n*m*m*n)
		//  = O(n^3*m^2)
		
		// NOTE: need to copy-on-write so vars get updated correctly
		// ie, return the new set of vars instead of overwriting existing ones
		// the immediate updates version actually does converge, but not monotonically,
		// so it breaks the guarantees that make early stopping sound
		GammaVars newGammas = new GammaVars(rcs);
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos2 == pos1) {
					continue;
				}
				
				for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
					
					double minVal = Double.POSITIVE_INFINITY;
					
					for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
						double theta = calcTheta(rcs, confIndex, pos1, rci1, pos2, rci2);
						double gamma = gammas.get(pos2, pos1, rci1);
						
						double sum = 0;
						for (int k=0; k<confIndex.getNumUndefined(); k++) {
							int pos3 = confIndex.getUndefinedPos()[k];
							if (pos3 != pos1) {
								sum += gammas.get(pos3, pos1, rci1);
							}
						}
						
						double val = theta - gamma + 2*sum/confIndex.getNumUndefined();
						
						minVal = Math.min(minVal, val);
					}
					
					newGammas.set(pos1, pos2, rci2, minVal);
				}
			}
		}
		
		return newGammas;
	}
}
