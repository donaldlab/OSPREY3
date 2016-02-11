/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.CreateMatrix;
import java.util.ArrayList;
import java.util.Collections;

/**
 * MPLP: Computes a bound on the best energy of a partial conformation Original
 * implementation: "Fixing Max Product: Convergent Message Passing Algorithms
 * for MAP LP-Relaxation" This implementation: "Tightening LP Relaxations for
 * MAP using Message Passing" The advantage of this implementation is that it
 * explicitly considers intra and pairwise energies
 *
 * @author hmn5, pablo (original author)
 *
 */
public class Mplp {

    //number of positions/nodes under consideration
    private int numPos;

    // counts of rots per positions
    int rotsPerPos[];

    // rotamers unpruned for each position
    ArrayList<ArrayList<Integer>> unprunedRotsPerPos;

    //For pairs pruning
    PruningMatrix pruneMat;

    EnergyMatrix emat;

    private boolean interactionGraph[][];

    int[] feasibleSolution;

    public Mplp(int aNumPos, EnergyMatrix aEmat, PruningMatrix aPruneMat) {
        init(aNumPos, aEmat, aPruneMat, 0.0);
    }

    public Mplp(int aNumPos, EnergyMatrix aEmat, PruningMatrix aPruneMat, double aEnergyCutOff) {
        init(aNumPos, aEmat, aPruneMat, aEnergyCutOff);
    }

    private void init(int aNumPos, EnergyMatrix aEmat, PruningMatrix aPruneMat, double aEnergyCutOff) {
        this.numPos = aNumPos;
        this.emat = aEmat;
        this.pruneMat = aPruneMat;

        this.unprunedRotsPerPos = new ArrayList<>();
        this.rotsPerPos = new int[numPos];
        for (int pos = 0; pos < numPos; pos++) {
            rotsPerPos[pos] = pruneMat.unprunedRCsAtPos(pos).size();
            unprunedRotsPerPos.add(pruneMat.unprunedRCsAtPos(pos));
        }

        this.interactionGraph = computePosPosInteractionGraph(emat, aEnergyCutOff);
    }

    /**
     * Returns a lower bound on the GMEC energy of a partial conformation Note
     * this implementation is based on Sontag et. al. "Tightening LP
     * Relaxations..." This implementation does not use higher-order clusters so
     * we only have edge to node messages Therefore the 2/3 and 1/3 in the above
     * paper are replaced by 1/2 and 1/2 However, it can be extended to allow
     * for triplet clusters and thus, to tighten the LP relaxation
     *
     * @param aPartialConf
     * @param iterations
     * @return
     */
    public double optimizeMPLP(int aPartialConf[], int iterations) {
        // The partial conf contains references to rotamers that were pruned, while our pairwise matrix removed those rotamers. Thus, we must creat
        //  a new partial conf that maps to our pruned matrix.
        int mappedPartialConf[] = new int[numPos];
        feasibleSolution = aPartialConf.clone();
        for (int pos = 0; pos < numPos; pos++) {
            if (aPartialConf[pos] == -1) {
                mappedPartialConf[pos] = -1;
            } else {
                mappedPartialConf[pos] = this.unprunedRotsPerPos.get(pos).indexOf(aPartialConf[pos]);
            }
        }

        // availableRots is just a list of lists with the list of rotamers available for the calculation at each rotamer; it is more convenient
        // 		than partialConf.  This code might be unnecessary but it makes the algorithm more elegant, IMO
        int availableRots[][] = new int[numPos][];
        for (int res = 0; res < numPos; res++) {
            if (mappedPartialConf[res] == -1) {
                availableRots[res] = new int[rotsPerPos[res]];
                for (int rot = 0; rot < rotsPerPos[res]; rot++) {
                    availableRots[res][rot] = rot;
                }
            } else { // residue has a rotamer already assigned in the partialConf.
                availableRots[res] = new int[1];
                availableRots[res][0] = mappedPartialConf[res];
            }
        }

        /**
         * Here comes the actual algorithm *
         */
        // lambda is the message matrix for MPLP; we set all messages to zero
        double lambda[][][] = CreateMatrix.create3DMsgMat(numPos, rotsPerPos, 0.0f);

        // the belief on each rotamer
        //b_i(x_i) = intrE(x_i) + sum_{k in N(i)} lambda_{k,i -> i}(x_i)
        double belief[][] = CreateMatrix.create2DRotMatrix(numPos, rotsPerPos, 0.0f);
        //initialize beliefs to intra energy term
        for (int pos = 0; pos < numPos; pos++) {
            for (int rot = 0; rot < this.rotsPerPos[pos]; rot++) {
                //need index of rot into original emat
                int rot_Emat = this.unprunedRotsPerPos.get(pos).get(rot);
                belief[pos][rot] += this.emat.getOneBody(pos, rot_Emat);
            }
        }

        // Ebound is the fscore
        double Ebound = Double.NEGATIVE_INFINITY;
        //To find out when we should stop, we use a delta factor to see the rate of change.
        //For now, our rate of change is fixed at 0.01kcal.
        double minRateOfChange = 0.01;

        //MPLP Algorithm
        for (int i = 0; i < iterations; i++) {
            for (int posI = 0; posI < numPos; posI++) {
                for (int posJ = posI + 1; posJ < numPos; posJ++) {
                    //Option to ignore pairs that are too far away to interact
                    if (interactionGraph[posI][posJ]) {
                        // We first update lambda[resJ][resI][rotIR] and immediately after we update lambda[resI][resJ][rotJS]
                        for (int rotIR : availableRots[posI]) {
                            //For edge to node update in "Tightening LP Relaxations..." (Figure 1)
                            //lambda_{i}^{-j}(x_i) + theta_i(x_i) = b_i(x_i) - lambda_{j,i -> i}(x_i)
                            //This is computationally more efficient
                            belief[posI][rotIR] -= lambda[posJ][posI][rotIR];// now b_i(x_i) = lambda_{i}^{-j}(x_i) + theta_i(x_i)

                            ArrayList<Double> msgsFromRotsAtJ_to_rotIR = new ArrayList<Double>();
                            for (int rotJS : availableRots[posJ]) {
                                //Create a tuple to check if pair is pruned
                                msgsFromRotsAtJ_to_rotIR.add(belief[posJ][rotJS] - lambda[posI][posJ][rotJS] + emat.getPairwise(posJ, unprunedRotsPerPos.get(posJ).get(rotJS), posI, unprunedRotsPerPos.get(posI).get(rotIR)));
                            }
                            if (availableRots[posJ].length == 0) {
                                System.out.println("NO ROTS MPLP CRASHING");
                                return Double.POSITIVE_INFINITY;
                            }
                            lambda[posJ][posI][rotIR] = -0.5 * belief[posI][rotIR] + 0.5 * Collections.min(msgsFromRotsAtJ_to_rotIR);
                            belief[posI][rotIR] += lambda[posJ][posI][rotIR];
                        }
                        //Now we update lambda[posI][posJ][rotJS]
                        for (int rotJS : availableRots[posJ]) {
                            belief[posJ][rotJS] -= lambda[posI][posJ][rotJS]; // now b_j(x_j) = lambda_(j}^{-i}(x_j) + theta_j(x_j)

                            ArrayList<Double> msgFromRotsAtI_to_rotJS = new ArrayList<>();
                            for (int rotIR : availableRots[posI]) {
                                msgFromRotsAtI_to_rotJS.add(belief[posI][rotIR] - lambda[posJ][posI][rotIR] + emat.getPairwise(posI, unprunedRotsPerPos.get(posI).get(rotIR), posJ, unprunedRotsPerPos.get(posJ).get(rotJS)));
                            }
                            if (availableRots[posI].length == 0) {
                                System.out.println("NO ROTS MPLP CRASHING");
                                return Double.POSITIVE_INFINITY;
                            }
                            lambda[posI][posJ][rotJS] = -0.5 * belief[posJ][rotJS] + 0.5 * Collections.min(msgFromRotsAtI_to_rotJS);
                            belief[posJ][rotJS] += lambda[posI][posJ][rotJS];
                        }
                    }
                }
            }
            double oldEbound = Ebound;
            Ebound = 0.0f;
            for (int posI = 0; posI < numPos; posI++) {
                // We must ignore other rotamers at position that were already assigned when computing the min.
                Ebound += computeMinBeliefInReducedAvailableRots(belief[posI], availableRots[posI]);
            }
            // If we already converged, then we can exit
            if (Math.abs(Ebound - oldEbound) < minRateOfChange) {
                break;
            }
        }

        getFeasibleSolution(belief, availableRots);
        return Ebound;
    }

    void getFeasibleSolution(double[][] belief, int[][] availableRots) {
        for (int pos = 0; pos < this.numPos; pos++) {
            if (feasibleSolution[pos] == -1) {
                //unprunedRotNum corresponds to to the rotamer number in the reduced unpruned set
                int unprunedRotNum = computeMinBelief(belief[pos], availableRots[pos]);
                //rotNum corresponds to the actual rot num (i.e. index into ematrix)
                int rotNum = this.pruneMat.unprunedRCsAtPos(pos).get(unprunedRotNum);
                feasibleSolution[pos] = rotNum;
            }
        }
    }

    int computeMinBelief(double[] belief, int[] availableRots) {
        double minValue = Double.POSITIVE_INFINITY;
        int minRot = -1;
        for (int rotIR : availableRots) {
            double score = belief[rotIR];
            if (score < minValue) {
                minValue = score;
                minRot = rotIR;
            }
        }
        return minRot;
    }

    double computeMinBeliefInReducedAvailableRots(double belief[], int availableRots[]) {
        double minValue = Double.POSITIVE_INFINITY;
        for (int rotIR : availableRots) {
            double score = belief[rotIR];
            if (score < minValue) {
                minValue = score;
            }
        }
        return minValue;
    }

    //computes an interaction graph between pairs of positions with an energy cutoff;
    private boolean[][] computePosPosInteractionGraph(EnergyMatrix emat, double energyCutOff) {
        int countTrueInteraction = 0;
        int countNoInteraction = 0;
        boolean[][] interactionGraph = new boolean[numPos][numPos];
        for (int xpos = 0; xpos < numPos; xpos++) {
            for (int ypos = 0; ypos < numPos; ypos++) {
                interactionGraph[xpos][ypos] = false;
            }
        }
        for (int xpos = 0; xpos < numPos; xpos++) {
            for (int ypos = xpos + 1; ypos < numPos; ypos++) {
                double maxInteraction_x_y = 0.0f;
                for (int xrot : this.unprunedRotsPerPos.get(xpos)) {
                    for (int yrot : this.unprunedRotsPerPos.get(ypos)) {
                        double pairInteraction = emat.getPairwise(xpos, xrot, ypos, yrot);
                        if (Math.abs(pairInteraction) > maxInteraction_x_y) {
                            maxInteraction_x_y = Math.abs(pairInteraction);
                        }
                    }
                }
                if (maxInteraction_x_y > energyCutOff) {
                    interactionGraph[xpos][ypos] = true;
                    interactionGraph[ypos][xpos] = true;
                    countTrueInteraction += 1;
                } else {
                    countNoInteraction += 1;
                }
            }
        }
        if (energyCutOff > 0.01) {
            System.out.println("Using an interaction cutoff of " + energyCutOff);
            System.out.println("Interaction between " + countTrueInteraction + " pairs");
            System.out.println("No interaction between " + countNoInteraction + " pairs");
        }
        return interactionGraph;
    }

    public void setInteractionGraph(boolean[][] newInteractionGraph) {
        if (newInteractionGraph.length != this.numPos) {
            throw new RuntimeException("ERROR: Cannot set interaction graph since its length != num positions");
        } else {
            this.interactionGraph = newInteractionGraph;
        }
    }

    public void setEmat(EnergyMatrix aEmat) {
        if (emat.oneBody.size() != numPos) {
            throw new RuntimeException("ERROR: Cannot set emat in MPLP because the size of one body terms != numPos");
        }
        this.emat = aEmat;
    }

    public void setEmatRecomputeInteractionGraph(EnergyMatrix aEmat, double eCutOff) {
        this.emat = aEmat;
        this.interactionGraph = computePosPosInteractionGraph(emat, eCutOff);
    }
}
