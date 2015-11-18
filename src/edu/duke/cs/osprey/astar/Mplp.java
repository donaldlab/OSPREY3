package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.tools.CreateMatrix;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * MPLP: computes a bound on the best energy of a partial conformation.
 *
 * @author pablo
 *
 */
public class Mplp {

    //number of residues under consideration
    private int numResidues;

    // counts of rots per position
    int rotsPerPos[];

    // rotamers unpruned for each residue (given by the assigned AA type)
    ArrayList<ArrayList<Integer>> unprunedRotsPerPos;

    //the offset in the array index for each level
    private int nodeIndexOffset[] = null;

    //the reduced min pairwise energy matrix
    private double[][][][] unifiedMinEnergyMatrix = null;

    private boolean interactionGraph[][];

    Mplp(int aNumResidues, ArrayList<ArrayList<Integer>> unprunedRotsPerPos, EnergyMatrix aPairwiseMinEnergyMatrix) {
        numResidues = aNumResidues;
        this.unprunedRotsPerPos = unprunedRotsPerPos;

        // A count of rots per position useful for the initialization of matrices\
        rotsPerPos = new int[numResidues];
        for (int res = 0; res < numResidues; res++) {
            rotsPerPos[res] = this.unprunedRotsPerPos.get(res).size();
        }

        interactionGraph = computeResidueResidueInteractionGraph(aPairwiseMinEnergyMatrix);
        unifiedMinEnergyMatrix = mergeIntraAndPairMats(aPairwiseMinEnergyMatrix, interactionGraph);
    }

    // Computes the low-energy bound using the EMPLP algorithm
    // availableRots is a "numRes*rots" matrix that tells us which 
    // rotamers at each position are available.  Only those rotamers
    // will be considered for the calculation.
    public double optimizeEMPLP(int aPartialConf[], int iterations) {

        // The partial conf contains references to rotamers that were pruned, while our pairwise matrix removed those rotamers. Thus, we must creat
        //  a new partial conf that maps to our pruned matrix.
        int mappedPartialConf[] = new int[numResidues];
        for (int res = 0; res < numResidues; res++) {
            if (aPartialConf[res] == -1) {
                mappedPartialConf[res] = -1;
            } else {
                mappedPartialConf[res] = this.unprunedRotsPerPos.get(res).indexOf(aPartialConf[res]);
            }
        }

        // availableRots is just a list of lists with the list of rotamers available for the calculation at each rotamer; it is more convenient
        // 		than partialConf.  This code might be unnecessary but it makes the algorithm more elegant, IMO
        int availableRots[][] = new int[numResidues][];
        for (int res = 0; res < numResidues; res++) {
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
        // lambda is the message matrix for EMPLP; we set all messages to zero
        double lambda[][][] = CreateMatrix.create3DMsgMat(numResidues, rotsPerPos, 0.0f);

        // The belief on each rotamer.
        double belief[][] = CreateMatrix.create2DRotMatrix(numResidues, rotsPerPos, 0.0f);

        // Ebound is the fscore.
        double Ebound = Double.NEGATIVE_INFINITY;
        // To find out when we should stop, we use a delta factor to see the rate of change.  For now, our rate of change is fixed at 0.01kcal.
        double minRateOfChange = 0.01;

        // EMPLP algorithm.  
        // The algorithm should loop until convergence.  For now we do it I times. 
        // Complexity of EMPLP: O(I*numRes*numRes*rotsPerRes*rotsPerRes*numRes)
        for (int i = 0; i < iterations; i++) {
            for (int resI = 0; resI < numResidues; resI++) {
                for (int resJ = resI + 1; resJ < numResidues; resJ++) {
                    // Ignore residue pairs that are too far away to interact.
                    if (interactionGraph[resI][resJ]) {
                        // We first update lambda[resJ][resI][rotIR] and immediately after we update lambda[resI][resJ][rotJS]
                        for (int rotIR : availableRots[resI]) {
                            belief[resI][rotIR] -= lambda[resJ][resI][rotIR];

                            ArrayList<Double> msgsFromRotsAtJ_to_rotIR = new ArrayList<Double>();
                            for (int rotJS : availableRots[resJ]) {
                                msgsFromRotsAtJ_to_rotIR.add(belief[resJ][rotJS] - lambda[resI][resJ][rotJS] + unifiedMinEnergyMatrix[resI][rotIR][resJ][rotJS]);
                            }
                            if (availableRots[resJ].length == 0) {
                                System.out.println("NO ROTS MPLP CRASHING");
                            }
                            lambda[resJ][resI][rotIR] = -0.5 * belief[resI][rotIR] + 0.5 * Collections.min(msgsFromRotsAtJ_to_rotIR);
                            belief[resI][rotIR] += lambda[resJ][resI][rotIR];
                        }
                        // Now we update lambda[resI][resJ][rotJS]
                        for (int rotJS : availableRots[resJ]) {
                            belief[resJ][rotJS] -= lambda[resI][resJ][rotJS];

                            ArrayList<Double> msgsFromRotsAtI_to_rotJS = new ArrayList<Double>();
                            for (int rotIR : availableRots[resI]) {
                                msgsFromRotsAtI_to_rotJS.add(belief[resI][rotIR] - lambda[resJ][resI][rotIR] + unifiedMinEnergyMatrix[resI][rotIR][resJ][rotJS]);
                            }

                            lambda[resI][resJ][rotJS] = -0.5 * belief[resJ][rotJS] + 0.5 * Collections.min(msgsFromRotsAtI_to_rotJS);
                            belief[resJ][rotJS] += lambda[resI][resJ][rotJS];
                        }
                    }
                }
            }
            double oldEbound = Ebound;
            Ebound = 0.0f;
            for (int resI = 0; resI < numResidues; resI++) {
                // It is important to ignore other rotamers at positions that were already assigned when computing the min.  
                Ebound += computeMinBeliefInReducedAvailableRots(belief[resI], availableRots[resI]);
            }
            // If we already converged, then we can exit. 
            if ((Ebound - oldEbound) < minRateOfChange) {
                break;
            }

        }
        return Ebound;
    }
    // Compute the minimum value in an array where only a subset of the rotamers are present.
    // Belief is an array for a residue, and availableRots is the indexes in the belief array for which to compute the search.	

    double computeMinBeliefInReducedAvailableRots(double belief[], int availableRots[]) {
        double minValue = Double.POSITIVE_INFINITY;
        for (int rotIR : availableRots) {
            if (belief[rotIR] < minValue) {
                minValue = belief[rotIR];
            }
        }
        return minValue;
    }

    /**
     * The original dual derivation of MPLP did not consider intra-energies.
     * Thus, I found that the easiest way to solve this would be to unify the
     * intra and pair energies by dividing the intra interaction between n-1 and
     * adding it to the edge. Note that the unified emat is 4D whereas the
     * standard matrices used by A* in OSPREY are 6D.
     *
     * @param numRes
     * @param rotsPerRes
     * @param twoDTo3D
     * @param emat
     * @param interactionGraph
     * @return
     */
    private double[][][][] mergeIntraAndPairMats(EnergyMatrix emat, boolean interactionGraph[][]) {
        double unifiedEmat[][][][] = CreateMatrix.create4DRotMatrix(numResidues, rotsPerPos, 0.0f);

        //Build Neighborhood
        int[] numNeighbors = new int[numResidues];
        for (int xres = 0; xres < numResidues; xres++) {
            numNeighbors[xres] = 0;
            for (int yres = 0; yres < numResidues; yres++) {
                if (interactionGraph[xres][yres] && xres != yres) {
                    numNeighbors[xres] += 1;
                }
            }
        }

        for (int resI = 0; resI < numResidues; resI++) {
            for (int resJ = resI + 1; resJ < numResidues; resJ++) {
                if (interactionGraph[resI][resJ]) {//If neighbors(resI,resJ)	
                    for (int rotIR_Ix = 0; rotIR_Ix < this.unprunedRotsPerPos.get(resI).size(); rotIR_Ix++) {
                        int rotIR_origMat = this.unprunedRotsPerPos.get(resI).get(rotIR_Ix);
                        for (int rotJS_Ix = 0; rotJS_Ix < this.unprunedRotsPerPos.get(resJ).size(); rotJS_Ix++) {
                            int rotJS_origMat = this.unprunedRotsPerPos.get(resJ).get(rotJS_Ix);
                            unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix] = emat.getPairwise(resI, rotIR_origMat, resJ, rotJS_origMat);
                            unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix] += emat.getOneBody(resI, rotIR_origMat) / ((double) numNeighbors[resI]);
                            unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix] += emat.getOneBody(resJ, rotJS_origMat) / ((double) numNeighbors[resJ]);
                            unifiedEmat[resJ][rotJS_Ix][resI][rotIR_Ix] = unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix];
                        }
                    }
                }
            }
        }
        return unifiedEmat;
    }

    // Computes an interaction graph with an energy cutoff; returns a pairs matrix between all pairs of residues. This is 
    private boolean[][] computeResidueResidueInteractionGraph(EnergyMatrix emat) {
        int countTrueInteraction = 0;
        int countNoInteraction = 0;
        double INTERACTION_CUTOFF = 0.0;
        boolean[][] interactionGraph = new boolean[numResidues][numResidues];
        for (int xres = 0; xres < numResidues; xres++) {
            for (int yres = 0; yres < numResidues; yres++) {
                interactionGraph[xres][yres] = false;
            }
        }
        for (int xres = 0; xres < numResidues; xres++) {
            for (int yres = xres + 1; yres < numResidues; yres++) {
                double maxInteraction_x_y = 0.0f;
                for (int xrot : this.unprunedRotsPerPos.get(xres)) {
                    for (int yrot : this.unprunedRotsPerPos.get(yres)) {
                        double pairInteraction = emat.getPairwise(xres, xrot, yres, yrot);
                        if (Math.abs(pairInteraction) > maxInteraction_x_y) {
                            maxInteraction_x_y = Math.abs(pairInteraction);
                        }
                    }
                }
                if (maxInteraction_x_y > INTERACTION_CUTOFF) {
                    interactionGraph[xres][yres] = true;
                    interactionGraph[yres][xres] = true;
                    countTrueInteraction += 1;
                } else {
                    countNoInteraction += 1;
                }
            }
        }
        System.out.println("Using an interaction cutoff of " + INTERACTION_CUTOFF);
        System.out.println("Interaction between " + countTrueInteraction + " pairs");
        System.out.println("No interaction between " + countNoInteraction + " pairs");
        return interactionGraph;

    }

    // Computes an interaction graph in which residues only interact if they are on different strands
    // or if the interactions are between the non-mutable strand 
    public void setCrossTermInteractionGraph(List<EnergyMatrix> boundResNumToUnboundEmat, List<Integer> boundResNumToUnboundResNum,
            List<Boolean> boundresNumToIsMutableStrand, boolean[][] belongToSameStrand) {
        boolean[][] interactionGraph = new boolean[numResidues][numResidues];
        int numInteractions = 0;
        for (int xres = 0; xres < numResidues; xres++) {
            for (int yres = 0; yres < numResidues; yres++) {
                interactionGraph[xres][yres] = false;
            }
        }
        for (int xres = 0; xres < numResidues; xres++) {
            for (int yres = xres + 1; yres < numResidues; yres++) {
                //If on the same strand (check if it is non-mutable)
                if (belongToSameStrand[xres][yres]) {
                    //xres and yres are not on the mutable strand, so create interaction
                    if (!boundresNumToIsMutableStrand.get(xres)) {
                        interactionGraph[xres][yres] = true;
                        interactionGraph[yres][xres] = true;
                        numInteractions++;
                    }
                } //else, they are not on the same strand so create interaction
                else {
                    interactionGraph[xres][yres] = true;
                    interactionGraph[yres][xres] = true;
                    numInteractions++;
                }
            }
        }
        System.out.println("MPLP: Cross-Term Interaction Graph Num Interaction: " + numInteractions);
        this.interactionGraph = interactionGraph;
    }

    public void createOnlyPairwiseMat(EnergyMatrix ematBound, List<EnergyMatrix> boundResNumToUnboundEmat, List<Integer> boundResNumToUnboundResNum,
            List<Boolean> boundresNumToIsMutableStrand, boolean[][] belongToSameStrand) {
        double unifiedEmat[][][][] = CreateMatrix.create4DRotMatrix(numResidues, rotsPerPos, 0.0f);

        //Build Neighborhood
        int[] numNeighbors = new int[numResidues];
        for (int xres = 0; xres < numResidues; xres++) {
            numNeighbors[xres] = 0;
            for (int yres = 0; yres < numResidues; yres++) {
                if (interactionGraph[xres][yres] && xres != yres) {
                    numNeighbors[xres] += 1;
                }
            }
        }

        for (int resI = 0; resI < numResidues; resI++) {
            for (int resJ = resI + 1; resJ < numResidues; resJ++) {
                //If they are interacting could be because: 1) They are on the same, non-mutable strand
                //                                          2) They are on opposite strands
                if (interactionGraph[resI][resJ]) {//If neighbors(resI,resJ)	
                    for (int rotIR_Ix = 0; rotIR_Ix < this.unprunedRotsPerPos.get(resI).size(); rotIR_Ix++) {
                        int rotIR_origMat = this.unprunedRotsPerPos.get(resI).get(rotIR_Ix);
                        for (int rotJS_Ix = 0; rotJS_Ix < this.unprunedRotsPerPos.get(resJ).size(); rotJS_Ix++) {
                            int rotJS_origMat = this.unprunedRotsPerPos.get(resJ).get(rotJS_Ix);
                            //If the interaction is not pairwise, then we will use only  the pairwise energies
                            //And the intraE that .getOneBody(resI, rotIR_origMat) + ematBound.getOneBody(resJ, rotJS_origMat);
                            EnergyMatrix ematUnbound1 = boundResNumToUnboundEmat.get(resI);
                            EnergyMatrix ematUnbound2 = boundResNumToUnboundEmat.get(resJ);
                            int resI_unbound = boundResNumToUnboundResNum.get(resI);
                            int resJ_unbound = boundResNumToUnboundResNum.get(resJ);
                            double pairwiseE = ematBound.getPairwise(resI, rotIR_origMat, resJ, rotJS_origMat);
                            double intraE_I = ematBound.getOneBody(resI, rotIR_origMat) / ((double) numNeighbors[resI]);
                            double intraE_J = ematBound.getOneBody(resJ, rotJS_origMat) / ((double) numNeighbors[resJ]);
                            double intraE_I_unbound = ematUnbound1.getOneBody(resI_unbound, rotIR_origMat) / ((double) numNeighbors[resI]);
                            double intraE_J_unbound = ematUnbound2.getOneBody(resJ_unbound, rotJS_origMat) / ((double) numNeighbors[resJ]);
                            double totalE = pairwiseE + intraE_I + intraE_J - intraE_I_unbound - intraE_J_unbound;
                            unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix] = totalE;
                            unifiedEmat[resJ][rotJS_Ix][resI][rotIR_Ix] = unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix];
                        }
                    }
                }
            }
        }
        this.unifiedMinEnergyMatrix = unifiedEmat;
    }

    public void createOnlyPairwiseMatWithLigand(EnergyMatrix ematBound, List<EnergyMatrix> boundResNumToUnboundEmat, List<Integer> boundResNumToUnboundResNum,
            List<Boolean> boundresNumToIsMutableStrand, boolean[][] belongToSameStrand) {
        double unifiedEmat[][][][] = CreateMatrix.create4DRotMatrix(numResidues, rotsPerPos, 0.0f);

        //Build Neighborhood
        int[] numNeighbors = new int[numResidues];
        //also keep track of neighbors on opposite strand
        int[] numNeighborsOppStrand = new int[numResidues];
        for (int xres = 0; xres < numResidues; xres++) {
            numNeighbors[xres] = 0;
            numNeighborsOppStrand[xres] = 0;
            for (int yres = 0; yres < numResidues; yres++) {
                if (interactionGraph[xres][yres] && xres != yres) {
                    numNeighbors[xres] += 1;
                    if (!belongToSameStrand[xres][yres]) {
                        numNeighborsOppStrand[xres] += 1;
                    }
                }
            }
        }

        for (int resI = 0; resI < numResidues; resI++) {
            for (int resJ = resI + 1; resJ < numResidues; resJ++) {
                //If they are interacting could be because: 1) They are on the same, non-mutable strand
                //                                          2) They are on opposite strands
                if (interactionGraph[resI][resJ]) {//If neighbors(resI,resJ)	
                    for (int rotIR_Ix = 0; rotIR_Ix < this.unprunedRotsPerPos.get(resI).size(); rotIR_Ix++) {
                        int rotIR_origMat = this.unprunedRotsPerPos.get(resI).get(rotIR_Ix);
                        for (int rotJS_Ix = 0; rotJS_Ix < this.unprunedRotsPerPos.get(resJ).size(); rotJS_Ix++) {
                            int rotJS_origMat = this.unprunedRotsPerPos.get(resJ).get(rotJS_Ix);
                            //If they belong to the same strand and interacting: they are both non-mutable
                            if (belongToSameStrand[resI][resJ]) {
                                //For the non-mutable non-mutable interactions we want to use the bound energies
                                //We will subtract out the unbound GMEC energy using the pre-computed const-term
                                //of the objective function for COMETS
                                double pairwiseE = ematBound.getPairwise(resI, rotIR_origMat, resJ, rotJS_origMat);
                                double intraE_I = ematBound.getOneBody(resI, rotIR_origMat) / ((double) numNeighbors[resI]);
                                double intraE_J = ematBound.getOneBody(resJ, rotJS_origMat) / ((double) numNeighbors[resJ]);
                                double totalE = pairwiseE + intraE_I + intraE_J;
                                unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix] = totalE;
                                unifiedEmat[resJ][rotJS_Ix][resI][rotIR_Ix] = unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix];
                            } else {
                                //If the interaction is not pairwise, then we will use only  the pairwise energies
                                //And the intraE that .getOneBody(resI, rotIR_origMat) + ematBound.getOneBody(resJ, rotJS_origMat);
                                EnergyMatrix ematUnbound1 = boundResNumToUnboundEmat.get(resI);
                                EnergyMatrix ematUnbound2 = boundResNumToUnboundEmat.get(resJ);
                                int resI_unbound = boundResNumToUnboundResNum.get(resI);
                                int resJ_unbound = boundResNumToUnboundResNum.get(resJ);
                                double pairwiseE = ematBound.getPairwise(resI, rotIR_origMat, resJ, rotJS_origMat);
                                double intraE_I = ematBound.getOneBody(resI, rotIR_origMat) / ((double) numNeighbors[resI]);
                                double intraE_J = ematBound.getOneBody(resJ, rotJS_origMat) / ((double) numNeighbors[resJ]);
                                double intraE_I_unbound = ematUnbound1.getOneBody(resI_unbound, rotIR_origMat) / ((double) numNeighborsOppStrand[resI]);
                                double intraE_J_unbound = ematUnbound2.getOneBody(resJ_unbound, rotJS_origMat) / ((double) numNeighborsOppStrand[resJ]);
                                double totalE = pairwiseE + intraE_I + intraE_J - intraE_I_unbound - intraE_J_unbound;
                                unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix] = totalE;
                                unifiedEmat[resJ][rotJS_Ix][resI][rotIR_Ix] = unifiedEmat[resI][rotIR_Ix][resJ][rotJS_Ix];
                            }
                        }
                    }
                }
            }
        }
        this.unifiedMinEnergyMatrix = unifiedEmat;
    }
}
