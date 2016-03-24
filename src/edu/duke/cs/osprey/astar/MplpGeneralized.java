/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.CreateMatrix;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import org.apache.commons.lang.ArrayUtils;

/**
 * This is a generalized MPLP algorithm to handle factors that are larger than pairs
 * It has not yet been rigorously tested, but I believe it works based on preliminary
 * tests.
 * @author hmn5
 */
public class MplpGeneralized {

    //number of positions/nodes under consideration
    private int numPos;

    // counts of rots per positions
    int rotsPerPos[];

    // rotamers unpruned for each position
    ArrayList<ArrayList<Integer>> unprunedRotsPerPos;

    //For pairs pruning
    PruningMatrix pruneMat;

    EnergyMatrix emat;

    public int[] feasibleSolution;

    double[][] belief;
    double[][][] delta;

    HashSet<Factor> factorSet;

    //We will keep a sorted list of the factors (sorted by size and then posNums)
    //in the future we may want to update our factors in a particular order
    //For now, I figured the pairwise factors should be updated first
    ArrayList<Factor> factorList = new ArrayList<>();

    int numFactors;

    public MplpGeneralized(int aNumPos, EnergyMatrix aEmat, PruningMatrix aPruneMat) {
        init(aNumPos, aEmat, aPruneMat);
    }

    private void init(int aNumPos, EnergyMatrix aEmat, PruningMatrix aPruneMat) {
        this.numPos = aNumPos;
        this.emat = aEmat;
        this.pruneMat = aPruneMat;

        this.unprunedRotsPerPos = new ArrayList<>();
        this.rotsPerPos = new int[numPos];
        for (int pos = 0; pos < numPos; pos++) {
            rotsPerPos[pos] = pruneMat.unprunedRCsAtPos(pos).size();
            unprunedRotsPerPos.add(pruneMat.unprunedRCsAtPos(pos));
        }
        this.factorSet = getFactors();
        for (Factor f : factorSet) {
            factorList.add(f);
        }
        Collections.sort(factorList, compare);
        this.numFactors = this.factorList.size();

        int[] partialConf = new int[this.numPos];
        Arrays.fill(partialConf, -1);
//        optimizeMPLP(partialConf, 1);
    }

    /**
     * Returns a lower bound on the GMEC energy of a partial conformation Note
     * this implementation is based on Sontag et. al. "Introduction to Dual
     * Decomposition for Inference" This implementation allows for factors of
     * arbitrary size, but only computes the "first-order" LP-relaxation (ie.
     * factors sum to the individual marginals) Higher-order relaxations can be
     * done by extending the method as in Sontag et. al. "Tightening LP
     * Relaxations..."
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
        // delta is the message matrix for MPLP; we set all messages to zero
        this.delta = get3DMsgMat(rotsPerPos, 0.0);

        // the belief on each rotamer
        //b_i(x_i) = intraE(x_i) + sum_{k in N(i)} lambda_{k,i -> i}(x_i)
        this.belief = CreateMatrix.create2DRotMatrix(numPos, rotsPerPos, 0.0f);
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
            //For every factor...
            for (int fNum = 0; fNum < this.numFactors; fNum++) {
                Factor f = this.factorList.get(fNum);
                //For every pos in factor...
                for (int posNum = 0; posNum < f.size; posNum++) {
                    int pos = f.posNums[posNum];
                    //For every rot in pos...
                    for (int rot : availableRots[pos]) {
                        //Now we update delta_{f,i}(x_i)

                        //delta_{i}^{-f} = intraE(x_i) + \sum_{fhat != f} delta_{fhat,i}(x_i)
                        //               = belief_{i}(x_i) - delta_{f,i}(x_i)
                        this.belief[pos][rot] -= this.delta[fNum][posNum][rot];

                        double minMarginal = Double.POSITIVE_INFINITY;
                        //TODO: If the number of configurations within a factor that we are computing a max marginal is large,
                        //we should use an A* variant to compute the max marginal
                        int[][] configsInFactor = getConfigurationsWithinFactor(f, pos, rot, availableRots);
                        for (int[] config : configsInFactor) {
                            double marginal = 0.0;
                            double configE = getFactorEnergy(config);
                            marginal += configE;
                            for (int posInFactor : f.posNums) {
                                int posIndex = ArrayUtils.indexOf(f.posNums, posInFactor);
                                int rotInFactor = config[posInFactor];
                                double deltaInverseFAtPosRot;
                                if (posInFactor != pos) {
                                    deltaInverseFAtPosRot = belief[posInFactor][rotInFactor] - this.delta[fNum][posIndex][rot];
                                } else {
                                    deltaInverseFAtPosRot = this.belief[pos][rot];
                                }
                                marginal += deltaInverseFAtPosRot;
                            }
                            minMarginal = Math.min(minMarginal, marginal);
                        }
                        //Update delta
                        this.delta[fNum][posNum][rot] = -this.belief[pos][rot] + (1.0 / ((double) f.size)) * minMarginal;
                        //Update belief with new delta
                        this.belief[pos][rot] += this.delta[fNum][posNum][rot];
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

    /**
     * Returns a list of configurations within a particular factor
     *
     * @param f
     * @param posFixed
     * @param rotFixed
     * @return
     */
    int[][] getConfigurationsWithinFactor(Factor f, int posFixed, int rotFixed, int[][] availableRotsPerPos) {
        int[][] configsInFactor = getConfigurationWithinFactorHelper(f.posNums[0], f, posFixed, rotFixed, availableRotsPerPos);
        return configsInFactor;
    }

    int[][] getConfigurationWithinFactorHelper(int posToGet, Factor f, int posFixed, int rotFixed, int[][] availableRotsPerPos) {
        int numPosInFactor = f.posNums.length;
        int[][] configsInFactor;

        int[] rcAtPos;
        if (posFixed == posToGet) {
            rcAtPos = new int[1];
            rcAtPos[0] = rotFixed;
        } else {
            rcAtPos = availableRotsPerPos[posToGet];
        }

        if (posToGet == f.posNums[numPosInFactor - 1]) {
            configsInFactor = new int[rcAtPos.length][this.numPos];
            for (int i = 0; i < rcAtPos.length; i++) {
                configsInFactor[i] = new int[this.numPos];
                Arrays.fill(configsInFactor[i], -1);
                configsInFactor[i][posToGet] = rcAtPos[i];
            }
        } else {
            int index = ArrayUtils.indexOf(f.posNums, posToGet);
            int[][] configsHelper = getConfigurationWithinFactorHelper(f.posNums[index + 1], f, posFixed, rotFixed, availableRotsPerPos);

            configsInFactor = new int[configsHelper.length * rcAtPos.length][this.numPos];
            int iter = 0;
            for (int i = 0; i < rcAtPos.length; i++) {
                for (int[] config : configsHelper) {
                    int[] newConfig = Arrays.copyOf(config, config.length);
                    newConfig[posToGet] = rcAtPos[i];
                    configsInFactor[iter] = newConfig;
                    iter++;
                }
            }
        }
        return configsInFactor;
    }

    /**
     * Gets the energy of a particular configuration within a factor
     * (theta_{f}(x_f)
     *
     * @param mappedPartialConf this is an array where each element is the rot
     * index, not the rot num
     * @return
     */
    double getFactorEnergy(int[] mappedPartialConf) {
        ArrayList<Integer> posNums = new ArrayList<>();
        ArrayList<Integer> rcAtPos = new ArrayList<>();
        for (int pos = 0; pos < this.numPos; pos++) {
            if (mappedPartialConf[pos] >= 0) {
                posNums.add(pos);
                rcAtPos.add(mappedPartialConf[pos]);
            }
        }
        int pos1 = posNums.get(0);
        int rc1 = rcAtPos.get(0);
        int pos2 = posNums.get(1);
        int rc2 = rcAtPos.get(1);

        double factorE;
        if (posNums.size() == 2) {
            factorE = emat.getPairwise(pos1, rc1, pos2, rc2);
        } else {//higher order factor
            HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1, rc1, pos2, rc2);
            if (htf == null) {
                factorE = 0.0;
            } else {
                ArrayList<Integer> remainingPosNums = new ArrayList<>();
                ArrayList<Integer> remainingRCAtPos = new ArrayList<>();
                for (int i = 2; i < posNums.size(); i++) {
                    remainingPosNums.add(posNums.get(i));
                    remainingRCAtPos.add(rcAtPos.get(i));
                }
                factorE = getFactorEnergyHigherOrder(htf, remainingRCAtPos, remainingRCAtPos);
            }
        }
        return factorE;
    }

    double getFactorEnergyHigherOrder(HigherTupleFinder<Double> htf, ArrayList<Integer> remainingPos, ArrayList<Integer> remainingRC) {
        if (remainingPos.size() == 1) {
            return htf.getInteraction(remainingPos.get(0), remainingRC.get(0));
        } else {
            HigherTupleFinder<Double> htf2 = htf.getHigherInteractions(remainingPos.get(0), remainingRC.get(0));
            if (htf2 == null) {
                return 0.0;
            } else {
                ArrayList<Integer> remainingPos2 = new ArrayList<>();
                ArrayList<Integer> remainingRC2 = new ArrayList<>();
                for (int i = 1; i < remainingPos.size(); i++) {
                    remainingPos2.add(remainingPos.get(i));
                    remainingRC2.add(remainingRC.get(i));
                }
                return getFactorEnergyHigherOrder(htf2, remainingPos2, remainingRC2);
            }
        }
    }

    /**
     * Get all the (non-unary) factors
     *
     * @return
     */
    HashSet<Factor> getFactors() {
        HashSet<Factor> setOfFactors = new HashSet<>();
        for (int posI = 0; posI < this.numPos; posI++) {
            for (int posJ = 0; posJ < posI; posJ++) {
                //create pairwise factor
                setOfFactors.add(new Factor(posJ, posI));
                HashSet<Factor> hoFactorSet = getHigherOrderFactors(posJ, posI);
                for (Factor f : hoFactorSet) {
                    setOfFactors.add(f);
                }
            }
        }
        return setOfFactors;
    }

    /**
     * HMN: This is a horribly written method to get higher order factors I will
     * try to make up for it by documenting my logic
     *
     * @param posI pos number (posI less than posJ)
     * @param posJ pos number (posJ greater than posI)
     * @return
     */
    HashSet<Factor> getHigherOrderFactors(int posI, int posJ) {
        //Make sure that posI < posJ
        if (posI > posJ) {
            throw new RuntimeException("ERROR: posI cannot be greater than posJ");
        }
        //Initialize a set 
        HashSet<Factor> hoFactorSet = new HashSet<>();
        //Iterate over RCs
        for (int rcI : this.pruneMat.unprunedRCsAtPos(posI)) {
            for (int rcJ : this.pruneMat.unprunedRCsAtPos(posJ)) {
                //get higher order tuple if it exists
                HigherTupleFinder<Double> htf = this.emat.getHigherOrderTerms(posI, rcI, posJ, rcJ);
                if (htf != null) {
                    //For interacting positions in this higher order tuple
                    //If the position comes before posI (and therefor posJ too)
                    ArrayList<Integer> interactingPos = htf.getInteractingPos();
                    for (int i = 0; i < interactingPos.size(); i++) {
                        int posK = interactingPos.get(i);
                        if (posK < posI) {
                            //Then we iterate over RCs
                            for (int rcK : this.pruneMat.unprunedRCsAtPos(posK)) {
                                //if there is some non-zero interaction, then these three positoins
                                //have a higher order interaction
                                if (htf.getInteraction(posK, rcK) != 0.0) {
                                    //So we add the factor to our set
                                    Factor hoFactor = new Factor(posI, posJ, posK);
                                    hoFactorSet.add(hoFactor);
                                }
                                //Now we must add all higher order interactions involving these
                                //three residues
                                HigherTupleFinder<Double> htf2 = htf.getHigherInteractions(posK, rcK);
                                if (htf2 != null) {
                                    int[] posNumsHO = {posK, posI, posJ};
                                    Arrays.sort(posNumsHO);
                                    //we recurse on this method to get all higher order factors 
                                    HashSet<Factor> hoFactorSet2 = getHigherOrderFactors(htf2, posNumsHO);
                                    //we add these to our set
                                    for (Factor hoFactor : hoFactorSet2) {
                                        hoFactorSet.add(hoFactor);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return hoFactorSet;
    }

    /**
     * HMN: This method gets all higher order factors that interact with the
     * list of positions given in posNumsHO So if the list contains [2,3,4] We
     * want all higher order factors, if they exist, such as [1,2,3,4] The
     * position that we add to this factor must be less than all other positions
     * This helps with speed (since the HashSet would simply remove the
     * duplicates anyway)
     *
     * @param htf the higher tuple finder
     * @param posNumsHO this list of positions that we are getting higher order
     * factors from
     * @return
     */
    HashSet<Factor> getHigherOrderFactors(HigherTupleFinder<Double> htf, int[] posNumsHO) {
        HashSet<Factor> hoFactorSet = new HashSet<>();
        ArrayList<Integer> interactingPos = htf.getInteractingPos();
        for (int i = 0; i < interactingPos.size(); i++) {
            int posK = interactingPos.get(i);
            //posNumsHO[0] is the smallest pos num, so we just make sure posK is 
            //smaller
            if (posK < posNumsHO[0]) {
                int[] posNumsInFactor = new int[posNumsHO.length + 1];
                posNumsInFactor[0] = posK;
                for (int j = 0; j < posNumsHO.length; j++) {
                    posNumsInFactor[j + 1] = posNumsHO[j];
                }
                for (int rcK : this.pruneMat.unprunedRCsAtPos(posK)) {
                    if (htf.getInteraction(posK, rcK) != 0.0) {
                        Factor hoFactor = new Factor(posNumsInFactor);
                        hoFactorSet.add(hoFactor);
                    }
                    HigherTupleFinder<Double> htf2 = htf.getHigherInteractions(posK, rcK);
                    if (htf2 != null) {
                        HashSet<Factor> hoFactorSet2 = getHigherOrderFactors(htf2, posNumsInFactor);
                        for (Factor hoFactor : hoFactorSet2) {
                            hoFactorSet.add(hoFactor);
                        }
                    }
                }
            }
        }
        return hoFactorSet;
    }

    class Factor {

        int[] posNums;
        int size;

        public Factor(int... posNums) {
            if (posNums.length == 0) {
                throw new RuntimeException("CANNOT CREATE FACTOR WITH 0 Pos Nums");
            } else {
                this.posNums = posNums;
                Arrays.sort(posNums);
                this.size = posNums.length;
            }
        }

        @Override
        public boolean equals(Object obj) {
            //return false if not an instance of Factor
            if (!(obj instanceof Factor)) {
                return false;
            }
            //if it is an instance, typecase to compare data
            Factor f = (Factor) obj;
            //Check to see if posNums are the same
            //First we see if the lengths are the same
            if (!(this.size == f.size)) {
                return false;
            }
            //If they are, then we iterate over one set of posNums and make sure
            //the other shares the same posNums
            for (int pos : this.posNums) {
                //see if factor f contains this posNum
                if (!(ArrayUtils.contains(f.posNums, pos))) {
                    return false;
                }
            }
            return true;
        }

        @Override
        public int hashCode() {
            int hash = 7;
            hash = 79 * hash + Arrays.hashCode(this.posNums);
            return hash;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(Integer.toString(this.posNums[0]));
            for (int i = 1; i < this.size; i++) {
                sb.append(" ");
                sb.append(Integer.toString(this.posNums[i]));
            }
            return sb.toString();
        }

    }

    /**
     * A lambda expression to compare two factors factorI is less than factorJ
     * is the size of factorI is less than factorJ If they are the same size
     * then we compare elements from the beginning to the end
     */
    Comparator<Factor> compare = (Factor o1, Factor o2) -> {
        //if the size of the factor is the same, then compare based on elements
        if (o1.size == o2.size) {
            for (int i = 0; i < o1.size; i++) {
                if (o1.posNums[i] < o2.posNums[i]) {
                    return -1;
                }
                if (o1.posNums[i] > o2.posNums[i]) {
                    return 1;
                }
            }
            //if they contain the same elements, they are equal
            return 0;
        } else {
            if (o1.size < o2.size) {
                return -1;
            }
            return 1;
        }
    };

    /**
     * Creates a 3D message matrix indexed by factor, receiving position,
     * receiving rot Thus we are holding the messages from a factor, to the
     * position/rot receiving the message
     *
     * @param rotsPerPos the number of rotamers per position
     * @param initVal the initial value of the messages (usually 0.0)
     * @return
     */
    double[][][] get3DMsgMat(int[] rotsPerPos, double initVal) {
        double[][][] msgMat = new double[this.numFactors][][];
        for (int factorNum = 0; factorNum < this.numFactors; factorNum++) {
            Factor f = this.factorList.get(factorNum);
            msgMat[factorNum] = new double[f.size][];
            for (int pos = 0; pos < f.size; pos++) {
                int posNum = f.posNums[pos];
                msgMat[factorNum][pos] = new double[rotsPerPos[posNum]];
                for (int rot = 0; rot < rotsPerPos[posNum]; rot++) {
                    msgMat[factorNum][pos][rot] = initVal;
                }
            }
        }
        return msgMat;
    }
}
