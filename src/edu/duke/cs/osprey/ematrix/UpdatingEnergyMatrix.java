/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupETrie;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.confspace.TupE;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;


public abstract class UpdatingEnergyMatrix<T extends TupE> extends ProxyEnergyMatrix {
    // Store the seen confs in a trie with wildcards.
    public static boolean debug = false;
    private TupETrie<T> corrections;
    private int numPos;
    
    //debug variable
    public final ConfEnergyCalculator sourceECalc;

    //Abstract constructors for the type of TupE we want to use
    protected abstract T makeT(RCTuple tup, double val);
    protected abstract TupETrie<T> makeTrie(List<SimpleConfSpace.Position> positions);

    public UpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target, ConfEnergyCalculator confECalc) {
        super(confSpace, target);
        this.corrections = makeTrie(confSpace.positions);
        this.numPos = confSpace.getNumPos();
        this.sourceECalc = confECalc;

    }

    public UpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
        this.numPos = confSpace.getNumPos();
        this.sourceECalc = null;
        this.corrections = makeTrie(confSpace.positions);
    }

    public List<T> getAllCorrections(){
        return corrections.getAllEntries();
    }

    public int getTrieSize(){
        return this.corrections.size();
    }

    /*Hack 1: Don't share residues*/
    @Override
    public boolean hasHigherOrderTerms() {
        return corrections.size() > 0;
    }

    @Override
    public boolean hasHigherOrderTermFor(RCTuple query) {
        return corrections.contains(query);
    }

    public String formatCorrections(List<T> corrections) {
        String out = "";
        for(T correction:corrections)
            out+=correction.tup.stringListing()+":"+correction.E+"\n";
            return out;
    }

    public double confE(int conf[], double parentCorrect){
        /**
         * Due to the greedy edge case, we need to ensure that our corrections never decrease
         */
        return getInternalEnergy(new RCTuple(conf), parentCorrect) + constTerm;
    }

    /*
    We don't want to use this like an energyMatrix, we want to just get corrections. This
    is essentially because we can't trust the greedy heuristic, so we can't feed it only
    the conformation and expect to get good results.
     */
    @Deprecated
    @Override
    public double confE(int[] conf){
        throw new UnsupportedOperationException();
    }
    /*
    We don't want to use this like an energyMatrix, we want to just get corrections. This
    is essentially because we can't trust the greedy heuristic, so we can't feed it only
    the conformation and expect to get good results.
     */
    @Deprecated
    public double getInternalEnergy(RCTuple tup, double parentCorrect){
        //internal energy of a tuple of residues when they're in the specified RCs

        // OPTIMIZATION: don't even check higher terms if the energy matrix doesn't have any
        // this does wonders to CPU cache performance!
        boolean useHigherOrderTerms = hasHigherOrderTerms();

        ArrayList<Integer> tuppos = tup.pos;
        ArrayList<Integer> tupRCs = tup.RCs;

        // OPTIMIZATION: split oneBody and pairwise energies into separate loops
        // to improve CPU cache performance

        int numPosInTuple = tup.pos.size();
        double energy = 0;

        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tuppos.get(indexInTuple);
            int RCNum = tupRCs.get(indexInTuple);

            energy += getOneBody(posNum,RCNum);
        }

        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tuppos.get(indexInTuple);
            int RCNum = tupRCs.get(indexInTuple);

            for(int index2=0; index2<indexInTuple; index2++){
                int pos2 = tuppos.get(index2);
                int rc2 = tupRCs.get(index2);

                energy += getPairwise(posNum,RCNum,pos2,rc2);

            }
        }
        if (useHigherOrderTerms) {
            // Since our corrections greediness can get us into trouble, don't accept corrections smaller than the parent
            double newCorrect = internalEHigherOrder(tup);
            energy+=Math.max(newCorrect, parentCorrect);
        }

        return energy;
    }

    @Deprecated
    public double getInternalEnergy(RCTuple tup){
        //internal energy of a tuple of residues when they're in the specified RCs

        // OPTIMIZATION: don't even check higher terms if the energy matrix doesn't have any
        // this does wonders to CPU cache performance!
        boolean useHigherOrderTerms = hasHigherOrderTerms();

        ArrayList<Integer> tuppos = tup.pos;
        ArrayList<Integer> tupRCs = tup.RCs;

        // OPTIMIZATION: split oneBody and pairwise energies into separate loops
        // to improve CPU cache performance

        int numPosInTuple = tup.pos.size();
        double energy = 0;

        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tuppos.get(indexInTuple);
            int RCNum = tupRCs.get(indexInTuple);

            energy += getOneBody(posNum,RCNum);
        }

        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tuppos.get(indexInTuple);
            int RCNum = tupRCs.get(indexInTuple);

            for(int index2=0; index2<indexInTuple; index2++){
                int pos2 = tuppos.get(index2);
                int rc2 = tupRCs.get(index2);

                energy += getPairwise(posNum,RCNum,pos2,rc2);

            }
        }
        if (useHigherOrderTerms) {
            energy+=internalEHigherOrder(tup);
        }

        return energy;
    }


    /*
    public double getPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2){
    //lookup by residue number, amino-acid index (for the residue's rotamer library),
    //and RC index for that residue and AA type (rotamer index for AA type if rigid backbone,
    //otherwise defined in the ConfSearchSpace)
    //RETURN ERROR IF RES1 AND RES2 ARE NOT SINGLE RESIDUES
    return getPairwise(res1, index1, res2, index2);
    }
     */
    //intra+shell similar...

    public double getCorrection(int[] conf){
        /**
         * Returns a Higher-Order Tuple energy correction for the input conf.
         *
         * Although ideally we would want this to be the *largest* possible HOT correction,
         * AFAIK that currently requires solving the maximum weighted set-packing problem.
         *
         */
       return internalEHigherOrder(new RCTuple(conf));
    }

    private double internalEHigherOrder(RCTuple tup){
        /**
         * Returns a Higher-Order Tuple energy correction for the input RCTuple.
         *
         * Although ideally we would want this to be the *largest* possible HOT correction,
         * AFAIK that currently requires solving the maximum weighted set-packing problem.
         *
         */
        double E = 0;
        List<T> confCorrections = corrections.getEntries(tup);
        if(confCorrections.size() > 0) {
            double corr = processCorrections(confCorrections);
            E += corr;
        }

        return E;
    }

    private double processCorrections(List<T> confCorrections) {
        /**
         * Returns a high-weight set packing of the list of corrections.
         *
         * Although ideally we would want this to be the *heaviest* possible set packing,
         * AFAIK that currently requires solving the maximum weighted set-packing problem,
         * which is in NP-hard.
         *
         * So, this uses a greedy strategy for now. This does not always return the optimal solution
         * in practice. Beware.
         *
         */
        Collections.sort(confCorrections, (a,b)->-Double.compare(a.E,b.E));
        double sum = 0;
        // Attempt 1: be greedy and start from the largest correction you
        // can get instead of trying to solve the NP-Complete problem.
        Set<Integer> usedPositions = new HashSet<>();
        List<T> usedCorrections = new ArrayList<>();
        int numApplied = 0;
        for(T correction: confCorrections) {
            if (usedPositions.size() >= numPos) {
                break;
            }
            Collection<Integer> positions = correction.tup.pos;
            boolean noIntersections = true;
            for(int position : positions) {
                if(usedPositions.contains(position)) {
                    noIntersections = false;
                    break;
                }
            }
            if(noIntersections) {
                usedPositions.addAll(correction.tup.pos);
                usedCorrections.add(correction);
                //System.out.println("Applying correction "+correction.tup.stringListing()+":"+correction.E);
                numApplied++;
                sum += correction.E;
            }
        }
        if(debug)
        {
            Set<Integer> positionCheck = new HashSet<>();
            for(T correction: usedCorrections) {
                for(int pos: correction.tup.pos) {
                    if (positionCheck.contains(pos))
                        System.err.println("REUSING POSITION "+pos);
                }
            }
        }
        return sum;
    }

    public void insertAll(List<T> corrList){
        for (T correction : corrList){
            corrections.insert(correction);
        }
    }

    @Override
    public void setHigherOrder(RCTuple tup, Double val) {
        /*
        if(tup.size() < 3)
        {
            System.err.println("Should not be trying to submit correction of lower-order term.");
            return;
        }
        */
        RCTuple orderedTup = tup.sorted();
        corrections.insert(makeT(orderedTup, val));
    }

    public void writeCorrectionsToFile(String filename){
        corrections.writeEntriesToFile(filename);
    }
    public void writeCorrectionsToFile(File f) {
        corrections.writeEntriesToFile(f);
    }

    public void readCorrectionsFromFile(String filename){
        corrections.readEntriesFromFile(filename);
    }
    public void readCorrectionsFromFile(File f){
        corrections.readEntriesFromFile(f);
    }

}
