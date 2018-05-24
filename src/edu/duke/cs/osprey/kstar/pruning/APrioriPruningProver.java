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

package edu.duke.cs.osprey.kstar.pruning;

import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.gmec.GMECFinder;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import java.util.HashMap;

/**
 *
 * This class selects an I0 and Ew for pruning so we can know ahead of time
 * that the error in partition functions due to skipping pruned conformations is < epsilon. 
 * Key use case: LUTE effectively zeroes out the Boltzmann weight for all pruned conformations,
 * and will likely have trouble if we try to unprune everything for a second round, 
 * so instead we can use this class to get a bounded error in K* in a single pass.  
 * 
 * 
 * @author mhall44
 */
public class APrioriPruningProver {
    
    KSConfigFileParser cfp;
    HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs;
    KSAbstract ks;
    boolean contSCFlex;
    
    /*public APrioriPruningProver(KSConfigFileParser cfp, 
            HashMap<Integer,KSAllowedSeqs> strand2AllowedSeqs){
        this.cfp = cfp;
        this.strand2AllowedSeqs = strand2AllowedSeqs;
    }*/
    public APrioriPruningProver(KSAbstract ks, KSConfigFileParser cfp, 
            HashMap<Integer,KSAllowedSeqs> strand2AllowedSeqs){
        this.ks = ks;
        this.cfp = cfp;
        this.strand2AllowedSeqs = strand2AllowedSeqs;
        contSCFlex = cfp.params.getBool("doMinimize", true);
    }
       
    
    public double calcI0(){
        //I0 is an upper bound, valid for all possibly desirable sequences,
        //on the gap between the lowest lower-bound and the GMEC energy
        //we can compute the lowest lower-bound over the whole sequence space here
        //and the user can provide a stability constraint, which upper bounds
        //how much higher any desirable sequence's GMEC can be compared to the "best" sequence GMEC
        
        double I0 = 0;
        if(!contSCFlex)//bounds will be tight in discrete case
            return 0;
        
        if(!cfp.params.getBool("iMinDEE"))
            throw new RuntimeException("ERROR: Trying to calc I0 but not using iMinDEE");
        
        System.out.println("Determining I0 to get provable K*");
        //OK to be provable I0 has to be the biggest among all strand
        
        for( int strand : new int[] {0, 1, 2} ){
            SearchProblem strandSP = ks.createPanSeqSP( contSCFlex, strand );//cfp.getSearchProblem(strand, strand2AllowedSeqs.get(strand));
            GMECFinder strandGF = new GMECFinder();
            strandGF.init(cfp, strandSP);
            double gmec = strandGF.calcGMECEnergy();
            double lowestBound = strandGF.getLowestBound();
            I0 = Math.max(I0, gmec-lowestBound);
        }
        
        double stabilityGap = cfp.params.getDouble("StabilityGap",10);
        //Add the stability constraint (how much less stable we can be than the most stable sequence)
        I0 += stabilityGap;
        
        System.out.println("I0 to get provable K*: "+I0);
        return I0;
    }
    
    
    public double calcEw(){
        //Let's demand that p* < (1-epsilon) * (q*)
        //p* <= (total_num_confs) exp(-(E_GMEC+Ew)/RT), while q* >=  exp(-E_GMEC/RT)
        //(these bounds are probably fairly loose but it's hard to do better a priori)
        //so let's set exp(-Ew/RT) * upper_bound(total_num_confs) = epsilon
        double epsilon = cfp.params.getDouble("epsilon");
        if(epsilon<=0 || epsilon >=1)
            throw new RuntimeException("ERROR: Can't calc Ew based on this epsilon: "+epsilon);
        
        double Ew = -PFAbstract.RT * Math.log(epsilon);
        
        //to get an upper bound on the total number of conformations in a sequence,
        //take product_i max_a num_RCs(i,a), where i is residue positions & a is amino-acid types
        //the complex will have the most possible conformations
        SearchProblem complexSP = ks.createPanSeqSP(contSCFlex, 2);//cfp.getSearchProblem(KS2, strand2AllowedSeqs.get(KS2));
        for(int pos=0; pos<complexSP.confSpace.numPos; pos++){
            PositionConfSpace pcs = complexSP.confSpace.posFlex.get(pos);
            Ew += PFAbstract.RT * Math.log(pcs.maxNumRCsPerResType());
        }
        
        System.out.println("Ew to get provable K*: "+Ew);
        return Ew;
    }
    
    
    
}
