/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ConfSpaceSuper;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class EnergyFunctionGenerator {
    //This is an object that generates an energy function (maps conformation-->energy) for a molecule or a portion thereof
    //it specifies settings for how this energy should be estimated
    
    //Start with AMBER/EEF1, add Poisson-Boltzmann (instead of EEF1), and ultimately QM and explicit water hopefully
    
    public ForcefieldParams ffParams;
    
    public double distCutoff;//distance cutoff for interactions (angstroms)
    
    public boolean usePoissonBoltzmann;//Use Poisson-Boltzmann energies for solvation
    //Note: these will not be included in the single-res and pair energy functions,
    //so we need to do only full-conf energies to get the Poisson-Boltzmann term
    
    
    public EnergyFunctionGenerator(ForcefieldParams fParams, double distC, boolean usePB){
        ffParams = fParams;
        distCutoff = distC;
        usePoissonBoltzmann = usePB;
    }
    
    
    public EnergyFunction resPairEnergy(Residue res1, Residue res2){
        //make an energy function estimating the interaction energy between two residues
        
        //for now definitely using ffParams
        return new ResPairEnergy(res1,res2,ffParams);
    }
    
    public EnergyFunction singleResEnergy(Residue res){
        //make an energy function estimating the energy of a single residue
        return new SingleResEnergy(res,ffParams);
    }
    
    
    public EnergyFunction intraAndShellEnergy(Residue res, ArrayList<Residue> shellResidues){
        //internal energy of res, plus its interactions with the residues in shellRes
        
        EnergyFunction intraE = singleResEnergy(res);
        
        MultiTermEnergyFunction intraAndShellE = new MultiTermEnergyFunction();
        intraAndShellE.addTerm(intraE);
        
        for(Residue shellRes : shellResidues){
            EnergyFunction pairE = resPairEnergy(res,shellRes);
            intraAndShellE.addTerm(pairE);
        }
        
        return intraAndShellE;
    }
    
    //want partial versions of the above too?
    
    
    
    public EnergyFunction fullConfEnergy(ConfSpace cSpace, ArrayList<Residue> shellResidues){
        //make an energy function estimating the energy of all flexible residues in cSpace,
        //plus any interactions with other (shell) residues that are within distCutoff of anything in cSpace
        //this can be used as the full energy for any conformation in cSpace
        
        ArrayList<Residue> flexibleResidues = new ArrayList<>();//array list for these so they stay in order
        for(PositionConfSpace pcs : cSpace.posFlex)
            flexibleResidues.add(pcs.res);
        
        //now we want the full energy (1-body + all interactions) of the flexibleResidues,
        //plus the interactions of the flexibleResidues with the shellResidues
        //(We don't count interactions among the shellResidues because they never change)
        MultiTermEnergyFunction fullEFunc = new MultiTermEnergyFunction();
        
        //start with interactions among flexible residues
        for(int flexResNum=0; flexResNum<flexibleResidues.size(); flexResNum++){
            
            Residue flexRes = flexibleResidues.get(flexResNum);
            EnergyFunction oneBodyE = singleResEnergy(flexRes);
            fullEFunc.addTerm(oneBodyE);
            
            for(int flexResNum2=0; flexResNum2<flexResNum; flexResNum2++){
                Residue flexRes2 = flexibleResidues.get(flexResNum2);
                EnergyFunction pairE = resPairEnergy(flexRes,flexRes2);
                fullEFunc.addTerm(pairE);
            }
        }
        
        //now flexible-to-shell interactions
        for(Residue flexRes : flexibleResidues){
            for(Residue shellRes : shellResidues){
                EnergyFunction pairE = resPairEnergy(flexRes,shellRes);
                fullEFunc.addTerm(pairE);
            }
        }
        
        //now add Poisson-Boltzmann energy, if applicable
        if(usePoissonBoltzmann){
            PoissonBoltzmannEnergy pbe = new PoissonBoltzmannEnergy(cSpace.m);
            fullEFunc.addTermWithCoeff(pbe, ffParams.getSolvScale());
        }
        
        return fullEFunc;
    }
    
  //HMN: Created new fullConfEnergy to hand ConfSpaceSuper input
    public EnergyFunction fullConfEnergy(ConfSpaceSuper cSpace, ArrayList<Residue> shellResidues){
        ArrayList<Residue> flexibleResidues = new ArrayList<>();
        for (PositionConfSpaceSuper pcs : cSpace.posFlex){
            for (Residue res : pcs.resList){
                flexibleResidues.add(res);
            }
        }
        
        MultiTermEnergyFunction fullEfun = new MultiTermEnergyFunction();
        
        //interaction among flexible residues
        for (int flexResNum=0; flexResNum<flexibleResidues.size(); flexResNum++){
            
            Residue flexres = flexibleResidues.get(flexResNum);
            EnergyFunction oneBodyE = singleResEnergy(flexres);
            fullEfun.addTerm(oneBodyE);
            
            for (int flexResNum2=0; flexResNum2<flexResNum; flexResNum2++){
                Residue flexRes2 = flexibleResidues.get(flexResNum2);
                EnergyFunction pairE = resPairEnergy(flexres, flexRes2);
                fullEfun.addTerm(pairE);
            }
        }
        
        //flexible-to-shell interaction
        for (Residue flexRes : flexibleResidues){
            for (Residue shellRes : shellResidues){
                EnergyFunction pairE = resPairEnergy(flexRes, shellRes);
                fullEfun.addTerm(pairE);
            }
        }
        
        return fullEfun;
    }   
    
    public EnergyFunction fullMolecEnergy(Molecule molec){
        //full energy of a molecule, with all residues interacting
        
        MultiTermEnergyFunction fullEFunc = new MultiTermEnergyFunction();
        
        //intra terms
        for(Residue res : molec.residues)
            fullEFunc.addTerm( singleResEnergy(res) );
        
        //pairwise terms
        for(int res1=0; res1<molec.residues.size(); res1++){
            for(int res2=0; res2<res1; res2++){
                fullEFunc.addTerm(resPairEnergy(molec.residues.get(res1),molec.residues.get(res2)));
            }
        }
        
        return fullEFunc;
    }
}
