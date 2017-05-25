/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;

/**
 *
 * @author mhall44
 */
public class EnergyFunctionGenerator {
    //This is an object that generates an energy function (maps conformation-->energy) for a molecule or a portion thereof
    //it specifies settings for how this energy should be estimated
    
    public ForcefieldParams ffParams;
    
    public EnergyFunctionGenerator(ForcefieldParams fParams) {
        ffParams = fParams;
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
    
    
    public EnergyFunction intraAndShellEnergy(Residue res, List<Residue> shellResidues){
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
    
    public EnergyFunction intraAndDistributedShellEnergy(Residue res, List<Residue> shellResidues, int numPos, double singleWeight) {
    	
    	if (singleWeight == 0) {
    		
    		return new SingleResEnergy(res, ffParams);
    		
    	} else {
    	
			MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
			efunc.addTerm(new SingleResEnergy(res, ffParams));
			for (Residue shellRes : shellResidues) {
				efunc.addTerm(scaleEfunc(new ResPairEnergy(res, shellRes, ffParams), singleWeight));
			}
			return efunc;
    	}
    }
    
    public EnergyFunction resPairAndDistributedShellEnergy(Residue res1, Residue res2, List<Residue> shellResidues, int numPos, double singleWeight) {
    	
    	if (singleWeight == 1) {
    		
    		return new ResPairEnergy(res1, res2, ffParams);
    		
    	} else {
    	
			double pairWeight = (1.0 - singleWeight)/(numPos - 1);
			
			MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
			efunc.addTerm(new ResPairEnergy(res1, res2, ffParams));
			for (Residue shellRes : shellResidues) {
				efunc.addTerm(scaleEfunc(new ResPairEnergy(res1, shellRes, ffParams), pairWeight));
				efunc.addTerm(scaleEfunc(new ResPairEnergy(res2, shellRes, ffParams), pairWeight));
			}
			return efunc;
    	}
    }
        
    private EnergyFunction scaleEfunc(EnergyFunction efunc, double scale) {
    	
    	// don't add the scaling function if the scale is just 1
    	if (scale == 1) {
    		return efunc;
    	}
    	
    	return new ScaledEnergyFunction(efunc, scale);
    }
    
    //want partial versions of the above too?
    
    
    
    public EnergyFunction fullConfEnergy(ConfSpace cSpace, List<Residue> shellResidues) {
    	return fullConfEnergy(cSpace, shellResidues, null);
    }
    
    public EnergyFunction fullConfEnergy(ConfSpace cSpace, List<Residue> shellResidues, Molecule mol) {
        //make an energy function estimating the energy of all flexible residues in cSpace,
        //plus any interactions with other (shell) residues that are within distCutoff of anything in cSpace
        //this can be used as the full energy for any conformation in cSpace
        
        List<Residue> flexibleResidues = new ArrayList<>();//array list for these so they stay in order
        for (PositionConfSpace pcs : cSpace.posFlex) {
            flexibleResidues.add(pcs.res);
        }
        
        // match residues to molecule if needed
        flexibleResidues = matchToMolecule(flexibleResidues, mol);
        shellResidues = matchToMolecule(shellResidues, mol);
        
        //now we want the full energy (1-body + all interactions) of the flexibleResidues,
        //plus the interactions of the flexibleResidues with the shellResidues
        //(We don't count interactions among the shellResidues because they never change)
        MultiTermEnergyFunction fullEFunc = new MultiTermEnergyFunction();
        
        //start with interactions among flexible residues
        for (int flexResNum=0; flexResNum<flexibleResidues.size(); flexResNum++) {
            
            Residue flexRes = flexibleResidues.get(flexResNum);
            fullEFunc.addTerm(singleResEnergy(flexRes));
            
            for (int flexResNum2=0; flexResNum2<flexResNum; flexResNum2++) {
                Residue flexRes2 = flexibleResidues.get(flexResNum2);
                fullEFunc.addTerm(resPairEnergy(flexRes, flexRes2));
            }
        }
        
        //now flexible-to-shell interactions
        for (Residue flexRes : flexibleResidues) {
            for (Residue shellRes : shellResidues) {
                fullEFunc.addTerm(resPairEnergy(flexRes, shellRes));
            }
        }
        
        //now add Poisson-Boltzmann energy, if applicable
        if (ffParams.solvationForcefield == SolvationForcefield.PoissonBoltzmann) {
            fullEFunc.addTermWithCoeff(new PoissonBoltzmannEnergy(cSpace.m), ffParams.solvScale);
        }
        
        return fullEFunc;
    }
        
    private List<Residue> matchToMolecule(List<Residue> residues, Molecule mol) {
        
        // no molecule, no matching needed
        if (mol == null) {
            return residues;
        }
        
        // otherwise, match all the residues
        List<Residue> matched = new ArrayList<>(residues.size());
        for (Residue res : residues) {
            matched.add(mol.residues.get(res.indexInMolecule));
        }
        return matched;
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
    
    public EnergyFunction interactionEnergy(ForcefieldInteractions interactions) {
        MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
        for (ForcefieldInteractions.AtomGroup[] groups : interactions) {
            // NOTE: this cast only works when the interactions are between residues
            // for now, that's always true, but maybe it won't be in the future?
            ForcefieldInteractions.ResidueAtomGroup group0 = (ForcefieldInteractions.ResidueAtomGroup)groups[0];
            ForcefieldInteractions.ResidueAtomGroup group1 = (ForcefieldInteractions.ResidueAtomGroup)groups[1];
            if (group0 == group1) {
                efunc.addTerm(singleResEnergy(group0.getResidue()));
            } else {
                efunc.addTerm(resPairEnergy(group0.getResidue(), group1.getResidue()));
            }
        }
        return efunc;
    }
    
    public EnergyFunction residueInteractionEnergy(Residues residues, ResidueInteractions interactions) {
    	
        MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
        
        for (ResidueInteractions.Pair pair : interactions) {
            Residue res1 = residues.getOrThrow(pair.resNum1);
            Residue res2 = residues.getOrThrow(pair.resNum2);
            
            EnergyFunction term;
            if (res1 == res2) {
                term = singleResEnergy(res1);
            } else {
                term = resPairEnergy(res1, res2);
            }
            
            // apply weight
            if (pair.weight != 1.0) {
            	term = new ScaledEnergyFunction(term, pair.weight);
            }
            
            // apply offset
            if (pair.offset != 0.0) {
            	final EnergyFunction fterm = term;
            	final double foffset = pair.offset;
            	term = new EnergyFunction() {

					private static final long serialVersionUID = 8343782653796072786L;

					@Override
					public double getEnergy() {
						return fterm.getEnergy() + foffset;
					}
            	};
            }
            
            efunc.addTerm(term);
        }
        return efunc;
    }
}
