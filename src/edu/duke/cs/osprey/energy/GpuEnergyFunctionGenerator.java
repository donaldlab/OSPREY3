package edu.duke.cs.osprey.energy;

import java.io.IOException;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.GpuInitException;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class GpuEnergyFunctionGenerator extends EnergyFunctionGenerator {
	
	public GpuEnergyFunctionGenerator(ForcefieldParams ffParams) {
		super(ffParams, Double.POSITIVE_INFINITY, false);
	}
	
	private GpuForcefieldEnergy makeGpuForcefield(ForcefieldInteractions interactions) {
		try {
			GpuForcefieldEnergy ff = new GpuForcefieldEnergy(ffParams, interactions);
			ff.initGpu();
			return ff;
		} catch (IOException ex) {
			throw new GpuInitException("can't init gpu forcefield", ex);
		}
	}
	
	@Override
	public GpuForcefieldEnergy singleResEnergy(Residue res) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		return makeGpuForcefield(interactions);
	}

	@Override
	public GpuForcefieldEnergy resPairEnergy(Residue res1, Residue res2) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResiduePair(res1, res2);
		return makeGpuForcefield(interactions);
	}
	
	@Override
	public GpuForcefieldEnergy intraAndShellEnergy(Residue res, ArrayList<Residue> shellResidues) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		for (Residue shellRes : shellResidues) {
			interactions.addResiduePair(res, shellRes);
		}
		return makeGpuForcefield(interactions);
	}
	
	@Override
	public GpuForcefieldEnergy intraAndDistributedShellEnergy(Residue res, ArrayList<Residue> shellResidues, int numPos, double singleWeight) {
		checkSingleWeight(singleWeight);
		return intraAndShellEnergy(res, shellResidues);
	}

	@Override
	public GpuForcefieldEnergy resPairAndDistributedShellEnergy(Residue res1, Residue res2, ArrayList<Residue> shellResidues, int numPos, double singleWeight) {
		checkSingleWeight(singleWeight);
		return resPairEnergy(res1, res2);
	}
	
	@Override
	public GpuForcefieldEnergy fullConfEnergy(ConfSpace confSpace, ArrayList<Residue> shellResidues) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (int pos1=0; pos1<confSpace.posFlex.size(); pos1++) {
			Residue res1 = confSpace.posFlex.get(pos1).res;
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (int pos2=0; pos2<pos1; pos2++) {
				Residue res2 = confSpace.posFlex.get(pos2).res;
				interactions.addResiduePair(res1, res2);
			}
			
			// shell energies
			for (Residue shellRes : shellResidues) {
				interactions.addResiduePair(res1, shellRes);
			}
		}
		
		return makeGpuForcefield(interactions);
	}
	
	@Override
	public GpuForcefieldEnergy fullMolecEnergy(Molecule mol) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (Residue res1 : mol.residues) {
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (Residue res2 : mol.residues) {
				interactions.addResiduePair(res1, res2);
			}
		}
		
		return makeGpuForcefield(interactions);
	}
	
	private void checkSingleWeight(double singleWeight) {
		
		// only single weight of 1 supported
		if (singleWeight != 1) {
			throw new UnsupportedOperationException("only all-on-pairs shell energy distribution supported in GPU forcefields so far");
		}
	}
}
