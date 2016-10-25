package edu.duke.cs.osprey.energy;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.ContextPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class GpuEnergyFunctionGenerator extends EnergyFunctionGenerator {
	
	public static class NoGpusException extends Exception {

		private static final long serialVersionUID = -8696449488422888399L;
	}
	
	public static GpuEnergyFunctionGenerator make(ForcefieldParams ffParams)
	throws NoGpusException {
		
		// if CUDA is supported, prefer that, since it's faster
		List<edu.duke.cs.osprey.gpu.cuda.Gpu> cudaGpus = edu.duke.cs.osprey.gpu.cuda.Gpus.get().getGpus();
		if (!cudaGpus.isEmpty()) {
			return new GpuEnergyFunctionGenerator(ffParams, new ContextPool());
		}
		
		// otherwise, try OpenCL
		List<edu.duke.cs.osprey.gpu.opencl.Gpu> openclGpus = edu.duke.cs.osprey.gpu.opencl.Gpus.get().getGpus();
		if (!openclGpus.isEmpty()) {
			return new GpuEnergyFunctionGenerator(ffParams, new GpuQueuePool());
		}
		
		throw new NoGpusException();
	}
	
	private GpuQueuePool openclQueues;
	private ContextPool cudaContexts;
	
	public GpuEnergyFunctionGenerator(ForcefieldParams ffParams, GpuQueuePool openclQueues) {
		super(ffParams, Double.POSITIVE_INFINITY, false);
		this.openclQueues = openclQueues;
		this.cudaContexts = null;
	}
	
	public GpuEnergyFunctionGenerator(ForcefieldParams ffParams, ContextPool cudaContexts) {
		super(ffParams, Double.POSITIVE_INFINITY, false);
		this.openclQueues = null;
		this.cudaContexts = cudaContexts;
	}
	
	public GpuQueuePool getOpenclQueuePool() {
		return openclQueues;
	}
	
	public ContextPool getCudaContexts() {
		return cudaContexts;
	}
	
	private GpuForcefieldEnergy makeGpuForcefield(ForcefieldInteractions interactions) {
		if (openclQueues != null) {
			return new GpuForcefieldEnergy(ffParams, interactions, openclQueues);
		} else if (cudaContexts != null) {
			return new GpuForcefieldEnergy(ffParams, interactions, cudaContexts);
		} else {
			throw new Error("bad gpu queue/context config, this is a bug");
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
	public GpuForcefieldEnergy intraAndShellEnergy(Residue res, List<Residue> shellResidues) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		for (Residue shellRes : shellResidues) {
			interactions.addResiduePair(res, shellRes);
		}
		return makeGpuForcefield(interactions);
	}
	
	@Override
	public GpuForcefieldEnergy intraAndDistributedShellEnergy(Residue res, List<Residue> shellResidues, int numPos, double singleWeight) {
		checkSingleWeight(singleWeight);
		return intraAndShellEnergy(res, shellResidues);
	}

	@Override
	public GpuForcefieldEnergy resPairAndDistributedShellEnergy(Residue res1, Residue res2, List<Residue> shellResidues, int numPos, double singleWeight) {
		checkSingleWeight(singleWeight);
		return resPairEnergy(res1, res2);
	}
	
	@Override
	public GpuForcefieldEnergy fullConfEnergy(ConfSpace confSpace, List<Residue> shellResidues) {
		return fullConfEnergy(confSpace, shellResidues, null);
	}
	
	@Override
	public GpuForcefieldEnergy fullConfEnergy(ConfSpace confSpace, List<Residue> shellResidues, Molecule mol) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (int pos1=0; pos1<confSpace.posFlex.size(); pos1++) {
			Residue res1 = confSpace.posFlex.get(pos1).res;
			res1 = matchResidue(res1, mol);
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (int pos2=0; pos2<pos1; pos2++) {
				Residue res2 = confSpace.posFlex.get(pos2).res;
				res2 = matchResidue(res2, mol);
				interactions.addResiduePair(res1, res2);
			}
			
			// shell energies
			for (Residue shellRes : shellResidues) {
				shellRes = matchResidue(shellRes, mol);
				interactions.addResiduePair(res1, shellRes);
			}
		}
		
		return makeGpuForcefield(interactions);
	}
	
	private Residue matchResidue(Residue res, Molecule mol) {
		if (mol == null) {
			return res;
		}
		return mol.residues.get(res.indexInMolecule);
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
	
	public void cleanup() {
		if (openclQueues != null) {
			openclQueues.cleanup();
		}
	}
}
