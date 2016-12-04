package edu.duke.cs.osprey.energy;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class GpuEnergyFunctionGenerator extends EnergyFunctionGenerator {
	
	private ForcefieldInteractionsGenerator intergen;
	private GpuQueuePool openclQueues;
	private GpuStreamPool cudaStreams;
	
	private GpuEnergyFunctionGenerator(ForcefieldParams ffParams) {
		super(ffParams, Double.POSITIVE_INFINITY, false);
		this.intergen = new ForcefieldInteractionsGenerator();
		this.openclQueues = null;
		this.cudaStreams = null;
	}
	
	public GpuEnergyFunctionGenerator(ForcefieldParams ffParams, GpuQueuePool openclQueues) {
		this(ffParams);
		this.openclQueues = openclQueues;
	}
	
	public GpuEnergyFunctionGenerator(ForcefieldParams ffParams, GpuStreamPool cudaContexts) {
		this(ffParams);
		this.cudaStreams = cudaContexts;
	}
	
	public GpuQueuePool getOpenclQueuePool() {
		return openclQueues;
	}
	
	public GpuStreamPool getCudaContexts() {
		return cudaStreams;
	}
	
	private GpuForcefieldEnergy makeGpuForcefield(ForcefieldInteractions interactions) {
		if (openclQueues != null) {
			return new GpuForcefieldEnergy(ffParams, interactions, openclQueues);
		} else if (cudaStreams != null) {
			return new GpuForcefieldEnergy(ffParams, interactions, cudaStreams);
		} else {
			throw new Error("bad gpu queue/context config, this is a bug");
		}
	}
	
	@Override
	public GpuForcefieldEnergy singleResEnergy(Residue res) {
		return makeGpuForcefield(intergen.makeSingleRes(res));
	}

	@Override
	public GpuForcefieldEnergy resPairEnergy(Residue res1, Residue res2) {
		return makeGpuForcefield(intergen.makeResPair(res1, res2));
	}
	
	@Override
	public GpuForcefieldEnergy intraAndShellEnergy(Residue res, List<Residue> shellResidues) {
		return makeGpuForcefield(intergen.makeIntraAndShell(res, shellResidues));
	}
	
	@Override
	public GpuForcefieldEnergy intraAndDistributedShellEnergy(Residue res, List<Residue> shellResidues, int numPos, double singleWeight) {
		return makeGpuForcefield(intergen.makeIntraAndDistributedShell(res, shellResidues, numPos, singleWeight));
	}

	@Override
	public GpuForcefieldEnergy resPairAndDistributedShellEnergy(Residue res1, Residue res2, List<Residue> shellResidues, int numPos, double singleWeight) {
		return makeGpuForcefield(intergen.makeResPairAndDistributedShell(res1, res2, shellResidues, numPos, singleWeight));
	}
	
	@Override
	public GpuForcefieldEnergy fullConfEnergy(ConfSpace confSpace, List<Residue> shellResidues) {
		return makeGpuForcefield(intergen.makeFullConf(confSpace, shellResidues));
	}
	
	@Override
	public GpuForcefieldEnergy fullConfEnergy(ConfSpace confSpace, List<Residue> shellResidues, Molecule mol) {
		return makeGpuForcefield(intergen.makeFullConf(confSpace, shellResidues, mol));
	}
	
	@Override
	public GpuForcefieldEnergy fullMolecEnergy(Molecule mol) {
		return makeGpuForcefield(intergen.makeFullMol(mol));
	}
	
	public void cleanup() {
		if (openclQueues != null) {
			openclQueues.cleanup();
		}
		if (cudaStreams != null) {
			cudaStreams.cleanup();
		}
	}
}
