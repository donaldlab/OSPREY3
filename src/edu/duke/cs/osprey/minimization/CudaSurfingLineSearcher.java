package edu.duke.cs.osprey.minimization;

import java.io.IOException;

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelCuda;
import edu.duke.cs.osprey.gpu.cuda.kernels.LineSearchKernel;

public class CudaSurfingLineSearcher implements LineSearcher.NeedsCleanup {
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs
	
	private ObjectiveFunction.OneDof f;
	private GpuForcefieldEnergy efunc;
	private ForcefieldKernelCuda ffKernel;
	private LineSearchKernel lsKernel;
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	@Override
	public void init(ObjectiveFunction.OneDof f) {
		
		this.f = f;
		
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
		
		// HACKHACK: break the interfaces to get the pieces we need
		// the runtime will complain if we can't do these casts, so we're still safe, but the compiler has no idea what's going on
		MoleculeModifierAndScorer mof = (MoleculeModifierAndScorer)f.getParent();
		FreeDihedral dof = (FreeDihedral)mof.getDOFs().get(f.getDimension());
		efunc = ((GpuForcefieldEnergy)mof.getEfunc(f.getDimension()));
		ffKernel = (ForcefieldKernelCuda)efunc.getKernel();
		
		// make the line search kernel
		try {
			lsKernel = new LineSearchKernel(ffKernel, dof);
		} catch (IOException ex) {
			throw new Error(ex);
		}
	}
	
	@Override
	public double search(double xd) {
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		double step;
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// TEMP
		// update forcefield subset
		boolean changed = ffKernel.setSubset(efunc.getSubset());
		if (changed) {
			System.out.println("updated subset!!");
		}
		
		// run the line search kernel
		// why on earth are we using degrees... oi
		lsKernel.runAsync(Math.toRadians(f.getXMin()), Math.toRadians(f.getXMax()), Math.toRadians(xd), Math.toRadians(step));
		
		// download the final outputs from the gpu
		lsKernel.downloadResultSync();
		double bestDihedral = Math.toDegrees(lsKernel.getXdstar());
		
		// update the step
		this.lastStep = bestDihedral - xd;
		if (iteration == 0) {
			firstStep = lastStep;
		}
		
		// pose the cpu-side protein too
		//f.setX(bestDihedral);
		
		iteration++;
		
		return bestDihedral;
	}
	
	@Override
	public void cleanup() {
		if (lsKernel != null) {
			lsKernel.cleanup();
		}
	}
}
