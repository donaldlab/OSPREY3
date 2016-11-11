package edu.duke.cs.osprey.minimization;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.CCDKernelCuda;

public class CudaCCDMinimizer implements Minimizer.NeedsCleanup {
	
	private GpuStreamPool streams;
	private MoleculeModifierAndScorer mof;
	private GpuStream stream;
	private CCDKernelCuda kernel;
	private ObjectiveFunction.DofBounds dofBounds;
	private DoubleMatrix1D x;

	public CudaCCDMinimizer(GpuStreamPool streams, MoleculeModifierAndScorer mof)
	throws IOException {
		
		this.streams = streams;
		this.mof = mof;
		this.stream = streams.checkout();
		
		// collects the dofs
		List<FreeDihedral> dofs = new ArrayList<>();
		for (DegreeOfFreedom dof : mof.getDOFs()) {
			
			// make sure it's a FreeDihedral. that's all we support on the gpu at the moment
			if (dof instanceof FreeDihedral) {
				dofs.add((FreeDihedral)dof);
			} else {
				throw new Error("degree-of-freedom type " + dof.getClass().getSimpleName() + " not yet supported by GPU CCD minimzer."
					+ " Use CPU minimizer with GPU energy function instead");
			}
		}
		
		// get the energy function
		EnergyFunction efunc = mof.getEfunc();
		BigForcefieldEnergy ffenergy;
		if (efunc instanceof BigForcefieldEnergy) {
			ffenergy = (BigForcefieldEnergy)efunc;
		} else {
			throw new Error("GPU CCD minimizer needs a " + BigForcefieldEnergy.class.getSimpleName() + ", not a " + efunc.getClass().getSimpleName() + ". this is a bug.");
		}
		
		kernel = new CCDKernelCuda(this.stream, ffenergy, dofs);
		
		// init x to the center of the bounds
		dofBounds = new ObjectiveFunction.DofBounds(mof.getConstraints());
		x = DoubleFactory1D.dense.make(dofBounds.size());
		dofBounds.getCenter(x);
	}
	
	@Override
	public Minimizer.Result minimize() {
		
		// do the minimization
		mof.setDOFs(x);
		kernel.uploadCoordsAsync();
		Minimizer.Result result = kernel.runSync(x, dofBounds);
		
		// update the CPU-side molecule
		mof.setDOFs(result.dofValues);
		
		return result;
	}

	@Override
	public void cleanup() {
		kernel.cleanup();
		streams.release(stream);
	}
}
