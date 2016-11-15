package edu.duke.cs.osprey.minimization;

import java.io.IOException;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
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

	public CudaCCDMinimizer(GpuStreamPool streams, MoleculeModifierAndScorer mof) {
		
		this.streams = streams;
		this.mof = mof;
		this.stream = streams.checkout();
		
		// make the kernel
		try {
			kernel = new CCDKernelCuda(this.stream, mof);
		} catch (IOException ex) {
			throw new Error("can't make CCD kernel", ex);
		}
		
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
