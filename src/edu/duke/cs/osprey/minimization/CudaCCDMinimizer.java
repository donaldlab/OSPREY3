package edu.duke.cs.osprey.minimization;

import java.io.IOException;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.CCDKernelCuda;

public class CudaCCDMinimizer implements Minimizer.NeedsCleanup, Minimizer.Reusable {
	
	private GpuStreamPool streams;
	private MoleculeModifierAndScorer mof;
	private GpuStream stream;
	private CCDKernelCuda kernel;
	private ObjectiveFunction.DofBounds dofBounds;
	private DoubleMatrix1D x;

	public CudaCCDMinimizer(GpuStreamPool streams) {
		this.streams = streams;
	}
	
	public CudaCCDMinimizer(GpuStreamPool streams, ObjectiveFunction f) {
		this(streams);
		init(f);
	}
	
	@Override
	public void init(ObjectiveFunction f) {
		
		// get the molecule objective function
		if (f instanceof MoleculeModifierAndScorer) {
			mof = (MoleculeModifierAndScorer)f;
		} else {
			throw new Error("objective function should be a " + MoleculeModifierAndScorer.class.getSimpleName() + ", not a " + f.getClass().getSimpleName() + ". this is a bug");
		}
		
		if (kernel == null) {
			// make the kernel
			try {
				stream = streams.checkout();
				kernel = new CCDKernelCuda(this.stream);
			} catch (IOException ex) {
				throw new Error("can't make CCD kernel", ex);
			}
		}
		kernel.init(mof);
		
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
