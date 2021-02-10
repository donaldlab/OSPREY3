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

package edu.duke.cs.osprey.minimization;

import java.io.IOException;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.CCDKernelCuda;

public class CudaCCDMinimizer implements Minimizer.NeedsCleanup, Minimizer.Reusable {
	
	private GpuStreamPool streams;
	private ObjectiveFunction f;
	private GpuStream stream;
	private CCDKernelCuda kernel;
	private ObjectiveFunction.DofBounds dofBounds;

	public CudaCCDMinimizer(GpuStreamPool streams) {
		this.streams = streams;
	}
	
	public CudaCCDMinimizer(GpuStreamPool streams, ObjectiveFunction f) {
		this(streams);
		init(f);
	}
	
	@Override
	public void init(ObjectiveFunction f) {
		
		this.f = f;
		
		if (kernel == null) {
			// make the kernel
			try {
				stream = streams.checkout();
				kernel = new CCDKernelCuda(stream);
			} catch (IOException ex) {
				streams.release(stream);
				stream = null;
				throw new Error("can't make CCD kernel", ex);
			}
		}
		
		// get the molecule objective function
		if (f instanceof MoleculeModifierAndScorer) {
			kernel.init(new MoleculeObjectiveFunction((MoleculeModifierAndScorer)f));
		} else if (f instanceof MoleculeObjectiveFunction) {
			kernel.init((MoleculeObjectiveFunction)f);
		} else {
			throw new Error("objective function should be a " + MoleculeModifierAndScorer.class.getSimpleName() + ", not a " + f.getClass().getSimpleName() + ". this is a bug");
		}

		dofBounds = new ObjectiveFunction.DofBounds(f.getConstraints());
	}

	@Override
	public Minimizer.Result minimizeFromCenter() {

		DoubleMatrix1D x = DoubleFactory1D.dense.make(dofBounds.size());
		dofBounds.getCenter(x);

		return minimizeFrom(x);
	}

	@Override
	public Minimizer.Result minimizeFrom(DoubleMatrix1D x) {
		
		// do the minimization
		f.setDOFs(x);
		kernel.uploadCoordsAsync();
		Minimizer.Result result = kernel.runSync(x, dofBounds);
		
		// update the CPU-side molecule
		f.setDOFs(result.dofValues);
		
		return result;
	}

	@Override
	public void clean() {
		if (kernel != null) {
			kernel.cleanup();
			kernel = null;
		}
		streams.release(stream);
	}
}
