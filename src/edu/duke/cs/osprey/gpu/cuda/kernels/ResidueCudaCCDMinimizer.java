package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache.ResPair;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.MathTools;
import jcuda.Pointer;

public class ResidueCudaCCDMinimizer extends Kernel implements Minimizer.NeedsCleanup {
	
	private class Dihedral {
		
		public final int d;
		public final Residue res;
		public final int[] dihedralIndices;
		public final int[] rotatedIndices;
		
		public int resIndex;
		public int[] resPairIndices;
		public double xdmin;
		public double xdmax;
		
		public Dihedral(int d, FreeDihedral diehedral) {
			this.d = d;
			this.res = diehedral.getResidue();
			this.dihedralIndices = res.template.getDihedralDefiningAtoms(diehedral.getDihedralNumber());
			this.rotatedIndices = toArray(res.template.getDihedralRotatedAtoms(diehedral.getDihedralNumber()));
		}
		
		private int[] toArray(List<Integer> list) {
			int[] array = new int[list.size()];
			for (int i=0; i<list.size(); i++) {
				array[i] = list.get(i);
			}
			return array;
		}
	}
	
	public final GpuStreamPool streams;
	
	private Kernel.Function func;
	private static AtomicInteger blockThreads = new AtomicInteger(-1);
	
	private MoleculeObjectiveFunction mof;
	private ResidueForcefieldEnergy efunc;
	
	private List<Dihedral> dihedrals;
	
	/* buffer layout:
	 * NOTE: try to use 8-byte alignments to be happy on 64-bit machines
	 * 
	 * int flags (useHEs, useHvdW, distDepDielect, useEEF1)
	 * int numDihedrals
	 * int numResPairs
	 * int maxNumAtoms
	 * double coulombFactor
	 * double scaledCoulombFactor
	 *  
	 * for each dihedral:
	 *    long offset
	 * 
	 * for each residue pair:
	 *    long offset
	 *    
	 * for each dihedral:
	 *    int resIndex
	 *    int numAtoms;
	 *    long atomsOffset;
	 *    short[4] angleAtomOffsets // in index space of residue
	 *    int numRotatedAtoms
	 *    int numResPairs
	 *    double xdmin
	 *    double xdmax
	 *    
	 *    for each rotated atom: // word-aligned to 8 bytes
	 *       short atomOffset; // in index space of residue
	 *    
	 *    for each res pair: // word-aligned to 8 bytes
	 *       int index
	 *    
	 * for each residue pair:
	 *    long numAtomPairs
	 *    long atomsOffset1
	 *    long atomsOffset2
	 *    double weight
	 *    double offset
	 *    
	 *    // NOTE: use struct-of-arrays here so GPU memory accesses are coalesced
	 *    for each atom pair:
	 *       long flags  (bit isHeavyPair, bit is14Bonded, 6 bits space, 3 byte space, short atomOffset1, short atomOffset2)
	 *    for each atom pair:
	 *       double charge
	 *    for each atom pair:
	 *       double Aij
	 *    for each atom pair:
	 *       double Bij
	 *    for each atom pair:
	 *       double radius1
	 *    for each atom pair:
	 *       double lambda1
	 *    for each atom pair:
	 *       double alpha1
	 *    for each atom pair:
	 *       double radius2
	 *    for each atom pair:
	 *       double lambda2
	 *    for each atom pair:
	 *       double alpha2
	 */
	private CUBuffer<ByteBuffer> data;
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<DoubleBuffer> xin;
	private CUBuffer<DoubleBuffer> out;
	
	private static final int HeaderBytes = Integer.BYTES*4 + Double.BYTES*2;
	private static final int DihedralBytes = Long.BYTES + Integer.BYTES*4 + Short.BYTES*4 + Double.BYTES*2;
	private static final int ResPairBytes = Long.BYTES*3 + Double.BYTES*2;
	private static final int AtomPairBytes = Short.BYTES*2 + Integer.BYTES + Double.BYTES*9;
	
	
	public ResidueCudaCCDMinimizer(GpuStreamPool streams, ObjectiveFunction f) {
		super(streams.checkout(), "residueCcd");
		
		this.streams = streams;
		GpuStream stream = getStream();
		
		// get the molecule objective function
		if (f instanceof MoleculeObjectiveFunction) {
			mof = (MoleculeObjectiveFunction)f;
		} else {
			throw new Error("objective function should be a " + MoleculeObjectiveFunction.class.getSimpleName() + ", not a " + f.getClass().getSimpleName() + ". this is a bug");
		}
		
		// get the energy function
		if (mof.efunc instanceof ResidueForcefieldEnergy) {
			efunc = (ResidueForcefieldEnergy)mof.efunc;
		} else {
			throw new Error("energy function should be a " + ResidueForcefieldEnergy.class.getSimpleName() + ", not a " + mof.efunc.getClass().getSimpleName() + ". this is a bug");
		}

		// no need to waste time on broken conformations
		if (efunc.isBroken) {
			return;
		}
		
		// make sure this thread can use the cuda context
		stream.getContext().attachCurrentThread();
		
		// count atoms and offsets for each residue
		int[] atomOffsetsByResIndex = new int[efunc.residues.size()];
		Arrays.fill(atomOffsetsByResIndex, -1);
		int atomOffset = 0;
		int numAtoms = 0;
		for (int i=0; i<efunc.residues.size(); i++) {
			Residue res = efunc.residues.get(i);
			atomOffsetsByResIndex[i] = atomOffset;
			atomOffset += 3*res.atoms.size();
			numAtoms += res.atoms.size();
		}
		
		// make the coords buffer
		coords = stream.doubleBuffers.checkout(numAtoms*3);
		
		// collect the dihedral angles
		dihedrals = new ArrayList<>();
		ObjectiveFunction.DofBounds dofBounds = new ObjectiveFunction.DofBounds(mof.getConstraints());
		Map<Residue,int[]> cache = new IdentityHashMap<>();
		int maxNumAtoms = 0;
		for (int d=0; d<mof.getNumDOFs(); d++) {
			
			// make sure the DoF is a FreeDihedral. that's all we support on the gpu at the moment
			DegreeOfFreedom dof = mof.pmol.dofs.get(d);
			if (!(dof instanceof FreeDihedral)) {
				throw new Error("degree-of-freedom type " + dof.getClass().getSimpleName() + " not yet supported by CCD kernel."
					+ " Use CPU minimizer with GPU energy function instead");
			}
			
			// make the dihedral
			Dihedral dihedral = new Dihedral(d, (FreeDihedral)dof);
			
			// find the residue, if any
			dihedral.resIndex = efunc.residues.findIndex(dihedral.res);
			if (dihedral.resIndex < 0) {
				continue;
			}
			
			// get the residue pairs
			dihedral.resPairIndices = cache.get(dihedral.res);
			if (dihedral.resPairIndices == null) {
				dihedral.resPairIndices = efunc.makeResPairIndicesSubset(dof.getResidue());
				cache.put(dihedral.res, dihedral.resPairIndices);
			}
			
			dihedral.xdmin = dofBounds.getMin(d);
			dihedral.xdmax = dofBounds.getMax(d);
			
			maxNumAtoms = Math.max(maxNumAtoms, dihedral.res.atoms.size());
			
			dihedrals.add(dihedral);
		}
		
		// allocate the data buffer
		int numRotatedAtoms = 0;
		int numResPairs = 0;
		for (Dihedral dihedral : dihedrals) {
			numRotatedAtoms += MathTools.roundUpToMultiple(dihedral.rotatedIndices.length, Long.BYTES/Short.BYTES);
			numResPairs += MathTools.roundUpToMultiple(dihedral.resPairIndices.length, Long.BYTES/Integer.BYTES);
		}
		int numAtomPairs = 0;
		for (ResPair resPair : efunc.resPairs) {
			numAtomPairs += resPair.info.numAtomPairs;
		}
		data = stream.byteBuffers.checkout(
			HeaderBytes
			+ Long.BYTES*dihedrals.size()
			+ Long.BYTES*efunc.resPairs.length
			+ DihedralBytes*dihedrals.size()
				+ Short.BYTES*numRotatedAtoms
				+ Integer.BYTES*numResPairs
			+ ResPairBytes*efunc.resPairs.length
			+ AtomPairBytes*numAtomPairs
		);
		ByteBuffer databuf = data.getHostBuffer();
		
		ForcefieldParams ffparams = efunc.resPairCache.ffparams;
		
		// put the data header
		int flags = ffparams.hElect ? 1 : 0;
		flags <<= 1;
		flags |= ffparams.hVDW ? 1 : 0;
		flags <<= 1;
		flags |= ffparams.distDepDielect ? 1 : 0;
		flags <<= 1;
		flags |= ffparams.solvationForcefield == SolvationForcefield.EEF1 ? 1 : 0;
		databuf.putInt(flags);
		databuf.putInt(dihedrals.size());
		databuf.putInt(efunc.resPairs.length);
		databuf.putInt(maxNumAtoms);
		double coulombFactor = ForcefieldParams.coulombConstant/ffparams.dielectric;
		double scaledCoulombFactor = coulombFactor*ffparams.forcefld.coulombScaling;
		databuf.putDouble(coulombFactor);
		databuf.putDouble(scaledCoulombFactor);
		
		// put space for the dihedral offsets
		int dihedralOffsetsPos = databuf.position();
		for (int d=0; d<dihedrals.size(); d++) {
			databuf.putLong(0);
		}
		
		// put space for the res pair offsets
		int resPairsOffsetsPos = databuf.position();
		for (int d=0; d<efunc.resPairs.length; d++) {
			databuf.putLong(0);
		}
		
		// put the dihedrals
		for (int d=0; d<dihedrals.size(); d++) {
			Dihedral dihedral = dihedrals.get(d);
			
			// update the offset
			databuf.putLong(dihedralOffsetsPos + d*Long.BYTES, databuf.position());
			
			// put the dihedral header
			databuf.putInt(dihedral.resIndex);
			databuf.putInt(dihedral.res.atoms.size());
			databuf.putLong(atomOffsetsByResIndex[dihedral.resIndex]);
			for (int j=0; j<dihedral.dihedralIndices.length; j++) {
				databuf.putShort((short)(dihedral.dihedralIndices[j]*3));
			}
			databuf.putInt(dihedral.rotatedIndices.length);
			databuf.putInt(dihedral.resPairIndices.length);
			databuf.putDouble(Math.toRadians(dihedral.xdmin));
			databuf.putDouble(Math.toRadians(dihedral.xdmax));
			
			// put the rotated offsets (pad up to word boundary)
			int n = MathTools.roundUpToMultiple(dihedral.rotatedIndices.length, Long.BYTES/Short.BYTES);
			for (int i=0; i<n; i++) {
				if (i < dihedral.rotatedIndices.length) {
					databuf.putShort((short)(dihedral.rotatedIndices[i]*3)); 
				} else {
					databuf.putShort((short)0);
				}
			}
			
			// put the res pairs (pad up to word boundary)
			n = MathTools.roundUpToMultiple(dihedral.resPairIndices.length, Long.BYTES/Integer.BYTES);
			for (int i=0; i<n; i++) {
				if (i < dihedral.resPairIndices.length) {
					databuf.putInt(dihedral.resPairIndices[i]);
				} else {
					databuf.putInt(0);
				}
			}
		}
		
		// put the res pairs and atom pairs
		for (int i=0; i<efunc.resPairs.length; i++) {
			ResPair resPair = efunc.resPairs[i];
			
			// update the offset
			databuf.putLong(resPairsOffsetsPos + i*Long.BYTES, databuf.position());
			
			// put the res pair
			databuf.putLong(resPair.info.numAtomPairs);
			databuf.putLong(atomOffsetsByResIndex[efunc.residues.findIndexOrThrow(resPair.res1)]);
			databuf.putLong(atomOffsetsByResIndex[efunc.residues.findIndexOrThrow(resPair.res2)]);
			databuf.putDouble(resPair.weight);
			databuf.putDouble(resPair.offset + resPair.solvEnergy);
			
			// put the atom pairs
			// NOTE: use struct-of-arrays here, not array-of-structs
			// so the GPU can coalesce memory accesses
			for (int j=0; j<resPair.info.numAtomPairs; j++) {
				databuf.putLong(resPair.info.flags[j]);
			}
			for (int k=0; k<resPair.info.numPrecomputedPerAtomPair; k++) {
				for (int j=0; j<resPair.info.numAtomPairs; j++) {
					databuf.putDouble(resPair.info.precomputed[j*resPair.info.numPrecomputedPerAtomPair + k]);
				}
			}
		}
		
		databuf.flip();
		data.uploadAsync();

		// allocate the x input
		xin = stream.doubleBuffers.checkout(dihedrals.size());
		
		// allocate the output
		out = stream.doubleBuffers.checkout(dihedrals.size() + 1);
		
		func = makeFunction("ccd");
		func.numBlocks = 1;
		func.sharedMemCalc = (int blockThreads) -> {
			return blockThreads*Double.BYTES // energy reduction
				+ dihedrals.size()*Double.BYTES // nextx
				+ dihedrals.size()*Double.BYTES // firstSteps
				+ dihedrals.size()*Double.BYTES; // lastSteps
		};
		func.setArgs(Pointer.to(
			data.getDevicePointer(),
			coords.getDevicePointer(),
			xin.getDevicePointer(),
			out.getDevicePointer()
		));
		
		// calc the best number of block threads for this machine
		func.blockThreads = func.getBestBlockThreads(blockThreads);
	}
	
	@Override
	public Minimizer.Result minimizeFromCenter() {

		DoubleMatrix1D x = DoubleFactory1D.dense.make(dihedrals.size());
		for (Dihedral dihedral : dihedrals) {
			x.set(dihedral.d, (dihedral.xdmax + dihedral.xdmin)/2);
		}

		return minimizeFrom(x);
	}

	@Override
	public Minimizer.Result minimizeFrom(DoubleMatrix1D x) {

		// no need to waste time on broken conformations
		if (efunc.isBroken) {
			return new Minimizer.Result(null, Double.POSITIVE_INFINITY);
		}
		
		// put the molecule in the staring state
		// TODO: do initial pose on the GPU?
		for (Dihedral dihedral : dihedrals) {
			mof.setDOF(dihedral.d, x.get(dihedral.d));
		}
		
		// pack and upload the coords
		DoubleBuffer coordsbuf = coords.getHostBuffer();
		coordsbuf.clear();
		for (Residue res : efunc.residues) {
			coordsbuf.put(res.coords);
		}
		coordsbuf.clear();
		coords.uploadAsync();

		// upload x
		DoubleBuffer xinbuf = xin.getHostBuffer();
		xinbuf.clear();
		for (Dihedral dihedral : dihedrals) {
			xinbuf.put(Math.toRadians(x.get(dihedral.d)));
		}
		xinbuf.clear();
		xin.uploadAsync();
		
		// launch the kernel and wait
		func.runAsync();
		DoubleBuffer outbuf = out.downloadSync();
		
		// read the result
		outbuf.rewind();
		x = DoubleFactory1D.dense.make(mof.getNumDOFs());
		for (Dihedral dihedral : dihedrals) {
			x.set(dihedral.d, Math.toDegrees(outbuf.get()));
		}
		Minimizer.Result result = new Minimizer.Result(x, outbuf.get());
		
		// update the CPU-side molecule
		mof.setDOFs(result.dofValues);
		// TODO: or download the coords, test which is faster?
		
		return result;
	}

	@Override
	public void clean() {
		GpuStream stream = getStream();
		if (data != null) {
			stream.byteBuffers.release(data);
			data = null;
		}
		if (coords != null) {
			stream.doubleBuffers.release(coords);
			coords = null;
		}
		if (xin != null) {
			stream.doubleBuffers.release(xin);
			xin = null;
		}
		if (out != null) {
			stream.doubleBuffers.release(out);
			out = null;
		}
		streams.release(stream);
	}
}
