package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.structure.Residue;
import jcuda.Pointer;

public class ForcefieldKernelOneBlockCuda extends Kernel {
	
	private class DofInfo {
		
		public Residue res;
		public BigForcefieldEnergy.Subset subset;
		public int[] dihedralIndices;
		public List<Integer> rotatedIndices;
		public int atomOffset;
		public int numAtoms;
		
		public DofInfo(FreeDihedral dof) {
			
			this.res = dof.getResidue();
			this.subset = ffenergy.new Subset(ffenergy.getInteractions().makeSubsetByResidue(res));
			
			this.dihedralIndices = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
			this.rotatedIndices = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
			this.atomOffset = ffenergy.getAtomOffset(res);
			this.numAtoms = res.atoms.size();
		}
	}
	
	public static class LinesearchResult {
		
		public double dihedralRadians;
		public double energy;
		
		public LinesearchResult() {
			dihedralRadians = 0;
			energy = 0;
		}
		
		public LinesearchResult(LinesearchResult other) {
			this.dihedralRadians = other.dihedralRadians;
			this.energy = other.energy;
		}
		
		public void set(DoubleBuffer buf, double internalSolvationEnergy) {
			dihedralRadians = buf.get(0);
			energy = buf.get(1) + internalSolvationEnergy;
		}
	}
	
	private Kernel.Function funcFull;
	private Kernel.Function funcDof;
	private Kernel.Function funcPoseDof;
	private Kernel.Function funcLinesearch;
	
	private CUBuffer<DoubleBuffer> coords;
	private CUBuffer<IntBuffer> atomFlags;
	private CUBuffer<DoubleBuffer> precomputed;
	private CUBuffer<ByteBuffer> ffargs;
	
	private CUBuffer<ByteBuffer> dofargs;
	private CUBuffer<IntBuffer> subsetTables;
	private CUBuffer<IntBuffer> dihedralIndices;
	private CUBuffer<IntBuffer> rotatedIndices;
	
	private CUBuffer<DoubleBuffer> energies;
	private CUBuffer<DoubleBuffer> linesearchOut;
	private LinesearchResult linesearchResult;
	
	private List<DofInfo> dofInfos;
	
	private int[] dofArg;
	private double[] xdArg;
	private double[] xdminArg;
	private double[] xdmaxArg;
	private double[] stepArg;
	
	private int blockThreads;
	private int maxNumBlocks;
	
	private BigForcefieldEnergy ffenergy;
	
	public ForcefieldKernelOneBlockCuda(GpuStream stream, BigForcefieldEnergy ffenergy, List<FreeDihedral> dofs)
	throws IOException {
		super(stream, "forcefieldOneBlock");
		
		// OPTIMIZATION: 1024 threads seems to work best on modern cards, including GTX 1070 and Tesla K80
		blockThreads = 512;
		
		// TODO: make arg
		maxNumBlocks = 1;
		
		this.ffenergy = ffenergy;
		
		// wrap the incoming buffers
		coords = stream.makeBuffer(ffenergy.getCoords());
		atomFlags = stream.makeBuffer(ffenergy.getAtomFlags());
		precomputed = stream.makeBuffer(ffenergy.getPrecomputed());
		
		// make the args buffer
		ffargs = stream.makeByteBuffer(40);
		ByteBuffer argsBuf = ffargs.getHostBuffer();
		argsBuf.rewind();
		argsBuf.putInt(ffenergy.getFullSubset().getNumAtomPairs());
		argsBuf.putInt(ffenergy.getFullSubset().getNum14AtomPairs());
		argsBuf.putDouble(ffenergy.getCoulombFactor());
		argsBuf.putDouble(ffenergy.getScaledCoulombFactor());
		argsBuf.putDouble(ffenergy.getSolvationCutoff2());
		argsBuf.put((byte)(ffenergy.useDistDependentDielectric() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHElectrostatics() ? 1 : 0));
		argsBuf.put((byte)(ffenergy.useHVdw() ? 1 : 0));
		
		// make the energy buffer for the full subset
		int allPairs = ffenergy.getFullSubset().getNum14AtomPairs();
		energies = stream.makeDoubleBuffer(calcNumBlocks(allPairs));
		
		// upload static forcefield info
		atomFlags.uploadAsync();
		precomputed.uploadAsync();
		ffargs.uploadAsync();
		
		// get info about the dofs
		dofInfos = new ArrayList<>();
		int subsetsSize = 0;
		int numRotatedAtoms = 0;
		for (FreeDihedral dof : dofs) {
			
			DofInfo dofInfo = new DofInfo(dof);
			dofInfos.add(dofInfo);
			
			// update counts
			subsetsSize += dofInfo.subset.getNumAtomPairs();
			numRotatedAtoms += dofInfo.rotatedIndices.size();
		}
		
		// allocate dof-related buffers
		subsetTables = stream.makeIntBuffer(subsetsSize);
		IntBuffer subsetTablesBuf = subsetTables.getHostBuffer();
		subsetTablesBuf.clear();
		
		dofargs = stream.makeByteBuffer(dofs.size()*32);
		ByteBuffer dofargsBuf = dofargs.getHostBuffer();
		dofargsBuf.clear();
		
		dihedralIndices = stream.makeIntBuffer(dofs.size()*4);
		IntBuffer dihedralIndicesBuf = dihedralIndices.getHostBuffer();
		dihedralIndicesBuf.clear();
		
		rotatedIndices = stream.makeIntBuffer(numRotatedAtoms);
		IntBuffer rotatedIndicesBuf = rotatedIndices.getHostBuffer();
		rotatedIndicesBuf.clear();
		
		// populate dof buffers
		for (int d=0; d<dofs.size(); d++) {
			DofInfo dofInfo = dofInfos.get(d);
			
			int firstCoord = dofInfo.atomOffset*3;
			int lastCoord = (dofInfo.atomOffset + dofInfo.numAtoms)*3 - 1;
			
			// update dofargs
			dofargsBuf.putInt(subsetTablesBuf.position());
			dofargsBuf.putInt(dofInfo.subset.getNumAtomPairs());
			dofargsBuf.putInt(dofInfo.subset.getNum14AtomPairs());
			dofargsBuf.putInt(rotatedIndicesBuf.position());
			dofargsBuf.putInt(dofInfo.rotatedIndices.size());
			dofargsBuf.putInt(firstCoord);
			dofargsBuf.putInt(lastCoord);
			dofargsBuf.putInt(0); // spacer for word alignment
			
			// append subset table
			dofInfo.subset.getSubsetTable().rewind();
			subsetTablesBuf.put(dofInfo.subset.getSubsetTable());
			
			// collect dihedral indices
			for (int i=0; i<dofInfo.dihedralIndices.length; i++) {
				dihedralIndicesBuf.put(dofInfo.dihedralIndices[i]);
			}
			
			// collect rotated indices
			for (int i=0; i<dofInfo.rotatedIndices.size(); i++) {
				rotatedIndicesBuf.put(dofInfo.rotatedIndices.get(i));
			}
		}
		
		// upload more bufs
		dofargs.uploadAsync();
		subsetTables.uploadAsync();
		dihedralIndices.uploadAsync();
		rotatedIndices.uploadAsync();
		
		// init other kernel function inputs,outputs
		dofArg = new int[1];
		xdArg = new double[1];
		xdminArg = new double[1];
		xdmaxArg = new double[1];
		stepArg = new double[1];
		
		linesearchOut = stream.makeDoubleBuffer(2);
		linesearchResult = new LinesearchResult();
		
		// init the kernel functions
		
		funcFull = makeFunction("calcEnergy");
		funcFull.numBlocks = calcNumBlocks(allPairs);
		funcFull.blockThreads = blockThreads;
		funcFull.sharedMemBytes = blockThreads*Double.BYTES;
		funcFull.setArgs(Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			ffargs.makeDevicePointer(),
			energies.makeDevicePointer()
		));
		
		funcDof = makeFunction("calcDofEnergy");
		funcDof.setArgs(Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			ffargs.makeDevicePointer(),
			subsetTables.makeDevicePointer(),
			dofargs.makeDevicePointer(),
			Pointer.to(dofArg),
			energies.makeDevicePointer()
		));
		
		funcPoseDof = makeFunction("poseAndCalcDofEnergy");
		funcPoseDof.setArgs(Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			ffargs.makeDevicePointer(),
			subsetTables.makeDevicePointer(),
			dihedralIndices.makeDevicePointer(),
			rotatedIndices.makeDevicePointer(),
			dofargs.makeDevicePointer(),
			Pointer.to(dofArg),
			Pointer.to(xdArg),
			energies.makeDevicePointer()
		));
		
		funcLinesearch = makeFunction("linesearch");
		funcLinesearch.numBlocks = 1;
		funcLinesearch.blockThreads = blockThreads;
		funcLinesearch.setArgs(Pointer.to(
			coords.makeDevicePointer(),
			atomFlags.makeDevicePointer(),
			precomputed.makeDevicePointer(),
			ffargs.makeDevicePointer(),
			subsetTables.makeDevicePointer(),
			dihedralIndices.makeDevicePointer(),
			rotatedIndices.makeDevicePointer(),
			dofargs.makeDevicePointer(),
			Pointer.to(dofArg),
			Pointer.to(xdArg),
			Pointer.to(xdminArg),
			Pointer.to(xdmaxArg),
			Pointer.to(stepArg),
			linesearchOut.makeDevicePointer()
		));
	}
	
	public void uploadCoordsAsync() {
		
		// tell the forcefield to gather updated coords
		ffenergy.updateCoords();
		
		coords.uploadAsync();
	}
	
	public double calcEnergySync() {
		funcFull.runAsync();
		return downloadEnergySync(funcFull.numBlocks) + ffenergy.getFullSubset().getInternalSolvationEnergy();
	}
	
	public double calcEnergyDofSync(int d) {
		DofInfo dofInfo = dofInfos.get(d);
		funcDof.numBlocks = calcNumBlocks(dofInfo.subset.getNumAtomPairs());
		funcDof.blockThreads = blockThreads;
		funcDof.sharedMemBytes = blockThreads*Double.BYTES;
		dofArg[0] = d;
		funcDof.runAsync();
		return downloadEnergySync(funcDof.numBlocks) + dofInfo.subset.getInternalSolvationEnergy();
	}
	
	public double poseAndCalcEnergyDofSync(int d, double xdRadians) {
		DofInfo dofInfo = dofInfos.get(d);
		funcPoseDof.numBlocks = calcNumBlocks(dofInfo.subset.getNumAtomPairs());
		funcPoseDof.blockThreads = blockThreads;
		funcPoseDof.sharedMemBytes = blockThreads*Double.BYTES + dofInfo.numAtoms*3*Double.BYTES;
		dofArg[0] = d;
		xdArg[0] = xdRadians;
		funcPoseDof.runAsync();
		return downloadEnergySync(funcPoseDof.numBlocks) + dofInfo.subset.getInternalSolvationEnergy();
	}
	
	public LinesearchResult linesearchSync(int d, double xdRadians, double xdminRadians, double xdmaxRadians, double stepRadians) {
		DofInfo dofInfo = dofInfos.get(d);
		funcLinesearch.sharedMemBytes = blockThreads*Double.BYTES + dofInfo.numAtoms*3*Double.BYTES;
		dofArg[0] = d;
		xdArg[0] = xdRadians;
		xdminArg[0] = xdminRadians;
		xdmaxArg[0] = xdmaxRadians;
		stepArg[0] = stepRadians;
		funcLinesearch.runAsync();
		linesearchResult.set(linesearchOut.downloadSync(), dofInfo.subset.getInternalSolvationEnergy());
		return linesearchResult;
	}
	
	public void cleanup() {
		coords.cleanup();
		atomFlags.cleanup();
		precomputed.cleanup();
		ffargs.cleanup();
		dofargs.cleanup();
		subsetTables.cleanup();
		dihedralIndices.cleanup();
		rotatedIndices.cleanup();
		energies.cleanup();
		linesearchOut.cleanup();
	}
	
	private int calcNumBlocks(int numPairs) {
		return Math.min(maxNumBlocks, divUp(numPairs, blockThreads));
	}
	
	private double downloadEnergySync(int numBlocks) {
		
		DoubleBuffer buf = energies.downloadSync();
		buf.rewind();
		
		// do the final reduction across blocks on the cpu
		double energy = 0;
		for (int i=0; i<numBlocks; i++) {
			energy += buf.get();
		}
		return energy;
	}
}
