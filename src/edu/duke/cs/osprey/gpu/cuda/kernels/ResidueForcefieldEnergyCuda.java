package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvPairParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.VdwParams;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairs;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class ResidueForcefieldEnergyCuda extends Kernel implements EnergyFunction.DecomposableByDof, EnergyFunction.NeedsCleanup {
	
	private static final long serialVersionUID = 4015880661919715967L;
	
	public static class ResPair {
		
		public final int resIndex1;
		public final Residue res1;
		public final int resIndex2;
		public final Residue res2;
		public final ResidueInteractions.Pair pair;
		public final AtomPairs atomPairs;
		public final int numAtomPairs;
		
		public ResPair(Residue[] ress, ResidueInteractions.Pair pair, AtomConnectivity connectivity) {
			
			resIndex1 = findIndex(ress, pair.resNum1);
			res1 = ress[resIndex1];
			resIndex2 = findIndex(ress, pair.resNum2);
			res2 = ress[resIndex2];
			this.pair = pair;
			
			atomPairs = connectivity.getAtomPairs(res1, res2);
			numAtomPairs = atomPairs.getNumPairs(AtomNeighbors.Type.BONDED14)
				+ atomPairs.getNumPairs(AtomNeighbors.Type.NONBONDED);
		}
	}
	
	private static Residue findRes(Collection<Residue> residues, String resNum) {
		for (Residue res : residues) {
			if (res.getPDBResNumber().equals(resNum)) {
				return res;
			}
		}
		throw new NoSuchElementException("no residue " + resNum + " found in " + residues);
	}
	
	private static int findIndex(Residue[] ress, String resNum) {
		for (int i=0; i<ress.length; i++) {
			if (ress[i].getPDBResNumber().equals(resNum)) {
				return i;
			}
		}
		throw new NoSuchElementException("no residue " + resNum + " found in " + Arrays.toString(ress));
	}
	
	public final ForcefieldParams params;
	public final ResidueInteractions inters;
	public final AtomConnectivity connectivity;
	
	private boolean isBroken;
	
	/* buffer layout:
	 * NOTE: try to use 8-byte alignments to be happy on 64-bit machines
	 * 
	 * long flags (useHEs, useHvdW, distDepDielect, useEEF1)
	 * double coulombFactor
	 * double scaledCoulombFactor
	 * 
	 * for each residue pair:
	 *    long offset
	 * 
	 * for each residue pair:
	 *    long numAtomPairs
	 *    long atomsOffset1
	 *    long atomsOffset2
	 *    
	 *    for each atom pair:
	 *       short atomOffset1
	 *       short atomOffset2
	 *       int flags  (isHeavyPair, is14Bonded)
	 *       double charge
	 *       double Aij
	 *       double Bij
	 *       double radius1
	 *       double lambda1
	 *       double alpha1
	 *       double radius2
	 *       double lambda2
	 *       double alpha2
	 *       
	 *    double weight
	 *    double offset
	 */
	private ByteBuffer data;
	private DoubleBuffer coords;
	private IntBuffer allIndices;
	
	private static final int HeaderBytes = Long.BYTES + Double.BYTES*2;
	private static final int ResPairBytes = Long.BYTES*3 + Double.BYTES*2;
	private static final int AtomPairBytes = Short.BYTES*2 + Integer.BYTES + Double.BYTES*9;
	
	private Residue[] ress;
	private ResPair[] resPairs;
	
	public ResidueForcefieldEnergyCuda(GpuStream stream, ForcefieldParams params, ResidueInteractions inters, Molecule mol, AtomConnectivity connectivity)
	throws IOException {
		this(stream, params, inters, mol.residues, connectivity);
	}
	
	public ResidueForcefieldEnergyCuda(GpuStream stream, ForcefieldParams params, ResidueInteractions inters, Collection<Residue> residues, AtomConnectivity connectivity)
	throws IOException {
		super(stream, "residueForcefield");
		
		this.params = params;
		this.inters = inters;
		this.connectivity = connectivity;
		
		// map the residue numbers to residues
		ress = new Residue[inters.getResidueNumbers().size()];
		int index = 0;
		for (String resNum : inters.getResidueNumbers()) {
			ress[index++] = findRes(residues, resNum);
		}
		resPairs = new ResPair[inters.size()];
		index = 0;
		for (ResidueInteractions.Pair pair : inters) {
			resPairs[index++] = new ResPair(ress, pair, connectivity);
		}
		
		// is this a broken conformation?
		isBroken = false;
		for (ResPair pair : resPairs) {
			if (pair.res1.confProblems.size() + pair.res2.confProblems.size() > 0) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		
		// count atoms and offsets for each residue
		int[] atomOffsetsByResIndex = new int[ress.length];
		Arrays.fill(atomOffsetsByResIndex, -1);
		int atomOffset = 0;
		int numAtoms = 0;
		for (int i=0; i<ress.length; i++) {
			Residue res = ress[i];
			atomOffsetsByResIndex[i] = atomOffset;
			atomOffset += 3*res.atoms.size();
			numAtoms += res.atoms.size();
		}
		
		// make the coords buffer
		coords = BufferTools.Type.Direct.makeDouble(numAtoms*3);
		
		// allocate the data
		int totalNumAtomPairs = 0;
		for (int i=0; i<resPairs.length; i++) {
			totalNumAtomPairs += resPairs[i].numAtomPairs;
		}
		data = BufferTools.Type.Direct.makeByte(HeaderBytes + (Long.BYTES + ResPairBytes)*resPairs.length + AtomPairBytes*totalNumAtomPairs);
		
		// put the data header
		long flags = params.hElect ? 1 : 0;
		flags <<= 1;
		flags |= params.hVDW ? 1 : 0;
		flags <<= 1;
		flags |= params.distDepDielect ? 1 : 0;
		flags <<= 1;
		flags |= params.solvationForcefield == SolvationForcefield.EEF1 ? 1 : 0;
		data.putLong(flags);
		double coulombFactor = ForcefieldParams.coulombConstant/params.dielectric;
		double scaledCoulombFactor = coulombFactor*params.forcefld.coulombScaling;
		data.putDouble(coulombFactor);
		data.putDouble(scaledCoulombFactor);
		
		// put the res pair offsets
		long offset = HeaderBytes + Long.BYTES*resPairs.length;
		for (int i=0; i<resPairs.length; i++) {
			data.putLong(offset);
			offset += ResPairBytes + AtomPairBytes*resPairs[i].numAtomPairs;
		}
		
		VdwParams vdwparams = new VdwParams();
		SolvPairParams solvparams = new SolvPairParams();
		
		// put the res pairs and atom pairs
		for (int i=0; i<resPairs.length; i++) {
			ResPair resPair = resPairs[i];
			
			// put the res pair
			data.putLong(resPair.numAtomPairs);
			data.putLong(atomOffsetsByResIndex[resPair.resIndex1]);
			data.putLong(atomOffsetsByResIndex[resPair.resIndex2]);
			
			for (AtomNeighbors.Type type : Arrays.asList(AtomNeighbors.Type.BONDED14, AtomNeighbors.Type.NONBONDED)) {
				for (int[] atomPair : resPair.atomPairs.getPairs(type)) {
			
					// put the atom offsets
					data.putShort((short)(atomPair[0]*3));
					data.putShort((short)(atomPair[1]*3));
					
					Atom atom1 = resPair.res1.atoms.get(atomPair[0]);
					Atom atom2 = resPair.res2.atoms.get(atomPair[1]);
					
					// pack the flags
					boolean isHeavyPair = !atom1.isHydrogen() && !atom2.isHydrogen();
					boolean is14Bonded = type == AtomNeighbors.Type.BONDED14;
					int atomFlags = is14Bonded ? 1 : 0;
					atomFlags <<= 1;
					atomFlags |= isHeavyPair ? 1 : 0;
					data.putInt(atomFlags);
					
					// calc electrostatics params
					data.putDouble(atom1.charge*atom2.charge);
					
					// calc vdw params
					params.getVdwParams(atom1, atom2, type, vdwparams);
					data.putDouble(vdwparams.Aij);
					data.putDouble(vdwparams.Bij);
					
					// compute solvation params if needed
					if (isHeavyPair) {
						switch (params.solvationForcefield) {
							
							case EEF1:
								params.eef1parms.getSolvationPairParams(atom1, atom2, params.solvScale, solvparams);
								data.putDouble(solvparams.radius1);
								data.putDouble(solvparams.lambda1);
								data.putDouble(solvparams.alpha1);
								data.putDouble(solvparams.radius2);
								data.putDouble(solvparams.lambda2);
								data.putDouble(solvparams.alpha2);
							break;
							
							default:
								data.position(data.position() + Double.BYTES*6);
						}
							
					} else {
						
						data.position(data.position() + Double.BYTES*6);
					}
				}
			}
			
			// put the res pair weight and offset
			data.putDouble(resPair.pair.weight);
			
			double resPairOffset = resPair.pair.offset;
			if (resPair.res1 == resPair.res2) {
				// update the pair offset with internal solvation energy
				switch (params.solvationForcefield) {
					
					case EEF1:
						resPairOffset += params.eef1parms.getInternalEnergy(resPair.res1)*params.solvScale*resPair.pair.weight;
					break;
					
					default:
						// do nothing
				}
			}
			data.putDouble(resPairOffset);
		}
		
		assert (data.position() == data.capacity());
		
		data.flip();
		
		// make the indices
		allIndices = BufferTools.Type.Direct.makeInt(resPairs.length);
		for (int i=0; i<resPairs.length; i++) {
			allIndices.put(i);
		}
		allIndices.flip();
	}
	
	@Override
	public void cleanup() {
		// TODO
	}
	
	public void copyCoords() {
		coords.clear();
		for (Residue res : ress) {
			coords.put(res.coords);
		}
		coords.clear();
	}

	@Override
	public double getEnergy() {
		
		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}
		
		// capture the current molecule state
		copyCoords();
		
		// TODO: launch CUDA kernel
		return getEnergy(allIndices);
	}
	
	// TEMP: just for cpu-side testing
	private double getEnergy(IntBuffer indices) {
		
		ByteBuffer data = this.data;
		DoubleBuffer coords = this.coords;
		
		// read the data header
		data.rewind();
		long flags = data.getLong();
		boolean useEEF1 = (flags & 0x1) == 0x1;
		flags >>= 1;
		boolean distDepDielect = (flags & 0x1) == 0x1;
		flags >>= 1;
		boolean useHvdW = (flags & 0x1) == 0x1;
		flags >>= 1;
		boolean useHEs = (flags & 0x1) == 0x1;
		double coulombFactor = data.getDouble();
		double scaledCoulombFactor = data.getDouble();
		
		double energy = 0;
		
		indices.rewind();
		int numIndices = indices.capacity();
		for (int i=0; i<numIndices; i++) {
			
			// seek to the res pair
			data.position(HeaderBytes + i*Long.BYTES);
			data.position((int)data.getLong());
			
			long numAtomPairs = data.getLong();
			int atomsOffset1 = (int)data.getLong();
			int atomsOffset2 = (int)data.getLong();
			
			double resPairEnergy = 0;
			
			// for each atom pair...
			for (int j=0; j<numAtomPairs; j++) {
				
				int atomOffset1 = data.getShort();
				int atomOffset2 = data.getShort();
				
				// get the radius
				double r2;
				double r;
				{
					// read atom coords
					coords.position(atomsOffset1 + atomOffset1);
					double x1 = coords.get();
					double y1 = coords.get();
					double z1 = coords.get();
					coords.position(atomsOffset2 + atomOffset2);
					double x2 = coords.get();
					double y2 = coords.get();
					double z2 = coords.get();
					
					double d;
					
					d = x1 - x2;
					r2 = d*d;
					d = y1 - y2;
					r2 += d*d;
					d = z1 - z2;
					r2 += d*d;
					r = Math.sqrt(r2);
				}
				
				// read the flags
				// NOTE: this is efficient, but destructive to the val
				flags = data.getInt();
				boolean isHeavyPair = (flags & 0x1) == 0x1;
				flags >>= 1;
				boolean is14Bonded = (flags & 0x1) == 0x1;
				
				// electrostatics
				if (isHeavyPair || useHEs) {
					double charge = data.getDouble();
					if (is14Bonded) {
						if (distDepDielect) {
							resPairEnergy += scaledCoulombFactor*charge/r2;
						} else {
							resPairEnergy += scaledCoulombFactor*charge/r;
						}
					} else {
						if (distDepDielect) {
							resPairEnergy += coulombFactor*charge/r2;
						} else {
							resPairEnergy += coulombFactor*charge/r;
						}
					}
				} else {
					data.position(data.position() + Double.BYTES);
				}
				
				// van der Waals
				if (isHeavyPair || useHvdW) {
					
					double Aij = data.getDouble();
					double Bij = data.getDouble();
					
					// compute vdw
					double r6 = r2*r2*r2;
					double r12 = r6*r6;
					resPairEnergy += Aij/r12 - Bij/r6;
					
				} else {
					data.position(data.position() + Double.BYTES*2);
				}
				
				// solvation
				if (useEEF1 && isHeavyPair && r2 < ForcefieldParams.solvCutoff2) {
							
					double radius1 = data.getDouble();
					double lambda1 = data.getDouble();
					double alpha1 = data.getDouble();
					double radius2 = data.getDouble();
					double lambda2 = data.getDouble();
					double alpha2 = data.getDouble();
					
					// compute solvation energy
					double Xij = (r - radius1)/lambda1;
					double Xji = (r - radius2)/lambda2;
					resPairEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
					
				} else {
					data.position(data.position() + Double.BYTES*6);
				}
			}
			
			// apply weights and offsets
			double weight = data.getDouble();
			double offset = data.getDouble();
			energy += resPairEnergy*weight;
			energy += offset;
		}
		
		return energy;
	}
	
	private class Subset implements EnergyFunction {
		
		private static final long serialVersionUID = -1749739381007657718L;
		
		private IntBuffer indices;
		
		public Subset(Residue res) {
			
			// pass 1: count
			int num = 0;
			for (int i=0; i<resPairs.length; i++) {
				ResPair resPair = resPairs[i];
				if (resPair.res1 == res || resPair.res2 == res) {
					num++;
				}
			}
			
			// pass 2: collect
			indices = BufferTools.Type.Direct.makeInt(num);
			for (int i=0; i<resPairs.length; i++) {
				ResPair resPair = resPairs[i];
				if (resPair.res1 == res || resPair.res2 == res) {
					indices.put(i);
				}
			}
		}
		
		@Override
		public double getEnergy() {
			return ResidueForcefieldEnergyCuda.this.getEnergy(indices);
		}
	}
	
	@Override
	public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {
		
		Map<Residue,Subset> cache = new HashMap<>();
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		for (DegreeOfFreedom dof : dofs) {
			Residue res = dof.getResidue();
			
			if (res == null) {
				
				// no res, just use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// make a subset energy function
				Subset subset = cache.get(res);
				if (subset == null) {
					subset = new Subset(res);
					cache.put(res, subset);
				}
				efuncs.add(subset);
			}
		}
		
		return efuncs;
	}
}
