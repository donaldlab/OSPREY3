package edu.duke.cs.osprey.ematrix;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.Gpus;
import edu.duke.cs.osprey.minimization.CudaCCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.Minimizer.Result;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

@Deprecated
public abstract class SimpleEnergyCalculator {
	
	public final ForcefieldParams ffparams;
	public ConfSpace confSpace;
	public List<Residue> shellResidues;
	
	protected SimpleEnergyCalculator(ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues) {
		this.ffparams = ffparams;
		this.confSpace = confSpace;
		this.shellResidues = shellResidues;
	}
	
	public abstract EnergyFunctionGenerator getEnergyFunctionGenerator();
	
	public EnergyFunction makeSingleEfunc(int pos) {
		return makeSingleEfunc(pos, null);
	}
	
	public abstract EnergyFunction makeSingleEfunc(int pos, Molecule mol);
	
	public Minimizer.Result calcSingle(int pos, int rc) {
		return calcSingle(pos, rc, ParameterizedMoleculeCopy.makeNoCopy(confSpace));
	}
	
	public abstract Minimizer.Result calcSingle(int pos, int rc, ParameterizedMoleculeCopy pmol);
	
	public EnergyFunction makePairEfunc(int pos1, int pos2) {
		return makePairEfunc(pos1, pos2, null);
	}
	
	public abstract EnergyFunction makePairEfunc(int pos1, int pos2, Molecule mol);
	
	public Minimizer.Result calcPair(int pos1, int rc1, int pos2, int rc2) {
		return calcPair(pos1, rc1, pos2, rc2, ParameterizedMoleculeCopy.makeNoCopy(confSpace));
	}
	
	public abstract Minimizer.Result calcPair(int pos1, int rc1, int pos2, int rc2, ParameterizedMoleculeCopy pmol);
	
	
	public static class Cpu extends SimpleEnergyCalculator {

		private EnergyFunctionGenerator efuncs;
		
		public Cpu(ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues) {
			super(ffparams, confSpace, shellResidues);
			efuncs = new EnergyFunctionGenerator(ffparams);
		}
		
		@Override
		public EnergyFunctionGenerator getEnergyFunctionGenerator() {
			return efuncs;
		}
		
		@Override
		public EnergyFunction makeSingleEfunc(int pos, Molecule mol) {
			return efuncs.interactionEnergy(FFInterGen.makeIntraAndShell(confSpace, pos, shellResidues, mol));
		}

		@Override
		public EnergyFunction makePairEfunc(int pos1, int pos2, Molecule mol) {
			return efuncs.interactionEnergy(FFInterGen.makeResPair(confSpace, pos1, pos2, mol));
		}

		@Override
		public Result calcSingle(int pos, int rc, ParameterizedMoleculeCopy pmol) {
			return calc(makeSingleEfunc(pos, pmol.getCopiedMolecule()), new RCTuple(pos, rc), pmol);
		}

		@Override
		public Result calcPair(int pos1, int rc1, int pos2, int rc2, ParameterizedMoleculeCopy pmol) {
			return calc(makePairEfunc(pos1, pos2, pmol.getCopiedMolecule()), new RCTuple(pos1, rc1, pos2, rc2), pmol);
		}
		
		public Minimizer.Result calc(EnergyFunction efunc, RCTuple tuple, ParameterizedMoleculeCopy pmol) {
			
			// put molecule in correct conformation for rcs
			MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple, pmol);
			
			// optimize the degrees of freedom, if needed
			if (mof.getNumDOFs() > 0) {
				return new SimpleCCDMinimizer(mof).minimize();
			}

			// otherwise, just evaluate the energy function
			return new Minimizer.Result(null, efunc.getEnergy());
		}
	}
	
	
	/**
	 * This is pretty slow compared to the CPU,
	 * so don't actually use it in the real world.
	 * Maybe someday it could be faster...
	 * Use the multi-threaded CPU calculator instead
	 */
	@Deprecated
	public static class Cuda extends SimpleEnergyCalculator {
		
		public static boolean isSupported() {
			return !Gpus.get().getGpus().isEmpty();
		}
		
		private GpuStreamPool pool;
		private GpuEnergyFunctionGenerator efuncs;
		
		public Cuda(GpuStreamPool pool, ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues) {
			super(ffparams, confSpace, shellResidues);
			this.pool = pool;
			this.efuncs = new GpuEnergyFunctionGenerator(ffparams, pool);
		}
		
		private BufferTools.Type getBufType(boolean isMinimizing) {
			// if we're not minimizing, don't incur the overhead of directly-allocated buffers
			if (isMinimizing) {
				return BufferTools.Type.Direct;
			} else {
				return BufferTools.Type.Normal;
			}
		}
		
		@Override
		public EnergyFunctionGenerator getEnergyFunctionGenerator() {
			return efuncs;
		}

		@Override
		public BigForcefieldEnergy makeSingleEfunc(int pos, Molecule mol) {
			// this function is typically not called during minimization,
			// but rather by code that's inspecting subsets of the larger energy function,
			return makeSingleEfunc(pos, mol, getBufType(false));
		}
		
		public BigForcefieldEnergy makeSingleEfunc(int pos, Molecule mol, BufferTools.Type bufType) {
			ForcefieldInteractions ffinteractions = FFInterGen.makeIntraAndShell(confSpace, pos, shellResidues, mol);
			return new BigForcefieldEnergy(ffparams, ffinteractions, bufType);
		}

		@Override
		public BigForcefieldEnergy makePairEfunc(int pos1, int pos2, Molecule mol) {
			return makePairEfunc(pos1, pos2, mol, getBufType(false));
		}

		private BigForcefieldEnergy makePairEfunc(int pos1, int pos2, Molecule mol, BufferTools.Type bufType) {
			ForcefieldInteractions ffinteractions = FFInterGen.makeResPair(confSpace, pos1, pos2, mol);
			return new BigForcefieldEnergy(ffparams, ffinteractions, bufType);
		}

		@Override
		public Result calcSingle(int pos, int rc, ParameterizedMoleculeCopy pmol) {
			RCTuple tuple = new RCTuple(pos, rc);
			boolean isMinimizing = MoleculeModifierAndScorer.hasMinimizableDofs(confSpace, tuple);
			BigForcefieldEnergy efunc = makeSingleEfunc(pos, pmol.getCopiedMolecule(), getBufType(isMinimizing));
			return calc(efunc, tuple, pmol);
		}

		@Override
		public Result calcPair(int pos1, int rc1, int pos2, int rc2, ParameterizedMoleculeCopy pmol) {
			RCTuple tuple = new RCTuple(pos1, rc1, pos2, rc2);
			boolean isMinimizing = MoleculeModifierAndScorer.hasMinimizableDofs(confSpace, tuple);
			BigForcefieldEnergy efunc = makePairEfunc(pos1, pos2, pmol.getCopiedMolecule(), getBufType(isMinimizing));
			return calc(efunc, tuple, pmol);
		}
		
		public Minimizer.Result calc(EnergyFunction efunc, RCTuple tuple, ParameterizedMoleculeCopy pmol) {
			
			// put molecule in correct conformation for rcs
			MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple, pmol);
			
			// optimize the degrees of freedom, if needed
			if (mof.getNumDOFs() > 0) {
				try (CudaCCDMinimizer minimizer = new CudaCCDMinimizer(pool, mof)) {
					return minimizer.minimize();
				}
			}

			// otherwise, just evaluate the energy function
			return new Minimizer.Result(null, efunc.getEnergy());
		}
	}
}
