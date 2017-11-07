package edu.duke.cs.osprey.energy.forcefield;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ResidueForcefieldEnergyCuda;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.HashCalculator;

public class TestForcefieldEnergy extends TestBase {
	
	private static final double EnergyEpsilon = 1e-12;
	
	private static Strand strand = null;
	
	@BeforeClass
	public static void before() {
		strand = new Strand.Builder(PDBIO.readFile("examples/DAGK/2KDC.P.forOsprey.pdb"))
			.setErrorOnNonTemplateResidues(true)
			.build();
	}
	
	public static class TestResidues {

		public final Molecule mol;
		public final Residue gly06, gly15, ser17, trp18, trp25, arg22, ala24, ile26, phe31, arg32, glu34,
			val36, leu39, trp47, leu48, ile53, arg55, val56, leu57, ile59, val62, leu64, val65, met66;
		
		public TestResidues() {
			mol = new Molecule(strand.mol);
			gly06 = mol.getResByPDBResNumber("6");
			gly15 = mol.getResByPDBResNumber("15");
			ser17 = mol.getResByPDBResNumber("17");
			trp18 = mol.getResByPDBResNumber("18");
			trp25 = mol.getResByPDBResNumber("25");
			arg22 = mol.getResByPDBResNumber("22");
			ala24 = mol.getResByPDBResNumber("24");
			ile26 = mol.getResByPDBResNumber("26");
			phe31 = mol.getResByPDBResNumber("31");
			arg32 = mol.getResByPDBResNumber("32");
			glu34 = mol.getResByPDBResNumber("34");
			val36 = mol.getResByPDBResNumber("36");
			leu39 = mol.getResByPDBResNumber("39");
			trp47 = mol.getResByPDBResNumber("47");
			leu48 = mol.getResByPDBResNumber("48");
			ile53 = mol.getResByPDBResNumber("53");
			arg55 = mol.getResByPDBResNumber("55");
			val56 = mol.getResByPDBResNumber("56");
			leu57 = mol.getResByPDBResNumber("57");
			ile59 = mol.getResByPDBResNumber("59");
			val62 = mol.getResByPDBResNumber("62");
			leu64 = mol.getResByPDBResNumber("64");
			val65 = mol.getResByPDBResNumber("65");
			met66 = mol.getResByPDBResNumber("66");
			
			// clear conf problems
			for (Residue res : mol.residues) {
				res.confProblems.clear();
			}
		}
	}
	
	public static enum IntersType {
		
		AllPairs {
			@Override
			public ResidueInteractions makeInters(Residues residues) {
				ResidueInteractions inters = new ResidueInteractions();
				for (int pos1=0; pos1<residues.size(); pos1++) {
					inters.addSingle(residues.get(pos1).getPDBResNumber());
					for (int pos2=0; pos2<pos1; pos2++) {
						inters.addPair(residues.get(pos1).getPDBResNumber(), residues.get(pos2).getPDBResNumber());
					}
				}
				return inters;
			}
		},
		SingleAndShell {
			@Override
			public ResidueInteractions makeInters(Residues residues) {
				ResidueInteractions inters = new ResidueInteractions();
				inters.addSingle(residues.get(0).getPDBResNumber());
				for (int pos1=1; pos1<residues.size(); pos1++) {
					inters.addPair(residues.get(0).getPDBResNumber(), residues.get(pos1).getPDBResNumber());
				}
				return inters;
			}
		};
		
		public abstract ResidueInteractions makeInters(Residues residues);
	}
	
	public static enum FFType {
		
		NoSolv {
			@Override
			public ForcefieldParams makeFFParams() {
				ForcefieldParams ffparams = new ForcefieldParams();
				ffparams.solvationForcefield = null;
				return ffparams;
			}
		},
		EEF1 {
			@Override
			public ForcefieldParams makeFFParams() {
				ForcefieldParams ffparams = new ForcefieldParams();
				ffparams.solvationForcefield = ForcefieldParams.SolvationForcefield.EEF1;
				return ffparams;
			}
		};
		
		public abstract ForcefieldParams makeFFParams();
	}
	
	
	// stupid glue code...
	private static interface EfuncGen {
		
		EnergyFunction make(Residues residues, ResidueInteractions inters, ForcefieldParams ffparams);
		default void init() {}
		default void cleanup() {}
	
		static class FFInters implements EfuncGen {
			
			public static interface SubGen {
				EnergyFunction make(ForcefieldParams ffparams, ForcefieldInteractions inters);
				default void init() {}
				default void cleanup() {}
			}
			
			static class Cuda implements SubGen {
				
				public static interface SubSubGen {
					EnergyFunction make(GpuStreamPool streams, ForcefieldParams ffparams, ForcefieldInteractions inters);
				}
				
				GpuStreamPool streams;
				SubSubGen subGen;
				
				public Cuda(SubSubGen subGen) {
					this.subGen = subGen;
				}
				
				@Override
				public void init() {
					streams = new GpuStreamPool(1, 1);
				}
				
				@Override
				public EnergyFunction make(ForcefieldParams ffparams, ForcefieldInteractions inters) {
					return subGen.make(streams, ffparams, inters);
				}
				
				@Override
				public void cleanup() {
					streams.cleanup();
				}
			}
			
			static class OpenCL implements SubGen {
				
				public static interface SubSubGen {
					EnergyFunction make(GpuQueuePool queues, ForcefieldParams ffparams, ForcefieldInteractions inters);
				}
				
				GpuQueuePool queues;
				SubSubGen subGen;
				
				public OpenCL(SubSubGen subGen) {
					this.subGen = subGen;
				}
				
				@Override
				public void init() {
					queues = new GpuQueuePool(1, 1);
				}
				
				@Override
				public EnergyFunction make(ForcefieldParams ffparams, ForcefieldInteractions inters) {
					return subGen.make(queues, ffparams, inters);
				}
				
				@Override
				public void cleanup() {
					queues.cleanup();
				}
			}
			
			private SubGen subGen;
			
			public FFInters(SubGen subGen) {
				this.subGen = subGen;
			}
			
			@Override
			public void init() {
				subGen.init();
			}
			
			@Override
			public EnergyFunction make(Residues residues, ResidueInteractions inters, ForcefieldParams ffparams) {
				return subGen.make(ffparams, new ForcefieldInteractions(inters, new Molecule(residues)));
			}
			
			@Override
			public void cleanup() {
				subGen.cleanup();
			}
		}
		
		static class Cuda implements EfuncGen {
			
			public static interface SubGen {
				EnergyFunction make(GpuStreamPool streams, Residues residues, ResidueInteractions inters, ForcefieldParams ffparams);
			}
			
			GpuStreamPool streams;
			SubGen subGen;
			
			public Cuda(SubGen subGen) {
				this.subGen = subGen;
			}
			
			@Override
			public void init() {
				streams = new GpuStreamPool(1, 1);
			}
			
			@Override
			public EnergyFunction make(Residues residues, ResidueInteractions inters, ForcefieldParams ffparams) {
				return subGen.make(streams, residues, inters, ffparams);
			}
			
			@Override
			public void cleanup() {
				streams.cleanup();
			}
		}
	}
	
	private static ResPairCache makeResPairCache(Residues residues, ForcefieldParams ffparams) {
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(residues)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		return new ResPairCache(ffparams, connectivity);
	}
	
	// define all the different forcefields and how to make them
	private static EfuncGen efuncsCpu = (residues, inters, ffparams) -> new EnergyFunctionGenerator(ffparams).residueInteractionEnergy(residues, inters);
	private static EfuncGen efuncsBigCpu = new EfuncGen.FFInters((ffparams, inters) -> new BigForcefieldEnergy(ffparams, inters));
	private static EfuncGen efuncsResidueCpu = (residues, inters, ffparams) -> new ResidueForcefieldEnergy(makeResPairCache(residues, ffparams), inters, residues);
	private static EfuncGen efuncsOpenCL = new EfuncGen.FFInters(new EfuncGen.FFInters.OpenCL((queues, ffparams, inters) -> new GpuForcefieldEnergy(ffparams, inters, queues)));
	private static EfuncGen efuncsCuda = new EfuncGen.FFInters(new EfuncGen.FFInters.Cuda((streams, ffparams, inters) -> new GpuForcefieldEnergy(ffparams, inters, streams)));
	private static EfuncGen efuncsResidueCuda = new EfuncGen.Cuda((streams, residues, inters, ffparams) -> new ResidueForcefieldEnergyCuda(streams, makeResPairCache(residues, ffparams), inters, residues));
	
	
	private void checkEnergies(EfuncGen efuncs, Residues residues, TestParams ... params) {
		efuncs.init();
		try {
			for (TestParams p : params) {
				EnergyFunction efunc = efuncs.make(residues, p.intersType.makeInters(residues), p.ffType.makeFFParams());
				try {
					assertThat(p.toString(), efunc.getEnergy(), isAbsolutely(p.expectedEnergy, EnergyEpsilon));
				} finally {
					EnergyFunction.Tools.cleanIfNeeded(efunc);
				}
			}
		} finally {
			efuncs.cleanup();
		}
	}
	
	
	private static class TestParams {
		
		public final IntersType intersType;
		public final FFType ffType;
		public final double expectedEnergy;
		
		public TestParams(IntersType intersType, FFType ffType, double expectedEnergy) {
			this.intersType = intersType;
			this.ffType = ffType;
			this.expectedEnergy = expectedEnergy;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(intersType.hashCode(), ffType.hashCode());
		}
		
		@Override
		public boolean equals(Object obj) {
			TestParams other = (TestParams)obj;
			return this.intersType == other.intersType
				&& this.ffType == other.ffType;
		}
	}


	public void singleGly(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -4.572136255843063),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -4.572136255843063),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, 0.7914065324002029),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 0.7914065324002029)
		);
	}
	@Test public void singleGlyCpu()         { singleGly(efuncsCpu); }
	@Test public void singleGlyBigCpu()      { singleGly(efuncsBigCpu); }
	@Test public void singleGlyResidueCpu()  { singleGly(efuncsResidueCpu); }
	@Test public void singleGlyOpenCL()      { singleGly(efuncsOpenCL); }
	@Test public void singleGlyCuda()        { singleGly(efuncsCuda); }
	@Test public void singleGlyResidueCuda() { singleGly(efuncsResidueCuda); }
	
	
	public void glyPair(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly06, r.gly15);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -9.17380398335906),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -4.601667727515996),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, 1.5576201272528598),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 0.7662135948526569)
		);
	}
	@Test public void glyPairCpu()         { glyPair(efuncsCpu); }
	@Test public void glyPairBigCpu()      { glyPair(efuncsBigCpu); }
	@Test public void glyPairResidueCpu()  { glyPair(efuncsResidueCpu); }
	@Test public void glyPairOpenCL()      { glyPair(efuncsOpenCL); }
	@Test public void glyPairCuda()        { glyPair(efuncsCuda); }
	@Test public void glyPairResidueCuda() { glyPair(efuncsResidueCuda); }
	
	
	public void glySerPair(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -9.48559560659799),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -2.6911081922156552),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, 3.0264132059679643),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 1.899770277423194)
		);
	}
	@Test public void glySerPairCpu()         { glySerPair(efuncsCpu); }
	@Test public void glySerPairBigCpu()      { glySerPair(efuncsBigCpu); }
	@Test public void glySerPairResidueCpu()  { glySerPair(efuncsResidueCpu); }
	@Test public void glySerPairOpenCL()      { glySerPair(efuncsOpenCL); }
	@Test public void glySerPairCuda()        { glySerPair(efuncsCuda); }
	@Test public void glySerPairResidueCuda() { glySerPair(efuncsResidueCuda); }
	
	
	public void trpPair(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.trp18, r.trp25);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -12.625574526252965),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -6.218018599252964),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, 3.6482861076947595),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 1.8445932583378253)
		);
	}
	@Test public void trpPairCpu()         { trpPair(efuncsCpu); }
	@Test public void trpPairBigCpu()      { trpPair(efuncsBigCpu); }
	@Test public void trpPairResidueCpu()  { trpPair(efuncsResidueCpu); }
	@Test public void trpPairOpenCL()      { trpPair(efuncsOpenCL); }
	@Test public void trpPairCuda()        { trpPair(efuncsCuda); }
	@Test public void trpPairResidueCuda() { trpPair(efuncsResidueCuda); }
	
	
	public void the4Residues(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -23.31199205572296),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -2.756905624257449),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, 2.9701822745081854),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 1.768239618380071)
		);
	}
	@Test public void the4ResiduesCpu()         { the4Residues(efuncsCpu); }
	@Test public void the4ResiduesBigCpu()      { the4Residues(efuncsBigCpu); }
	@Test public void the4ResiduesResidueCpu()  { the4Residues(efuncsResidueCpu); }
	@Test public void the4ResiduesOpenCL()      { the4Residues(efuncsOpenCL); }
	@Test public void the4ResiduesCuda()        { the4Residues(efuncsCuda); }
	@Test public void the4ResiduesResidueCuda() { the4Residues(efuncsResidueCuda); }
	
	
	public void the6Residues(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -52.31617653073317),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -2.7906943839799343),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, -9.856849417475418),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 1.7344508586575857)
		);
	}
	@Test public void the6ResiduesCpu()         { the6Residues(efuncsCpu); }
	@Test public void the6ResiduesBigCpu()      { the6Residues(efuncsBigCpu); }
	@Test public void the6ResiduesResidueCpu()  { the6Residues(efuncsResidueCpu); }
	@Test public void the6ResiduesOpenCL()      { the6Residues(efuncsOpenCL); }
	@Test public void the6ResiduesCuda()        { the6Residues(efuncsCuda); }
	@Test public void the6ResiduesResidueCuda() { the6Residues(efuncsResidueCuda); }
	
	
	public void the10Residues(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -93.3333779512777),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -2.7991581273906516),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, -17.302367217342816),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 1.7259871152468682)
		);
	}
	@Test public void the10ResiduesCpu()         { the10Residues(efuncsCpu); }
	@Test public void the10ResiduesBigCpu()      { the10Residues(efuncsBigCpu); }
	@Test public void the10ResiduesResidueCpu()  { the10Residues(efuncsResidueCpu); }
	@Test public void the10ResiduesOpenCL()      { the10Residues(efuncsOpenCL); }
	@Test public void the10ResiduesCuda()        { the10Residues(efuncsCuda); }
	@Test public void the10ResiduesResidueCuda() { the10Residues(efuncsResidueCuda); }
	
	
	public void the14Residues(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -112.44246575304817),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -2.799964297346741),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, -21.693144259934765),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 1.7251809452907791)
		);
	}
	@Test public void the14ResiduesCpu()         { the14Residues(efuncsCpu); }
	@Test public void the14ResiduesBigCpu()      { the14Residues(efuncsBigCpu); }
	@Test public void the14ResiduesResidueCpu()  { the14Residues(efuncsResidueCpu); }
	@Test public void the14ResiduesOpenCL()      { the14Residues(efuncsOpenCL); }
	@Test public void the14ResiduesCuda()        { the14Residues(efuncsCuda); }
	@Test public void the14ResiduesResidueCuda() { the14Residues(efuncsResidueCuda); }
	
	
	public void the24Residues(EfuncGen efuncs) {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		);
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   -163.74206898485193),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   -4.612537058185951),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, -34.03236066473183),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, 0.7553442641827014)
		);
	}
	@Test public void the24ResiduesCpu()         { the24Residues(efuncsCpu); }
	@Test public void the24ResiduesBigCpu()      { the24Residues(efuncsBigCpu); }
	@Test public void the24ResiduesResidueCpu()  { the24Residues(efuncsResidueCpu); }
	@Test public void the24ResiduesOpenCL()      { the24Residues(efuncsOpenCL); }
	@Test public void the24ResiduesCuda()        { the24Residues(efuncsCuda); }
	@Test public void the24ResiduesResidueCuda() { the24Residues(efuncsResidueCuda); }
	
	
	public void brokenProline(EfuncGen efuncs) {
		
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15);
		
		// mutate to a proline, which will be broken at this pos
		Residue res = r.gly15;
		res.pucker = new ProlinePucker(strand.templateLib, res);
		ResidueTypeDOF.switchToTemplate(strand.templateLib, res, strand.templateLib.getTemplate("PRO"), true);
		
		assertThat(res.confProblems.size(), is(1));
		checkEnergies(efuncs, residues,
			new TestParams(IntersType.AllPairs,       FFType.EEF1,   Double.POSITIVE_INFINITY),
			new TestParams(IntersType.SingleAndShell, FFType.EEF1,   Double.POSITIVE_INFINITY),
			new TestParams(IntersType.AllPairs,       FFType.NoSolv, Double.POSITIVE_INFINITY),
			new TestParams(IntersType.SingleAndShell, FFType.NoSolv, Double.POSITIVE_INFINITY)
		);
	}
	@Test public void brokenProlineCpu()         { brokenProline(efuncsCpu); }
	@Test public void brokenProlineBigCpu()      { brokenProline(efuncsBigCpu); }
	@Test public void brokenProlineResidueCpu()  { brokenProline(efuncsResidueCpu); }
	@Test public void brokenProlineOpenCL()      { brokenProline(efuncsOpenCL); }
	@Test public void brokenProlineCuda()        { brokenProline(efuncsCuda); }
	@Test public void brokenProlineResidueCuda() { brokenProline(efuncsResidueCuda); }
	
	
	private double calcEnergy(EfuncGen efuncs, Residues residues, ResidueInteractions inters) {
		EnergyFunction efunc = efuncs.make(residues, inters, FFType.EEF1.makeFFParams());
		try {
			return efunc.getEnergy();
		} finally {
			EnergyFunction.Tools.cleanIfNeeded(efunc);
		}
	}
	
	public void oneIntraWeight(EfuncGen efuncs) {
		
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15);
		String resNum = residues.get(0).getPDBResNumber();
		
		efuncs.init();
		
		// check base value
		final double baseEnergy = -4.5721362558430645;
		ResidueInteractions inters = new ResidueInteractions();
		inters.addSingle(resNum);
		assertThat("base", calcEnergy(efuncs, residues, inters), isAbsolutely(baseEnergy, EnergyEpsilon));
		
		// test weight
		for (double weight : Arrays.asList(-0.5, -2.0, -1.0, 0.0, 0.5, 2.0)) {
			inters = new ResidueInteractions();
			inters.addSingle(resNum, weight, 0);
			assertThat("weight: " + weight, calcEnergy(efuncs, residues, inters), isAbsolutely(baseEnergy*weight, EnergyEpsilon));
		}
		
		efuncs.cleanup();
		
	}
	@Test public void oneIntraWeightCpu()         { oneIntraWeight(efuncsCpu); }
	@Test public void oneIntraWeightResidueCpu()  { oneIntraWeight(efuncsResidueCpu); }
	@Test public void oneIntraWeightResidueCuda() { oneIntraWeight(efuncsResidueCuda); }
	
	
	public void oneIntraOffset(EfuncGen efuncs) {
		
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15);
		String resNum = residues.get(0).getPDBResNumber();
		
		efuncs.init();
		
		// check base value
		final double baseEnergy = -4.5721362558430645;
		ResidueInteractions inters = new ResidueInteractions();
		inters.addSingle(resNum);
		assertThat("base", calcEnergy(efuncs, residues, inters), isAbsolutely(baseEnergy, EnergyEpsilon));
		
		// test offset
		for (double offset : Arrays.asList(-0.5, -2.0, -1.0, 0.0, 0.5, 2.0)) {
			inters = new ResidueInteractions();
			inters.addSingle(resNum, 1, offset);
			assertThat("offset: " + offset, calcEnergy(efuncs, residues, inters), isAbsolutely(baseEnergy + offset, EnergyEpsilon));
		}
		
		efuncs.cleanup();
	}
	@Test public void oneIntraOffsetCpu()         { oneIntraOffset(efuncsCpu); }
	@Test public void oneIntraOffsetResidueCpu()  { oneIntraOffset(efuncsResidueCpu); }
	@Test public void oneIntraOffsetResidueCuda() { oneIntraOffset(efuncsResidueCuda); }
}
