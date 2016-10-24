package edu.duke.cs.osprey.gpu.cuda;

import java.io.IOException;
import java.nio.DoubleBuffer;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelCuda;
import edu.duke.cs.osprey.gpu.cuda.kernels.TestKernel;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class CudaPlayground extends TestBase {
	
	private static final int NumElements = 20;
	private static final int NumRuns = 1000;
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// NOTE: samples and such here:
		// https://github.com/jcuda/jcuda-samples/tree/master/JCudaSamples/src/main/java/jcuda
		
		// info on dynamic parallelism:
		// http://docs.nvidia.com/cuda/cuda-c-programming-guide/#cuda-dynamic-parallelism
		
		// init CUDA
		Context context = new Context(Gpus.get().getGpus().get(0));
		
		forcefield(context);
		//hostLoop(context);
		//deviceLoop(context);
	
		context.cleanup();
	}
	
	private static void hostLoop(Context context)
	throws IOException {
		
		TestKernel kernel = new TestKernel(context, NumElements);
		
		// init host buffers
		DoubleBuffer a = kernel.getA();
		DoubleBuffer b = kernel.getB();
		a.clear();
		b.clear();
		for (int i=0; i<NumElements; i++) {
			a.put((double)i);
			b.put((double)i);
		}
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		for (int i=0; i<NumRuns; i++) {
			kernel.uploadAsync();
			kernel.runAsync();
			kernel.downloadSync();
		}
		stopwatch.stop();
		System.out.println("total time: " + TimeFormatter.format(stopwatch.getTimeNs(), 1));
		System.out.println("avg time per op: " + TimeFormatter.format(stopwatch.getTimeNs()/NumRuns, 1));
		
		System.out.println("sum");
		DoubleBuffer out = kernel.getOut();
		out.rewind();
		for (int i=0; i<10; i++) {
			System.out.println(out.get());
		}
		
		kernel.cleanup();
	}
	
	private static void deviceLoop(Context context)
	throws IOException {
		
		// get the kernel func
		CUmodule module = context.getKernel("testDP");
		CUfunction func = new CUfunction();
		JCudaDriver.cuModuleGetFunction(func, module, "loop");
		
		// init host buffers
		float[] a = new float[NumElements];
		float[] b = new float[NumElements];
		for (int i=0; i<NumElements; i++) {
			a[i] = (float)i;
			b[i] = (float)i;
		}
		float[] out = new float[NumElements];
		
		// allocate buffers
		CUdeviceptr pdA = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdA, NumElements*Sizeof.FLOAT);
		CUdeviceptr pdB = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdB, NumElements*Sizeof.FLOAT);
		CUdeviceptr pdOut = new CUdeviceptr();
		JCudaDriver.cuMemAlloc(pdOut, NumElements*Sizeof.FLOAT);
		
		// prep kernel args
		Pointer pParams = Pointer.to(
			Pointer.to(new int[] {NumElements}),
			Pointer.to(pdA),
			Pointer.to(pdB),
			Pointer.to(pdOut)
		);
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		// upload inputs
		JCudaDriver.cuMemcpyHtoDAsync(pdA, Pointer.to(a), NumElements*Sizeof.FLOAT, null);
		JCudaDriver.cuMemcpyHtoDAsync(pdB, Pointer.to(b), NumElements*Sizeof.FLOAT, null);
		
		// launch the kernel
		JCudaDriver.cuLaunchKernel(
			func,
			1, 1, 1,
			1, 1, 1,
			0,
			null,
			pParams,
			null
		);
		
		// download the output
		JCudaDriver.cuMemcpyDtoH(Pointer.to(out), pdOut, NumElements*Sizeof.FLOAT);
		
		stopwatch.stop();
		System.out.println("total time: " + TimeFormatter.format(stopwatch.getTimeNs(), 1));
		System.out.println("avg time per op: " + TimeFormatter.format(stopwatch.getTimeNs()/NumRuns, 1));
		
		System.out.println("sum");
		for (int i=0; i<10; i++) {
			System.out.println(out[i]);
		}
		
		// cleanup
		JCudaDriver.cuMemFree(pdA);
		JCudaDriver.cuMemFree(pdB);
		JCudaDriver.cuMemFree(pdOut);
	}
	
	private static void forcefield(Context context)
	throws IOException {
		
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		
		// read a protein and get some arbitrary residues
		Molecule mol = PDBFileReader.readPDBFile("test/DAGK/2KDC.P.forOsprey.pdb");
		Residue gly06 = mol.getResByPDBResNumber("6");
		Residue gly15 = mol.getResByPDBResNumber("15");
		Residue ser17 = mol.getResByPDBResNumber("17");
		Residue trp18 = mol.getResByPDBResNumber("18");
		Residue trp25 = mol.getResByPDBResNumber("25");
		Residue arg22 = mol.getResByPDBResNumber("22");
		Residue ala24 = mol.getResByPDBResNumber("24");
		Residue ile26 = mol.getResByPDBResNumber("26");
		Residue phe31 = mol.getResByPDBResNumber("31");
		Residue arg32 = mol.getResByPDBResNumber("32");
		Residue glu34 = mol.getResByPDBResNumber("34");
		Residue val36 = mol.getResByPDBResNumber("36");
		Residue leu39 = mol.getResByPDBResNumber("39");
		Residue trp47 = mol.getResByPDBResNumber("47");
		Residue leu48 = mol.getResByPDBResNumber("48");
		Residue ile53 = mol.getResByPDBResNumber("53");
		Residue arg55 = mol.getResByPDBResNumber("55");
		Residue val56 = mol.getResByPDBResNumber("56");
		Residue leu57 = mol.getResByPDBResNumber("57");
		Residue ile59 = mol.getResByPDBResNumber("59");
		Residue val62 = mol.getResByPDBResNumber("62");
		Residue leu64 = mol.getResByPDBResNumber("64");
		Residue val65 = mol.getResByPDBResNumber("65");
		Residue met66 = mol.getResByPDBResNumber("66");
		
		//Residue[] residues = { gly06 };
		//Residue[] residues = { gly06, gly15 };
		//Residue[] residues = { trp18 };
		//Residue[] residues = { trp18, trp25 };
		Residue[] residues = {
			gly06, gly15, ser17, trp18, trp25, arg22, ala24, ile26, phe31, arg32, glu34, val36,
			leu39, trp47, leu48, ile53, arg55, val56, leu57, ile59, val62, leu64, val65, met66
		};
		
		// make the all pairs energy functions
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		for (int pos1=0; pos1<residues.length; pos1++) {
			interactions.addResidue(residues[pos1]);
			for (int pos2=0; pos2<pos1; pos2++) {
				interactions.addResiduePair(residues[pos1], residues[pos2]);
			}
		}
		
		// test the cpu energy
		BigForcefieldEnergy bigFfenergy = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		System.out.println("cpu energy: " + bigFfenergy.getEnergy());
		
		// make the kernel
		ForcefieldKernelCuda kernel = new ForcefieldKernelCuda(context);
		kernel.setForcefield(bigFfenergy);
		
		// TODO
		kernel.uploadCoordsAsync();
		kernel.runAsync();
		System.out.println("gpu energy: " + kernel.downloadEnergySync());
		
		kernel.cleanup();
	}
}
