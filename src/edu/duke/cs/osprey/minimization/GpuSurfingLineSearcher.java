package edu.duke.cs.osprey.minimization;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.List;

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.GpuQueue;
import edu.duke.cs.osprey.gpu.kernels.DihedralMinimizer;
import edu.duke.cs.osprey.gpu.kernels.ForceFieldKernel;
import edu.duke.cs.osprey.structure.Residue;

public class GpuSurfingLineSearcher implements LineSearcher.NeedsCleanup {
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs
	
	private ObjectiveFunction.OneDof f;
	private GpuForcefieldEnergy efunc;
	private int coordsOffset;
	
	private DihedralMinimizer.PoseKernel poseKernel;
	private DihedralMinimizer.SearchKernel[] searchKernels;
	private DihedralMinimizer.SurfKernel surfKernel;
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	// DEBUG
	//private GpuStyleSurfingLineSearcher cpuSide;
	
	@Override
	public void init(ObjectiveFunction.OneDof f) {
		
		this.f = f;
		
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
		
		// HACKHACK: break the interfaces to get the pieces we need
		// the runtime will complain if we can't do these casts, so we're still safe, but the compiler has no idea what's going on
		MoleculeModifierAndScorer mof = (MoleculeModifierAndScorer)f.getParent();
		FreeDihedral dof = (FreeDihedral)mof.getDOFs().get(f.getDimension());
		efunc = ((GpuForcefieldEnergy)mof.getEfunc(f.getDimension()));
		
		ForceFieldKernel ffKernel = efunc.getKernel();
		
		// translate residue atom indices to global atom indices
		Residue res = dof.getResidue();
		ForcefieldInteractions.ResidueAtomGroup group = ffKernel.getForcefield().getInteractions().getResidueAtomGroup(res);
		coordsOffset = ffKernel.getForcefield().getAtomOffset(group);
		
		int[] dihedralAtomIndicesSrc = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
		int[] dihedralAtomIndices = new int[dihedralAtomIndicesSrc.length];
		for (int i=0; i<dihedralAtomIndicesSrc.length; i++) {
			dihedralAtomIndices[i] = dihedralAtomIndicesSrc[i] + coordsOffset;
		}
		
		List<Integer> rotatedAtomIndicesSrc = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
		int[] rotatedAtomIndices = new int[rotatedAtomIndicesSrc.size()];
		for (int i=0; i<rotatedAtomIndicesSrc.size(); i++) {
			rotatedAtomIndices[i] = rotatedAtomIndicesSrc.get(i) + coordsOffset;
		}
		
		// init kernels
		try {
			
			GpuQueue queue = ffKernel.getQueue();
			
			poseKernel = new DihedralMinimizer.PoseKernel(queue);
			searchKernels = new DihedralMinimizer.SearchKernel[9];
			for (int i=0; i<7; i++) {
				searchKernels[i] = new DihedralMinimizer.SearchKernel(queue, i);
			}
			surfKernel = new DihedralMinimizer.SurfKernel(queue);
			
		} catch (IOException ex) {
			throw new Error("can't init dihedral minimizer kernel", ex);
		}
		
		int numEnergies = ffKernel.getEnergySize(efunc.getSubset());
		poseKernel.initForcefield(ffKernel.getCoords(), dihedralAtomIndices, rotatedAtomIndices, ffKernel.getArgs());
		for (int i=0; i<7; i++) {
			searchKernels[i].init(ffKernel.getEnergies(), numEnergies, ffKernel.getArgs(), poseKernel.getArgs());
		}
		surfKernel.init(ffKernel.getEnergies(), numEnergies, ffKernel.getArgs(), poseKernel.getArgs());
		
		/* DEBUG: check against the cpu implementation
		// make sure we use the cpu efunc and not the gpu efunc though
		cpuSide = new GpuStyleSurfingLineSearcher();
		cpuSide.init(new ObjectiveFunction.OneDof(
			new MoleculeModifierAndScorer(ffKernel.getForcefield(), mof.getConstraints(), mof.getMolec(), mof.DOFs),
			f.getDimension()
		));
		*/
	}
	
	@Override
	public double search(double xd) {
		
		/* PROFILING
		ProfilingEvents events = new ProfilingEvents(100);
		efunc.getKernel().setProfilingEvents(events);
		poseKernel.setProfilingEvents(events);
		for (int i=0; i<7; i++) {
			searchKernels[i].setProfilingEvents(events);
		}
		surfKernel.setProfilingEvents(events);
		*/
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		double step;
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// PROFILING
		//Profiler p = new Profiler("args");
		
		ForceFieldKernel ffKernel = efunc.getKernel();
		GpuQueue queue = efunc.getKernel().getQueue();
		
		// init pose/search/surf kernel args
		poseKernel.initArgs(step, f.getXMin(), f.getXMax(), xd, efunc.getSubset().getInternalSolvationEnergy());
		
		// PROFILING
		//queue.waitForGpu();
		//p.start("ff");
		
		// prep the forcefield kernel
		boolean subsetChanged = ffKernel.setSubset(efunc.getSubset());
		if (!subsetChanged) {
			ffKernel.setDoEnergy(true);
		}
		
		/* DEBUG
		cpuSide.preSearch(xd);
		check("pre");
		*/
		
		// PROFILING
		//queue.waitForGpu();
		//p.start("kernels");
		
		// pipeline all the kernels
		for (int i=0; i<=3; i++) {
			
			poseKernel.runAsync();
			ffKernel.runAsync();
			searchKernels[i].runAsync();
			
			// make sure all the kernels get sent to the gpu right away
			queue.flush();
			
			/* DEBUG
			cpuSide.pose();
			cpuSide.energy();
			cpuSide.search(i);
			check("search " + i);
			*/
		}
		
		for (int i=0; i<100; i++) {
			
			surfKernel.runAsync();
			
			/* DEBUG
			cpuSide.surf();
			check("surf " + i + " a");
			*/
			
			// should we stop surfing?
			// NOTE: most line searches stop after the first surf step
			// so it's worth it to pay a download cost now to see if we should
			// otherwise, only download and check every other surf
			if (i % 2 == 0 && !ffKernel.isDoEnergy()) {
				break;
			}
			
			poseKernel.runAsync();
			ffKernel.runAsync();
			
			/* DEBUG
			cpuSide.pose();
			cpuSide.energy();
			check("surf " + i + " b");
			*/
		}
		
		searchKernels[4].runAsync();
		
		/* DEBUG
		cpuSide.search(4);
		check("search 4");
		*/
		
		for (int i=5; i<=6; i++) {
			
			poseKernel.runAsync();
			ffKernel.runAsync();
			searchKernels[i].runAsync();
			
			// make sure all the kernels get sent to the gpu right away
			queue.flush();
			
			/* DEBUG
			cpuSide.pose();
			cpuSide.energy();
			cpuSide.search(i);
			check("search " + i);
			*/
		}
		
		poseKernel.runAsync();
		
		/* DEBUG
		cpuSide.pose();
		check("final pose");
		*/
		
		// PROFILING
		//queue.waitForGpu();
		//p.start("download");
		
		// download the final outputs from the gpu
		ByteBuffer args = poseKernel.downloadArgsSync();
		double bestDihedral = args.getDouble(80);
		double lastStep = args.getDouble(40);
		
		// PROFILING
		//p.start("post");
		
		// update the step
		this.lastStep = lastStep;
		if (iteration == 0) {
			firstStep = lastStep;
		}
		
		/* DEBUG
		cpuSide.postSearch();
		check("post");
		*/
		
		// pose the cpu-side protein too
		f.setX(bestDihedral);
		
		/* PROFILING
		p.stop();
		System.out.println(String.format("i:%d, d:%2d, gpu kernels:%s, %s",
			iteration, f.getDimension(), TimeFormatter.format(events.getTotalNs(), TimeUnit.MICROSECONDS), p.makeReport(TimeUnit.MICROSECONDS)
		));
		events.cleanup();
		efunc.getKernel().setProfilingEvents(null);
		poseKernel.setProfilingEvents(null);
		for (int i=0; i<7; i++) {
			searchKernels[i].setProfilingEvents(null);
		}
		surfKernel.setProfilingEvents(null);
		*/
		
		iteration++;
		
		return bestDihedral;
	}
	
	/* DEBUG
	private void check(String name) {
		cpuSide.check(
			name,
			efunc.getKernel().downloadCoordsSync(), coordsOffset,
			poseKernel.downloadArgsSync(),
			efunc.getSubset().getInternalSolvationEnergy(), firstStep, lastStep, iteration
		);
	}
	*/

	@Override
	public void cleanup() {
		if (poseKernel != null) {
			poseKernel.cleanup();
		}
		for (int i=0; i<9; i++) {
			if (searchKernels[i] != null) {
				searchKernels[i].cleanup();
			}
		}
		if (surfKernel != null) {
			surfKernel.cleanup();
		}
	}
}
