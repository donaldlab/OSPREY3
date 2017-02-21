package edu.duke.cs.osprey.confspace;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.Forcefield;
import edu.duke.cs.osprey.kstar.KSTermini;
import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

/**
 * A molecule with associated residue flexibility information.
 */
public class Strand {
	
	public static final String WildType = "__WT__";
	
	public static class Builder {
		
		private Molecule mol;
		private String firstResNum;
		private String lastResNum;
		private GenericResidueTemplateLibrary templateLib;
		private boolean errorOnNonTemplateResidues;
		
		public Builder(Molecule mol) {
			this.mol = mol;
			this.firstResNum = mol.residues.get(0).getPDBResNumber();
			this.lastResNum = mol.residues.get(mol.residues.size() - 1).getPDBResNumber();
			this.templateLib = new GenericResidueTemplateLibrary.Builder().build();
			this.errorOnNonTemplateResidues = false;
		}
		
		public Builder setResidues(int firstResNum, int lastResNum) {
			return setResidues(Integer.toString(firstResNum), Integer.toString(lastResNum));
		}
		
		public Builder setResidues(String firstResNum, String lastResNum) {
			this.firstResNum = firstResNum;
			this.lastResNum = lastResNum;
			return this;
		}
		
		/**
		 * temporary glue code to support old KSTermini, but KSTermini will eventually be removed in favor of this Strand class
		 */
		@Deprecated
		public Builder setResidues(KSTermini termini) {
			if (termini != null) {
				return setResidues(termini.getTerminusBegin(), termini.getTerminusEnd());
			}
			return this;
		}
		
		public Builder setTemplateLibrary(GenericResidueTemplateLibrary val) {
			this.templateLib = val;
			return this;
		}
		
		public Builder setDefaultTemplateLibrary(Forcefield forcefield) {
			return setTemplateLibrary(new GenericResidueTemplateLibrary.Builder()
				.setForcefield(forcefield)
				.build());
		}
		
		public Builder setLovellTemplateLibrary() {
			return setTemplateLibrary(new GenericResidueTemplateLibrary.Builder()
				.setLovellRotamers()
				.build());
		}
		
		public Builder setLovellTemplateLibrary(Forcefield forcefield) {
			return setTemplateLibrary(new GenericResidueTemplateLibrary.Builder()
				.setForcefield(forcefield)
				.setLovellRotamers()
				.build());
		}
		
		public Builder setErrorOnNonTemplateResidues(boolean val) {
			this.errorOnNonTemplateResidues = val;
			return this;
		}
		
		public Strand build() {
			return new Strand(mol, firstResNum, lastResNum, templateLib, errorOnNonTemplateResidues);
		}
	}
	
	/**
	 * configured flexibility for one residue
	 */
	public static class ResidueFlex {
		
		public final String wildType;
		
		public VoxelShape voxelShape;
		public Set<String> resTypes;
		public boolean addWildTypeRotamers;
		
		protected ResidueFlex(String wildType) {
			
			this.wildType = wildType;
			
			this.resTypes = new LinkedHashSet<>();
			
			// default flexibility
			setDiscrete();
			setNoRotamers();
		}
		
		public boolean isFlexible() {
			return !resTypes.isEmpty() || addWildTypeRotamers == true;
		}
		
		public ResidueFlex setNoRotamers() {
			resTypes.clear();
			addWildTypeRotamers = false;
			return this;
		}
		
		/**
		 * Add the wild type conformation as a rotamer, including any alternate conformations.
		 */
		public ResidueFlex addWildTypeRotamers() {
			addWildTypeRotamers = true;
			return this;
		}
		
		/**
		 * Set the mutatable residue types (e.g., amino acids) for this residue.
		 * Use {@link WildType} to represent the wild type residue type.
		 * If no residue types are passed, the wild type will be used. 
		 */
		public ResidueFlex setLibraryRotamers(String ... resTypes) {
			return setLibraryRotamers(Arrays.asList(resTypes));
		}
		
		public ResidueFlex setLibraryRotamers(List<String> resTypes) {
			
			this.resTypes.clear();
			
			if (resTypes.isEmpty()) {
				
				// no mutations explicitly chosen, assume wild type only
				this.resTypes.add(wildType);
				
			} else {
				
				for (String resType : resTypes) {
					
					// replace WT with actual amino acid
					if (resType.equalsIgnoreCase(WildType)) {
						resType = wildType;
					}
					
					this.resTypes.add(resType);
				}
			}
			
			return this;
		}
		
		public ResidueFlex setDiscrete() {
			return setVoxelShape(new VoxelShape.Point());
		}
		
		public ResidueFlex setContinuous() {
			return setVoxelShape(new VoxelShape.Rect());
		}
		
		public ResidueFlex setContinuous(double voxelWidth) {
			return setVoxelShape(new VoxelShape.Rect(voxelWidth));
		}
		
		public ResidueFlex setContinuousEllipses() {
			return setVoxelShape(new VoxelShape.Ellipse());
		}
		
		public ResidueFlex setContinuousEllipses(double voxelWidth) {
			return setVoxelShape(new VoxelShape.Ellipse(voxelWidth));
		}
		
		public ResidueFlex setVoxelShape(VoxelShape voxelShape) {
			this.voxelShape = voxelShape;
			return this;
		}
	}
	
	public class Flexibility {
		
		private Map<String,ResidueFlex> residues;
		
		public Flexibility(List<Residue> residues) {
			this.residues = new HashMap<>();
			for (Residue res : residues) {
				this.residues.put(res.getPDBResNumber(), new ResidueFlex(res.template.name));
			}
		}
		
		public ResidueFlex get(int resNum) {
			return get(Integer.toString(resNum));
		}
		
		public ResidueFlex get(String resNum) {
			return residues.get(resNum);
		}
		
		public List<String> getFlexibleResidueNumbers() {
			return residues.entrySet().stream()
				.filter((entry) -> {
					ResidueFlex resFlex = entry.getValue();
					return resFlex.isFlexible();
				})
				.map((entry) -> {
					String resNum = entry.getKey();
					return resNum;
				})
				.collect(Collectors.toList());
		}
		
		public List<String> getStaticResidueNumbers() {
			return residues.entrySet().stream()
				.filter((entry) -> {
					ResidueFlex resFlex = entry.getValue();
					return !resFlex.isFlexible();
				})
				.map((entry) -> {
					String resNum = entry.getKey();
					return resNum;
				})
				.collect(Collectors.toList());
		}
	}
	
	public final Molecule mol;
	public final GenericResidueTemplateLibrary templateLib;
	public final Set<String> nonTemplateResNames;
	public final Flexibility flexibility;
	
	/**
	 * it's probably easier to use the {@link Builder} rather than call this constructor directly.
	 */
	public Strand(Molecule mol, String firstResNumber, String lastResNumber, GenericResidueTemplateLibrary templateLib, boolean errorOnNonTemplateResidues) {
		
		// build our molecule copy from the residue subset
		this.mol = new Molecule();
		for (Residue res : mol.getResRangeByPDBResNumber(firstResNumber, lastResNumber)) {
			
			// copy the main res
			Residue newRes = new Residue(res);
			this.mol.appendResidue(newRes);
			
			// copy the alternates, if any
			for (Residue altRes : mol.getAlternates(res.indexInMolecule)) {
				this.mol.addAlternate(newRes.indexInMolecule, new Residue(altRes));
			}
			
			assert (this.mol.getAlternates(newRes.indexInMolecule).size() == mol.getAlternates(res.indexInMolecule).size());
		}
		
		// assign templates and mark intra-residue bonds
		this.templateLib = templateLib;
		nonTemplateResNames = new LinkedHashSet<>();
		for (Residue res : this.mol.residues) {
			
			// We accept D-amino acid named using the usual L names,
			// but must change them here so the right template name is used
			DAminoAcidHandler.tryRenamingAsD(res);
			for (Residue altRes : this.mol.getAlternates(res.indexInMolecule)) {
				DAminoAcidHandler.tryRenamingAsD(altRes);
			}
			
			// try to assign the template
			boolean templateAssigned = res.assignTemplate(templateLib);
			if (templateAssigned) {
				
				// assign the alternates too
				Iterator<Residue> altIter = this.mol.getAlternates(res.indexInMolecule).iterator();
				while (altIter.hasNext()) {
					Residue altRes = altIter.next();
					boolean altTemplateAssigned = altRes.assignTemplate(templateLib);
					if (!altTemplateAssigned) {
						
						// sometimes alts have fewer atoms than the main residue and we can't do template assignment
						// just delete the alt
						altIter.remove();
					
						nonTemplateResNames.add(altRes.fullName + " alternate");
					}
				}
				
			} else {
				nonTemplateResNames.add(res.fullName);
			}
		}
		
		// delete non template residues if needed
		if (!nonTemplateResNames.isEmpty()) {
			
			// get a list of residues
			String resNames = String.join("\n", nonTemplateResNames);
			
			if (errorOnNonTemplateResidues) {
				throw new Error("ERROR: " + nonTemplateResNames.size() + " Strand residue(s) could not be matched to templates:\n" + resNames);
			}
			
			this.mol.deleteResidues(nonTemplateResNames);
			
			// warn user about deleted residues
			// TODO: write to special log?
			System.out.println("WARNING: " + nonTemplateResNames.size() + " Strand residue(s) could not be matched to templates and were automatically deleted:\n" + resNames);
		}
		
		// assigning templates marks intra-res bonds; we can now mark inter-res too
		HardCodedResidueInfo.markInterResBonds(this.mol);
		
		// init flexibility
		flexibility = new Flexibility(this.mol.residues);
	}
}
