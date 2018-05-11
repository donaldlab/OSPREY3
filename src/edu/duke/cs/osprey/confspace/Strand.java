package edu.duke.cs.osprey.confspace;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;

import java.io.Serializable;

/**
 * A molecule with associated residue flexibility information.
 */
public class Strand implements Serializable {
	
	/** Magic value that represents the wild-type residue type. Used by {@link ResidueFlex#setLibraryRotamers} */
	public static final String WildType = "__WT__";
	
	public static class Builder {
		
		private Molecule mol;
		private String firstResNum;
		private String lastResNum;

		/**
		 * The template library to use for this strand.
		 *
		 * If no template library is specified, a default template library is created.
		 **/
		private ResidueTemplateLibrary templateLib;
		private boolean errorOnNonTemplateResidues;
		
		public Builder(Molecule mol) {
			this.mol = mol;
			this.firstResNum = mol.residues.get(0).getPDBResNumber();
			this.lastResNum = mol.residues.get(mol.residues.size() - 1).getPDBResNumber();
			this.templateLib = new ResidueTemplateLibrary.Builder().build();
			this.errorOnNonTemplateResidues = false;
		}

		@Deprecated
		public Builder setResidues(int firstResNum, int lastResNum) {
			return setResidues(stringResNumForMolec(firstResNum,mol),
					stringResNumForMolec(lastResNum,mol));
		}
		
		public Builder setResidues(String firstResNum, String lastResNum) {
			this.firstResNum = firstResNum;
			this.lastResNum = lastResNum;
			return this;
		}
		
		/**
		 * temporary glue code to support old ResidueTermini, but ResidueTermini will eventually be removed in favor of this Strand class
		 */
		@Deprecated
		public Builder setResidues(ResidueTermini termini) {
			if (termini != null) {
				return setResidues(termini.lBound, termini.uBound);
			}
			return this;
		}
		
		public Builder setTemplateLibrary(ResidueTemplateLibrary val) {
			this.templateLib = val;
			return this;
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
	public static class ResidueFlex implements Serializable {
		
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
		 * Use {@link #WildType} to represent the wild type residue type.
		 * If no residue types are passed, the wild type will be used. 
		 */
		public ResidueFlex setLibraryRotamers(String ... resTypes) {
			return setLibraryRotamers(Arrays.asList(resTypes));
		}
		
		public ResidueFlex setLibraryRotamers(List<String> resTypes) {
			
			this.resTypes.clear();
			
			if (resTypes.isEmpty()) {
				
				// no mutations explicitly chosen, don't assume anything, this is an error
				throw new IllegalArgumentException("no residue types chosen");
				
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

		/** Uses a continuous voxel with default half width of {@link VoxelShape#DefaultHalfWidthDegrees} */
		public ResidueFlex setContinuous() {
			return setVoxelShape(new VoxelShape.Rect());
		}

		public ResidueFlex setContinuous(double voxelHalfWidth) {
			return setVoxelShape(new VoxelShape.Rect(voxelHalfWidth));
		}
		
		public ResidueFlex setVoxelShape(VoxelShape voxelShape) {
			this.voxelShape = voxelShape;
			return this;
		}
	}
	
	public class Flexibility implements Serializable {
		
		private Map<String,ResidueFlex> residues;
		
		public Flexibility(List<Residue> residues) {
			this.residues = new LinkedHashMap<>();
			for (Residue res : residues) {
				this.residues.put(
					Residues.normalizeResNum(res.getPDBResNumber()),
					new ResidueFlex(res.template.name)
				);
			}
		}

		@Deprecated
		public ResidueFlex get(int resNum) {
			return get(stringResNumForMolec(resNum, mol));
		}
		
		public ResidueFlex get(String resNum) {
			resNum = Residues.normalizeResNum(resNum);
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
	
	/** The molecule this strand represents */
	public final Molecule mol;
	
	/** The template library used to pick templates for this strand */
	public final ResidueTemplateLibrary templateLib;
	
	/** Names of residues that couldn't be matched to templates */
	public final Set<String> nonTemplateResNames;
	
	/** Flexibility parameters for this strand */
	public final Flexibility flexibility;
		
	private Strand(Molecule mol, String firstResNumber, String lastResNumber, ResidueTemplateLibrary templateLib, boolean errorOnNonTemplateResidues) {
		
		// make sure the mol has these residues, otherwise the ranges won't work correctly
		mol.residues.getOrThrow(firstResNumber);
		mol.residues.getOrThrow(lastResNumber);
		
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
		nonTemplateResNames = tryAssigningTemplates(this.mol, templateLib);
		
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
		this.mol.markInterResBonds();
		
		// init flexibility
		flexibility = new Flexibility(this.mol.residues);
	}
        
        
    public static LinkedHashSet<String> tryAssigningTemplates(Molecule mol, ResidueTemplateLibrary templateLib){
        LinkedHashSet<String> nonTemplateResNames = new LinkedHashSet<>();
        for (Residue res : mol.residues) {

                // We accept D-amino acid named using the usual L names,
                // but must change them here so the right template name is used
                DAminoAcidHandler.tryRenamingAsD(res);
                for (Residue altRes : mol.getAlternates(res.indexInMolecule)) {
                        DAminoAcidHandler.tryRenamingAsD(altRes);
                }

                // try to assign the template
                boolean templateAssigned = res.assignTemplate(templateLib);
                if (templateAssigned) {

                        // assign the alternates too
                        Iterator<Residue> altIter = mol.getAlternates(res.indexInMolecule).iterator();
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
        return nonTemplateResNames;
    }
	
	private static String stringResNumForMolec(int resNumInt, Molecule molec) {
		//given an integer residue number, find the full (with chain ID) residue number
		//in this molec that matches it.  Error if more than one.  
		String pureNumber = Integer.toString(resNumInt);
		return molec.getResByPDBResNumber(pureNumber).getPDBResNumber();
	}
}
