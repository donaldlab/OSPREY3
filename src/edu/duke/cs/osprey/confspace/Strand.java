package edu.duke.cs.osprey.confspace;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

import edu.duke.cs.osprey.control.Defaults;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KSTermini;
import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class Strand {
	
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
			this.templateLib = Defaults.genericTemplateLibrary;
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
		
		public Builder setDefaultTemplateLibrary(ForcefieldParams ffparams) {
			return setTemplateLibrary(GenericResidueTemplateLibrary.builder()
				.setForcefield(ffparams)
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
	
	public static Builder builder(Molecule mol) {
		return new Builder(mol);
	}
	
	public final Molecule mol;
	public final GenericResidueTemplateLibrary templateLib;
	public final Set<String> nonTemplateResNames;
	
	/**
	 * it's probably easier to use the Builder rather than call this constructor directly.
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
		
		// DEBUG: make sure various structure the pointers are all set correctly
		for (int i=0; i<this.mol.residues.size(); i++) {
			
			Residue res = this.mol.residues.get(i);
			assert (res.molec == this.mol);
			
			// make sure all the bonds are marked
			assert (res.intraResBondsMarked);
			assert (res.interResBondsMarked);
			
			// make sure all the atoms point to the right residues
			for (Atom atom : res.atoms) {
				assert (atom.res == res);
			}
			
			// check the alternates too
			for (Residue altRes : this.mol.getAlternates(i)) {
				for (Atom atom : altRes.atoms) {
					assert (atom.res == altRes);
				}
			}
		}
	}
}
