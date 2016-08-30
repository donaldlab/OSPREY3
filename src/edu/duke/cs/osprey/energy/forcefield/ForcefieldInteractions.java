package edu.duke.cs.osprey.energy.forcefield;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;

public class ForcefieldInteractions extends ArrayList<AtomGroup[]> {
	
	private static final long serialVersionUID = -2052346528544825763L;
	
	public static interface AtomGroup {
		
		int getId();
		List<Atom> getAtoms();
		double[] getCoords();
		boolean hasChemicalChange();
		void ackChemicalChange();
	}
	
	public static class ResidueAtomGroup implements AtomGroup {
		
		private Residue res;
		private ResidueTemplate template;
		
		public ResidueAtomGroup(Residue res) {
			this.res = res;
		}
		
		@Override
		public int getId() {
			return res.indexInMolecule;
		}
		
		@Override
		public List<Atom> getAtoms() {
			return res.atoms;
		}
		
		@Override
		public double[] getCoords() {
			return res.coords;
		}
		
		@Override
		public boolean hasChemicalChange() {
			return template != res.template;
		}
		
		@Override
		public void ackChemicalChange() {
			template = res.template;
		}
	}
	
	public void addResidue(Residue res) {
		AtomGroup group = new ResidueAtomGroup(res);
		add(new AtomGroup[] { group, group });
	}
	
	public void addResiduePair(Residue res1, Residue res2) {
		add(new AtomGroup[] { new ResidueAtomGroup(res1), new ResidueAtomGroup(res2) });
	}
	
	public ForcefieldInteractions makeSubsetByResidue(Residue res) {
		ForcefieldInteractions subInteractions = new ForcefieldInteractions();
		for (AtomGroup[] groupPair : this) {
			if (groupPair[0].getId() == res.indexInMolecule || groupPair[1].getId() == res.indexInMolecule) {
				subInteractions.add(groupPair);
			}
		}
		return subInteractions;
	}
}
