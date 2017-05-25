package edu.duke.cs.osprey.energy.forcefield;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

/**
 * Use ResidueInteractions and ResidueForcefieldEnergy instead
 */
@Deprecated
public class ForcefieldInteractions extends ArrayList<AtomGroup[]> {
	
	private static final long serialVersionUID = -2052346528544825763L;
	
	public static interface AtomGroup {
		
		int getId();
		List<Atom> getAtoms();
		double[] getCoords();
		int getSequenceNumber();
		boolean isBroken();
	}
	
	public static class ResidueAtomGroup implements AtomGroup {
		
		private Residue res;
		private ResidueTemplate template;
		private int sequenceNumber;
		
		public ResidueAtomGroup(Residue res) {
			this.res = res;
			this.template = res.template;
			this.sequenceNumber = 0;
		}
		
		public Residue getResidue() {
			return res;
		}
		
		public static int getId(Residue res) {
			return res.indexInMolecule;
		}
		
		@Override
		public int getId() {
			return getId(res);
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
		public int getSequenceNumber() {
			
			// did the template change?
			if (this.template != res.template) {
				this.template = res.template;
				sequenceNumber++;
			}
			
			return sequenceNumber;
		}
		
		@Override
		public boolean isBroken() {
			return !res.confProblems.isEmpty();
		}
	}
	
	private Map<Integer,AtomGroup> groupsById;
	
	public ForcefieldInteractions() {
		groupsById = new HashMap<>();
	}
	
	public ForcefieldInteractions(ResidueInteractions inters, Molecule mol) {
		this();
		
		for (ResidueInteractions.Pair pair : inters) {
			
			if (pair.weight != ResidueInteractions.Pair.IdentityWeight || pair.offset != ResidueInteractions.Pair.IdentityOffset) {
				throw new UnsupportedOperationException("weights and offsets not supported by ForcefieldInteractions");
			}
			
			// make the first atom group
			AtomGroup group1 = makeResidueAtomGroup(mol.getResByPDBResNumber(pair.resNum1));
			groupsById.put(group1.getId(), group1);
			
			// make the second atom group if it's different
			AtomGroup group2;
			if (pair.resNum1.equals(pair.resNum2)) {
				group2 = group1;
			} else {
				group2 = makeResidueAtomGroup(mol.getResByPDBResNumber(pair.resNum2));
				groupsById.put(group2.getId(), group2);
			}
			
			add(new AtomGroup[] { group1, group2 });
		}
	}
	
	private ResidueAtomGroup makeResidueAtomGroup(Residue res) {
		int id = ResidueAtomGroup.getId(res);
		ResidueAtomGroup group = (ResidueAtomGroup)groupsById.get(id);
		if (group == null) {
			group = new ResidueAtomGroup(res);
			groupsById.put(id, group);
		}
		return group;
	}
	
	public void addResidue(Residue res) {
		AtomGroup group = makeResidueAtomGroup(res);
		groupsById.put(group.getId(), group);
		add(new AtomGroup[] { group, group });
	}
	
	public void addResiduePair(Residue res1, Residue res2) {
		AtomGroup group1 = makeResidueAtomGroup(res1);
		AtomGroup group2 = makeResidueAtomGroup(res2);
		groupsById.put(group1.getId(), group1);
		groupsById.put(group2.getId(), group2);
		add(new AtomGroup[] { group1, group2 });
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
	
	public ResidueAtomGroup getResidueAtomGroup(Residue res) {
		AtomGroup group = groupsById.get(res.indexInMolecule);
		if (group instanceof ResidueAtomGroup) {
			return (ResidueAtomGroup)group;
		}
		return null;
	}
}
