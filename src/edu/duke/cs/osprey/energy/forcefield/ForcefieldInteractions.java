package edu.duke.cs.osprey.energy.forcefield;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;

public class ForcefieldInteractions extends ArrayList<AtomGroup[]> {
	
	private static final long serialVersionUID = -2052346528544825763L;

	public static class AtomGroup {
		
		private int id;
		private List<Atom> atoms;
		private double[] coords;
		private boolean isDynamic;
		
		public AtomGroup(int id, List<Atom> atoms, double[] coords, boolean isDynamic) {
			this.id = id;
			this.atoms = atoms;
			this.coords = coords;
			this.isDynamic = isDynamic;
		}
		
		public int getId() {
			return id;
		}
		
		public List<Atom> getAtoms() {
			return atoms;
		}
		
		public double[] getCoords() {
			return coords;
		}
		
		public boolean isDynamic() {
			return isDynamic;
		}
	}
	
	public void addResidue(Residue res, boolean isDynamic) {
		AtomGroup group = new AtomGroup(res.indexInMolecule, res.atoms, res.coords, isDynamic);
		add(new AtomGroup[] { group, group });
	}
	
	public void addResiduePair(Residue res1, boolean isDynamic1, Residue res2, boolean isDynamic2) {
		add(new AtomGroup[] {
			new AtomGroup(res1.indexInMolecule, res1.atoms, res1.coords, isDynamic1),
			new AtomGroup(res2.indexInMolecule, res2.atoms, res2.coords, isDynamic2)
		});
	}
}
