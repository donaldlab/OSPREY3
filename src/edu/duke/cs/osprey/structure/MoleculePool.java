package edu.duke.cs.osprey.structure;

import java.util.ArrayDeque;
import java.util.Deque;

public class MoleculePool {
	
	private Molecule mol;
	private Deque<Molecule> mols;
	
	public MoleculePool(Molecule mol) {
		this.mol = mol;
		this.mols = new ArrayDeque<>();
	}
	
	public synchronized Molecule checkout() {
		if (mols.isEmpty()) {
			mols.add(new Molecule(mol));
		}
		return mols.removeFirst();
	}
	
	public synchronized void release(Molecule mol) {
		mols.add(mol);
	}
}
