package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectPool;

public class MoleculePool extends ObjectPool<Molecule> {
	
	public MoleculePool(final Molecule mol) {
		super(new Factory<Molecule,Void>() {
			@Override
			public Molecule make(Void context) {
				return new Molecule(mol);
			}
		});
	}
}
