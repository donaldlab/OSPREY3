package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectPool;

public class ParameterizedMoleculePool extends ObjectPool<ParameterizedMoleculeCopy> {
	
	public ParameterizedMoleculePool(final ConfSpace confSpace) {
		super(new Factory<ParameterizedMoleculeCopy,Void>() {
			@Override
			public ParameterizedMoleculeCopy make(Void context) {
				return new ParameterizedMoleculeCopy(confSpace);
			}
		});
	}
}
