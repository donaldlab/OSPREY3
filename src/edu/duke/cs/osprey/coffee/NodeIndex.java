package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.HashCalculator;

import java.io.File;
import java.util.Arrays;
import java.util.function.Consumer;


public class NodeIndex extends FixedIndex<BigExp,NodeIndex.StateConf> {

	public static class StateConf {

		final int statei;
		final int[] conf;

		public StateConf(int statei, int[] conf) {
			this.statei = statei;
			this.conf = conf;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				Integer.hashCode(statei),
				Arrays.hashCode(conf)
			);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof StateConf && equals((StateConf)other);
		}

		public boolean equals(StateConf other) {
			return this.statei == other.statei
				&& Arrays.equals(this.conf, other.conf);
		}
	}

	public NodeIndex(File file, long maxBytes, MultiStateConfSpace confSpace, Consumer<StateConf> evictionListener) {
		super(
			file,
			maxBytes,
			Serializers.BigExp,
			Serializers.stateConf(confSpace),
			evictionListener
		);
	}
}
