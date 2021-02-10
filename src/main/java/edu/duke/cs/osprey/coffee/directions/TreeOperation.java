package edu.duke.cs.osprey.coffee.directions;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;

import java.io.IOException;


public class TreeOperation extends Operation {

	private int statei;
	private NodeTree tree;

	@SuppressWarnings("unused") // used by hazelcast
	public TreeOperation() {
		statei = -1;
		tree = null;
	}

	public TreeOperation(int statei, NodeTree tree) {
		this.statei = statei;
		this.tree = tree;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return Directions.ServiceName;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		out.writeInt(statei);
		assert (tree.rcs.getPruneMat() == null);
		out.writeInt(tree.rcs.getNumPos());
		for (int posi=0; posi<tree.rcs.getNumPos(); posi++) {
			int[] confs = tree.rcs.get(posi);
			out.writeInt(confs.length);
			for (int conf : confs) {
				out.writeInt(conf);
			}
		}
		if (tree.maxSimultaneousMutations != null) {
			out.writeInt(tree.maxSimultaneousMutations);
		} else {
			out.writeInt(-1);
		}
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		statei = in.readInt();
		var rcsAtPos = new int[in.readInt()][];
		for (int posi=0; posi<rcsAtPos.length; posi++) {
			var rcs = new int[in.readInt()];
			rcsAtPos[posi] = rcs;
			for (int confi=0; confi<rcs.length; confi++) {
				rcs[confi] = in.readInt();
			}
		}
		Integer maxSimultaneousMutations = in.readInt();
		if (maxSimultaneousMutations <= 0) {
			maxSimultaneousMutations = null;
		}
		tree = new NodeTree(new RCs(rcsAtPos), maxSimultaneousMutations);
	}

	@Override
	public final void run() {
		Directions directions = getService();
		directions.receiveTree(statei, tree);
	}
}
