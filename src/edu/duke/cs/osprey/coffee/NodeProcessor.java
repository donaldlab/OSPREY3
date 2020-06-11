package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.parallelism.TaskExecutor;


public class NodeProcessor {

	public final TaskExecutor tasks;
	public final SeqDB seqdb;
	public final NodeDB nodedb;
	public final StateInfo[] stateInfos;

	public NodeProcessor(TaskExecutor tasks, SeqDB seqdb, NodeDB nodedb, StateInfo[] stateInfos) {
		this.tasks = tasks;
		this.seqdb = seqdb;
		this.nodedb = nodedb;
		this.stateInfos = stateInfos;
	}

	private void log(String msg, Object ... args) {
		nodedb.member.log(msg, args);
	}

	public enum Result {

		ProcessedNode(true),
		NeedInfo(false),
		NoNode(false);

		public final boolean processedNode;

		Result(boolean processedNode) {
			this.processedNode = processedNode;
		}
	}

	public Sequence makeSeq(int statei, int[] conf) {
		return seqdb.confSpace.seqSpace.makeSequence(stateInfos[statei].confSpace, conf);
	}

	public Result process(Directions directions) {

		// get the currently focused state
		int statei = directions.getFocusedStatei();
		if (statei < 0) {
			return Result.NeedInfo;
		}
		var stateInfo = stateInfos[statei];

		// get the tree for this state
		RCs tree = directions.getTree(statei);
		if (tree == null) {
			return Result.NeedInfo;
		}

		// get the next node from that state, or wait for one to happen
		NodeIndex.Node node = nodedb.removeHigh(statei);
		if (node == null) {
			return Result.NoNode;
		}

		var confIndex = stateInfo.makeConfIndex();
		Conf.index(node.conf, confIndex);

		var seqdbBatch = seqdb.batch();

		// is this a leaf node?
		if (confIndex.isFullyDefined()) {

			// compute the full conf weighted energy
			var z = stateInfo.zPath(node.conf);

			seqdbBatch.addZConf(
				stateInfo.config.state,
				makeSeq(statei, node.conf),
				z,
				node.score
			);

		} else {

			// remove the old bound from SeqDB
			seqdbBatch.subZSumUpper(
				stateInfo.config.state,
				makeSeq(statei, node.conf),
				node.score
			);

			// pick the next position to expand
			// TODO: allow permuting the positions
			int posi = confIndex.numDefined;

			// expand the node
			for (int confi : tree.get(posi)) {
				confIndex.assignInPlace(posi, confi);

				// compute an upper bound for the assignment
				var zSumUpper = stateInfo.zSumUpper(confIndex, tree).normalize(true);

				// update nodedb
				var childNode = new NodeIndex.Node(statei, Conf.make(confIndex), zSumUpper);
				nodedb.addLocal(childNode);

				// add the new bound
				seqdbBatch.addZSumUpper(
					stateInfo.config.state,
					makeSeq(statei, childNode.conf),
					childNode.score
				);

				confIndex.unassignInPlace(posi);
			}
		}

		seqdbBatch.save();

		return Result.ProcessedNode;
	}
}
