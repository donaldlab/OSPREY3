package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.parallelism.TaskExecutor;


public class NodeProcessor {

	public final TaskExecutor tasks;
	public final SeqDB seqdb;
	public final NodeDB nodedb;

	public NodeProcessor(TaskExecutor tasks, SeqDB seqdb, NodeDB nodedb) {
		this.tasks = tasks;
		this.seqdb = seqdb;
		this.nodedb = nodedb;
	}

	private void log(String msg, Object ... args) {
		nodedb.member.log(msg, args);
	}

	public void process(NodeIndex.Node node) {

		// TEMP
		log("process node: %s", node);

		// TODO: expand the node, compute bounds, update seqdb,nodedb
	}
}
