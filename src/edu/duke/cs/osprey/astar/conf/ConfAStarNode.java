package edu.duke.cs.osprey.astar.conf;

public interface ConfAStarNode extends Comparable<ConfAStarNode> {

	ConfAStarNode assign(int pos, int rc);
	double getGScore();
	void setGScore(double val);
	double getHScore();
	void setHScore(double val);
	int getLevel();
	void getConf(int[] conf);
	void index(ConfIndex index);
	
	default int[] makeConf(int numPos) {
		int[] conf = new int[numPos];
		getConf(conf);
		return conf;
	}
	
	default double getScore() {
		return getGScore() + getHScore();
	}
	
	@Override
	default int compareTo(ConfAStarNode other) {
		return Double.compare(getScore(), other.getScore());
	}
}
