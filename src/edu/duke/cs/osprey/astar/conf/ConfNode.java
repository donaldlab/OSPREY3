package edu.duke.cs.osprey.astar.conf;

public interface ConfNode {
    ConfNode assign(int pos, int rc);

    void getConf(int[] conf);
    void index(ConfIndex index);
    int getLevel();

    default int[] makeConf(int numPos) {
        int[] conf = new int[numPos];
        getConf(conf);
        return conf;
    }
}
