package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.tools.MathTools;
import org.apache.commons.lang3.builder.ToStringBuilder;

import java.util.List;

/**
 * Result class used for benchmarking, passing stats through to python interface
 */
public class Result {
    public double runtimeS; //runtime in s
    public long startNs; // stats start ns
    public long stopNs; // stats end ns
    public long minimized; // number of nodes expanded
    public long expanded; // number of nodes expanded
    public long rescored; // number of nodes rescored
    public long finished; // number of nodes finished
    public List<SeqFreeEnergies> seqs; // list of finished sequences
    public MathTools.DoubleBounds targetFreeEnergy; // target free energy

    public Result(double runtimeS){
        this.runtimeS=runtimeS;
    }

    public Result(double runtimeS, NodeStats.Report stats, List<SeqFreeEnergies> seqs, MathTools.DoubleBounds targetFreeEnergy){
        this.runtimeS=runtimeS;
        this.startNs=stats.startNs;
        this.stopNs=stats.stopNs;
        this.minimized=stats.values.minimized;
        this.expanded=stats.values.expanded;
        this.rescored=stats.values.rescored;
        this.finished=stats.values.finished;
        this.seqs=seqs;
        this.targetFreeEnergy=targetFreeEnergy;
    }

    public String toString(){
        return ToStringBuilder.reflectionToString(this);
    }
}
