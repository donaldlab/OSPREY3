/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.collections4.queue.CircularFifoQueue;

import edu.duke.cs.osprey.gmec.ConsoleConfPrinter;
import java.util.HashMap;


//This is a general interface for things that search conformational space
//like an A* trees, a BWM* search, or a WCSP solver
//For each of these, we instantiate the ConfSearch, 
//call nextConf to get the GMEC for the model being searched,
//and then successive calls to nextConf() return the next conformation in the gap-free list
/**
 * A generic interface for A* searches over a conformatio space.
 * 
 * Each search defines an order for all conformations in the conformation space,
 * and enumerates those conformations order of increasing score.
 * 
 * @author mhall44
 */
public interface ConfSearch {
    
    /**
     * Get the conformation in the conformation space with the next lowest score.
     */
    ScoredConf nextConf();
    
    /**
     * Get the total number of conformations in the conformation space.
     * This can be an astronomically large number.
     */
    default BigInteger getNumConformations() {
    	throw new UnsupportedOperationException();
    }
    
    /**
     * Get the next conformations in the conformation space with scores up to maxEnergy.
     * @param foo cows are tasty
     * @param bar cheese is too
     */
    default List<ScoredConf> nextConfs(double maxEnergy) {
		List<ScoredConf> nodes = new ArrayList<>();
		while (true) {
			
			ScoredConf conf = nextConf();
			if (conf == null) {
				break;
			}
			
			nodes.add(conf);
			
			if (conf.getScore() >= maxEnergy) {
				break;
			}
		}
		return nodes;
    }
    
    /**
     * A conformation from a conformation space with an associated score.
     */
    public static class ScoredConf {
        
        private int[] assignments;
        private double score;
        
        /** Make a new scored conformation */
        public ScoredConf(int[] assignments, double score) {
            this.assignments = assignments;
            this.score = score;
        }
        
        /** Make a copy of an existing scored conformation */
        public ScoredConf(ScoredConf other) {
        	this.assignments = other.assignments.clone();
        	this.score = other.score;
        }
        
        /**
         * Get the design position assignments of the conformation.
         * 
         * Design positions are listed in order of increasing residue number.
         * Each integer in the array describes the index of the assigned residue
         * conformation at that design position.
         * */
        public int[] getAssignments() {
            return assignments;
        }
        
        /** Get the score for the conformation */
        public double getScore() {
            return score;
        }
        
        /** Sets the score for this conformation */
        public void setScore(double val) {
        	score = val;
        }
        
        /** Adds val to the score for this conformation */
        public void offsetScore(double val) {
        	score += val;
        }

        @Override
		public boolean equals(Object other) {
        	return other instanceof ScoredConf && equals((ScoredConf)other);
		}

		public boolean equals(ScoredConf other) {
        	return Arrays.equals(this.assignments, other.assignments)
				&& Double.compare(this.score, other.score) == 0;
		}
    }
    
    /**
     * A conformation from a conformation space with an associated score,
     * and an associated energy. The definition of "energy" in this context
     * depends on the conformation energy calculator used in the design.
     */
    public static class EnergiedConf extends ScoredConf {
        
        private double energy;
        
        public EnergiedConf(ScoredConf conf, double energy) {
        	super(conf);
            this.energy = energy;
        }
        
        public EnergiedConf(int[] conf, double score, double energy) {
        	super(conf, score);
        	this.energy = energy;
        }
        
        /** Get the energy for this conformation */
        public double getEnergy() {
            return energy;
        }
        
        /** Set the energy for this conformation */
        public void setEnergy(double val) {
        	energy = val;
        }
        
        /** Adds val to the energy for this conformation */
        public void offsetEnergy(double val) {
        	energy += val;
        }
        
        @Override
        public String toString() {
        	return toString(null);
        }
        
        /**
         * Generates a nicely formatted description of this conformation
         * with its score and energy.
         * 
         * For example:
         * <pre>
         * Residue Conf Ids       1  20  27
         * Residue types        GLY ARG LYS
         * Rotamer numbers        L  L9  W0
         * Energy               -21.788082
         * Score                -22.675928 (gap: 0.887846)
         * </pre>
         */
        public String toString(SimpleConfSpace confSpace) {
        	return ConsoleConfPrinter.makeReport(this, confSpace, null);
        }
        
        public HashMap<String,List<String>> getProperties(SimpleConfSpace confSpace){
            return ConsoleConfPrinter.makeReportMap(this, confSpace, null);
        }

		@Override
		public boolean equals(Object other) {
			return other instanceof EnergiedConf && equals((EnergiedConf)other);
		}

		public boolean equals(EnergiedConf other) {
			return super.equals(other)
				&& Double.compare(this.energy, other.energy) == 0;
		}
	}
    
	/**
	 * Lets multiple consumers read confs from the stream regardless of order of reads.
	 */
	public static class Splitter {
		
		public class Stream implements ConfSearch {
			
			private long index;
			
			private Stream() {
				index = 0;
			}
			
			public ConfSearch getSource() {
				return confs;
			}
			
			@Override
			public ScoredConf nextConf() {
			
				// where in the buffer should we read?
				long pos = index - firstIndex;
				
				// just in case
				assert (pos >= 0);
				assert (pos <= buf.size());
				
				ScoredConf conf = null;
				
				// off the end?
				if (pos == buf.size()) {
					
					// read a new conf from the tree
					conf = confs.nextConf();
					
					if (conf != null) {
						
						// if we're out of space, make a bigger queue
						if (buf.isAtFullCapacity()) {
							CircularFifoQueue<ScoredConf> newBuf = new CircularFifoQueue<>(buf.maxSize()*2);
							newBuf.addAll(buf);
							buf = newBuf;
						}
						
						buf.add(conf);
						
					} else {
						// end of the tree
					}
					
				} else {
	
					// read the conf from the buffer
					if (pos > Integer.MAX_VALUE) {
						throw new Error("Integer overflow! Conf buffer grew too large: " + pos);
					}
					conf = buf.get((int)pos);
				}
				
				if (conf != null) {
					index++;
					
					assert (index - firstIndex <= buf.size() + 1);
				}
				
				// prune the first conf from the buffer if everyone's read it
				if (pos == 0 && everyoneHasReadFirst()) {
					buf.remove();
					firstIndex++;
				}
				
				return conf;
			}
				
			private boolean everyoneHasReadFirst() {
				for (Stream other : streams) {
					if (other.index == firstIndex) {
						return false;
					}
				}
				return true;
			}

			@Override
			public BigInteger getNumConformations() {
				return confs.getNumConformations();
			}
			
			public void close() {
				
				// remove our stream from the splitter
				streams.remove(this);
				
				// what's the earliest remaining stream index?
				long minIndex = Long.MAX_VALUE;
				for (Stream stream : streams) {
					minIndex = Math.min(minIndex, stream.index);
				}
				
				// prune the buffer
				while (index < minIndex) {
					buf.remove();
					index++;
					firstIndex++;
				}
			}
		}
		
		private ConfSearch confs;
		private CircularFifoQueue<ScoredConf> buf;
		private long firstIndex;
		private List<Stream> streams;
		
		/**
		 * Create a splitter for a conformation search
		 */
		public Splitter(ConfSearch confs) {
			
			this.confs = confs;
			
			buf = new CircularFifoQueue<>();
			firstIndex = 0;
			streams = new ArrayList<>();
		}
		
		/**
		 * Make a new stream for this conf search that
		 * starts at the current position of the search.
		 */
		public Stream makeStream() {
			
			// don't start a new stream if we've already started reading
			if (firstIndex > 0) {
				throw new IllegalStateException("can't start new stream after first read");
			}
			
			Stream stream = new Stream();
			streams.add(stream);
			return stream;
		}
		
		public int getBufferSize() {
			return buf.size();
		}
	}
}
