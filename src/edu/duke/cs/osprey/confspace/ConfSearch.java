/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.confspace;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.externalMemory.ScoredConfFIFOSerializer;
import org.apache.commons.collections4.queue.CircularFifoQueue;

import edu.duke.cs.osprey.gmec.ConsoleConfPrinter;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicBoolean;


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
	 * Get the next `num` conformations in the conformation space
	 */
    default List<ScoredConf> nextConfs(int num) {
		List<ScoredConf> confs = new ArrayList<>();
		for (int i=0; i<num; i++) {

			ScoredConf conf = nextConf();
			if (conf == null) {
				break;
			}

			confs.add(conf);
		}
		return confs;
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

		@Override
		public String toString() {
        	return String.format("%s %.4f", Conf.toString(assignments), score);
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
	public static class MultiSplitter {
		
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
		public MultiSplitter(ConfSearch confs) {
			
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


	/**
	 * Lets exactly two consumers read confs from the stream, where one consumer
	 * always reads before the other.
	 *
	 * Supports external memory for the conformation buffer
	 */
	public static class Splitter {

		public static class OutOfOrderException extends RuntimeException {
			public OutOfOrderException() {
				super("second reader tried to read confs before first reader");
			}
		}

		public final ConfSearch confs;
		public final ConfSearch first;
		public final ConfSearch second;

		private Queue.FIFO<ScoredConf> buf;

		public Splitter(ConfSearch confs) {
			this(confs, false, null);
		}

		public Splitter(ConfSearch confs, boolean useExternalMemory, RCs rcs) {

			this.confs = confs;

			if (useExternalMemory) {
				buf = Queue.ExternalFIFOFactory.of(new ScoredConfFIFOSerializer(rcs));
			} else {
				buf = Queue.FIFOFactory.of();
			}

			AtomicBoolean exhausted = new AtomicBoolean(false);

			first = new ConfSearch() {

				@Override
				public ScoredConf nextConf() {

					// read from the ConfSearch
					ScoredConf conf = confs.nextConf();
					if (conf == null) {

						// I am le tired
						exhausted.set(true);

						return null;
					}

					// add to the buffer
					buf.push(conf);

					return conf;
				}

				@Override
				public BigInteger getNumConformations() {
					return confs.getNumConformations();
				}
			};

			second = new ConfSearch() {

				@Override
				public ScoredConf nextConf() {

					// read from the buffer
					ScoredConf conf = buf.poll();

					if (conf == null && !exhausted.get()) {
						throw new OutOfOrderException();
					}

					return conf;
				}

				@Override
				public BigInteger getNumConformations() {
					return confs.getNumConformations();
				}
			};
		}
	}
}
