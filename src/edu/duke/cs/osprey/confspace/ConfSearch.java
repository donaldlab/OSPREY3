/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections4.queue.CircularFifoQueue;

//This is a general interface for things that search conformational space
//like an A* trees, a BWM* search, or a WCSP solver
//For each of these, we instantiate the ConfSearch, 
//call nextConf to get the GMEC for the model being searched,
//and then successive calls to nextConf() return the next conformation in the gap-free list
/**
 *
 * @author mhall44
 */
public interface ConfSearch {
    
    BigInteger getNumConformations();
    ScoredConf nextConf();
    List<ScoredConf> nextConfs(double maxEnergy);
    
    public static class ScoredConf {
        
        private int[] assignments;
        private double score;
        
        public ScoredConf(int[] assignments, double score) {
            this.assignments = assignments;
            this.score = score;
        }
        
        public ScoredConf(ScoredConf other) {
        	this.assignments = other.assignments.clone();
        	this.score = other.score;
        }
        
        public int[] getAssignments() {
            return assignments;
        }
        
        public double getScore() {
            return score;
        }
        public void setScore(double val) {
        	score = val;
        }
        public void offsetScore(double val) {
        	score += val;
        }
    }
    
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
        
        public double getEnergy() {
            return energy;
        }
        public void setEnergy(double val) {
        	energy = val;
        }
        public void offsetEnergy(double val) {
        	energy += val;
        }
    }
    
	// lets multiple consumers read confs from the stream regardless of order of reads
	public static class Splitter {
		
		public class Stream {
			
			private int index;
			
			private Stream() {
				index = 0;
			}
			
			public ScoredConf next() {
			
				// where in the buffer should we read?
				int pos = index - firstIndex;
				
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
					conf = buf.get(pos);
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
				
				// TEMP
				for (Stream other : streams) {
					assert (other.index >= firstIndex);
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
		}
		
		private ConfSearch confs;
		private CircularFifoQueue<ScoredConf> buf;
		private int firstIndex;
		private List<Stream> streams;
		
		public Splitter(ConfSearch confs) {
			
			this.confs = confs;
			
			buf = new CircularFifoQueue<>();
			firstIndex = 0;
			streams = new ArrayList<>();
		}
		
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
