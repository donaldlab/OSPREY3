package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Objects;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Function;

public abstract class AbstractTupleMatrix<T> implements TupleMatrix<T>, Serializable {

	private static final long serialVersionUID = 1654821458320722522L;
	
    //We will need "matrices" of quantities defined
    //for example, the energy matrix (T=Double) stores single, pairwise, and higher-order energies
    //and we'll also have pruning (T=Boolean) and EPIC (T=EPoly) matrices
    //we'll store things as ArrayLists to make it easier to merge residues, partition RCs, etc.
    //and also to facilitate generics
	
	private int numPos; // eg residues
	private int[] numConfAtPos; // eg RCs at each residue
    
    // index the arrays
    private int[] oneBodyOffsets;
    private int[] pairwiseOffsets;
    private int numPairwiseTerms;
    
    // TODO: do we really need a pruning interval to make an energy matrix?
    // or can we simplify the code overall by storing the pruning interval somewhere else?
    private double pruningInterval;//This matrix needs to hold entries for all RCs
    //that cannot be pruned with the specified pruning interval (Ew + Ival)
    //i.e. the matrix must describe all conformations within pruningInterval 
    //of the lowest pairwise lower bound
    
	private ArrayList<HigherTupleFinder<T>> higherTerms; // indices: same as pairwise, can be null if no interactions
    private T defaultHigherInteraction;//We only mark sparse higher interactions;
    //if unmarked we assume this value (e.g., 0 for energy, false for pruning)
    
    
    protected AbstractTupleMatrix() {
    	// do nothing
    	// so subclasses can avoid allocation if they want
    }
    
    protected AbstractTupleMatrix(ConfSpace cSpace, double pruningInterval, T defaultHigherInteraction) {
        //allocate the matrix based on the provided conformational space
    	this(cSpace.numPos, cSpace.getNumRCsAtPos(), pruningInterval, defaultHigherInteraction);
    }
    
    protected AbstractTupleMatrix(SimpleConfSpace confSpace, double pruningInterval, T defaultHigherInteraction) {
        this(confSpace.positions.size(), confSpace.getNumResConfsByPos(), pruningInterval, defaultHigherInteraction);
    }
    
    protected AbstractTupleMatrix(int numPos, int[] numConfAtPos, double pruningInterval, T defaultHigherInteraction) {
        //allocate the matrix based on the provided conformational space size
        //also specify what pruningInterval it's valid up to
    	
    	this.pruningInterval = pruningInterval;
    	this.defaultHigherInteraction = defaultHigherInteraction;
    	
    	this.numPos = numPos;
    	this.numConfAtPos = numConfAtPos;
        
        // compute the indices and allocate space
        int oneBodyOffset = 0;
        
        // first one-body offsets
        oneBodyOffsets = new int[numPos];
        for (int res1=0; res1<numPos; res1++) {
        	oneBodyOffsets[res1] = oneBodyOffset;
        	oneBodyOffset += numConfAtPos[res1];
        }
        
        // then pairwise offsets
        pairwiseOffsets = new int[numPos*(numPos - 1)/2];
        int pairwiseOffset = 0;
        int pairwiseIndex = 0;
        for (int res1=0; res1<numPos; res1++) {
        	for (int res2=0; res2<res1; res2++) {
        		pairwiseOffsets[pairwiseIndex++] = pairwiseOffset;
        		pairwiseOffset += numConfAtPos[res1]*numConfAtPos[res2];
        	}
        }
        numPairwiseTerms = pairwiseOffset;
        assert (pairwiseIndex == pairwiseOffsets.length);
        
        allocate(oneBodyOffset, numPairwiseTerms);
        
    	// don't allocate space for higher terms right now
        // wait till we write something
        higherTerms = null;
    }
    
    protected AbstractTupleMatrix(AbstractTupleMatrix<T> other) {
    	this.numPos = other.numPos;
    	this.numConfAtPos = other.numConfAtPos.clone();
    	this.oneBodyOffsets = other.oneBodyOffsets.clone();
    	this.pairwiseOffsets = other.pairwiseOffsets.clone();
    	this.numPairwiseTerms = other.numPairwiseTerms;
    	this.pruningInterval = other.pruningInterval;
    	if (other.higherTerms != null) {
    		throw new UnsupportedOperationException("copying higher order terms isn't implemented yet");
    	}
    	this.higherTerms = null;
    	this.defaultHigherInteraction = null;
    }
    
    protected abstract void allocate(int numOneBody, int numPairwise);
    
    public double getPruningInterval() {
        return pruningInterval;
    }
    
    public void setPruningInterval(double pruningInterval) {
    	this.pruningInterval = pruningInterval;
    }
    
    public T getDefaultHigherInteraction() {
    	return defaultHigherInteraction;
    }
    
    public int getNumPos() {
    	return numPos;
    }
    
    public int getNumConfAtPos(int pos) {
    	return numConfAtPos[pos];
    }

	public int[] getNumConfAtPos() {//get em all
		return numConfAtPos;
	}
    
    protected int getOneBodyIndex(int res, int conf) {
    	return oneBodyOffsets[res] + conf;
    }
    
    private int getPairwiseIndexNoCheck(int res1, int res2) {
    	return res1*(res1 - 1)/2 + res2;
    }
    
    protected int getPairwiseIndex(int res1, int res2) {
    	
    	// res2 should be strictly less than res1
    	if (res2 > res1) {
    		int swap = res1;
    		res1 = res2;
    		res2 = swap;
    	} else if (res1 == res2) {
    		throw new Error("Can't pair residue " + res1 + " with itself");
    	}
    	
    	return getPairwiseIndexNoCheck(res1, res2);
    }
    
    protected int getPairwiseIndex(int res1, int conf1, int res2, int conf2) {
    	
    	// res2 should be strictly less than res1
    	if (res2 > res1) {
    		int swap = res1;
    		res1 = res2;
    		res2 = swap;
    		swap = conf1;
    		conf1 = conf2;
    		conf2 = swap;
    	} else if (res1 == res2) {
    		throw new Error("Can't pair residue " + res1 + " with itself");
    	}
    	
    	return pairwiseOffsets[getPairwiseIndexNoCheck(res1, res2)] + numConfAtPos[res2]*conf1 + conf2;
    }
    
    @Override
    public void fill(T val) {
		for (int res1=0; res1<getNumPos(); res1++) {
			int m1 = getNumConfAtPos(res1);
			for (int i1=0; i1<m1; i1++) {
				setOneBody(res1, i1, val);
				for (int res2=0; res2<res1; res2++) {
					int m2 = getNumConfAtPos(res2);
					for (int i2=0; i2<m2; i2++) {
						setPairwise(res1, i1, res2, i2, val);
					}
				}
			}
		}
    }
    
    @Override
    public void fill(Iterator<T> val) {
		for (int res1=0; res1<getNumPos(); res1++) {
			int m1 = getNumConfAtPos(res1);
			for (int i1=0; i1<m1; i1++) {
				setOneBody(res1, i1, val.next());
				for (int res2=0; res2<res1; res2++) {
					int m2 = getNumConfAtPos(res2);
					for (int i2=0; i2<m2; i2++) {
						setPairwise(res1, i1, res2, i2, val.next());
					}
				}
			}
		}
    }
    
    @Override
    public boolean hasHigherOrderTerms() {
    	return higherTerms != null;
    }
    
    @Override
    public void setTupleValue(RCTuple tup, T val){
        //assign the given value to the specified RC tuple
        int tupSize = tup.pos.size();
        
        if(tupSize==1)//just a one-body quantity
            setOneBody( tup.pos.get(0), tup.RCs.get(0), val);
        else if(tupSize==2)//two-body
            setPairwise( tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), val );
        else if(tupSize>2){//higher-order
            setHigherOrder(tup,val);
        }
        else
            throw new UnsupportedOperationException( "ERROR: Not supporting tuple size " + tupSize );
    }
    
    public void setHigherOrder(RCTuple tup, T val){
        //set a higher-order term
        //we need all pairs contained in tup to know about it
        
        //loop over pairs
        for(int index1=0; index1<tup.pos.size(); index1++){
            for(int index2=0; index2<index1; index2++){
                
                int pos1 = tup.pos.get(index1);
                int rc1 = tup.RCs.get(index1);
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
                
                //put tup into the HigherTupleFinder for this pair
                HigherTupleFinder<T> htf = getHigherOrderTerms(pos1,rc1,pos2,rc2);
                
                //create a HigherTupleFinder if there is none yet
                if(htf==null){
                    htf = new HigherTupleFinder<>(defaultHigherInteraction);
                    setHigherOrderTerms(pos1, rc1, pos2, rc2, htf);
                }
                
                RCTuple subTup = tup.subtractMember(index1).subtractMember(index2);
                
                htf.setInteraction(subTup,val);
            }
        }
        
    }
    
    @Override
    public HigherTupleFinder<T> getHigherOrderTerms(int res1, int conf1, int res2, int conf2) {
    	if (higherTerms != null) {
    		return higherTerms.get(getPairwiseIndex(res1, conf1, res2, conf2));
    	}
    	return null;
    }
    
    @Override
    public void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<T> val) {
    	if (val != null && higherTerms == null) {
    		
    		// lazy allocation
			higherTerms = new ArrayList<>(numPairwiseTerms);
			for (int i=0; i<numPairwiseTerms; i++) {
				higherTerms.add(null);
			}
    	}
    	if (higherTerms != null) {
    		higherTerms.set(getPairwiseIndex(res1, conf1, res2, conf2), val);
    	}
    }

	public String toString(int cellWidth, Function<T,String> formatter) {

		StringBuilder buf = new StringBuilder();

		// make sure the cell width is at least 6
		final int fCellWidth = Math.max(cellWidth, 6);

		final String spacer = "  ";
		Consumer<String> center = (val) -> {
			int length = val.length();

			if (length > fCellWidth) {
				buf.append(val.substring(0, fCellWidth));
				return;
			}

			int pad = fCellWidth - length;
			int halfPad = pad/2;
			if (pad % 2 == 0) {
				for (int i=0; i<halfPad; i++) {
					buf.append(" ");
				}
				buf.append(val);
				for (int i=0; i<halfPad; i++) {
					buf.append(" ");
				}
			} else {
				for (int i=0; i<halfPad; i++) {
					buf.append(" ");
				}
				buf.append(val);
				for (int i=0; i<halfPad+1; i++) {
					buf.append(" ");
				}
			}
		};
		Consumer<Integer> posPrinter = (pos) -> {
			buf.append(spacer);
			center.accept(String.format("%d", pos));
		};
		Consumer<Integer> rcPrinter = (rc) -> {
			buf.append(spacer);
			center.accept(String.format("%3d", rc));
		};
		BiConsumer<Integer,Integer> labelPrinter = (pos, rc) -> {
			String label = String.format("%2d:%-3d", pos, rc);
			buf.append(String.format("%s%-" + fCellWidth + "s", spacer, label));
		};
		Consumer<T> valuePrinter = (energy) -> {
			String value = formatter.apply(energy);
			buf.append(spacer);
			buf.append(String.format("%" + fCellWidth + "s", value));
		};
		Runnable blankPrinter = () -> {
			buf.append(spacer);
			for (int i=0; i<fCellWidth; i++) {
				buf.append(" ");
			}
		};

		// find the position with the most singles
		int maxNumRcs = 0;
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			maxNumRcs = Math.max(maxNumRcs, getNumConfAtPos(pos1));
		}

		// singles
		buf.append("singles:\n");
		blankPrinter.run();
		for (int rc=0; rc<maxNumRcs; rc++) {
			rcPrinter.accept(rc);
		}
		buf.append("\n");
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			int n1 = getNumConfAtPos(pos1);

			posPrinter.accept(pos1);
			for (int rc1=0; rc1<n1; rc1++) {
				valuePrinter.accept(getOneBody(pos1, rc1));
			}
			buf.append("\n");
		}

		// pairs
		buf.append("pairs:\n");
		blankPrinter.run();
		for (int pos1=0; pos1<getNumPos() - 1; pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {
				labelPrinter.accept(pos1, rc1);
			}
		}
		buf.append("\n");
		for (int pos1=1; pos1<getNumPos(); pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {

				labelPrinter.accept(pos1, rc1);

				for (int pos2=0; pos2<pos1; pos2++) {
					int n2 = getNumConfAtPos(pos2);
					for (int rc2=0; rc2<n2; rc2++) {
						valuePrinter.accept(getPairwise(pos1, rc1, pos2, rc2));
					}
				}

				buf.append("\n");
			}
		}

		return buf.toString();
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof AbstractTupleMatrix
			&& equals((AbstractTupleMatrix<?>) other);
	}

	public boolean equals(AbstractTupleMatrix<?> other) {

		if (this.getNumPos() != other.getNumPos()) {
			return false;
		}
		int n = this.getNumPos();

		for (int pos1=0; pos1<n; pos1++) {

			if (this.getNumConfAtPos(pos1) != other.getNumConfAtPos(pos1)) {
				return false;
			}
			int n1 = this.getNumConfAtPos(pos1);

			for (int rc1=0; rc1<n1; rc1++) {

				if (!Objects.equals(this.getOneBody(pos1, rc1), other.getOneBody(pos1, rc1))) {
					return false;
				}

				for (int pos2=0; pos2<pos1; pos2++) {

					if (this.getNumConfAtPos(pos2) != other.getNumConfAtPos(pos2)) {
						return false;
					}
					int n2 = this.getNumConfAtPos(pos2);

					for (int rc2=0; rc2<n2; rc2++) {

						if (!Objects.equals(this.getPairwise(pos1, rc1, pos2, rc2), other.getPairwise(pos1, rc1, pos2, rc2))) {
							return false;
						}
					}
				}
			}
		}

		return true;
	}
}
