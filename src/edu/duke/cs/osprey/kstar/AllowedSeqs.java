package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;

public class AllowedSeqs {

	private DEEPerSettings dset;
	private ArrayList<String[]> moveableStrandTermini;
	private ArrayList<String[]> freeBBZoneTermini;
	private ArrayList<String> flexRes;
	private ArrayList<ArrayList<String>> allowedAAs;
	private ArrayList<String> wt;
	private int dist;
	private int strand;
	private int maxSequences = (int)Math.pow(2, 27);
	ArrayList<ArrayList<String>> allowedSeqs = null;
	ArrayList<ArrayList<ArrayList<String>>> allowedSubSeqs = null;
	HashMap<ArrayList<String>, ArrayList<String>> seq2FlexRes = new HashMap<>();

	public AllowedSeqs( int strand, DEEPerSettings dset, 
			ArrayList<String[]> freeBBZoneTermini,
			ArrayList<String[]> moveableStrandTermini,
			ArrayList<String> flexRes, 
			ArrayList<ArrayList<String>> allowedAAs, 
			ArrayList<String> wt, int dist ) {

		this.strand = strand;
		this.dset = dset;
		this.freeBBZoneTermini = freeBBZoneTermini;
		this.moveableStrandTermini = moveableStrandTermini;
		this.flexRes = flexRes;
		this.allowedAAs = allowedAAs;
		this.wt = wt;
		this.dist = dist;
		this.allowedSeqs = generateSequences();
	}


	public AllowedSeqs( int strand, DEEPerSettings dset, 
			ArrayList<String[]> freeBBZoneTermini,
			ArrayList<String[]> moveableStrandTermini,
			ArrayList<String> flexRes, AllowedSeqs in, 
			ArrayList<ArrayList<String>> allowedAAs, int lb, int ub ) {

		this.strand = strand;
		this.dset = dset;
		this.freeBBZoneTermini = freeBBZoneTermini;
		this.moveableStrandTermini = moveableStrandTermini;
		this.flexRes = flexRes;
		this.allowedAAs = allowedAAs;
		this.wt = new ArrayList<String>( in.wt.subList(lb, ub) );
		this.dist = in.dist;
		this.allowedSeqs = new ArrayList<ArrayList<String>>();

		// filter allowedSequences by index limits
		for(ArrayList<String> seq : in.allowedSeqs) {
			this.allowedSeqs.add( new ArrayList<String>( seq.subList(lb, ub) ) );
		}
	}


	public int getStrand() {
		return strand;
	}


	public int getDistFromWT( ArrayList<String> seq ) {

		if( seq.size() > wt.size() )
			throw new RuntimeException("ERROR: sequence length > wt sequence length");

		int dist = 0;
		for( int i = 0; i < seq.size(); ++i ) {
			if( seq.get(i).compareTo(wt.get(i)) != 0 ) {
				dist++;
			}
		}

		return dist;
	}

	public DEEPerSettings getDEEPerSettings() {
		return dset;
	}

	public ArrayList<String[]> getMoveableStrandTermini() {
		return moveableStrandTermini;
	}


	public ArrayList<String[]> getFreeBBZoneTermini() {
		return freeBBZoneTermini;
	}


	public ArrayList<String> getFlexRes( ArrayList<String> seq ) {

		ArrayList<String> ans = null;

		if( seq2FlexRes.get(seq) == null ) {
			
			ans = new ArrayList<>();

			for( int i = 0; i < seq.size(); ++i ) {

				if( !allowedAAs.get(i).contains(seq.get(i)) )
					throw new RuntimeException("ERROR: amino acid " + seq.get(i) + 
							" is not allowed at position " + flexRes.get(i));

				ans.add( flexRes.get(i) );
			}

			seq2FlexRes.put(seq, ans);
		}

		else ans = seq2FlexRes.get(seq);

		return ans;
	}


	public ArrayList<String> getFlexRes() {
		return flexRes;
	}


	public int getSequenceLength() {
		return flexRes.size();
	}


	public ArrayList<ArrayList<String>> getAllowedAAs() {
		return allowedAAs;
	}


	public ArrayList<ArrayList<ArrayList<String>>> getStrandSubSeqList() {

		if( allowedSubSeqs == null) {

			allowedSubSeqs = new ArrayList<>();

			// create arraylist for all depths except the last
			// depth 0 is empty
			for( int depth = 0; depth < getSequenceLength(); ++depth ) {
				allowedSubSeqs.add( new ArrayList<ArrayList<String>>() );
			}

			for( ArrayList<String> seq : allowedSeqs ) {

				// make all subsequences of seq
				for( int i = 1; i < getSequenceLength(); ++i ) {

					ArrayList<String> tmp = new ArrayList<>(seq.subList(0, i));

					//tmp.addAll( seq.subList(0, i) );

					if( !allowedSubSeqs.get(i).contains(tmp) )
						allowedSubSeqs.get(i).add(tmp);
				}

			}

			// add sequences at final depth
			// since we are creating a successor function, we cannot allow duplicates
			ArrayList<ArrayList<String>> fullyDefSeqs = new ArrayList<>();
			for( ArrayList<String> seq : getStrandSeqList() ) {
				if( !fullyDefSeqs.contains(seq) )
					fullyDefSeqs.add(seq);
			}

			allowedSubSeqs.add( fullyDefSeqs );
		}

		return allowedSubSeqs;
	}


	public ArrayList<ArrayList<ArrayList<String>>> getStrandSubSeqList( 
			AllowedSeqs p, AllowedSeqs l ) {

		if( strand != Strand.COMPLEX )
			throw new RuntimeException("ERROR: this version of the method "
					+ "should only be called for the COMPLEX strand");

		if( allowedSubSeqs == null ) {

			allowedSubSeqs = new ArrayList<>();

			// create arraylist for all depths, including the last. this is because
			// p and l already contain the last depth
			// depth 0 is empty
			for( int depth = 0; depth <= getSequenceLength(); ++depth ) {
				allowedSubSeqs.add( new ArrayList<ArrayList<String>>() );
			}

			for( int depthP = 1; depthP <= p.getStrandSubSeqsMaxDepth(); ++depthP ) {

				for( ArrayList<String> subSeqP : p.getStrandSubSeqsAtDepth(depthP) ) {

					for( int depthL = 1; depthL <= l.getStrandSubSeqsMaxDepth(); ++depthL ) {

						for( ArrayList<String> subSeqL : l.getStrandSubSeqsAtDepth(depthL) ) {

							if( p.getDistFromWT(subSeqP) + l.getDistFromWT(subSeqL) <= dist ) {

								ArrayList<String> tmpSubSeq = new ArrayList<>();
								ArrayList<String> tmpFlexRes = new ArrayList<>();

								tmpSubSeq.addAll(subSeqP);
								tmpSubSeq.addAll(subSeqL);

								tmpFlexRes.addAll(p.getFlexRes(subSeqP));
								tmpFlexRes.addAll(l.getFlexRes(subSeqL));

								// add complex subsequence
								if( !allowedSubSeqs.get(tmpSubSeq.size()).contains(tmpSubSeq) ) 
									allowedSubSeqs.get(tmpSubSeq.size()).add(tmpSubSeq);
								
								// add complex flexible res positions for subsequence
								seq2FlexRes.put(tmpSubSeq, tmpFlexRes);
							}
						}		
					}
				}
			}

			/*
			// i believe all fully defined sequences should have been added from the above loop.
			// if not, then run the commented code
			ArrayList<ArrayList<String>> finalDepth = allowedSubSeqs.get(allowedSubSeqs.size()-1);
			for( ArrayList<String> seq : getStrandSeqList() ) {
				if( !finalDepth.contains(seq) )
					finalDepth.add(seq);
			}
			 */
		}

		return allowedSubSeqs;
	}


	public ArrayList<ArrayList<String>> getStrandSubSeqsAtDepth( int depth, AllowedSeqs p, AllowedSeqs l ) {

		if( strand != Strand.COMPLEX )
			throw new RuntimeException("ERROR: this version of the method "
					+ "should only be called for the COMPLEX strand");

		if( allowedSubSeqs == null )
			getStrandSubSeqList(p, l);

		if( depth < 0 || depth > allowedSubSeqs.size()-1 )
			throw new RuntimeException("ERROR: the requested depth " + depth + 
					" is not within the valid range [0," + (allowedSubSeqs.size()-1) + "]");

		return allowedSubSeqs.get( depth );
	}


	public ArrayList<ArrayList<String>> getStrandSubSeqsAtDepth( int depth ) {

		if( allowedSubSeqs == null )
			getStrandSubSeqList();

		if( depth < 0 || depth > allowedSubSeqs.size()-1 )
			throw new RuntimeException("ERROR: the requested depth " + depth + 
					" is not within the valid range [0," + (allowedSubSeqs.size()-1) + "]");

		return allowedSubSeqs.get( depth );
	}


	public int getStrandSubSeqsMaxDepth() {
		return getFlexRes().size();
	}


	public void truncateAllowedAAs() {

		ArrayList<ArrayList<String>> newAllowedAAs = new ArrayList<>();

		for(int i = 0; i < getSequenceLength(); ++i) {

			ArrayList<String> aasAtPos = new ArrayList<>();

			for(ArrayList<String> al : getStrandSeqList()) {
				String aa = al.get(i);
				if(!aasAtPos.contains(aa))
					aasAtPos.add(aa);
			}

			newAllowedAAs.add(aasAtPos);
		}

		allowedAAs = newAllowedAAs;
	}


	public int getNumSeqs() {
		return allowedSeqs.size();
	}


	public ArrayList<ArrayList<String>> getStrandSeqList() {
		return allowedSeqs;
	}


	public ArrayList<String> getStrandSeqAtPos( int index ) {
		if(index > -1 && index < allowedSeqs.size()) 
			return allowedSeqs.get(index);
		return null;
	}


	public ArrayList<String> getStrandSubSeq( int index, int begin, int end ) {
		if(begin < 0 || end > getSequenceLength()) {
			throw new RuntimeException("ERROR: begin and end indexes [" + begin + "," + end + "] are out of range."
					+ " Valid range is [0," + getSequenceLength() + "].");
		}

		ArrayList<String> seq = getStrandSeqAtPos(index);
		if(seq == null) return null;

		ArrayList<String> ans = new ArrayList<>();
		for(int i = begin; i < end; ++i) {
			ans.add(seq.get(i));
		}

		return ans;
	}


	public ArrayList<String> removeStrandSeq(int index) {
		if(index > -1 && index < allowedSeqs.size()) 
			return allowedSeqs.remove(index);
		else
			throw new RuntimeException("ERROR: index "+ index + " is invalid. "
					+ "Valid range is [0," + (getSequenceLength()-1) +"]");
	}


	private ArrayList<ArrayList<String>> generateSequences() {
		return generateAllSequencesWithDist( allowedAAs );
	}


	/**
	 * Generates all sequences that differ from WT in exactly dist positions.
	 * WT sequence is assumed to be in the first column of the input.
	 * @param input
	 * @param dist
	 * @return
	 */
	private ArrayList<ArrayList<String>> generateAllSequencesWithDist ( 
			ArrayList<ArrayList<String>> input ) {

		// pre-allocate buffer and wt
		ArrayList<String> buffer = new ArrayList<>();
		for ( int it = 0; it < input.size(); it++ ) {
			buffer.add("");
		}
		buffer.trimToSize();

		ArrayList<ArrayList<String>> output = new ArrayList<>();
		output.ensureCapacity(maxSequences);

		generatePermutations( input, output, buffer, 0, 0 );

		// remove objects that differ from wt by more than dist elements
		int numRemoved = 0;
		for( int it = 0; it < output.size(); ) {
			if( !isSpecifiedDist(wt, output.get(it)) ) {
				output.remove(it);
				numRemoved++;
			}
			else ++it;
		}

		if( numRemoved > 0 ) {
			System.out.println("Error check on created mutants removed " + numRemoved + " sequences");
		}

		System.out.println("Number of sequences with " + this.dist + 
				" mutation(s) from wild type: " + output.size());

		output.trimToSize();
		return output;
	}


	private void generatePermutations( ArrayList<ArrayList<String>> input, 
			ArrayList<ArrayList<String>> output, ArrayList<String> current, 
			int depth, int diff ) {

		if( output.size() >= maxSequences )
			throw new RuntimeException("ERROR: the number of requested sequence "
					+ "combinations exceeds " + maxSequences + ". Reduce the value of the NUMMUTATIONS parameter.");

		if(diff == dist) {
			// if we have the max num of mutants but have not assigned all residue
			// positions, then take wt for the rest of the positions
			if(depth != input.size()) {
				for( ; depth < wt.size(); depth++) {
					current.set(depth, wt.get(depth));
				}
			}

			if(!output.contains(current))
				output.add(new ArrayList<String>(current));

			return;
		}

		if(depth == input.size()) {
			if(diff == dist) {
				output.add(new ArrayList<String>(current));
			}
			return;
		}

		for( int it = 0; it < input.get(depth).size(); ++it ) {

			current.set(depth, input.get(depth).get(it));

			int tmpDiff = current.get(depth).equalsIgnoreCase(wt.get(depth)) ? diff : diff + 1;
			if( tmpDiff > dist ) continue;

			generatePermutations( input, output, current, depth + 1, tmpDiff );
		}
	}


	/**
	 * Returns true iff the two sequences (of the same length) differ by exactly
	 * inDist number of elements
	 * @param wtSeq
	 * @param seq
	 * @param dist
	 * @return
	 */
	private boolean isSpecifiedDist( ArrayList<String> s1, ArrayList<String> s2 ) {

		if( s1.size() != s2.size() )
			throw new RuntimeException("Error: input strings " + s1 + " and " 
					+ s2 + " are not the same length.");

		int dist = 0;
		for( int it = 0 ; it < s1.size(); ++it ) {
			if( !s1.get(it).equalsIgnoreCase(s2.get(it)) ) ++dist;
			if( dist > this.dist ) return false;
		}

		return dist == this.dist ? true : false;
	}


	public boolean containsWT() {
		boolean ans = allowedSeqs.contains(wt);
		return ans;
	}


	public void addWT() {
		if( containsWT() )
			allowedSeqs.remove(wt);

		allowedSeqs.add(0, wt);
	}
}
