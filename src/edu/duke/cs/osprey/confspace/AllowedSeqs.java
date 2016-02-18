package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;

import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;

public class AllowedSeqs {

	private DEEPerSettings dset;
	private ArrayList<String[]> moveableStrandTermini;
	private ArrayList<String[]> freeBBZoneTermini;
	private ArrayList<String> flexRes;
	private ArrayList<ArrayList<String>> allowedAAs;
	private ArrayList<String> wt;
	private int dist;
	private int maxSequences = (int)Math.pow(2, 27);
	ArrayList<ArrayList<String>> allowedSeqs;

	public AllowedSeqs( DEEPerSettings dset, ArrayList<String[]> freeBBZoneTermini,
			ArrayList<String[]> moveableStrandTermini,
			ArrayList<String> flexRes, 
			ArrayList<ArrayList<String>> allowedAAs, 
			ArrayList<String> wt, int dist ) {

		this.dset = dset;
		this.freeBBZoneTermini = freeBBZoneTermini;
		this.moveableStrandTermini = moveableStrandTermini;
		this.flexRes = flexRes;
		this.allowedAAs = allowedAAs;
		this.wt = wt;
		this.dist = dist;
		this.allowedSeqs = generateSequences();
	}


	public AllowedSeqs( DEEPerSettings dset, ArrayList<String[]> freeBBZoneTermini,
			ArrayList<String[]> moveableStrandTermini,
			ArrayList<String> flexRes, AllowedSeqs in, 
			ArrayList<ArrayList<String>> allowedAAs, int lb, int ub ) {

		this.dset = dset;
		this.freeBBZoneTermini = freeBBZoneTermini;
		this.moveableStrandTermini = moveableStrandTermini;
		this.flexRes = flexRes;
		this.allowedAAs = allowedAAs;
		this.wt = in.wt;
		this.dist = in.dist;
		this.allowedSeqs = new ArrayList<ArrayList<String>>();

		// filter allowedSequences by index limits
		for(ArrayList<String> seq : in.allowedSeqs) {
			this.allowedSeqs.add( new ArrayList<String>( seq.subList(lb, ub) ) );
		}
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


	public ArrayList<String> getFlexRes() {
		return flexRes;
	}
	
	
	public int getSequenceLength() {
		return flexRes.size();
	}


	public ArrayList<ArrayList<String>> getAllowedAAs() {
		return allowedAAs;
	}


	public int getNumSeqs() {
		return allowedSeqs.size();
	}


	public ArrayList<ArrayList<String>> getStrandSeqList() {
		return allowedSeqs;
	}
	
	
	public ArrayList<String> getStrandSeq(int index) {
		if(index > -1 && index < allowedSeqs.size()) 
			return allowedSeqs.get(index);
		return null;
	}
	
	
	public void removeStrandSeq(int index) {
		if(index > -1 && index < allowedSeqs.size()) 
			allowedSeqs.remove(index);
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
