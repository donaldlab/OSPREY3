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

package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.tools.ObjectIO;

public class KSAllowedSeqs {

	private DEEPerSettings dset;
	private ArrayList<String[]> moveableStrandTermini;
	private ArrayList<String[]> freeBBZoneTermini;
	private ArrayList<String> flexRes;
	private ArrayList<ArrayList<String>> allowedAAs;
	private ArrayList<String> wt;
	public boolean addWT;
	private int dist;
        private boolean allowLessMut;
	private int strand;
	private ResidueTermini limits;
	private int maxSequences = (int)Math.pow(2, 28);
	private ArrayList<ArrayList<String>> allowedSeqs = null;
	private LinkedHashMap<ArrayList<String>, Integer> allowedSeq2Index = null;
	ArrayList<HashSet<ArrayList<String>>> allowedSubSeqs = null;

	public KSAllowedSeqs( int strand, ResidueTermini limits, DEEPerSettings dset, 
			ArrayList<String[]> freeBBZoneTermini,
			ArrayList<String[]> moveableStrandTermini,
			ArrayList<String> flexRes, 
			ArrayList<ArrayList<String>> allowedAAs, 
			ArrayList<String> wt, boolean addWT, int dist, boolean allowLessMut ) {

		this.strand = strand;
		this.limits = limits;
		this.dset = dset;
		this.freeBBZoneTermini = freeBBZoneTermini;
		this.moveableStrandTermini = moveableStrandTermini;
		this.flexRes = flexRes;
		this.allowedAAs = addPosToAllowedAAs(allowedAAs, flexRes);		
		this.wt = addPosToSeq(wt, flexRes);
		this.addWT = addWT;
		this.dist = dist;
                this.allowLessMut = allowLessMut;
		this.allowedSeqs = generateSequences();
	}


	public KSAllowedSeqs( int strand, ResidueTermini limits, DEEPerSettings dset, 
			ArrayList<String[]> freeBBZoneTermini,
			ArrayList<String[]> moveableStrandTermini,
			ArrayList<String> flexRes, KSAllowedSeqs in, 
			ArrayList<ArrayList<String>> allowedAAs, int lb, int ub ) {

		this.strand = strand;
		this.limits = limits;
		this.dset = dset;
		this.freeBBZoneTermini = freeBBZoneTermini;
		this.moveableStrandTermini = moveableStrandTermini;
		this.flexRes = flexRes;
		this.allowedAAs = allowedAAs;
		this.wt = new ArrayList<String>( in.wt.subList(lb, ub) );
		this.addWT = in.addWT;
		this.dist = in.dist;
		this.allowedSeqs = new ArrayList<ArrayList<String>>();

		// filter allowedSequences by index limits
		for(ArrayList<String> seq : in.allowedSeqs) {
			this.allowedSeqs.add( new ArrayList<String>( seq.subList(lb, ub) ) );
		}
	}
	
	
	public ResidueTermini getStrandLimits() {
		return limits;
	}
	

	public static ArrayList<ArrayList<String>> removePosFromAllowedAAs(ArrayList<ArrayList<String>> in) {
		@SuppressWarnings("unchecked")
		ArrayList<ArrayList<String>> ans = (ArrayList<ArrayList<String>>) ObjectIO.deepCopy(in);

		for(ArrayList<String> list : ans) {
			for(int i = 0; i < list.size(); ++i) {
				list.set(i, list.get(i).split("-")[0]);
			}
		}

		return ans;
	}


	public static ArrayList<ArrayList<String>> addPosToAllowedAAs(ArrayList<ArrayList<String>> in, ArrayList<String> flexRes) {

		if(in.size() > flexRes.size())
			throw new RuntimeException("ERROR: cannot assign positions to all AAs in list.");

		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		for(int i = 0; i < in.size(); ++i) {
			ans.add(addSinglePosToSeq(in.get(i), flexRes.get(i)));
		}

		ans.trimToSize();
		return ans;
	}


	private static ArrayList<String> addSinglePosToSeq(ArrayList<String> in, String pos) {

		ArrayList<String> ans = new ArrayList<>(); 
		for(int i = 0; i < in.size(); ++i) ans.add(in.get(i) + "-" + pos);

		ans.trimToSize();
		return ans;
	}


	public static ArrayList<String> addPosToSeq(ArrayList<String> in, ArrayList<String> flexRes) {

		if(in.size() > flexRes.size())
			throw new RuntimeException("ERROR: cannot assign positions to all AAs in list.");

		ArrayList<String> ans = new ArrayList<>(); 
		for(int i = 0; i < in.size(); ++i) ans.add(in.get(i) + "-" + flexRes.get(i));

		return ans;
	}


	public static ArrayList<String> getAAsFromSeq(ArrayList<String> seq) {

		ArrayList<String> ans = new ArrayList<>();

		for(String res : seq) ans.add(res.split("-")[0]);

		ans.trimToSize();
		return ans;
	}


	public static ArrayList<String> getFlexResFromSeq(ArrayList<String> seq) {

		ArrayList<String> ans = new ArrayList<>();

		for(String res : seq) ans.add(res.split("-")[1]);

		ans.trimToSize();
		return ans;
	}


	public ArrayList<Integer> getFlexResIndexesFromSeq(ArrayList<String> seq) {

		ArrayList<String> seqFlexRes = getFlexResFromSeq(seq);
		ArrayList<Integer> ans = new ArrayList<>();

		for(String pos : seqFlexRes) {
			
			for(int j = 0; j < flexRes.size(); ++j) {
				if(pos.equalsIgnoreCase(flexRes.get(j))) {
					ans.add(j);
					break;
				}
			}
		}

		ans.trimToSize();
		return ans;
	}

	
	public int getFlexPosIndex( String res ) {
		
		String flexPos = res.split("-")[1];
		
		for(int index = 0; index < flexRes.size(); ++index) {		
			if(flexPos.equalsIgnoreCase(flexRes.get(index))) return index;
		}
		return Integer.MIN_VALUE;
	}
	
	
	public boolean isAllowed( String res ) {
		int pos;
		if((pos = getFlexPosIndex(res)) == Integer.MIN_VALUE) return false;
		
		return allowedAAs.get(pos).contains(res);
	}
	
	
	public boolean isAllowed( ArrayList<String> seq ) {
		for(String res : seq) {
			if(!isAllowed(res)) return false;
		}
		return true;
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


	public ArrayList<String> getFlexRes() {
		return flexRes;
	}


	public int getSequenceLength() {
		return flexRes.size();
	}


	public ArrayList<ArrayList<String>> getAllowedAAs() {
		return allowedAAs;
	}


	public ArrayList<HashSet<ArrayList<String>>> getStrandSubSeqList() {

		if( allowedSubSeqs == null) {

			allowedSubSeqs = new ArrayList<>();

			// create arraylist for all depths except the last
			// depth 0 is empty
			for( int depth = 0; depth < getSequenceLength(); ++depth ) {
				allowedSubSeqs.add( new HashSet<ArrayList<String>>() );
			}

			for( ArrayList<String> seq : allowedSeqs ) {

				// make all subsequences of seq
				for( int i = 1; i < getSequenceLength(); ++i ) {

					ArrayList<String> tmp = new ArrayList<>(seq.subList(0, i));

					if( !allowedSubSeqs.get(i).contains(tmp) )
						allowedSubSeqs.get(i).add(tmp);
				}

			}

			// add fully defined sequences at final depth
			HashSet<ArrayList<String>> fullyDefSeqs = new HashSet<>(getStrandSeqList());

			allowedSubSeqs.add( fullyDefSeqs );
		}

		return allowedSubSeqs;
	}
	
	
	public ArrayList<HashSet<ArrayList<String>>> getStrandSubSeqList2( 
			KSAllowedSeqs p, KSAllowedSeqs l ) {

		if( strand != 2 )
			throw new RuntimeException("ERROR: this version of the method "
					+ "should only be called for the COMPLEX strand");

		if( allowedSubSeqs == null ) {

			allowedSubSeqs = new ArrayList<>();

			// create arraylist for all depths, including the last. this is because
			// p and l already contain the last depth
			// depth 0 is empty
			for( int depth = 0; depth <= getSequenceLength(); ++depth ) {
				allowedSubSeqs.add( new HashSet<ArrayList<String>>() );
			}

			for( int depth = 1; depth <= Math.max(p.getStrandSubSeqsMaxDepth(), l.getStrandSubSeqsMaxDepth()); ++depth ) {
				// adjust indices to prevent out of bounds error
				int depthP = Math.min(depth, p.getStrandSubSeqsMaxDepth());
				int depthL = Math.min(depth, l.getStrandSubSeqsMaxDepth());
				
				for( ArrayList<String> subSeqP : p.getStrandSubSeqsAtDepth(depthP) ) {
					
					for( ArrayList<String> subSeqL : l.getStrandSubSeqsAtDepth(depthL) ) {
						
						if( p.getDistFromWT(subSeqP) + l.getDistFromWT(subSeqL) <= dist ) {

							ArrayList<String> tmpSubSeq = new ArrayList<>();

							tmpSubSeq.addAll(subSeqP);
							tmpSubSeq.addAll(subSeqL);

							if( tmpSubSeq.size() == getSequenceLength() && getDistFromWT(tmpSubSeq) != dist ) continue;
							
							tmpSubSeq.trimToSize();
							
							// add complex subsequence
							if( !allowedSubSeqs.get(tmpSubSeq.size()).contains(tmpSubSeq) ) 
								allowedSubSeqs.get(tmpSubSeq.size()).add(tmpSubSeq);
						}
					}
				}
			}
		}
		
		return allowedSubSeqs;
	}


	public HashSet<ArrayList<String>> getStrandSubSeqsAtDepth( int depth, KSAllowedSeqs p, KSAllowedSeqs l ) {

		if( strand != 2 )
			throw new RuntimeException("ERROR: this version of the method "
					+ "should only be called for the COMPLEX strand");

		if( allowedSubSeqs == null )
			getStrandSubSeqList2(p, l);

		if( depth < 0 || depth > allowedSubSeqs.size()-1 )
			throw new RuntimeException("ERROR: the requested depth " + depth + 
					" is not within the valid range [0," + (allowedSubSeqs.size()-1) + "]");

		return allowedSubSeqs.get( depth );
	}


	public HashSet<ArrayList<String>> getStrandSubSeqsAtDepth( int depth ) {
		
		if( allowedSubSeqs == null ) {
			
			if(strand == 2) 
				throw new RuntimeException("ERROR: sub-sequences of the COMPLEX "
						+ "strand cannot be initialized using this method");
			
			getStrandSubSeqList();
		}

		if( depth < 0 || depth > allowedSubSeqs.size()-1 )
			throw new RuntimeException("ERROR: the requested depth " + depth + 
					" is not within the valid range [0," + (allowedSubSeqs.size()-1) + "]");

		return allowedSubSeqs.get( depth );
	}
	
	
	public static void deleteFromSet( ArrayList<String> item, HashSet<ArrayList<String>> set ) {
		// delete item from any matching element of set
		for( Iterator<ArrayList<String>> iterator = set.iterator(); iterator.hasNext(); ) {
			
			ArrayList<String> element = iterator.next();
			
			if(element.containsAll(item))
				iterator.remove();
		}
	}

	
	public int getNumSubSeqs() {
		int ans = 0;
		for(int depth = 0; depth <= getStrandSubSeqsMaxDepth(); ++depth)
			ans += getStrandSubSeqsAtDepth(depth).size();
		
		return ans;
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

			aasAtPos.trimToSize();
			newAllowedAAs.add(aasAtPos);
		}

		newAllowedAAs.trimToSize();
		allowedAAs = newAllowedAAs;
	}


	public int getNumSeqs() {
		return allowedSeqs.size();
	}


	public ArrayList<ArrayList<String>> getStrandSeqList() {
		return allowedSeqs;
	}
	
	
	public int getPosOfSeq( ArrayList<String> seq ) {
		
		if( allowedSeq2Index == null ) {
			allowedSeq2Index = new LinkedHashMap<>();
			
			for( int index = 0; index < allowedSeqs.size(); ++index )
				allowedSeq2Index.put(allowedSeqs.get(index), index);
		}
			
		return allowedSeq2Index.get(seq);
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


	private ArrayList<ArrayList<String>> generateAllSequencesWithDist ( 
			ArrayList<ArrayList<String>> input ) {

		// pre-allocate buffer and wt
		ArrayList<String> buffer = new ArrayList<>();
		for ( int it = 0; it < input.size(); it++ ) buffer.add("");

		buffer.trimToSize();

		// linked hashset to preserve order
		LinkedHashSet<ArrayList<String>> output = new LinkedHashSet<>();

		generatePermutations( input, output, buffer, 0, 0 );

		// remove wt, if present
		boolean wtIsPresent = output.remove(wt);

		ArrayList<ArrayList<String>> ans = new ArrayList<ArrayList<String>>(output); 
		if( addWT || wtIsPresent ) ans.add(0, wt);
		
		int size = dist == 0 ? ans.size() : ans.size()-1;
		System.out.println("\nNumber of sequences with " + dist + " mutation(s) from wild type: " + size + "\n");

		ans.trimToSize();
		return ans;
	}


	private void generatePermutations( ArrayList<ArrayList<String>> input, 
			LinkedHashSet<ArrayList<String>> output, ArrayList<String> current, 
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

			if(!output.contains(current)) {
				ArrayList<String> ans = new ArrayList<String>(current);
				ans.trimToSize();
				output.add(ans);
			}

			return;
		}

		if(depth == input.size()) {
			if( (diff==dist||allowLessMut) && !output.contains(current)) {
				ArrayList<String> ans = new ArrayList<String>(current);
				ans.trimToSize();
				output.add(ans);
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
	protected boolean isSpecifiedDist( ArrayList<String> s1, ArrayList<String> s2 ) {

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


	public ArrayList<String> getWTSeq() {
		return wt;
	}
	
	
	public boolean containsWTSeq() {
		if( allowedSeqs == null || allowedSeqs.size() < 1 ) 
			return false;
		
		return allowedSeqs.get(0).equals(getWTSeq());
	}
}
