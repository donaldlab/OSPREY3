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

package edu.duke.cs.osprey.structure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.tools.HashCalculator;

/**
 * High-performance lookups of atom pairs by bond chain length
 */
public class AtomConnectivity {

	public static class Builder {
		
		/** should 15 bonded atoms, where at least one is H, be treated as non-bonded? */
		private boolean treat15HasNonBonded = true;

		public Builder set15HasNonBonded(boolean val) {
			treat15HasNonBonded = val;
			return this;
		}
		
		public AtomConnectivity build() {
			return new AtomConnectivity(treat15HasNonBonded);
		}
	}
	
	public static class AtomPairs {
		
		public final ResidueTemplate templ1;
		public final ResidueTemplate templ2;
		
		public AtomPairs(ResidueTemplate templ1, ResidueTemplate templ2) {
			this.templ1 = templ1;
			this.templ2 = templ2;
		}
		
		private int[][][] pairsByType = new int[AtomNeighbors.Type.values().length][][];
		
		public int[][] getPairs(AtomNeighbors.Type type) {
			return pairsByType[type.ordinal()];
		}
		
		public int getNumPairs(AtomNeighbors.Type type) {
			return getPairs(type).length;
		}
		
		public String dumpPairs(AtomNeighbors.Type type) {
			StringBuilder buf = new StringBuilder();
			buf.append("[");
			for (int i=0; i<getNumPairs(type); i++) {
				if (buf.length() > 1) {
					buf.append(", ");
				}
				buf.append(Arrays.toString(getPairs(type)[i]));
			}
			buf.append("]");
			return buf.toString();
		}
	}
	
	private static class Key1 {
		
		private ResidueTemplate templ1;
		
		public Key1(ResidueTemplate templ1) {
			this.templ1 = templ1;
		}
		
		@Override
		public int hashCode() {
			return System.identityHashCode(templ1);
		}
		
		@Override
		public boolean equals(Object obj) {
			Key1 other = (Key1)obj;
			return this.templ1 == other.templ1;
		}
	}
	
	private static class Key2 {
		
		private ResidueTemplate templ1;
		private ResidueTemplate templ2;
		private boolean isForward;
		
		public Key2(ResidueTemplate templ1, ResidueTemplate templ2, boolean isForward) {
			this.templ1 = templ1;
			this.templ2 = templ2;
			this.isForward = isForward;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				System.identityHashCode(templ1),
				System.identityHashCode(templ2),
				Boolean.hashCode(isForward)
			);
		}
		
		@Override
		public boolean equals(Object obj) {
			Key2 other = (Key2)obj;
			return this.templ1 == other.templ1
				&& this.templ2 == other.templ2
				&& this.isForward == other.isForward;
		}
	}
	
	private static class KeySeparate {
		
		private ResidueTemplate templ1;
		private ResidueTemplate templ2;
		
		public KeySeparate(ResidueTemplate templ1, ResidueTemplate templ2) {
			this.templ1 = templ1;
			this.templ2 = templ2;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				System.identityHashCode(templ1),
				System.identityHashCode(templ2)
			);
		}
		
		@Override
		public boolean equals(Object obj) {
			KeySeparate other = (KeySeparate)obj;
			return this.templ1 == other.templ1 && this.templ2 == other.templ2;
		}
	}


	public final boolean treat15HasNonBonded;
	
	private final Map<Key1,AtomPairs> atomPairs1 = new HashMap<>();
	private final Map<Key2,AtomPairs> atomPairs2 = new HashMap<>();
	private final Map<KeySeparate,AtomPairs> atomPairsSeparate = new HashMap<>();
	
	private AtomConnectivity(boolean treat15HasNonBonded) {
		this.treat15HasNonBonded = treat15HasNonBonded;
	}

	public AtomPairs getAtomPairs(Residue res1, Residue res2) {

		// NOTE: remember, don't safe references to these residues in the caches anywhere

		// NOTE: also, this function gets hammered A LOT by multiple threads,
		// so be thread-safe when doing writes to shared data structures
		
		// do we want intra pairs?
		if (res1 == res2) {
			synchronized (atomPairs1) {
				return atomPairs1.computeIfAbsent(new Key1(res1.template), key ->
					makeAtomPairs(res1, res1)
				);
			}
		}

		/* TEMP
		// are they bonded together?
		if (isInterResBondedForward(res1, res2)) {
			// yup, in forward order
			synchronized (atomPairs2) {
				return atomPairs2.computeIfAbsent(new Key2(res1.template, res2.template, true), key ->
					makeAtomPairs(res1, res2)
				);
			}
		} else if (isInterResBondedForward(res2, res1)) {
			// yup, in reverse order
			synchronized (atomPairs2) {
				return atomPairs2.computeIfAbsent(new Key2(res2.template, res1.template, false), key ->
					makeAtomPairs(res1, res2) // keep res1,res2 order here, to match inputs
				);
			}
		} else {
			// res1 and res2 are not bonded
			synchronized (atomPairsSeparate) {
				return atomPairsSeparate.computeIfAbsent(new KeySeparate(res1.template, res2.template), key ->
					makeAtomPairs(res1, res2)
				);
			}
		}
		*/

		// TEMP: don't try to cache the atom pairs, since this design is near a disulfide bond
		return makeAtomPairs(res1, res2);
	}

	private boolean isInterResBondedForward(Residue res1, Residue res2) {
		return res1.template.interResBonding.getClass() == res2.template.interResBonding.getClass()
			&& res1.template.interResBonding.isInterResBondedForward(res1, res2);
	}

	private AtomPairs makeAtomPairs(Residue res1, Residue res2) {
		
		Map<AtomNeighbors.Type,List<int[]>> pairsByType = new EnumMap<>(AtomNeighbors.Type.class);
		for (AtomNeighbors.Type type : AtomNeighbors.Type.values()) {
			pairsByType.put(type, new ArrayList<>());
		}
		
		// collect all the atom pairs by type
		for (int i=0; i<res1.atoms.size(); i++) {
			Atom atom1 = res1.atoms.get(i);
			
			AtomNeighbors neighbors;
			if (treat15HasNonBonded) {
				neighbors = new AtomNeighbors(atom1);
			} else {
				neighbors = new ProbeAtomNeighbors(atom1);
			}
			
			// for self residue pairs, skip self atom pairs and atom pairs in the other direction
			int n = i;
			if (res1 != res2) {
				n = res2.atoms.size();
			}
			
			for (int j=0; j<n; j++) {
				Atom atom2 = res2.atoms.get(j);
				
				AtomNeighbors.Type type = neighbors.classifyAtom(atom2);
				pairsByType.get(type).add(new int[] { i, j });
			}
		}

		// make the atom pairs
		AtomPairs pairs = new AtomPairs(res1.template, res2.template);
		for (Map.Entry<AtomNeighbors.Type,List<int[]>> entry : pairsByType.entrySet()) {
			AtomNeighbors.Type type = entry.getKey();
			List<int[]> atomPairs = entry.getValue();
			pairs.pairsByType[type.ordinal()] = new int[atomPairs.size()][2];
			atomPairs.toArray(pairs.pairsByType[type.ordinal()]);
		}
		
		return pairs;
	}
}
