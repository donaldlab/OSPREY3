package edu.duke.cs.osprey.structure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.HashCalculator;

/**
 * High-performance lookups of atom pairs by bond distance
 */
public class AtomConnectivity {
	
	public static Set<ResidueTemplate> collectTemplates(SimpleConfSpace confSpace) {
		Set<ResidueTemplate> templates = new HashSet<>();
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
				templates.add(rc.template);
			}
		}
		for (Strand strand : confSpace.strands) {
			for (Residue res : strand.mol.residues) {
				templates.add(res.template);
			}
		}
		return templates;
	}
	
	public static Set<ResidueTemplate> collectTemplates(Residues residues) {
		Set<ResidueTemplate> templates = new HashSet<>();
		for (Residue res : residues) {
			templates.add(res.template);
		}
		return templates;
	}
	
	public static class Builder {
		
		private Set<ResidueTemplate> templates = new HashSet<>();
		private Parallelism parallelism = Parallelism.makeCpu(1);

		public Builder addTemplates(Collection<ResidueTemplate> val) {
			templates.addAll(val);
			return this;
		}

		public Builder addTemplates(ResidueTemplateLibrary val) {
			addTemplates(val.templates);
			return this;
		}

		public Builder addTemplates(SimpleConfSpace confSpace) {
			addTemplates(collectTemplates(confSpace));
			return this;
		}
		
		public Builder addTemplates(Residues residues) {
			addTemplates(collectTemplates(residues));
			return this;
		}
		
		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}
		
		public AtomConnectivity build() {
			return new AtomConnectivity(new ArrayList<>(templates), parallelism);
		}
	}
	
	public static class AtomPairs {
		
		public final Residue res1;
		public final Residue res2;
		
		public AtomPairs(Residue res1, Residue res2) {
			this.res1 = res1;
			this.res2 = res2;
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
	
	private Map<Key1,AtomPairs> atomPairs1;
	private Map<Key2,AtomPairs> atomPairs2;
	private Map<KeySeparate,AtomPairs> atomPairsSeparate;
	
	private AtomConnectivity(List<ResidueTemplate> templates, Parallelism parallelism) {
		
		// make sure we have residue templates
		if (templates == null || templates.isEmpty()) {
			throw new IllegalArgumentException("templates cannot be empty. Try adding templates from a molecule or a conf space");
		}
		
		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
		
			/* DEBUG: show template info
			for (ResidueTemplate template : templates) {
				if (template.name.equals("ALA") || template.name.equals("GLU")) {
					System.out.println(template);
					for (int i=0; i<template.templateRes.atoms.size(); i++) {
						System.out.println("\t" + template.templateRes.atoms.get(i).name);
					}
				}
			}
			*/
			
			// do singles
			atomPairs1 = new HashMap<>();
			for (int i=0; i<templates.size(); i++) {
				ResidueTemplate templ1  = templates.get(i);
				tasks.submit(
					() -> makeSingle(templ1),
					(AtomPairs pairs) -> {
						boolean wasAdded = atomPairs1.put(new Key1(templ1), pairs) == null;
						assert (wasAdded);
					}
				);
			}
			tasks.waitForFinish();
			
			// do doubles
			atomPairs2 = new HashMap<>();
			for (int i=0; i<templates.size(); i++) {
				ResidueTemplate templ1  = templates.get(i);
				for (int j=0; j<templates.size(); j++) {
					ResidueTemplate templ2  = templates.get(j);
					
					// make the forward order
					tasks.submit(
						() -> makeDouble(templ1, templ2),
						(AtomPairs pairs) -> {
							boolean wasAdded = atomPairs2.put(new Key2(templ1, templ2, true), pairs) == null;
							assert (wasAdded);
							
							// make the reverse order
							AtomPairs swappedPairs = makeSwappedPairs(pairs);
							boolean wasSwappedAdded = atomPairs2.put(new Key2(templ1, templ2, false), swappedPairs) == null;
							assert (wasSwappedAdded);
						}
					);
				}
			}
			tasks.waitForFinish();
			
			// do separates
			atomPairsSeparate = new HashMap<>();
			for (int i=0; i<templates.size(); i++) {
				ResidueTemplate templ1  = templates.get(i);
				for (int j=0; j<=i; j++) {
					ResidueTemplate templ2  = templates.get(j);
					
					// make the forward order
					tasks.submit(
						() -> makeSeparate(templ1, templ2),
						(AtomPairs pairs) -> {
							boolean wasAdded = atomPairsSeparate.put(new KeySeparate(templ1, templ2), pairs) == null;
							assert (wasAdded);
							
							// make the reverse order if needed
							if (templ1 != templ2) {
								AtomPairs swappedPairs = makeSwappedPairs(pairs);
								boolean wasSwappedAdded = atomPairsSeparate.put(new KeySeparate(templ2, templ1), swappedPairs) == null;
								assert (wasSwappedAdded);
							}
						}
					);
				}
			}
			tasks.waitForFinish();
		}
	}
	
	public AtomPairs getAtomPairs(Residue res1, Residue res2) {
		
		// do we want intra pairs?
		if (res1 == res2) {
			return atomPairs1.get(new Key1(res1.template));
		}
	
		// are they bonded together?
		if (isDipeptide(res1, res2)) {
			// yup, in forward order
			return atomPairs2.get(new Key2(res1.template, res2.template, true));
		} else if (isDipeptide(res2, res1)) {
			// yup, in reverse order
			return atomPairs2.get(new Key2(res2.template, res1.template, false));
		} else {
			// res1 and res2 are not bonded
			return atomPairsSeparate.get(new KeySeparate(res1.template, res2.template));
		}
	}
	
	private AtomPairs makeSingle(ResidueTemplate templ) {
		Residue res = makeResidue(templ);
		return makeAtomPairs(res, res);
	}
	
	private AtomPairs makeDouble(ResidueTemplate templ1, ResidueTemplate templ2) {
		
		Residue res1 = makeResidue(templ1);
		Residue res2 = makeResidue(templ2);
		if (!makePeptideBond(res1, res2)) {
			return null;
		}
		
		return makeAtomPairs(res1, res2);
	}
	
	private AtomPairs makeSeparate(ResidueTemplate templ1, ResidueTemplate templ2) {
		Residue res1 = makeResidue(templ1);
		Residue res2 = makeResidue(templ2);
		return makeAtomPairs(res1, res2);
	}
	
	private Residue makeResidue(ResidueTemplate templ) {
		Residue res = new Residue(Residue.copyAtoms(templ.templateRes.atoms), (double[])null, null, null);
		res.copyIntraBondsFrom(templ.templateRes);
		res.template = templ;
		return res;
	}
	
	private boolean makePeptideBond(Residue res1, Residue res2) {
		
		Atom C = res1.getAtomByName("C");
		Atom N = res2.getAtomByName("N");
		
		// no way to make a peptide bond? then we don't care about this sequence of templates
		// TODO: what about non-protein chains?
		if (C == null || N == null) {
			return false;
		}
		
		C.addBond(N);
		
		return true;
	}
	
	private boolean isDipeptide(Residue res1, Residue res2) {
		
		Atom C = res1.getAtomByName("C");
		Atom N = res2.getAtomByName("N");
		
		if (C == null || N == null) {
			return false;
		}
		
		return C.bonds.contains(N);
	}
	
	private AtomPairs makeAtomPairs(Residue res1, Residue res2) {
		
		Map<AtomNeighbors.Type,List<int[]>> pairsByType = new EnumMap<>(AtomNeighbors.Type.class);
		for (AtomNeighbors.Type type : AtomNeighbors.Type.values()) {
			pairsByType.put(type, new ArrayList<>());
		}
		
		// collect all the atom pairs by type
		for (int i=0; i<res1.atoms.size(); i++) {
			Atom atom1 = res1.atoms.get(i);
			
			AtomNeighbors neighbors = new AtomNeighbors(atom1);
			
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
		AtomPairs pairs = new AtomPairs(res1, res2);
		for (Map.Entry<AtomNeighbors.Type,List<int[]>> entry : pairsByType.entrySet()) {
			AtomNeighbors.Type type = entry.getKey();
			List<int[]> atomPairs = entry.getValue();
			pairs.pairsByType[type.ordinal()] = new int[atomPairs.size()][2];
			atomPairs.toArray(pairs.pairsByType[type.ordinal()]);
		}
		
		return pairs;
	}
	
	private AtomPairs makeSwappedPairs(AtomPairs pairs) {

		if (pairs == null) {
			return null;
		}
		
		AtomPairs swapped = new AtomPairs(pairs.res2, pairs.res1);
		
		for (AtomNeighbors.Type type : AtomNeighbors.Type.values()) {
			
			int n = pairs.getNumPairs(type);
			
			int[][] swappedTypedPairs = new int[n][2];
			int[][] typedPairs = pairs.getPairs(type);
			
			for (int i=0; i<n; i++) {
				swappedTypedPairs[i][0] = typedPairs[i][1];
				swappedTypedPairs[i][1] = typedPairs[i][0];
			}
			
			swapped.pairsByType[type.ordinal()] = swappedTypedPairs;
		}
		
		return swapped;
	}
}
