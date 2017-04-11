package edu.duke.cs.osprey.structure;

import java.util.ArrayList;
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
import edu.duke.cs.osprey.tools.HashCalculator;

/**
 * High-performance lookups of atom pairs by bond distance
 */
public class AtomConnectivity {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace = null;
		private Molecule mol = null;
		private Parallelism parallelism = Parallelism.makeCpu(1);
	
		public Builder setConfSpace(SimpleConfSpace val) {
			confSpace = val;
			return this;
		}
		
		public Builder setMolecule(Molecule val) {
			mol = val;
			return this;
		}
		
		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}
		
		public AtomConnectivity build() {
			return new AtomConnectivity(confSpace, mol, parallelism);
		}
	}
	
	public static class AtomPair {
		
		public final int index1;
		public final int index2;
		public final AtomNeighbors.Type type;
		
		public AtomPair(int index1, int index2, AtomNeighbors.Type type) {
			this.index1 = index1;
			this.index2 = index2;
			this.type = type;
		}
		
		@Override
		public String toString() {
			return String.format("%d - %d - %s", index1, index2, type);
		}
		
		public String toString(List<Atom> atoms1, List<Atom> atoms2) {
			return String.format("%d:%s - %d:%s - %s", index1, atoms1.get(index1).name, index2, atoms2.get(index2).name, type);
		}
	}
	
	public static class AtomPairList {
		
		public final ResidueTemplate templa;
		public final ResidueTemplate templb;
		public final List<AtomPair> pairs;
		
		private int[] numAtomPairsByConnectivity;
		
		public AtomPairList(ResidueTemplate templa, ResidueTemplate templb) {
			assert (templa != null);
			assert (templb != null);
			this.templa = templa;
			this.templb = templb;
			this.pairs = new ArrayList<>();
			numAtomPairsByConnectivity = null;
		}
		
		public int size() {
			return pairs.size();
		}
		
		// residues in the cache might be in a different order and than the query
		// so these are little helper methods to swap the residue order when needed
		
		public int getIndex1(Residue res1, Residue res2, int i) {
			if (res1.template == templa && res2.template == templb) {
				return pairs.get(i).index1;
			} else if (res1.template == templb && res2.template == templa) {
				return pairs.get(i).index2;
			} else {
				throw new IllegalArgumentException(String.format("Residue templates %s, %d don't match atom pair templates: %s, %s",
					res1.template, res2.template, templa, templb
				));
			}
		}
		
		public int getIndex2(Residue res1, Residue res2, int i) {
			if (res1.template == templa && res2.template == templb) {
				return pairs.get(i).index2;
			} else if (res1.template == templb && res2.template == templa) {
				return pairs.get(i).index1;
			} else {
				throw new IllegalArgumentException(String.format("Residue templates %s, %d don't match atom pair templates: %s, %s",
					res1.template, res2.template, templa, templb
				));
			}
		}
		
		public AtomNeighbors.Type getType(int i) {
			return pairs.get(i).type;
		}
		
		private void updateCounts() {
			numAtomPairsByConnectivity = new int[pairs.size()];
			for (AtomNeighbors.Type type : AtomNeighbors.Type.values()) {
				numAtomPairsByConnectivity[type.ordinal()] = (int)pairs.stream()
					.filter((AtomPair pair) -> pair.type == type)
					.count();
			}
		}
		
		public int getNumPairs(AtomNeighbors.Type type) {
			return numAtomPairsByConnectivity[type.ordinal()];
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
		
		public Key2(ResidueTemplate templ1, ResidueTemplate templ2) {
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
			Key2 other = (Key2)obj;
			return this.templ1 == other.templ1
				&& this.templ2 == other.templ2;
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
			return HashCalculator.combineHashesCommutative(
				System.identityHashCode(templ1),
				System.identityHashCode(templ2)
			);
		}
		
		@Override
		public boolean equals(Object obj) {
			KeySeparate other = (KeySeparate)obj;
			return
				(this.templ1 == other.templ1 && this.templ2 == other.templ2)
				|| (this.templ1 == other.templ2 && this.templ2 == other.templ1);
		}
	}
	
	private TaskExecutor tasks;
	private Map<Key1,AtomPairList> atomPairs1;
	private Map<Key2,AtomPairList> atomPairs2;
	private Map<KeySeparate,AtomPairList> atomPairsSeparate;
	
	private AtomConnectivity(SimpleConfSpace confSpace, Molecule mol, Parallelism parallelism) {
		
		tasks = parallelism.makeTaskExecutor();
		
		// collect all the templates
		Set<ResidueTemplate> templatesSet = new HashSet<>();
		if (confSpace != null) {
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
					templatesSet.add(rc.template);
				}
			}
			for (Strand strand : confSpace.strands) {
				for (Residue res : strand.mol.residues) {
					templatesSet.add(res.template);
				}
			}
		}
		if (mol != null) {
			for (Residue res : mol.residues) {
				templatesSet.add(res.template);
			}
		}
		List<ResidueTemplate> templates = new ArrayList<>(templatesSet);
		
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
				(AtomPairList pairs) -> {
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
				tasks.submit(
					() -> makeDouble(templ1, templ2),
					(AtomPairList pairs) -> {
						boolean wasAdded = atomPairs2.put(new Key2(templ1, templ2), pairs) == null;
						assert (wasAdded);
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
				tasks.submit(
					() -> makeSeparate(templ1, templ2),
					(AtomPairList pairs) -> {
						boolean wasAdded = atomPairsSeparate.put(new KeySeparate(templ1, templ2), pairs) == null;
						assert (wasAdded);
					}
				);
			}
		}
		tasks.waitForFinish();
	}
	
	public AtomPairList getAtomPairs(Molecule mol, Residue res1, Residue res2) {
		
		// do we want intra pairs?
		if (res1 == res2) {
			return atomPairs1.get(new Key1(res1.template));
		}
	
		// are they bonded together?
		if (isDipeptide(res1, res2)) {
			// res1 first, then res2 adjacent
			return atomPairs2.get(new Key2(res1.template, res2.template));
		} else if (isDipeptide(res2, res1)) {
			// res2 first, then res1 adjacent
			return atomPairs2.get(new Key2(res2.template, res1.template));
		} else {
			// res1 and res2 are separated by at least one residue, or not bonded
			return atomPairsSeparate.get(new KeySeparate(res1.template, res2.template));
		}
	}
	
	private AtomPairList makeSingle(ResidueTemplate templ) {
		return makeAtomPairs(templ, templ.templateRes, templ, templ.templateRes);
	}
	
	private AtomPairList makeDouble(ResidueTemplate templ1, ResidueTemplate templ2) {
		
		Residue res1 = makeResidue(templ1);
		Residue res2 = makeResidue(templ2);
		if (!makePeptideBond(res1, res2)) {
			return null;
		}
		
		return makeAtomPairs(templ1, res1, templ2, res2);
	}
	
	private AtomPairList makeSeparate(ResidueTemplate templ1, ResidueTemplate templ2) {
		
		Residue res1 = makeResidue(templ1);
		Residue res2 = makeResidue(templ2);
		
		return makeAtomPairs(templ1, res1, templ2, res2);
	}
	
	private Residue makeResidue(ResidueTemplate templ) {
		Residue res = new Residue(Residue.copyAtoms(templ.templateRes.atoms), (double[])null, null, null);
		res.copyIntraBondsFrom(templ.templateRes);
		return res;
	}
	
	private boolean makePeptideBond(Residue res1, Residue res2) {
		
		int cIndex = res1.getAtomIndexByName("C");
		int nIndex = res2.getAtomIndexByName("N");
		
		// no way to make a peptide bond? then we don't care about this sequence of templates
		// TODO: what about non-protein chains?
		if (cIndex < 0 || nIndex < 0) {
			return false;
		}
		
		Atom C = res1.atoms.get(cIndex);
		Atom N = res2.atoms.get(nIndex);
		C.addBond(N);
		
		return true;
	}
	
	private boolean isDipeptide(Residue res1, Residue res2) {
		
		int cIndex = res1.getAtomIndexByName("C");
		int nIndex = res2.getAtomIndexByName("N");
		
		if (cIndex < 0 || nIndex < 0) {
			return false;
		}
		
		Atom C = res1.atoms.get(cIndex);
		Atom N = res2.atoms.get(nIndex);
		return C.bonds.contains(N);
	}
	
	private AtomPairList makeAtomPairs(ResidueTemplate templ1, Residue res1, ResidueTemplate templ2, Residue res2) {
		
		AtomPairList pairs = new AtomPairList(templ1, templ2);
		
		// just in case...
		assert (templ1.templateRes.atoms.size() == res1.atoms.size());
		assert (templ2.templateRes.atoms.size() == res2.atoms.size());
		
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
				
				pairs.pairs.add(new AtomPair(i, j, neighbors.classifyAtom(atom2)));
			}
		}
		
		// just in case...
		if (res1 == res2) {
			assert (pairs.pairs.size() == res1.atoms.size()*(res1.atoms.size() - 1)/2);
		} else {
			assert (pairs.pairs.size() == res1.atoms.size()*res2.atoms.size());
		}
		
		pairs.updateCounts();
		
		return pairs;
	}
}
