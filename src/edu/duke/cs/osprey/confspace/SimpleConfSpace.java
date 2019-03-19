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
import java.util.*;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.dof.*;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.Serializable;
import java.util.stream.Collectors;


/**
 * Maintains the design positions and residue conformations for a list of strands.
 * 
 * Also creates molecules in conformations on-demand.
 */
public class SimpleConfSpace implements Serializable {

	public static class Builder {
		
		private SimpleConfSpace confSpace;
		
		/**
		 * A residue is included in the steric shell if any of its atoms lies within
		 * shellDist angstroms of any atom in any flexible residue.
		 */
		private double shellDist = Double.POSITIVE_INFINITY;

		private List<Strand> strands = new ArrayList<>();
		private Map<Strand,List<StrandFlex>> strandFlex = new IdentityHashMap<>();
		
		public Builder addStrand(Strand strand, StrandFlex ... flexTypes) {
			addStrand(strand, Arrays.asList(flexTypes));
			return this;
		}
		
		public Builder addStrand(Strand strand, List<StrandFlex> flexTypes) {
			strands.add(strand);
			strandFlex.put(strand, flexTypes);
			return this;
		}
		
		public Builder addStrands(Strand ... strands) {
			for (Strand strand : strands) {
				addStrand(strand);
			}
			return this;
		}

		public Builder addStrands(Iterable<Strand> strands) {
			for (Strand strand : strands) {
				addStrand(strand);
			}
			return this;
		}
		
		public Builder setShellDistance(double val) {
			shellDist = val;
			return this;
		}
		
		public SimpleConfSpace build() {
			return new SimpleConfSpace(strands, strandFlex, shellDist);
		}
	}

	public static class NonMutablePositionException extends IllegalStateException {
		public NonMutablePositionException(SimpleConfSpace.Position pos) {
			super(String.format("Position %s was not mutable", pos));
		}
	}

	public static class Position implements Serializable {

		/** index in the positions list */
		public final int index;
		/** index in the mutable positions list, or -1 if not there */
		public final int mindex;
		public final Strand strand;
		public final String resNum;
		public final Strand.ResidueFlex resFlex;
		public final List<ResidueConf> resConfs;
		public final List<String> resTypes;

		public SeqSpace.Position seqPos = null;

		public Position(int index, int mindex, Strand strand, Residue res, List<String> resTypes) {
			this.index = index;
			this.mindex = mindex;
			this.strand = strand;
			this.resNum = res.getPDBResNumber();
			this.resFlex = strand.flexibility.get(resNum);
			this.resConfs = new ArrayList<>();
			this.resTypes = resTypes;
		}

		@Override
		public String toString() {
			return resNum;
		}

		public String formatConfPos(ConfSearch.ScoredConf conf) {
			return formatConfPos(conf.getAssignments());
		}

		public String formatConfPos(int[] conf) {
			SimpleConfSpace.ResidueConf rc = resConfs.get(conf[index]);
			return String.format("%s-%s-%d",
				rc.template.name.toUpperCase(),
				rc.getRotamerCode(),
				conf[index]
			);
		}

		public boolean hasMutations() {
			return mindex >= 0;
		}

		public int mutableIndexOrThrow() {
			if (mindex < 0) {
				throw new NonMutablePositionException(this);
			}
			return mindex;
		}
	}
	
	public static class ResidueConf implements Serializable {
		
		public static enum Type {
			
			Library('L'),
			WildType('W');
			
			public final char letter;
			
			private Type(char letter) {
				this.letter = letter;
			}
			
			public static Type getFromLetter(char letter) {
				for (Type type : values()) {
					if (type.letter == letter) {
						return type;
					}
				}
				return null;
			}
		}
		
		public static interface PostTemplateModifier extends Serializable {
			void modify(Residue res);
		}
		
		public final int index;
		public final ResidueTemplate template;
		public final Type type;
		public final Integer rotamerIndex;
		
		public PostTemplateModifier postTemplateModifier;

		public HashMap<String,double[]> dofBounds;//for DEEPer need this, else it's just a rotamer
		//for EPIC it is necessary to have some kind of persistent representation of the DOFs.  Using names here
		//so can dynamically create dihedral DOFs, etc. without getting confused

		public ResidueConf(Position pos, int index, ResidueTemplate template, Type type, HashMap<String,double[]> bbVoxel) {
			this(pos, index, template, type, null, bbVoxel);
		}
		
		public ResidueConf(Position pos, int index, ResidueTemplate template, Type type, Integer rotamerIndex, HashMap<String,double[]> bbVoxel) {
			
			if (template.getNumRotamers() > 0 && rotamerIndex == null) {
				throw new IllegalArgumentException("template " + template.name + " has rotamers, so rotamer index is required");
			}
			
			this.index = index;
			this.template = template;
			this.type = type;
			this.rotamerIndex = rotamerIndex;
			
			this.postTemplateModifier = null;

			this.dofBounds = sidechainDOFBounds(pos,template,rotamerIndex);
			this.dofBounds.putAll(bbVoxel);
		}

		public void updateResidue(ResidueTemplateLibrary templateLib, Residue res, MutAlignmentCache mutAlignmentCache) {

			// HACKHACK: make sure prolines have puckers
			boolean toProline = template.name.equalsIgnoreCase("PRO");
			if (toProline) {
				res.pucker = new ProlinePucker(templateLib, res);
			}

			ResidueTypeDOF.switchToTemplate(templateLib, res, template, false, mutAlignmentCache, false);

			// NOTE: currently only set for prolines, will idealize sidechain and check for problems
			if (postTemplateModifier != null) {
				postTemplateModifier.modify(res);
			}
		}
		
		public String getRotamerCode() {
			StringBuilder buf = new StringBuilder();
			buf.append(type.letter);
			if (rotamerIndex != null) {
				buf.append(rotamerIndex);
			}
			for(String dofName:dofBounds.keySet()) {
				if(dofName.contains("PERT")) {
					double[] bounds = dofBounds.get(dofName);
					buf.append("P"+dofBounds.get(dofName)[0]);
				}
			}
			return buf.toString();
		}
		
		@Override
		public String toString() {
			return template.name + " " + getRotamerCode();
		}
                
                public boolean isParametricallyIncompatibleWith(ResidueConf rc2){
                    //Two RCs are parametrically incompatible if and only if there is a DOF that they share
                    //for which they have different intervals
                    //Thus, a conformation has a well-defined voxel if and only if it contains no parametrically
                    //incompatible pairs of RCs

                    final double tol = 1e-8;
                    for(String dofName1 : dofBounds.keySet()){
                        for(String dofName2 : rc2.dofBounds.keySet()){
                            if(dofName1.equalsIgnoreCase(dofName2)){
                                double bounds1[] = dofBounds.get(dofName1);
                                double bounds2[] = dofBounds.get(dofName2);
                                for(int a=0; a<2; a++){//look for min or max not matching
                                    if( Math.abs( bounds1[a]-bounds2[a] ) > tol )
                                        return true;
                                }
                            }
                        }
                    }

                    //no incompatibility found!
                    return false;
                }
	}
	
	/**
	 * What kind of degree of freedom?
	 * 
	 * Used mainly to decide if we can use GPU CCD implementations or not,
	 * which only support dihedral angles for now.
	 */
	public static enum DofTypes {
		
		// NOTE: order is important for combine()
		None,
		OnlyDihedrals,
		Any;
		
		public static DofTypes combine(DofTypes a, DofTypes b) {
			if (a == b) {
				return a;
			} else if (a == null && b != null) {
				return b;
			} else if (a != null && b == null) {
				return a;
			} else {
				if (a.ordinal() > b.ordinal()) {
					return a;
				} else {
					return b;
				}
			}
		}

		public static DofTypes combine(Iterable<DofTypes> types) {
			Iterator<DofTypes> iter = types.iterator();
			if (!iter.hasNext()) {
				return null;
			}
			DofTypes out = iter.next();
			if (iter.hasNext()) {
				out = DofTypes.combine(out, iter.next());
			}
			return out;
		}
	}
	
	
	/** The strands */
	public final List<Strand> strands;

	/** The flexibility of each strand */
	public final Map<Strand,List<StrandFlex>> strandFlex;

	/** Distance from flexible residues (in Angstroms) within which to include shell residues */
	public final double shellDist;

	/** The design positions */
	public final List<Position> positions;

	/** The subset of design positions with allowed mutations */
	public final List<Position> mutablePositions;

	/** The subset of design positions with no mutations */
	public final List<Position> immutablePositions;

	/** The sequence space for this conformation space */
	public final SeqSpace seqSpace;

	/** The design positions, indexed by residue number */
	private final Map<String,Position> positionsByResNum;
	
	/** The residue numbers of the shell residues */
	public final Set<String> shellResNumbers;

	private final int[] numResConfsByPos;
	private final Molecule molTemplate;

	private final MutAlignmentCache mutAlignmentCache = new MutAlignmentCache();

	public SimpleConfSpace(List<Strand> strands, Map<Strand,List<StrandFlex>> strandFlex, double shellDist) {

		this.strands = strands;
		this.strandFlex = strandFlex;
		this.shellDist = shellDist;

		// make sure all residues are unique
		Set<String> resNums = new HashSet<>();
		for (Strand strand : strands) {
			for (Residue res : strand.mol.residues) {
				String resNum = res.getPDBResNumber();
				boolean isUnique = resNums.add(resNum);
				if (!isUnique) {
					throw new IllegalArgumentException("residue " + resNum + " appears more than once");
				}
			}
		}

		// build the design positions
		positions = new ArrayList<>();
		mutablePositions = new ArrayList<>();
		immutablePositions = new ArrayList<>();
		for (Strand strand : strands) {

			for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
				Residue res = strand.mol.getResByPDBResNumber(resNum);
				Strand.ResidueFlex resFlex = strand.flexibility.get(resNum);

				// make the pos
				int index = positions.size();
				int mindex = resFlex.isMutable() ? mutablePositions.size() : -1;
				List<String> resTypes = new ArrayList<>(resFlex.getAllResTypes());
				Position pos = new Position(index, mindex, strand, res, resTypes);
				positions.add(pos);
				if (mindex >= 0) {
					mutablePositions.add(pos);
				} else {
					immutablePositions.add(pos);
				}

				// make residue confs from library rotamers
				for (String resType : resFlex.resTypes) {
					makeResidueConfsFromTemplate(pos, strand.templateLib.getTemplateOrThrow(resType, true), ResidueConf.Type.Library);
				}

				// make residue confs from wild type rotamers
				if (resFlex.addWildTypeRotamers) {
					makeResidueConfsFromTemplate(pos, strand.templateLib.getOrMakeWildTypeTemplate(res), ResidueConf.Type.WildType);
				}
			}
		}

		// make sure we have some design positions
		if (positions.isEmpty()) {
			throw new IllegalArgumentException("ConfSpace has no design positions, try adding some strand flexibility");
		}

		// index the positions
		positionsByResNum = new HashMap<>();
		for (Position pos : positions) {
			positionsByResNum.put(Residues.normalizeResNum(pos.resNum), pos);
		}

		// collect all the static,flexible residues
		List<Residue> staticResidues = new ArrayList<>();
		List<Residue> flexibleResidues = new ArrayList<>();
		for (Strand strand : strands) {
			for (String resNum : strand.flexibility.getStaticResidueNumbers()) {
				staticResidues.add(strand.mol.getResByPDBResNumber(resNum));
			}
			for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
				flexibleResidues.add(strand.mol.getResByPDBResNumber(resNum));
			}
		}
		assert (flexibleResidues.size() == positions.size());
		assert (staticResidues.size() == resNums.size() - positions.size());

		// make the shell
		shellResNumbers = new HashSet<>();
		for (Residue staticRes : staticResidues) {

			// see if a flexible residue is nearby
			for (Residue flexibleRes : flexibleResidues) {
				if (staticRes.distanceTo(flexibleRes) <= shellDist) {
					shellResNumbers.add(Residues.normalizeResNum(staticRes.getPDBResNumber()));
					break;
				}
			}
		}
		assert (shellResNumbers.size() <= staticResidues.size());
		assert (shellDist != Double.POSITIVE_INFINITY || shellResNumbers.size() == staticResidues.size());

		// count the residue conformations
		numResConfsByPos = new int[positions.size()];
		for (int i=0; i<positions.size(); i++) {
			numResConfsByPos[i] = positions.get(i).resConfs.size();
		}

		// make the sequence space
		seqSpace = new SeqSpace(this);
		for (SeqSpace.Position seqPos : seqSpace.positions) {
			positionsByResNum.get(seqPos.resNum).seqPos = seqPos;
		}

		// make the molecule template, to speed up makeMolecule()
		molTemplate = new Molecule();
		for (Strand strand : strands) {
			for (Residue res : strand.mol.residues) {
				res = new Residue(res);
				res.molec = molTemplate;
				res.indexInMolecule = molTemplate.residues.size();
				molTemplate.residues.add(res);
			}
		}
	}

	private void makeResidueConfsFromTemplate(Position pos, ResidueTemplate template, ResidueConf.Type type) {
		
		//make one RC for each (backbone voxel, rotamer) pair
		for (HashMap<String,double[]> bbState : listBackboneVoxels(pos)) {

			if (template.name.equalsIgnoreCase("PRO") && type == ResidueConf.Type.Library) {

				// HACKHACK: add one conf for each proline pucker
				for (ProlinePucker.Direction dir : ProlinePucker.Direction.values()) {
					ResidueConf resConf = new ResidueConf(
						pos,
						pos.resConfs.size(),
						template,
						type,
						dir.ordinal(),
						bbState
					);
					resConf.postTemplateModifier = (res) -> res.pucker.apply(dir);
					pos.resConfs.add(resConf);
				}

			} else if (template.getNumRotamers() <= 0) {

				// make one template for the library template
				pos.resConfs.add(new ResidueConf(
					pos,
					pos.resConfs.size(),
					template,
					type,
					bbState
				));

			} else {

				// make a template for each rotamer
				for (int rotamerIndex=0; rotamerIndex<template.getNumRotamers(); rotamerIndex++) {
					pos.resConfs.add(new ResidueConf(
						pos,
						pos.resConfs.size(),
						template,
						type,
						rotamerIndex,
						bbState
					));
				}
			}
		}
	}

	/** Gets a design position by residue number, or returns null */
	public Position getPositionOrNull(String resNum) {
		return positionsByResNum.get(Residues.normalizeResNum(resNum));
	}

	/** Gets a design position by residue number, or throws an exception */
	public Position getPositionOrThrow(String resNum) {
		SimpleConfSpace.Position pos = getPositionOrNull(resNum);
		if (pos != null) {
			return pos;
		}
		throw new NoSuchElementException(
			"no residue found with number " + resNum + " among " + positions.stream()
				.map((p) -> p.resNum)
				.collect(Collectors.toList())
		);
	}

	/** Gets the number of residue confs per position */
	public int[] getNumResConfsByPos() {
		return numResConfsByPos;
	}
	
	/** Gets the number of residue confs for a position */
	public int getNumResConfs(int pos) {
		return positions.get(pos).resConfs.size();
	}
	
	/** Gets the total number of residue confs for all positions */
	public int getNumResConfs() {
		int count = 0;
		for (int pos=0; pos<positions.size(); pos++) {
			count += numResConfsByPos[pos];
		}
		return count;
	}
	
	/** Gets the total number of residue conf pairs for all positions */
	public int getNumResConfPairs() {
		int count = 0;
		for (int pos1=1; pos1<positions.size(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				count += numResConfsByPos[pos1]*numResConfsByPos[pos2];
			}
		}
		return count;
	}

	/** Gets the total number of residue conf triples for all positions */
	public int getNumResConfTriples() {
		int count = 0;
		for (int pos1=2; pos1<positions.size(); pos1++) {
			for (int pos2=1; pos2<pos1; pos2++) {
				for (int pos3=0; pos3<pos2; pos3++) {
					count += numResConfsByPos[pos1]*numResConfsByPos[pos2]*numResConfsByPos[pos3];
				}
			}
		}
		return count;
	}

	/** Gets the total number of residue conf quads for all positions */
	public int getNumResConfQuads() {
		int count = 0;
		for (int pos1=2; pos1<positions.size(); pos1++) {
			for (int pos2=1; pos2<pos1; pos2++) {
				for (int pos3=0; pos3<pos2; pos3++) {
					for (int pos4=0; pos4<pos3; pos4++) {
						count += numResConfsByPos[pos1]*numResConfsByPos[pos2]*numResConfsByPos[pos3]*numResConfsByPos[pos4];
					}
				}
			}
		}
		return count;
	}

	/** @see #makeMolecule(RCTuple) */
	public ParametricMolecule makeMolecule(ScoredConf conf) {
		return makeMolecule(conf.getAssignments());
	}

	/** @see #makeMolecule(RCTuple) */
	public ParametricMolecule makeMolecule(int[] conf) {
		return makeMolecule(new RCTuple(conf));
	}

	/**
	 * create a new {@link Molecule} in the specified conformation
	 * but without any continuous flexibility
	 */
	public Molecule makeDiscreteMolecule(RCTuple conf) {

		Molecule mol = new Molecule();
		for (Residue res : molTemplate.residues) {

			// is this a conformation residue?
			Position pos = getPositionOrNull(res.getPDBResNumber());
			if (pos != null) {

				// is this residue part of the given conformation?
				int index = conf.pos.indexOf(pos.index);
				if (index >= 0) {

					// yup, get the RC
					ResidueConf rc = pos.resConfs.get(conf.RCs.get(index));

					// build the residue from the RC
					Residue newRes = res.copyToMol(mol, false);
					rc.updateResidue(pos.strand.templateLib, newRes, mutAlignmentCache);

					continue;
				}
			}

			// otherise, just copy from the template molecule
			res.copyToMol(mol, true);
		}

		mol.markInterResBonds();

		return mol;
	}
	
	/**
	 * create a new {@link ParametricMolecule} in the specified conformation
	 * for analysis (e.g., energy calculation, minimization)
	 * 
	 * To increase stability of analysis, each analysis should be conducted
	 * with a new molecule instance. this completely prevents roundoff error
	 * from accumulating across separate analyses. 
	 */
	public ParametricMolecule makeMolecule(RCTuple conf) {

		// copy the molecule from the template and make mutations as needed
		Molecule mol = makeDiscreteMolecule(conf);

		// figure out what conformational DOFs are specified by the conf
		HashSet<String> confDOFNames = new HashSet<>();
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			confDOFNames.addAll(resConf.dofBounds.keySet());
		}

		// OK now apply all DOF vals including puckers

		// make all the DOFs
		List<DegreeOfFreedom> dofs = new ArrayList<>();
		
		// first, backbone flexibility: strand DOFs and DEEPer/CATS DOFs
		for (Strand strand : getConfStrands(conf)) {
			for (StrandFlex flex : strandFlex.get(strand)) {
	 			for (DegreeOfFreedom dof : flex.makeDofs(strand, mol)) {
					if (confDOFNames.contains(dof.getName())) {
						//DEEPer and CATS DOFS may not involve all positions in a strand
						dofs.add(dof);
					}
				}
			}
		}

		// then, residue conf DOFs
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			Residue res = mol.getResByPDBResNumber(pos.resNum);
			
			// make the residue DOFs
			Strand.ResidueFlex resFlex = pos.strand.flexibility.get(pos.resNum);
			List<DegreeOfFreedom> contDihedralDOFs = resFlex.voxelShape.makeDihedralDOFs(res);
			dofs.addAll(contDihedralDOFs);

			//so we apply them here instead
			if (contDihedralDOFs.isEmpty()) {
				// pose the residue to match the rotamer
				List<DegreeOfFreedom> dihedralDofs = new VoxelShape.Rect().makeDihedralDOFs(res);
				for (int d=0; d<resConf.template.numDihedrals; d++) {
					dihedralDofs.get(d).apply(resConf.template.getRotamericDihedrals(resConf.rotamerIndex, d));
				}
			}
		}

		//Figure out the bounds and apply the DOF values in the middle of the voxel
		DofBounds dofBounds = new DofBounds(dofs.size());//the bounds to use
		HashMap<String,Integer> name2Index = DegreeOfFreedom.nameToIndexMap(dofs);//map DOF names to their indices in DOFs
		HashSet<String> dofsAdded = new HashSet<>();

		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			HashMap<String,double[]> rcDOFBounds = resConf.dofBounds;
			for(String DOFName : rcDOFBounds.keySet()) {
				int dofIndex = name2Index.get(DOFName);
				double curDOFBounds[] = rcDOFBounds.get(DOFName);
				if (dofsAdded.contains(DOFName)) {//we have bounds on this DOF from another position...check consistency
					if (Math.abs(dofBounds.getMax(dofIndex)-curDOFBounds[1])>1e-10
							|| Math.abs(dofBounds.getMin(dofIndex)-curDOFBounds[0])>1e-10) {
						throw new RuntimeException("ERROR: Conformation has inconsistent DOF bounds"
								+ " between different positions' RCs");
					}
				} else { // record the bounds and apply the middle value
					dofBounds.set(dofIndex, curDOFBounds[0], curDOFBounds[1]);
					dofs.get(dofIndex).apply(0.5*(curDOFBounds[0]+curDOFBounds[1]));
					dofsAdded.add(DOFName);
				}
			}
		}

		return new ParametricMolecule(mol, dofs, dofBounds);
	}
	
	/** @see #isContinuouslyFlexible(RCTuple) */
	public boolean isContinuouslyFlexible(int[] conf) {
		return isContinuouslyFlexible(new RCTuple(conf));
	}
	
	/** Return true if the conformation has continuous degrees of freedom, false otherwise. */
	public boolean isContinuouslyFlexible(RCTuple conf) {
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			for (double[] bounds : resConf.dofBounds.values()) {
				if (bounds[1]>bounds[0]) {
					return true;
				}
			}
		}
		return false;
	}

	private Set<Strand> getConfStrands(RCTuple conf) {
		Set<Strand> confStrands = new HashSet<>();
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			confStrands.add(pos.strand);
		}
		return confStrands;
	}

	public DofTypes getDofTypes() {
		
		DofTypes dofTypes = null;
		
		// check strand flex
		for (Strand strand : strands) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				dofTypes = DofTypes.combine(dofTypes, flex.getDofTypes());
			}
		}
		
		// check residue flex
		for (Position pos : positions) {
			dofTypes = DofTypes.combine(dofTypes, pos.resFlex.voxelShape.getDofTypes());
		}
		
		return dofTypes;
	}

	public DofInfo makeDofInfo(RCTuple conf) {

		// NOTE: this method should describe dofs in the same order that makeMolecule() does

		DofInfo info = new DofInfo(conf);

		// backbone flexibility
		for (Strand strand : getConfStrands(conf)) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				info.addStrand(strand, flex.countDofs(strand));
			}
		}

		// residue flexibility
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf rc = pos.resConfs.get(conf.RCs.get(i));
			Strand.ResidueFlex resFlex = pos.strand.flexibility.get(pos.resNum);
			info.addPos(pos, rc, resFlex.voxelShape.countDihedralDOFs(rc.template));
		}

		return info;
	}
        
        public int getNumPos(){
            return positions.size();
        }
        
        private static HashMap<String,double[]> sidechainDOFBounds(Position pos, ResidueTemplate template, Integer rotamerIndex){
            //get bounds on the sidechain dihedrals associated with a particular rotamer
            if(rotamerIndex==null)//no dihedrals
                return new HashMap<>();
            
            DofBounds bounds = pos.resFlex.voxelShape.makeDihedralBounds(template, rotamerIndex);
            Residue res = pos.strand.mol.getResByPDBResNumber(pos.resNum);

            HashMap<String,double[]> ans = new HashMap<>();
            for(int d=0; d<bounds.size(); d++){//in rigid case will not make bounds
                FreeDihedral dih = new FreeDihedral(res,d);
                ans.put(dih.getName(), new double[]{bounds.getMin(d),bounds.getMax(d)});
            }
            
            return ans;
        }
        
        
        private ArrayList<HashMap<String,double[]>> listBackboneVoxels(Position pos){
            ArrayList<HashMap<String,double[]> > ans = new ArrayList<>();
            
            for(StrandFlex flexType : strandFlex.get(pos.strand)){
                ArrayList<HashMap<String,double[]> > flexTypeVox = flexType.listBackboneVoxels(pos);
                if( (!flexTypeVox.isEmpty()) && (!ans.isEmpty()) ){
                    throw new RuntimeException("ERROR: Can't have multiple types of backbone flexibility for the same residue");
                    //not supported, because current backbone DOF implementations depend on the DOF block
                    //keeping track of the unperturbed backbone conformation, so if another type of motion
                    //changes that unperturbed bb conf, then there will be errors
                    //Mutations and sidechain dihedrals don't move the backbone so no issue there
                }
                else
                    ans = flexTypeVox;
            }
            
            if(ans.isEmpty())//no backbone flexibility
                return new ArrayList(Arrays.asList(new HashMap<>()));
            else
                return ans;
        }

	public SimpleConfSpace makeSubspace(Strand strand) {
		return new Builder()
			.addStrand(strand, strandFlex.get(strand))
			.setShellDistance(shellDist)
			.build();
	}

	public Map<Position,Position> mapPositionsTo(SimpleConfSpace other) {

		Map<Position,Position> matching = new IdentityHashMap<>();

		for (Position pos : positions) {

			// try to match to the other pos by res num
			Position matchedPos = other.getPositionOrNull(pos.resNum);

			// just in case...
			if (matchedPos != null) {
				if (pos.resConfs.size() != matchedPos.resConfs.size()) {
					throw new IllegalStateException("residue conformations for residue " + pos.resNum + " don't match across conformation spaces");
				}
				assert (pos.resNum.equals(matchedPos.resNum));
			}

			matching.put(pos, matchedPos);
		}

		return matching;
	}

	public List<ResidueConf> getResidueConfs(ConfSearch.ScoredConf conf) {
		return getResidueConfs(conf.getAssignments());
	}

	public List<ResidueConf> getResidueConfs(int[] conf) {
		return positions.stream()
			.map((pos) -> pos.resConfs.get(conf[pos.index]))
			.collect(Collectors.toList());
	}

	public String formatConf(ConfSearch.ScoredConf conf) {
		return formatConf(conf.getAssignments());
	}

	public String formatConf(int[] conf) {
		return String.join(" ", positions.stream()
			.map((pos) -> String.format("%-12s", pos.formatConfPos(conf)))
			.collect(Collectors.toList())
		);
	}

	public static String formatConfRCs(ConfSearch.ScoredConf conf) {
		return formatConfRCs(conf.getAssignments());
	}

	public static String formatConfRCs(int[] conf) {
		return String.join(" ", Arrays.stream(conf)
			.mapToObj((int rc) -> String.format("%-5d", rc))
			.collect(Collectors.toList())
		);
	}

	public String formatConfSequence(ConfSearch.ScoredConf conf) {
		return formatConfSequence(conf.getAssignments());
	}

	public String formatConfSequence(int[] conf) {
		return String.join(" ", getResidueConfs(conf).stream()
			.map((resConf) -> String.format("%-5s", resConf.template.name))
			.collect(Collectors.toList())
		);
	}

	public String formatConfRotamers(ConfSearch.ScoredConf conf) {
		return formatConfRotamers(conf.getAssignments());
	}

	public String formatConfRotamers(int[] conf) {
		return String.join(" ", getResidueConfs(conf).stream()
			.map((resConf) -> String.format("%-5s", resConf.getRotamerCode()))
			.collect(Collectors.toList())
		);
	}

	public String formatConfRotamersWithResidueNumbers(int[] conf) {
		return String.join(",",positions.stream()
				.map((pos) -> {
					String rcString = pos.resNum+":*";
				    if(conf[pos.index] > -1) {
				    	ResidueConf rc = pos.resConfs.get(conf[pos.index]);
						rcString = pos.resNum+":"+rc.template.name+"-"+rc.getRotamerCode();
					}
					return rcString;
				})
				.collect(Collectors.toList()));
	}

	public String formatResidueNumbers() {
		return String.join(" ", positions.stream()
			.map((pos) -> String.format("%-5s", pos.resNum))
			.collect(Collectors.toList())
		);
	}

	public Sequence makeUnassignedSequence() {
		return seqSpace.makeUnassignedSequence();
	}

	public Sequence makeWildTypeSequence() {
		return seqSpace.makeWildTypeSequence();
	}

	public Sequence makeSequenceFromAssignments(int[] assignments) {
		Sequence seq = seqSpace.makeUnassignedSequence();
		for (Position pos : positions) {
			seq.set(pos.resNum, pos.resConfs.get(assignments[pos.index]).template.name);
		}
		return seq;
	}

	public Sequence makeSequenceFromConf(ConfSearch.ScoredConf conf) {
		return makeSequenceFromAssignments(conf.getAssignments());
	}

	public BigInteger calcNumSequences() {
		BigInteger numSequences = BigInteger.ZERO;
		for (SimpleConfSpace.Position pos : positions) {
			BigInteger numResTypes = BigInteger.valueOf(pos.resFlex.resTypes.size());
			if (MathTools.isZero(numSequences)) {
				numSequences = numResTypes;
			} else {
				numSequences = numSequences.multiply(numResTypes);
			}
		}
		return numSequences;
	}
}
