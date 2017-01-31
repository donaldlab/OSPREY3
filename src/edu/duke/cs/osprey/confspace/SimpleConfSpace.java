package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

/**
 * maintains the design positions and residue conformations for a list of strands
 * 
 * also creates molecules in conformations on-demand
 */
public class SimpleConfSpace {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace;
		private double shellDist;
		
		public Builder() {
			confSpace = new SimpleConfSpace();
			shellDist = Double.POSITIVE_INFINITY;
		}
		
		public Builder addStrand(Strand strand, StrandFlex ... flexType) {
			confSpace.addStrand(strand, flexType);
			return this;
		}
		
		public Builder setShellDistance(double val) {
			shellDist = val;
			return this;
		}
		
		public SimpleConfSpace build() {
			
			// make sure we have some design positions
			if (confSpace.positions.isEmpty()) {
				throw new IllegalStateException("ConfSpace has no design positions, try adding some strand flexibility");
			}
			
			confSpace.makeShell(shellDist);
			confSpace.countResConfs();
			return confSpace;
		}
	}
	
	public static Builder builder() {
		return new Builder();
	}
	
	public static SimpleConfSpace build(Strand strand) {
		return builder().addStrand(strand).build();
	}
	
	public static class Position {
		
		public final int index;
		public final Strand strand;
		public final String resNum;
		public final Strand.ResidueFlex resFlex;
		public final List<ResidueConf> resConfs;
		
		public Position(int index, Strand strand, Residue res) {
			this.index = index;
			this.strand = strand;
			this.resNum = res.getPDBResNumber();
			this.resFlex = strand.flexibility.get(resNum);
			this.resConfs = new ArrayList<>();
		}
	}
	
	public static class ResidueConf {
		
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
		
		public static interface PostTemplateModifier {
			void modify(Residue res);
		}
		
		public final int index;
		public final ResidueTemplate template;
		public final Type type;
		public final Integer rotamerIndex;
		
		public PostTemplateModifier postTemplateModifier;
		
		public ResidueConf(int index, ResidueTemplate template, Type type) {
			this(index, template, type, null);
		}
		
		public ResidueConf(int index, ResidueTemplate template, Type type, Integer rotamerIndex) {
			
			if (template.getNumRotamers() > 0 && rotamerIndex == null) {
				throw new IllegalArgumentException("template " + template.name + " has rotamers, so rotamer index is required");
			}
			
			this.index = index;
			this.template = template;
			this.type = type;
			this.rotamerIndex = rotamerIndex;
			
			this.postTemplateModifier = null;
		}

		public void updateResidue(GenericResidueTemplateLibrary templateLib, Residue res) {
			
			// HACKHACK: make sure prolines have puckers
			if (res.template.name.equalsIgnoreCase("PRO") || template.name.equalsIgnoreCase("PRO")) {
				res.pucker = new ProlinePucker(templateLib, res);
			}
			
			if (postTemplateModifier != null) {
				postTemplateModifier.modify(res);
			}
			
			ResidueTypeDOF.switchToTemplate(templateLib, res, template, false);
		}
		
		public String getRotamerCode() {
			StringBuilder buf = new StringBuilder();
			buf.append(type.letter);
			if (rotamerIndex != null) {
				buf.append(rotamerIndex);
			}
			return buf.toString();
		}
		
		@Override
		public String toString() {
			return template.name + " " + getRotamerCode();
		}
	}
	
	public final List<Strand> strands;
	public final List<Position> positions;
	public final Set<String> shellResNumbers;
	public final Map<Strand,List<StrandFlex>> strandFlex; // yeah, map on instance, not identity
	
	private Map<Strand,Set<String>> shellResNumbersByStrand;
	private int[] numResConfsByPos;
	
	private SimpleConfSpace() {
		strands = new ArrayList<>();
		positions = new ArrayList<>();
		shellResNumbers = new HashSet<>();
		strandFlex = new HashMap<>();
		
		shellResNumbersByStrand = new HashMap<>();
		numResConfsByPos = null;
	}
	
	private void addStrand(Strand strand, StrandFlex ... flexType) {
		
		// add the strand
		strands.add(strand);
		strandFlex.put(strand, Arrays.asList(flexType));
		
		// make the positions
		for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
			Residue res = strand.mol.getResByPDBResNumber(resNum);
			Strand.ResidueFlex resFlex = strand.flexibility.get(resNum);
			
			// make the pos
			Position pos = new Position(positions.size(), strand, res);
			positions.add(pos);
			
			// make residue confs from library rotamers
			for (String resType : resFlex.resTypes) {
				makeResidueConfsFromTemplate(pos, strand.templateLib.getTemplate(resType), ResidueConf.Type.Library);
			}
			
			// make residue confs from wild type rotamers
			if (resFlex.addWildTypeRotamers) {
				
				List<Residue> residues = new ArrayList<>();
				residues.add(res);
				residues.addAll(strand.mol.getAlternates(res.indexInMolecule));
				
				makeResidueConfsFromTemplate(pos, ResidueTemplate.makeFromResidueConfs(residues), ResidueConf.Type.WildType);
			}
		}
		
		// TODO: make positions for non-residue discrete confs (eg blocks, DEEPer)
	}
	
	private void makeResidueConfsFromTemplate(Position pos, ResidueTemplate template, ResidueConf.Type type) {
		
		if (template.name.equalsIgnoreCase("PRO")) {
			
			// HACKHACK: add one cone for each proline pucker
			for (ProlinePucker.Direction dir : ProlinePucker.Direction.values()) {
				ResidueConf resConf = new ResidueConf(
					pos.resConfs.size(),
					template,
					type
				);
				resConf.postTemplateModifier = (res) -> res.pucker.apply(dir);
				pos.resConfs.add(resConf);
			}
			
		} else if (template.getNumRotamers() <= 0) {
			
			// make one template for the library template
			pos.resConfs.add(new ResidueConf(
				pos.resConfs.size(),
				template,
				type
			));
			
		} else {
			
			// make a template for each rotamer
			for (int rotamerIndex=0; rotamerIndex<template.getNumRotamers(); rotamerIndex++) {
				pos.resConfs.add(new ResidueConf(
					pos.resConfs.size(),
					template,
					type,
					rotamerIndex
				));
			}
		}
	}
	
	private void makeShell(double shellDist) {
		
		shellResNumbers.clear();
		shellResNumbersByStrand.clear();
	
		for (Strand strand : strands) {
			
			Set<String> resNumbers = new HashSet<>();
			
			for (String staticResNum : strand.flexibility.getStaticResidueNumbers()) {
				Residue staticRes = strand.mol.getResByPDBResNumber(staticResNum);
				
				for (String flexResNum : strand.flexibility.getFlexibleResidueNumbers()) {
					Residue flexRes = strand.mol.getResByPDBResNumber(flexResNum);
					
					if (staticRes.distanceTo(flexRes) <= shellDist) {
						resNumbers.add(staticResNum);
						break;
					}
				}
			}
			
			shellResNumbers.addAll(resNumbers);
			shellResNumbersByStrand.put(strand, resNumbers);
		}
	}
	
	private void countResConfs() {
		numResConfsByPos = new int[positions.size()];
		for (int i=0; i<positions.size(); i++) {
			numResConfsByPos[i] = positions.get(i).resConfs.size();
		}
	}
	
	public int[] getNumResConfsByPos() {
		return numResConfsByPos;
	}
	
	public int getNumResConfs(int pos) {
		return positions.get(pos).resConfs.size();
	}
	
	public int getNumResConfs() {
		int count = 0;
		for (int pos=0; pos<positions.size(); pos++) {
			count += numResConfsByPos[pos];
		}
		return count;
	}
	
	public int getNumResConfPairs() {
		int count = 0;
		for (int pos1=0; pos1<positions.size(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				for (int rc1=0; rc1<numResConfsByPos[pos2]; rc1++) {
					count += numResConfsByPos[pos2];
				}
			}
		}
		return count;
	}
	
	public ParametricMolecule makeMolecule(int[] conf) {
		return makeMolecule(new RCTuple(conf));
	}
	
	/**
	 * create a new {@link ParametricMolecule} in the specified conformation
	 * for analysis (e.g., minimization)
	 * 
	 * To increase stability of analysis, each analysis should be conducted
	 * with a new molecule instance. this completely prevents roundoff error
	 * from accumulating across separate analyses. 
	 */
	public ParametricMolecule makeMolecule(RCTuple conf) {
		
		// make the molecule from the strands (ignore alternates)
		Molecule mol = new Molecule();
		for (Strand strand : strands) {
			for (Residue res : strand.mol.residues) {
				res = new Residue(res);
				res.molec = mol;
				mol.residues.add(res);
			}
		}
		HardCodedResidueInfo.markInterResBonds(mol);
		
		// mutate to the conf templates
		for (int i=0; i<conf.size(); i++) {
			
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			Residue res = mol.getResByPDBResNumber(pos.resNum);
			
			resConf.updateResidue(pos.strand.templateLib, res);
		}
		
		// make all the DOFs
		List<DegreeOfFreedom> dofs = new ArrayList<>();
		
		// first, strand DOFs
		for (Strand strand : getConfStrands(conf)) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				dofs.addAll(flex.makeDofs(strand));
			}
			// NOTE: strand DOFs are "centered" on initial molecule pos, so there's no need to set anything here to match the conf
		}
		
		// then, residue conf DOFs
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			Residue res = mol.getResByPDBResNumber(pos.resNum);
			
			// make the residue DOFs
			Strand.ResidueFlex resFlex = pos.strand.flexibility.get(pos.resNum);
			dofs.addAll(resFlex.voxelShape.makeDihedralDOFs(res));
			
			// pose the residue to match the rotamer
			List<DegreeOfFreedom> dihedralDofs = new VoxelShape.Rect().makeDihedralDOFs(res);
			for (int d=0; d<resConf.template.numDihedrals; d++) {
				dihedralDofs.get(d).apply(resConf.template.getRotamericDihedrals(resConf.rotamerIndex, d));
			}
		}
		
		// TODO: make sure prolines have puckers
		// TODO: make DOFs for non-residue continuous motions (eg, blocks, DEEPer)
		
		return new ParametricMolecule(mol, dofs);
	}
	
	public boolean isContinuouslyFlexible(int[] conf) {
		return isContinuouslyFlexible(new RCTuple(conf));
	}
	
	public boolean isContinuouslyFlexible(RCTuple conf) {
		return makeBounds(conf).size() > 0;
	}
	
	public DofBounds makeBounds(int[] conf) {
		return makeBounds(new RCTuple(conf));
	}
	
	public DofBounds makeBounds(RCTuple conf) {
		
		// gather all the DOF bounds
		List<DofBounds> bounds = new ArrayList<>();
		
		// first, strand bounds
		for (Strand strand : getConfStrands(conf)) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				bounds.add(flex.makeBounds(strand));
			}
		}
		
		// then, residue conf bounds
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			ResidueConf resConf = pos.resConfs.get(conf.RCs.get(i));
			if (resConf.rotamerIndex != null) {
				bounds.add(pos.resFlex.voxelShape.makeDihedralBounds(resConf.template, resConf.rotamerIndex));
			}
		}
		
		// TODO: make DOF bounds for non-residue continuous motions (eg, blocks, DEEPer)
		
		return DofBounds.concatenate(bounds);
	}
	
	private Set<Strand> getConfStrands(RCTuple conf) {
		Set<Strand> confStrands = new HashSet<>();
		for (int i=0; i<conf.size(); i++) {
			Position pos = positions.get(conf.pos.get(i));
			confStrands.add(pos.strand);
		}
		return confStrands;
	}

	public boolean isGpuCcdSupported() {
		
		// check strand flex
		for (Strand strand : strands) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				if (!flex.isGpuCcdSupported()) {
					return false;
				}
			}
		}
		
		// check residue flex
		for (Position pos : positions) {
			if (!pos.resFlex.voxelShape.isGpuCcdSupported()) {
				return false;
			}
		}
		
		return true;
	}
}
