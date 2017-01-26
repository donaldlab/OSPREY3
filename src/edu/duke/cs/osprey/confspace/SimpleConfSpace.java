package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
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
	
	public static class Position {
		
		public final int index;
		public final Strand strand;
		public final Residue res;
		public final List<ResidueConf> resConfs;
		
		public Position(int index, Strand strand, Residue residue) {
			this.index = index;
			this.strand = strand;
			this.res = residue;
			this.resConfs = new ArrayList<>();
		}
	}
	
	public static class ResidueConf {
		
		public static enum Type {
			Library,
			WildType;
		}
		
		public final int index;
		public final ResidueTemplate template;
		public final Type type;
		public final Integer rotamerIndex;
		
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
		}
	}
	
	public final List<Strand> strands;
	public final List<Position> positions;
	
	private Map<Strand,List<StrandFlex>> strandFlex; // yeah, map on instance, not identity
	
	public SimpleConfSpace() {
		strands = new ArrayList<>();
		positions = new ArrayList<>();
		
		strandFlex = new HashMap<>();
	}
	
	public SimpleConfSpace(Strand ... strands) {
		this();
		for (Strand strand : strands) {
			addStrand(strand);
		}
	}
	
	public void addStrand(Strand strand, StrandFlex ... flexType) {
		
		// add the strand
		strands.add(strand);
		strandFlex.put(strand, Arrays.asList(flexType));
		
		// make the positions
		for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
			Residue res = strand.mol.getResByPDBResNumber(resNum);
			
			// make the pos
			Position pos = new Position(positions.size(), strand, res);
			positions.add(pos);
			
			// make the residue confs
			Strand.ResidueFlex resFlex = strand.flexibility.get(resNum);
			
			for (String resType : resFlex.resTypes) {
				makeResidueConfsFromTemplate(pos, strand.templateLib.getTemplate(resType), ResidueConf.Type.Library);
			}
			
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
		
		if (template.getNumRotamers() <= 0) {
			
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

	/**
	 * create a new {@link ParametricMolecule} in the specified conformation
	 * for analysis (e.g., minimization)
	 * 
	 * To increase stability of analysis, each analysis should be conducted
	 * with a new molecule instance. this completely prevents roundoff error
	 * from accumulating across separate analyses. 
	 */
	public ParametricMolecule makeMolecule(int[] conf) {
		
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
		for (Position pos : positions) {
			ResidueConf resConf = pos.resConfs.get(conf[pos.index]);
			Residue res = mol.getResByPDBResNumber(pos.res.getPDBResNumber());
			ResidueTypeDOF.switchToTemplate(pos.strand.templateLib, res, resConf.template, false);
		}
		
		// make all the DOFs
		List<DegreeOfFreedom> dofs = new ArrayList<>();
		for (Strand strand : strands) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				dofs.addAll(flex.makeDofs(strand));
			}
			// NOTE: strand DOFs are "centered" on initial molecule pos, so there's no need to set anything here to match the conf
		}
		for (Position pos : positions) {
			ResidueConf resConf = pos.resConfs.get(conf[pos.index]);
			Residue res = mol.getResByPDBResNumber(pos.res.getPDBResNumber());
			
			// make the residue DOFs
			Strand.ResidueFlex resFlex = pos.strand.flexibility.get(res.getPDBResNumber());
			dofs.addAll(resFlex.voxelShape.makeDihedralDOFs(res));
			
			// pose the residue to match the rotamer
			List<DegreeOfFreedom> dihedralDofs = new VoxelShape.Rect().makeDihedralDOFs(res);
			for (int d=0; d<resConf.template.numDihedrals; d++) {
				dihedralDofs.get(d).apply(resConf.template.getRotamericDihedrals(resConf.rotamerIndex, d));
			}
		}
		
		// TODO: make DOFs for non-residue continuous motions (eg, blocks, DEEPer)
		
		return new ParametricMolecule(mol, dofs);
	}
	
	public DofBounds makeBounds(int[] conf) {
		
		// gather all the DOF bounds
		List<DofBounds> bounds = new ArrayList<>();
		for (Strand strand : strands) {
			for (StrandFlex flex : strandFlex.get(strand)) {
				bounds.add(flex.makeBounds(strand));
			}
		}
		for (Position pos : positions) {
			ResidueConf resConf = pos.resConfs.get(conf[pos.index]);
			if (resConf.rotamerIndex != null) {
				Strand.ResidueFlex resFlex = pos.strand.flexibility.get(pos.res.getPDBResNumber());
				bounds.add(resFlex.voxelShape.makeDihedralBounds(resConf.template, resConf.rotamerIndex));
			}
		}
		
		// TODO: make DOF bounds for non-residue continuous motions (eg, blocks, DEEPer)
		
		return DofBounds.concatenate(bounds);
	}
}
