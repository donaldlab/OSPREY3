package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.SimpleConfSpace.DofTypes;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;

public abstract class VoxelShape implements Serializable {
	
	public static final double DefaultHalfWidthDegrees = 9;
	
	/**
	 * make degrees of freedom for the residue in its current state (ie, ignore possible mutations)
	 */
	public abstract List<DegreeOfFreedom> makeDihedralDOFs(Residue res);
	
	/**
	 * make bounds for degrees of freedom for the residue in its current state (ie, ignore possible mutations)
	 */
	public abstract DofBounds makeDihedralBounds(ResidueTemplate template, int rotamerIndex);
	
	public abstract DofTypes getDofTypes();

	public static class Point extends VoxelShape {

		@Override
		public List<DegreeOfFreedom> makeDihedralDOFs(Residue res) {
			
			// no continuous degrees of freedom
			return new ArrayList<>();
		}

		@Override
		public DofBounds makeDihedralBounds(ResidueTemplate template, int rotamerIndex) {
			
			// no continuous degrees of freedom
			return new DofBounds(0);
		}

		@Override
		public DofTypes getDofTypes() {
			return DofTypes.None;
		} 
	}
	
	public static class Rect extends VoxelShape {
		
		public final double halfWidth;
		
		public Rect() {
			this(DefaultHalfWidthDegrees);
		}
		
		public Rect(double halfWidth) {
			this.halfWidth = halfWidth;
		}

		@Override
		public List<DegreeOfFreedom> makeDihedralDOFs(Residue res) {
			List<DegreeOfFreedom> dihedrals = new ArrayList<>();
			for (int i=0; i<res.template.numDihedrals; i++) {
				dihedrals.add(new FreeDihedral(res, i));
			}
			return dihedrals;
		}

		@Override
		public DofBounds makeDihedralBounds(ResidueTemplate template, int rotamerIndex) {
			DofBounds bounds = new DofBounds(template.numDihedrals);
			for (int d=0; d<template.numDihedrals; d++) {
				double chi = template.getRotamericDihedrals(rotamerIndex, d);
				bounds.set(d, chi - halfWidth, chi + halfWidth);
			}
			return bounds;
		}

		@Override
		public DofTypes getDofTypes() {
			return DofTypes.OnlyDihedrals;
		}
	}
}
