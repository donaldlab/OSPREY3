package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.SimpleConfSpace.DofTypes;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Residue;

public abstract class VoxelShape {
	
	public static final double DefaultWidthDegrees = 9;
	
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
		
		public final double width;
		
		public Rect() {
			this(DefaultWidthDegrees);
		}
		
		public Rect(double width) {
			this.width = width;
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
				bounds.set(d, chi - width, chi + width);
			}
			return bounds;
		}

		@Override
		public DofTypes getDofTypes() {
			return DofTypes.OnlyDihedrals;
		}
	}
	
	public static class Ellipse extends VoxelShape {

		public final double width;
		
		public Ellipse() {
			this(DefaultWidthDegrees);
		}
		
		public Ellipse(double width) {
			this.width = width;
		}
		
		@Override
		public List<DegreeOfFreedom> makeDihedralDOFs(Residue res) {
			
			Rect rectVoxel = new Rect(width);
			List<DegreeOfFreedom> dihedrals = rectVoxel.makeDihedralDOFs(res);
			
			// TODO: build elliptical DOFs
			throw new UnsupportedOperationException("not implement yet");
		}

		@Override
		public DofBounds makeDihedralBounds(ResidueTemplate template, int rotamerIndex) {
			// TODO: implement me
			throw new UnsupportedOperationException("not implement yet");
		}

		@Override
		public DofTypes getDofTypes() {
			// TODO: implement me
			throw new UnsupportedOperationException("not implement yet");
		} 
	}
}
