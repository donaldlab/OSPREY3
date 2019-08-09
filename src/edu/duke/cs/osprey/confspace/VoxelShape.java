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

	public abstract int countDihedralDOFs(ResidueTemplate template);
	
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
		public int countDihedralDOFs(ResidueTemplate template) {
			return 0;
		}

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
		public int countDihedralDOFs(ResidueTemplate template) {
			return template.numDihedrals;
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
