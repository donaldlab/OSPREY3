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
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;

public abstract class StrandFlex implements Serializable {
	
	public abstract List<? extends DegreeOfFreedom> makeDofs(Strand strand, Molecule mol);//DOF will control the conformation of mol
	public abstract DofBounds makeBounds(Strand strand);
        
        public abstract ArrayList<HashMap<String,double[]> > listBackboneVoxels(Position pos);
        //list the backbone voxels (in DOF name --> bounds form) for a particular position
	
	public DofTypes getDofTypes() {
		return DofTypes.Any;
	}
        
        protected HashMap<String, double[]> defaultBackboneVoxel(Strand strand){
            List<? extends DegreeOfFreedom> dofs = makeDofs(strand, strand.mol);
            //can make DOFs in strand molec because just need their names
            DofBounds bounds = makeBounds(strand);
            HashMap<String,double[]> vox = new HashMap<>();
            for(int d=0; d<dofs.size(); d++)
                vox.put(dofs.get(d).getName(), new double[]{bounds.getMin(d),bounds.getMax(d)});
            return vox;
        }
        	
	public static class None extends StrandFlex {
		
		@Override
		public List<DegreeOfFreedom> makeDofs(Strand strand, Molecule mol) {
			// no movement, no DOFs
			return new ArrayList<>();
		}

		@Override
		public DofBounds makeBounds(Strand strand) {
			// no movement, no DOFs
			return new DofBounds(0);
		}

                @Override
                public ArrayList<HashMap<String, double[]>> listBackboneVoxels(Position pos) {
                        return new ArrayList<>();
                }
	}
	
	public static class TranslateRotate extends StrandFlex {
		
		public static final double DefaultMaxRotDegrees = 5;
		public static final double DefaultMaxTranslation = 1.2;
		
		private double maxRotDegrees;
		private double maxTranslation;
		
		public TranslateRotate() {
			this(DefaultMaxRotDegrees, DefaultMaxTranslation);
		}
		
		public TranslateRotate(double maxRotDegrees, double maxTranslation) {
			this.maxRotDegrees = maxRotDegrees;
			this.maxTranslation = maxTranslation;
		}
		
		@Override
		public List<DegreeOfFreedom> makeDofs(Strand strand, Molecule mol) {
                    //Need to do move the new (copied) molecule rather than the one in the strand
                    ArrayList<Residue> movingResidues = new ArrayList<>();
                    for(Residue origRes : strand.mol.residues)
                        movingResidues.add(origRes.equivalentInMolec(mol));
                    return new MoveableStrand(movingResidues).getDOFs();
		}

		@Override
		public DofBounds makeBounds(Strand strand) {
			DofBounds bounds = new DofBounds(6);
			bounds.set(0, -maxRotDegrees, maxRotDegrees);
			bounds.set(1, -maxRotDegrees, maxRotDegrees);
			bounds.set(2, -maxRotDegrees, maxRotDegrees);
			bounds.set(3, -maxTranslation, maxTranslation);
			bounds.set(4, -maxTranslation, maxTranslation);
			bounds.set(5, -maxTranslation, maxTranslation);
			return bounds;
		}
                
                @Override
                public ArrayList<HashMap<String, double[]>> listBackboneVoxels(Position pos) {
                        //just one voxel
                        return new ArrayList(Arrays.asList(defaultBackboneVoxel(pos.strand)));
                }
	}
}
