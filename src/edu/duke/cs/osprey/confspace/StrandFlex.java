package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;

public abstract class StrandFlex {
	
	public abstract List<DegreeOfFreedom> makeDofs(Strand strand);
	public abstract DofBounds makeBounds(Strand strand);
	
	public static class None extends StrandFlex {
		
		@Override
		public List<DegreeOfFreedom> makeDofs(Strand strand) {
			// no movement, no DOFs
			return new ArrayList<>();
		}

		@Override
		public DofBounds makeBounds(Strand strand) {
			// no movement, no DOFs
			return new DofBounds(0);
		}
	}
	
	public static class TranslateRotate extends StrandFlex {
		
		public static final double DefaultMaxRotDegrees = 5;
		public static final double DefaultMaxTranslation = 1.2;
		
		@Override
		public List<DegreeOfFreedom> makeDofs(Strand strand) {
			return new MoveableStrand(strand.mol.residues).getDOFs();
		}

		@Override
		public DofBounds makeBounds(Strand strand) {
			DofBounds bounds = new DofBounds(6);
			bounds.set(0, -DefaultMaxRotDegrees, DefaultMaxRotDegrees);
			bounds.set(1, -DefaultMaxRotDegrees, DefaultMaxRotDegrees);
			bounds.set(2, -DefaultMaxRotDegrees, DefaultMaxRotDegrees);
			bounds.set(3, -DefaultMaxTranslation, DefaultMaxTranslation);
			bounds.set(4, -DefaultMaxTranslation, DefaultMaxTranslation);
			bounds.set(5, -DefaultMaxTranslation, DefaultMaxTranslation);
			return bounds;
		}
	}
}
