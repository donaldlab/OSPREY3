package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.SimpleConfSpace.DofTypes;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;

public abstract class StrandFlex {
	
	public abstract List<DegreeOfFreedom> makeDofs(Strand strand);
	public abstract DofBounds makeBounds(Strand strand);
	
	public DofTypes getDofTypes() {
		return DofTypes.Any;
	}
	
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
		public List<DegreeOfFreedom> makeDofs(Strand strand) {
			return new MoveableStrand(strand.mol.residues).getDOFs();
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
	}
}
