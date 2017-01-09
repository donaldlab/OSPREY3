package edu.duke.cs.osprey.minimization;

public class IntervalLineSearcher implements LineSearcher {
	
	private static final double Tau = (Math.sqrt(5) + 1)/2;
	private static final double MinIntervalWidth = 1e-3; // TODO: calibrate for different DOFs
	
	private ObjectiveFunction.OneDof f;

	@Override
	public void init(ObjectiveFunction.OneDof f) {
		this.f = f;
	}

	@Override
	public double search(double xd) {
		
		// use golden interval minimization
		// it won't necessarily pick the closest minimum,
		// but it will converge to one of them arbitrarily
		// https://en.wikipedia.org/wiki/Golden_section_search
		
		// start with the interval equal to the full range
		double a = f.getXMin();
		double b = f.getXMax();
		
		while (Math.abs(a - b) > MinIntervalWidth) {
			
			// use the golden ratio to sample points in the interior of the interval
			double step = (b - a)/Tau;
			double c = b - step;
			double d = a + step;
			
			// update interval bounds based on which side is lower
			double fc = f.getValue(d);
			double fd = f.getValue(d);
			if (fc < fd) {
				b = d;
			} else {
				a = c;
			}
		}
		
		// the interval is small now, pick the center as the final value
		double x = (a + b)/2;
		f.setX(x);
		return x;
	}
}
