package edu.duke.cs.osprey;

import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;

public class TestBase {
	
	private static final double Epsilon = 1e-14;

	protected double getRelativeError(double expected, double observed) {
		return Math.abs(expected - observed)/Math.abs(observed);
	}
	
	protected Matcher<Double> isRelatively(final double expected) {
		return new BaseMatcher<Double>() {

			@Override
			public boolean matches(Object obj) {
				double observed = ((Double)obj).doubleValue();
				return getRelativeError(expected, observed) <= Epsilon;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				double observed = ((Double)obj).doubleValue();
				double relErr = getRelativeError(expected, observed);
				desc.appendText("value ").appendValue(observed)
					.appendText(" has relative err ").appendValue(relErr)
					.appendText(" that's greater than epsilon ").appendValue(Epsilon);
			}
		};
	}
}
