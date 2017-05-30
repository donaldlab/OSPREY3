package edu.duke.cs.osprey.astar;

import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;

public class Matchers {
	
	public static Matcher<int[]> startsWith(int ... expected) {
		return new BaseMatcher<int[]>() {

			@Override
			public boolean matches(Object obj) {
				int[] observed = (int[]) obj;
				for (int i=0; i<expected.length; i++) {
					if (expected[i] != observed[i]) {
						return false;
					}
				}
				return true;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("starts with ")
					.appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				int[] observed = (int[]) obj;
				desc.appendText("array is ")
					.appendValue(observed);
			}
		};
	}
}
