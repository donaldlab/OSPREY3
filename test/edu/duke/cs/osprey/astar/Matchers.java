/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

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
