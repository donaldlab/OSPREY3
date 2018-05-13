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

package edu.duke.cs.osprey.tools;

import org.junit.Test;

import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;

public class TestExpFunction {

	@Test
	public void exp() {
		final double Epsilon = 1e-8;
		ExpFunction exp = new ExpFunction();
		assertThat(exp.exp(-10.0).doubleValue(), isAbsolutely(Math.exp(-10.0), Epsilon));
		assertThat(exp.exp( -5.0).doubleValue(), isAbsolutely(Math.exp( -5.0), Epsilon));
		assertThat(exp.exp( -1.0).doubleValue(), isAbsolutely(Math.exp( -1.0), Epsilon));
		assertThat(exp.exp(  0.0).doubleValue(), isAbsolutely(Math.exp(  0.0), Epsilon));
		assertThat(exp.exp(  1.0).doubleValue(), isAbsolutely(Math.exp(  1.0), Epsilon));
		assertThat(exp.exp(  5.0).doubleValue(), isAbsolutely(Math.exp(  5.0), Epsilon));
		assertThat(exp.exp( 10.0).doubleValue(), isAbsolutely(Math.exp( 10.0), Epsilon));
	}
}
