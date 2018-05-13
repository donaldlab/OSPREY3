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

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

public class TestTimeTools {

	@Test
	public void test()
	throws Exception {

		long ms = TimeTools.getTimestampMs();
		long us = TimeTools.getTimestampUs();
		long ns = TimeTools.getTimestampNs();

		assertThat(us/1000L - ms, lessThanOrEqualTo(1L));
		assertThat(ns/1000000L - ms, lessThanOrEqualTo(1L));

		Thread.sleep(37);

		ms = TimeTools.getTimestampMs();
		us = TimeTools.getTimestampUs();
		ns = TimeTools.getTimestampNs();

		assertThat(us/1000L - ms, lessThanOrEqualTo(1L));
		assertThat(ns/1000000L - ms, lessThanOrEqualTo(1L));

		Thread.sleep(379);

		ms = TimeTools.getTimestampMs();
		us = TimeTools.getTimestampUs();
		ns = TimeTools.getTimestampNs();

		assertThat(us/1000L - ms, lessThanOrEqualTo(1L));
		assertThat(ns/1000000L - ms, lessThanOrEqualTo(1L));
	}
}
