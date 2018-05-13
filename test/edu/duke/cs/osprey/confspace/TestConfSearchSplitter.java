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

package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.hamcrest.Matchers.is;
import static org.hamcrest.Matchers.nullValue;
import static org.hamcrest.core.IsNull.notNullValue;
import static org.junit.Assert.assertThat;

public class TestConfSearchSplitter {

	private static SimpleConfSpace confSpace;
	private static EnergyMatrix emat;
	private static List<ConfSearch.ScoredConf> expectedConfs;

	@BeforeClass
	public static void beforeClass() {

		// make a conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7")) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}

		// get all of the expected conformations (there are 27)
		expectedConfs = makeSearch().nextConfs(Double.POSITIVE_INFINITY);
	}

	private static ConfSearch makeSearch() {

		// make any conf search, doesn't matter
		return new ConfAStarTree.Builder(emat, confSpace).build();
	}

	private static class Checker {

		int index = 0;
		ConfSearch confs;

		Checker(ConfSearch confs) {
			this.confs = confs;
		}

		void assertConfs(int n) {
			for (int i=0; i<n; i++) {
				ConfSearch.ScoredConf conf = confs.nextConf();
				assertThat(conf, is(notNullValue()));
				assertThat(conf, is(expectedConfs.get(index++)));
			}
		}

		void assertEnd() {
			assertThat(confs.nextConf(), is(nullValue()));
		}
	}

	@Test
	public void someFirstSomeSecond() {

		ConfSearch.Splitter splitter = new ConfSearch.Splitter(makeSearch());

		Checker first = new Checker(splitter.first);
		Checker second = new Checker(splitter.second);

		first.assertConfs(10);
		second.assertConfs(10);

		first.assertConfs(10);
		second.assertConfs(10);

		first.assertConfs(7);
		first.assertEnd();
		second.assertConfs(7);
		second.assertEnd();
	}

	@Test
	public void allFirstAllSecond() {

		ConfSearch.Splitter splitter = new ConfSearch.Splitter(makeSearch());

		Checker first = new Checker(splitter.first);
		first.assertConfs(27);
		first.assertEnd();

		Checker second = new Checker(splitter.second);
		second.assertConfs(27);
		second.assertEnd();
	}

	@Test(expected=ConfSearch.Splitter.OutOfOrderException.class)
	public void outOfOrderDirect() {

		ConfSearch.Splitter splitter = new ConfSearch.Splitter(makeSearch());

		Checker second = new Checker(splitter.second);
		second.assertConfs(1);
	}

	@Test(expected=ConfSearch.Splitter.OutOfOrderException.class)
	public void outOfOrderDelayed() {

		ConfSearch.Splitter splitter = new ConfSearch.Splitter(makeSearch());

		Checker first = new Checker(splitter.first);
		first.assertConfs(10);

		Checker second = new Checker(splitter.second);
		second.assertConfs(10);
		second.assertConfs(1);
	}

	@Test
	public void externalMemory() {
		ExternalMemory.use(16, () -> {

			ConfSearch.Splitter splitter = new ConfSearch.Splitter(
				makeSearch(),
				true,
				new RCs(confSpace)
			);

			Checker first = new Checker(splitter.first);
			first.assertConfs(10);

			Checker second = new Checker(splitter.second);
			second.assertConfs(10);
		});
	}

	@Test
	public void externalMemoryExhaust() {
		ExternalMemory.use(16, () -> {

			ConfSearch.Splitter splitter = new ConfSearch.Splitter(
				makeSearch(),
				true,
				new RCs(confSpace)
			);

			Checker first = new Checker(splitter.first);
			first.assertConfs(27);

			Checker second = new Checker(splitter.second);
			second.assertConfs(27);

			assertThat(splitter.first.nextConf(), is(nullValue()));
			assertThat(splitter.second.nextConf(), is(nullValue()));
		});
	}
}
