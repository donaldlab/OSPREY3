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

package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;

@Deprecated
public interface GMECConfEnergyCalculator {
	
	EnergiedConf calcEnergy(ScoredConf conf);

	// use asynchronous techniques so we can parallelize conformation evaluation
	public static interface Async extends GMECConfEnergyCalculator {

		void calcEnergyAsync(ScoredConf conf, Listener listener);
		TaskExecutor getTasks();

		public static interface Listener extends TaskListener<EnergiedConf> {
			// nothing else to do
		}

		public static class Adapter implements Async {

			private GMECConfEnergyCalculator calc;
			private TaskExecutor tasks;

			public Adapter(GMECConfEnergyCalculator calc) {
				this.calc = calc;
				this.tasks = new TaskExecutor();
			}

			@Override
			public EnergiedConf calcEnergy(ScoredConf conf) {
				return calc.calcEnergy(conf);
			}

			@Override
			public void calcEnergyAsync(ScoredConf conf, Listener listener) {
				listener.onFinished(calc.calcEnergy(conf));
			}
			
			@Override
			public TaskExecutor getTasks() {
				return tasks;
			}
		}
	}
}
