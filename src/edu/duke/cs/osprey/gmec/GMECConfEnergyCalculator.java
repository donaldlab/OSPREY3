/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
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
