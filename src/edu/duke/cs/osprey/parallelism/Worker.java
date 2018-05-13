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

package edu.duke.cs.osprey.parallelism;

public abstract class Worker extends WorkThread {
	
	private WorkCrew<Worker> crew;
	
	// protected so the WorkCrew can use it
	protected boolean hasWork;
	
	@SuppressWarnings("unchecked")
	public Worker(WorkCrew<? extends Worker> crew) {
		super(crew.getName() + "-" + crew.workers.size());
		this.crew = (WorkCrew<Worker>)crew;
		this.crew.workers.add(this);
		hasWork = false;
	}
	
	@Override
	public void doWork()
	throws InterruptedException {
		
		// is there any work to do?
		if (hasWork) {
			
			workIt();
			hasWork = false;
			
			crew.finishedWork();
		}
		
		// wait until we get more work
		// but check the isRunning flag every second or so
		crew.waitForWork(this, 1000);
	}
	
	protected abstract void workIt();
}
