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
