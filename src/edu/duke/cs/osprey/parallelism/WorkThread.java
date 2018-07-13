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

import edu.duke.cs.tpie.Cleaner.Cleanable;

public abstract class WorkThread extends Thread implements Cleanable {
	
	// this flag is hit from multiple threads concurrently, so make it volatile
	private volatile boolean isRunning;
	
	protected WorkThread(String name) {
		super(name);
		setDaemon(true);
		isRunning = false;
	}
	
	public boolean isRunning() {
		return isRunning();
	}
	
	public void askToStop() {
		isRunning = false;
	}
	
	public void askToStopAndWait()
	throws InterruptedException {
		askToStop();
		join();
	}
	
	@Override
	public void run() {
		
		isRunning = true;
		
		while (isRunning) {
			try {
				doWork();
			} catch (InterruptedException ex) {
				break;
			}
		}
		
		isRunning = false;
	}
	
	protected abstract void doWork() throws InterruptedException;
	
	@Override
	public void clean() {
		askToStop();
	}
}
