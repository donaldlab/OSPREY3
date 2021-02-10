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

import java.util.List;

import edu.duke.cs.osprey.tools.Stopwatch;

public abstract class TimingThread extends Thread {
	
	public Signal warmupSignal;
	public Signal goSignal;
	public Signal doneSignal;
	
	public TimingThread(String name) {
		super(name);
		setDaemon(true);
		
		warmupSignal = new Signal();
		goSignal = new Signal();
		doneSignal = new Signal();
	}
	
	protected void init() throws Exception {
		// nothing to do by default
	}
	
	protected abstract void warmup();
	protected abstract void time();
	
	protected void cleanup() {
		// nothing to do by default
	}
	
	@Override
	public void run() {
		try {
			init();
		} catch (Exception ex) {
			throw new Error(ex);
		}
		warmup();
		warmupSignal.sendSignal();
		goSignal.waitForSignal();
		time();
		doneSignal.sendSignal();
		cleanup();
	}
	
	public void waitForFinish() {
		try {
			join();
		} catch (InterruptedException ex) {
			throw new Error(ex);
		}
	}
	
	public static Stopwatch timeThreads(List<TimingThread> threads) {
	
		// start the threads and do the warmup
		for (TimingThread thread : threads) {
			thread.start();
		}
		
		// wait for the warmup to finish
		for (TimingThread thread : threads) {
			thread.warmupSignal.waitForSignal();
		}
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		// start the main event
		for (TimingThread thread : threads) {
			thread.goSignal.sendSignal();
		}
		
		// wait for timing to finish
		for (TimingThread thread : threads) {
			thread.doneSignal.waitForSignal();
		}
		
		stopwatch.stop();
		
		// wait for cleanup
		for (TimingThread thread : threads) {
			thread.waitForFinish();
		}
		
		return stopwatch;
	}
}
