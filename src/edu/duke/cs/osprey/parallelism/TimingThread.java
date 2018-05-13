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
