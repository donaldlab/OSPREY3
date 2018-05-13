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

import java.util.ArrayList;
import java.util.Iterator;

public class WorkCrew<T extends Worker> implements Iterable<T> {
	
	public static @interface CalledByWorker {}
	public static @interface CalledByBoss {}
	
	// protected so Worker can access it
	protected ArrayList<T> workers;
	
	private String name;
	private Object workSync;
	private Object resultsSync;
	private volatile int counter;
	
	@CalledByBoss
	public WorkCrew(String name) {
		this.name = name;
		workers = new ArrayList<>();
		workSync = new Object();
		resultsSync = new Object();
		counter = 0;
	}
	
	public String getName() {
		return name;
	}
	
	public ArrayList<T> getWorkers() {
		return workers;
	}
	
	@Override
	public Iterator<T> iterator() {
		return workers.iterator();
	}
	
	public void start() {
		for (Worker worker : workers) {
			worker.start();
		}
	}
	
	public void askToStop() {
		for (Worker worker : workers) {
			worker.askToStop();
		}
	}
	
	@SuppressWarnings("deprecation")
	public void killThreads() {
		for (Worker worker : workers) {
			worker.stop();
		}
	}

	@CalledByWorker
	public void waitForWork(Worker processor, int timeoutMs)
	throws InterruptedException {
		if (processor.hasWork) {
			return;
		}
		synchronized (workSync) {
			// NOTE: after we waited to synchronize this, we might have been given more work
			// if that's the case, don't wait, just exit and do the work
			if (processor.hasWork) {
				return;
			}
			workSync.wait(timeoutMs);
		}
	}
	
	@CalledByBoss
	public void sendWork() {
		counter = 0;
		for (Worker worker : workers) {
			worker.hasWork = true;
		}
		synchronized (workSync) {
			workSync.notifyAll();
		}
	}
	
	@CalledByBoss
	public boolean waitForResults(int timeoutMs)
	throws InterruptedException {
		if (counter == workers.size()) {
			return true;
		}
		synchronized (resultsSync) {
			// NOTE: after waiting to synchronize this, the processors might have finished
			// if that's the case, don't wait
			if (counter == workers.size()) {
				return true;
			}
			resultsSync.wait(timeoutMs);
		}
		return counter == workers.size();
	}
	
	@CalledByWorker
	public void finishedWork() {
		int c;
		synchronized (this) {
			counter++;
			c = counter;
		}
		if (c == workers.size()) {
			synchronized (resultsSync) {
				resultsSync.notify();
			}
		}
	}
}
