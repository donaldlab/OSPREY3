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
