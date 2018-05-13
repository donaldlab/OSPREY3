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

public class CounterSignal {
	
	public static interface SignalCondition {
		boolean shouldSignal(int counter);
	}
	
	private int count;
	private SignalCondition condition;
	
	public CounterSignal(int initialCount, SignalCondition condition) {
		this.count = initialCount;
		this.condition = condition;
	}
	
	public synchronized void waitForSignal() {
		waitForSignal(0);
	}
	
	public synchronized void waitForSignal(long timeoutMs) {
		if (!condition.shouldSignal(count)) {
			try {
				wait(timeoutMs);
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}
	
	public synchronized void offset(int delta) {
		count += delta;
		if (condition.shouldSignal(count)) {
			notifyAll();
		}
	}
}
