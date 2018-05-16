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
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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




