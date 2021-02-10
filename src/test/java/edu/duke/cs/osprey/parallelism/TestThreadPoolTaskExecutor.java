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

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskException;

public class TestThreadPoolTaskExecutor {
	
	@Test
	public void countToTen() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		int[] count = { 0 };
		
		for (int i=0; i<10; i++) {
			tasks.submit(
				() -> {
					// no work to do
					return null;
				},
				(Void ignore) -> {
					// increment the counter on the listener thread
					count[0]++;
				}
			);
		}
		tasks.waitForFinish();
		
		assertThat(count[0], is(10));
	}
	
	@Test
	public void countLotsOfTimes() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(4);
		
		for (int r=0; r<1000; r++) {
			
			int[] count = { 0 };
			
			for (int i=0; i<4; i++) {
				tasks.submit(
					() -> {
						// on worker thread: no work to do
						return null;
					},
					(Void ignore) -> {
						// on listener thread: increment counter
						count[0]++;
					}
				);
			}
			tasks.waitForFinish();
		
			assertThat(count[0], is(4));
		}
	}
	
	@Test
	public void handleTaskExceptionsGracefully() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		for (int r=0; r<100; r++) {
		
			try {
				for (int i=0; i<10; i++) {
					tasks.submit(
						() -> {
							// crash in the task
							throw new Error("Oh No! a Bad Thing has happened");
						},
						(Void ignore) -> {
							fail("task should not finish");
						}
					);
				}
				tasks.waitForFinish();
				
				fail("should have thrown Error");
				
			} catch (TaskException ex) {
				
				assertThat(tasks.getNumRunningTasks(), is(0L));
				
				// all is well
				continue;
			}
		}
	}
	
	@Test
	public void handleListenerExceptionsGracefully() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		for (int r=0; r<100; r++) {
		
			try {
				for (int i=0; i<10; i++) {
					tasks.submit(
						() -> {
							// easiest task ever!
							return null;
						},
						(Void ignore) -> {
							// crash in the listener
							throw new Error("Oh No! a Bad Thing has happened");
						}
					);
				}
				tasks.waitForFinish();
				
				fail("should have thrown error");
				
			} catch (TaskException ex) {
				
				assertThat(tasks.getNumRunningTasks(), is(0L));
				
				// all is well
				continue;
			}
		}
	}
}
