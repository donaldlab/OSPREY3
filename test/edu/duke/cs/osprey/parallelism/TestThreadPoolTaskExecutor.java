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
	
	@Test(expected = TaskException.class)
	public void handleTaskExceptionsGracefully() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
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
	}
	
	@Test(expected = TaskException.class)
	public void handleListenerExceptionsGracefully() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		for (int i=0; i<10; i++) {

			tasks.submit(
				() -> {
					// we can do it!
					return null;
				},
				(Void ignore) -> {
					// crash in the listener
					throw new Error("Oh No! a Bad Thing has happened");
				}
			);
		}
		
		tasks.waitForFinish();
	}
}
