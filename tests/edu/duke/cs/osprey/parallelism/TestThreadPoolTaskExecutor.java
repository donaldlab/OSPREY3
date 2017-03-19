package edu.duke.cs.osprey.parallelism;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

public class TestThreadPoolTaskExecutor {
	
	@Test
	public void countToTen() {
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		int[] count = { 0 };
		
		for (int i=0; i<10; i++) {
			tasks.submit(() -> {
				// no work to do
			}, (task) -> {
				// increment the counter on the listener thread
				count[0]++;
			});
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
				tasks.submit(() -> {
					// on worker thread: no work to do
				}, (task) -> {
					// on listener thread: increment counter
					count[0]++;
				});
			}
			tasks.waitForFinish();
		
			assertThat(count[0], is(4));
		}
	}
}
