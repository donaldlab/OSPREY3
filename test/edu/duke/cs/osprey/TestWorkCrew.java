package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;

import edu.duke.cs.osprey.parallelism.WorkCrew;
import edu.duke.cs.osprey.parallelism.Worker;

public class TestWorkCrew {
	
	private static class IncrementWorker extends Worker {
		
		public int val;
		
		protected long numWork;

		public IncrementWorker(WorkCrew<IncrementWorker> crew, long numWork) {
			super(crew);
			this.val = 0;
			this.numWork = numWork;
		}

		@Override
		protected void workIt() {
			
			// really simple workload
			
			// NOTE: for some reason an int-type loop variable always takes
			// less than 10 ms for any size workload that fits in an int,
			// so use long-type for loop variable
			
			// typical times for workloads:
			// 10000000000L: 5 s
			// 1000000000L:  774 ms
			// 100000000L:   88 ms
			// 10000000L:    12 ms
			// 1000000L:     4 ms
			// 100000L:      1.5 ms
			// 10000L:       150 us
			// 1000L:        15 us
			
			for (long i=0; i<numWork; i++);
			
			val++;
		}
	}
	
	private static class RandomIncrementWorker extends IncrementWorker {
		
		private int[] workloads;
		private int workIndex;

		public RandomIncrementWorker(WorkCrew<IncrementWorker> crew, int maxWork) {
			super(crew, 0);
			
			workloads = new int[10000];
			Random rand = new Random(12345 * crew.getWorkers().size());
			for (int i=0; i<workloads.length; i++) {
				workloads[i] = rand.nextInt(maxWork);
			}
			workIndex = 0;
		}
		
		@Override
		protected void workIt() {
			
			// do a random amount of work
			workIndex = (workIndex + 1) % workloads.length;
			this.numWork = workloads[workIndex];
			
			super.workIt();
		}
	}
		
	
	private void testWorkCrew(WorkCrew<IncrementWorker> crew, int numRounds)
	throws Exception {
		
		crew.start();
		
		// everyone should start at 0
		for (IncrementWorker worker : crew) {
			assertThat(worker.val, is(0));
		}
		
		for (int i=0; i<numRounds; i++) {
			
			// do the parallel work
			crew.sendWork();
			boolean isFinished = crew.waitForResults(1000);
			assertThat(isFinished, is(true));
			
			// check the result
			for (IncrementWorker worker : crew) {
				assertThat(worker.val, is(i + 1));
			}
		}
		
		crew.askToStop();
	}
	
	@Test
	public void testOneRound()
	throws Exception {
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new IncrementWorker(crew, 0);
		}
		testWorkCrew(crew, 1);
	}
	
	@Test
	public void testManyRoundsLightWorkload()
	throws Exception {
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new IncrementWorker(crew, 10);
		}
		testWorkCrew(crew, 20000);
	}
	
	@Test
	public void testManyRoundsMediumWorkload()
	throws Exception {
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new IncrementWorker(crew, 10000);
		}
		testWorkCrew(crew, 20000);
	}
	
	@Test
	public void testManyRoundsHeavyWorkload()
	throws Exception {
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new IncrementWorker(crew, 1000000);
		}
		testWorkCrew(crew, 1000);
	}
	
	@Test
	public void testManyRoundsRandomWorkload()
	throws Exception {
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new RandomIncrementWorker(crew, 100000);
		}
		testWorkCrew(crew, 10000);
	}
	
	@Test
	public void testTimeout()
	throws Exception {
		
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new IncrementWorker(crew, Long.MAX_VALUE);
		}
		
		crew.start();
		crew.sendWork();
		boolean isFinished = crew.waitForResults(1000);
		crew.askToStop();
		assertThat(isFinished, is(false));
		
		// force the threads to stop
		crew.killThreads();
	}
	
	// uncomment this test for deep debugging of hard-to-reproduce concurrency issues
	// it takes way too long for typical regression testing though
	//@Test
	public void testManyManyManyRoundsRandomWorkload()
	throws Exception {
		WorkCrew<IncrementWorker> crew = new WorkCrew<>("Test");
		for (int i=0; i<4; i++) {
			new RandomIncrementWorker(crew, 100000000);
		}
		testWorkCrew(crew, 1000);
	}
}
