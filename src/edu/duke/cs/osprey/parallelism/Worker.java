package edu.duke.cs.osprey.parallelism;

public abstract class Worker extends WorkThread {
	
	private WorkCrew<Worker> crew;
	
	// protected so the WorkCrew can use it
	protected boolean hasWork;
	
	@SuppressWarnings("unchecked")
	public Worker(WorkCrew<? extends Worker> crew) {
		super(crew.getName() + "-" + crew.workers.size());
		this.crew = (WorkCrew<Worker>)crew;
		this.crew.workers.add(this);
		hasWork = false;
	}
	
	@Override
	public void doWork()
	throws InterruptedException {
		
		// is there any work to do?
		if (hasWork) {
			
			workIt();
			hasWork = false;
			
			crew.finishedWork();
		}
		
		// wait until we get more work
		// but check the isRunning flag every second or so
		crew.waitForWork(this, 1000);
	}
	
	protected abstract void workIt();
}
