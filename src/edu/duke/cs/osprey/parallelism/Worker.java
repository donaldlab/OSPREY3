package edu.duke.cs.osprey.parallelism;

public abstract class Worker extends Thread {
	
	// this flag is hit from multiple threads concurrently, so make it volatile
	private volatile boolean isRunning;
	
	private WorkCrew<Worker> crew;
	
	// protected so the WorkCrew can use it
	protected boolean hasWork;
	
	@SuppressWarnings("unchecked")
	public Worker(WorkCrew<? extends Worker> crew) {
		super(crew.getName() + "-" + crew.workers.size());
		this.crew = (WorkCrew<Worker>)crew;
		this.crew.workers.add(this);
		setDaemon(true);
		hasWork = false;
	}
	
	@Override
	public void run() {
		
		isRunning = true;
		while (isRunning) {
			
			// is there any work to do?
			if (hasWork) {
				
				workIt();
				hasWork = false;
				
				crew.finishedWork();
			}
			
			// wait until we get more work
			// but check the isRunning flag every second or so
			try {
				crew.waitForWork(this, 1000);
			} catch (InterruptedException ex) {
				// something wants us to stop, so exit this thread
				break;
			}
		}
	}
	
	public void askToStop() {
		isRunning = false;
	}
	
	protected abstract void workIt();
}
