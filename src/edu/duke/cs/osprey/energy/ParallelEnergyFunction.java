package edu.duke.cs.osprey.energy;

import java.util.ArrayList;

public class ParallelEnergyFunction implements EnergyFunction {
	
	
	private static class Sync {
		
		private Object workSync;
		private Object resultsSync;
		private volatile int counter;
		
		public Sync() {
			workSync = new Object();
			resultsSync = new Object();
		}
	
		public void waitForWork(int timeoutMs)
		throws InterruptedException {
			synchronized (workSync) {
				workSync.wait(timeoutMs);
			}
		}
		
		public void sendWork() {
			counter = 0;
			for (int i=0; i<processors.length; i++) {
				processors[i].hasWork = true;
			}
			synchronized (workSync) {
				workSync.notifyAll();
			}
		}
		
		public void waitForResults(int timeoutMs)
		throws InterruptedException {
			synchronized (resultsSync) {
				resultsSync.wait(timeoutMs);
			}
		}
		
		public void finishedWork() {
			int c;
			synchronized (this) {
				counter++;
				c = counter;
			}
			if (c == processors.length) {
				synchronized (resultsSync) {
					resultsSync.notify();
				}
			}
		}
	};
	
	private static class Processor extends Thread {
		
		// this flag is hit from multiple threads concurrently, so make it volatile
		private volatile boolean isRunning;
		
		private int startIndex;
		private int stopIndex;
		private boolean hasWork;
		private double energy;
		
		public Processor(int index) {
			super("EnergyProcessor-" + index);
			setDaemon(true);
			startIndex = 0;
			stopIndex = 0;
			hasWork = false;
			energy = 0;
		}
		
		@Override
		public void run() {
			
			double energy;
			ArrayList<EnergyFunction> terms;
			ArrayList<Double> coeffs;
			
			isRunning = true;
			while (isRunning) {
				
				// is there any work to do?
				if (hasWork) {
					
					// copy some references/values to the stack for that little extra bit of speed
					terms = efunc.terms;
					coeffs = efunc.coeffs;
					
					// do it!
					energy = 0;
					for (int i=startIndex; i<=stopIndex; i++) {
						energy += terms.get(i).getEnergy()*coeffs.get(i);
					}
					
					// save results to 'this' memory
					this.energy = energy;
					hasWork = false;
					
					sync.finishedWork();
					
					// NOTE: after syncing, someone could have given us more work, so we have to check again
					if (hasWork) {
						continue;
					}
				}
				
				// wait until we get more work
				// but check the isRunning flag every second or so
				try {
					sync.waitForWork(1000);
				} catch (InterruptedException ex) {
					// something wants us to stop, so exit this thread
					break;
				}
			}
		}
		
		public void askToStop() {
			isRunning = false;
		}
	}
	
	private static Sync sync = new Sync();
	private static Processor[] processors;
	private static ParallelEnergyFunction efunc;
	
	static {
		processors = null;
		efunc = null;
	}
	
	public static boolean areProcessorsStarted() {
		return processors != null;
	}
	
	public static void startProcessors(int numThreads) {
		processors = new Processor[numThreads];
		for (int i=0; i<numThreads; i++) {
			processors[i] = new Processor(i);
			processors[i].start();
		}
	}
	
	public static void stopProcessors() {
		if (processors == null) {
			return;
		}
		for (int i=0; i<processors.length; i++) {
			processors[i].askToStop();
		}
		processors = null;
	}
	
	public static void setEFunc(ParallelEnergyFunction val) {
		
		if (efunc == val) {
			return;
		}
		
		efunc = val;
		
		// set up partition
		// partition terms among processors
		int numTerms = efunc.terms.size();
		int numProcs = processors.length;
		int width = (numTerms + numProcs - 1)/numProcs;
		int startIndex = 0;
		int stopIndex = startIndex + width - 1;
		for (Processor proc : processors) {
			proc.startIndex = startIndex;
			proc.stopIndex = stopIndex;
			startIndex += width;
			stopIndex = Math.min(stopIndex + width, numTerms - 1);
		}
	}

	private static final long serialVersionUID = -2789380428939629566L;
	
	ArrayList<EnergyFunction> terms;
	ArrayList<Double> coeffs;
	
	public ParallelEnergyFunction(ArrayList<EnergyFunction> terms, ArrayList<Double> coeffs) {
		this.terms = terms;
		this.coeffs = coeffs;
	}
	
	public ArrayList<EnergyFunction> getTerms() {
		return terms;
	}

	public ArrayList<Double> getCoeffs() {
		return coeffs;
	}
	
	@Override
	public double getEnergy() {
		
		setEFunc(this);
		
		// start all the processors
		sync.sendWork();
		
		// wait for the work to finish
		try {
			sync.waitForResults(10000);
			if (sync.counter != processors.length) {
				throw new Error("Timed our waiting 10 seconds for energy calculations to finish!"
					+ "\nEnergy calculation shouldn't take more than 10 seconds, right?");
			}
		} catch (InterruptedException ex) {
			// something wanted us to stop, so stop, then forward the exception
			throw new Error(ex);
		}
		
		// collect the results
		double energy = 0;
		for (int i=0; i<processors.length; i++) {
			energy += processors[i].energy;
		}
		return energy;
	}
}
