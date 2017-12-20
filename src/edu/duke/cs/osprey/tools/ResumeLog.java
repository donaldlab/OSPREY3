package edu.duke.cs.osprey.tools;

import edu.duke.cs.osprey.externalMemory.Queue;

import java.io.*;
import java.util.*;
import java.util.function.Consumer;

/**
 * The resume log is designed to help resume long-running operations where
 * a list of things are transformed from "in" things to "out" things.
 * Assuming this transformation is expensive (e.g., energy minimization).
 * the resume log saves a list of the "out" things in the order of the
 * corresponding "in" things. Using the log, the transformation list can be
 * resumed close to the position of the most recently transformed "in" thing,
 * by reading the saved "out" things and matching them to a list of expected
 * "in" things.
 *
 * @param <I> The type of the "in" thing
 * @param <O> The type of the "out" thing
 */
public class ResumeLog<I,O> {

	public static interface Equator<I,O> {
		boolean isEqual(I inThing, O outThing);
	}

	public static interface OutThingSerializer<O> {
		void write(O outThing, DataOutputStream stream) throws IOException;
		O read(DataInputStream stream) throws IOException;
	}

	public final File logFile;
	public final int flushEveryNThings;
	public final Equator<I,I> inEquator;
	public final Equator<I,O> inOutEquator;
	public final OutThingSerializer<O> outSerializer;

	public final Deque<I> inThings = new ArrayDeque<>();
	private List<O> outThings = new ArrayList<>();

	public ResumeLog(File logFile, int flushEveryNThings, Equator<I,I> inEquator, Equator<I,O> inOutEquator, OutThingSerializer<O> outSerializer) {
		this.logFile = logFile;
		this.flushEveryNThings = flushEveryNThings;
		this.inEquator = inEquator;
		this.inOutEquator = inOutEquator;
		this.outSerializer = outSerializer;
	}

	public boolean isLogging() {
		return logFile != null;
	}

	public boolean hasLog() {
		return isLogging() && logFile.exists();
	}

	public void appendInThing(I inThing) {
		inThings.add(inThing);
	}

	private int findIndex(I queryThing) {

		// linear search should be plenty fast enough here
		// the expected list should always be pretty small
		int i = 0;
		for (I inThing : inThings) {
			if (inEquator.isEqual(inThing, queryThing)) {
				return i;
			}
			i++;
		}

		throw new NoSuchElementException("expected " + queryThing + " not found");
	}

	public void matchOutThing(I inThing, O outThing) {

		int i = findIndex(inThing);

		// expand the out buffer to make space if needed
		while (outThings.size() <= i) {
			outThings.add(null);
		}

		// save the out thing in position
		outThings.set(i, outThing);

		// flush buffer as needed
		if (isOutThingsFull()) {
			flushOutThings();
		}
	}

	private boolean isOutThingsFull() {

		// do we have enough in things?
		if (inThings.size() < flushEveryNThings) {
			return false;
		}

		// do we have enough consecutive out things?
		if (outThings.size() < flushEveryNThings) {
			return false;
		}
		for (int i = 0; i< flushEveryNThings; i++) {
			if (outThings.get(i) == null) {
				return false;
			}
		}

		return true;
	}

	private void flushOutThings() {

		try (DataOutputStream out = new DataOutputStream(new FileOutputStream(logFile, true))) {

			for (int i=0; i<outThings.size(); i++) {

				// get the next energy to save
				O outThing = outThings.get(i);
				if (outThing == null) {
					break;
				}

				// pop the corresponding in thing from the queue
				inThings.pop();

				// save it to the file
				outSerializer.write(outThing, out);
			}

		} catch (IOException ex) {
			throw new RuntimeException("Can't write to resume log: " + logFile, ex);
		}
	}

	public void readAll(Queue<I> inThings, Consumer<O> onReadOutThing) {

		try (FileInputStream fin = new FileInputStream(logFile)) {
			DataInputStream in = new DataInputStream(fin);

			while (true) {

				// are we at the end of the file?
				if (fin.getChannel().position()  == fin.getChannel().size()) {

					// end of the log, we're done reading
					break;
				}

				// read the next log conf
				O outThing = outSerializer.read(in);

				// compare to the next in thing
				if (inThings.isEmpty()) {
					throw new IllegalStateException(String.format("resume log expected %s, but didn't find anything.", outThing));
				}
				I inThing = inThings.poll();

				if (!inOutEquator.isEqual(inThing, outThing)) {
					throw new IllegalStateException(String.format("resume log expected %s, but but found %s instead.", outThing, inThing));
				}

				// all is well, pass up the thing
				onReadOutThing.accept(outThing);
			}
		} catch (IOException ex) {
			throw new RuntimeException("Can't read from resume log: " + logFile, ex);
		}
	}

	public void delete() {
		logFile.delete();
	}
}
