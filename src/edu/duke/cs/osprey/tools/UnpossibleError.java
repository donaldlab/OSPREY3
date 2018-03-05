package edu.duke.cs.osprey.tools;

public class UnpossibleError extends Error {

	public UnpossibleError() {
		super("A programmer thought this code path was un-possible, even though the compiler didn't."
			+ "\nLooks like the compiler was right. This is definitely a bug.");
	}
}
