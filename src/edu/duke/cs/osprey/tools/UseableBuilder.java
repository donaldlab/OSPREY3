package edu.duke.cs.osprey.tools;

public interface UseableBuilder<T extends AutoCleanable> {
	
	public static interface Block<T> {
		void run(T thing) throws Exception;
	}
	
	T build();
	
	default void use(Block<T> block) {
		try (T thing = build()) {
			block.run(thing);
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
	}
}
