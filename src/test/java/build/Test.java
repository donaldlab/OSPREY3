package build;

import java.util.List;
import java.util.Map;

/**
 * Javadoc comment for Test
 */
public class Test {

	/** javadoc for i */
	int i = 1;

	/** javadoc for getRandonNumber */
	int getRandomNumber() {
		// chosen by fair die roll, guaranteed to be random
		return 4;
	}

	/** javadoc for stuff */
	void stuff(int a, float b) {}

	/** high-precision stuff */
	void stuff(double c) {}

	/** javadoc for Foo */
	class Foo {

		/** javadoc for j */
		float j = 42f;

		/** javadoc for bars */
		List<Bar> bars;

		/** javadoc for barIndex */
		Map<String,Bar> barIndex;

		/** javadoc for Bar */
		class Bar {

			/** javadoc for k */
			String k = "cheese";
		}
	}

	/** javadoc for container class */
	class Container<T> {}

	/** javadoc for container field */
	Container<String> container;

	/**
	 * This method has params
	 * @param a the a
	 * @param b the b
	 */
	void params(int a, int b) {}

	/**
	 * This javadoc has links!
	 * {@link #container}
	 * {@link List}
	 * {@link Foo a foo}
	 * {@link #stuff}
	 * {@link #stuff(double)}
	 */
	int j = 5;

	/**
	 * Here are some non-standard javadoc tags
	 * {@note This is a note.}
	 * {@warning This is a warning.}
	 * {@cite This is a citation.
	 * It can have multiple lines.
	 * Sometimes even three.}
	 */
	int k = 7;

	/** hello {@note hi} world */
	int m = 9;

	/**
	 * This one has block comments.
	 *
	 * @note the note
	 */
	int blocks;

	/** This is an eenumm! */
	enum Eenumm {
		Val1,
		Val2
	}
}
