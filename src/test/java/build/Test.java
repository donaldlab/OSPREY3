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
}
