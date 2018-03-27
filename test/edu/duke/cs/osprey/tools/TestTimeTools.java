package edu.duke.cs.osprey.tools;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

public class TestTimeTools {

	@Test
	public void test()
	throws Exception {

		long ms = TimeTools.getTimestampMs();
		long us = TimeTools.getTimestampUs();
		long ns = TimeTools.getTimestampNs();

		assertThat(us/1000L - ms, lessThanOrEqualTo(1L));
		assertThat(ns/1000000L - ms, lessThanOrEqualTo(1L));

		Thread.sleep(37);

		ms = TimeTools.getTimestampMs();
		us = TimeTools.getTimestampUs();
		ns = TimeTools.getTimestampNs();

		assertThat(us/1000L - ms, lessThanOrEqualTo(1L));
		assertThat(ns/1000000L - ms, lessThanOrEqualTo(1L));

		Thread.sleep(379);

		ms = TimeTools.getTimestampMs();
		us = TimeTools.getTimestampUs();
		ns = TimeTools.getTimestampNs();

		assertThat(us/1000L - ms, lessThanOrEqualTo(1L));
		assertThat(ns/1000000L - ms, lessThanOrEqualTo(1L));
	}
}
