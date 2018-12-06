package edu.duke.cs.osprey.sofea;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.tools.MathTools;
import org.junit.Test;

import java.io.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

public class TestBigDecimalIO {

	private static class Buf {

		ByteArrayOutputStream buf = new ByteArrayOutputStream();
		DataOutput out = new DataOutputStream(buf);
		DataInput in = null;

		int size() {
			return buf.size();
		}

		void swap() {
			in = new DataInputStream(new ByteArrayInputStream(buf.toByteArray()));
		}
	}

	@Test
	public void variable()
	throws IOException {

		BigDecimalIO io = new BigDecimalIO.Variable();

		Buf buf = new Buf();
		io.write(buf.out, BigDecimal.ONE);
		io.write(buf.out, BigDecimal.ZERO);
		io.write(buf.out, new BigDecimal("2.0"));
		io.write(buf.out, new BigDecimal("4.2"));
		io.write(buf.out, new BigDecimal("1.94538e-4893"));
		io.write(buf.out, MathTools.BigNaN);
		io.write(buf.out, MathTools.BigNegativeInfinity);
		io.write(buf.out, MathTools.BigPositiveInfinity);

		buf.swap();

		assertThat(io.read(buf.in), is(BigDecimal.ONE));
		assertThat(io.read(buf.in), is(BigDecimal.ZERO));
		assertThat(io.read(buf.in), is(new BigDecimal("2.0")));
		assertThat(io.read(buf.in), is(new BigDecimal("4.2")));
		assertThat(io.read(buf.in), is(new BigDecimal("1.94538e-4893")));
		assertThat(io.read(buf.in), is(MathTools.BigNaN));
		assertThat(io.read(buf.in), is(MathTools.BigNegativeInfinity));
		assertThat(io.read(buf.in), is(MathTools.BigPositiveInfinity));
	}

	@Test
	public void fixed()
		throws IOException {

		BigDecimalIO.Fixed io = new BigDecimalIO.Fixed(new MathContext(16, RoundingMode.HALF_UP));

		Buf buf = new Buf();
		io.write(buf.out, BigDecimal.ONE);
		io.write(buf.out, BigDecimal.ZERO);
		io.write(buf.out, new BigDecimal("2.0"));
		io.write(buf.out, new BigDecimal("4.2"));
		io.write(buf.out, new BigDecimal("1.94538e-4893"));
		io.write(buf.out, MathTools.BigNaN);
		io.write(buf.out, MathTools.BigNegativeInfinity);
		io.write(buf.out, MathTools.BigPositiveInfinity);

		assertThat(buf.size(), is(io.numBytes*8));

		buf.swap();

		assertThat(io.read(buf.in), is(BigDecimal.ONE));
		assertThat(io.read(buf.in), is(BigDecimal.ZERO));
		assertThat(io.read(buf.in), is(new BigDecimal("2.0")));
		assertThat(io.read(buf.in), is(new BigDecimal("4.2")));
		assertThat(io.read(buf.in), is(new BigDecimal("1.94538e-4893")));
		assertThat(io.read(buf.in), is(MathTools.BigNaN));
		assertThat(io.read(buf.in), is(MathTools.BigNegativeInfinity));
		assertThat(io.read(buf.in), is(MathTools.BigPositiveInfinity));
	}
}
