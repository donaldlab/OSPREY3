package edu.duke.cs.osprey.tools;

import java.math.BigDecimal;
import java.math.BigInteger;

public class Log {

	public static void log(String format, Object ... args) {
		System.out.println(String.format(format, args));
	}

	public static void log(StringBuilder buf, String format, Object ... args) {
		buf.append(String.format(format + "\n", args));
	}

	public static String formatBig(BigInteger i) {
		if (i.compareTo(BigInteger.valueOf(1000000)) < 0) {
			return String.format("%s", i);
		} else {
			return String.format("%e", i.doubleValue());
		}
	}

	public static String formatBig(BigDecimal f) {
		return String.format("%e (%.2f)", f.doubleValue(), MathTools.log10p1(f));
	}
}
