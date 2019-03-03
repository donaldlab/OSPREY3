/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.tools;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;


public class Log {

	/** write a fragment to the log (no newline is appended) */
	public static void logf(String format, Object ... args) {
		System.out.print(String.format(format, args));
	}

	/** write a line to the log (a new line is appended) */
	public static void log(String format, Object ... args) {
		System.out.println(String.format(format, args));
	}

	public static void log(StringBuilder buf, String format, Object ... args) {
		buf.append(String.format(format + "\n", args));
	}

	public static String formatBig(BigInteger i) {
		if (i == null) {
			return "null";
		} else if (i.compareTo(BigInteger.valueOf(1000000)) < 0) {
			return String.format("%s", i);
		} else {
			return String.format("%e", i.doubleValue());
		}
	}

	public static String formatBig(BigDecimal f) {
		if (f == null) {
			return "null";
		} else {
			return String.format("%e (%.2f)", f.doubleValue(), MathTools.log10p1(f));
		}
	}

	public static String formatBigLn(BigDecimal f) {
		if (f == null) {
			return "null";
		} else if (MathTools.isZero(f)) {
			return "0";
		} else {
			MathContext mathContext = new MathContext(16, RoundingMode.HALF_UP);
			BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);
			return String.format("%9.4f", bcalc.ln1p(f));
		}
	}

	public static String formatBigLn(MathTools.BigDecimalBounds b) {
		if (b == null) {
			return "null";
		} else {
			return String.format("[%-9s,%9s]",
				formatBigLn(b.lower),
				formatBigLn(b.upper)
			);
		}
	}
}
