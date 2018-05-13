/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.tools;

import java.math.BigDecimal;
import java.math.BigInteger;

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
