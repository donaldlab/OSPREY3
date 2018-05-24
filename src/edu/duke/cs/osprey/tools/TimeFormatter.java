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

import java.util.concurrent.TimeUnit;

public class TimeFormatter {
	
	public static final long NSpUS = 1000;
	public static final long NSpMS = NSpUS * 1000;
	public static final long NSpS = NSpMS * 1000;
	public static final long NSpM = NSpS * 60;
	public static final long NSpH = NSpM * 60;

	public static double getTimeUs(long ns) {
		return (double)ns/NSpUS;
	}
	
	public static double getTimeMs(long ns) {
		return (double)ns/NSpMS;
	}
	
	public static double getTimeS(long ns) {
		return (double)ns/NSpS;
	}
	
	public static double getTimeM(long ns) {
		return (double)ns/NSpM;
	}
	
	public static double getTimeH(long ns) {
		return (double)ns/NSpH;
	}
	
	public static String format(long ns) {
		return format(ns, 0);
	}
	
	public static String format(long ns, int decimals) {
		long nsabs = Math.abs(ns);
		if (getTimeH(nsabs) > 1) {
			return format(ns, TimeUnit.HOURS, decimals);
		} else if (getTimeM(nsabs) > 1) {
			return format(ns, TimeUnit.MINUTES, decimals);
		} else if (getTimeS(nsabs) > 1) {
			return format(ns, TimeUnit.SECONDS, decimals);
		} else if (getTimeMs(nsabs) > 1) {
			return format(ns, TimeUnit.MILLISECONDS, decimals);
		} else if (getTimeUs(nsabs) > 1) {
			return format(ns, TimeUnit.MICROSECONDS, decimals);
		} else {
			return format(ns, TimeUnit.NANOSECONDS, decimals);
		}
	}
	
	public static String format(long ns, TimeUnit unit) {
		return format(ns, unit, 0);
	}
	
	public static String format(long ns, TimeUnit unit, int decimals) {
		
		String formatSpec;
		if (unit == TimeUnit.NANOSECONDS) {
			formatSpec = "%d ";
		} else {
			formatSpec = "%." + decimals + "f ";
		}
		
		switch (unit) {
			case NANOSECONDS: return String.format(formatSpec + "ns", ns);
			case MICROSECONDS: return String.format(formatSpec + "us", getTimeUs(ns));
			case MILLISECONDS: return String.format(formatSpec + "ms", getTimeMs(ns));
			case SECONDS: return String.format(formatSpec + "s", getTimeS(ns));
			case MINUTES: return String.format(formatSpec + "m", getTimeM(ns));
			case HOURS: return String.format(formatSpec + "h", getTimeH(ns));
			default:
				throw new Error("Unsupported time unit: " + unit);
		}
	}
}
