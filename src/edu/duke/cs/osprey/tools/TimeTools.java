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

public class TimeTools {

	private static final long baseMs = System.currentTimeMillis();
	private static final long baseNs = System.nanoTime();

	/**
	 * Returns the number of milliseconds elapsed since midnight, January 1, 1970 UTC.
	 *
	 * Won't overflow a signed long until about the year 292473178. Hopefully by then we'll still a have a sun.
	 */
	public static long getTimestampMs() {
		return System.currentTimeMillis();
	}

	/**
	 * Returns approximately the number of microseconds elapsed since midnight, January 1, 1970 UTC.
	 *
	 * Won't overflow a signed long until about the year 294441. Hopefully by then we won't be using this code anymore.
	 */
	public static long getTimestampUs() {
		return baseMs*1000L + (System.nanoTime() - baseNs)/1000;
	}

	/**
	 * Returns approximately the number of nanoseconds elapsed since midnight, January 1, 1970 UTC.
	 *
	 * Won't overflow a signed long until about the year 2262. Hopefully by then we'll have better computers.
	 * And flying cars.
	 */
	public static long getTimestampNs() {
		return baseMs*1000000L + System.nanoTime() - baseNs;
	}
}
