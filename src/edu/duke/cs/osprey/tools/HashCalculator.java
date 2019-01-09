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

public class HashCalculator {
	
	public static int combineHashes(long ... nums) {
		long hashCode = 1;
		for (long i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return (int)((hashCode >> 32) ^ hashCode);
	}
	
	public static int combineHashes(int ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return hashCode;
	}
	
	public static int combineHashes(short ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return hashCode;
	}
	
	public static int combineHashes(byte ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return hashCode;
	}
	
	public static int combineHashesCommutative(int ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode += i;
		}
		return hashCode;
	}
	
	public static int hashIds(int ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 37 ^ (i + 1);
		}
		return hashCode;
	}

	@SafeVarargs
	public static <T> int combineObjHashes(T ... objs) {
		int hashCode = 1;
		for (T obj : objs) {
			hashCode = hashCode * 31 + obj.hashCode();
		}
		return hashCode;
	}
}
