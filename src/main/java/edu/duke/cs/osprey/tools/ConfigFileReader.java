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

import java.util.Iterator;
import java.util.function.BiConsumer;


/**
 * tool for reading text files that have sections
 *
 * section headers are lines that look like:
 * [section name]
 */
public class ConfigFileReader {

	private Iterator<String> iter;
	private String line;

	public ConfigFileReader(String text) {
		iter =  FileTools.parseLines(text).iterator();
	}

	private String normalizeLine(String line) {

		// trim off comments
		int commentPos = line.indexOf('#');
		if (commentPos >= 0) {
			line = line.substring(0, commentPos);
		}

		// trim off surrounding whitespace
		line = line.trim();

		return line;
	}

	/** returns the current line, or null if end of file */
	public String getLine() {
		return line;
	}

	/** returns the section name of the current line, if any */
	public String getSectionName() {
		if (line != null && line.startsWith("[") && line.endsWith("]")) {
			return line.substring(1, line.length() - 1);
		}
		return null;
	}

	/** parses a line like Name=Value and returns Name and Value in a callback */
	public void getAssignment(BiConsumer<String,String> callback) {

		String[] parts = line.split("=");
		for (int i=0; i<parts.length; i++) {
			parts[i] = parts[i].trim();
		}

		if (parts.length < 2) {
			throw new IllegalArgumentException("line is not an assignment: " + line);
		}

		callback.accept(parts[0], parts[1]);
	}

	/** set the current line to the next line */
	public void advance() {
		if (iter.hasNext()) {
			line = normalizeLine(iter.next());
		} else {
			line = null;
		}
	}

	/** set the current line to the next non-empty line */
	public void advanceToNonEmptyLine() {
		while (iter.hasNext()) {
			line = normalizeLine(iter.next());
			if (line.isEmpty()) {
				line = null;
			} else {
				break;
			}
		}
	}
}
