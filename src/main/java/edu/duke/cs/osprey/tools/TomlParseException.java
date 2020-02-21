package edu.duke.cs.osprey.tools;

import org.tomlj.TomlPosition;


public class TomlParseException extends RuntimeException {

	public final TomlPosition pos;

	public TomlParseException(String msg, TomlPosition pos) {
		super(msg + (pos != null ? " at " + pos.toString() : ""));
		this.pos = pos;
	}

	public TomlParseException(String msg) {
		this(msg, null);
	}
}
