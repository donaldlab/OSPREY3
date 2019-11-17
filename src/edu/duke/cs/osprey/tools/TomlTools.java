package edu.duke.cs.osprey.tools;


import org.tomlj.*;


public class TomlTools {

	public static TomlTable parseOrThrow(String toml) {
		TomlParseResult doc = Toml.parse(toml);
		if (doc.hasErrors()) {
			throw new TomlParseException("TOML parsing failure:\n"
				+ Streams.joinToString(doc.errors(), "\n", err -> err.toString())
			);
		}
		return doc;
	}

	private static String formatPos(TomlPosition pos) {
		if (pos != null) {
			return String.format(" @ %s", pos);
		} else {
			return "";
		}
	}

	private static Integer longToIntOrNull(Long longval) {
		int intval = longval.intValue();
		// check for overflow
		if (intval != longval) {
			return null;
		}
		return intval;
	}


	// TABLE TOOLS

	public static Object getOrThrow(TomlTable table, String key) {
		return getOrThrow(table, key, null);
	}
	public static Object getOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = table.get(key);
		if (value == null) {
			throw new TomlParseException(String.format("missing %s in table%s",
				key,
				formatPos(tablePos)
			));
		}
		return value;
	}

	private static TomlParseException invalidTableType(String key, String desiredType, Object value, TomlPosition tablePos) {
		return new TomlParseException(String.format("%s in table%s should be a %s, not a(n) %s",
			key,
			desiredType,
			value.getClass().getSimpleName(),
			formatPos(tablePos)
		));
	}

	public static TomlTable getTableOrThrow(TomlTable table, String key) {
		return getTableOrThrow(table, key, null);
	}
	public static TomlTable getTableOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = getOrThrow(table, key, tablePos);
		if (value instanceof TomlTable) {
			return (TomlTable)value;
		}
		throw invalidTableType(key, "Table", value, tablePos);
	}

	public static TomlArray getArrayOrThrow(TomlTable table, String key) {
		return getArrayOrThrow(table, key, null);
	}
	public static TomlArray getArrayOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = getOrThrow(table, key, tablePos);
		if (value instanceof TomlArray) {
			return (TomlArray)value;
		}
		throw invalidTableType(key, "Array", value, tablePos);
	}

	public static String getStringOrThrow(TomlTable table, String key) {
		return getStringOrThrow(table, key, null);
	}
	public static String getStringOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = getOrThrow(table, key, tablePos);
		if (value instanceof String) {
			return (String)value;
		}
		throw invalidTableType(key, "String", value, tablePos);
	}

	public static long getLongOrThrow(TomlTable table, String key) {
		return getLongOrThrow(table, key, null);
	}
	public static long getLongOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = getOrThrow(table, key, tablePos);
		if (value instanceof Long) {
			return (Long)value;
		}
		throw invalidTableType(key, "Long", value, tablePos);
	}

	public static int getIntOrThrow(TomlTable table, String key) {
		return getIntOrThrow(table, key, null);
	}
	public static int getIntOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = getOrThrow(table, key, tablePos);
		if (value instanceof Long) {
			Integer intval = longToIntOrNull((Long)value);
			if (intval == null) {
				throw new TomlParseException(String.format("%s = %d in table%s overflows int value",
					key,
					value,
					tablePos
				));
			}
			return intval;
		}
		throw invalidTableType(key, "Int", value, tablePos);
	}

	public static double getDoubleOrThrow(TomlTable table, String key) {
		return getDoubleOrThrow(table, key, null);
	}
	public static double getDoubleOrThrow(TomlTable table, String key, TomlPosition tablePos) {
		Object value = getOrThrow(table, key, tablePos);
		if (value instanceof Double) {
			return (Double)value;
		}
		throw invalidTableType(key, "Double", value, tablePos);
	}


	// ARRAY TOOLS

	public static Object getOrThrow(TomlArray array, int i) {
		return getOrThrow(array, i, null);
	}
	public static Object getOrThrow(TomlArray array, int i, TomlPosition tablePos) {
		Object value = null;
		if (i >= 0 && i < array.size()) {
			value = array.get(i);
		}
		if (value == null) {
			throw new TomlParseException(String.format("missing index %d in array%s",
				i,
				formatPos(tablePos)
			));
		}
		return value;
	}

	private static TomlParseException invalidArrayType(int i, String desiredType, Object value, TomlPosition iPos) {
		return new TomlParseException(String.format("index %d%s in array should be a %s, not a(n) %s",
			i,
			formatPos(iPos),
			desiredType,
			value.getClass().getSimpleName()
		));
	}

	public static TomlTable getTableOrThrow(TomlArray array, int i) {
		return getTableOrThrow(array, i, null);
	}
	public static TomlTable getTableOrThrow(TomlArray array, int i, TomlPosition arrayPos) {
		Object value = getOrThrow(array, i, arrayPos);
		if (value instanceof TomlTable) {
			return (TomlTable)value;
		}
		throw invalidArrayType(i, "Table", value, array.inputPositionOf(i));
	}

	public static TomlArray getArrayOrThrow(TomlArray array, int i) {
		return getArrayOrThrow(array, i, null);
	}
	public static TomlArray getArrayOrThrow(TomlArray array, int i, TomlPosition arrayPos) {
		Object value = getOrThrow(array, i, arrayPos);
		if (value instanceof TomlArray) {
			return (TomlArray)value;
		}
		throw invalidArrayType(i, "Array", value, array.inputPositionOf(i));
	}

	public static String getStringOrThrow(TomlArray array, int i) {
		return getStringOrThrow(array, i, null);
	}
	public static String getStringOrThrow(TomlArray array, int i, TomlPosition arrayPos) {
		Object value = getOrThrow(array, i, arrayPos);
		if (value instanceof String) {
			return (String)value;
		}
		throw invalidArrayType(i, "String", value, array.inputPositionOf(i));
	}

	public static long getLongOrThrow(TomlArray array, int i) {
		return getLongOrThrow(array, i, null);
	}
	public static long getLongOrThrow(TomlArray array, int i, TomlPosition arrayPos) {
		Object value = getOrThrow(array, i, arrayPos);
		if (value instanceof Long) {
			return (Long)value;
		}
		throw invalidArrayType(i, "Long", value, array.inputPositionOf(i));
	}

	public static int getIntOrThrow(TomlArray array, int i) {
		return getIntOrThrow(array, i, null);
	}
	public static int getIntOrThrow(TomlArray array, int i, TomlPosition arrayPos) {
		Object value = getOrThrow(array, i, arrayPos);
		if (value instanceof Long) {
			Integer intval = longToIntOrNull((Long)value);
			if (intval == null) {
				throw new TomlParseException(String.format("index %d%s in array = %d overflows int value",
					i,
					formatPos(array.inputPositionOf(i)),
					value
				));
			}
			return intval;
		}
		throw invalidArrayType(i, "Int", value, array.inputPositionOf(i));
	}

	public static double getDoubleOrThrow(TomlArray array, int i) {
		return getDoubleOrThrow(array, i, null);
	}
	public static double getDoubleOrThrow(TomlArray array, int i, TomlPosition arrayPos) {
		Object value = getOrThrow(array, i, arrayPos);
		if (value instanceof Double) {
			return (Double)value;
		}
		throw invalidArrayType(i, "Double", value, array.inputPositionOf(i));
	}
}
