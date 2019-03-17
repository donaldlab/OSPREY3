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

package edu.duke.cs.osprey.tools.resultdoc;


import com.google.common.base.Charsets;
import edu.duke.cs.osprey.tools.Streams;
import org.apache.commons.io.IOUtils;

import java.io.*;
import java.util.Map;


public class GNUPlot {

	public static class Cmd {

		private final BufferedWriter out;

		private Cmd(BufferedWriter out) {
			this.out = out;
		}

		public void command(String text)
		throws IOException {
			out.write(text);
			out.write("\n");
		}

		public void command(String format, Object ... args)
		throws IOException {
			command(String.format(format, args));
		}

		public void flush()
		throws IOException {
			out.flush();
		}

		public void datablock(String name, String data)
		throws IOException {
			command("$" + name + " << EOD");
			command(data);
			command("EOD");
		}

		/** escape text values and box with " chars */
		public static String renderText(String val) {
			return "\""
				+ val
					.replace("\"", "\\\"")
					.replace("\n", "\\n")
				+ "\"";
		}

		private abstract class Prop {

			private final String set;
			private final String unset;

			public Prop(String set, String unset) {
				this.set = set;
				this.unset = unset;
			}

			public Prop(String name) {
				this("set " + name + " %s", "unset " + name);
			}

			protected void set(String val)
			throws IOException {
				if (val == null) {
					command(unset);
				} else {
					command(String.format(set, val));
				}
			}
		}

		public class NumberProp extends Prop {

			private Number val = null;

			public NumberProp(String set, String unset) {
				super(set, unset);
			}

			public NumberProp(String name) {
				super(name);
			}

			public void set(Number val)
			throws IOException {
				this.val = val;
				if (val != null) {
					super.set(val.toString());
				} else {
					super.set(null);
				}
			}

			public Number get() {
				return val;
			}
		}

		public class LiteralProp extends Prop {

			private String val = null;

			public LiteralProp(String set, String unset) {
				super(set, unset);
			}

			public LiteralProp(String name) {
				super(name);
			}

			public void set(String val)
			throws IOException {
				this.val = val;
				super.set(val);
			}

			public void set(String format, Object ... args)
			throws IOException {
				set(String.format(format, args));
			}

			public String get() {
				return val;
			}
		}

		public class TextProp extends Prop {

			private String val = null;

			public TextProp(String set, String unset) {
				super(set, unset);
			}

			public TextProp(String name) {
				super(name);
			}

			public void set(String val)
			throws IOException {
				this.val = val;
				if (val != null) {
					val = renderText(val);
				}
				super.set(val);
			}

			public void set(String format, Object ... args)
			throws IOException {
				set(String.format(format, args));
			}

			public String get() {
				return val;
			}
		}

		public class Tics {

			private final String name;

			public final TextProp format;
			public final NumberProp scale;
			public final NumberProp rotate;

			public Tics(String name) {

				this.name = name;

				format = new TextProp(name + " format");
				scale = new NumberProp(name + " scale");
				rotate = new NumberProp("set " + name + " rotate by %s", "set " + name + " norotate");
			}

			/** hide the axis entirely **/
			public void hide()
			throws IOException {
				command("unset " + name);
			}

			public void set(Map<Number,String> ticks)
			throws IOException {
				// eg, set xtics ("low" 0, "medium" 50, "high" 100)
				command("set " + name + " ("
					+ Streams.joinToString(ticks, ", ", (val, label) -> renderText(label) + " " + val)
					+ ")"
				);
			}

			public void add(Number val, String label)
			throws IOException {
				command("set " + name + " add (" + renderText(label) + " " + val + ")");
			}
		}

		public final Tics xtics = new Tics("xtics");
		public final Tics ytics = new Tics("ytics");

		public class Style {
			public final LiteralProp data = new LiteralProp("style data");
		}
		public final Style style = new Style();

		public final LiteralProp terminal = new LiteralProp("terminal");
		public final TextProp xlabel = new TextProp("xlabel");
		public final TextProp ylabel = new TextProp("ylabel");
		public final LiteralProp key = new LiteralProp("key");
	}

	public static interface Plotter {
		void plot(Cmd cmd) throws IOException;
	}

	public static byte[] plotPng(int width, int height, Plotter plotter) {

		// wrap the plotter with SVG commands
		Plotter pngPlotter = c -> {
			c.terminal.set("png size %s, %s", width, height);
			plotter.plot(c);
			c.flush();
		};

		// plot it and grab the binary result
		try (InputStream result = plot(pngPlotter)) {
			return IOUtils.toByteArray(result);
		} catch (IOException ex) {
			throw new RuntimeException("can't read PNG plot", ex);
		}
	}

	public static String plotSvg(int width, int height, Plotter plotter) {

		// wrap the plotter with SVG commands
		Plotter svgPlotter = c -> {
			c.terminal.set("svg size %s, %s", width, height);
			plotter.plot(c);
			c.flush();
		};

		// plot it and grab the string result
		try (InputStream result = plot(svgPlotter)) {
			return IOUtils.toString(result, Charsets.UTF_8);
		} catch (IOException ex) {
			throw new RuntimeException("can't read SVG plot", ex);
		}
	}

	private static InputStream plot(Plotter plotter) {
		try {

			// start the GNUPlot process
			Process process = new ProcessBuilder().command("gnuplot").start();

			try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()))) {
				plotter.plot(new Cmd(writer));
			}

			// dump any errors from stderr
			for (String line : IOUtils.readLines(process.getErrorStream(), Charsets.UTF_8)) {
				System.err.println("GNUPlot: " + line);
			}

			// read GNUPlot response from stdout
			return process.getInputStream();

		} catch (IOException ex) {
			throw new RuntimeException("can't run GNUPlot", ex);
		}
	}
}
