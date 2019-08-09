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

import edu.duke.cs.osprey.tools.Streams;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Base64;


public class ResultDoc implements AutoCloseable {

	private final OutputStream out;
	private final Writer writer;

	public ResultDoc(OutputStream out) {
		this.out = out;
		this.writer = new BufferedWriter(new OutputStreamWriter(out));
	}

	@Override
	public void close() {
		try {
			out.close();
		} catch (IOException ex) {
			// don't care
		}
	}

	public ResultDoc() {
		this(System.out);
	}

	private static OutputStream pathToStream(Path path) {
		try {
			return Files.newOutputStream(path);
		} catch (IOException ex) {
			throw new RuntimeException("can't open Path: " + path, ex);
		}
	}

	public ResultDoc(Path path) {
		this(pathToStream(path));
	}

	public ResultDoc(java.io.File file) {
		this(file.toPath());
	}

	public ResultDoc(String path) {
		this(Paths.get(path));
	}

	public void print(String text) {
		try {
			writer.write(text);
			writer.flush();
		} catch (IOException ex) {
			throw new RuntimeException("can't write to ResultDoc", ex);
		}
	}

	public void println() {
		print("\n");
	}

	public void println(String text) {
		print(text);
		println();
	}

	public void h1(String text) {
		print("# ");
		println(text);
	}

	public void h2(String text) {
		print("## ");
		println(text);
	}

	public void h3(String text) {
		print("### ");
		println(text);
	}

	public void h4(String text) {
		print("#### ");
		println(text);
	}

	public void h5(String text) {
		print("##### ");
		println(text);
	}

	public void h6(String text) {
		print("###### ");
		println(text);
	}

	/** unordered list */
	public void ul(String text) {
		print(" * ");
		println(text);
	}

	/** ordered list */
	public void ol(String text) {
		print(" 1. ");
		println(text);
	}

	private static String urlEncodeBinary(byte[] data) {
		return Base64.getEncoder().encodeToString(data)
			// keep everything on one line for URIs
			.replace("\n", "");
	}

	private static String urlEncodeText(String data) {
		return data
			// need to escape some chars to keep from confusing markdown parsers
			.replace("%", "%25")
			.replace("\n", "%0A")
			.replace("\t", "%09")
			.replace(" ", "%20")
			.replace("\"", "%22")
			.replace("(", "%28")
			.replace(")", "%29")
			.replace("#", "%23")
			.replace("|", "%7C")
			.replace("[", "%5B")
			.replace("]", "%5D");
	}

	private static String escapeText(String text) {
		return text
			.replace("(", "\\(")
			.replace(")", "\\)")
			.replace("[", "\\[")
			.replace("]", "\\]");
	}


	public class Image {

		/** binary image */
		public void embed(byte[] data, String type) {
			print("![](data:");
			print(type);
			print(";base64,");
			print(urlEncodeBinary(data));
			println(")");
		}

		/** text image */
		public void embed(String data, String type) {
			print("![](data:");
			print(type);
			print(";utf8,");
			print(urlEncodeText(data));
			println(")");
		}

		public void png(byte[] data) {
			embed(data, "image/png");
		}

		public void svg(String data) {
			embed(data, "image/svg+xml");
		}
	}
	public final Image image = new Image();


	public class File {

		/** binary file */
		public void embed(byte[] data, String type, String link) {
			print("[");
			print(escapeText(link));
			print("](data:");
			print(type);
			print(";base64,");
			print(urlEncodeBinary(data));
			println(")");
		}

		/** text file */
		public void embed(String data, String type, String link) {
			print("[");
			print(escapeText(link));
			print("](data:");
			print(type);
			print(";utf8,");
			print(urlEncodeText(data));
			println(")");
		}

		public void csv(String data, String link) {
			embed(data, "text/csv", link);
		}

		public void tsv(String data, String link) {
			embed(data, "text/tab-separated-values", link);
		}
	}
	public final File file = new File();


	private final String tableSeparator = " | ";

	public void tableHeader(String ... text) {
		println(Streams.joinToString(text, tableSeparator));
		println(Streams.joinToString(text, tableSeparator, t -> "---"));
	}

	public void tableRow(String ... text) {
		println(Streams.joinToString(text, tableSeparator));
	}

	/** horizontal rule */
	public void hr() {
		println("---");
	}

	public static enum PlotType {

		PNG {
			@Override
			public void plot(ResultDoc doc, Plot plot) {
				doc.image.png(plot.renderPng());
			}
		},
		SVG {
			@Override
			public void plot(ResultDoc doc, Plot plot) {
				doc.image.svg(plot.renderSvg());
			}
		};

		abstract void plot(ResultDoc doc, Plot plot);
	}

	public PlotType defaultPlotType = PlotType.PNG;

	public void plot(Plot plot, PlotType plotType) {
		plotType.plot(this, plot);
	}

	public void plot(Plot plot) {
		plot(plot, defaultPlotType);
	}

	// TODO: support LaTeX via MathJax?
}
