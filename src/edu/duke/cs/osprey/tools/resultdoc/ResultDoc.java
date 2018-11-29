package edu.duke.cs.osprey.tools.resultdoc;


import com.sun.org.apache.xml.internal.security.utils.Base64;
import edu.duke.cs.osprey.tools.Streams;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

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

	public ResultDoc(File file) {
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

	public class Image {

		/** binary image */
		public void embed(byte[] data, String type) {
			print("![](data:" + type + ";base64,");
			print(Base64.encode(data)
				// keep everything on one line for URIs
				.replace("\n", "")
			);
			println(")");
		}

		/** text image */
		public void embed(String data, String type) {
			print("![](data:" + type + ";utf8,");
			print(data
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
			);
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


	private final String tableSeparator = " | ";
	/** NOTE: some markdown parsers don't support this table syntax =( */
	public void tableHeader(String ... text) {
		println(Streams.joinToString(text, tableSeparator));
		println(Streams.joinToString(text, tableSeparator, t -> "..."));
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
}
