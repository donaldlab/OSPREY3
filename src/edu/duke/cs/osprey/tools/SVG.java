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

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.dom.GenericDocument;
import org.apache.batik.svggen.DOMGroupManager;
import org.apache.batik.svggen.SVGGeneratorContext;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.svggen.SVGSyntax;
import org.apache.batik.util.CSSConstants;
import org.apache.batik.util.SVGConstants;
import org.w3c.dom.*;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.*;
import java.lang.reflect.Field;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class SVG {

	public static enum LengthUnit {

		em("em"),
		ex("ex"),
		px("px"),
		in("in"),
		cm("cm"),
		mm("mm"),
		pt("pt"),
		pc("pc"),
		percent("%");

		public final String token;

		private LengthUnit(String token) {
			this.token = token;
		}

		public String render(double value) {
			return String.format("%f%s", value, token);
		}
	}

	public static class StyleClass {

		public final String name;
		public final Map<String,String> values = new HashMap<>();

		public int priority = 0;

		public StyleClass(String name) {
			this.name = name;
		}

		private String colorToString(int color) {
			return String.format("#%6s", Integer.toHexString(color)).replace(' ', '0');
		}

		public void set(String prop, String value) {
			values.put(prop, value);
		}

		public void setNoStroke() {
			set(CSSConstants.CSS_STROKE_PROPERTY, "none");
		}

		public void setStrokeColor(int color) {
			set(CSSConstants.CSS_STROKE_PROPERTY, colorToString(color));
		}

		public void setStrokeWidth(double size) {
			set(CSSConstants.CSS_STROKE_WIDTH_PROPERTY, Double.toString(size));
		}

		public void setStrokeOpacity(double opacity) {
			set(CSSConstants.CSS_STROKE_OPACITY_PROPERTY, Double.toString(opacity));
		}

		public void setStrokeDashArray(int ... array) {
			String arrayString = String.join(",",
				Arrays.stream(array)
					.mapToObj((i) -> Integer.toString(i))
					.collect(Collectors.toList())
			);
			set(CSSConstants.CSS_STROKE_DASHARRAY_PROPERTY, arrayString);
		}

		public void setNoFill() {
			set(CSSConstants.CSS_FILL_PROPERTY, "none");
		}

		public void setFillColor(int color) {
			set(CSSConstants.CSS_FILL_PROPERTY, colorToString(color));
		}

		public void setFontSize(double size, LengthUnit unit) {
			set(CSSConstants.CSS_FONT_SIZE_PROPERTY, unit.render(size));
		}

		public static enum TextAnchor {
			Start,
			Middle,
			End
		}

		public void setTextAnchor(TextAnchor val) {
			set(CSSConstants.CSS_TEXT_ANCHOR_PROPERTY, val.name().toLowerCase());
		}

		public static enum DominantBaseline {

			Auto("auto"),
			UseScript("use-script"),
			NoChange("no-change"),
			ResetSize("reset-size"),
			Ideographic("ideographic"),
			Alphabetic("alphabetic"),
			Hanging("hanging"),
			Mathematical("mathematical"),
			Central("central"),
			Middle("middle"),
			TextAfterEdge("text-after-edge"),
			TextBeforeEdge("text-before-edge"),
			Inherit("inherit");

			public final String token;

			private DominantBaseline(String token) {
				this.token = token;
			}
		}

		public void setDominantBaseline(DominantBaseline val) {
			set(CSSConstants.CSS_DOMINANT_BASELINE_PROPERTY, val.token);
		}

		public static String renderNames(List<StyleClass> currentStyleClasses) {
			if (currentStyleClasses == null || currentStyleClasses.isEmpty()) {
				return "";
			} if (currentStyleClasses.size() == 1) {
				// handle the simple case efficiently
				return currentStyleClasses.get(0).name;
			} else {
				return String.join(" ", currentStyleClasses.stream()
					.filter((style) -> style != null)
					.map((style) -> style.name)
					.collect(Collectors.toList())
				);
			}
		}
	}

	private static interface Drawer {
		void draw();
	}

	public abstract class Drawable {

		private List<StyleClass> styleClasses = null;
		private String id = null;

		public Drawable setStyleClasses(StyleClass ... vals) {
			if (styleClasses == null) {
				styleClasses = new ArrayList<>();
			}
			styleClasses.clear();
			for (StyleClass val : vals) {
				if (val != null) {
					styleClasses.add(val);
				}
			}
			return this;
		}

		// TODO: "id" attribute doesn't seem to display for all SVG elements in inkscape?
		public Drawable setId(String val) {
			this.id = val;
			return this;
		}

		public void draw() {
			currentStyleClasses = styleClasses;
			currentId = id;
			reallyDraw();
		}

		protected abstract void reallyDraw();
	}

	public class ShapeDrawable extends Drawable {

		public final Shape shape;

		public ShapeDrawable(Shape shape) {
			this.shape = shape;
		}

		@Override
		protected void reallyDraw() {
			g.draw(shape);
		}
	}

	public class TextDrawable extends Drawable {

		public final String text;

		private double x = 0.0;
		private double y = 0.0;

		private String dy = null;

		public TextDrawable(String text) {
			this.text = text;
		}

		public TextDrawable setPos(double x, double y) {
			this.x = x;
			this.y = y;
			return this;
		}

		public TextDrawable setDY(double dy, LengthUnit unit) {
			// NOTE: this seems to apply before the y-axis flip transformation,
			// so flip the dy value here
			this.dy = unit.render(-dy);
			return this;
		}

		@Override
		protected void reallyDraw() {

			// NOTE: don't use g.drawString here, since AWT's text API is terrible
			// use the raw SVG api here instead

			// make the text SVG element
			Element elem = doc.createElementNS(SVGConstants.SVG_NAMESPACE_URI, SVGConstants.SVG_TEXT_TAG);
			elem.appendChild(doc.createTextNode(text));
			elem.setAttribute(SVGConstants.SVG_X_ATTRIBUTE, ctx.doubleString(x));
			elem.setAttribute(SVGConstants.SVG_Y_ATTRIBUTE, ctx.doubleString(y));
			elem.setAttributeNS(SVGConstants.XML_NAMESPACE_URI, SVGConstants.XML_SPACE_QNAME, SVGConstants.XML_PRESERVE_VALUE);

			// set optional attributes
			if (currentStyleClasses != null && !currentStyleClasses.isEmpty()) {
				elem.setAttribute(SVGConstants.SVG_CLASS_ATTRIBUTE, StyleClass.renderNames(currentStyleClasses));
			}
			if (currentId != null) {
				elem.setAttribute(SVGConstants.SVG_ID_ATTRIBUTE, currentId);
			}
			if (dy != null) {
				elem.setAttribute(SVGConstants.SVG_DY_ATTRIBUTE, dy);
			}

			groups.addElement(elem, DOMGroupManager.FILL);
		}
	}

	private final GenericDocument doc;
	private final SVGGeneratorContext ctx;
	private final SVGGraphics2D g;
	private final DOMGroupManager groups;
	private final Map<String,StyleClass> styleClasses = new HashMap<>();

	private List<StyleClass> currentStyleClasses = null;
	private String currentId = null;

	private Rectangle2D.Double bounds;

	public SVG() {

		// make the SVG DOM
		doc = (GenericDocument)GenericDOMImplementation.getDOMImplementation().createDocument(
			"http://www.w3.org/2000/svg",
			"svg",
			null
		);

		// get the graphics object that makes SVG DOM elements
		ctx = SVGGeneratorContext.createDefault(doc);
		ctx.setStyleHandler((Element elem, Map styleMap, SVGGeneratorContext ctxAgain) -> {

			// don't modify <g> elements
			if (elem.getTagName().equalsIgnoreCase("g")) {
				return;
			}

			// set the current style classes, if any
			if (currentStyleClasses != null && !currentStyleClasses.isEmpty()) {
				elem.setAttribute(SVGConstants.SVG_CLASS_ATTRIBUTE, StyleClass.renderNames(currentStyleClasses));
			}

			// set the current id, if any
			if (currentId != null) {
				elem.setAttribute(SVGConstants.SVG_ID_ATTRIBUTE, currentId);
			}
		});
		g = new SVGGraphics2D(ctx, false);

		// HACKHACK: get some internal stuff using reflection so we can add SVG elements directly
		try {
			Field groupsField = g.getClass().getDeclaredField("domGroupManager");
			groupsField.setAccessible(true);
			groups = (DOMGroupManager)groupsField.get(g);
		} catch (Throwable t) {
			throw new Error("Reflection hacks failed, maybe the Batik version changed?", t);
		}

		// set the default bonuds
		setBounds(0, 100, 0, 100);
	}

	public void setBounds(double x1, double x2, double y1, double y2) {
		bounds = new Rectangle2D.Double(
			x1, y1,
			x2 - x1,
			y2 - y1
		);
	}

	public StyleClass makeStyleClass(String name) {
		StyleClass styleClass = new StyleClass(name);
		putStyleClasses(styleClass);
		return styleClass;
	}

	public void putStyleClasses(StyleClass ... styleClasses) {
		for (StyleClass styleClass : styleClasses) {
			this.styleClasses.put(styleClass.name, styleClass);
		}
	}

	public ShapeDrawable makeLine(double x1, double y1, double x2, double y2) {
		return new ShapeDrawable(new Line2D.Double(
			x1, y1,
			x2, y2
		));
	}

	public ShapeDrawable makeRect(double x1, double x2, double y1, double y2) {
		return new ShapeDrawable(new Rectangle2D.Double(
			x1, y1,
			x2 - x1,
			y2 - y1
		));
	}

	public ShapeDrawable makeRectWH(double x, double y, double w, double h) {
		return new ShapeDrawable(new Rectangle2D.Double(
			x, y,
			w, h
		));
	}

	public ShapeDrawable makePoint(double x, double y, double radius) {
		double size = radius*2;
		return new ShapeDrawable(new Ellipse2D.Double(x - radius, y - radius, size, size));
	}

	public TextDrawable makeText(String text) {
		return new TextDrawable(text);
	}

	public Finished finish() {
		return new Finished();
	}

	public class Finished {

		private final Element svgElem;

		private Finished() {

			// move the SVG DOM inside g out to an external element
			svgElem = g.getRoot();
			Element topGroup = (Element)svgElem.getElementsByTagName("g").item(0);

			// apply the bounds and flip the y axis
			svgElem.setAttribute(SVGConstants.SVG_VIEW_BOX_ATTRIBUTE, String.format("%f %f %f %f", bounds.x, -(bounds.y + bounds.height), bounds.width, bounds.height));
			topGroup.setAttribute(SVGConstants.SVG_TRANSFORM_ATTRIBUTE, "scale(1,-1)");

			// flip all the text elements back
			for (Element textElem : elements(svgElem.getElementsByTagName("text"))) {
				try {
					double y = Double.parseDouble(textElem.getAttribute(SVGConstants.SVG_Y_ATTRIBUTE));
					textElem.setAttribute(SVGConstants.SVG_TRANSFORM_ATTRIBUTE, String.format("scale(1,-1) translate(0, %f) translate(0, %f)", -y, -y));
				} catch (NumberFormatException ex) {
					// something went wrong, don't flip the text
				}
			}

			// build the stylesheet from the classes
			CDATASection stylesheet = doc.createCDATASection("");

			// sort the style classes by priority
			List<StyleClass> sortedStyleClasses = new ArrayList<>(styleClasses.values());
			sortedStyleClasses.sort(Comparator.comparingInt(a -> a.priority));

			for (StyleClass styleClass : sortedStyleClasses) {
				stylesheet.appendData(".");
				stylesheet.appendData(styleClass.name);
				stylesheet.appendData(" { ");
				styleClass.values.forEach((key, value) -> {
					stylesheet.appendData(key);
					stylesheet.appendData(":");
					stylesheet.appendData(value);
					stylesheet.appendData("; ");
				});
				stylesheet.appendData("}\n");
			}

			// flip back all the text elements
			//stylesheet.appendData("g text { transform:scaleY(-1); }\n");

			// add the stylesheet to the svg dom
			Element style = doc.createElementNS(SVGSyntax.SVG_NAMESPACE_URI, SVGSyntax.SVG_STYLE_TAG);
			style.setAttributeNS(null, SVGSyntax.SVG_TYPE_ATTRIBUTE, "text/css");
			style.appendChild(stylesheet);
			doc.getChildElementById(svgElem, SVGSyntax.ID_PREFIX_GENERIC_DEFS).appendChild(style);
		}

		private Iterable<? extends Element> elements(NodeList nodes) {
			return () -> new Iterator<Element>() {

				int i = 0;
				int n = nodes.getLength();

				@Override
				public boolean hasNext() {
					return i < n;
				}

				@Override
				public Element next() {
					return (Element)nodes.item(i++);
				}
			};
		}

		public void write(Writer out) {
			try {
				g.stream(svgElem, out, true, false);
			} catch (IOException ex) {
				// checked exceptions are dumb...
				throw new RuntimeException(ex);
			}
		}

		public void write(File file) {
			try (FileWriter out = new FileWriter(file)) {
				write(out);
			} catch (IOException ex) {
				// checked exceptions are dumb...
				throw new RuntimeException(ex);
			}
		}
	}
}
