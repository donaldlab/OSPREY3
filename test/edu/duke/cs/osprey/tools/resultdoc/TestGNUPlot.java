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


import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;


public class TestGNUPlot {

	private static void writePlot(Plot plot, String filename) {

		// no idea how to automate test checking here...
		// I just look at the SVG file manually
		FileTools.writeFile(plot.renderSvg(), new File(filename + ".svg"));
		FileTools.writeFileBytes(plot.renderPng(), new File(filename + ".png"));
	}

	@Test
	public void pointsX() {

		Plot plot = new Plot();
		plot.xlabel = "Value";
		plot.ylabels = Arrays.asList("A", "B", "C", "D");

		Plot.PointsX data1 = plot.new PointsX();
		data1.name = "Points 1";
		data1.data = Arrays.asList(0.0, 1.0, 5.0, 6.2);
		data1.size = 2.0;

		Plot.PointsX data2 = plot.new PointsX();
		data2.name = "Points 2";
		data2.data = Arrays.asList(3.3, 4.6, 9.3, 1.5);
		data2.size = 2.0;

		writePlot(plot, "pointsx");
	}

	@Test
	public void linesX() {

		Plot plot = new Plot();
		plot.xlabel = "Value";
		plot.ylabels = Arrays.asList("A", "B", "C", "D");

		Plot.LinesX data1 = plot.new LinesX();
		data1.name = "Lines 1";
		data1.data = Arrays.asList(0.0, 1.0, 5.0, 6.2);

		Plot.LinesX data2 = plot.new LinesX();
		data2.name = "Lines 2";
		data2.data = Arrays.asList(3.3, 4.6, 9.3, 1.5);

		writePlot(plot, "linesx");
	}

	@Test
	public void intervalsX() {

		Plot plot = new Plot();
		plot.xlabel = "Value";
		plot.ylabels = Arrays.asList("A", "B", "C");

		Plot.IntervalsX data = plot.new IntervalsX();
		data.name = "Intervals";
		data.data = Arrays.asList(
			new Plot.Interval(1.0, 2.5),
			new Plot.Interval(5.3, 7.2),
			new Plot.Interval(2.0, 4.0)
		);

		writePlot(plot, "intervalsx");
	}
}
