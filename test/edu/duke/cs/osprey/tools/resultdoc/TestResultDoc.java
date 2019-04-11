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


import org.junit.Test;

import java.util.Arrays;

public class TestResultDoc {

	@Test
	public void test() {

		try (ResultDoc doc = new ResultDoc("test.md")) {

			doc.h1("Hello World");
			doc.println("Look at this `awesome` result!");

			doc.println();

			Plot plot1 = new Plot();
			Plot.LinesY lines1 = plot1.new LinesY();
			lines1.name = "A";
			lines1.data = Arrays.asList(3.0, 1.0, 2.0, 4.0);
			Plot.LinesY lines2 = plot1.new LinesY();
			lines2.data = Arrays.asList(5.0, 7.0, 6.0);
			doc.plot(plot1);

			doc.println();

			Plot plot2 = new Plot();
			Plot.IntervalsX intervals = plot2.new IntervalsX();
			intervals.data = Arrays.asList(
				new Plot.Interval(1.0, 2.0),
				new Plot.Interval(2.0, 3.0),
				new Plot.Interval(3.0, 4.0)
			);
			Plot.PointsX points = plot2.new PointsX();
			points.data = Arrays.asList(1.4, 2.6, 3.9);
			points.type = Plot.PointType.TriangleUp;
			plot2.ylabels = Arrays.asList("Abracadabra", "Doofingus", "Ghabrolamabama");
			doc.plot(plot2);

			doc.println();
			doc.hr();
			doc.println();

			doc.h2("Foods!");

			doc.tableHeader("Name", "Type", "Tasty?");
			doc.tableRow("Provolone", "Cheese", "Yes");
			doc.tableRow("Gouda", "Cheese", "Yes");
			doc.tableRow("Apple", "Fruit", "Yes");
			doc.tableRow("Banana", "Fruit", "No");

			doc.println();

			StringBuilder tsv = new StringBuilder();
			tsv.append("Name\tType\tTasty?\n");
			tsv.append("Provolone\tCheese\tYes\n");
			tsv.append("Gouda\tCheese\tYes\n");
			tsv.append("Apple\tFruit\tYes\n");
			tsv.append("Banana\tFruit\tNo\n");
			doc.file.tsv(tsv.toString(), "Download foods.tsv");

			doc.println();

			doc.println("And that's the whole story.");
		}

		// no idea how to automate test checking here...
		// I just look at the MD file manually
	}
}
