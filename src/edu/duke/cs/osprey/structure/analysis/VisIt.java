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

package edu.duke.cs.osprey.structure.analysis;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/** tool to write VisIt files for various data sets */
public class VisIt {

	public static void writeVoxels(List<SmallAngleVoxel> voxels, int d1, int d2, File file) {

		// convert the voxels to corner points
		List<double[]> corners = new ArrayList<>();

		for (SmallAngleVoxel voxel : voxels) {
			corners.add(new double[] { voxel.intervals[d1].min(), voxel.intervals[d2].min() });
			corners.add(new double[] { voxel.intervals[d1].min(), voxel.intervals[d2].max() });
			corners.add(new double[] { voxel.intervals[d1].max(), voxel.intervals[d2].min() });
			corners.add(new double[] { voxel.intervals[d1].max(), voxel.intervals[d2].max() });
		}

		writeAngles2D(corners, d1, d2, file);
	}

	public static void writeAngles2D(Collection<double[]> points, int d1, int d2, File file) {

		try (Writer out = new FileWriter(file)) {

			/* write out the vtk file, e.g.:
				# vtk DataFile Version 3.0
				any text here
				ASCII
				DATASET POLYDATA
				POINTS 4 float
				0 0 0
				1 0 0
				1.1 1.1 0
				0 1 0

				see: http://www.visitusers.org/index.php?title=ASCII_VTK_Files
				and: https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
			*/

			// write headers
			out.write("# vtk DataFile Version 3.0\n");
			out.write("whatever\n");
			out.write("ASCII\n");
			out.write("DATASET POLYDATA\n");
			out.write(String.format("POINTS %d float\n", points.size()));

			for (double[] p : points) {

				for (int i=0; i<3; i++) {
					if (i > 0) {
						out.write(" ");
					}
					if (i < 2) {

						// map to [0,360] for visualization
						double pd = p[i == 0 ? d1 : d2];
						while (pd < 0) {
							pd += 360;
						}
						while (pd >= 360) {
							pd -= 360;
						}

						assert (Double.isFinite(pd));

						out.write(String.format("%.4f", pd));
					} else {
						out.write("0");
					}
				}
				out.write("\n");
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	/** values stored in y-major order */
	public static void writeGrid2D(double[] x, double[] y, double[][] values, File file) {

		try (Writer out = new FileWriter(file)) {

			/* write out the vtk file, e.g.:
				# vtk DataFile Version 3.0
				beer is super awesome!
				ASCII
				DATASET RECTILINEAR_GRID
				DIMENSIONS 3 4 1
				X_COORDINATES 3 float
				0 2 4
				Y_COORDINATES 4 float
				1 2 3 4
				Z_COORDINATES 1 float
				0
				POINT_DATA 12
				FIELD FieldData 1
				cellscalar 1 12 float
				0 0 0
				1 1 1
				2 2 2
				3 3 3

				see: http://www.visitusers.org/index.php?title=ASCII_VTK_Files
				and: https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
			*/

			// write headers
			out.write("# vtk DataFile Version 3.0\n");
			out.write("whatever\n");
			out.write("ASCII\n");
			out.write("DATASET RECTILINEAR_GRID\n");
			out.write(String.format("DIMENSIONS %d %d %d\n", x.length, y.length, 1));

			out.write(String.format("X_COORDINATES %d float\n", x.length));
			for (int i=0; i<x.length; i++) {
				if (i % 10 > 0) {
					out.write(" ");
				} else if (i > 0) {
					out.write("\n");
				}
				out.write(String.format("%.4f", x[i]));
			}
			out.write("\n");

			out.write(String.format("Y_COORDINATES %d float\n", y.length));
			for (int i=0; i<y.length; i++) {
				if (i % 10 > 0) {
					out.write(" ");
				} else if (i > 0) {
					out.write("\n");
				}
				out.write(String.format("%.4f", y[i]));
			}
			out.write("\n");

			out.write("Z_COORDINATES 1 float\n");
			out.write("0\n");

			out.write(String.format("POINT_DATA %d\n", x.length*y.length));
			out.write("FIELD fieldDelta 1\n");
			out.write(String.format("delta 1 %d float\n", x.length*y.length));
			for (int iy=0; iy<y.length; iy++) {
				for (int ix=0; ix<x.length; ix++) {
					if (ix % 10 > 0) {
						out.write(" ");
					} else if (ix > 0) {
						out.write("\n");
					}
					out.write(String.format("%.4f", values[iy][ix]));
				}
				out.write("\n");
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}
}
