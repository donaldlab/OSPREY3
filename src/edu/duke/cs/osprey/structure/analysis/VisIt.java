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

						out.write(String.format("%.2f", pd));
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
}
