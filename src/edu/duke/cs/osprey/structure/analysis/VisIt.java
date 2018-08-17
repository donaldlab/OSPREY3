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

	public static void writeVoxels(List<SmallAngleVoxel> voxels, File file) {

		// convert the voxels to corner points
		List<double[]> corners = new ArrayList<>();
		for (SmallAngleVoxel voxel : voxels) {

			assert (voxel.intervals.length == 2);

			corners.add(new double[] { voxel.intervals[0].min(), voxel.intervals[1].min() });
			corners.add(new double[] { voxel.intervals[0].min(), voxel.intervals[1].max() });
			corners.add(new double[] { voxel.intervals[0].max(), voxel.intervals[1].min() });
			corners.add(new double[] { voxel.intervals[0].max(), voxel.intervals[1].max() });
		}

		writeDihedrals(corners, file);
	}

	public static void writeDihedrals(Collection<double[]> points, File file) {

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

				assert (p.length == 2);

				for (int d=0; d<3; d++) {
					if (d > 0) {
						out.write(" ");
					}
					if (d < p.length) {

						// map to [0,360] for visualization
						double pd = p[d];
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
