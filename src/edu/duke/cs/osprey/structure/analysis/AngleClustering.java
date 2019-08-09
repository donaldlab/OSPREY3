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


import de.lmu.ifi.dbs.elki.algorithm.clustering.hierarchical.AnderbergHierarchicalClustering;
import de.lmu.ifi.dbs.elki.algorithm.clustering.hierarchical.CentroidLinkageMethod;
import de.lmu.ifi.dbs.elki.algorithm.clustering.hierarchical.extraction.ExtractFlatClusteringFromHierarchy;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.DendrogramModel;
import de.lmu.ifi.dbs.elki.data.spatial.SpatialComparable;
import de.lmu.ifi.dbs.elki.data.type.TypeUtil;
import de.lmu.ifi.dbs.elki.data.type.VectorFieldTypeInformation;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.database.ids.DBIDIter;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.bundle.MultipleObjectsBundle;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.LPIntegerNormDistanceFunction;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Protractor;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class AngleClustering {

	/** cluster points in a toroidal space */
	public static List<List<double[]>> cluster(Collection<double[]> points, int numDimensions, double distCutoff) {

		// prep the database for ELKI
		Database db = new StaticArrayDatabase(() -> {

			List<DoubleVector> vecs = new ArrayList<>(points.size());
			for (double[] p : points) {
				vecs.add(new DoubleVector(p));
			}

			MultipleObjectsBundle b = new MultipleObjectsBundle();
			b.appendColumn(new VectorFieldTypeInformation<>(
				DoubleVector.FACTORY, numDimensions, numDimensions, DoubleVector.FACTORY.getDefaultSerializer()
			), vecs);
			return b;

		}, null);
		db.initialize();
		Relation<NumberVector> vectors = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);


		// use agglomerative clustering so we can use a custom distance metric
		LPIntegerNormDistanceFunction dist = new LPIntegerNormDistanceFunction(1) {

			@Override
			public double distance(NumberVector a, NumberVector b) {
				double dist = 0.0;
				for (int d=0; d<numDimensions; d++) {
					double dd = Protractor.getDistDegrees(a.doubleValue(d), b.doubleValue(d));
					dist += dd*dd;
				}
				return dist;
			}

			@Override
			public double norm(NumberVector v) {
				// hope we don't need this
				throw new UnsupportedOperationException();
			}

			@Override
			public double minDist(SpatialComparable a, SpatialComparable b) {
				// hope we don't need this
				throw new UnsupportedOperationException();
			}

			@Override
			public String toString() {
				return "ToroidalSquaredDistance";
			}
		};

		ExtractFlatClusteringFromHierarchy hac = new ExtractFlatClusteringFromHierarchy(
			new AnderbergHierarchicalClustering<>(
				dist,
				CentroidLinkageMethod.STATIC
			),
			distCutoff*distCutoff, // square the cutuff, since we're using a squared distance metric
			false,
			false
		);

		log("running clustering...");
		Clustering<DendrogramModel> result = hac.run(db);
		log("clustering done: %s clusters", result.getToplevelClusters().size());

		// collect the points for each cluster
		List<List<double[]>> clusters = new ArrayList<>();
		for (Cluster<DendrogramModel> cluster : result.getToplevelClusters()) {

			List<double[]> clusterPoints = new ArrayList<>(cluster.size());
			for (DBIDIter id1 = cluster.getIDs().iter(); id1.valid(); id1.advance()) {
				clusterPoints.add(vectors.get(id1).getColumnVector().getArrayCopy());
			}
			clusters.add(clusterPoints);

			log("\tcluster: n=%d", cluster.size());
		}

		return clusters;
	}

	public static double[] calcMedoid(Collection<double[]> points) {

		// calculating centroids in curved spaces is hard, let's try medoids instead
		// this is a navie O(n2) algorithm, but it's fast enough for this amount of data
		double[] medoid = null;
		double minScore = Double.POSITIVE_INFINITY;
		for (double[] p1 : points) {

			double score = 0.0;
			for (double[] p2 : points) {

				// add the squared toroidal distance to the score
				for (int d=0; d<p1.length; d++) {
					double dd = Protractor.getDistDegrees(p1[d], p2[d]);
					score += dd*dd;
				}
			}

			// is this the best one?
			if (medoid == null || score < minScore) {
				minScore = score;
				medoid = p1;
			}
		}

		return medoid;
	}

	public static SmallAngleVoxel calcVoxel(List<double[]> points) {
		SmallAngleVoxel voxel = new SmallAngleVoxel(calcMedoid(points));
		for (double[] p : points) {
			voxel.expand(p);
		}
		return voxel;
	}

	public static SmallAngleVoxel fitFixedVoxel(Collection<double[]> points, SmallAngleVoxel voxel, double radiusDegrees) {

		int n = voxel.intervals.length;

		int[] sizes = new int[n];
		for (int d=0; d<n; d++) {
			sizes[d] = Math.max(1, (int)Math.ceil(voxel.intervals[d].size() - radiusDegrees*2));
		}

		double[] center = new double[n];
		SmallAngleVoxel fixedVoxel = new SmallAngleVoxel(center);
		for (int d=0; d<n; d++) {
			fixedVoxel.intervals[d].less = 0;
			fixedVoxel.intervals[d].more = radiusDegrees*2;
		}

		long maxCount = 0;
		int[] bestIndices = null;
		for (int[] indices : new MathTools.GridIterable(sizes)) {

			// set the center
			for (int d=0; d<n; d++) {
				fixedVoxel.intervals[d].center = Protractor.normalizeDegrees(voxel.intervals[d].min() + indices[d]);
			}

			long count = points.stream()
				.filter(p -> fixedVoxel.contains(p))
				.count();

			if (count > maxCount) {
				maxCount = count;
				bestIndices = indices.clone();
			}
		}
		assert (bestIndices != null);

		// reposition fixed voxel to the best spot
		for (int d=0; d<n; d++) {
			fixedVoxel.intervals[d].center = Protractor.normalizeDegrees(voxel.intervals[d].min() + bestIndices[d]);
		}

		return fixedVoxel;
	}
}
