package edu.duke.cs.osprey.triplesBounds;


import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.tools.Progress;

import java.io.*;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * efficient storage for three-body energies
 */
public class TriplesEnergy {

	public final SimpleConfSpace confSpace;

	private int[] numRCsAtPos;
	private int[] offsets;
	private double[] energies;

	public TriplesEnergy(SimpleConfSpace confSpace) {

		this.confSpace = confSpace;

		// just in case...
		if (confSpace.positions.size() < 3) {
			throw new IllegalArgumentException("conformation space is too small for triples energies");
		}


		int numPos = confSpace.positions.size();
		numRCsAtPos = confSpace.getNumResConfsByPos();

		int numTriples = numPos*(numPos - 1)*(numPos - 2)/6;

		// compute the rc block offsets
		offsets = new int[numTriples];
		numTriples = 0;
		int offset = 0;
		for (int pos1=2; pos1<numPos; pos1++) {
			int numRCs1 = numRCsAtPos[pos1];
			for (int pos2=1; pos2<pos1; pos2++) {
				int numRCs2 = numRCsAtPos[pos2];
				for (int pos3=0; pos3<pos2; pos3++) {
					int numRCs3 = numRCsAtPos[pos3];

					offsets[numTriples] = offset;
					offset += numRCs1*numRCs2*numRCs3;
					numTriples++;
				}
			}
		}

		// allocate the space for the energies
		energies = new double[offset];
		fill(Double.NaN);
	}

	public void fill(double val) {
		Arrays.fill(energies, val);
	}

	public boolean isFullyDefined() {
		for (double energy : energies) {
			if (Double.isNaN(energy)) {
				return false;
			}
		}
		return true;
	}

	public int size() {
		return energies.length;
	}

	public double get(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		return energies[index(pos1, rc1, pos2, rc2, pos3, rc3)];
	}

	public double set(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, double val) {
		return energies[index(pos1, rc1, pos2, rc2, pos3, rc3)] = val;
	}

	private int index(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {

		// enforce pos1 > pos2 > pos3
		// can do it with a 3-step sort tree

		if (pos2 > pos1) {

			int swap = pos1;
			pos1 = pos2;
			pos2 = swap;

			swap = rc1;
			rc1 = rc2;
			rc2 = swap;

		} else if (pos2 == pos1) {
			throw new RuntimeException("Can't pair pos " + pos1 + " with itself");
		}

		if (pos3 > pos2) {

			int swap = pos2;
			pos2 = pos3;
			pos3 = swap;

			swap = rc2;
			rc2 = rc3;
			rc3 = swap;

		} else if (pos3 == pos2) {
			throw new RuntimeException("Can't pair pos " + pos2 + " with itself");
		}

		if (pos2 > pos1) {

			int swap = pos1;
			pos1 = pos2;
			pos2 = swap;

			swap = rc1;
			rc1 = rc2;
			rc2 = swap;

		} else if (pos2 == pos1) {
			throw new RuntimeException("Can't pair pos " + pos1 + " with itself");
		}

		//assert (pos1 > pos2); // the IntelliJ thinks this is always true, guess we don't need a runtime check
		assert (pos2 > pos3);

		return offsets[pos1*(pos1 - 1)*(pos1 - 2)/6 + pos2*(pos2 - 1)/2 + pos3]
			+ numRCsAtPos[pos3]*numRCsAtPos[pos2]*rc1 + numRCsAtPos[pos3]*rc2 + rc3;
	}

	public void calculateOrCache(EnergyCalculator ecalc, File file) {

		if (file.exists()) {
			try {
				read(file);
				log("read triples energies from file: " + file.getAbsolutePath());
				return;
			} catch (IllegalArgumentException ex) {
				// fall through to compute it again
			}
		}

		calculate(ecalc);

		write(file);
		log("wrote triples energies to file: " + file.getAbsolutePath());
	}

	public void calculate(EnergyCalculator ecalc) {

		int numPos = confSpace.positions.size();
		double weightSingle = 2.0/(numPos - 1)/(numPos - 2); // 1/(n-1 choose 2)
		double weightPair = 1.0/(numPos - 2); // 1/(n-2 choose 1)

		Progress progress = new Progress(size());
		log("calculating %d triple energies...", progress.getTotalWork());

		for (int pos1=2; pos1<numPos; pos1++) {
			final int fpos1 = pos1;
			int numRCs1 = numRCsAtPos[pos1];
			for (int pos2=1; pos2<pos1; pos2++) {
				final int fpos2 = pos2;
				int numRCs2 = numRCsAtPos[pos2];
				for (int pos3=0; pos3<pos2; pos3++) {
					final int fpos3 = pos3;
					int numRCs3 = numRCsAtPos[pos3];

					for (int rc1=0; rc1<numRCs1; rc1++) {
						final int frc1 = rc1;
						for (int rc2=0; rc2<numRCs2; rc2++) {
							final int frc2 = rc2;
							for (int rc3=0; rc3<numRCs3; rc3++) {
								final int frc3 = rc3;

								ecalc.tasks.submit(
									() -> {
										ResidueInteractions inters = ResInterGen.of(confSpace)
											.addIntra(fpos1, weightSingle, 0.0)
											.addIntra(fpos2, weightSingle, 0.0)
											.addIntra(fpos3, weightSingle, 0.0)
											.addShell(fpos1, weightSingle, 0.0)
											.addShell(fpos2, weightSingle, 0.0)
											.addShell(fpos3, weightSingle, 0.0)
											.addInter(fpos1, fpos2, weightPair, 0.0)
											.addInter(fpos1, fpos3, weightPair, 0.0)
											.addInter(fpos2, fpos3, weightPair, 0.0)
											.make();

										RCTuple frag = new RCTuple(fpos1, frc1, fpos2, frc2, fpos3, frc3);
										ParametricMolecule pmol = confSpace.makeMolecule(frag);
										return ecalc.calcEnergy(pmol, inters).energy;
									},
									(energy) -> {
										set(fpos1, frc1, fpos2, frc2, fpos3, frc3, energy);
										progress.incrementProgress();
									}
								);
							}
						}
					}
				}
			}
		}

		ecalc.tasks.waitForFinish();
	}

	public class GScorer implements AStarScorer {

		@Override
		public AStarScorer make() {
			return new GScorer();
		}

		@Override
		public double calc(ConfIndex index, RCs rcs) {

			double score = 0.0;

			for (int i1=2; i1<index.numDefined; i1++) {
				int pos1 = index.definedPos[i1];
				int rc1 = index.definedRCs[i1];

				for (int i2=1; i2<i1; i2++) {
					int pos2 = index.definedPos[i2];
					int rc2 = index.definedRCs[i2];

					for (int i3=0; i3<i2; i3++) {
						int pos3 = index.definedPos[i3];
						int rc3 = index.definedRCs[i3];

						score += get(pos1, rc1, pos2, rc2, pos3, rc3);
					}
				}
			}

			return score;
		}
	}

	public class HScorer implements AStarScorer {

		@Override
		public AStarScorer make() {
			return new HScorer();
		}

		@Override
		public double calc(ConfIndex index, RCs rcs) {

			double hscore = 0;

			// get the score for each undefined position
			for (int u1=0; u1<index.numUndefined; u1++) {
				int pos1 = index.undefinedPos[u1];

				// min over possible assignments to pos1
				double pos1Energy = Double.POSITIVE_INFINITY;
				for (int rc1 : rcs.get(pos1)) {

					double rc1Energy = 0.0;

					// add triples with defined positions
					for (int d2=0; d2<index.numDefined; d2++) {
						int pos2 = index.definedPos[d2];
						int rc2 = index.definedRCs[d2];

						for (int d3=0; d3<d2; d3++) {
							int pos3 = index.definedPos[d3];
							int rc3 = index.definedRCs[d3];

							rc1Energy += get(pos1, rc1, pos2, rc2, pos3, rc3);
						}
					}

					// add triples with undefined positions
					for (int u2=0; u2<u1; u2++) {
						int pos2 = index.undefinedPos[u2];

						// min over possible assignments to pos2
						double pos2Energy = Double.POSITIVE_INFINITY;
						for (int rc2 : rcs.get(pos2)) {

							double rc2Energy = 0.0;

							// add triples with defined positions
							for (int d3=0; d3<index.numDefined; d3++) {
								int pos3 = index.definedPos[d3];
								int rc3 = index.definedRCs[d3];
								rc2Energy += get(pos1, rc1, pos2, rc2, pos3, rc3);
							}

							// add triples with undefined positions
							for (int u3=0; u3<u2; u3++) {
								int pos3 = index.undefinedPos[u3];

								// min over possible assignments to pos3
								double pos3Energy = Double.POSITIVE_INFINITY;
								for (int rc3 : rcs.get(pos3)) {
									double rc3Energy = get(pos1, rc1, pos2, rc2, pos3, rc3);
									pos3Energy = Math.min(pos3Energy, rc3Energy);
								}

								rc2Energy += pos3Energy;
							}

							pos2Energy = Math.min(pos2Energy, rc2Energy);
						}

						rc1Energy += pos2Energy;
					}

					pos1Energy = Math.min(pos1Energy, rc1Energy);
				}

				hscore += pos1Energy;
			}

			return hscore;
		}
	}

	public void write(File file) {
		try (FileOutputStream out = new FileOutputStream(file)) {
			write(out);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	public void write(OutputStream out) {
		try {

			GZIPOutputStream zip = new GZIPOutputStream(out);
			DataOutputStream dout = new DataOutputStream(zip);

			dout.writeInt(energies.length);
			for (double energy : energies) {
				dout.writeDouble(energy);
			}

			zip.finish();

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	public void read(File file) {
		try (FileInputStream in = new FileInputStream(file)) {
			read(in);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	public void read(InputStream in) {
		try {

			GZIPInputStream zip = new GZIPInputStream(in);
			DataInputStream din = new DataInputStream(zip);

			if (din.readInt() != energies.length) {
				throw new IllegalArgumentException("stream does not match conf space");
			}

			for (int i=0; i<energies.length; i++) {
				energies[i] = din.readDouble();
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}
}
