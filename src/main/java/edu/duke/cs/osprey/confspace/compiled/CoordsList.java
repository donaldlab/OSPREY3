package edu.duke.cs.osprey.confspace.compiled;

import org.joml.Vector3d;
import org.joml.Vector3dc;


/**
 * A supremely efficient representation for atomic coordinates.
 *
 * Optimized for lookup speed and read/writes to coordinates.
 *
 * Does not provide any information whatsoever about atom names,
 * atom elements, molecules, residues, bonds, etc.
 */
public class CoordsList {

	public final int size;

	private final double[] coords;

	public CoordsList(int size) {
		this.size = size;
		coords = new double[size*3];
	}

	public CoordsList(CoordsList other) {
		this(other.size);
		copyFrom(other, 0);
	}

	public double x(int i) {
		return coords[i*3];
	}
	public void setX(int i, double val) {
		coords[i*3] = val;
	}

	public double y(int i) {
		return coords[i*3 + 1];
	}
	public void setY(int i, double val) {
		coords[i*3 + 1] = val;
	}

	public double z(int i) {
		return coords[i*3 + 2];
	}
	public void setZ(int i, double val) {
		coords[i*3 + 2] = val;
	}

	public void get(int i, Vector3d out) {
		int o = i*3;
		out.x = coords[o];
		out.y = coords[++o];
		out.z = coords[++o];
	}

	public void set(int i, double x, double y, double z) {
		int o = i*3;
		coords[o] = x;
		coords[++o] = y;
		coords[++o] = z;
	}

	public void set(int i, Vector3dc in) {
		int o = i*3;
		coords[o] = in.x();
		coords[++o] = in.y();
		coords[++o] = in.z();
	}

	public void copyFrom(CoordsList src, int destIndex) {
		System.arraycopy(src.coords, 0, coords, destIndex*3, src.size*3);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof CoordsList && equals((CoordsList)other);
	}

	public boolean equals(CoordsList other) {

		if (this.size != other.size) {
			return false;
		}

		for (int i=0; i<coords.length; i++) {
			if (this.coords[i] != other.coords[i]) {
				return false;
			}
		}

		return true;
	}

	@Override
	public String toString() {
		return String.format("CoordsList[n=%d]", size);
	}

	public void print() {
		for (int i = 0; i < size; i++) {
			System.out.printf("%f %f %f%n", x(i), y(i), z(i));
		}
	}
}
