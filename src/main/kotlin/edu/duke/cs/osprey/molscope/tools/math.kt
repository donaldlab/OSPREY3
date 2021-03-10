package edu.duke.cs.osprey.molscope.tools

import org.joml.AABBf
import org.joml.Vector3f
import kotlin.math.PI


fun Float.toStringAngstroms() = "%.3f".format(this)

fun Vector3f.toStringAngstroms() = "(%.3f,%.3f,%.3f)".format(x, y, z)

fun AABBf.toStringAngstroms() = "[%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f]".format(minX, maxX, minY, maxY, minZ, maxZ)

const val TwoPI = PI*2

/**
 * Normalizes the radians into the range (-pi,pi]
 */
fun Double.normalizeMinusPIToPI(): Double {
	var v = this
	while (v <= -PI) {
		v += TwoPI
	}
	while (v > PI) {
		v -= TwoPI
	}
	return v
}

/**
 * Normalizes the radians into the range [0,2pi)
 */
fun Double.normalizeZeroToTwoPI(): Double {
	var v = this
	while (v < 0) {
		v += TwoPI
	}
	while (v >= TwoPI) {
		v -= TwoPI
	}
	return v
}
