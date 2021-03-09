package edu.duke.cs.osprey.service

import kotlinx.serialization.Serializable


@Serializable
data class Point3d(val x: Double, val y: Double, val z: Double)
