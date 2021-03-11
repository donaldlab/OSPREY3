package edu.duke.cs.osprey.gui.io

import cuchaz.kludge.tools.x
import cuchaz.kludge.tools.y
import cuchaz.kludge.tools.z
import edu.duke.cs.osprey.gui.compiler.AtomPairs
import edu.duke.cs.osprey.gui.compiler.CompiledConfSpace
import edu.duke.cs.osprey.gui.tools.UnsupportedClassException
import org.joml.Vector3dc
import java.io.ByteArrayOutputStream
import java.io.DataOutput
import java.io.DataOutputStream


/**
 * Compiled conf spaces are actually big enough that parsing the TOML
 * file can run the JVM out of heap space! So, we need a more efficient format...
 *
 * This is an attempt at encoding a compiled conf space with
 * an efficient binary format.
 *
 * NOTE: This does not produce human-readable output!
 */
fun CompiledConfSpace.toBytes(): ByteArray {

	val buf = ByteArrayOutputStream()
	val out = DataOutputStream(buf)

	// start with an 8-byte magic number
	// _ospccs_ = [Osp]rey [C]ompiled [C]onformation [S]pace
	out.writeBytes("_ospccs_")

	// then write a version number
	out.writeInt(1)

	// write out conf space properties
	out.writeUTF(name)
	out.writeInt(forcefields.size)
	for (forcefield in forcefields) {
		out.writeUTF(forcefield.name)
		out.writeUTF(forcefield.ospreyImplementation)
	}

	// write out forcefield settings
	for (ff in forcefields) {
		for ((key, value) in ff.settings) {
			when (value) {
				is String -> out.writeUTF(value)
				is Int -> out.writeInt(value)
				is Double -> out.writeDouble(value)
				is Boolean -> out.writeBoolean(value)
				else -> throw Error("don't know how to serialize forcefield setting: $key, which is a ${value::class}")
			}
		}
	}

	// write out the molecule infos
	out.writeInt(molInfos.size)
	for (molInfo in molInfos) {
		out.writeUTF(molInfo.name)
		if (molInfo.type != null) {
			out.writeByte(1)
			out.writeUTF(molInfo.type)
		} else {
			out.writeByte(0)
		}

		// write the motions, if any
		out.writeInt(molInfo.motions.size)
		for (motion in molInfo.motions) {
			when (motion) {
				is CompiledConfSpace.MotionInfo.DihedralAngle -> out.write(motion)
				is CompiledConfSpace.MotionInfo.TranslationRotation -> out.write(motion)
				else -> throw UnsupportedClassException("motion not supported on molecules", motion)
			}
		}
	}

	// write out the resdidue infos
	out.writeInt(resInfos.size)
	for (resInfo in resInfos) {
		out.writeUTF(resInfo.chainId)
		out.writeUTF(resInfo.id)
		out.writeUTF(resInfo.type)
		out.writeInt(resInfo.indexInChain)
	}

	fun DataOutput.write(atom: CompiledConfSpace.AtomInfo) {
		write(atom.pos)
		writeInt(atom.molIndex)
		writeInt(atom.resIndex)
		writeUTF(atom.name)
	}

	// write out the static atoms in order
	out.writeInt(staticAtoms.size)
	for (atom in staticAtoms) {
		out.write(atom)
	}

	// write out the internal energies for the static atoms
	for (energy in staticEnergies) {
		out.writeDouble(energy)
	}

	// write out the design positions
	out.writeInt(positions.size)
	for (pos in positions) {

		// write out the design position properties
		out.writeUTF(pos.name)
		out.writeUTF(pos.wildType)

		// write out the fragments
		out.writeInt(pos.fragments.size)

		// write out the conformations
		out.writeInt(pos.confs.size)
		for (conf in pos.confs) {

			// write out the conformation properties
			out.writeUTF(conf.id)
			out.writeUTF(conf.type)
			out.writeInt(conf.fragIndex)

			// write out the atoms for this conformation
			out.writeInt(conf.atoms.size)
			for (atom in conf.atoms) {
				out.write(atom)
			}

			// write the continuous motions, if needed
			out.writeInt(conf.motions.size)
			for (motion in conf.motions) {
				when (motion) {
					is CompiledConfSpace.MotionInfo.DihedralAngle -> out.write(motion)
					else -> throw UnsupportedClassException("motion not supported on conformations", motion)
				}
			}

			// write the internal energies
			for (ffi in forcefields.indices) {
				out.writeDouble(conf.internalEnergies[ffi])
			}
		}
	}

	// write out atom pairs
	for (ffi in atomPairs.indices) {
		out.write(atomPairs[ffi].static)
	}
	for ((posi1, pos1) in positions.withIndex()) {
		for (fragi1 in pos1.fragments.indices) {
			for (ffi in atomPairs.indices) {
				out.write(atomPairs[ffi].pos[posi1, fragi1])
				out.write(atomPairs[ffi].posStatic[posi1, fragi1])
			}
		}
	}

	// write out more atom pairs
	for ((posi1, pos1) in positions.withIndex()) {
		for (posi2 in 0 until posi1) {
			val pos2 = positions[posi2]

			for (fragi1 in pos1.fragments.indices) {
				for (fragi2 in pos2.fragments.indices) {

					// write the pos-pos atom pairs
					for (ffi in atomPairs.indices) {
						val pairs = atomPairs[ffi].posPos[posi1, fragi1, posi2, fragi2]
						out.writeInt(pairs.size)
						for (atomPair in pairs) {
							out.write(atomPair)
						}
					}
				}
			}
		}
	}

	// write the cached forcefield params
	for (atomPairs in atomPairs) {
		out.writeInt(atomPairs.paramsCache.size)
		for (params in atomPairs.paramsCache) {
			out.writeInt(params.size)
			for (param in params) {
				out.writeDouble(param)
			}
		}
	}

	// write the magic bytes again at the end, for error checking
	out.writeBytes("_ospccs_")

	return buf.toByteArray()
}

private fun DataOutput.write(v: Vector3dc) {
	writeDouble(v.x)
	writeDouble(v.y)
	writeDouble(v.z)
}

private fun DataOutput.write(atomPair: AtomPairs.AtomPair) {
	writeInt(atomPair.atomi1)
	writeInt(atomPair.atomi2)
	writeInt(atomPair.paramsi)
}

private fun DataOutput.write(atomPairs: List<AtomPairs.AtomPair>) {
	writeInt(atomPairs.size)
	for (atomPair in atomPairs) {
		write(atomPair)
	}
}

private fun DataOutput.write(motion: CompiledConfSpace.MotionInfo.DihedralAngle) {
	writeUTF("dihedralAngle")
	writeDouble(motion.minDegrees)
	writeDouble(motion.maxDegrees)
	writeInt(motion.abcd[0])
	writeInt(motion.abcd[1])
	writeInt(motion.abcd[2])
	writeInt(motion.abcd[3])
	writeInt(motion.rotated.size)
	for (atomi in motion.rotated) {
		writeInt(atomi)
	}
}

private fun DataOutput.write(motion: CompiledConfSpace.MotionInfo.TranslationRotation) {
	writeUTF("translationRotation")
	writeDouble(motion.maxTranslationDistance)
	writeDouble(motion.maxRotationRadians)
	write(motion.centroid)
}
