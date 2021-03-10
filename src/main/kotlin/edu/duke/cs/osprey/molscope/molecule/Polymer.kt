package edu.duke.cs.osprey.molscope.molecule

import edu.duke.cs.osprey.molscope.tools.Bijection


/**
 * An extension to Molecule that supports residues in a chain topology,
 * and separates mainchain from sidechain
 */
class Polymer(
	name: String
) : Molecule(name) {

	class Residue(
		/** A unique id for this residue in the chain, often a 4-character sequence number. */
		val id: String,
		/** A machine-readable description of this residue, often 3 characters. */
		var type: String
	) {
		val atoms: MutableList<Atom> = ArrayList()

		constructor(id: String, type: String, atoms: List<Atom>) : this(id, type) {
			this.atoms.addAll(atoms)
		}

		fun findAtomOrThrow(name: String) =
			atoms.find { it.name == name }
				?: throw NoSuchElementException("no atom with name=$name")

		override fun toString() = "$id:$type"
	}

	class Chain(
		/** A unique id for this chain, often a single character. */
		val id: String
	) {
		val residues: MutableList<Residue> = ArrayList()

		fun findResidueOrThrow(id: String) =
			residues
				.find { it.id == id }
				?: throw NoSuchElementException("no residue with id=$id")

		override fun toString() = id
	}

	val chains: MutableList<Chain> = ArrayList()

	override fun copy() =
		copyWithMaps().first

	override fun copyWithMaps(): Pair<Polymer,MoleculeMaps> {
		val src = this
		val dst = Polymer(name)
		val maps = src.copyTo(dst)
		return dst to maps
	}

	fun copyTo(dst: Polymer): PolymerMaps {
		val src = this

		val maps = src.copyTo(dst as Molecule)

		// copy all the chains
		val chainMap = ChainMap()
		val resMap = ResidueMap()
		for (srcChain in src.chains) {
			val dstChain = Chain(srcChain.id)
			chainMap.add(srcChain, dstChain)
			for (srcRes in srcChain.residues) {
				val dstRes = Residue(
					srcRes.id,
					srcRes.type,
					srcRes.atoms.map { maps.atoms.getBOrThrow(it) }
				)
				resMap.add(srcRes, dstRes)
				dstChain.residues.add(dstRes)
			}
			dst.chains.add(dstChain)
		}

		return PolymerMaps(maps, chainMap, resMap)
	}

	fun findChainOrThrow(id: String) =
		chains
			.find { it.id == id }
			?: throw NoSuchElementException("no chain with id=$id")

	fun findResidue(atom: Atom) =
		chains
			.flatMap { it.residues }
			.find { atom in it.atoms }

	fun findResidueOrThrow(atom: Atom) = findResidue(atom)
		?: throw NoSuchElementException("atom ${atom.name} was not found in any residue")

	fun findChainAndResidue(atom: Atom): Pair<Chain,Residue>? {
		for (chain in chains) {
			for (res in chain.residues) {
				if (atom in res.atoms) {
					return chain to res
				}
			}
		}
		return null
	}
}

class ChainMap : Bijection<Polymer.Chain>()
class ResidueMap : Bijection<Polymer.Residue>()
class PolymerMaps(
	maps: MoleculeMaps,
	val chains: ChainMap,
	val residues: ResidueMap
) : MoleculeMaps(maps)
