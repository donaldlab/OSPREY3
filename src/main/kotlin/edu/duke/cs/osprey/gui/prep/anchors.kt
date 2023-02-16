package edu.duke.cs.osprey.gui.prep

import cuchaz.kludge.tools.isFinite
import cuchaz.kludge.tools.toString
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.tools.normalizeZeroToTwoPI
import edu.duke.cs.osprey.gui.io.ConfLib
import org.joml.Quaterniond
import org.joml.Vector3d
import kotlin.math.abs
import kotlin.math.atan2


sealed class Anchor(val mol: Molecule) {

	abstract val anchorAtoms: List<Atom>

	/**
	 * Filters the given atoms down to just the ones connected to the anchor atoms.
	 */
	abstract fun getConnectedAtoms(atoms: Set<Atom>): MutableSet<Atom>
	abstract fun connectedAtomsIsSingleComponent(atoms: Set<Atom>): Boolean

	abstract fun isAnchorCompatible(anchor: ConfLib.Anchor): Boolean
	abstract fun bondToAnchors(anchor: ConfLib.Anchor, bondFunc: (List<ConfLib.AtomInfo>, Atom) -> Unit)

	/**
	 * Transforms the given atoms from the ConfLib coordinate system
	 * into the molecular coordinate system.
	 */
	abstract fun align(atoms: Set<Atom>, anchorCoords: ConfLib.AnchorCoords)

	/**
	 * Translate the molecule anchor atoms into a ConfLib anchor definition
	 */
	abstract fun makeLibraryAnchor(atoms: Set<Atom>, id: Int, getAtomInfo: (Atom) -> ConfLib.AtomInfo): ConfLib.Anchor
	abstract fun makeLibraryCoords(): ConfLib.AnchorCoords

	abstract fun copyToMol(mol: Molecule, oldToNew: AtomMap): Anchor

	/**
	 * Returns the residue of the first anchor atom, if any
	 */
	fun findResidue(): Polymer.Residue? =
		(mol as? Polymer)?.findResidue(anchorAtoms.first())

	// a few convenience functions for transforming a list of atoms
	protected fun Iterable<Atom>.makeLibraryAnchor(t: Vector3d) {
		for (atom in this) {
			atom.pos.add(t)
		}
	}
	protected fun Iterable<Atom>.rotate(q: Quaterniond) {
		for (atom in this) {
			atom.pos.rotate(q)
		}
	}

	protected fun Atom.getConnectedAtoms(included: Set<Atom>, blockers: Set<Atom> = emptySet()) =
		mol.bonds.bondedAtoms(this)
			.filter { it in included }
			.flatMap { sourceNeighbor ->
				mol
					.bfs(
						sourceNeighbor,
						visitSource = true,
						shouldVisit = { _, dst, _ -> dst in included && dst !in blockers }
					)
					.map { it.atom }
			}
			.toIdentitySet()

	class Single(
		mol: Molecule,
		val a: Atom,
		val b: Atom,
		val c: Atom
	) : Anchor(mol) {

		fun copy(mol: Molecule = this.mol, a: Atom = this.a, b: Atom = this.b, c: Atom = this.c) =
			Single(mol, a, b, c)

		override val anchorAtoms = listOf(a, b, c)

		override fun getConnectedAtoms(atoms: Set<Atom>) =
			a.getConnectedAtoms(included = atoms)

		override fun connectedAtomsIsSingleComponent(atoms: Set<Atom>) = true // this is trivially true, since there's only one anchor

		override fun isAnchorCompatible(anchor: ConfLib.Anchor) = anchor is ConfLib.Anchor.Single

		private fun ConfLib.Anchor.cast() = this as ConfLib.Anchor.Single
		private fun ConfLib.AnchorCoords.cast() = this as ConfLib.AnchorCoords.Single

		override fun bondToAnchors(anchor: ConfLib.Anchor, bondFunc: (List<ConfLib.AtomInfo>, Atom) -> Unit) {
			bondFunc(anchor.cast().bonds, a)
		}

		/**
		 * Aligns the attached atoms such that:
		 *   anchors a coincide
		 *   anchor vectors a-b are parallel and in the same direction
		 *   anchor planes a,b,c are parallel
		 * This alignment method exactly aligns the two coordinate systems without approximation.
		 */
		override fun align(atoms: Set<Atom>, anchorCoords: ConfLib.AnchorCoords) {

			val connectedAtoms = getConnectedAtoms(atoms)
			val coords = anchorCoords.cast()

			// center coords on a
			connectedAtoms.makeLibraryAnchor(
				Vector3d(coords.a)
					.negate()
			)

			// rotate so anchor a->b vectors are parallel, and the a,b,c planes are parallel
			connectedAtoms.rotate(
				Quaterniond()
					.lookAlong(
						Vector3d(coords.b).sub(coords.a),
						Vector3d(coords.c).sub(coords.a)
					)
					.apply {
						if (!isFinite()) {
							throw IllegalAlignmentException("conformation anchors a,b,c must not be co-linear, or nearly co-linear:\n\t" +
								listOf("a" to coords.a, "b" to coords.b, "c" to coords.c)
									.joinToString("\n\t") { (name, pos) -> "$name = ${pos.toString(2)}" }
							)
						}
					}
					.premul(
						Quaterniond()
							.lookAlong(
								Vector3d(b.pos).sub(a.pos),
								Vector3d(c.pos).sub(a.pos)
							)
							.apply {
								if (!isFinite()) {
									throw IllegalAlignmentException("design position anchor atoms a,b,c must not be co-linear, or nearly co-linear:\n\t" +
										listOf("a" to a, "b" to b, "c" to c)
											.joinToString("\n\t") { (name, atom) -> "$name = ${atom.name} ${atom.pos.toString(2)}" }
									)
								}
							}
							.conjugate()
					)
			)

			// translate back to a
			connectedAtoms.makeLibraryAnchor(a.pos)
		}

		override fun makeLibraryAnchor(atoms: Set<Atom>, id: Int, getAtomInfo: (Atom) -> ConfLib.AtomInfo) =
			ConfLib.Anchor.Single(
				id,
				// use sorted atoms here so we get a deterministic order
				bonds = mol.bonds.bondedAtomsSorted(a)
					.filter { it in atoms }
					.map { getAtomInfo(it) }
			)

		override fun makeLibraryCoords() =
			ConfLib.AnchorCoords.Single(
				Vector3d(a.pos),
				Vector3d(b.pos),
				Vector3d(c.pos)
			)

		override fun copyToMol(mol: Molecule, oldToNew: AtomMap) =
			Single(
				mol,
				a = oldToNew.getBOrThrow(a),
				b = oldToNew.getBOrThrow(b),
				c = oldToNew.getBOrThrow(c)
			)
	}

	class Double(
		mol: Molecule,
		val a: Atom,
		val b: Atom,
		val c: Atom,
		val d: Atom
	) : Anchor(mol) {

		fun copy(mol: Molecule = this.mol, a: Atom = this.a, b: Atom = this.b, c: Atom = this.c, d: Atom = this.d) =
			Double(mol, a, b, c, d)

		override val anchorAtoms = listOf(a, b, c, d)

		override fun getConnectedAtoms(atoms: Set<Atom>) =
			listOf(
				a.getConnectedAtoms(included = atoms),
				b.getConnectedAtoms(included = atoms)
			).union()

		override fun connectedAtomsIsSingleComponent(atoms: Set<Atom>) =
			listOf(
				a.getConnectedAtoms(included = atoms),
				b.getConnectedAtoms(included = atoms)
			).intersection().isNotEmpty()

		override fun isAnchorCompatible(anchor: ConfLib.Anchor) = anchor is ConfLib.Anchor.Double

		private fun ConfLib.Anchor.cast() = this as ConfLib.Anchor.Double
		private fun ConfLib.AnchorCoords.cast() = this as ConfLib.AnchorCoords.Double

		override fun bondToAnchors(anchor: ConfLib.Anchor, bondFunc: (List<ConfLib.AtomInfo>, Atom) -> Unit) {
			bondFunc(anchor.cast().bondsa, a)
			bondFunc(anchor.cast().bondsb, b)
		}

		/**
		 * Aligns the attached atoms such that:
		 *   anchors line segments a->b are in the same direction and have midpoints coincident at m.
		 *   anchor wedges c<-m->d are rotated about the a-b axis so their center directions are parallel.
		 * Since the different anchors can have different lengths of the a-b line segments,
		 * and different dihedral angles c->a->b->d, this method is an approximate method for
		 * aligning the two coordinate systems that tries to keep the error from accumulating
		 * all in one linear or angular distance.
		 */
		override fun align(atoms: Set<Atom>, anchorCoords: ConfLib.AnchorCoords) {

			val connectedAtoms = getConnectedAtoms(atoms)
			val coords = anchorCoords.cast()

			// let m = midpoint between a and b
			val coordsm = Vector3d(coords.a)
				.add(coords.b)
				.mul(0.5)
			val mpos = Vector3d(a.pos)
				.add(b.pos)
				.mul(0.5)

			/* DEBUG: show linear error
			val coordsALen = Vector3d(coords.a).sub(coordsm).length()
			val posALen = Vector3d(a.pos).sub(mpos).length()
			println("delta distance: ${abs(coordsALen - posALen)}")
			*/

			// center coords on m
			connectedAtoms.makeLibraryAnchor(
				Vector3d(coordsm)
					.negate()
			)

			// rotate into a coordinate system where:
			//   a->b becomes +z,
			//   a->c lies in the y,z plane
			val coordsQ = Quaterniond()
				.lookAlong(
					Vector3d(coords.b).sub(coords.a),
					Vector3d(coords.c).sub(coords.a)
				)
				.apply {
					if (!isFinite()) {
						throw IllegalAlignmentException("conformation anchors a,b,c must not be co-linear, or nearly co-linear:\n\t" +
							listOf("a" to coords.a, "b" to coords.b, "c" to coords.c)
								.joinToString("\n\t") { (name, pos) -> "$name = ${pos.toString(2)}" }
						)
					}
				}
			val posQ = Quaterniond()
				.lookAlong(
					Vector3d(b.pos).sub(a.pos),
					Vector3d(c.pos).sub(a.pos)
				)
				.apply {
					if (!isFinite()) {
						throw IllegalAlignmentException("design position anchor atoms a,b,c must not be co-linear, or nearly co-linear:\n\t" +
							listOf("a" to a, "b" to b, "c" to c)
								.joinToString("\n\t") { (name, atom) -> "$name = ${atom.name} ${atom.pos.toString(2)}" }
						)
					}
				}
			connectedAtoms.rotate(coordsQ)

			// measure the c->a->b->d dihedral angles in [0,2pi)
			val coordsDihedralRadians = abs(
				Vector3d(coords.d)
					.sub(coordsm)
					.rotate(coordsQ)
					.let { d -> atan2(d.y, d.x) } -
				Vector3d(coords.c)
					.sub(coordsm)
					.rotate(coordsQ)
					.let { c -> atan2(c.y, c.x) }
			).normalizeZeroToTwoPI()
			val posDihedralRadians = abs(
				Vector3d(d.pos)
					.sub(mpos)
					.rotate(posQ)
					.let { d -> atan2(d.y, d.x) } -
				Vector3d(c.pos)
					.sub(mpos)
					.rotate(posQ)
					.let { c -> atan2(c.y, c.x) }
			).normalizeZeroToTwoPI()

			/* DEBUG: show angular error
			println("delta angle: ${abs(coordsDihedralRadians - posDihedralRadians).toDegrees()/2.0}")
			*/

			// rotate about +z half the difference in dihedral angles
			connectedAtoms.rotate(
				Quaterniond()
					.rotationZ((posDihedralRadians - coordsDihedralRadians)/2.0)
			)

			// rotate back to the molecular frame where:
			//   +z becomes a->b
			//   +y lies in the a,b,c plane
			connectedAtoms.rotate(
				Quaterniond(posQ)
					.conjugate()
			)

			// translate back to m
			connectedAtoms.makeLibraryAnchor(mpos)
		}

		override fun makeLibraryAnchor(atoms: Set<Atom>, id: Int, getAtomInfo: (Atom) -> ConfLib.AtomInfo) =
			ConfLib.Anchor.Double(
				id,
				bondsa = mol.bonds.bondedAtoms(a)
					.filter { it in atoms }
					.map { getAtomInfo(it) },
				bondsb = mol.bonds.bondedAtoms(b)
					.filter { it in atoms }
					.map { getAtomInfo(it) }
			)

		override fun makeLibraryCoords() =
			ConfLib.AnchorCoords.Double(
				Vector3d(a.pos),
				Vector3d(b.pos),
				Vector3d(c.pos),
				Vector3d(d.pos)
			)

		override fun copyToMol(mol: Molecule, oldToNew: AtomMap) =
			Double(
				mol,
				a = oldToNew.getBOrThrow(a),
				b = oldToNew.getBOrThrow(b),
				c = oldToNew.getBOrThrow(c),
				d = oldToNew.getBOrThrow(d)
			)
	}
}

class IllegalAlignmentException(msg: String)
	: IllegalArgumentException("Can't align conformation to anchor:\n$msg")
