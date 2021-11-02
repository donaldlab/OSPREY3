package edu.duke.cs.osprey.gui.io

import cuchaz.kludge.tools.x
import cuchaz.kludge.tools.y
import cuchaz.kludge.tools.z
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.tools.identityHashSet
import org.joml.Vector3d
import org.joml.Vector3dc
import org.tomlj.*
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap


/**
 * Conformation Library
 */
class ConfLib(
	/** should be a valid TOML key */
	val id: String,
	/** any human-readable short name */
	val name: String,
	val fragments: Map<String,Fragment>,
	val description: String? = null,
	val citation: String? = null
) {

	/**
	 * A globally unique id for each library, assigned at runtime.
	 * This id is not intrinsic to the library itself, so should not be persisted anywhere.
	 * Its only purpose is to allow assigning globally unique ids to fragments, conformations, etc at runtime.
	 */
	val runtimeId = "library-${nextId++}"

	/**
	 * A globally unique id for each fragment, assigned at runtime.
	 * This id is not intrinsic to the library itself, so should not be persisted anywhere.
	 */
	fun fragRuntimeId(frag: Fragment) =
		"$runtimeId.${frag.id}"

	/**
	 * A globally unique id for each conformation, assigned at runtime.
	 * This id is not intrinsic to the library itself, so should not be persisted anywhere.
	 */
	fun confRuntimeId(frag: Fragment, conf: Conf) =
		"$runtimeId.${frag.id}.${conf.id}"

	fun getFragmentOrThrow(fragId: String) =
		fragments[fragId] ?: throw NoSuchElementException("no fragment found with id $fragId")

	override fun toString() = "$name ($runtimeId)"

	data class AtomInfo(
		val id: Int,
		val name: String,
		val element: Element
	) : AtomPointer {

		override fun matchIn(atoms: List<AtomInfo>, anchors: List<Anchor>) =
			atoms.find { it.name == name }

		override fun resolveCoords(conf: Conf) =
			conf.coords[this]
	}

	data class Bond(
		val a: AtomInfo,
		val b: AtomInfo
	)

	sealed class Anchor {

		abstract val id: Int

		abstract fun matchIn(atoms: List<AtomInfo>, anchors: List<Anchor>): Anchor?

		/**
		 * Find all atoms bonded to this anchor.
		 */
		abstract fun findAtoms(frag: Fragment): List<AtomInfo>

		data class Single(
			override val id: Int,
			val bonds: List<AtomInfo>
		) : Anchor() {

			override fun matchIn(atoms: List<AtomInfo>, anchors: List<Anchor>) =
				anchors
					.filterIsInstance<Single>()
					.find { it.bonds.all { it.matchIn(atoms, anchors) != null } }

			override fun findAtoms(frag: Fragment) =
				this.bonds
					.flatMap { findAtoms(frag, it) }
					.toCollection(identityHashSet())
					.toList()
		}

		data class Double(
			override val id: Int,
			val bondsa: List<AtomInfo>,
			val bondsb: List<AtomInfo>
		) : Anchor() {

			override fun matchIn(atoms: List<AtomInfo>, anchors: List<Anchor>) =
				anchors
					.filterIsInstance<Double>()
					.find {
						it.bondsa.all { it.matchIn(atoms, anchors) != null }
						&& it.bondsb.all { it.matchIn(atoms, anchors) != null }
					}

			override fun findAtoms(frag: Fragment) =
				listOf(bondsa, bondsb)
					.flatMap { bonds ->
						bonds.flatMap { findAtoms(frag, it) }
					}
					.toCollection(identityHashSet())
					.toList()
		}

		companion object {

			fun findAtoms(frag: Fragment, source: AtomInfo): List<AtomInfo> {

				// do BFS in the bond graph
				val toVisit = ArrayDeque<AtomInfo>()
				val visitScheduled = identityHashSet<AtomInfo>()

				fun scheduleVisit(atom: AtomInfo) {
					toVisit.add(atom)
					visitScheduled.add(atom)
				}

				// start with the source
				scheduleVisit(source)

				val out = ArrayList<AtomInfo>()

				while (true) {

					// visit the next atom
					val atom = toVisit.pollFirst() ?: break
					out.add(atom)

					// schedule visits to neighbors
					for (bondedAtom in frag.bondedAtoms(atom)) {
						if (bondedAtom !in visitScheduled) {
							scheduleVisit(bondedAtom)
						}
					}
				}

				return out
			}
		}
	}

	sealed class AnchorCoords {

		abstract fun getCoords(index: Int): Vector3d?

		data class Single(
			val a: Vector3d,
			val b: Vector3d,
			val c: Vector3d
		) : AnchorCoords() {

			override fun getCoords(index: Int) =
				when (index) {
					0 -> a
					1 -> b
					2 -> c
					else -> null
				}
		}

		data class Double(
			val a: Vector3d,
			val b: Vector3d,
			val c: Vector3d,
			val d: Vector3d
		) : AnchorCoords() {

			override fun getCoords(index: Int) =
				when (index) {
					0 -> a
					1 -> b
					2 -> c
					3 -> d
					else -> null
				}
		}
	}

	data class Conf(
		val id: String,
		val name: String,
		val description: String?,
		val coords: Map<AtomInfo,Vector3d>,
		val anchorCoords: Map<Anchor,AnchorCoords>
	) {

		fun uniqueId(conflib: ConfLib, frag: Fragment) = "${conflib.id}/${frag.id}/$id"
	}

	interface AtomPointer {

		/**
		 * Returns a new atom pointer that pointes to an equivalent same atom in a specified fragment.
		 */
		fun matchIn(atoms: List<AtomInfo>, anchors: List<Anchor>): AtomPointer?

		fun matchInOrThrow(atoms: List<AtomInfo>, anchors: List<Anchor>): AtomPointer =
			matchIn(atoms, anchors) ?: throw NoSuchElementException("no match found for $this")

		/**
		 * Looks up the atom coordinates in the conformation,
		 * whether the pointer be to a conformation atom or an anchor atom.
		 */
		fun resolveCoords(conf: Conf): Vector3d?

		fun resolveCoordsOrThrow(conf: Conf) =
			resolveCoords(conf) ?: throw NoSuchElementException("no coords found for atom pointer $this")

		companion object {
			// defined so we can extend it
		}

		interface Resolver {

			fun resolve(p: AtomPointer): Atom?

			fun resolveOrThrow(p: AtomPointer) =
				resolve(p) ?: throw NoAtomException(p)
		}

		class NoAtomException(val p: AtomPointer)
			: RuntimeException("can't find atom from pointer: $p")
	}

	data class AnchorAtomPointer(
		val anchor: Anchor,
		val index: Int
	) : AtomPointer {

		@Suppress("UselessCallOnCollection") // IntelliJ warning on filterIsInstance() is incorrect here T_T
		override fun matchIn(atoms: List<AtomInfo>, anchors: List<Anchor>) =
			anchor.matchIn(atoms, anchors)
				?.let {
					AnchorAtomPointer(it, index)
				}

		override fun resolveCoords(conf: Conf) =
			conf.anchorCoords[anchor]?.getCoords(index)


		override fun toString() = "${this::class.simpleName}[${anchor.id},$index]"
	}

	sealed class ContinuousMotion {

		abstract val id: Int
		abstract fun affectedAtoms(frag: Fragment): List<AtomInfo>

		/**
		 * Make a copy of this degree of freedom, but pointing to the given atoms and anchors,
		 * instead of the atoms and anchors it already has.
		 */
		abstract fun copyTo(atoms: List<AtomInfo>, anchors: List<Anchor>, id: Int): ContinuousMotion

		data class DihedralAngle(
			override val id: Int,
			val a: AtomPointer,
			val b: AtomPointer,
			val c: AtomPointer,
			val d: AtomPointer
		) : ContinuousMotion() {

			override fun affectedAtoms(frag: Fragment): List<AtomInfo> {

				// b might be an anchor atom
				val b = b as? AtomInfo

				// but if c is an anchor atom, that means all the sidechain atoms move, right?
				val c = c as? AtomInfo
					?: return frag.atoms

				// grab all the atoms connected to c not through b-c
				return frag.bfs(
					source = c,
					visitSource = false,
					shouldVisit = { from, to -> to !== b }
				).toList()
			}

			override fun copyTo(atoms: List<AtomInfo>, anchors: List<Anchor>, id: Int) =
				DihedralAngle(
					id,
					a.matchInOrThrow(atoms, anchors),
					b.matchInOrThrow(atoms, anchors),
					c.matchInOrThrow(atoms, anchors),
					d.matchInOrThrow(atoms, anchors)
				)
		}

		// TODO: other motions?
	}

	data class Fragment(
		val id: String,
		val name: String,
		val type: String,
		val atoms: List<AtomInfo>,
		val bonds: List<Bond>,
		val anchors: List<Anchor>,
		val confs: Map<String,Conf>,
		val motions: List<ContinuousMotion>
	) {

		fun uniqueId(conflib: ConfLib) = "${conflib.id}/$id"

		fun bondedAtoms(atom: AtomInfo): List<AtomInfo> =
			bonds
				.mapNotNull { (a, b) ->
					when {
						atom === a -> b
						atom === b -> a
						else -> null
					}
				}

		private val atomsByAnchor = IdentityHashMap<Anchor,List<AtomInfo>>()

		fun getAtomsFor(anchor: Anchor) =
			atomsByAnchor.getOrPut(anchor) {
				anchor.findAtoms(this)
			}

		/**
		 * Walk the bond graph using breadth-first search.
		 */
		fun bfs(
			source: AtomInfo,
			visitSource: Boolean = false,
			shouldVisit: (fromAtom: AtomInfo, toAtom: AtomInfo) -> Boolean = { _, _ -> true }
		) = object : Iterable<AtomInfo> {

			override fun iterator() = object : Iterator<AtomInfo> {

				// track the atom visitation schedule
				val toVisit = ArrayDeque<AtomInfo>()
				val visitScheduled = identityHashSet<AtomInfo>()

				fun scheduleVisit(atom: AtomInfo) {
					toVisit.add(atom)
					visitScheduled.add(atom)
				}

				override fun hasNext() = toVisit.isNotEmpty()

				override fun next(): AtomInfo {

					// take the next step
					val atom = toVisit.pollFirst()
						?: throw NoSuchElementException("no more atoms to visit")

					// schedule visits to neighbors
					for (neighbor in bondedAtoms(atom)) {
						if (neighbor !in visitScheduled && shouldVisit(atom, neighbor)) {
							scheduleVisit(neighbor)
						}
					}

					return atom
				}

				init {

					// seed with the source atom
					scheduleVisit(source)

					// skip the source atom, if need
					if (!visitSource) {
						next()
					}
				}
			}
		}

		fun getConfs(vararg id: String): MutableSet<Conf> =
			id
				.map { confs.getValue(it) }
				.toCollection(identityHashSet())
	}

	companion object {

		private var nextId = 0

		fun from(toml: String): ConfLib {

			// parse the TOML
			val doc = Toml.parse(toml)
			if (doc.hasErrors()) {
				throw TomlParseException("TOML parsing failure:\n${doc.errors().joinToString("\n")}")
			}

			return from(doc)
		}

		fun from(doc: TomlTable): ConfLib {

			// read the header
			val libId = doc.getStringOrThrow("id")
			val libName = doc.getStringOrThrow("name")
			val libDesc = doc.getString("description")
			val citation = doc.getString("citation")

			// read the fragments
			val frags = fragmentsFrom(doc)

			return ConfLib(libId, libName, frags, libDesc, citation)
		}

		fun fragmentsFrom(toml: String): Map<String,Fragment> {

			// parse the TOML
			val doc = Toml.parse(toml)
			if (doc.hasErrors()) {
				throw TomlParseException("TOML parsing failure:\n${doc.errors().joinToString("\n")}")
			}

			return fragmentsFrom(doc)
		}

		fun fragmentsFrom(doc: TomlTable): Map<String,Fragment> {

			val frags = HashMap<String,Fragment>()

			val fragsTable = doc.getTable("frag") ?: return frags

			for (fragId in fragsTable.keySet()) {
				val fragTable = fragsTable.getTableOrThrow(fragId)
				val fragPos = fragsTable.inputPositionOf(fragId)

				val fragName = fragTable.getString("name") ?: fragId
				val fragType = fragTable.getStringOrThrow("type")

				// read the atoms
				val atoms = HashMap<Int,AtomInfo>()
				val atomsArray = fragTable.getArrayOrThrow("atoms", fragPos)
				for (i in 0 until atomsArray.size()) {
					val atomTable = atomsArray.getTable(i)
					val pos = atomsArray.inputPositionOf(i)

					val id = atomTable.getIntOrThrow("id", pos)
					atoms[id] = AtomInfo(
						id,
						atomTable.getStringOrThrow("name", pos),
						Element[atomTable.getStringOrThrow("elem", pos)]
					)
				}

				fun getAtom(id: Int, pos: TomlPosition? = null) =
					atoms[id] ?: throw TomlParseException("no atom with id $id in fragment $fragId", pos)

				// read the bonds
				val bonds = ArrayList<Bond>()
				val bondsArray = fragTable.getArrayOrThrow("bonds", fragPos)
				for (i in 0 until bondsArray.size()) {
					val bondArray = bondsArray.getArray(i)
					val pos = bondsArray.inputPositionOf(i)

					bonds.add(Bond(
						getAtom(bondArray.getInt(0), pos),
						getAtom(bondArray.getInt(1), pos)
					))
				}

				fun TomlArray.toBonds(pos: TomlPosition): List<AtomInfo> =
					(0 until size())
						.map { i -> getAtom(getInt(i), pos) }

				// read the anchors
				val anchors = HashMap<Int,Anchor>()
				val anchorArray = fragTable.getArrayOrThrow("anchors", fragPos)
				for (i in 0 until anchorArray.size()) {
					val anchorTable = anchorArray.getTable(i)
					val pos = anchorArray.inputPositionOf(i)

					val id = anchorTable.getIntOrThrow("id", pos)
					val type = anchorTable.getString("type")?.toLowerCase()
					anchors[id] = when (type) {
						"single" -> {
							Anchor.Single(
								id,
								anchorTable.getArrayOrThrow("bonds", pos).toBonds(pos)
							)
						}
						"double" -> {
							Anchor.Double(
								id,
								anchorTable.getArrayOrThrow("bondsa", pos).toBonds(pos),
								anchorTable.getArrayOrThrow("bondsb", pos).toBonds(pos)
							)
						}
						else -> throw TomlParseException("unrecognized anchor type: $type", pos)
					}
				}

				fun getAnchor(id: Int, pos: TomlPosition? = null) =
					anchors[id] ?: throw TomlParseException("no anchor with id $id in fragment $fragId", pos)

				fun TomlArray.toVector3d() =
					Vector3d(
						getDouble(0),
						getDouble(1),
						getDouble(2)
					)

				// read the confs
				val confs = HashMap<String,Conf>()
				val confsTable = fragTable.getTableOrThrow("conf", fragPos)
				for (confId in confsTable.keySet()) {
					val confTable = confsTable.getTableOrThrow(confId)
					val confPos = confsTable.inputPositionOf(confId)

					val name = confTable.getString("name") ?: confId
					val desc = confTable.getString("description")
					val coords = IdentityHashMap<AtomInfo,Vector3d>()
					val anchorCoords = IdentityHashMap<Anchor,AnchorCoords>()

					val coordsArray = confTable.getArrayOrThrow("coords", confPos)
					for (i in 0 until coordsArray.size()) {
						val coordTable = coordsArray.getTable(i)
						val pos = coordsArray.inputPositionOf(i)

						val atomInfo = getAtom(coordTable.getIntOrThrow("id", pos))
						coords[atomInfo] = coordTable.getArrayOrThrow("xyz", pos).toVector3d()
					}

					val anchorCoordsArray = confTable.getArrayOrThrow("anchorCoords", confPos)
					for (i in 0 until anchorCoordsArray.size()) {
						val coordTable = anchorCoordsArray.getTable(i)
						val pos = anchorCoordsArray.inputPositionOf(i)

						val anchor = getAnchor(coordTable.getIntOrThrow("id", pos))
						anchorCoords[anchor] = when (anchor) {
							is Anchor.Single -> AnchorCoords.Single(
								a = coordTable.getArrayOrThrow("a", pos).toVector3d(),
								b = coordTable.getArrayOrThrow("b", pos).toVector3d(),
								c = coordTable.getArrayOrThrow("c", pos).toVector3d()
							)
							is Anchor.Double -> AnchorCoords.Double(
								a = coordTable.getArrayOrThrow("a", pos).toVector3d(),
								b = coordTable.getArrayOrThrow("b", pos).toVector3d(),
								c = coordTable.getArrayOrThrow("c", pos).toVector3d(),
								d = coordTable.getArrayOrThrow("d", pos).toVector3d()
							)
						}
					}

					confs[confId] = Conf(confId, name, desc, coords, anchorCoords)
				}

				fun makeAtomPointer(key: Any?, pos: TomlPosition): AtomPointer =
					when (key) {
						is Long -> {
							// look up the atom info by index
							getAtom(key.toInt(), pos)
						}
						is TomlArray -> {
							// make an anchor atom pointer from the two indices
							if (key.size() == 2 && key.containsLongs()) {
								val anchorIndex = key.getInt(0)
								val atomIndex = key.getInt(1)
								AnchorAtomPointer(
									getAnchor(anchorIndex),
									atomIndex - 1
								)
							} else {
								throw TomlParseException("atom pointer doesn't look like an anchor atom pointer: $key", pos)
							}
						}
						else -> throw TomlParseException("unrecognized atom pointer: $key", pos)
					}

				// read the motions
				val motions = HashMap<Int,ContinuousMotion>()
				val motionsArray = fragTable.getArrayOrThrow("motions", fragPos)
				for (i in 0 until motionsArray.size()) {
					val motionTable = motionsArray.getTable(i)
					val pos = motionsArray.inputPositionOf(i)

					val id = motionTable.getIntOrThrow("id", pos)
					val type = motionTable.getStringOrThrow("type", pos)

					when (type) {
						"dihedral" -> {
							motions[id] = ContinuousMotion.DihedralAngle(
								id,
								makeAtomPointer(motionTable.get("a"), pos),
								makeAtomPointer(motionTable.get("b"), pos),
								makeAtomPointer(motionTable.get("c"), pos),
								makeAtomPointer(motionTable.get("d"), pos)
							)
						}
						else -> throw TomlParseException("unrecognized degree of freedom type: $type", pos)
					}
				}

				frags[fragId] = Fragment(
					fragId,
					fragName,
					fragType,
					atoms.map { (_, atom) -> atom },
					bonds,
					anchors.map { (_, anchor) -> anchor },
					confs,
					motions.map { (_, motion) -> motion }
				)
			}

			return frags
		}
	}
}

fun ConfLib.toToml(table: String? = null): String {

	val buf = StringBuilder()
	fun write(str: String, vararg args: Any) = buf.append(String.format(str, *args))

	// write the header
	if (table != null) {
		write("[$table]\n")
	}
	write("id = %s\n", id.quote())
	write("name = %s\n", name.quote())
	description?.let {
		write("description = %s\n", it.multilineQuote())
	}
	citation?.let {
		write("citation = %s\n", it.multilineQuote())
	}

	// write the fragments
	write(fragments.values.sortedBy { it.id }.toToml(table = table).toml)

	return buf.toString()
}


data class FragmentsTOML(
	val toml: String,
	val idsByFrag: Map<ConfLib.Fragment,String>
)

/**
 * Writes out a list of fragments to TOML
 */
fun List<ConfLib.Fragment>.toToml(
	/**
	 * If true, appends sequence numbers to fragment ids to avoid collisions.
	 * If false, throws an exception when an id collision is found.
	 */
	resolveIdCollisions: Boolean = false,
	/**
	 * If given, add fragments to the table, instead of the root document.
	 */
	table: String? = null
): FragmentsTOML {

	val buf = StringBuilder()
	fun write(str: String, vararg args: Any) = buf.append(String.format(str, *args))

	val fragIds = HashSet<String>()
	val idsByFrag = IdentityHashMap<ConfLib.Fragment,String>()

	for (frag in this) {

		// find the id, resolving collisions if desired
		val id = when {
			frag.id !in fragIds -> frag.id
			resolveIdCollisions -> {
				val num = fragIds
					.filter { it.startsWith(frag.id) }
					.map {
						it.substring(frag.id.length)
							.filter { it.isDigit() }
							.toIntOrNull()
							?: 1
					}
					.maxOrNull()
					?: 1
				"${frag.id}${num + 1}"
			}
			else -> throw IllegalArgumentException("multiple fragments with id: ${frag.id}")
		}
		fragIds.add(id)
		idsByFrag[frag] = id

		val tablePrefix = if (table != null) {
			"$table."
		} else {
			""
		}

		write("\n")
		write("[${tablePrefix}frag.$id]\n")
		write("name = %s\n", frag.name.quote())
		write("type = %s\n", frag.type.quote())

		// write the atoms
		write("atoms = [\n")
		for (atom in frag.atoms) {
			write("\t{ id = %2d, name = %6s, elem = %s },\n",
				atom.id,
				atom.name.quote(),
				atom.element.symbol.quote()
			)
		}
		write("]\n")

		// write the bonds
		write("bonds = [\n")
		for (bond in frag.bonds) {
			write("\t[ %2d, %2d ], # %4s - %-4s\n",
				bond.a.id,
				bond.b.id,
				bond.a.name,
				bond.b.name
			)
		}
		write("]\n")

		// write the anchors
		write("anchors = [\n")
		for (anchor in frag.anchors) {
			when (anchor) {
				is ConfLib.Anchor.Single -> {
					write("\t{ id = %2d, type = %s, bonds = [ %s ] }, # %s\n",
						anchor.id,
						"single".quote(),
						anchor.bonds.joinToString(", ") { it.id.toString() },
						anchor.bonds.joinToString(", ") { it.name }
					)
				}
				is ConfLib.Anchor.Double -> {
					write("\t{ id = %2d, type = %s, bondsa=[ %s ], bondsb=[ %s ] }, # %s; %s\n",
						anchor.id,
						"double".quote(),
						anchor.bondsa.joinToString(", ") { it.id.toString() },
						anchor.bondsb.joinToString(", ") { it.id.toString() },
						anchor.bondsa.joinToString(", ") { it.name },
						anchor.bondsb.joinToString(", ") { it.name }

					)
				}
			}
		}
		write("]\n")

		fun ConfLib.AtomPointer.name() =
			when (this) {
				is ConfLib.AtomInfo -> name
				is ConfLib.AnchorAtomPointer -> "A$index"
				else -> "?"
			}

		// write the motions
		write("motions = [\n")
		for (motion in frag.motions) {
			when (motion) {
				is ConfLib.ContinuousMotion.DihedralAngle ->
					write("\t{ id = %2d, type = %s, a = %s, b = %s, c = %s, d = %s }, # %s, %s, %s, %s\n",
						motion.id,
						"dihedral".quote(),
						motion.a.toToml(),
						motion.b.toToml(),
						motion.c.toToml(),
						motion.d.toToml(),
						motion.a.name(),
						motion.b.name(),
						motion.c.name(),
						motion.d.name()
					)
			}
		}
		write("]\n")

		// write the confs
		for ((confId, conf) in frag.confs) {

			write("[${tablePrefix}frag.$id.conf.$confId]\n")
			write("name = %s\n", conf.name.quote())
			conf.description?.let { write("description = %s\n", it.quote()) }

			write("coords = [\n")
			for ((atom, pos) in conf.coords.entries.sortedBy { (atom, _) -> atom.id }) {
				write("\t{ id = %2d, xyz = %s }, # %s\n",
					atom.id,
					pos.toToml(),
					atom.name
				)
			}
			write("]\n")

			write("anchorCoords = [\n")
			for ((anchor, coords) in conf.anchorCoords.entries.sortedBy { (anchor, _) -> anchor.id }) {
				when (coords) {
					is ConfLib.AnchorCoords.Single ->
						write("\t{ id = %2d, a = %s, b = %s, c = %s },\n",
							anchor.id,
							coords.a.toToml(),
							coords.b.toToml(),
							coords.c.toToml()
						)
					is ConfLib.AnchorCoords.Double ->
						write("\t{ id = %2d, a = %s, b = %s, c = %s, d = %s },\n",
							anchor.id,
							coords.a.toToml(),
							coords.b.toToml(),
							coords.c.toToml(),
							coords.d.toToml()
						)
				}
			}
			write("]\n")
		}
	}

	return FragmentsTOML(
		buf.toString(),
		idsByFrag
	)
}

private fun String.quote() = "'$this'"

private fun String.multilineQuote() = "'''\n$this'''"

private fun Vector3dc.toToml() =
	"[ %12.6f, %12.6f, %12.6f ]".format(x, y, z)

private fun ConfLib.AtomPointer.toToml() =
	when (this) {
		is ConfLib.AtomInfo -> id.toString()
		is ConfLib.AnchorAtomPointer -> "[ %d, %d ]".format(anchor.id, index + 1)
		else -> throw IllegalArgumentException("unkown atom pointer: $this")
	}
