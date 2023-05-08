package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.compiler.AtomPairs
import edu.duke.cs.osprey.gui.compiler.CompiledConfSpace.AtomInfo
import edu.duke.cs.osprey.gui.compiler.CompiledConfSpace.ResInfo
import edu.duke.cs.osprey.gui.compiler.ConfSpaceCompiler
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.motions.TranslationRotation
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.Proteins
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.tools.identityHashMapOf
import edu.duke.cs.osprey.tools.FileTools
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.Matcher
import io.kotest.matchers.MatcherResult
import io.kotest.matchers.collections.shouldContainExactly
import io.kotest.matchers.collections.shouldContainExactlyInAnyOrder
import io.kotest.matchers.should
import io.kotest.matchers.shouldBe
import org.joml.Vector3d


/**
 * Check each and every forcefield atom pair in a whole conf space,
 * to make sure none of them are missing, or no unexpected ones are there either
 */
@Suppress("LocalVariableName")
class TestCompiledAtomPairs : FunSpec({

	// load the amino acid conf lib
	val conflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))

	// make Vector3d.toString() do something reasonable instead of hiding all the precision
	System.setProperty("joml.format", "false")
	org.joml.internal.Options.NUMBER_FORMAT.maximumFractionDigits = 8

	fun ConfSpace.compile() = withService {
		ConfSpaceCompiler(this).run {

			// just the amber forcefield for this test, since it has more atom pairs than EEF1
			forcefields.add(Forcefield.Amber96)

			compile().run {
				waitForFinish()
				report!!.run {
					compiled
						?: error?.let { throw Error("compilation failed", it) }
						?: throw Error("no compiled")
				}
			}
		}
	}

	fun haveAtoms(vararg expected: AtomInfo) = object : Matcher<List<AtomInfo>> {
		override fun test(value: List<AtomInfo>): MatcherResult {

			val passed = value.size == expected.size && value.indices.all { i ->
				val exp = expected[i]
				val obs = value[i]
				exp.name == obs.name
				&& exp.molIndex == obs.molIndex
				&& exp.resIndex == obs.resIndex
				&& exp.pos.distance(obs.pos) <= 1e-3
			}

			fun Double.format() =
				String.format("%.3f", this)
			fun Collection<AtomInfo>.format() =
				joinToString("\n\t") { "${it.molIndex}:${it.resIndex}:${it.name}   (${it.pos.x.format()}, ${it.pos.y.format()}, ${it.pos.z.format()})" }

			return MatcherResult(
				passed,
				{
					"""
						|List should contain:
						|	${expected.toList().format()}
						|but was:
						|	${value.format()}
					""".trimMargin()
				},
				{
					"nope"
				}
			)
		}
	}

	fun List<AtomInfo>.indexOfAtom(q: AtomInfo) =
		indexOfFirst { atom ->
			atom.molIndex == q.molIndex
			&& atom.resIndex == q.resIndex
			&& atom.name == q.name
			&& atom.pos.distance(q.pos) <= 1e-9
		}
		.takeIf { it >= 0 }
		?: throw NoSuchElementException("atom info $q was not found in:\n\t${joinToString("\n\t")}")

	data class AtomKey(
		val moli: Int,
		val resi: Int,
		val atomName: String
	)

	fun List<AtomInfo>.indexOfAtom(key: AtomKey) =
		indexOfFirst { atom ->
			atom.molIndex == key.moli
			&& atom.resIndex == key.resi
			&& atom.name == key.atomName
		}
		.takeIf { it >= 0 }
		?: throw NoSuchElementException("atom $key was not found in:\n\t${joinToString("\n\t")}")

	fun List<AtomInfo>.indexOfAtom(key: Any) =
		when (key) {
			is AtomInfo -> indexOfAtom(key)
			is AtomKey -> indexOfAtom(key)
			else -> throw IllegalArgumentException("unrecognized atom key type: ${key::class}")
		}

	fun AtomInfo.format() =
		"$molIndex:$resIndex:$name"

	fun Pair<AtomInfo,AtomInfo>.format(): String {
		val (a1, a2) = this
		return "${a1.format()} <-> ${a2.format()}"
	}

	fun List<Pair<AtomInfo,AtomInfo>>.format() =
		if (isNotEmpty()) {
			joinToString("\n\t") { it.format() }
		} else {
			"(none)"
		}

	fun AtomPairs.AtomPair.matchesSameOrder(i1: Int, i2: Int) =
		(atomi1 == i1 && atomi2 == i2)

	fun AtomPairs.AtomPair.matchesEitherOrder(i1: Int, i2: Int) =
		(atomi1 == i1 && atomi2 == i2) || (atomi1 == i2 && atomi2 == i1)

	fun AtomPairs.AtomPair.matches(atoms1: List<AtomInfo>, i1: Int, atoms2: List<AtomInfo>, i2: Int) =
		if (atoms1 === atoms2) {
			matchesEitherOrder(i1, i2)
		} else {
			matchesSameOrder(i1, i2)
		}

	fun haveOnlyAtomPairs(
		atoms1: List<AtomInfo>,
		atoms2: List<AtomInfo>,
		vararg expected: Pair<Any,Any>
	) = object : Matcher<List<AtomPairs.AtomPair>> {
		override fun test(value: List<AtomPairs.AtomPair>): MatcherResult {

			// give the variable a better name for this context
			val observed = value

			// expand the expected list into a concrete list of atom pairs
			fun Any.toList(): List<Any> =
				(this as? List<Any?> ?: listOf(this)).filterNotNull()
			val expectedAtoms: List<Pair<AtomInfo,AtomInfo>> = expected
				.flatMap { (any1, any2) ->
					any1.toList().flatMap { exp1 ->
						val atom1 = atoms1[atoms1.indexOfAtom(exp1)]
						any2.toList().map { exp2 ->
							val atom2 = atoms2[atoms2.indexOfAtom(exp2)]
							atom1 to atom2
						}
					}
				}

			// gather the atoms pair we expected but didn't find
			val missing: List<Pair<AtomInfo,AtomInfo>> = expectedAtoms
				.filter { (exp1, exp2) ->
					val atomi1 = atoms1.indexOfAtom(exp1)
					val atomi2 = atoms2.indexOfAtom(exp2)
					observed.none { obs ->
						obs.matches(atoms1, atomi1, atoms2, atomi2)
					}
				}

			// gather the atoms pairs we didn't expect but found anyway
			val unexpected: List<Pair<AtomInfo,AtomInfo>> = observed
				.filter { obs ->
					expectedAtoms.none { (exp1, exp2) ->
						val atomi1 = atoms1.indexOfAtom(exp1)
						val atomi2 = atoms2.indexOfAtom(exp2)
						obs.matches(atoms1, atomi1, atoms2, atomi2)
					}
				}
				.map { pair ->
					atoms1[pair.atomi1] to atoms2[pair.atomi2]
				}

			return MatcherResult(
				missing.isEmpty() && unexpected.isEmpty(),
				{
					"""
						|List should contain:
						|	${expectedAtoms.format()}
						|but these pairs were missing:
						|	${missing.format()}
						|and these pairs weren't expected:
						|	${unexpected.format()}
					""".trimMargin()
				},
				{ "nope" }
			)
		}
	}

	fun AtomPairs.haveParams() = object : Matcher<List<AtomPairs.AtomPair>> {
		override fun test(value: List<AtomPairs.AtomPair>): MatcherResult {

			val uncachedPairs = value.filter {
				it.paramsi !in paramsCache.indices
			}

			return MatcherResult(
				uncachedPairs.isEmpty(),
				{ "Atom pairs list has ${uncachedPairs.size} entries without any forcefield params" },
				{ "nope" }
			)
		}
	}


	context("amber96") {

		test("gly-glyala") {

			// load some tiny tiny molecules
			// small enough that we can track every atom pair without too much trouble
			val (gly, glyala) = Molecule.fromOMOL(FileTools.readResource("/confSpaces/gly-glyala.omol"))
			gly as Polymer
			glyala as Polymer

			// make sure we got all the residues
			val A1 = gly.findResidueOrThrow("A", "1")
			val B2 = glyala.findResidueOrThrow("B", "2")
			val B3 = glyala.findResidueOrThrow("B", "3")
			
			// these should be all the atoms in the molecules
			val A1N   = AtomInfo(  "N", Vector3d(4.400000, -0.515000,  7.533000), 0, 0)
			val A1CA  = AtomInfo( "CA", Vector3d(5.791000, -0.751000,  7.871000), 0, 0)
			val A1C   = AtomInfo(  "C", Vector3d(6.672000,  0.451000,  7.612000), 0, 0)
			val A1O   = AtomInfo(  "O", Vector3d(7.716000,  0.674000,  8.236000), 0, 0)
			val A1OXT = AtomInfo("OXT", Vector3d(6.370663,  1.270451,  6.747051), 0, 0)
			val A1H   = AtomInfo(  "H", Vector3d(4.142885,  0.391196,  7.168582), 0, 0)
			val A1HA2 = AtomInfo("HA2", Vector3d(5.846169, -1.009957,  8.928354), 0, 0)
			val A1HA3 = AtomInfo("HA3", Vector3d(6.149443, -1.588911,  7.273068), 0, 0)

			val B2N   = AtomInfo(  "N", Vector3d( 6.124000, 2.886000, 10.277000), 1, 1)
			val B2CA  = AtomInfo( "CA", Vector3d( 6.834000, 2.327000, 11.416000), 1, 1)
			val B2C   = AtomInfo(  "C", Vector3d( 8.328000, 2.590000, 11.400000), 1, 1)
			val B2O   = AtomInfo(  "O", Vector3d( 8.953000, 2.898000, 12.415000), 1, 1)
			val B2H1  = AtomInfo( "H1", Vector3d( 6.797028, 3.352733,  9.685989), 1, 1)
			val B2H2  = AtomInfo( "H2", Vector3d( 5.498584, 3.577379, 10.665523), 1, 1)
			val B2H3  = AtomInfo( "H3", Vector3d( 5.750795, 2.090706,  9.778678), 1, 1)
			val B2HA2 = AtomInfo("HA2", Vector3d( 6.406782, 2.762838, 12.319122), 1, 1)
			val B2HA3 = AtomInfo("HA3", Vector3d( 6.663359, 1.250441, 11.416928), 1, 1)

			val B3N   = AtomInfo(  "N", Vector3d( 8.953000, 2.496000, 10.228000), 1, 2)
			val B3CA  = AtomInfo( "CA", Vector3d(10.391000, 2.781000, 10.119000), 1, 2)
			val B3C   = AtomInfo(  "C", Vector3d(10.720000, 4.212000, 10.519000), 1, 2)
			val B3O   = AtomInfo(  "O", Vector3d(11.712000, 4.502000, 11.209000), 1, 2)
			val B3CB  = AtomInfo( "CB", Vector3d(10.850000, 2.513000,  8.682000), 1, 2)
			val B3OXT = AtomInfo("OXT", Vector3d( 9.996575, 5.140503, 10.165478), 1, 2)
			val B3H   = AtomInfo(  "H", Vector3d( 8.460974, 2.192115,  9.399952), 1, 2)
			val B3HA  = AtomInfo( "HA", Vector3d(10.939808, 2.119261, 10.789082), 1, 2)
			val B3HB1 = AtomInfo("HB1", Vector3d(10.315832, 3.174595,  8.000049), 1, 2)
			val B3HB2 = AtomInfo("HB2", Vector3d(11.921269, 2.697540,  8.601822), 1, 2)
			val B3HB3 = AtomInfo("HB3", Vector3d(10.639905, 1.475934,  8.420363), 1, 2)

			// make sure we got all the atoms
			A1.atoms.map { it.name } shouldContainExactlyInAnyOrder listOf(
				A1N, A1CA, A1C, A1O, A1OXT, A1H, A1HA2, A1HA3
			).map { it.name }
			B2.atoms.map { it.name } shouldContainExactlyInAnyOrder listOf(
				B2N, B2CA, B2C, B2O, B2H1, B2H2, B2H3, B2HA2, B2HA3
			).map { it.name }
			B3.atoms.map { it.name } shouldContainExactlyInAnyOrder listOf(
				B3N, B3CA, B3C, B3O, B3CB, B3OXT, B3H, B3HA, B3HB1, B3HB2, B3HB3
			).map { it.name }

			// build a conf space
			val confSpace = ConfSpace(listOf(
				MoleculeType.Protein to gly,
				MoleculeType.Protein to glyala,
			)).apply {

				name = "test"

				conflibs.add(conflib)

				// define three design positions, one for each residue
				val posA1 = Proteins.makeDesignPosition(gly, A1, "A1 GLY")
				val posB2 = Proteins.makeDesignPosition(glyala, B2, "B2 GLY")
				val posB3 = Proteins.makeDesignPosition(glyala, B3, "B3 ALA")
				designPositionsByMol[gly] = mutableListOf(posA1)
				designPositionsByMol[glyala] = mutableListOf(posB2, posB3)

				// add wild-type conformations
				for (pos in listOf(posA1, posB2, posB3)) {
					addWildTypeConformation(pos)
				}

				// let A1 mutate to ALA
				addConformationsFromLibraries(posA1, "ALA")

				// let B2 mutate to ALA
				addConformationsFromLibraries(posB2, "ALA")

				// keep B3 at ALA, allow flexibility
				addConformationsFromLibraries(posB3, "ALA")

				// let the free glycine translate and rotate
				molMotions[gly] = mutableListOf(
					TranslationRotation.MolDescription(
						gly,
						2.0,
						5.0
					)
				)
			}

			// compile it
			val compiled = confSpace.compile()
			
			// check the easy things
			compiled.name shouldBe confSpace.name
			compiled.molInfos.map { it.name } shouldBe listOf(gly.name, glyala.name)

			// make some tools to refer to atoms easily			
			val molisByRes = identityHashMapOf(
				A1 to 0,
				B2 to 1,
				B3 to 1
			)
			val resisByRes = identityHashMapOf(
				A1 to 0,
				B2 to 1,
				B3 to 2
			)
			operator fun Polymer.Residue.invoke(a: String) =
				AtomKey(
					molisByRes.getValue(this),
					resisByRes.getValue(this),
					a
				)

			// check the residues
			compiled.resInfos.size shouldBe 3
			compiled.resInfos[0] shouldBe ResInfo("A", "1", "GLY", 0)
			compiled.resInfos[1] shouldBe ResInfo("B", "2", "GLY", 0)
			compiled.resInfos[2] shouldBe ResInfo("B", "3", "ALA", 1)

			// check the static atoms
			// ie, atoms that are fixed (don't get changed by conformations),
			// and whose forcefield params aren't affected by conformation changes
			compiled.staticAtoms.shouldContainExactly(listOf(

				// free gly should be just NH group
				// amber96 uses different charges for ALA and GLY on the rest of the backbone
				A1N, A1H,

				// the gly in glyala should be just the backbone, except CA
				// amber96 uses different charges for the NH3 group
				B2C, B2O,

				// the ala in glyala should be the whole backbone, since nothing mutates here
				B3N, B3CA, B3C, B3O, B3OXT
			))

			// check the positions
			compiled.positions.size shouldBe 3
			compiled.positions[0].apply {

				name shouldBe "A1 GLY"
				wildType shouldBe "GLY"

				val dynamicFixedAtoms = listOf(
					A1C, A1CA, A1O, A1OXT
				).map { it.name }

				// check the fragments
				fragments.size shouldBe 2
				fragments[0].apply {
					name shouldBe "ALAn"
					atomNames shouldContainExactlyInAnyOrder dynamicFixedAtoms + listOf(
						"CB", "HA", "HB1", "HB2", "HB3" // sidechain atoms
					)
				}
				fragments[1].apply {
					name shouldBe "wt-Chain_A-A1_GLY"
					atomNames shouldContainExactlyInAnyOrder dynamicFixedAtoms + listOf(
						"HA2", "HA3" // sidechain atoms
					)
				}

				// check the confs
				confs.size shouldBe 2
				confs[0].apply {
					id shouldBe "ALAn:ALA"
					atoms should haveAtoms(
						A1CA, A1C, A1O, A1OXT,
						AtomInfo( "HA", Vector3d( 5.849, -0.971, 8.947), 0, 0),
						AtomInfo( "CB", Vector3d( 6.327, -1.945, 7.090), 0, 0),
						AtomInfo("HB1", Vector3d( 7.381, -2.114, 7.355), 0, 0),
						AtomInfo("HB2", Vector3d( 5.740, -2.841, 7.339), 0, 0),
						AtomInfo("HB3", Vector3d( 6.248, -1.744, 6.011), 0, 0)
					)
					fragments.getOrNull(fragIndex)?.name shouldBe "ALAn"
				}
				confs[1].apply {
					id shouldBe "wt-Chain_A-A1_GLY:conf1"
					atoms should haveAtoms(
						A1CA, A1C, A1O, A1OXT,
						AtomInfo( "HA2", Vector3d( 5.846, -1.010, 8.928), 0, 0),
						AtomInfo( "HA3", Vector3d( 6.149, -1.589, 7.273), 0, 0)
					)
					fragments.getOrNull(fragIndex)?.name shouldBe "wt-Chain_A-A1_GLY"
				}
			}
			compiled.positions[1].apply {

				name shouldBe "B2 GLY"
				wildType shouldBe "GLY"

				val dynamicFixedAtoms = listOf(
					B2N, B2H1, B2H2, B2H3, B2CA
				).map { it.name }

				// check the fragments
				fragments.size shouldBe 2
				fragments[0].apply {
					name shouldBe "ALAn"
					atomNames shouldContainExactlyInAnyOrder dynamicFixedAtoms + listOf(
						"CB", "HA", "HB1", "HB2", "HB3" // sidechain atoms
					)
				}
				fragments[1].apply {
					name shouldBe "wt-Chain_B-B2_GLY"
					atomNames shouldContainExactlyInAnyOrder dynamicFixedAtoms + listOf(
						"HA2", "HA3" // sidechain atoms
					)
				}

				// check the confs
				confs.size shouldBe 2
				confs[0].apply {
					id shouldBe "ALAn:ALA"
					atoms should haveAtoms(
						B2N, B2CA, B2H1, B2H2, B2H3,
						AtomInfo( "HA", Vector3d( 6.440, 2.798, 12.329), 1, 1),
						AtomInfo( "CB", Vector3d( 6.604, 0.822, 11.493), 1, 1),
						AtomInfo("HB1", Vector3d( 7.145, 0.412, 12.358), 1, 1),
						AtomInfo("HB2", Vector3d( 5.529, 0.620, 11.605), 1, 1),
						AtomInfo("HB3", Vector3d( 6.971, 0.347, 10.571), 1, 1)
					)
					fragments.getOrNull(fragIndex)?.name shouldBe "ALAn"
				}
				confs[1].apply {
					id shouldBe "wt-Chain_B-B2_GLY:conf1"
					atoms should haveAtoms(
						B2N, B2CA, B2H1, B2H2, B2H3,
						AtomInfo( "HA2", Vector3d( 6.407, 2.763, 12.319), 1, 1),
						AtomInfo( "HA3", Vector3d( 6.663, 1.250, 11.417), 1, 1)
					)
					fragments.getOrNull(fragIndex)?.name shouldBe "wt-Chain_B-B2_GLY"
				}
			}
			compiled.positions[2].apply {

				name shouldBe "B3 ALA"
				wildType shouldBe "ALA"

				val dynamicFixedAtoms = listOf(
					B3H
				).map { it.name }

				// check the fragments
				fragments.size shouldBe 2
				fragments[0].apply {
					name shouldBe "ALA"
					atomNames shouldContainExactlyInAnyOrder dynamicFixedAtoms + listOf(
						"CB", "HA", "HB1", "HB2", "HB3" // sidechain atoms
					)
				}
				fragments[1].apply {
					name shouldBe "wt-Chain_B-B3_ALA"
					atomNames shouldContainExactlyInAnyOrder dynamicFixedAtoms + listOf(
						"CB", "HA", "HB1", "HB2", "HB3" // sidechain atoms
					)
				}

				// check the confs
				confs.size shouldBe 2
				confs[0].apply {
					id shouldBe "ALA:ALA"
					atoms should haveAtoms(
						AtomInfo(  "H", Vector3d(  8.335, 2.379,  9.441), 1, 2), // this ALA conf moved the H too
						AtomInfo( "HA", Vector3d( 10.926, 2.111, 10.808), 1, 2),
						AtomInfo( "CB", Vector3d( 10.873, 2.520,  8.697), 1, 2),
						AtomInfo("HB1", Vector3d( 11.949, 2.736,  8.628), 1, 2),
						AtomInfo("HB2", Vector3d( 10.695, 1.466,  8.436), 1, 2),
						AtomInfo("HB3", Vector3d( 10.324, 3.168,  7.998), 1, 2)
					)
					fragments.getOrNull(fragIndex)?.name shouldBe "ALA"
				}
				confs[1].apply {
					id shouldBe "wt-Chain_B-B3_ALA:conf1"
					atoms should haveAtoms(
						AtomInfo( "CB", Vector3d(10.850, 2.513,  8.682), 1, 2),
						B3H,
						AtomInfo( "HA", Vector3d(10.940, 2.119, 10.789), 1, 2),
						AtomInfo("HB1", Vector3d(10.316, 3.175,  8.000), 1, 2),
						AtomInfo("HB2", Vector3d(11.921, 2.698,  8.602), 1, 2),
						AtomInfo("HB3", Vector3d(10.640, 1.476,  8.420), 1, 2)
					)
					fragments.getOrNull(fragIndex)?.name shouldBe "wt-Chain_B-B3_ALA"
				}
			}

			for (pos in compiled.positions) {
				for (conf in pos.confs) {

					// none of the confs should have internal energies, since we're only using Amber96
					conf.internalEnergies.size shouldBe 1
					conf.internalEnergies[0] shouldBe 0.0

					// or motions either, since it's all ala and gly
					conf.motions.size shouldBe 0
				}
			}

			// TODO: make sure all the present atom pairs actually have params in the cache

			// check the atom pairs
			compiled.atomPairs.size shouldBe 1
			compiled.atomPairs[0].apply { // amber96

				// static atom pairs, look for just:
				// * atom pairs within A
				// * between A and B
				// but not within B, since those interactions get counted up in the static energy
				static should haveOnlyAtomPairs(
					atoms1 = compiled.staticAtoms,
					atoms2 = compiled.staticAtoms,

					listOf(A1N, A1H) to listOf(B2C, B2O, B3N, B3CA, B3C, B3O, B3OXT)
				)
				static should haveParams()

				// A1 GLY -> ALAn
				posStatic[0, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[0].confs[0].atoms,
					atoms2 = compiled.staticAtoms,

					// all 1-4 and farther pairs within A
					listOf(
						A1("HA"),
						A1("CB"), A1("C")
					) to A1H,
					listOf(
						A1("HB1"), A1("HB2"), A1("HB3"),
						A1("O"), A1("OXT")
					) to listOf(A1N, A1H),

					// all pairs between A and B, since they're all nonbonded
					listOf(
						A1("CA"), A1("C"), A1("O"), A1("OXT"), A1("HA"),
						A1("CB"), A1("HB1"), A1("HB2"), A1("HB3")
					) to listOf(B2C, B2O, B3N, B3CA, B3C, B3O, B3OXT)
				)
				posStatic[0, 0] should haveParams()
				pos[0, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[0].confs[0].atoms,
					atoms2 = compiled.positions[0].confs[0].atoms,

					// all 1-4 and farther pairs within A
					listOf(
						A1("HA"), A1("CB")
					) to listOf(
						A1("O"), A1("OXT")
					),
					listOf(
						A1("HB1"), A1("HB2"), A1("HB3")
					) to listOf(
						A1("HA"), A1("C"), A1("O"), A1("OXT")
					)
				)
				pos[0, 0] should haveParams()

				// A1 GLY -> wt
				posStatic[0, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[0].confs[1].atoms,
					atoms2 = compiled.staticAtoms,

					// all 1-4 and farther pairs within A
					listOf(A1HA2, A1HA3, A1C) to A1H,
					listOf(A1O, A1OXT) to listOf(A1H, A1N),

					// all pairs between A and B, since they're all nonbonded
					listOf(
						A1HA2, A1HA3, A1CA, A1C, A1O, A1OXT
					) to listOf(
						B2C, B2O, B3N, B3CA, B3C, B3O, B3OXT
					)
				)
				posStatic[0, 1] should haveParams()
				pos[0, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[0].confs[1].atoms,
					atoms2 = compiled.positions[0].confs[1].atoms,

					// all 1-4 and farther pairs within A
					listOf(A1HA2, A1HA3) to listOf(A1O, A1OXT)
				)
				pos[0, 1] should haveParams()

				// B2 GLY -> ALAn
				posStatic[1, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[0].atoms,
					atoms2 = compiled.staticAtoms,

					listOf(
						B2("N"), B2("H1"), B2("H2"), B2("H3"), B2("HA"),
						B2("CA"), B2("CB"), B2("HB1"), B2("HB2"), B2("HB3")
					) to listOf(A1N, A1H),

					listOf(
						B2("H1"), B2("H2"), B2("H3"),
						B2("HB1"), B2("HB2"), B2("HB3")
					) to listOf(B2C, B2O, B3N, B3CA, B3C, B3O, B3OXT),

					listOf(
						B2("N"), B2("CB"), B2("HA")
					) to listOf(B2O, B3N, B3CA, B3C, B3O, B3OXT),

					B2("CA") to listOf(B3CA, B3C, B3O, B3OXT)
				)
				posStatic[1, 0] should haveParams()
				pos[1, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[0].atoms,
					atoms2 = compiled.positions[1].confs[0].atoms,

					listOf(
						B2("H1"), B2("H2"), B2("H3")
					) to listOf(
						B2("CB"), B2("HB1"), B2("HB2"), B2("HB3"), B2("HA")
					),

					listOf(
						B2("N"), B2("HA")
					) to listOf(
						B2("HB1"), B2("HB2"), B2("HB3")
					)
				)
				pos[1, 0] should haveParams()

				// B2 GLY -> wt
				posStatic[1, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[1].atoms,
					atoms2 = compiled.staticAtoms,

					listOf(B2N, B2H1, B2H2, B2H3, B2CA, B2HA2, B2HA3) to listOf(A1N, A1H),
					listOf(B2H1, B2H2, B2H3) to listOf(B2C, B2O, B3N, B3CA, B3C, B3O, B3OXT),
					listOf(B2N, B2HA2, B2HA3) to listOf(B2O, B3N, B3CA, B3C, B3O, B3OXT),
					listOf(B2CA) to listOf(B3CA, B3C, B3O, B3OXT)
				)
				posStatic[1, 1] should haveParams()
				pos[1, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[1].atoms,
					atoms2 = compiled.positions[1].confs[1].atoms,

					listOf(B2H1, B2H2, B2H3) to listOf(B2HA2, B2HA3)
				)
				pos[1, 1] should haveParams()

				// B3 ALA -> ALA
				posStatic[2, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[0].atoms,
					atoms2 = compiled.staticAtoms,

					listOf(
						B3("H"), B3("HA"),
						B3("CB"), B3("HB1"), B3("HB2"), B3("HB3")
					) to listOf(A1N, A1H),

					listOf(
						B3("H")
					) to listOf(B2O, B3C, B3O, B3OXT),

					listOf(
						B3("CB"), B3("HA")
					) to listOf(B2C, B2O, B3O, B3OXT),

					listOf(
						B3("HB1"), B3("HB2"), B3("HB3")
					) to listOf(B2C, B2O, B3N, B3C, B3O, B3OXT)
				)
				posStatic[2, 0] should haveParams()
				pos[2, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[0].atoms,
					atoms2 = compiled.positions[2].confs[0].atoms,

					listOf(
						B3("H")
					) to listOf(
						B3("HA"), B3("CB"), B3("HB1"), B3("HB2"), B3("HB3")
					),

					listOf(
						B3("HA")
					) to listOf(
						B3("HB1"), B3("HB2"), B3("HB3")
					)
				)
				pos[2, 0] should haveParams()

				// B2 GLY > ALAn <-> A1 GLY > ALAn
				posPos[1, 0, 0, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[0].atoms,
					atoms2 = compiled.positions[0].confs[0].atoms,

					// different molecules, just all pairs
					listOf(
						B2("N"), B2("H1"), B2("H2"), B2("H3"),
						B2("CA"), B2("HA"),
						B2("CB"), B2("HB1"), B2("HB2"), B2("HB3")
					) to listOf(
						A1("CA"), A1("C"), A1("O"), A1("OXT"),
						A1("HA"), A1("CB"), A1("HB1"), A1("HB2"), A1("HB3")
					)
				)
				posPos[1, 0, 0, 0] should haveParams()

				// B2 GLY > ALAn <-> A1 GLY > wt
				posPos[1, 0, 0, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[0].atoms,
					atoms2 = compiled.positions[0].confs[1].atoms,

					// different molecules, just all pairs
					listOf(
						B2("N"), B2("H1"), B2("H2"), B2("H3"),
						B2("CA"), B2("HA"),
						B2("CB"), B2("HB1"), B2("HB2"), B2("HB3")
					) to listOf(A1CA, A1C, A1O, A1OXT, A1HA2, A1HA3)
				)
				posPos[1, 0, 0, 1] should haveParams()

				// B2 GLY > wt <-> A1 GLY > ALA
				posPos[1, 1, 0, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[1].atoms,
					atoms2 = compiled.positions[0].confs[0].atoms,

					// different molecules, just all pairs
					listOf(B2N, B2H1, B2H2, B2H3, B2CA, B2HA2, B2HA3) to listOf(
						A1("CA"), A1("C"), A1("O"), A1("OXT"),
						A1("HA"), A1("CB"), A1("HB1"), A1("HB2"), A1("HB3")
					)
				)
				posPos[1, 1, 0, 0] should haveParams()

				// B2 GLY > wt <-> A1 GLY > wt
				posPos[1, 1, 0, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[1].confs[1].atoms,
					atoms2 = compiled.positions[0].confs[1].atoms,

					// different molecules, just all pairs
					listOf(B2N, B2H1, B2H2, B2H3, B2CA, B2HA2, B2HA3) to listOf(A1CA, A1C, A1O, A1OXT, A1HA2, A1HA3)
				)
				posPos[1, 1, 0, 1] should haveParams()

				// B3 ALA > ALA <-> A1 GLY -> ALAn
				posPos[2, 0, 0, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[0].atoms,
					atoms2 = compiled.positions[0].confs[0].atoms,

					// different molecules, just all pairs
					listOf(
						B3("H"), B3("HA"), B3("CB"), B3("HB1"), B3("HB2"), B3("HB3")
					) to listOf(
						A1("CA"), A1("C"), A1("O"), A1("OXT"),
						A1("HA"), A1("CB"), A1("HB1"), A1("HB2"), A1("HB3")
					)
				)
				posPos[2, 0, 0, 0] should haveParams()

				// B3 ALA > ALA <-> A1 GLY -> wt
				posPos[2, 0, 0, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[0].atoms,
					atoms2 = compiled.positions[0].confs[1].atoms,

					// different molecules, just all pairs
					listOf(
						B3("H"), B3("HA"), B3("CB"), B3("HB1"), B3("HB2"), B3("HB3")
					) to listOf(A1CA, A1C, A1O, A1OXT, A1HA2, A1HA3)
				)
				posPos[2, 0, 0, 1] should haveParams()

				// B3 ALA > wt <-> A1 GLY -> ALAn
				posPos[2, 1, 0, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[1].atoms,
					atoms2 = compiled.positions[0].confs[0].atoms,

					// different molecules, just all pairs
					listOf(B3H, B3HA, B3CB, B3HB1, B3HB2, B3HB3) to listOf(
						A1("CA"), A1("C"), A1("O"), A1("OXT"),
						A1("HA"), A1("CB"), A1("HB1"), A1("HB2"), A1("HB3")
					)
				)
				posPos[2, 1, 0, 0] should haveParams()

				// B3 ALA > wt <-> A1 GLY -> wt
				posPos[2, 1, 0, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[1].atoms,
					atoms2 = compiled.positions[0].confs[1].atoms,

					// different molecules, just all pairs
					listOf(B3H, B3HA, B3CB, B3HB1, B3HB2, B3HB3) to listOf(A1CA, A1C, A1O, A1OXT, A1HA2, A1HA3)
				)
				posPos[2, 1, 0, 1] should haveParams()

				// B3 ALA > ALA <-> B2 GLY > ALAn
				posPos[2, 0, 1, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[0].atoms,
					atoms2 = compiled.positions[1].confs[0].atoms,

					// just far enough apart, all pairs
					listOf(
						B3("H"), B3("HA"), B3("CB"), B3("HB1"), B3("HB2"), B3("HB3")
					) to listOf(
						B2("N"), B2("H1"), B2("H2"), B2("H3"),
						B2("CA"), B2("HA"), B2("CB"), B2("HB1"), B2("HB2"), B2("HB3")
					)
				)
				posPos[2, 0, 1, 0] should haveParams()

				// B3 ALA > ALA <-> B2 GLY > wt
				posPos[2, 0, 1, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[0].atoms,
					atoms2 = compiled.positions[1].confs[1].atoms,

					// just far enough apart, all pairs
					listOf(
						B3("H"), B3("HA"), B3("CB"), B3("HB1"), B3("HB2"), B3("HB3")
					) to listOf(B2N, B2H1, B2H2, B2H3, B2CA, B2HA2, B2HA3)
				)
				posPos[2, 0, 1, 1] should haveParams()

				// B3 ALA > wt <-> B2 GLY > ALAn
				posPos[2, 1, 1, 0] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[1].atoms,
					atoms2 = compiled.positions[1].confs[0].atoms,

					// just far enough apart, all pairs
					listOf(B3H, B3HA, B3CB, B3HB1, B3HB2, B3HB3) to listOf(
						B2("N"), B2("H1"), B2("H2"), B2("H3"),
						B2("CA"), B2("HA"), B2("CB"), B2("HB1"), B2("HB2"), B2("HB3")
					)
				)
				posPos[2, 1, 1, 0] should haveParams()

				// B3 ALA > wt <-> B2 GLY > wt
				posPos[2, 1, 1, 1] should haveOnlyAtomPairs(
					atoms1 = compiled.positions[2].confs[1].atoms,
					atoms2 = compiled.positions[1].confs[1].atoms,

					// just far enough apart, all pairs
					listOf(B3H, B3HA, B3CB, B3HB1, B3HB2, B3HB3) to listOf(B2N, B2H1, B2H2, B2H3, B2CA, B2HA2, B2HA3)
				)
				posPos[2, 1, 1, 1] should haveParams()
			}
		}
	}
})
