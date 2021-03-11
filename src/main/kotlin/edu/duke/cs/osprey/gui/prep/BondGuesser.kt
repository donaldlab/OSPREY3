package edu.duke.cs.osprey.gui.prep

import cuchaz.kludge.tools.square
import edu.duke.cs.osprey.molscope.molecule.*
import java.util.*
import kotlin.collections.ArrayList


/**
 * Attempts to guess where the bonds should be based on
 * known elemental valences and covalent bond lengths.
 *
 * Presents guessed bonds to the caller (rather than automatically
 * adding them to the molecule), so the caller can decide what to do next.
 */
fun Molecule.guessBonds(
	guesser: BondGuesser = BondGuesser(),
	searchDist: Double = 5.0
): List<BondGuesser.Bond> {

	fun Atom.props() = guesser.getProps(this)
	val searchDistSq = searchDist.square()

	// build a tree of the atoms so we can do nearest-neighbor queries efficiently
	val tree = atoms.toTree()
	
	val guessedBonds = ArrayList<BondGuesser.Bond>()

	for (atom1 in atoms) {
		val props1 = atom1.props()

		// find the closest neighbors, up to the number allowed by the max valence
		val nearbyAtoms = tree.nearest(atom1)
		while (bonds.bondedAtoms(atom1).size < props1.maxValence && nearbyAtoms.hasNext()) {

			val (atom2, distSq) = nearbyAtoms.next()
			val props2 = atom2.props()

			// if we get too far away, stop searching
			if (distSq > searchDistSq) {
				break
			}

			// if the atom is out of range, skip it
			if (distSq < guesser.minDistSq) {
				continue
			}
			val covDistSq = (props1.bondRadius + props2.bondRadius + guesser.tolerance).square()
			if (distSq > covDistSq) {
				continue
			}

			// atom is in range, make the bond
			guessedBonds.add(BondGuesser.Bond(atom1, atom2))
		}
	}

	return guessedBonds
}

fun Atom.inCovalentRange(
	a2: Atom,
	guesser: BondGuesser = BondGuesser()
): Boolean {

	val a1 = this
	val distSq = a1.pos.distanceSquared(a2.pos)

	// check the lower bound first
	if (distSq < guesser.minDistSq) {
		return false
	}

	fun Atom.props() = guesser.getProps(this)

	// check check the upper bound
	return distSq <= (a1.props().bondRadius + a2.props().bondRadius + guesser.tolerance).square()
}

fun Atom.minCovalentDist(
	guesser: BondGuesser = BondGuesser()
): Double =
	guesser.minDist

fun Atom.maxCovalentDist(
	a2: Atom,
	guesser: BondGuesser = BondGuesser()
): Double {

	val a1 = this
	fun Atom.props() = guesser.getProps(this)

	return a1.props().bondRadius + a2.props().bondRadius + guesser.tolerance
}

fun Atom.covalentRange(
	a2: Atom,
	guesser: BondGuesser = BondGuesser()
): ClosedFloatingPointRange<Double> =
	minCovalentDist() .. maxCovalentDist(a2)


class BondGuesser(
	val tolerance: Double = defaultTolerance,
	val minDist: Double = defaultMinDist,
	val elementProps: EnumMap<Element,ElementProps> = defaultElementProps
) {

	data class Bond(val a1: Atom, val a2: Atom)

	data class ElementProps(val bondRadius: Double, val maxValence: Int)

	fun getProps(element: Element) =
		elementProps[element]
		?: throw NoSuchElementException("no properties for element $element")

	fun getProps(atom: Atom) = getProps(atom.element)

	val minDistSq = minDist.square()

	companion object {

		const val defaultMinDist = 0.4
		const val defaultTolerance = 0.45

		/* covalent radii from:
			Covalent radii revisited
			Beatriz Cordero, et al
			Dalton Transactions, 2008
			https://doi.org/10.1039/b801115j
		*/
		val defaultElementProps =
			EnumMap<Element,ElementProps>(Element::class.java).apply {
				this[Element.Hydrogen] = ElementProps(0.31, 1)
				this[Element.Helium] = ElementProps(0.28, 0)
				this[Element.Lithium] = ElementProps(1.28, 1)
				this[Element.Beryllium] = ElementProps(0.96, 2)
				this[Element.Boron] = ElementProps(0.84, 4)
				this[Element.Carbon] = ElementProps(0.76, 4)
				this[Element.Nitrogen] = ElementProps(0.71, 4)
				this[Element.Oxygen] = ElementProps(0.66, 2)
				this[Element.Fluorine] = ElementProps(0.57, 1)
				this[Element.Neon] = ElementProps(0.58, 0)
				this[Element.Sodium] = ElementProps(1.66, 1)
				this[Element.Magnesium] = ElementProps(1.41, 2)
				this[Element.Aluminium] = ElementProps(1.21, 6)
				this[Element.Silicon] = ElementProps(1.11, 6)
				this[Element.Phosphorus] = ElementProps(1.07, 6)
				this[Element.Sulfur] = ElementProps(1.05, 6)
				this[Element.Chlorine] = ElementProps(1.02, 1)
				this[Element.Argon] = ElementProps(1.06, 0)
				this[Element.Potassium] = ElementProps(2.03, 1)
				this[Element.Calcium] = ElementProps(1.76, 2)
				this[Element.Scandium] = ElementProps(1.7, 6)
				this[Element.Titanium] = ElementProps(1.6, 6)
				this[Element.Vanadium] = ElementProps(1.53, 6)
				this[Element.Chromium] = ElementProps(1.39, 6)
				this[Element.Manganese] = ElementProps(1.39, 8)
				this[Element.Iron] = ElementProps(1.32, 6)
				this[Element.Cobalt] = ElementProps(1.26, 6)
				this[Element.Nickel] = ElementProps(1.24, 6)
				this[Element.Copper] = ElementProps(1.32, 6)
				this[Element.Zinc] = ElementProps(1.22, 6)
				this[Element.Gallium] = ElementProps(1.22, 3)
				this[Element.Germanium] = ElementProps(1.2, 4)
				this[Element.Arsenic] = ElementProps(1.19, 3)
				this[Element.Selenium] = ElementProps(1.2, 2)
				this[Element.Bromine] = ElementProps(1.2, 1)
				this[Element.Krypton] = ElementProps(1.16, 0)
				this[Element.Rubidium] = ElementProps(2.2, 1)
				this[Element.Strontium] = ElementProps(1.95, 2)
				this[Element.Yttrium] = ElementProps(1.9, 6)
				this[Element.Zirconium] = ElementProps(1.75, 6)
				this[Element.Niobium] = ElementProps(1.64, 6)
				this[Element.Molybdenum] = ElementProps(1.54, 6)
				this[Element.Technetium] = ElementProps(1.47, 6)
				this[Element.Ruthenium] = ElementProps(1.46, 6)
				this[Element.Rhodium] = ElementProps(1.42, 6)
				this[Element.Palladium] = ElementProps(1.39, 6)
				this[Element.Silver] = ElementProps(1.45, 6)
				this[Element.Cadmium] = ElementProps(1.44, 6)
				this[Element.Indium] = ElementProps(1.42, 3)
				this[Element.Tin] = ElementProps(1.39, 4)
				this[Element.Antimony] = ElementProps(1.39, 3)
				this[Element.Tellurium] = ElementProps(1.38, 2)
				this[Element.Iodine] = ElementProps(1.39, 1)
				this[Element.Xenon] = ElementProps(1.4, 0)
				this[Element.Caesium] = ElementProps(2.44, 1)
				this[Element.Barium] = ElementProps(2.15, 2)
				this[Element.Lanthanum] = ElementProps(2.07, 12)
				this[Element.Cerium] = ElementProps(2.04, 6)
				this[Element.Praseodymium] = ElementProps(2.03, 6)
				this[Element.Neodymium] = ElementProps(2.01, 6)
				this[Element.Promethium] = ElementProps(1.99, 6)
				this[Element.Samarium] = ElementProps(1.98, 6)
				this[Element.Europium] = ElementProps(1.98, 6)
				this[Element.Gadolinium] = ElementProps(1.96, 6)
				this[Element.Terbium] = ElementProps(1.94, 6)
				this[Element.Dysprosium] = ElementProps(1.92, 6)
				this[Element.Holmium] = ElementProps(1.92, 6)
				this[Element.Erbium] = ElementProps(1.89, 6)
				this[Element.Thulium] = ElementProps(1.9, 6)
				this[Element.Ytterbium] = ElementProps(1.87, 6)
				this[Element.Lutetium] = ElementProps(1.87, 6)
				this[Element.Hafnium] = ElementProps(1.75, 6)
				this[Element.Tantalum] = ElementProps(1.7, 6)
				this[Element.Tungsten] = ElementProps(1.62, 6)
				this[Element.Rhenium] = ElementProps(1.51, 6)
				this[Element.Osmium] = ElementProps(1.44, 6)
				this[Element.Iridium] = ElementProps(1.41, 6)
				this[Element.Platinum] = ElementProps(1.36, 6)
				this[Element.Gold] = ElementProps(1.36, 6)
				this[Element.Mercury] = ElementProps(1.32, 6)
				this[Element.Thallium] = ElementProps(1.45, 3)
				this[Element.Lead] = ElementProps(1.46, 4)
				this[Element.Bismuth] = ElementProps(1.48, 3)
				this[Element.Polonium] = ElementProps(1.4, 2)
				this[Element.Astatine] = ElementProps(1.5, 1)
				this[Element.Radon] = ElementProps(1.5, 0)
				this[Element.Francium] = ElementProps(2.6, 1)
				this[Element.Radium] = ElementProps(2.21, 2)
				this[Element.Actinium] = ElementProps(2.15, 6)
				this[Element.Thorium] = ElementProps(2.06, 6)
				this[Element.Protactinium] = ElementProps(2.0, 6)
				this[Element.Uranium] = ElementProps(1.96, 6)
				this[Element.Neptunium] = ElementProps(1.9, 6)
				this[Element.Plutonium] = ElementProps(1.87, 6)
				this[Element.Americium] = ElementProps(1.8, 6)
				this[Element.Curium] = ElementProps(1.69, 6)
			}
	}
}
