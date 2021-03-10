package edu.duke.cs.osprey.molscope.molecule

import edu.duke.cs.osprey.SharedSpec
import io.kotlintest.shouldBe


class TestMoleculeSearch : SharedSpec({

	val c1 = Atom(Element.Carbon, "C1", 6.778, 10.510, 20.665)
	val c2 = Atom(Element.Carbon, "C2", 5.994, 9.710, 19.842)
	val c3 = Atom(Element.Carbon, "C3", 6.562, 9.055, 18.751)
	val c4 = Atom(Element.Carbon, "C4", 7.916, 9.259, 18.444)
	val c5 = Atom(Element.Carbon, "C5", 8.711, 10.120, 19.210)
	val c6 = Atom(Element.Carbon, "C6", 8.128, 10.734, 20.335)
	val c = Atom(Element.Carbon, "C", 6.244, 11.152, 21.851)
	val n1 = Atom(Element.Nitrogen, "N1", 7.014, 11.257, 22.910)
	val n2 = Atom(Element.Nitrogen, "N2", 4.965, 11.590, 21.821)

	val benzamidine = Molecule("Benzamidine").apply {

		atoms.addAll(listOf(c1, c2, c3, c4, c5, c6, c, n1, n2))

		bonds.add(c1, c2)
		bonds.add(c1, c6)
		bonds.add(c1, c)
		bonds.add(c2, c3)
		bonds.add(c3, c4)
		bonds.add(c4, c5)
		bonds.add(c5, c6)
		bonds.add(c, n1)
		bonds.add(c, n2)
	}

	fun Iterator<Molecule.Searched>.shouldHaveNext(atom: Atom, dist: Int) {
		hasNext() shouldBe true
		next() shouldBe Molecule.Searched(atom, dist)
	}

	fun Iterator<Molecule.Searched>.shouldBeDone() {
		hasNext() shouldBe false
	}

	group("dfs") {

		test("full order") {

			benzamidine.dfs(
				source = c1,
				visitSource = true
			).iterator().apply {

				shouldHaveNext(c1, 0)
				shouldHaveNext(c, 1)
				shouldHaveNext(n1, 2)
				shouldHaveNext(n2, 2)
				shouldHaveNext(c2, 1)
				shouldHaveNext(c3, 2)
				shouldHaveNext(c4, 3)
				shouldHaveNext(c5, 4)
				shouldHaveNext(c6, 1)
				shouldBeDone()
			}
		}

		test("don't visit source") {

			benzamidine.dfs(
				source = c1,
				visitSource = false
			).iterator().apply {

				shouldHaveNext(c, 1)
				shouldHaveNext(n1, 2)
				shouldHaveNext(n2, 2)
				shouldHaveNext(c2, 1)
				shouldHaveNext(c3, 2)
				shouldHaveNext(c4, 3)
				shouldHaveNext(c5, 4)
				shouldHaveNext(c6, 1)
				shouldBeDone()
			}
		}

		test("within 3 bonds") {

			benzamidine.dfs(
				source = n2,
				visitSource = false,
				shouldVisit = { _, _, dist -> dist <= 3 }
			).iterator().apply {

				shouldHaveNext(c, 1)
				shouldHaveNext(c1, 2)
				shouldHaveNext(c2, 3)
				shouldHaveNext(c6, 3)
				shouldHaveNext(n1, 2)
				shouldBeDone()
			}
		}

		test("only carbons") {

			benzamidine.dfs(
				source = c,
				visitSource = true,
				shouldVisit = { _, dst, _ -> dst.element == Element.Carbon }
			).iterator().apply {

				shouldHaveNext(c, 0)
				shouldHaveNext(c1, 1)
				shouldHaveNext(c2, 2)
				shouldHaveNext(c3, 3)
				shouldHaveNext(c4, 4)
				shouldHaveNext(c5, 5)
				shouldHaveNext(c6, 2)
				shouldBeDone()
			}
		}
	}

	group("bfs") {

		test("full order") {

			benzamidine.bfs(
				source = c1,
				visitSource = true
			).iterator().apply {

				shouldHaveNext(c1, 0)
				shouldHaveNext(c, 1)
				shouldHaveNext(c2, 1)
				shouldHaveNext(c6, 1)
				shouldHaveNext(n1, 2)
				shouldHaveNext(n2, 2)
				shouldHaveNext(c3, 2)
				shouldHaveNext(c5, 2)
				shouldHaveNext(c4, 3)
				shouldBeDone()
			}
		}

		test("don't visit source") {

			benzamidine.bfs(
				source = c1,
				visitSource = false
			).iterator().apply {

				shouldHaveNext(c, 1)
				shouldHaveNext(c2, 1)
				shouldHaveNext(c6, 1)
				shouldHaveNext(n1, 2)
				shouldHaveNext(n2, 2)
				shouldHaveNext(c3, 2)
				shouldHaveNext(c5, 2)
				shouldHaveNext(c4, 3)
				shouldBeDone()
			}
		}

		test("within 3 bonds") {

			benzamidine.bfs(
				source = n2,
				visitSource = false,
				shouldVisit = { _, _, dist -> dist <= 3 }
			).iterator().apply {

				shouldHaveNext(c, 1)
				shouldHaveNext(c1, 2)
				shouldHaveNext(n1, 2)
				shouldHaveNext(c2, 3)
				shouldHaveNext(c6, 3)
				shouldBeDone()
			}
		}

		test("only carbons") {

			benzamidine.bfs(
				source = c,
				visitSource = true,
				shouldVisit = { _, dst, _ -> dst.element == Element.Carbon }
			).iterator().apply {

				shouldHaveNext(c, 0)
				shouldHaveNext(c1, 1)
				shouldHaveNext(c2, 2)
				shouldHaveNext(c6, 2)
				shouldHaveNext(c3, 3)
				shouldHaveNext(c5, 3)
				shouldHaveNext(c4, 4)
				shouldBeDone()
			}
		}
	}
})
