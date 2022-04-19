+++
title = "Conformation Energy"
weight = 3
+++


## Task 2: Evaluate a conformation using physics

### Forcefields in more detail

Once we have a conformation, we want to score it to see how biophysically likely it
is to occur in some assumed environment, like the aqueous interior of a cell.

Osprey uses two different forcefields to do this:
 * [AMBER](https://en.wikipedia.org/wiki/AMBER)
 * [EEF1 or Effective Energy Function 1](http://charmm.sunhwanj.com/html/eef1.html)
 
For the AMBER forcefield, we only use the non-bonded terms modeling van der Waals
and electrostatic interactions. EEF1 models the interactions between the molecule
and an implicit bulk solvent (eg water).

The forcefields are essentially functions that map a set of atoms into a real
scalar value called *energy*. Generally speaking, the lower the energy of a conformation,
the more likey we are to observe that conformation occurring in nature.

These particular forcefields work by treating the conformation (which could be a single
molecule by itself, or multiple molecules interacting with each other) as a list of pairs
of atoms, where each pair of atoms is associated with *forcefield parameters*.

For example, in AMBER, every 1-4 bonded or non-bonded pair of atoms is associated with three values:
`vdwA`, `vdwB`, and `esQ`.

The energy of that atom pair is represented by an equation similar to:
```
E = eqQ/(r^2) + vdwA/(r^12) - vdwB/(r^6)
```
where `r` is the radius (aka distance) between the two atom positions.
The total energy for the conformation is the sum of all the energies of all
the atom pairs.

Therefore, computing the energy of a conformation consists of enumerating all
the atom pairs, computing the atom pair energies, and then computing the sum of all
the atom pair energies.

The conformation space also contains the lists of atom pairs for each interaction
between design positions, and the forcefield parameters for each atom pair.


### Energetic interactions between design positions

Since atoms in the molecule change depending on which conformation we're looking at,
the list of atom pairs and forcefield parameters must change too. Groups of atom pairs
in the conformation space are organized by pairs of design positions, and the conformations
those design positions can have.

For example, `posi1=0,posi2=1` represents all the atom pairs between the atoms in
position `0` and position `1` in a conformation. If the conformation in this case is `[0,0]`,
then we'll lookup the atom pairs for this position interaction in the conformation space
at position `0` conformation `0`, and position `1` conformation `0`.

Positions can also be paired with themselves. So `posi1=0,posi2=0` refers to all the atom
pairs *within* position `0`.

Positions can also be paired with the *static* position, which is all of the atoms outside
of any design position. eg `posi1=0,posi2=static`.

Finally, there's the interaction within the static position,
eg `posi1=static,posi2=static`. This represents the internal energy of all the atoms
not in any design positions.

The `static` position was so named because its atoms aren't changed by any design position.
However, other mechanisms can move the `static` atoms, such as dihedral angles defined
for a whole molecule, or translation and rotation motions for a whole molecule.
So it's not correct to assume that the atom coordinates of the `static` position are
constant.

When we're evaluating a conformation, we'll want the complete set of all position
interactions, including interactions within design positions, interactions between all
pairs of design positions, interactions between all positions and the static position,
and interactions within the static position.

Osprey's design algorithm that try to find conformations to evaluate sometimes only
care about interactions between some of the positions. For example, a design algorithm
in Osprey may want to compute the pairwise energy between positions `0` and `1`
independently of the other interactions. So our energy calculators need to support
compute energies on arbitrary subsets of the position interactions.


### Computing minimized energies for conformations 

Computing the energy from the atom coordinates exactly as described by the conformation
space is a crude model of molecular conformational flexibility. Often, small overlaps
between atoms that have just been mutated would lead to extremely high forces between
these atoms (recall the `r^12` terms in the van der Waals equations!) and hence extremely
poor energies. Osprey improves its biophysical modeling by allowing the atom coordinates
to "relax" before recording the conformation energy. That is, the atoms are allowed to
move a little closer to an equilibrium state where the atomic forces are more in balance.
As long as the atoms don't move too far away from their original positions, it's reasonable
to say the "relaxed" conformation is a better physical model of the original conformation
before relaxation.

Osprey performs this relaxation by minimizing the energy function over continuous motions
of the atoms, like dihedral rotations around rotatable bonds. The parameters for these
dihedral rotations are also defined in the conformation space. In general, there are
many different kinds of continuous motions we need to model for molecular flexiblity,
and dihedral angles are one of these. More continuous motions are planned for the future,
but for now, dihedral angle motions are the only motions implemented in the newest
minimizer code because they're the ones we use most often by far.

Our minimization algorithm of choice in Osprey is
[Cyclic Coordinate Descent (CCD)](https://en.wikipedia.org/wiki/Coordinate_descent) due
to its simplicity and speed. All the continuous motions are translated into scalar
[*degrees of freedom*](https://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)).
For dihedral angle rotations, there's only one degree of freedom:
the rotation angle. But in general, a continuous motion can map to many degrees of freedom.
For example, translation and rotation motions map to 6 degees of freedom.

Each degree of freedom is bounded to lie within an interval. These bounding intervals
are also described in the conformation space.

Since the forcefield can be broken down by position interaction, the line search
subroutine inside CCD can focus only on the atom pairs affected by the degree of freedom
being searched. This single optimization causes Osprey's CCD implementation to consistently
out-perform other gradient-based minimization algorithms we've tried in the past, like
[Newtonian methods](https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization).
