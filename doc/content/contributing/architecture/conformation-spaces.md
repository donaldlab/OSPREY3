+++
title = "Conformations and Conformation Spaces"
menuTitle = "Conformation Spaces"
weight = 1
+++


A *conformation space* is the set of all mutations we can make to a molecule,
and all the ways we allow the atom coordinates to change within a mutation.

A conformation space describes its mutations to a molecule by dividing up the atoms
into non-overlapping groups called *design positions*. A design position is a small
subset of atoms in the molecule that we can replace with other sets of atoms.

In a protein, a design position might include all the atoms at, say, residue 42.
That design position might also allow residue 42 to mutate from its wild-type Alanine
amino acid to a Valine amino acid. The design position might also include three
different sets of atom coordinates for the Valine amino acid, representing three
different conformational states for the residue.

Here's a description of a conformation space that includes the design
position we described above that contains the mutation from Alanine to Valine,
and also another design position at residue 17 that only has discrete flexibility
for the wild-type amino acid, Tyrosine:

* Residue 42
    * Alanine (wild-type)
        * conformation 1
    * Valine (mutation)
        * conformation 1
        * conformation 2
        * conformation 3
* Residue 17
    * Tyrosine (wild-type)
        * conformation 1
        * conformation 2
        * ...
        * conformation 8

{{% notice note %}}
The term *conformation* is slightly overloaded here. OSPREY uses several different
senses of the word conformation, depending on the context. A conformation can be:

1: An entire molecule where sets of atoms have been swapped out with others or moved from their original positions.

2: The mathematical representation of such a molecule in a conformation space.

3: A molecular fragment that can be exchanged with the existing atoms at a design position in a conformation space.
In the specific case of protein molecules, these are sometimes called *residue conformations*.
{{% /notice %}}

In this example conformation space, there are a total of 4 * 8 = 32 conformations,
since Residue 42 can be one of 4 different options, Reisdue 17 can be one of 8 different
options, and both of these choices are completely independent.

Also in this example conformation, there are a total of 2*1 = 1 sequences. The sequences
are (in no particular order):
* Alanine @ Residue 42, Tyrosine @ Residue 17
    * (this sequence has 1*8 = 8 conformations)
* Valine @ Residue 42, Tyrosine @ Residue 17
    * (this sequence has 3*8 = 24 conformations)

To describe these conformations compactly, we index everything in the conformation space
with integers. So the example conformation space becomes:

* `posi=0` Residue 42
    * Alanine (wild-type)
        * `confi=0` conformation 1
    * Valine (mutation)
        * `confi=1` conformation 1
        * `confi=2` conformation 2
        * `confi=3` conformation 3
* `posi=1` Residue 17
    * Tyrosine (wild-type)
        * `confi=0` conformation 1
        * `confi=1` conformation 2
        * ...
        * `confi=7` conformation 8

Now we can describe an entire conformation with a simple list of integers. For example,
the conformation `[0,0]` represents the molecule where Residue 42 is Alanine,
conformation 1 and Residue 17 is Tyrosine, conformation 1.
Another conformation is `[3,7]`.

Sometimes a conformation will only be partially assigned. In that case, the conformation
index at that design position will be `-1`. So the conformation `[-1,-1]` describes only
the atoms in the molecule that lie outside any design positions. The state of the atoms
inside the design positions is undefined. Or rather, unassigned.
