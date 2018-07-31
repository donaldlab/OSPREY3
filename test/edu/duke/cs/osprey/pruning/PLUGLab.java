package edu.duke.cs.osprey.pruning;


import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.*;

import java.util.Arrays;

import static edu.duke.cs.osprey.tools.Log.log;


public class PLUGLab {

	public static void main(String[] args) {

		// load a protein
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A3").setLibraryRotamers(Strand.WildType);
		//strand.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		//strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType);
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// load probe
		Probe probe = new Probe();
		probe.matchTemplates(strand.templateLib);

		// TEMP
		for (Residue res : strand.mol.residues) {
			Probe.Template template = probe.getTemplate(res);
			log("%s -> %s", res.getPDBResNumber(), template.id);
		}

		// TEMP
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.build();

		Residue res1 = strand.mol.residues.getOrThrow("A3"); // glu
		Residue res2 = strand.mol.residues.getOrThrow("A5"); // lys
		AtomConnectivity.AtomPairs atomPairs = connectivity.getAtomPairs(res1, res2);
		for (AtomNeighbors.Type dist : Arrays.asList(AtomNeighbors.Type.NONBONDED, AtomNeighbors.Type.BONDED15H)) {
			for (int[] atomPair : atomPairs.getPairs(dist)) {
				Atom a1 = res1.atoms.get(atomPair[0]);
				Atom a2 = res2.atoms.get(atomPair[1]);
				Probe.Interaction inter = probe.getInteraction(a1, a2);
				log("%s:%-4s <-> %s:%-4s = %s", res1.getPDBResNumber(), a1.name, res2.getPDBResNumber(), a2.name, inter);
			}
		}
	}
}
