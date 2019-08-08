/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class TestUpdatingEnergyMatrix {
    public static final int NUM_FLEX = 10;
    public static final int NUM_TUPS = 100;

    @Test
    public void testTupleTrieRandom()
    {
        SimpleConfSpace confSpace = make1GUASmall(NUM_FLEX);
        UpdatingEnergyMatrix.TupleTrie trie = new UpdatingEnergyMatrix.TupleTrie(confSpace.positions);
        runManual(trie);
        for(int i = 0; i < NUM_TUPS; i++)
        {
            TupE tupE = makeRandomTupE(confSpace);
            System.out.println("Inserting "+tupE.tup.stringListing()+":"+tupE.E);
            trie.insert(tupE);
            assert(trie.contains(tupE.tup));
            List<TupE> corrections = trie.getCorrections(tupE.tup);
            for(TupE tupEout: corrections)
                System.out.println("Retrieved "+tupEout.tup.stringListing()+":"+tupEout.E);

        }
    }

    @Test
    public void testTupleTrieManual () {
        SimpleConfSpace confSpace = make1GUASmall(4);
        UpdatingEnergyMatrix.TupleTrie trie = new UpdatingEnergyMatrix.TupleTrie(confSpace.positions);
        runManual(trie);
    }

    private void runManual(UpdatingEnergyMatrix.TupleTrie trie) {
        for(TupE tupE : makeManualTupE()) {
            System.out.println("Inserting "+tupE.tup.stringListing()+":"+tupE.E);
            trie.insert(tupE);
            List<TupE> corrections = trie.getCorrections(tupE.tup);
            for(TupE tupEout: corrections)
                System.out.println("Retrieved "+tupEout.tup.stringListing()+":"+tupEout.E);
        }
    }

    private List<TupE> makeManualTupE() {
        List<TupE> out = new ArrayList<>();
        RCTuple tup4 = new RCTuple();
        tup4 = tup4.addRC(3, 20);
        double energy4 = Math.random() * -40;
        out.add(new TupE(tup4, energy4));
        RCTuple tup2 = new RCTuple();
        tup2 = tup2.addRC(2, 15);
        tup2 = tup2.addRC(3, 20);
        double energy2 = Math.random() * -40;
        out.add(new TupE(tup2, energy2));
        RCTuple tup3 = new RCTuple();
        tup3 = tup3.addRC(0, 15);
        tup3 = tup3.addRC(3, 20);
        double energy3 = Math.random() * -40;
        out.add(new TupE(tup3, energy3));
        RCTuple tup = new RCTuple();
        tup = tup.addRC(0, 15);
        tup = tup.addRC(2, 15);
        tup = tup.addRC(3, 20);
        double energy = Math.random() * -40;
        out.add(new TupE(tup, energy));
        return out;
    }

    private TupE makeRandomTupE(SimpleConfSpace space) {
        RCTuple tup = new RCTuple();
        for(SimpleConfSpace.Position pos: space.positions) {
            int index = pos.index;
            if(Math.random() >= 0.5 || (pos.index == 3 && tup.size() < 1))
            tup = tup.addRC(index, (int) (Math.random() * 20));
        }
        double energy = Math.random()*-40;
        return new TupE(tup, energy);
    }

    private static SimpleConfSpace make1GUASmall(int numFlex) {


        // configure the forcefield
        ForcefieldParams ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
                .addMoleculeForWildTypeRotamers(mol)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("1", "180")
                .build();
        int start = 21;
        for(int i = start; i < start+numFlex; i++) {
            protein.flexibility.get(i+"").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        }

        SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();

        return proteinConfSpace;
    }
}
