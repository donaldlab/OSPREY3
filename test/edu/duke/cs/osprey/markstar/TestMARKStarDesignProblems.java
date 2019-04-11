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

package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.confspace.CATSStrandFlex;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.markstar.TestMARKStar.Result;
import edu.duke.cs.osprey.kstar.TestKStar.ConfSpaces;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.util.List;

import static edu.duke.cs.osprey.kstar.TestBBKStar.runBBKStar;
import static edu.duke.cs.osprey.markstar.TestMARKStar.*;
import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

public class TestMARKStarDesignProblems {
    @Test
    public void test5UCE () {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/5uce.cfs");
            runMARKStar(confSpaces, 0.01);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test5UCF () {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/5ucf.cfs");
            runMARKStar(confSpaces, 0.01);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test5UCFMut () {
        try {
            ConfSpaces confSpaces = loadFromCFS("examples/python.KStar/5ucf.cfs");
            runMARKStar(confSpaces, 0.01);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test1A0RBBKStar() {
        ConfSpaces confSpaces = make1A0RBBKStar();
        int numSequences = 2;
        double epsilon = 0.99;
        String kstartime = "(Not run)";
        Stopwatch watch = new Stopwatch();
        boolean runkstar = true;
        if(runkstar) {
            watch.start();
            runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
            watch.stop();
            kstartime = watch.getTime(2);
            watch.reset();
        }
        watch.start();
        runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
        watch.stop();
        String bbkstartime = watch.getTime(2);
        System.out.println("MARK*: "+bbkstartime+", Traditional K*: "+kstartime);
    }

    private static ConfSpaces make1A0RBBKStar() {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.readFile("examples/python.KStar/1a0r_prepped.pdb");

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("P13", "P230")
                .build();
        protein.flexibility.get("P219").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("P220").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("P222").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("P223").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR",
                "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIP", "HIE", "HID", "ASP", "GLU", "ASN", "GLN",
                "GLY").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("P224").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("B2", "B340")
                .build();
        ligand.flexibility.get("B46").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("B47").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // make the conf spaces ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }

    @Test
    public void test2RL0Python() {
        ConfSpaces confSpaces = make2RL0Python();
        List<MARKStar.ScoredSequence> markstarResult = runMARKStar(confSpaces, 0.1).scores;
        boolean runkstar = true;
        List<KStar.ScoredSequence> kstarResult = null;
        if(runkstar) {
            kstarResult = runKStar(confSpaces, 0.1);
        }
        for(MARKStar.ScoredSequence seq: markstarResult) {
            System.out.println(seq);
        }
        if(runkstar)
            for (KStar.ScoredSequence seq : kstarResult) {
                System.out.println(seq);
            }
    }

    public static ConfSpaces make2RL0Python() {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.readFile("examples/python.KStar/2RL0.min.reduce.pdb");

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addMoleculeForWildTypeRotamers(mol)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // make the conf spaces ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }

    public static ConfSpaces make2RL0() {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.readFile("examples/python.KStar/2RL0.min.reduce.pdb");

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addMoleculeForWildTypeRotamers(mol)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("G648", "G654")
                .build();
        protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A155", "A194")
                .build();
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

        // make the conf spaces ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }

    @Test
    public void test2RL0() {

        double epsilon = 0.95;
        TestMARKStar.Result result = runMARKStar(make2RL0(), epsilon);

        // check the results (values collected with e = 0.1 and 64 digits precision)
        // NOTE: these values don't match the ones in the TestKSImplLinear test because the conf spaces are slightly different
        // also, the new K* code has been updated to be more precise
        assertSequence(result,   0, "PHE ASP GLU THR PHE LYS ILE THR", 4.300422e+04, 4.347270e+30, 4.201039e+50, epsilon); // K* = 15.351629 in [15.312533,15.396058] (log10)
        assertSequence(result,   1, "PHE ASP GLU THR PHE LYS ILE SER", 4.300422e+04, 1.076556e+30, 4.045744e+50, epsilon); // K* = 15.941451 in [15.878562,15.986656] (log10)
        assertSequence(result,   2, "PHE ASP GLU THR PHE LYS ILE ASN", 4.300422e+04, 4.650623e+29, 1.854792e+49, epsilon); // K* = 14.967273 in [14.920727,15.011237] (log10)
        assertSequence(result,   3, "PHE ASP GLU THR PHE LYS ALA THR", 4.300422e+04, 1.545055e+27, 7.003938e+45, epsilon); // K* = 14.022887 in [13.984844,14.066296] (log10)
        assertSequence(result,   4, "PHE ASP GLU THR PHE LYS VAL THR", 4.300422e+04, 5.694044e+28, 9.854022e+47, epsilon); // K* = 14.604682 in [14.558265,14.649295] (log10)
        //assertSequence(result,   5, "PHE ASP GLU THR PHE LYS LEU THR", 4.300422e+04, 3.683508e-11, 3.644143e+08, epsilon); // K* = 14.361823 in [14.324171,14.405292] (log10)
        //assertSequence(result,   6, "PHE ASP GLU THR PHE LYS PHE THR", 4.300422e+04, 2.820863e+24, null        , epsilon); // K* = none      in [-Infinity,-Infinity] (log10)
        //assertSequence(result,   7, "PHE ASP GLU THR PHE LYS TYR THR", 4.300422e+04, 1.418587e+26, null        , epsilon); // K* = none      in [-Infinity,-Infinity] (log10)
        assertSequence(result,   8, "PHE ASP GLU THR PHE ASP ILE THR", 4.300422e+04, 4.294128e+20, 1.252820e+36, epsilon); // K* = 10.831503 in [10.802450,10.873479] (log10)
        assertSequence(result,   9, "PHE ASP GLU THR PHE GLU ILE THR", 4.300422e+04, 4.831904e+20, 2.273475e+35, epsilon); // K* = 10.039061 in [10.009341,10.081194] (log10)
        assertSequence(result,  10, "PHE ASP GLU THR TYR LYS ILE THR", 4.300422e+04, 4.583429e+30, 2.294671e+50, epsilon); // K* = 15.066019 in [15.026963,15.110683] (log10)
        assertSequence(result,  11, "PHE ASP GLU THR ALA LYS ILE THR", 4.300422e+04, 3.310340e+28, 2.171286e+47, epsilon); // K* = 14.183333 in [14.156555,14.226141] (log10)
        assertSequence(result,  12, "PHE ASP GLU THR VAL LYS ILE THR", 4.300422e+04, 9.004068e+29, 1.866542e+49, epsilon); // K* = 14.683088 in [14.652599,14.725905] (log10)
        assertSequence(result,  13, "PHE ASP GLU THR ILE LYS ILE THR", 4.300422e+04, 3.398648e+30, 1.598348e+50, epsilon); // K* = 15.038854 in [15.002827,15.082763] (log10)
        assertSequence(result,  14, "PHE ASP GLU THR LEU LYS ILE THR", 4.300422e+04, 5.296285e+27, 1.045234e+47, epsilon); // K* = 14.661731 in [14.616778,14.705997] (log10)
        assertSequence(result,  15, "PHE ASP GLU SER PHE LYS ILE THR", 3.153051e+06, 4.347270e+30, 1.477525e+53, epsilon); // K* = 16.032587 in [15.970427,16.078237] (log10)
        assertSequence(result,  16, "PHE ASP GLU ASN PHE LYS ILE THR", 1.484782e+06, 4.347270e+30, 9.591789e+52, epsilon); // K* = 16.172020 in [16.114119,16.217413] (log10)
        assertSequence(result,  17, "PHE ASP GLU GLN PHE LYS ILE THR", 2.531411e+06, 4.347270e+30, 3.346589e+53, epsilon); // K* = 16.483023 in [16.424659,16.528464] (log10)
        assertSequence(result,  18, "PHE ASP ASP THR PHE LYS ILE THR", 1.216972e+01, 4.347270e+30, 1.469949e+45, epsilon); // K* = 13.443805 in [13.413093,13.487750] (log10)
        assertSequence(result,  19, "PHE GLU GLU THR PHE LYS ILE THR", 1.986991e+05, 4.347270e+30, 1.097189e+50, epsilon); // K* = 14.103869 in [14.056796,14.148018] (log10)
        assertSequence(result,  20, "TYR ASP GLU THR PHE LYS ILE THR", 1.666243e+04, 4.347270e+30, 2.814673e+46, epsilon); // K* = 11.589473 in [11.550098,11.633104] (log10)
        assertSequence(result,  21, "ALA ASP GLU THR PHE LYS ILE THR", 6.100779e+02, 4.347270e+30, 1.671418e+45, epsilon); // K* = 11.799483 in [11.777675,11.841368] (log10)
        assertSequence(result,  22, "VAL ASP GLU THR PHE LYS ILE THR", 1.271497e+02, 4.347270e+30, 2.380877e+45, epsilon); // K* = 12.634205 in [12.613207,12.676919] (log10)
        assertSequence(result,  23, "ILE ASP GLU THR PHE LYS ILE THR", 5.942890e+02, 4.347270e+30, 2.012605e+46, epsilon); // K* = 12.891544 in [12.844167,12.936569] (log10)
        assertSequence(result,  24, "LEU ASP GLU THR PHE LYS ILE THR", 4.614233e+00, 4.347270e+30, 4.735376e+43, epsilon); // K* = 12.373038 in [12.339795,12.417250] (log10)
    }

    public static void assertSequence(TestMARKStar.Result result, int sequenceIndex, String sequence, Double proteinQStar, Double ligandQStar, Double complexQStar, double epsilon) {

        MARKStar.ScoredSequence scoredSequence = result.scores.get(sequenceIndex);

        // check the sequence
        assertThat(scoredSequence.sequence.toString(Sequence.Renderer.ResType), is(sequence));

        // check q* values and epsilon
        assertResult(scoredSequence.score.protein, proteinQStar, epsilon);
        assertResult(scoredSequence.score.ligand, ligandQStar, epsilon);
        assertResult(scoredSequence.score.complex, complexQStar, epsilon);
    }

    @Test
    public void test2RL0MARKStar() {
        ConfSpaces confSpaces = null;
        try {
            confSpaces = loadSSFromCFS("examples/python.KStar/2rl0_A_13res_ss_6.837E+28.cfs");
            int numSequences = 2;
            double epsilon = 0.99999;
            String kstartime = "(Not run)";
            Stopwatch watch = new Stopwatch();
            boolean runkstar = false;
            if(runkstar) {
                watch.start();
                runKStar(confSpaces, epsilon);
                watch.stop();
                kstartime = watch.getTime(2);
                watch.reset();
            }
            watch.start();
            runMARKStar(confSpaces, epsilon);
            watch.stop();
            String bbkstartime = watch.getTime(2);
            System.out.println("MARK*: "+bbkstartime+", Traditional K*: "+kstartime);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        int numSequences = 2;
        double epsilon = 0.99;
        String kstartime = "(Not run)";
        Stopwatch watch = new Stopwatch();
        boolean runkstar = false;
        TestComparison(confSpaces, epsilon, runkstar);
    }

    public void test1GUA11() {

        double epsilon = 0.999;
        Stopwatch runtime = new Stopwatch().start();
        Result result = runMARKStar(make1GUA11(), epsilon);
        runtime.stop();

        for (int index = 0; index <6; index++){
            printSequence(result, index);
        }
        // check the results (values collected with e = 0.1 and 64 digits precision)
        assertSequence(result,   0, "HIE VAL", 1.194026e+42, 2.932628e+07, 1.121625e+66, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [7.467257 , 7.467257] (log10)                    complex [66.049848,66.051195] (log10)                    K* = 16.505577 in [16.505576,16.506925] (log10)
        assertSequence(result,   1, "HIE HID", 1.194026e+42, 5.738568e+07, 3.346334e+66, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [7.758803 , 7.758803] (log10)                    complex [66.524569,66.543073] (log10)                    K* = 16.688752 in [16.688752,16.707256] (log10)
        assertSequence(result,   2, "HIE HIE", 1.194026e+42, 6.339230e+06, 5.544100e+65, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [6.802036 , 6.802036] (log10)                    complex [65.743831,65.769366] (log10)                    K* = 16.864781 in [16.864780,16.890316] (log10)
        assertSequence(result,   3, "HIE LYS", 1.194026e+42, 6.624443e+04, 3.315130e+63, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [4.821149 , 4.826752] (log10)                    complex [63.520501,63.563549] (log10)                    K* = 16.622337 in [16.616735,16.665386] (log10)
        assertSequence(result,   4, "HIE ARG", 1.194026e+42, 1.196619e+05, 5.375633e+64, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [5.077956 , 5.087238] (log10)                    complex [64.730430,64.774106] (log10)                    K* = 17.575460 in [17.566178,17.619136] (log10)
        assertSequence(result,   5, "HID VAL", 9.813429e+41, 2.932628e+07, 2.680104e+66, epsilon); // protein [41.991821,41.992159] (log10)                    ligand [7.467257 , 7.467257] (log10)                    complex [66.428152,66.446408] (log10)                    K* = 16.969074 in [16.968735,16.987330] (log10)
        System.out.println("Total time: "+runtime.getTime(2));
    }

    @Test
    public void test4KT6 () {
        ConfSpaces confSpaces = make4KT6();
        final double epsilon = 0.99;
        boolean runkstar = false;
        TestComparison(confSpaces, epsilon, runkstar);
    }

    private ConfSpaces make4KT6() {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.readFile("examples/python.KStar/4kt6_prepped.pdb");

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addMoleculeForWildTypeRotamers(mol)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("C193", "C446")
                .build();
        protein.flexibility.get("C290").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("D1", "D161")
                .build();
        ligand.flexibility.get("D148").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("D149").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("D151").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("D152").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("D153").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("D155").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("D156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // make the conf spaces ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }

    public void test2RL0BBKStar() {
        ConfSpaces confSpaces = null;
        try {
            confSpaces = loadFromCFS("examples/python.KStar/2rl0_A_13res_6.837E+28.cfs");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        int numSequences = 2;
        double epsilon = 0.99999;
        String kstartime = "(Not run)";
        Stopwatch watch = new Stopwatch();
        boolean runkstar = false;
        if(runkstar) {
            watch.start();
            runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
            watch.stop();
            kstartime = watch.getTime(2);
            watch.reset();
        }
        watch.start();
        runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
        watch.stop();
        String bbkstartime = watch.getTime(2);
        System.out.println("MARK*: "+bbkstartime+", Traditional K*: "+kstartime);
    }

    @Test
    public void test2XXM() {
        ConfSpaces confSpaces = make2XXM();
        final double epsilon = 0.95;
        String kstartime = "(not run)";
        boolean runkstar = true;
        Stopwatch runtime = new Stopwatch().start();
        List<KStar.ScoredSequence> kStarSeqs = null;
        if(runkstar) {
            kStarSeqs = runKStar(confSpaces, epsilon);
            runtime.stop();
            kstartime = runtime.getTime(2);
            runtime.reset();
            runtime.start();
        }
        Result result = runMARKStar(confSpaces, epsilon);
        runtime.stop();
        String markstartime = runtime.getTime(2);
        System.out.println("MARK* time: "+markstartime+", K* time: "+kstartime);
        for(MARKStar.ScoredSequence seq: result.scores)
            printMARKStarComputationStats(seq);
        if(runkstar)
            for(KStar.ScoredSequence seq: kStarSeqs)
                printKStarComputationStats(seq);
    }

    @Test
    public void test3U7Y() {
        try {
            //runBBKStar(loadSSFromCFS("3u7y_L_15res_1.326E+48.cfs"),
             //       5, 0.99, null, 2, true);
            runMARKStar(loadSSFromCFS("3u7y_L_15res_1.326E+48.cfs"),0.99);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test1a0r() {

        try {
            runBBKStar(loadFromCFS("examples/python.KStar/1a0r_B_10res_8.034E+41.cfs"),
                   5, 0.68, null, 2, true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private ConfSpaces make2XXM() {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.readFile("examples/python.KStar/2xxm_prepped.pdb");

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("A146", "A218")
                .build();
        protein.flexibility.get("A177").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A179").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A180").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("A181").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("B4", "B113")
                .build();
        ligand.flexibility.get("B58").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("B60").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("B61").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("B64").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // make the conf spaces ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }

    public static ConfSpaces make1A0R() {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1A0R_prepped.pdb"));

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
                .addMoleculeForWildTypeRotamers(mol)
                .build();

        // define the protein strand
        Strand protein = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("039", "0339")
                .build();
        protein.flexibility.get("0313").setLibraryRotamers(Strand.WildType, "ASN", "SER", "THR", "GLN", "HID", "ALA", "VAL", "ILE", "LEU", "GLY", "ALA", "VAL", "GLY").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("0311").setLibraryRotamers(Strand.WildType, "HID", "HIE", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "ALA", "VAL", "ILE", "LEU", "GLY").addWildTypeRotamers().setContinuous();
        protein.flexibility.get("0332").setLibraryRotamers(Strand.WildType, "TRP", "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "MET", "SER", "THR", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();

        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("0520", "0729")
                .build();
        ligand.flexibility.get("0605").setLibraryRotamers(Strand.WildType, "LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("0696").setLibraryRotamers(Strand.WildType, "GLU", "ASP", "PHE", "TYR", "ALA", "VAL", "ILE", "LEU", "HIE", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("0697").setLibraryRotamers(Strand.WildType, "LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("0698").setLibraryRotamers(Strand.WildType, "LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("0601").setLibraryRotamers(Strand.WildType, "MET", "ILE", "ALA", "VAL", "LEU", "PHE", "TYR", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("0729").setLibraryRotamers(Strand.WildType, "GLU", "ASP", "PHE", "TYR", "ALA", "VAL", "ILE", "LEU", "HIE", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();

        // make the complex conf space ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrands(protein, ligand)
                .build();

        return confSpaces;
    }

    public static ConfSpaces make1GUASmallCATS(int numFlex) {

        ConfSpaces confSpaces = new ConfSpaces();

        // configure the forcefield
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        // make sure all strands share the same template library
        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
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
        CATSStrandFlex bbflex = new CATSStrandFlex(protein, "22", "25");



        // define the ligand strand
        Strand ligand = new Strand.Builder(mol)
                .setTemplateLibrary(templateLib)
                .setResidues("181", "215")
                .build();
        ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();

        // make the complex conf space ("complex" SimpleConfSpace, har har!)
        confSpaces.protein = new SimpleConfSpace.Builder()
                .addStrand(protein, bbflex)
                .build();
        confSpaces.ligand = new SimpleConfSpace.Builder()
                .addStrand(ligand)
                .build();
        confSpaces.complex = new SimpleConfSpace.Builder()
                .addStrand(protein, bbflex)
                .addStrand(ligand)
                .build();

        return confSpaces;
    }

    @Test
    public void testk4hemBig () {
        ConfSpaces confSpaces = makeConfSpaces(
                "examples/python.KStar/4hem_prepped.pdb",
                new String[]{"A34", "A163"},
                new String[]{"A67", "A68", "A78", "A76", "A88", "A89", "A157"},
                new String[]{"G1", "G123"},
                new String[]{"G101", "G102", "G103", "G57", "G104"});
        boolean runkstar = false;
        final double epsilon = 0.9999;
        TestComparison(confSpaces, epsilon, runkstar);
    }

    @Test
    public void test1a0rSlowness() {
        ConfSpaces confSpaces = makeConfSpaces(
                "examples/python.KStar/1a0r_prepped.pdb",
                new String[]{"B2", "B340"},
                new String[]{"B246", "B188", "B275", "B274", "B231", "B230", "B290"},
                new String[]{"P13", "P230"},
                new String[]{"P20", "P23", "P17"}
        );
        boolean runkstar = true;
        final double epsilon = 0.99999;
        TestComparison(confSpaces, epsilon, runkstar);

    }

    @Test
    public void test3K3Q () {

        ConfSpaces confSpaces = makeConfSpaces(
                "examples/python.KStar/3k3q_prepped.pdb",
                new String[]{"C253", "C417"},
                new String[]{"C412", "C404", "C409", "C391", "C390", "C384", "C392", "C383", "C397", "C378", "C376", "C413"},
                new String[]{"B2","B250"},
                new String[]{"B193", "B216", "B214"}
        );
        final double epsilon = 0.99999;
        boolean runkstar = false;
        TestComparison(confSpaces, epsilon, runkstar);
    }

    @Test
    public void test4WWI() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/4wwi_F_12res_1.780E+30.cfs");
            TestComparison(confSpaces, 0.99, true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test1A0R19ResE52() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/1a0r_B_19res_9.186E+52.cfs");
            TestComparison(confSpaces, 0.99, false);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test4WWI27ResE78() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/4wwi_B_27res_1.845E+78.cfs");
            TestComparison(confSpaces, 0.999999999999, false);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test1B6C() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/1b6c_E_19res_1.010E+52.cfs");
            TestComparison(confSpaces, 0.99, true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test4HEM12ResE66() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/4hem_A_12res_4.754E+66.cfs");
            TestComparison(confSpaces, 0.99, true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test3MA28ResE24() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/3ma2_A_8res_1.509E+24.cfs");
            TestComparison(confSpaces, 0.01, true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test3CAL3ResE20() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/3cal_C_3res_7.245E+20.cfs");
            TestComparison(confSpaces, 0.99, false);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test3K3QHuge() {
        try {
            ConfSpaces confSpaces = loadSSFromCFS("examples/python.KStar/3k3q_C_26res_1.223E+76.cfs");
            TestComparison(confSpaces, 0.99, true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
