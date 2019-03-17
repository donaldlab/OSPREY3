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

///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package edu.duke.cs.osprey.tests;
//
//import edu.duke.cs.osprey.astar.ewakstar.NewEWAKStarDoer;
//import edu.duke.cs.osprey.confspace.Sequence;
//import edu.duke.cs.osprey.confspace.SimpleConfSpace;
//import edu.duke.cs.osprey.confspace.Strand;
//import edu.duke.cs.osprey.ematrix.EnergyMatrix;
//import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
//import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
//import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
//import edu.duke.cs.osprey.energy.EnergyCalculator;
//import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
//import edu.duke.cs.osprey.parallelism.Parallelism;
//import edu.duke.cs.osprey.structure.Molecule;
//import edu.duke.cs.osprey.structure.PDBIO;
//
//import java.io.File;
//import java.util.ArrayList;
//import java.util.Arrays;
//
///**
// *
// * @author lowegard, based on NewCOMETS
// */
//public class TestEWAKStar {
//
//    public static void main(String[] args) {
//
//        Integer[] pos = new Integer[]{0, 1, 2, 3, 4, 5, 6, 7};
//        Integer[] posL = new Integer[]{4, 5, 6, 7};
//        Integer[] posP = new Integer[]{0, 1, 2, 3};
//
//        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
//                new String[]{"PHE"},
//                new String[]{"LYS"},
//                new String[]{"ILE"},
//                new String[]{"THR"},
//                new String[]{"PHE", "ALA", "VAL", "ILE", "LEU", "TYR"},
//                new String[]{"ASP"},
//                new String[]{"GLU"},
//                new String[]{"THR"}
//        );
//
//
//        int numCPUs = 4;
//        String PLmatrixName = "ewak.*";
//        String mutableType = "exact"; //can be "exact", "max", or "all"
//        int numMutable = 1;
//        int numFilteredSeqs = 10000;
//        double orderOfMag = 10.0;
//        double unboundEw = 30.0;
//        double boundEw = 30.0;
//        double ewakstarEw = 1.0;
//        double Ival = 0.0;
//        double epsilon = 0.01;
//        int maxPFConfs = 5000;
//        int numTopSeqs = 6;
//        boolean seqFilterOnly = false;
//        boolean wtBenchmark = false;
//        String startResL = "G648";
//        String endResL = "G654";
//        String startResP = "A155";
//        String endResP = "A194";
//        String pdbFile = "examples/python.KStar/2RL0.min.reduce.pdb";
//        String[] resNumsPL = new String[]{"A156", "A172", "A192", "A193", "G649", "G650", "G651", "G654"};
//        String[] resNumsL = new String[]{"G649", "G650", "G651", "G654"};
//        String[] resNumsP = new String[]{"A156", "A172", "A192", "A193"};
////
////        Integer[] pos = new Integer[]{0, 1, 2, 3, 4, 5, 6, 7, 8};
////        Integer[] posL = new Integer[]{3, 4, 5, 6, 7, 8};
////        Integer[] posP = new Integer[]{0, 1, 2};
//
////        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
////                new String[]{"HID", "HIE", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "ALA", "VAL", "ILE", "LEU", "GLY"},
////                new String[]{"ASN", "SER", "THR", "GLN", "HID", "VAL", "ILE", "LEU", "GLY", "ALA"},
////                new String[]{"TRP", "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "MET", "SER", "THR", "ASN", "GLN", "GLY"},
////                new String[]{"MET", "ILE", "ALA", "VAL", "LEU", "PHE", "TYR", "GLU", "ASP", "HID", "ASN", "GLN", "GLY"},
////                new String[]{"LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY"},
////                new String[]{"GLU", "ASP", "PHE", "TYR", "ALA", "VAL", "ILE", "LEU", "HIE", "HID", "ASN", "GLN", "GLY"},
////                new String[]{"LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY"},
////                new String[]{"LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY"},
////                new String[]{"GLU", "ASP", "PHE", "TYR", "ALA", "VAL", "ILE", "LEU", "HIE", "HID", "ASN", "GLN", "GLY"}
////        );
//
////        ArrayList<ArrayList<String>> AATypeOptions = toDoubleList(
////                new String[]{"HID", "HIE"},
////                new String[]{"ASN", "SER"},
////                new String[]{"TRP", "ALA"},
////                new String[]{"MET", "ILE"},
////                new String[]{"LEU", "ILE"},
////                new String[]{"GLU", "ASP"},
////                new String[]{"LEU", "ILE"},
////                new String[]{"LEU", "ILE"},
////                new String[]{"GLU", "ASP"}
////        );
////
////        String mutableType = "exact"; //can be "exact", "max", or "all"
////        int numMutable = 1;
////        int numFilteredSeqs = 10000;
////        double orderOfMag = 5.0;
////        double unboundEw = 8.0;
////        double boundEw = 8.0;
////        double ewakstarEw = 1.0;
////        double Ival = 0.0;
////        double epsilon = 0.68;
////        int maxPFConfs = 500;
////        int numTopSeqs = 5;
////        boolean seqFilterOnly = false;
////        boolean wtBenchmark = false;
////        String startResL = "0520";
////        String endResL = "0729";
////        String startResP = "039";
////        String endResP = "0339";
////        String pdbFile = "examples/python.EWAKStar/1A0R/1A0R.b.shell.pdb";
////        String[] resNumsPL = new String[]{"0311", "0313", "0332", "0601", "0605", "0696", "0697", "0698", "0729"};
////        String[] resNumsL = new String[]{"0601", "0605", "0696", "0697", "0698", "0729"};
////        String[] resNumsP = new String[]{"0311", "0313", "0332"};
//
//        Molecule mol = PDBIO.readFile(pdbFile);
//        Strand strandL = new Strand.Builder(mol).setResidues(startResL, endResL).build();
//        Strand strandP = new Strand.Builder(mol).setResidues(startResP, endResP).build();
//
//        for (Integer p : posL) {
//            strandL.flexibility.get(resNumsPL[p]).setLibraryRotamers(AATypeOptions.get(p)).addWildTypeRotamers().setContinuous();
//        }
//        for (Integer p : posP) {
//            strandP.flexibility.get(resNumsPL[p]).setLibraryRotamers(AATypeOptions.get(p)).addWildTypeRotamers().setContinuous();
//        }
//
//        Parallelism parallelism = Parallelism.makeCpu(numCPUs);
//
//        SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(strandP, strandL).build();
//        SimpleConfSpace confSpaceP = new SimpleConfSpace.Builder().addStrand(strandP).build();
//        SimpleConfSpace confSpaceL = new SimpleConfSpace.Builder().addStrand(strandL).build();
//
//        ForcefieldParams ffparams = new ForcefieldParams();
//        EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams).setParallelism(parallelism).build();
//        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(ecalc).setIsMinimizing(false).build();
//
//        SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpace, ecalc).build();
//        SimpleReferenceEnergies rigidEref = new SimpleReferenceEnergies.Builder(confSpace, rigidEcalc).build();
//
//        ConfEnergyCalculator confECalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
//                .setReferenceEnergies(eref)
//                .build();
//        ConfEnergyCalculator confRigidECalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
//                .setReferenceEnergies(rigidEref)
//                .build();
//
//        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc)
//                .setCacheFile(new File("ewak.PL.emat"))
//                .build()
//                .calcEnergyMatrix();
//
//        NewEWAKStarDoer ewakstar = new NewEWAKStarDoer(mutableType, numMutable,
//                numCPUs, wtBenchmark, seqFilterOnly, numTopSeqs, maxPFConfs, epsilon, confRigidECalc, confECalc, emat,
//                ecalc, confSpace, confSpaceL, confSpaceP, pos, posL, posP, numFilteredSeqs, orderOfMag, unboundEw,
//                boundEw, ewakstarEw, startResL, endResL, startResP, endResP, mol, resNumsPL, resNumsL, resNumsP, Ival,
//                PLmatrixName);
//
//
//        ArrayList<Sequence> bestSequences = ewakstar.run();
//
////        for (int i=0; i<bestSequences.size(); i++) {
////            System.out.println(bestSequences.get(i));
////        }
//
//
//    }
//
//
//    private static ArrayList<ArrayList<String>> toDoubleList(String[]... arr){
//        ArrayList<ArrayList<String>> ans = new ArrayList<>();
//        for(String[] a : arr){
//            ans.add( new ArrayList(Arrays.asList(a)) );
//        }
//        return ans;
//    }
//
//}
