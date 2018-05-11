package edu.duke.cs.osprey;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.StrandFlex;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;



public class DfrG {

    public static void main(String[] args) {




        Molecule mol = PDBIO.read(FileTools.readFile("examples/dhfr/3frd-1k-substitute.GMEC.min.reduce.nuclear.pdb"));

        ForcefieldParams ffparams = new ForcefieldParams();

        ResidueTemplateLibrary custom_templatelib = new ResidueTemplateLibrary.Builder()
                .addTemplates(FileTools.readFile("examples/dhfr/datafiles/all_nuc94_and_gr.in"))
                .addTemplateCoords(FileTools.readFile("examples/dhfr/datafiles/all_amino_coords.in"))
                .addRotamers(FileTools.readFile("examples/dhfr/datafiles/GenericRotamers-ndp-fol.dat"))
                .build();

        Strand protein = new Strand.Builder(mol).setTemplateLibrary(custom_templatelib).setResidues("1", "158").build();

        protein.flexibility.get("31").setLibraryRotamers(protein.WildType,"LEU").setContinuous();
        protein.flexibility.get("54").setLibraryRotamers(protein.WildType).setContinuous();

        Strand ligand = new Strand.Builder(mol).setTemplateLibrary(custom_templatelib).setResidues("159", "159").build();

        ligand.flexibility.get("159").setLibraryRotamers(protein.WildType).setContinuous();

        StrandFlex ligandRotTrans = new StrandFlex.TranslateRotate();

        SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder().addStrand(protein).build();
        SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder().addStrand(ligand).build();
        SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder().addStrand(protein).addStrand(ligand,ligandRotTrans).build();

        Parallelism parallelism = Parallelism.makeCpu(4);


        new EnergyCalculator.Builder(complexConfSpace, ffparams)
                .setParallelism(parallelism)
                .use((ecalc) -> {

                    // how should we define energies of conformations?
                    KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
                        return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                                        .build()
                                        .calcReferenceEnergies()
                                ).build();
                    };


                    // how should confs be ordered and searched?
                    KStar.ConfSearchFactory confSearchFactory = (emat, pmat) -> {
                        return new ConfAStarTree.Builder(emat, pmat)
                                .setTraditional()
                                .build();
                    };


                    // run K*
                    TestKStar.Result result = new TestKStar.Result();

                    KStar.Settings settings = new KStar.Settings.Builder()
                            .setEpsilon(0.68)
                            .setStabilityThreshold(null)
                            .setMaxNumConfs(20000)
                            //.addScoreConsoleWriter(testFormatter)
                            .setConfDBPattern("kstar.*.conf.db")
                            .setShowPfuncProgress(true)
                            .build();
                    result.kstar = new KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalc, confEcalcFactory, confSearchFactory, settings);
                    result.scores = result.kstar.run();


                });


    }

}
