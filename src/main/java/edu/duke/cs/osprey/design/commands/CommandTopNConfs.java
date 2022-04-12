package edu.duke.cs.osprey.design.commands;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.design.Main;
import edu.duke.cs.osprey.design.models.MoleculeDesign;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.*;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

@Parameters(commandDescription = CommandTopNConfs.CommandDescription)
public class CommandTopNConfs extends RunDesignCommand {

    public static final String CommandName = "topn";
    static final String CommandDescription = "Calculate the energies of the top (n) confs";

    @SuppressWarnings("MismatchedQueryAndUpdateOfCollection") // Ignored because updated by command line arguments
    @Parameter(names = "--energy", description = "Analyze the energy of conformation(s).")
    private List<Integer> captureEnergies = new ArrayList<>();

    @Parameter(names = "--n", description = "Number of conformations to retrieve energy for.")
	private int numConfs = 10;

    @Parameter(names = "--m", description = "Number of residues for each conformation to get atom-specific breakdowns for.")
    private int topNRes = 5;

    @Override
    public int run(JCommander commander, String[] args) {
        var retVal = processHelpAndNoArgs(commander, args);

        if (retVal.isPresent()) {
            return retVal.get();
        }

        var designOpt = parseDesignSpec(MoleculeDesign.class);
        if (designOpt.isEmpty()) { // couldn't parse design
            return Main.Failure;
        }

        var design = designOpt.get();
        return runTopNConfs(design);
    }

    @Override
    public String getCommandName() {
        return CommandName;
    }

    @Override
    public String getCommandDescription() {
        return CommandDescription;
    }

    private int runTopNConfs(MoleculeDesign design) {
        /* This reads parm96a.dat, which contains the energy parameters of DNA, RNA, and protein residues */
        var ffParams = new ForcefieldParams();

        /* Maintains flexibility information with the molecule, and can use that to make new molecules */
        var confSpace = delegate.createConfSpace(design.molecule, ffParams);

        if (delegate.verifyInput) {
            return Main.Success;
        }

        EnergyMatrix emat;
        try (var ecalc = new EnergyCalculator.Builder(confSpace, ffParams).setParallelism(delegate.getParallelism()).build()) {
            emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                .build()
                .calcEnergyMatrix();

            var confECalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

            /* Contains the confSpace and a pruning matrix */
            var astarTree = new ConfAStarTree.Builder(emat, confSpace).setTraditional().build();

            for (int i = 0; i < numConfs; i++) {
                var sConf = astarTree.nextConf();
                var eConf = confECalc.calcEnergy(sConf);
                var interactions = confECalc.makeFragInters(new RCTuple(sConf.getAssignments()));
                var parametricMolecule = confSpace.makeMolecule(sConf.getAssignments());
                printConf(i, sConf, eConf);

                var efunc = new ResidueForcefieldEnergy(ecalc.resPairCache, interactions, parametricMolecule.mol);
                var energyContributions = efunc.getEnergyContributions(efunc.resPairs);
                var residueEnergies = energyContributions.stream()
                        .sorted(Comparator.comparingDouble(ResPairEnergyContribution::getEnergy).reversed())
                        .collect(Collectors.toList());

                printTopAtomContributions(residueEnergies);
            }

            return Main.Success;
        }
    }

    private void printTopAtomContributions(List<ResPairEnergyContribution> residueEnergies) {
    	var n = 3;
        for (int i = 0; i < n; i++) {
            for(var contribution : residueEnergies.get(i).getTopEnergyContributions(n)) {
                var builder = new StringBuilder();
            	if (contribution instanceof AtomPairElectrostaticContribution) {
            		builder.append("Electrostatic:\t");
                } else if (contribution instanceof AtomPairVanDerWaalsContribution) {
                    builder.append("van der Waals:\t");
                } else if (contribution instanceof AtomPairSolvationContribution) {
                    builder.append("Solvation:\t");
                }

            	builder.append(contribution.getEnergy());
            	builder.append("\t");
                builder.append(contribution.getAtom1().toString()).append("\t");
                builder.append(contribution.getAtom2().toString()).append("\t");
                System.out.println(builder.toString());
            }
        }
    }

    private void printConf(int confIdx, ConfSearch.ScoredConf sConf, ConfSearch.EnergiedConf eConf) {
        var score = sConf.getScore();
        var energy = eConf.getEnergy();

        System.out.println(String.format("Conf %d: Score:\t%f\t\tEnergy:\t%f", confIdx, score, energy));
    }
}
