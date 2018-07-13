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

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.SVG;
import edu.duke.cs.osprey.tools.SVGPlot;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class KStarRunner {

	public static void main(String[] args)
	throws Exception {

		// try JJ's design
		Molecule mol = PDBIO.readFile("/home/jeff/dlab/osprey test cases/jj-serialization/gp120SRVRC26.09SR.pdb");

		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder()
			.addMoleculeForWildTypeRotamers(mol)
			.addTemplates(FileTools.readFile("/home/jeff/dlab/osprey test cases/jj-serialization/all_nuc94_and_gr.in"))
			.addTemplateCoords(FileTools.readFile("/home/jeff/dlab/osprey test cases/jj-serialization/all_amino_coords.in"))
			.addRotamers(FileTools.readFile("/home/jeff/dlab/osprey test cases/jj-serialization/GenericRotamers.dat"))
			.build();

		Strand ligand = new Strand.Builder(mol)
			.setResidues("H1792", "L2250")
			.setTemplateLibrary(templateLib)
			.build();
		ligand.flexibility.get("H1901").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1904").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1905").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIS", "ASP", "GLU", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1906").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1907").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("H1908").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		Strand target = new Strand.Builder(mol)
			.setResidues("F379", "J1791")
			.setTemplateLibrary(templateLib)
			.build();
		target.flexibility.get("G973").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G977").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G978").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G979").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("G980").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		target.flexibility.get("J1448").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.addStrand(target)
			.build();

		SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();

		SimpleConfSpace targetConfSpace = new SimpleConfSpace.Builder()
			.addStrand(target)
			.build();

		//runKStar(targetConfSpace, ligandConfSpace, complexConfSpace);
		analyze(targetConfSpace, ligandConfSpace, complexConfSpace);
	}

	private static void runKStar(SimpleConfSpace targetConfSpace, SimpleConfSpace ligandConfSpace, SimpleConfSpace complexConfSpace) {

		KStar.Settings settings = new KStar.Settings.Builder()
			.setEpsilon(0.9999)
			.setStabilityThreshold(null)
			.setMaxSimultaneousMutations(3)
			.addScoreConsoleWriter()
			.addScoreFileWriter(new File("kstar.txt"))
			.setShowPfuncProgress(true)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complexConfSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			//.setParallelism(Parallelism.make(4, 1, 1))
			.build()
		) {
			KStar kstar = new KStar(targetConfSpace, ligandConfSpace, complexConfSpace, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, ecalc).build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.setCacheFile(new File(String.format("kstar.emat.%s.dat", info.id)))
					.build()
					.calcEnergyMatrix();

				info.confSearchFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();
			}
			kstar.run();
		}
	}

	private static void analyze(SimpleConfSpace targetConfSpace, SimpleConfSpace ligandConfSpace, SimpleConfSpace complexConfSpace)
	throws Exception {

		class Entry {
			public Sequence sequence;
			public double ligandPfuncLB;
			public double ligandPfuncUB;
			public double complexPfuncLB;
			public double complexPfuncUB;
			public double lowerBinding;
			public double upperBinding;
		}
		List<Entry> entries = new ArrayList<>();

		// read the kstar data into entries
		List<String> lines = Files.readAllLines(Paths.get("kstar.better.txt"));
		for (String line : lines) {

			// skip lines that don't start with "sequence"
			if (!line.startsWith("sequence ")) {
				continue;
			}

			// e.g.:
			// 0                                                                                                   1                                                                                                   2
			// 0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8
			// 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
			// sequence    1/ 216   A5=lys A6=hie A7=tyr A11=val A12=val A13=met   K*(log10): 23.997236 in [23.825229,24.188758] (log10)   protein: [36.043239,36.043240] (log10)                      ligand: [34.175654,34.347660] (log10)                      complex: [94.216129,94.407650] (log10)
			// sequence   22/ 361   H1901=asp H1904=VAL H1905=LEU H1906=gly H1907=lys H1908=gln G973=thr G977=arg G978=asp G979=lys G980=lys J1448=arg   K*(log10): none      in [-4.173163,32.036108] (log10)

			//int sequenceEnd = 65;
			int sequenceEnd = 135;

			Entry entry = new Entry();

			// re-construct the sequence
			entry.sequence = complexConfSpace.makeUnassignedSequence();
			for (String assignment : line.substring(21, sequenceEnd).split(" ")) {
				String parts2[] = assignment.split("=");
				SimpleConfSpace.Position pos = complexConfSpace.getPositionOrThrow(parts2[0]);
				entry.sequence.set(pos.resNum, parts2[1]);
			}

			// read the pfuncs
			entry.ligandPfuncLB = Math.pow(10, Double.parseDouble(line.substring(193 - 65 + sequenceEnd, 202 - 65 + sequenceEnd)));
			entry.ligandPfuncUB = Math.pow(10, Double.parseDouble(line.substring(203 - 65 + sequenceEnd, 212 - 65 + sequenceEnd)));
			entry.complexPfuncLB = Math.pow(10, Double.parseDouble(line.substring(253 - 65 + sequenceEnd, 262 - 65 + sequenceEnd)));
			entry.complexPfuncUB = Math.pow(10, Double.parseDouble(line.substring(263 - 65 + sequenceEnd, 272 - 65 + sequenceEnd)));

			// compute the binding score
			entry.lowerBinding = entry.complexPfuncLB/entry.ligandPfuncUB;
			entry.upperBinding = entry.complexPfuncUB/entry.ligandPfuncLB;

			entries.add(entry);
		}
		log("%d sequences parsed", entries.size());

		entries.sort(Comparator.comparing((Entry entry) -> entry.upperBinding).reversed());

		// plot the pfunc bounds
		{
			SVG svg = new SVG();

			SVG.StyleClass wildtypeBoxStyle = svg.makeStyleClass("wildtype-box-style");
			wildtypeBoxStyle.priority = 10;
			wildtypeBoxStyle.setStrokeColor(0x66cc55);
			wildtypeBoxStyle.setStrokeWidth(0.5);

			SVGPlot.Boxes boxes = new SVGPlot.Boxes();
			boxes.boxStyle.setNoFill();
			for (Entry entry : entries) {

				SVGPlot.Boxes.Box box = boxes.addBox(
					MathTools.log10p1(entry.ligandPfuncLB),
					MathTools.log10p1(entry.ligandPfuncUB),
					MathTools.log10p1(entry.complexPfuncLB),
					MathTools.log10p1(entry.complexPfuncUB)
				);
				box.id = String.format("%s: ligand:[%.2f,%.2f] complex:[%.2f,%.2f]",
					entry.sequence.toString(),
					MathTools.log10p1(entry.ligandPfuncLB),
					MathTools.log10p1(entry.ligandPfuncUB),
					MathTools.log10p1(entry.complexPfuncLB),
					MathTools.log10p1(entry.complexPfuncUB)
				);
				if (entry.sequence.isWildType()) {
					box.extraStyle = wildtypeBoxStyle;
				}
			}
			boxes.xaxis = boxes.makeXAxis();
			boxes.xaxis.tickFormat = "%.0f";
			boxes.yaxis = boxes.makeYAxis();
			boxes.yaxis.tickFormat = "%.0f";
			boxes.draw(svg);
			boxes.setBounds(svg, 10, 16);

			svg.finish().write(new File("pfuncs.kstar.svg"));
		}

		// plot the binding scores
		{
			SVG svg = new SVG();

			SVG.StyleClass wildtypeIntervalStyle = svg.makeStyleClass("wildtype-interval-style");
			wildtypeIntervalStyle.priority = 10;
			wildtypeIntervalStyle.setStrokeColor(0x66cc55);
			wildtypeIntervalStyle.setStrokeWidth(0.5);

			SVGPlot.Intervals intervals = new SVGPlot.Intervals();
			intervals.intervalWidth = 2.0;

			for (Entry entry : entries) {

				// skip entries that have no binding score
				if (Double.isNaN(entry.lowerBinding) || Double.isNaN(entry.upperBinding)) {
					log("sequence %s has no binding score, ligand:[%.1f,%.1f] complex:[%.1f,%.1f]",
						entry.sequence, entry.ligandPfuncLB, entry.ligandPfuncUB, entry.complexPfuncLB, entry.complexPfuncUB
					);
					continue;
				}

				SVGPlot.Intervals.Interval interval = intervals.addInterval(
					MathTools.log10p1(entry.lowerBinding),
					MathTools.log10p1(entry.upperBinding)
				);
				interval.id = String.format("%s: [%.2f,%.2f]",
					entry.sequence.toString(),
					MathTools.log10p1(entry.lowerBinding),
					MathTools.log10p1(entry.upperBinding)
				);
				if (entry.sequence.isWildType()) {
					interval.extraStyle = wildtypeIntervalStyle;
				}
			}
			intervals.axis = intervals.makeAxis();
			intervals.axis.tickFormat = "%.0f";
			intervals.draw(svg);
			intervals.setBounds(svg, 10, 16);

			svg.finish().write(new File("binding.kstar.svg"));
		}
	}
}
