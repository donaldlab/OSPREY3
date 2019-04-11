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

package edu.duke.cs.osprey.structure;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Streams;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;


public class TestProbe {

	@Test
	public void matchProbeContacts() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();

		Probe probe = new Probe();
		probe.matchTemplates(strand.mol.residues);

		// normalize after matching templates to our probe
		normalizeMol(strand.mol);

		// load all the overlaps from real probe output
		List<Overlap> overlaps = parseLines(strand.mol, readLines("/probe/1CC8.probe"));
		assertThat(overlaps.size(), is(3234));

		// match them to our probe results
		for (Overlap overlap : overlaps) {
			Probe.AtomPair pair = probe.new AtomPair(overlap.a1, overlap.a2);
			Probe.AtomPair.Interaction interaction = pair.getInteraction();

			try {
				assertThat(pair.getDist(), isAbsolutely(overlap.dist, 1e-3));
				assertThat(pair.infoa.vdwRadius, isAbsolutely(overlap.radius1, 1e-3));
				assertThat(pair.infob.vdwRadius, isAbsolutely(overlap.radius2, 1e-3));
				assertThat(interaction.overlap, isAbsolutely(overlap.overlap, 1e-3));
				assertThat(interaction.contact, is(overlap.contact));
			} catch (AssertionError e) {
				String msg = String.format("%s   osprey: r1=%5.3f r2=%5.3f g=%6.3f %-12s   probe: r1=%5.3f r2=%5.3f g=%6.3f %-12s",
					pair,
					pair.infoa.vdwRadius, pair.infob.vdwRadius, interaction.overlap, interaction.contact,
					overlap.radius1, overlap.radius2, overlap.overlap, overlap.contact
				);
				throw new AssertionError(msg, e);
			}
		}
	}

	@Test
	public void findAllContacts() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A23").addWildTypeRotamers();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		Probe probe = new Probe();
		probe.matchTemplates(strand.mol.residues);

		// normalize after matching templates to our probe
		normalizeMol(strand.mol);

		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.set15HasNonBonded(false)
			.build();

		// get all the probe interactions, and filter down to just the contacts
		ResidueInteractions inters = new ResidueInteractions();
		inters.addComplete(strand.mol.residues);
		List<Probe.AtomPair.Interaction> interactions = probe.getInteractions(strand.mol.residues, inters, connectivity).stream()
			.filter(i -> i.contact.isContact)
			.collect(Collectors.toList());

		// our hacked probe version actually outputs each contact twice
		assertThat(interactions.size()*2, is(3234));

		// load all the overlaps from real probe output
		List<Overlap> overlaps = parseLines(strand.mol, readLines("/probe/1CC8.probe"));
		assertThat(overlaps.size(), is(3234));

		Function<Probe.AtomPair,Boolean> findOverlap = (pair) -> {
			for (Overlap overlap : overlaps) {
				// don't need to check other direction, since the probe outputs have both directions
				if (pair.a == overlap.a1 && pair.b == overlap.a2) {
					return true;
				}
			}
			return false;
		};

		// check we matched the same atom pairs
		for (Probe.AtomPair.Interaction interaction : interactions) {
			assertThat(findOverlap.apply(interaction.atomPair), is(true));
		}
	}

	@Test
	public void regularOverlaps() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();

		Probe probe = new Probe();
		probe.matchTemplates(strand.mol.residues);

		// normalize after matching templates to our probe
		normalizeMol(strand.mol);

		// get a regular (non-hbond) atom pair
		Probe.AtomPair pair = getAtomPair(probe, strand.mol.residues, "A6 C A8 C");

		double r = 1.65;
		assertThat(pair.infoa.vdwRadius, isAbsolutely(r));
		assertThat(pair.infob.vdwRadius, isAbsolutely(r));

		// make the atoms at "easy" positions
		setAtomPos(pair.a, 0, 0, 0);
		setAtomPos(pair.b, 0, 0, 10);

		// should be no contact at this range
		assertThat(pair.getInteraction().contact, is(Probe.Contact.NoContact));
		//assertThat(probeContact(pair), is(Probe.Contact.NoContact));

		// should be right at edge of wide contact
		setAtomPos(pair.b,0, 0, 2*r + 0.5 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.NoContact));
		//assertThat(probeContact(pair), is(Probe.Contact.NoContact));
		setAtomPos(pair.b,0, 0, 2*r + 0.5 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.WideContact));
		//assertThat(probeContact(pair), is(Probe.Contact.WideContact));

		// should be right at edge of close contact
		setAtomPos(pair.b,0, 0, 2*r + 0.25 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.WideContact));
		//assertThat(probeContact(pair), is(Probe.Contact.WideContact));
		setAtomPos(pair.b,0, 0, 2*r + 0.25 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.CloseContact));
		//assertThat(probeContact(pair), is(Probe.Contact.CloseContact));

		// should be right at edge of small clash
		setAtomPos(pair.b,0, 0, 2*r + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.CloseContact));
		//assertThat(probeContact(pair), is(Probe.Contact.CloseContact));
		setAtomPos(pair.b,0, 0, 2*r - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.SmallClash));
		//assertThat(probeContact(pair), is(Probe.Contact.SmallClash));

		// should be right at edge of bad clash
		setAtomPos(pair.b,0, 0, 2*r - 0.4 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.SmallClash));
		//assertThat(probeContact(pair), is(Probe.Contact.SmallClash));
		setAtomPos(pair.b,0, 0, 2*r - 0.4 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.BadClash));
		//assertThat(probeContact(pair), is(Probe.Contact.BadClash));
	}

	@Test
	public void hbondOverlaps() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();

		Probe probe = new Probe();
		probe.matchTemplates(strand.mol.residues);

		// normalize after matching templates to our probe
		normalizeMol(strand.mol);

		// get an hbond-able atom pair
		Probe.AtomPair pair = getAtomPair(probe, strand.mol.residues, "A6 O A8 H");

		double r1 = 1.40;
		double r2 = 1.00;
		assertThat(pair.infoa.vdwRadius, isAbsolutely(r1));
		assertThat(pair.infob.vdwRadius, isAbsolutely(r2));

		// make the atoms at "easy" positions
		setAtomPos(pair.a,0, 0, 0);
		setAtomPos(pair.b,0, 0, 10);

		// should be no contact at this range
		assertThat(pair.getInteraction().contact, is(Probe.Contact.NoContact));
		//assertThat(probeContact(pair), is(Probe.Contact.NoContact));

		// should be right at edge of wide contact
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.5 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.NoContact));
		//assertThat(probeContact(pair), is(Probe.Contact.NoContact));
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.5 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.WideContact));
		//assertThat(probeContact(pair), is(Probe.Contact.WideContact));

		// should be right at edge of close contact
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.25 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.WideContact));
		//assertThat(probeContact(pair), is(Probe.Contact.WideContact));
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.25 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.CloseContact));
		//assertThat(probeContact(pair), is(Probe.Contact.CloseContact));

		// should be right at edge of the bond
		setAtomPos(pair.b,0, 0, r1 + r2 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.CloseContact));
		//assertThat(probeContact(pair), is(Probe.Contact.CloseContact));
		setAtomPos(pair.b,0, 0, r1 + r2 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.Bonded));
		//assertThat(probeContact(pair), is(Probe.Contact.Bonded));

		// should be right at edge of small clash
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.6 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.Bonded));
		//assertThat(probeContact(pair), is(Probe.Contact.Bonded));
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.6 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.SmallClash));
		//assertThat(probeContact(pair), is(Probe.Contact.SmallClash));

		// should be right at edge of bad clash
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.6 - 0.4 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.SmallClash));
		//assertThat(probeContact(pair), is(Probe.Contact.SmallClash));
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.6 - 0.4 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.BadClash));
		//assertThat(probeContact(pair), is(Probe.Contact.BadClash));
	}

	@Test
	public void saltBridgeOverlaps() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();

		Probe probe = new Probe();
		probe.matchTemplates(strand.mol.residues);

		// normalize after matching templates to our probe
		normalizeMol(strand.mol);

		// get an salt bridge-able atom pair
		Probe.AtomPair pair = getAtomPair(probe, strand.mol.residues, "A37 OD2 A35 HZ2");

		double r1 = 1.40;
		double r2 = 1.00;
		assertThat(pair.infoa.vdwRadius, isAbsolutely(r1));
		assertThat(pair.infob.vdwRadius, isAbsolutely(r2));

		// make the atoms at "easy" positions
		setAtomPos(pair.a,0, 0, 0);
		setAtomPos(pair.b,0, 0, 10);

		// should be no contact at this range
		assertThat(pair.getInteraction().contact, is(Probe.Contact.NoContact));
		//assertThat(probeContact(pair), is(Probe.Contact.NoContact));

		// should be right at edge of wide contact
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.5 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.NoContact));
		//assertThat(probeContact(pair), is(Probe.Contact.NoContact));
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.5 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.WideContact));
		//assertThat(probeContact(pair), is(Probe.Contact.WideContact));

		// should be right at edge of close contact
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.25 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.WideContact));
		//assertThat(probeContact(pair), is(Probe.Contact.WideContact));
		setAtomPos(pair.b,0, 0, r1 + r2 + 0.25 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.CloseContact));
		//assertThat(probeContact(pair), is(Probe.Contact.CloseContact));

		// should be right at edge of the bond
		setAtomPos(pair.b,0, 0, r1 + r2 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.CloseContact));
		//assertThat(probeContact(pair), is(Probe.Contact.CloseContact));
		setAtomPos(pair.b,0, 0, r1 + r2 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.Bonded));
		//assertThat(probeContact(pair), is(Probe.Contact.Bonded));

		// should be right at edge of small clash
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.8 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.Bonded));
		//assertThat(probeContact(pair), is(Probe.Contact.Bonded));
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.8 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.SmallClash));
		//assertThat(probeContact(pair), is(Probe.Contact.SmallClash));

		// should be right at edge of bad clash
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.8 - 0.4 + 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.SmallClash));
		//assertThat(probeContact(pair), is(Probe.Contact.SmallClash));
		setAtomPos(pair.b,0, 0, r1 + r2 - 0.8 - 0.4 - 1e-3);
		assertThat(pair.getInteraction().contact, is(Probe.Contact.BadClash));
		//assertThat(probeContact(pair), is(Probe.Contact.BadClash));
	}

	public static void setAtomPos(Atom a, double x, double y, double z) {

		// translate the whole residue so the atom is at the pos
		// this preserves the relative position heavy atoms bonded to hydrogens
		// so probe can get the atom types correct

		double[] coords = a.res.coords;

		int i = a.indexInRes*3;
		double dx = x - coords[i++];
		double dy = y - coords[i++];
		double dz = z - coords[i];

		for (Atom atom : a.res.atoms) {
			i = atom.indexInRes*3;
			coords[i++] += dx;
			coords[i++] += dy;
			coords[i] += dz;
		}

		normalizeCoords(coords);
	}

	public static void normalizeMol(Molecule mol) {

		// rename HIE to HIS so we don't confuse the real probe
		mol.residues.getOrThrow("A6").fullName = "HIS A   6";

		// clamp all coords to PDB-resolution, so we see exactly what probe sees
		for (Residue res : mol.residues) {
			normalizeCoords(res.coords);
		}
	}

	public static void normalizeCoords(double[] coords) {
		for (int i=0; i<coords.length; i++) {
			double c = coords[i];
			c = Double.parseDouble(String.format("%.3f", c));
			coords[i] = c;
		}
	}

	public static Probe.AtomPair getAtomPair(Probe probe, Residues residues, String desc) {
		String[] parts = desc.split(" ");
		return probe.new AtomPair(
			residues.getOrThrow(parts[0]).getAtomByName(parts[1]),
			residues.getOrThrow(parts[2]).getAtomByName(parts[3])
		);
	}

	public static class Overlap {

		public final Atom a1;
		public final Atom a2;
		public final double radius1;
		public final double radius2;
		public final double dist;
		public final double overlap;
		public final Probe.Contact contact;
		public final boolean isHbond;
		public final boolean isTooCloseHbond;

		public Overlap(Atom a1, Atom a2, double radius1, double radius2, double dist, double overlap, Probe.Contact contact, boolean isHbond, boolean isTooCloseHbond) {
			this.a1 = a1;
			this.a2 = a2;
			this.radius1 = radius1;
			this.radius2 = radius2;
			this.dist = dist;
			this.overlap = overlap;
			this.contact = contact;
			this.isHbond = isHbond;
			this.isTooCloseHbond = isTooCloseHbond;
		}

		@Override
		public String toString() {
			return String.format("%4s:%-4s %5.3f <-> %4s:%-4s %5.3f   %6.3f   %s",
				a1.res.getPDBResNumber(), a1.name, radius1,
				a2.res.getPDBResNumber(), a2.name, radius2,
				overlap,
				contact
			);
		}
	}

	public static String atomPattern(Atom atom) {
		String chainId = Character.toString(atom.res.getChainId());
		int resNum = Integer.parseInt(atom.res.fullName.substring(6, 9).trim());
		return String.format("chain%s %d atom%4s",
			chainId,
			resNum,
			StringUtils.leftPad(StringUtils.rightPad(atom.name, 3, '_'), 4, '_')
		);
	}

	public static Probe.Contact probeContact(Probe.AtomPair pair) {
		Overlap overlap = probe(pair);
		if (overlap != null) {
			return overlap.contact;
		}
		return Probe.Contact.NoContact;
	}

	public static Overlap probe(Probe.AtomPair pair) {

		List<Overlap> overlaps = probe(
			pair.a.res.molec,
			atomPattern(pair.a),
			atomPattern(pair.b)
		);
		if (overlaps.size() > 1) {
			throw new IllegalStateException("found too many overlaps: " + overlaps.size());
		} else if (overlaps.size() == 1) {
			return overlaps.get(0);
		} else {
			return null;
		}
	}

	public static List<Overlap> probe(Molecule mol) {
		return probe(mol, "-self", "all");
	}

	public static List<Overlap> probe(Molecule mol, String source, String target) {

		if (true) throw new Error("this function requires a customized version of probe, which you probably don't have");

		List<String> lines;

		// write the mol to a temporary PDB file and call probe
		try (TestBase.TempFile file = new TestBase.TempFile("probe.pdb")) {
			PDBIO.writeFile(mol, file);

			List<String> args = new ArrayList<>(Arrays.asList(
				"/home/jeff/dlab/probe.2.13.110909.src/probe", // nope, you don't have this executable. sorry =(
				"-osprey", // use the special osprey output mode
				"-4H", "-mc", "-het", // probe defaults
				source, target,
				file.getAbsolutePath()
			));
			Process process = new ProcessBuilder(args).start();

			// read the stdout and stderr streams so the buffers don't fill up
			// but only read stdout
			StreamReader out = new StreamReader(process.getInputStream());
			StreamReader err = new StreamReader(process.getErrorStream());

			// wait for probe to finish
			lines = out.getLines();

			// DEBUG
			//for (String line : err.getLines()) System.err.println("PROBE: " + line);
			//for (String line : lines) System.out.println("PROBE: " + line);

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}

		// make sure the selections found something
		// parse the first two lines:
		// source atoms: 1165
		// target atoms: 1165
		int numSource = Integer.parseInt(lines.get(0).substring(14));
		int numTarget = Integer.parseInt(lines.get(1).substring(14));
		if (numSource <= 0) {
			throw new Error("atom selection is empty: " + source);
		}
		if (numTarget <= 0) {
			throw new Error("atom selection is empty: " + target);
		}

		return parseLines(mol, lines.subList(2, lines.size()));
	}

	public static class StreamReader {

		private final List<String> lines = new ArrayList<>();
		private final Thread thread;

		public StreamReader(InputStream in) {

			// start a thread to read the stream, so we can read multiple streams at once
			thread = new Thread(() -> {

				// save all the lines
				try (BufferedReader reader = new BufferedReader(new InputStreamReader(in))) {
					String line;
					while ((line = reader.readLine()) != null) {
						lines.add(line);
					}
				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}

			});
			thread.start();
		}

		public List<String> getLines() {
			try {
				thread.join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
			return lines;
		}
	}

	public static List<Overlap> parseLines(Molecule mol, List<String> lines) {

		Function<String,Probe.Contact> findContact = (code) -> {
			switch (code) {
				case "wc": return Probe.Contact.WideContact;
				case "cc": return Probe.Contact.CloseContact;
				case "so": return Probe.Contact.SmallClash;
				case "bo": return Probe.Contact.BadClash;
				case "hb": return Probe.Contact.Bonded;
				default: throw new IllegalArgumentException("unknown contact: " + code);
			}
		};

		// parse the lines to get the overlaps
		List<Overlap> overlaps = new ArrayList<>();
		for (String line : lines) {

			// parse the line, eg:
			// A  73   O   1.400   A  72   HG3 1.170  2.830  0.260  0.260 wc 0.000 N N
			// A  73   O   1.400   A   6   HB3 1.170  2.829  0.259  0.259 wc 0.000 N N
			// A  73   C   1.650   A   6   HB3 1.170  3.240  0.420  0.420 wc 0.000 N N

			StringTokenizer tok = new StringTokenizer(line);
			String resNum1 = tok.nextToken() + Integer.parseInt(tok.nextToken());
			String atomName1 = tok.nextToken();
			double radius1 = Double.parseDouble(tok.nextToken());
			String resNum2 = tok.nextToken() + Integer.parseInt(tok.nextToken());
			String atomName2 = tok.nextToken();
			double radius2 = Double.parseDouble(tok.nextToken());
			double dist = Double.parseDouble(tok.nextToken());
			double gap = Double.parseDouble(tok.nextToken());
			double effectivegap = Double.parseDouble(tok.nextToken());
			Probe.Contact contact = findContact.apply(tok.nextToken());
			double hbondCutoff = Double.parseDouble(tok.nextToken());
			boolean isHbond = tok.nextToken().equals("Y");
			boolean isTooCloseHbond = tok.nextToken().equals("Y");

			// find the atoms
			Atom a1 = mol.residues.getOrThrow(resNum1).getAtomByName(atomName1);
			Atom a2 = mol.residues.getOrThrow(resNum2).getAtomByName(atomName2);

			// build the overlap
			overlaps.add(new Overlap(a1, a2, radius1, radius2, dist, -gap, contact, isHbond, isTooCloseHbond));
		}

		return overlaps;
	}

	public static List<String> readLines(String path) {
		return Streams.of(FileTools.parseLines(FileTools.readResource(path)))
			.collect(Collectors.toList());
	}
}
