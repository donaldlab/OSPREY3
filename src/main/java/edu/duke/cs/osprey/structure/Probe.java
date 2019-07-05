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

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.ConfigFileReader;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Streams;

import java.util.*;
import java.util.stream.Collectors;


/**
 * implements probe-style rules for describing non-bonded atomic interactions,
 * eg close contacts, clashes, hydrogen bonds, salt bridges
 *
 * probe is software available at:
 * http://kinemage.biochem.duke.edu/software/probe.php
 */
public class Probe {

	public static enum AtomFlag {
		Donor,
		Acceptor,
		Positive,
		Negative
	}

	public static class AtomInfo {

		public final String name;
		public final double vdwRadius;
		public final EnumSet<AtomFlag> flags;

		public AtomInfo(String name, double vdwRadius, List<AtomFlag> flags) {

			this.name = name;
			this.vdwRadius = vdwRadius;

			this.flags = EnumSet.noneOf(AtomFlag.class);
			this.flags.addAll(flags);
		}
	}

	public static class Template {

		public final String id;
		public final String resType;
		public final String classifier;
		public final Map<String,AtomInfo> atoms;

		public Template(String id, String resType, String classifier, List<AtomInfo> atomInfos) {

			this.id = id;
			this.resType = resType;
			this.classifier = classifier;

			this.atoms = new HashMap<>();
			for (AtomInfo info : atomInfos) {
				this.atoms.put(info.name, info);
			}
		}
	}


	public double maxOverlapHBond = Double.NaN;
	public double maxOverlapSaltBridge = Double.NaN;
	public double minOverlapWideContact = Double.NaN;
	public double minOverlapCloseContact = Double.NaN;
	public double minOverlapBadClash = Double.NaN;

	private Map<String,List<Template>> templatesByResType = new HashMap<>();
	private Map<ResidueTemplate,Template> templateMap = new IdentityHashMap<>();

	public Probe() {
		this(true);
	}

	public Probe(boolean loadDefaults) {
		if (loadDefaults) {
			load(FileTools.readResource("/config/probe.cfg"));
		}
	}

	/** loads configuration info (eg, params, templates) into this instance */
	public void load(String configText) {

		// read the sections, eg params, templates
		ConfigFileReader reader = new ConfigFileReader(configText);
		reader.advanceToNonEmptyLine();
		while (reader.getLine() != null) {

			String section = reader.getSectionName();
			if (section.equalsIgnoreCase("params")) {
				loadParams(reader);
			} else if (section.equalsIgnoreCase("templates")) {
				loadTemplates(reader);
			} else {
				throw new IllegalArgumentException("unknown section: " + section);
			}
		}
	}

	private void loadParams(ConfigFileReader reader) {
		while (true) {

			reader.advanceToNonEmptyLine();
			if (reader.getLine() == null || reader.getSectionName() != null) {
				break;
			}

			try {
				reader.getAssignment((name, value) -> {
					if (name.equalsIgnoreCase("MaxOverlapHBond")) {
						maxOverlapHBond = Double.parseDouble(value);
					} else if (name.equalsIgnoreCase("MaxOverlapSaltBridge")) {
						maxOverlapSaltBridge = Double.parseDouble(value);
					} else if (name.equalsIgnoreCase("MinOverlapWideContact")) {
						minOverlapWideContact = Double.parseDouble(value);
					} else if (name.equalsIgnoreCase("MinOverlapCloseContact")) {
						minOverlapCloseContact = Double.parseDouble(value);
					} else if (name.equalsIgnoreCase("MinOverlapBadClash")) {
						minOverlapBadClash = Double.parseDouble(value);
					}
				});
			} catch (NumberFormatException ex) {
				throw new IllegalArgumentException("can't parse param value: " + reader.getLine());
			}
		}
	}

	private void loadTemplates(ConfigFileReader reader) {
		while (true) {

			reader.advanceToNonEmptyLine();
			if (reader.getLine() == null || reader.getSectionName() != null) {
				break;
			}

			// read the template, eg
			// ResType-chainPos
			//     AtomName vdwRadius flags
			// (blank line)

			// parse the id into res type and chain pos
			String id = reader.getLine();
			String[] parts = id.split("-");
			String resType = parts[0];
			String classifier = null;
			if (parts.length >= 2) {
				classifier = parts[1];
			}

			// read the atom infos
			List<AtomInfo> atomInfos = new ArrayList<>();
			while (true) {

				reader.advance();
				if (reader.getLine() == null || reader.getLine().isEmpty()) {
					break;
				}

				StringTokenizer tok = new StringTokenizer(reader.getLine(), " \t");
				try {

					String atomName = tok.nextToken();
					double vdwRadius = Double.parseDouble(tok.nextToken());

					// flags are optional
					List<AtomFlag> flags = new ArrayList<>();
					if (tok.hasMoreTokens()) {
						String flagsStr = tok.nextToken();

						// parse the flags
						if (flagsStr.startsWith("-")) {
							flags.add(AtomFlag.Negative);
							flagsStr = flagsStr.substring(1);
						} else if (flagsStr.startsWith("+")) {
							flags.add(AtomFlag.Positive);
							flagsStr = flagsStr.substring(1);
						}

						if (flagsStr.equalsIgnoreCase("donor")) {
							flags.add(AtomFlag.Donor);
						} else if (flagsStr.equalsIgnoreCase("acceptor")) {
							flags.add(AtomFlag.Acceptor);
						}
					}

					// make the atom info
					atomInfos.add(new AtomInfo(atomName, vdwRadius, flags));

				} catch (NoSuchElementException | NumberFormatException ex) {
					throw new IllegalArgumentException("can't parse atom record: " + reader.getLine());
				}

				// make the template
				templatesByResType
					.computeIfAbsent(resType, (resTypeAgain) -> new ArrayList<>())
					.add(new Template(id, resType, classifier, atomInfos));
			}
		}
	}

	/**
	 * finds a probe template for each template in the residues and remembers the mapping
	 */
	public boolean matchTemplates(Residues residues) {

		Set<ResidueTemplate> templates = new HashSet<>();
		for (Residue res : residues) {
			templates.add(res.template);
		}

		boolean allMatched = true;
		for (ResidueTemplate template : templates) {
			allMatched &= matchTemplate(template);
		}
		return allMatched;
	}

	/**
	 * finds a probe template for each template in the library and remembers the mapping
	 */
	public boolean matchTemplates(ResidueTemplateLibrary templateLib) {
		return matchTemplates(templateLib.templates)
			&& matchTemplates(templateLib.wildTypeTemplates.values());
	}

	/**
	 * finds a probe template for each template in the library and remembers the mapping
	 */
	public boolean matchTemplates(Collection<ResidueTemplate> templates) {
		boolean allMatched = true;
		for (ResidueTemplate templ : templates) {
			allMatched &= matchTemplate(templ);
		}
		return allMatched;
	}

	/**
	 * finds a probe template for each template in the conf space and remembers the mapping
	 */
	public boolean matchTemplates(SimpleConfSpace confSpace) {

		// collect all the templates in the conf space
		Set<ResidueTemplate> templates = new HashSet<>();
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
				 templates.add(rc.template);
			}
		}
		for (Strand strand : confSpace.strands) {
			for (Residue res : strand.mol.residues) {
				templates.add(res.template);
			}
		}

		boolean allMatched = true;
		for (ResidueTemplate template : templates) {
			allMatched &= matchTemplate(template);
		}
		return allMatched;
	}

	/**
	 * finds the templates that match the residue type,
	 * then tries to match a template to the atom names in a residue
	 * (should match eg N-terminal or C-terminal amino acids if needed)
	 */
	public boolean matchTemplate(ResidueTemplate resTemplate) {

		// look up templates by residue type
		List<Template> templates = templatesByResType.get(resTemplate.name);
		if (templates == null) {
			return false;
		}

		// make a quick set lookup for the residue atom names
		Set<String> resAtomNames = new HashSet<>();
		for (Atom atom : resTemplate.templateRes.atoms) {
			resAtomNames.add(atom.name);
		}

		// try to match the atom names to one of the templates for this res type
		for (Template template : templates) {
			if (resAtomNames.equals(template.atoms.keySet())) {

				// found a match! remember it
				templateMap.put(resTemplate, template);

				return true;
			}
		}

		return false;
	}

	public Template getTemplate(Residue res) {

		if (res.template == null) {
			throw new IllegalArgumentException("residue " + res + " does not have a template");
		}

		Template template = templateMap.get(res.template);
		if (template == null) {
			Set<String> atomNames = res.template.templateRes.atoms.stream()
				.map(atom -> atom.name)
				.collect(Collectors.toSet());
			throw new IllegalArgumentException(
				"no probe template for residue template " + res.template
				+ " with atoms " + Streams.joinToString(atomNames, ",")
			);
		}

		return template;
	}

	public AtomInfo getAtomInfo(Atom a) {

		if (a.res == null) {
			throw new IllegalArgumentException("atom " + a.name + " does not have a residue");
		}
		Template template = getTemplate(a.res);

		AtomInfo info = template.atoms.get(a.name);
		if (info == null) {
			throw new IllegalArgumentException("no probe atom info for " + a.name + " in template " + template.id);
		}

		return info;
	}

	public static enum Attraction {

		None(false),
		Hbond(true),
		SaltBridge(true);

		public final boolean isBond;

		private Attraction(boolean isBond) {
			this.isBond = isBond;
		}
	}

	public static enum Contact {

		NoContact(false, false, false),
		WideContact(true, false, false),
		CloseContact(true, false, false),
		SmallClash(true, true, false),
		BadClash(true, true, false),
		Bonded(true, false, true);

		public final boolean isContact;
		public final boolean isClash;
		public final boolean isBond;

		private Contact(boolean isContact, boolean isClash, boolean isBond) {
			this.isContact = isContact;
			this.isClash = isClash;
			this.isBond = isBond;
		}
	}

	public class AtomPair {

		public class Interaction {

			public final AtomPair atomPair = AtomPair.this;
			public final double dist;
			public final double overlap;
			public final Contact contact;

			public Interaction(double dist, double overlap, Contact contact) {
				this.dist = dist;
				this.overlap = overlap;
				this.contact = contact;
			}

			public double getViolation(double tolerance) {
				return probe.getViolation(overlap, maxOverlap, tolerance);
			}

			public boolean isClash(double tolerance) {
				return getViolation(tolerance) > 0;
			}

			@Override
			public String toString() {
				return String.format("%-10s overlap=%8.3f %s%5.3f%12s",
					attraction == Attraction.None ? "" : attraction.name(),
					overlap,
					overlap > maxOverlap ? " >" : "<=",
					maxOverlap,
					contact == Contact.NoContact ? "" : " " + contact.name()
				);
			}
		}

		public final Probe probe = Probe.this;

		public final Atom a;
		public final Atom b;

		public final AtomInfo infoa;
		public final AtomInfo infob;
		public final Attraction attraction;
		public final double maxOverlap;

		public AtomPair(Atom a, Atom b) {
			this(a, b, getAtomInfo(a), getAtomInfo(b));
		}

		public AtomPair(Atom a, Atom b, AtomInfo infoa, AtomInfo infob) {

			this.a = a;
			this.b = b;
			this.infoa = infoa;
			this.infob = infob;

			this.attraction = getAttraction(infoa, infob);
			this.maxOverlap = getMaxOverlap(attraction);
		}

		public double getDist() {
			return probe.getDist(a, b);
		}

		public Interaction getInteraction() {
			double dist = getDist();
			double overlap = getOverlap(infoa, infob, dist);
			Contact contact = getContact(overlap);
			return new Interaction(dist, overlap, contact);
		}

		public Contact getContact(double overlap) {
			return probe.getContact(overlap, maxOverlap, attraction);
		}

		/** return the overlap of the two atoms beyond the max allowed overlap, and beyond any extra tolerance */
		public double getViolation(double tolerance) {
			double dist = getDist();
			double overlap = getOverlap(infoa, infob, dist);
			return probe.getViolation(overlap, maxOverlap, tolerance);
		}

		@Override
		public String toString() {
			return String.format("%s:%-4s <-> %s:%-4s", a.res.getPDBResNumber(), a.name, b.res.getPDBResNumber(), b.name);
		}
	}

	public double getViolation(double overlap, double maxOverlap, double tolerance) {
		return overlap - maxOverlap - tolerance;
	}

	public double getDistSq(Atom a, Atom b) {

		// get the interatomic distance
		int i = a.indexInRes*3;
		double ax = a.res.coords[i];
		double ay = a.res.coords[i + 1];
		double az = a.res.coords[i + 2];
		i = b.indexInRes*3;
		double bx = b.res.coords[i];
		double by = b.res.coords[i + 1];
		double bz = b.res.coords[i + 2];

		// calc the squared distance
		double dx = ax - bx;
		double dy = ay - by;
		double dz = az - bz;
		return dx*dx + dy*dy + dz*dz;
	}

	public double getDist(Atom a, Atom b) {
		return Math.sqrt(getDistSq(a, b));
	}

	public Attraction getAttraction(AtomInfo a, AtomInfo b) {

		// is there any attraction to oppose the vdW forces?
		boolean donorAcceptorMatch =
			(a.flags.contains(AtomFlag.Donor) && b.flags.contains(AtomFlag.Acceptor))
				|| (a.flags.contains(AtomFlag.Acceptor) && b.flags.contains(AtomFlag.Donor));
		if (!donorAcceptorMatch) {

			// no donor,acceptor match, can't be an hbond
			return Attraction.None;

		} else {

			boolean bothCharged =
				(a.flags.contains(AtomFlag.Positive) || a.flags.contains(AtomFlag.Negative))
					&& (b.flags.contains(AtomFlag.Positive) || b.flags.contains(AtomFlag.Negative));
			if (bothCharged) {
				boolean chargeComplementarity =
					(a.flags.contains(AtomFlag.Positive) && b.flags.contains(AtomFlag.Negative))
						|| (a.flags.contains(AtomFlag.Negative) && b.flags.contains(AtomFlag.Positive));
				if (chargeComplementarity) {
					// donor,acceptor match + ionic attraction = salt bridge
					return Attraction.SaltBridge;
				} else {
					// donor,acceptor match + ionic repulsion != hbond
					return Attraction.None;
				}
			} else {
				// donor,acceptor match but no ionic effects = hbond
				return Attraction.Hbond;
			}
		}
	}

	public double getOverlap(AtomInfo a, AtomInfo b, double dist) {
		return a.vdwRadius + b.vdwRadius - dist;
	}

	public double getMaxOverlap(Attraction attraction) {
		if (attraction == Attraction.Hbond) {
			return maxOverlapHBond;
		} else if (attraction == Attraction.SaltBridge) {
			return maxOverlapSaltBridge;
		} else {
			return 0.0;
		}
	}

	public Contact getContact(double overlap, double maxOverlap, Attraction attraction) {

		if (overlap < minOverlapWideContact) {
			return Contact.NoContact;
		} else if (overlap < minOverlapCloseContact) {
			return Contact.WideContact;
		} else if (overlap <= 0) {
			return Contact.CloseContact;
		} else {
			// we have positive overlap, clash depends if hbond or not
			if (attraction.isBond) {
				// we have an hbond, is bonded or clashing?
				if (overlap <= maxOverlap) {
					return Contact.Bonded;
				} else {
					// adjust overlap by the allowed hbond amount
					if (overlap - maxOverlap < minOverlapBadClash) {
						return Contact.SmallClash;
					} else {
						return Contact.BadClash;
					}
				}
			} else {
				// non-bonded, treat as clash
				if (overlap < minOverlapBadClash) {
					return Contact.SmallClash;
				} else {
					return Contact.BadClash;
				}
			}
		}
	}

	/**
	 * gets probe interactions between all atom pairs defined by the residue pairs
	 *
	 * atom connectivity must have set15HasNonBonded(false) to match real probe reults
	 */
	public List<AtomPair.Interaction> getInteractions(Residues residues, ResidueInteractions inters, AtomConnectivity connectivity) {

		List<AtomPair.Interaction> interactions = new ArrayList<>();

		// for each res pair
		for (ResidueInteractions.Pair resPair : inters) {
			Residue res1 = residues.getOrThrow(resPair.resNum1);
			Residue res2 = residues.getOrThrow(resPair.resNum2);

			// for each atom pair
			for (int[] atomPair : connectivity.getAtomPairs(res1, res2).getPairs(AtomNeighbors.Type.NONBONDED)) {
				Atom a1 = res1.atoms.get(atomPair[0]);
				Atom a2 = res2.atoms.get(atomPair[1]);

				interactions.add(new AtomPair(a1, a2).getInteraction());
			}
		}

		return interactions;
	}

	/**
	 * gets probe interactions between all atom pairs within the residue
	 *
	 * atom connectivity must have set15HasNonBonded(false) to match real probe reults
	 */
	public List<AtomPair.Interaction> getInteractions(Residue res, AtomConnectivity connectivity) {

		List<AtomPair.Interaction> interactions = new ArrayList<>();

		// for each atom pair
		for (int[] atomPair : connectivity.getAtomPairs(res, res).getPairs(AtomNeighbors.Type.NONBONDED)) {
			Atom a1 = res.atoms.get(atomPair[0]);
			Atom a2 = res.atoms.get(atomPair[1]);

			interactions.add(new AtomPair(a1, a2).getInteraction());
		}

		return interactions;
	}
}
