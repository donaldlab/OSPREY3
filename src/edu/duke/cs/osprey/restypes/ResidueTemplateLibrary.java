package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.Serializable;
import java.util.*;

/**
 * Library of {@link ResidueTemplate} instances used by {@link Strand} instances to define discrete flexibilty.
 *
 * See Python script examples/python.GMEC/templateLibrary.py in your Osprey distribution for example usage.
 */
// NOTE: we should not be serializing this class, but DEEPer expects just about everything to be serializable
// objects referring to this class get serialized, taking the library with it
// when we deserialize said objects, we'll get two different template libraries:
// one from the deserialization, and one from the context of the code that called the deserializer
// ideally, we'd have context-aware deserialization, but that's a bit more work than I want to do today
// not to mention that stack overflow errors that Java's serialization system is prone to generate... oi
public class ResidueTemplateLibrary implements Serializable {

	public static class Builder {

		private static final String LovellRotamersPath = "/config/LovellRotamer.dat";

		@SafeVarargs
		private static <T> List<T> listOf(T ... vals) {
			List<T> list = new ArrayList<>();
			list.addAll(Arrays.asList(vals));
			return list;
		}

		/**
		 * The templates contain forcefield information like atom types,
		 * so they need to be matched with a specific forcefield.
		 */
		private ForcefieldParams.Forcefield forcefield = ForcefieldParams.Forcefield.AMBER;

		private final List<String> templatesTexts = listOf();
		private final List<String> templateCoordsTexts = listOf();
		private final List<String> rotamersTexts = listOf();
		private final List<String> backboneDependentRotamersTexts = listOf();
		private final List<String> entropyTexts = listOf();

		/**
		 * True to generate D amino acids for the template library
		 */
		private boolean makeDAminoAcidTemplates = true;

		private List<Molecule> molsForRotamers = listOf();

		public Builder() {
			initTexts();
		}

		public Builder(ForcefieldParams.Forcefield forcefield) {
			this.forcefield = forcefield;
			initTexts();
		}

		private void initTexts() {
			templatesTexts.add(FileTools.readResource(forcefield.aaPath));
			templatesTexts.add(FileTools.readResource(forcefield.aaNTPath));
			templatesTexts.add(FileTools.readResource(forcefield.aaCTPath));
			templatesTexts.add(FileTools.readResource(forcefield.grPath));
			templateCoordsTexts.add(FileTools.readResource("/config/all_amino_coords.in"));
			rotamersTexts.add(FileTools.readResource(LovellRotamersPath));
			entropyTexts.add(FileTools.readResource("/config/ResEntropy.dat"));
		}

		/**
		 * By default, templates for natual amino acids are added
		 */
		public Builder clearTemplates() {
			templatesTexts.clear();
			return this;
		}

		/**
		 * Templates to add to the library
		 */
		public Builder addTemplates(String text) {
			templatesTexts.add(text);
			return this;
		}

		/**
		 * By default, template coordinates for natual amino acids are added
		 */
		public Builder clearTemplateCoords() {
			templateCoordsTexts.clear();
			return this;
		}

		/**
		 * Template coordinates to add to the library
		 *
		 * @todo explain template coord file format. anyone know about this?
		 */
		public Builder addTemplateCoords(String text) {
			templateCoordsTexts.add(text);
			return this;
		}

		/**
		 * By default, the Lovell rotamer library is added.
		 *
		 * See {@link #addLovellRotamers()}
		 */
		public Builder clearRotamers() {
			rotamersTexts.clear();
			return this;
		}

		/**
		 * Rotamers to add to the library
		 *
		 * @todo explain rotamer file format. anyone know about this?
		 */
		public Builder addRotamers(String text) {
			rotamersTexts.add(text);
			return this;
		}

		/**
		 * Adds rotamers from the pentiultimate rotamer library
		 * {@cite Lovell2000 Lovell, Lovell, S.C., Word, J.M., Richardson, J.S. and Richardson, D.C., 2000.
		 * The penultimate rotamer library. Proteins: Structure, Function, and Bioinformatics, 40(3), pp.389-408.}
		 * is used.
		 *
		 * @todo explain rotamers file format. anyone know about this?
		 */
		public Builder addLovellRotamers() {
			return addRotamers(FileTools.readResource(LovellRotamersPath));
		}

		/**
		 * Backbone-dependent rotamers to add to the library
		 *
		 * @todo explain backbone-dependent rotamer file format. anyone know about this?
		 */
		public Builder addBackboneDependentRotamers(String text) {
			backboneDependentRotamersTexts.add(text);
			return this;
		}

		/**
		 * By default, residue entropies are added for natual amio acids
		 */
		public Builder clearResidueEntropies() {
			entropyTexts.clear();
			return this;
		}

		/**
		 * Residue entropies to add to the library
		 *
		 * @todo explain residue entropy file format. anyone know about this?
		 */
		public Builder addResidueEntropies(String text) {
			entropyTexts.add(text);
			return this;
		}

		public Builder setMakeDAminoAcidTemplates(boolean val) {
			makeDAminoAcidTemplates = val;
			return this;
		}

		/**
		 * Add a molecule to include its wild-type rotamers in the library
		 */
		public Builder addMoleculeForWildTypeRotamers(Molecule val) {
			molsForRotamers.add(val);
			return this;
		}

		public ResidueTemplateLibrary build() {
			return new ResidueTemplateLibrary(
					forcefield,
					templatesTexts,
					templateCoordsTexts,
					rotamersTexts,
					backboneDependentRotamersTexts,
					entropyTexts,
					makeDAminoAcidTemplates,
					molsForRotamers
			);
		}
	}

	// NAMING: We assume each distinct residue (AA or otherwise) has its own name
	// however, many residues will have multiple slightly different forms (N-terminal, etc.)
	// and these will share a name and a rotamer library entry

	public final ForcefieldParams ffparams;
	public final ArrayList<ResidueTemplate> templates = new ArrayList<>();
	public final Map<String,ResidueTemplate> wildTypeTemplates = new HashMap<>();
	public final ResidueEntropies residueEntropies = new ResidueEntropies();
	public int totalNumRotamers;//total number of rotamers read in from rotamer library file(s), starts at 0

	public ResidueTemplateLibrary(ForcefieldParams.Forcefield forcefield, List<String> templatesTexts, List<String> templateCoordTexts, List<String> rotamersTexts, List<String> backboneDependentRotamerTexts, List<String> resEntropyTexts, boolean makeDAminoAcids, List<Molecule> molsForRotamers) {

		this.ffparams = new ForcefieldParams(forcefield);

		// load templates
		TemplateParser templateParser = new TemplateParser(ffparams);
		for (String text : templatesTexts) {
			templates.addAll(templateParser.parse(text));
		}

		// load template coords
		TemplateCoordsParser templateCoordsParser = new TemplateCoordsParser(templates);
		for (String text : templateCoordTexts) {
			templateCoordsParser.parse(text);
		}

		// load rotamers
		int numRotamers = 0;
		for (String text : rotamersTexts) {
			numRotamers += RotamerLibraryReader.readRotLibrary(text, templates);
		}
		this.totalNumRotamers = numRotamers;
		for (String text : backboneDependentRotamerTexts) {
			RotamerLibraryReader.readDunbrackRotamerLibraryForResiduePosition(text, templates);
		}

		// make D amino acids if needed
		// NOTE: make D amino acids after loading all templates, coords, and rotamers
		if (makeDAminoAcids) {
			templates.addAll(DAminoAcidHandler.makeDTemplates(templates));
		}

		// load residue entropies
		ResidueEntropyParser entropyParser = new ResidueEntropyParser(residueEntropies);
		for (String text : resEntropyTexts) {
			entropyParser.parse(text);
		}

		// make wild type rotamers
		for (Molecule mol : molsForRotamers) {
			Strand strand = new Strand.Builder(mol)
				.setTemplateLibrary(this)
				.build();
			for (Residue res : strand.mol.residues) {
				getOrMakeWildTypeTemplate(res);
			}
		}
	}

	public ResidueTemplate getTemplate(String resType) {
		return getTemplate(resType, false);
	}

	public ResidueTemplate getTemplate(String resType, boolean requireCoords) {
		//get the first template with the given residue type
		//templates with the same type will all have the same rotamers, etc.
		//so this is good for looking up things like that
		for (ResidueTemplate templ : templates) {
			if (templ.name.equalsIgnoreCase(resType)) {
				if (!requireCoords || templ.templateRes.coords != null) {
					return templ;
				}
			}
		}
		return null;
	}

	public ResidueTemplate getTemplateOrThrow(String resType) {
		return getTemplateOrThrow(resType, false);
	}

	public ResidueTemplate getTemplateOrThrow(String resType, boolean requireCoords) {
		ResidueTemplate templ = getTemplate(resType, requireCoords);
		if (templ != null) {
			return templ;
		}
		throw new NoSuchElementException("no residue template with name: " + resType);
	}

	public ResidueTemplate getTemplateForMutation(String resType, Residue res) {
		try {
			return getTemplateOrThrow(resType, true);
		} catch (NoSuchElementException ex) {
			throw new NoSuchElementException("ERROR: Couldn't find a template for mutating " + res.fullName + " to " + resType);
		}
	}

	public double getResEntropy(String resType){
		return residueEntropies.get(resType);
	}

	public int numDihedralsForResType(String resType) {
		//number of free dihedrals (sidechain dihedrals for actual amino acids)
		return getTemplateOrThrow(resType).numDihedrals;
	}

	/**
	 * PGC 2015: 
	 * Returns the number of rotamers for the specified residue type, for backbone dependent or backbone independent rotamers.
	 * @param pos 
	 * @param resType in three letter amino acid type
	 * @param phi The backbone phi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @param psi The backbone psi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @return The number of rotamers.  
	 */
	public int numRotForResType(int pos, String resType, double phi, double psi) {
		return getTemplateOrThrow(resType).getNumRotamers(phi, psi);
	}

	/**
	 * PGC 2015:
	 * get ideal dihedral value for a particular rotamer of a particular residue type, for backbone dependent or backbone independent rotamers.
	 * 
	 * @param resType in three letter amino acid type
	 * @param phi The backbone phi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @param psi The backbone psi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @param rotNum The rotamer number within this residue type.
	 * @param dihedralNum The dihedral number within this rotamer.
	 * @return
	 */
	public double getDihedralForRotamer(int pos, String resType, double phi, double psi, int rotNum, int dihedralNum) {
		return getTemplateOrThrow(resType).getRotamericDihedrals(phi, psi, rotNum, dihedralNum);
	}
	public double getDihedralForRotamer(int pos, String resType, int rotNum, int dihedralNum) {
		return getTemplateOrThrow(resType).getRotamericDihedrals(0, 0, rotNum, dihedralNum);
	}
	public double getDihedralForRotamer(String resType, int rotNum, int dihedralNum) {
		return getTemplateOrThrow(resType).getRotamericDihedrals(0, 0, rotNum, dihedralNum);
	}

	public ResidueTemplate getOrMakeWildTypeTemplate(Residue res) {

		// check the library first
		String resNum = res.getPDBResNumber();
		ResidueTemplate template = wildTypeTemplates.get(resNum);
		if (template != null) {
			return template;
		}

		// not found, make a new one
		List<Residue> residues = new ArrayList<>();
		residues.add(res);
		residues.addAll(res.molec.getAlternates(res.indexInMolecule));

		template = ResidueTemplate.makeFromResidueConfs(residues);
		wildTypeTemplates.put(resNum, template);
		return template;
	}
	
    public boolean containsName(String name){
        //Is there a template with the specified name?
        for(ResidueTemplate templ : templates){
            if(templ.name.equalsIgnoreCase(name))
                return true;
        }
        return false;
    }

	public HashSet<String> templateNameSet(){
		HashSet<String> ans = new HashSet<>();
		for(ResidueTemplate templ : templates)
			ans.add(templ.name);
		return ans;
	}
}
