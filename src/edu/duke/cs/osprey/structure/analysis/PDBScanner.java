package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import static edu.duke.cs.osprey.tools.Log.log;

public class PDBScanner {

	public final File dir;
	public final Set<String> blacklist;
	public final File[] files;


	public PDBScanner(File dir, String ... blacklist) {
		this(dir, new HashSet<>(Arrays.asList(blacklist)));
	}

	public PDBScanner(File dir, Set<String> blacklist) {

		this.dir = dir;
		this.blacklist = blacklist;

		this.files = dir.listFiles((file, filename) ->
			filename.endsWith(".pdb") && !blacklist.contains(filename));
	}

	/** scans the folder for molecules */
	public void scan(BiConsumer<File,Molecule> callback) {

		Progress progress = new Progress(files.length);
		log("Reading %d PDB files...", files.length);

		for (File file : files) {

			String pdbText = FileTools.readFile(file);
			for (Molecule mol : PDBIO.readAll(pdbText)) {

				callback.accept(file, mol);
			}

			progress.incrementProgress();
		}

		log("Done reading PDB files!");
	}
}
