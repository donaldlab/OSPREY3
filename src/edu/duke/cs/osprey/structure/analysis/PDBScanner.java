package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;

import static edu.duke.cs.osprey.tools.Log.log;


public class PDBScanner {

	public final File dir;
	public final Set<String> blacklist;
	public final List<File> files;


	public PDBScanner(File dir, String ... blacklist) {
		this(dir, new HashSet<>(Arrays.asList(blacklist)));
	}

	public PDBScanner(File dir, Set<String> blacklist) {

		this.dir = dir;
		this.blacklist = blacklist;

		this.files = Arrays.asList(dir.listFiles((file, filename) ->
				filename.endsWith(".pdb") && !blacklist.contains(filename)
			));
	}

	/** scans the first N files in the folder for molecules */
	public void scan(int numFiles, BiConsumer<File,Molecule> callback) {
		scan(this.files.subList(0, numFiles), true, callback);
	}

	/** scans the folder for molecules */
	public void scan(BiConsumer<File,Molecule> callback) {
		scan(this.files, true, callback);
	}

	/** scan just the named file */
	public void scan(String filename, BiConsumer<File,Molecule> callback) {
		scan(Arrays.asList(new File(dir, filename)), false, callback);
	}

	private void scan(List<File> files, boolean showProgress, BiConsumer<File,Molecule> callback) {

		Progress progress = null;
		if (showProgress) {
			progress = new Progress(files.size());
			log("Reading %d PDB files...", files.size());
		}

		for (File file : files) {

			String pdbText = FileTools.readFile(file);

			// try to read the PDB file, or just skip it
			List<Molecule> mols;
			try {
				mols = PDBIO.readAll(pdbText);
			} catch (Exception ex) {
				System.err.println("error reading PDB file " + file.getName() + ", skipping it:\n\t" + ex.getMessage());
				continue;
			}

			// pass each model to the callback
			for (Molecule mol : mols) {
				callback.accept(file, mol);
			}

			if (showProgress) {
				progress.incrementProgress();
			}
		}

		if (showProgress) {
			log("Done reading PDB files!");
		}
	}
}
