package edu.duke.cs.osprey;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;

import org.jerkar.api.depmanagement.JkDependencies;
import org.jerkar.api.depmanagement.JkDependency;
import org.jerkar.api.depmanagement.JkDependency.JkFileDependency;
import org.jerkar.api.depmanagement.JkModuleDepFile;
import org.jerkar.api.depmanagement.JkModuleDependency;
import org.jerkar.api.depmanagement.JkModuleId;
import org.jerkar.api.depmanagement.JkResolveResult;
import org.jerkar.api.depmanagement.JkScopedDependency;
import org.jerkar.api.depmanagement.JkVersion;
import org.jerkar.api.file.JkFileTree;
import org.jerkar.api.file.JkFileTreeSet;
import org.jerkar.api.file.JkPathFilter;
import org.jerkar.api.java.JkJavaCompiler;
import org.jerkar.api.system.JkLog;
import org.jerkar.api.system.JkProcess;
import org.jerkar.api.utils.JkUtilsFile;
import org.jerkar.api.utils.JkUtilsIO;
import org.jerkar.api.utils.JkUtilsZip;
import org.jerkar.tool.JkDoc;
import org.jerkar.tool.builtins.eclipse.JkBuildPluginEclipse;
import org.jerkar.tool.builtins.javabuild.JkJavaBuild;
import org.jerkar.tool.builtins.javabuild.JkJavaPacker;

public class Build extends JkJavaBuild {
	
	@JkDoc("True (default) to make the docs, false to skip it")
	private boolean makeDocs = true;
	
	public Build() {
		
		// don't run unit tests in the build system
		// it takes forever
		tests.skip = true;
		
		manifest.mainClass = "edu.duke.cs.osprey.control.Main";
		
		// tell the eclipse plugin to stop using classpath vars
		JkBuildPluginEclipse eclipse = new JkBuildPluginEclipse();
		eclipse.useAbsolutePathsInClasspath = true;
		plugins.configure(eclipse);
	}

	@Override
	public JkModuleId moduleId() {
		return JkModuleId.of("edu.duke.cs", "osprey");
	}
	
	@Override
	public JkVersion version() {
		return JkVersion.name(JkUtilsFile.read(new File("resources/config/version")).trim());
	}

	@Override
	public String javaSourceVersion() {
		return JkJavaCompiler.V8;
	}
	
	@Override
	public JkDependencies dependencies() {
		return JkDependencies.builder()
			
			// test dependencies
			.on("org.hamcrest:hamcrest-all:1.3").scope(TEST)
			.on("junit:junit:4.12").scope(TEST)
			
			// compile dependencies (download from the repos)
			.on("colt:colt:1.2.0")
			.on("org.apache.commons:commons-math3:3.6.1")
			.on("org.apache.commons:commons-collections4:4.1")
			.on("commons-io:commons-io:2.5")
			.on("com.joptimizer:joptimizer:3.5.1")
			.on("org.ojalgo:ojalgo:41.0.0")
			.on("org.jogamp.gluegen:gluegen-rt:2.3.2")
			.on("org.jogamp.jocl:jocl:2.3.2")
			
			// apparently the JCuda repo is broken? this doesn't work:
			//.on("org.jcuda:jcuda:0.8.0")
			// just download the jar directly
			.on(url("https://search.maven.org/remotecontent?filepath=org/jcuda/jcuda/0.8.0/jcuda-0.8.0.jar"))
			
			// download TPIE-Java directly from github
			.on(url("https://github.com/donaldlab/TPIE-Java/releases/download/v1.0/edu.duke.cs.tpie-1.0.jar"))
			
			.build()
	
			// and the natives
			.and(nativeDependencies());
	}
	
	private JkDependencies nativeDependencies() {
		return JkDependencies.builder()
				
			.on("org.jogamp.gluegen:gluegen-rt:2.3.2:natives-linux-amd64").scope(RUNTIME)
			.on("org.jogamp.gluegen:gluegen-rt:2.3.2:natives-macosx-universal").scope(RUNTIME)
			.on("org.jogamp.gluegen:gluegen-rt:2.3.2:natives-windows-amd64").scope(RUNTIME)
			
			.on("org.jogamp.jocl:jocl:2.3.2:natives-linux-amd64").scope(RUNTIME)
			.on("org.jogamp.jocl:jocl:2.3.2:natives-macosx-universal").scope(RUNTIME)
			.on("org.jogamp.jocl:jocl:2.3.2:natives-windows-amd64").scope(RUNTIME)
			
			.on("org.jcuda:jcuda-natives:0.8.0:linux-x86_64").scope(RUNTIME)
			.on("org.jcuda:jcuda-natives:0.8.0:apple-x86_64").scope(RUNTIME)
			.on("org.jcuda:jcuda-natives:0.8.0:windows-x86_64").scope(RUNTIME)
			.build();
	}
	
	private File url(String url) {
		
		// get the filename
		String[] parts = url.split("/");
		String filename = parts[parts.length - 1];
		
		// pick the local file
		File libsDir = baseDir().file("build/libs");
		if (!libsDir.exists()) {
			libsDir.mkdirs();
		}
		File file = new File(libsDir, filename);
		
		// download it if needed
		if (!file.exists()) {
			try {
				JkLog.info("downloading " + url + " ...");
				JkUtilsIO.copyUrlToFile(new URL(url), file);
				JkLog.info("\t[SUCCESSFUL ]");
			} catch (MalformedURLException ex) {
				throw new RuntimeException(ex);
			}
		}
		
		return file;
	}
	
	@Override
	public JkFileTreeSet editedSources() {
		return JkFileTreeSet.of(file("src"));
	}
	
	@Override
	public JkFileTreeSet editedResources() {
		return JkFileTreeSet.of(file("resources"));
	}
	
	@Override
	public JkFileTreeSet unitTestEditedSources() {
		return JkFileTreeSet.of(file("test"));
	}
	
	@Override
	protected JkJavaPacker createPacker() {
		return JkJavaPacker.builder(this)
			.includeVersion(true)
			.doJar(false)
			.doSources(false)
			.doFatJar(true)
			.fullName(false)
			.fatJarSuffix(null)
			.fatJarExclusionFilter(JkJavaPacker.EXCLUDE_SIGNATURE_FILTER
				.and(JkPathFilter.exclude("meta-inf/license*").caseSensitive(false))
				.and(JkPathFilter.exclude("placeholder.txt"))
			)
			.extraFilesInJar(JkFileTreeSet.of(
				baseDir().include("LICENSE.txt"),
				baseDir().include("README.md")
			))
			.build();
	}
	
	public void doDist() {
		
		// make the osprey jar
		doPack();
	
		// clean the dist dir
		File dirDist = ouputDir("dist");
		dirDist.mkdirs();
		JkUtilsFile.deleteDirContent(dirDist);
		
		File dirOsprey = new File(dirDist, "osprey");
		
		// copy the python package
		JkFileTree.of(baseDir().file("python/osprey"))
			.exclude("**/*.pyc")
			.copyTo(dirOsprey);
		
		// copy install script
		JkUtilsFile.copyFileToDir(baseDir().file("python/setup.py"), dirDist);
		
		// can't exclude whole folders for some reason, so just delete it post-hoc
		JkUtilsFile.tryDeleteDir(new File(dirOsprey, "__pycache__"));
		
		// copy the osprey jar
		JkUtilsFile.copyFileToDir(packer().fatJarFile(), dirOsprey);
		
		if (makeDocs) {
			
			// make and copy the docs
			File dirDoc = baseDir().file("python/doc");
			JkProcess.of("make", "clean")
				.withWorkingDir(dirDoc)
				.failOnError(true)
				.runSync();
			JkProcess.of("make", "html")
				.withWorkingDir(dirDoc)
				.failOnError(true)
				.runSync();
			JkUtilsFile.copyDirContent(
				baseDir().file("python/doc/_build/html"),
				ouputDir("dist/doc"),
				false
			);
		}
		
		// extract the natives to a temp folder
		List<File> nativesArchives = new ArrayList<>();
		JkResolveResult resolved = dependencyResolver().resolve(RUNTIME);
		for (JkScopedDependency sdep : nativeDependencies()) {
			JkDependency dep = sdep.dependency();
			if (dep instanceof JkFileDependency) {
				JkFileDependency fileDep = (JkFileDependency)dep;
				nativesArchives.add(fileDep.files().iterator().next());
			} else if (dep instanceof JkModuleDependency) {
				JkModuleDependency modDep = (JkModuleDependency)dep;
				nativesArchives.addAll(resolved.filesOf(modDep.moduleId()));
			}
		}
		File dirTemp = ouputDir("temp");
		for (File file : nativesArchives) {
			JkUtilsZip.unzip(file, dirTemp);
		}
		
		// copy the files we want to the dist/natives folder
		File dirNatives = new File(dirOsprey, "natives");
		JkFileTree.of(dirTemp)
			.include("**/*.dll", "**/*.so", "**/*.*lib")
			.forEach((file) -> {
				JkUtilsFile.copyFileToDir(file, dirNatives);
			});
		
		// cleanup the temp folder
		JkUtilsFile.deleteDir(dirTemp);
		
		// update the osprey dev flag
		File initFile = new File(dirOsprey, "__init__.py");
		String initFileContent = JkUtilsFile.read(initFile);
		initFileContent = initFileContent.replaceFirst("_IS_DEV = True\\n", "_IS_DEV = False");
		JkUtilsFile.writeString(initFile, initFileContent, false);
		
		// copy the examples
		JkUtilsFile.copyDirContent(baseDir().file("examples"), ouputDir("dist/examples"), true);
		
		// copy text files
		JkUtilsFile.copyFileToDir(baseDir().file("LICENSE.txt"), dirDist);
		JkUtilsFile.copyFileToDir(baseDir().file("README.rst"), dirDist);
		JkUtilsFile.copyFileToDir(baseDir().file("CONTRIBUTING.rst"), dirDist);
		
		// create the distribution zip
		JkFileTree.of(dirDist).zip().to(ouputDir(packer().fatJarFile().getName().replaceAll(".jar", ".zip")));
	}
	
	public void doClasspath()
	throws IOException {
		
		// start with the bin dir
		LinkedHashSet<File> files = new LinkedHashSet<>();
		files.add(baseDir().file("bin"));
		
		// collect all the dependency files in order
		JkResolveResult result = dependencyResolver().resolve(RUNTIME);
		for (JkScopedDependency sdep : dependencies()) {
			JkDependency dep = sdep.dependency();
			if (dep instanceof JkFileDependency) {
				files.addAll(((JkFileDependency)dep).files());
			} else if (dep instanceof JkModuleDependency) {
				for (JkModuleDepFile depfile : result.moduleFiles()) {
					files.add(depfile.localFile());
				}
			}
		}
		
		// write out the runtime classpath for the python scripts dev mode
		try (FileWriter out = new FileWriter(ouputDir().file("classpath.txt"))) {
			for (File file : files) {
				out.write(file.getAbsolutePath());
				out.write("\n");
			}
		}
	}
	
	public void setupDevEnv()
	throws IOException {
		
		// build eclipse classpath
		JkBuildPluginEclipse eclipse = plugins.findInstanceOf(JkBuildPluginEclipse.class);
		eclipse.generateFiles();
		
		// build python classpath
		doClasspath();
	}
}
