package edu.duke.cs.osprey;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
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
			.on(url("https://github.com/donaldlab/TPIE-Java/releases/download/v1.1/edu.duke.cs.tpie-1.1.jar"))
			
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
	
		File dirOut = ouputDir(".");
		dirOut.mkdirs();
		
		// make the docs if needed
		File dirDocHtml = null;
		if (makeDocs) {
			File dirDoc = baseDir().file("python/doc");
			JkProcess.of("make", "clean")
				.withWorkingDir(dirDoc)
				.failOnError(true)
				.runSync();
			JkProcess.of("make", "html")
				.withWorkingDir(dirDoc)
				.failOnError(true)
				.runSync();
			dirDocHtml = new File(dirDoc, "_build/html");
		}
		
		// extract the natives
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
		
		// filter the natives dir to just the libs
		File dirNatives = ouputDir("natives");
		JkFileTree.of(dirTemp)
			.include("**/*.dll", "**/*.so", "**/*.*lib")
			.forEach((file) -> {
				JkUtilsFile.copyFileToDir(file, dirNatives);
			});
		
		JkUtilsFile.tryDeleteDir(dirTemp);
		
		// create the release name
		String name = "osprey-" + version();
		
		javaDist(name, dirOut, dirDocHtml, dirNatives);
		pythonDist(name, dirOut, dirDocHtml, dirNatives);
		
		// cleanup temp files
		JkUtilsFile.tryDeleteDir(dirNatives);
	}
	
	private static List<String> javaExamples = Arrays.asList(
		"1CC8", "1FSV", "2O9S_L2", "2RL0.kstar", "3K75.3LQC", "4HEM", "4NPD"
	);
	
	private void javaDist(String name, File dirOut, File dirDocHtml, File dirNatives) {
		
		File dir = new File(dirOut, "java");
		dir.mkdirs();
	
		// copy docs if needed
		if (dirDocHtml != null) {
			JkUtilsFile.copyDirContent(dirDocHtml, new File(dir, "doc"), false);
		}
		
		// copy the osprey jar
		JkUtilsFile.copyFileToDir(packer().fatJarFile(), dir);
		
		// copy the natives
		JkUtilsFile.copyDirContent(dirNatives, new File(dir, "natives"), false);
		
		// copy the examples
		copyExamples(javaExamples, new File(dir, "examples"));
		
		// copy text files
		JkUtilsFile.copyFileToDir(baseDir().file("LICENSE.txt"), dir);
		JkUtilsFile.copyFileToDir(baseDir().file("README.rst"), dir);
		JkUtilsFile.copyFileToDir(baseDir().file("CONTRIBUTING.rst"), dir);
		
		// create the distribution zip
		JkFileTree.of(dir).zip().to(new File(dirOut, name + ".java.zip"));
		
		// cleanup temp files
		JkUtilsFile.tryDeleteDir(dir);
	}
	
	private static List<String> pythonExamples = Arrays.asList(
		"1CC8.python", "gpu"
	);
	
	private void pythonDist(String name, File dirOut, File dirDocHtml, File dirNatives) {
		
		File dirPackage = new File(dirOut, "python-package");
		dirPackage.mkdirs();
		File dirDist = new File(dirOut, "python");
		dirDist.mkdirs();
		
		// copy the python package
		File dirOsprey = new File(dirPackage, "osprey");
		JkFileTree.of(baseDir().file("python/osprey"))
			.exclude("**/*.pyc")
			.copyTo(dirOsprey);
		
		// copy install script
		JkUtilsFile.copyFileToDir(baseDir().file("python/setup.py"), dirPackage);
		
		// can't exclude whole folders for some reason, so just delete it post-hoc
		JkUtilsFile.tryDeleteDir(new File(dirOsprey, "__pycache__"));
		
		// update the osprey dev flag
		File initFile = new File(dirOsprey, "__init__.py");
		String initFileContent = JkUtilsFile.read(initFile);
		initFileContent = initFileContent.replaceFirst("_IS_DEV = True\\n", "_IS_DEV = False\n");
		JkUtilsFile.writeString(initFile, initFileContent, false);
		
		// update the setup.py rootDir
		File setupFile = new File(dirPackage, "setup.py");
		String setupFileContent = JkUtilsFile.read(setupFile);
		setupFileContent = setupFileContent.replaceFirst("rootDir = '../'\\n", "rootDir = '../../../'\n");
		JkUtilsFile.writeString(setupFile, setupFileContent, false);
		
		// copy the osprey jar
		JkUtilsFile.copyFileToDir(packer().fatJarFile(), dirOsprey);
		
		// copy the natives
		JkUtilsFile.copyDirContent(dirNatives, new File(dirOsprey, "natives"), false);
		
		// copy text files
		JkUtilsFile.copyFileToDir(baseDir().file("LICENSE.txt"), dirOsprey);
		JkUtilsFile.copyFileToDir(baseDir().file("README.rst"), dirOsprey);
		JkUtilsFile.copyFileToDir(baseDir().file("CONTRIBUTING.rst"), dirOsprey);
		
		// build the python package
		JkProcess.of("python", "setup.py", "bdist_wheel", "--universal")
			.withWorkingDir(dirPackage)
			.failOnError(true)
			.runSync();
		File pack = new File(dirPackage, "dist").listFiles()[0];
		
		// copy the package to the dist 
		JkUtilsFile.copyFileToDir(pack, dirDist);
		
		// copy docs if needed
		if (dirDocHtml != null) {
			JkUtilsFile.copyDirContent(dirDocHtml, new File(dirDist, "doc"), false);
		}
		
		// copy text files
		JkUtilsFile.copyFileToDir(baseDir().file("LICENSE.txt"), dirDist);
		JkUtilsFile.copyFileToDir(baseDir().file("README.rst"), dirDist);
		JkUtilsFile.copyFileToDir(baseDir().file("CONTRIBUTING.rst"), dirDist);
		
		// create the distribution zip
		JkFileTree.of(dirDist).zip().to(new File(dirOut, name + ".python.zip"));
		
		// cleanup temp files
		JkUtilsFile.tryDeleteDir(dirPackage);
		JkUtilsFile.tryDeleteDir(dirDist);
	}
	
	private void copyExamples(List<String> names, File dir) {
		File srcDir = baseDir().file("examples");
		for (String name : names) {
			JkUtilsFile.copyDirContent(
				new File(srcDir, name),
				new File(dir, name),
				false
			);
		}
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
