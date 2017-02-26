package edu.duke.cs.osprey;

import java.io.File;

import org.jerkar.api.depmanagement.JkDependencies;
import org.jerkar.api.depmanagement.JkDependency;
import org.jerkar.api.depmanagement.JkDependency.JkFileDependency;
import org.jerkar.api.depmanagement.JkModuleId;
import org.jerkar.api.depmanagement.JkScope;
import org.jerkar.api.depmanagement.JkVersion;
import org.jerkar.api.file.JkFileTree;
import org.jerkar.api.file.JkFileTreeSet;
import org.jerkar.api.file.JkPathFilter;
import org.jerkar.api.java.JkJavaCompiler;
import org.jerkar.api.system.JkProcess;
import org.jerkar.api.utils.JkUtilsFile;
import org.jerkar.api.utils.JkUtilsZip;
import org.jerkar.tool.builtins.javabuild.JkJavaBuild;
import org.jerkar.tool.builtins.javabuild.JkJavaPacker;

public class Build extends JkJavaBuild {
	
	private static final JkScope NATIVES = JkScope.of("natives");
	
	public Build() {
		
		// don't run unit tests in the build system
		// it takes forever
		tests.skip = true;
		
		manifest.mainClass = "edu.duke.cs.osprey.control.Main";
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
			
			// for dependencies with natives, use local jars in lib/
			.on(new File("lib/jocl-2.3.2.jar"))
			.on(new File("lib/jocl-2.3.2-natives-linux-amd64.jar")).scope(NATIVES)
			.on(new File("lib/jocl-2.3.2-natives-macosx-universal.jar")).scope(NATIVES)
			.on(new File("lib/jocl-2.3.2-natives-windows-amd64.jar")).scope(NATIVES)
			.on(new File("lib/gluegen-rt-2.3.2.jar"))
			.on(new File("lib/gluegen-rt-2.3.2-natives-linux-amd64.jar")).scope(NATIVES)
			.on(new File("lib/gluegen-rt-2.3.2-natives-macosx-universal.jar")).scope(NATIVES)
			.on(new File("lib/gluegen-rt-2.3.2-natives-windows-amd64.jar")).scope(NATIVES)
			.on(new File("lib/jcuda-0.8.0RC.jar"))
			.on(new File("lib/jcuda-natives-0.8.0RC-linux-x86_64.jar")).scope(NATIVES)
			.on(new File("lib/jcuda-natives-0.8.0RC-windows-x86_64.jar")).scope(NATIVES)
			// no jcuda binaries released for mac =(
			// but we could try compiling from source if we wanted
			
			.build();
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
		
		// copy the docs
		File dirDoc = baseDir().file("python/doc");
		JkProcess.of("make", "clean")
			.withWorkingDir(dirDoc)
			.runSync();
		JkProcess.of("make", "html")
			.withWorkingDir(dirDoc)
			.runSync();
		JkUtilsFile.copyDirContent(
			baseDir().file("python/doc/_build/html"),
			ouputDir("dist/doc"),
			false
		);
		
		// extract the natives to a temp folder
		File dirTemp = ouputDir("temp");
		for (JkDependency dep : dependencies().dependenciesDeclaredWith(NATIVES)) {
			JkFileDependency fileDep = (JkFileDependency)dep;
			File depFile = fileDep.files().iterator().next();
			JkUtilsZip.unzip(depFile, dirTemp);
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
		
		// copy text files
		JkUtilsFile.copyFileToDir(baseDir().file("LICENSE.txt"), dirDist);
		JkUtilsFile.copyFileToDir(baseDir().file("README.rst"), dirDist);
		// TODO: do we need a different README (in RST format) for python packages?
		
		// create the distribution zip
		JkFileTree.of(dirDist).zip().to(ouputDir(packer().fatJarFile().getName().replaceAll(".jar", ".zip")));
	}
}
