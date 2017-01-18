package edu.duke.cs.osprey;

import java.io.File;

import org.jerkar.api.depmanagement.JkDependencies;
import org.jerkar.api.depmanagement.JkModuleId;
import org.jerkar.api.depmanagement.JkVersion;
import org.jerkar.api.file.JkFileTree;
import org.jerkar.api.file.JkFileTreeSet;
import org.jerkar.api.file.JkPathFilter;
import org.jerkar.api.java.JkJavaCompiler;
import org.jerkar.api.utils.JkUtilsFile;
import org.jerkar.tool.builtins.javabuild.JkJavaBuild;
import org.jerkar.tool.builtins.javabuild.JkJavaPacker;

public class Build extends JkJavaBuild {
	
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
		return JkVersion.name("3.0-alpha1");
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
			.on(new File("lib/jocl-2.3.2-natives-linux-amd64.jar"))
			.on(new File("lib/jocl-2.3.2-natives-macosx-universal.jar"))
			.on(new File("lib/jocl-2.3.2-natives-windows-amd64.jar"))
			.on(new File("lib/gluegen-rt-2.3.2.jar"))
			.on(new File("lib/gluegen-rt-2.3.2-natives-linux-amd64.jar"))
			.on(new File("lib/gluegen-rt-2.3.2-natives-macosx-universal.jar"))
			.on(new File("lib/gluegen-rt-2.3.2-natives-windows-amd64.jar"))
			.on(new File("lib/jcuda-0.8.0RC.jar"))
			.on(new File("lib/jcuda-natives-0.8.0RC-linux-x86_64.jar"))
			.on(new File("lib/jcuda-natives-0.8.0RC-windows-x86_64.jar"))
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
		
		// copy the python package
		JkFileTree.of(baseDir().file("python"))
			.exclude("**/*.pyc")
			.copyTo(dirDist);
		
		// can't exclude whole folders for some reason, so just delete it post-hoc
		JkUtilsFile.deleteDir(new File(dirDist, "osprey.egg-info"));
		
		// copy the osprey jar
		File dirOsprey = new File(dirDist, "osprey");
		File jarFile = packer().fatJarFile();
		JkUtilsFile.copyFileToDir(jarFile, dirOsprey);
		
		// update the osprey jar path
		File initFile = new File(dirOsprey, "__init__.py");
		String initFileContent = JkUtilsFile.read(initFile);
		initFileContent = initFileContent.replace("_ospreyPath = '../../bin'", "_ospreyPath = '" + jarFile.getName() + "'");
		JkUtilsFile.writeString(initFile, initFileContent, false);
		
		// copy text files
		JkUtilsFile.copyFileToDir(baseDir().file("LICENSE.txt"), dirDist);
		JkUtilsFile.copyFileToDir(baseDir().file("README.md"), dirDist);
		// TODO: do we need a different README (in RST format) for python packages?
	}
}
