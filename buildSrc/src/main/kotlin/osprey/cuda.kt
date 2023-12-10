package osprey

import org.gradle.api.Project
import org.gradle.api.tasks.Exec
import org.gradle.api.tasks.JavaExec
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue


fun Project.makeCudaTasks() {

	val compileCuda_residueForcefield by tasks.creating(Exec::class) {
		nvcc(this, "residueForcefield")
	}

	val compileCuda_residueCcd by tasks.creating(Exec::class) {
		nvcc(this, "residueCcd", maxRegisters=64)
	}

	val benchmarkCuda by tasks.creating(JavaExec::class) {
		group = "Execution"
		description = "Run Benchmark Energies"
		classpath = sourceSets.test.runtimeClasspath
		mainClass.set("edu.duke.cs.osprey.BenchmarkEnergies")
		jvmArgs = listOf("--enable-preview")
	}

	@Suppress("UNUSED_VARIABLE")
	val compileCuda by tasks.creating {
		description = "Compile cuda kernels"
		dependsOn(
			compileCuda_residueForcefield,
			compileCuda_residueCcd
		)
	}
}

fun Project.nvcc(exec: Exec, kernelName: String, maxRegisters: Int? = null, profile: Boolean = false) {

	val args = mutableListOf("nvcc")

	if (profile) {
		// if profiling, compile for one arch with profiling/debug info
		// NOTE: change this to your GPU's arch
		args.addAll(listOf("-cubin", "-gencode=arch=compute_61,code=sm_61", "-lineinfo", "--ptxas-options=-v"))
	} else {
		// otherwise, compile for all archs
		// see Maxwell compatibility guide:
		// http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
		args.addAll(listOf("-fatbin",
			"-gencode=arch=compute_20,code=sm_20",
			"-gencode=arch=compute_30,code=sm_30",
			"-gencode=arch=compute_35,code=sm_35",
			"-gencode=arch=compute_50,code=sm_50",
			"-gencode=arch=compute_52,code=sm_52",
			"-gencode=arch=compute_60,code=sm_60",
			"-gencode=arch=compute_61,code=sm_61",
			"-gencode=arch=compute_62,code=sm_62",
			"-gencode=arch=compute_62,code=compute_62"
		))
	}

	if (maxRegisters != null) {
		args.addAll(listOf("-maxrregcount", "$maxRegisters"))
	}

	args.addAll(listOf("$kernelName.cu", "-o", "$kernelName.bin"))

	exec.workingDir = file("src/main/resources/gpuKernels/cuda")
	exec.commandLine(args)
}
