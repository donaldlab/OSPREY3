package edu.duke.cs.osprey.gpu.opencl;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLProgram;

import edu.duke.cs.osprey.tools.FileTools;

public class Gpu {
	
	// NOTE: the opencl driver caches compiled kernels
	// so you won't get a log if the driver pulls the binary from the cache
	// make some trivial change to the source to force a recompile
	private static final boolean DumpCompilerLog = false;
	
	private CLDevice device;
	private Map<String,CLProgram> programs;
	
	public Gpu(CLDevice device) {
		this.device = device;
		programs = new HashMap<>();
	}
	
	public CLDevice getDevice() {
		return device;
	}
	
	public boolean supportsDoubles() {
		Set<String> extensions = device.getExtensions();
		return extensions.contains("cl_khr_fp64") || extensions.contains("cl_amd_fp64");
	}
	
	@Override
	public String toString() {
		return device.getName();
	}
	
	public CLProgram getProgram(String filename)
	throws IOException {
		CLProgram program = programs.get(filename);
		if (program != null) {
			return program;
		}
		return compileProgram(filename);
	}
	
	public CLProgram compileProgram(String filename)
	throws IOException {
		
		try (InputStream in = FileTools.openResource("/gpuKernels/opencl/" + filename, getClass())) {
			
			CLProgram program;
			if (DumpCompilerLog) {
				program = device.getContext().createProgram(in).build("-cl-nv-verbose");
				System.out.println(program.getBuildLog());
			} else {
				program = device.getContext().createProgram(in).build();
			}
			
			programs.put(filename, program);
			return program;
			
		} catch (IOException ex) {
			throw new IOException(String.format("error compiling gpu program file: %s", filename), ex);
		}
	}
}
