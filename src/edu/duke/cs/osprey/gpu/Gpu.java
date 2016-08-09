package edu.duke.cs.osprey.gpu;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLProgram;

public class Gpu {
	
	private CLDevice device;
	private Map<String,CLProgram> programs;
	
	public Gpu(CLDevice device) {
		this.device = device;
		programs = new HashMap<>();
	}
	
	public CLDevice getDevice() {
		return device;
	}
	
	public GpuQueue makeQueue() {
		return makeQueue(false);
	}
	
	public GpuQueue makeQueue(boolean useProfiling) {
		CLCommandQueue queue;
		if (useProfiling) {
			queue = device.createCommandQueue(CLCommandQueue.Mode.PROFILING_MODE);
		} else {
			queue = device.createCommandQueue();
		}
		return new GpuQueue(this, queue);
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
		URL url = getClass().getResource("kernelSource/" + filename);
		if (url == null) {
			throw new IOException(String.format("can't find gpu program file: %s\nlooked in folder: %s", filename, getClass().getResource("kernels").toString()));
		}
		try (InputStream in = url.openStream()) {
			
			CLProgram program = device.getContext().createProgram(in).build();
			
			// DEBUG: see compiler output, if there is any
			//CLProgram program = context.createProgram(in).build("-cl-nv-verbose");
			//System.out.println(program.getBuildLog());
			
			programs.put(filename, program);
			return program;
			
		} catch (IOException ex) {
			throw new IOException(String.format("error compiling gpu program file: %s", filename), ex);
		}
	}
}
