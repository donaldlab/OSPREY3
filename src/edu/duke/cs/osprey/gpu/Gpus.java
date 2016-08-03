package edu.duke.cs.osprey.gpu;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLProgram;

public class Gpus {
	
	public static boolean useProfiling;
	
	private static Gpus instance;
	
	static {
		useProfiling = false;
		instance = null;
	}
	
	public static Gpus get() {
		if (instance == null) {
			instance = new Gpus(useProfiling);
		}
		return instance;
	}
	
	private CLContext context;
	private List<Gpu> gpus;
	private Map<String,CLProgram> programs;
	
	private Gpus(boolean useProfiling) {
		
		assert (context == null);
		context = CLContext.create();
		
		// get the gpus that support doubles
		gpus = new ArrayList<>();
		for (CLDevice device : context.getDevices()) {
			if (device.getType() == CLDevice.Type.GPU) {
				Gpu gpu = new Gpu(device, useProfiling);
				if (gpu.supportsDoubles()) {
					gpus.add(gpu);
				}
			}
		}
		
		// sort gpus by flops
		Collections.sort(gpus, new Comparator<Gpu>() {
			
			@Override
			public int compare(Gpu a, Gpu b) {
				// NOTE: compare in reverse to get descending sort
				return Long.compare(getSpeed(b), getSpeed(a));
			}
			
			private long getSpeed(Gpu gpu) {
				return gpu.getDevice().getMaxComputeUnits() * gpu.getDevice().getMaxClockFrequency();
			}
		});
		
		// other init
		programs = new HashMap<>();
	}
	
	public CLContext getContext() {
		return context;
	}
	
	public List<Gpu> getGpus() {
		return Collections.unmodifiableList(gpus);
	}
	
	public Gpu getBestGpu() {
		if (gpus.isEmpty()) {
			return null;
		}
		return gpus.get(0);
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
			
			CLProgram program = context.createProgram(in).build();
			
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
