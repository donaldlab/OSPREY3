package edu.duke.cs.osprey.gpu.cuda;

import jcuda.driver.CUdevice;
import jcuda.driver.CUdevice_attribute;
import jcuda.driver.JCudaDriver;

public class Diagnostics {
	
	public static void main(String[] args) {
	
		try {
			// get the osprey-usable GPUs first (inits the CUDA libraries)
			Gpus gpus = Gpus.get();
			
			// how many gpus are there?
			int[] ints = new int[1];
			JCudaDriver.cuDeviceGetCount(ints);
			int count = ints[0];
			
			// show devices
			System.out.println("All CUDA devices:");
			for (int i=0; i<count; i++) {
				
				CUdevice device = new CUdevice();
				JCudaDriver.cuDeviceGet(device, i);
				Gpu gpu = new Gpu(device);
				
				System.out.println("\n\t" + gpu.getName());
				System.out.println(String.format("\t\t%-30s %d MiB", "memory:", gpu.getTotalMemory()/1024/1024));
				System.out.println(String.format("\t\t%-30s %d.%d", "compute version:", gpu.getComputeVersion()[0], gpu.getComputeVersion()[1]));
				
				// get some device attributes
				System.out.println(String.format("\t\t%-30s %s", "integrated or discrete:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_INTEGRATED) == 1 ? "integrate" : "discrete"));
				System.out.println(String.format("\t\t%-30s %d MHz", "clock rate:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_CLOCK_RATE)/1000));
				System.out.println(String.format("\t\t%-30s %d", "multiprocessors:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)));
				System.out.println(String.format("\t\t%-30s %d", "threads per multiprocessor:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR)));
				System.out.println(String.format("\t\t%-30s %d", "warp size:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_WARP_SIZE)));
				System.out.println(String.format("\t\t%-30s %d", "max block threads:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X)));
				System.out.println(String.format("\t\t%-30s %d", "max grid blocks:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X)));
				System.out.println(String.format("\t\t%-30s %d", "registers per block (32-bit):", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK)));
				System.out.println(String.format("\t\t%-30s %d KiB", "memory per block:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK)/1024));
				System.out.println(String.format("\t\t%-30s %d KiB", "constant memory:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY)/1024));
				System.out.println(String.format("\t\t%-30s %d", "Async engine count:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_ASYNC_ENGINE_COUNT)));
				
				System.out.println();
				
				System.out.println(String.format("\t\t%-30s %b", "64-bit double support:", gpu.supportsDoubles()));
				System.out.println(String.format("\t\t%-30s %b", "dynamic parallelism support:", gpu.supportsDynamicParallelism()));
				System.out.println(String.format("\t\t%-30s %b", "concurrent kernels:", gpu.getAttribute(CUdevice_attribute.CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS)));
			}
			
			// show the usable GPUs
			System.out.println("\nCUDA GPUs usable by Osprey:");
			for (Gpu gpu : gpus.getGpus()) {
				System.out.println(String.format("%30s:     memory free: %.1f%%", gpu, (100.0f*gpu.getFreeMemory()/gpu.getTotalMemory())));
			}
			
		} catch (UnsatisfiedLinkError ex) {
			System.out.println("CUDA does not appear to be installed. No GPUs to show.");
			return;
		}
	}
}
