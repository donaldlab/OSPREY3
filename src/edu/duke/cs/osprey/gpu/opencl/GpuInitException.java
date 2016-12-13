package edu.duke.cs.osprey.gpu.opencl;

public class GpuInitException extends RuntimeException {

	private static final long serialVersionUID = -7618850500199218177L;
	
	public GpuInitException(String msg, Throwable cause) {
		super(msg, cause);
	}
}
