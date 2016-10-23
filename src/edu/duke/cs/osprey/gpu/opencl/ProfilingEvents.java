package edu.duke.cs.osprey.gpu.opencl;

import java.lang.reflect.Field;
import java.util.concurrent.TimeUnit;

import com.jogamp.common.nio.PointerBuffer;
import com.jogamp.opencl.CLEvent;
import com.jogamp.opencl.CLEvent.CommandType;
import com.jogamp.opencl.CLEvent.ProfilingCommand;
import com.jogamp.opencl.CLEventList;
import com.jogamp.opencl.CLException;

import edu.duke.cs.osprey.tools.TimeFormatter;

public class ProfilingEvents {
	
	private static Field idsField;
	
	static {
		try {
			idsField = CLEventList.class.getDeclaredField("IDs");
			idsField.setAccessible(true);
		} catch (Exception ex) {
			throw new Error("can't hack the events ids, did the JOCL library change?", ex);
		}
	}
	
	private CLEventList events;
	private PointerBuffer ids;
	
	public ProfilingEvents(int numEvents) {
		
		events = new CLEventList(numEvents);
		
		// HACKHACK: get the events ids from the private field for low-level commands
		if (events != null) {
			try {
				ids = (PointerBuffer)idsField.get(events);
			} catch (Exception ex) {
				throw new Error("can't hack the events ids, did the JOCL library change?", ex);
			}
		}
	}
	
	public CLEventList getCLEvents() {
		return events;
	}
	
	public PointerBuffer getCLIds() {
		return ids;
	}
	
	public String makeReport() {
		
		StringBuilder buf = new StringBuilder();
		
		buf.append("GPU Profile: ");
		buf.append(events.size());
		buf.append(" events");
		
		long sumNs = 0;
		for (CLEvent event : events) {
			
			CommandType type = event.getType();
			
			try {
				
				long startNs = event.getProfilingInfo(ProfilingCommand.START);
				long endNs = event.getProfilingInfo(ProfilingCommand.END);
				long diffNs = endNs - startNs;
				sumNs += diffNs;
				
				buf.append(String.format("\n\t%16s %s", type.name(), TimeFormatter.format(diffNs, TimeUnit.MICROSECONDS)));
				
			} catch (CLException.CLProfilingInfoNotAvailableException ex) {
				
				buf.append(String.format("\n\t%16s (unknown timing)", type.name()));
			}
		}
		
		buf.append(String.format("\n\tTotal: %s", TimeFormatter.format(sumNs, TimeUnit.MICROSECONDS)));
		
		return buf.toString();
	}
	
	public long getTotalNs() {
		long sumNs = 0;
		for (CLEvent event : events) {
			try {
				
				long startNs = event.getProfilingInfo(ProfilingCommand.START);
				long endNs = event.getProfilingInfo(ProfilingCommand.END);
				long diffNs = endNs - startNs;
				sumNs += diffNs;
				
			} catch (CLException.CLProfilingInfoNotAvailableException ex) {
				System.err.println("can't get timing for profiling event: " + event);
			}
		}
		return sumNs;
	}
	
	public void cleanup() {
		events.release();
	}
}
