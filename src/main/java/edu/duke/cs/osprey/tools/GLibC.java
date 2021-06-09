package edu.duke.cs.osprey.tools;

import com.sun.jna.Native;
import com.sun.jna.NativeLong;
import com.sun.jna.Structure;

import java.util.EnumMap;


/**
 * Access GLibC functions from Java
 */
public class GLibC {

	private static class NativeLib {

		static {
			Native.register("c");
		}

		/*
		struct rlimit {
			__kernel_ulong_t	rlim_cur;
			__kernel_ulong_t	rlim_max;
		};
		*/
		@Structure.FieldOrder({"rlim_curr", "rlim_max"})
		public static class RLimit extends Structure {
			public NativeLong rlim_curr;
			public NativeLong rlim_max;
		}

		// https://man7.org/linux/man-pages/man2/getrlimit.2.html
		public static native int getrlimit(int resource, RLimit rlim);
	}

	/*
	enum __rlimit_resource {
		RLIMIT_CPU = 0,
		RLIMIT_FSIZE = 1,
		RLIMIT_DATA = 2,
		RLIMIT_STACK = 3,
		RLIMIT_CORE = 4,
		RLIMIT_RSS = 5,
		RLIMIT_NPROC = 6,
		RLIMIT_NOFILE = 7,
		RLIMIT_MEMLOCK = 8,
		RLIMIT_AS = 9,
		RLIMIT_LOCKS = 10,
		RLIMIT_SIGPENDING = 11,
		RLIMIT_MSGQUEUE = 12,
		RLIMIT_NICE = 13,
		RLIMIT_RTPRIO = 14,
		RLIMIT_RTTIME = 15,
	};
	*/
	enum RLimitResource {
		CPU,
		FSIZE,
		DATA,
		STACK,
		CORE,
		RSS,
		NPROC,
		NOFILE,
		MEMLOCK,
		AS,
		LOCKS,
		SIGPENDING,
		MSGQUEUE,
		NICE,
		RTPRIO,
		RTTIME
	}

	public static class RLimit {

		long curr;
		long max;

		private RLimit(NativeLib.RLimit n) {
			curr = n.rlim_curr.longValue();
			max = n.rlim_max.longValue();
		}
	}

	public static RLimit getrlimit(RLimitResource resource) {
		var rlimit = new NativeLib.RLimit();
		NativeLib.getrlimit(resource.ordinal(), rlimit);
		return new RLimit(rlimit);
	}

	public static EnumMap<RLimitResource,RLimit> getrlimits() {
		var map = new EnumMap<RLimitResource,RLimit>(RLimitResource.class);
		for (var resource : RLimitResource.values()) {
			map.put(resource, getrlimit(resource));
		}
		return map;
	}
}
