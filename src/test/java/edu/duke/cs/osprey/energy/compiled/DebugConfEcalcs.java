package edu.duke.cs.osprey.energy.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;
import jdk.incubator.foreign.MemoryHandles;
import jdk.incubator.foreign.MemorySegment;
import org.hamcrest.MatcherAssert;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.lang.invoke.VarHandle;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Supplier;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * A simple debugging harness to replay issues with conf energy calculators
 */
public class DebugConfEcalcs {

	public static void main(String[] args) {

		log("reading conf space ...");
		var file = new File("/path/to/confspace.ccsx");
		var confSpace = ConfSpace.fromBytes(FileTools.readFileBytes(file));
		log("done!");

		// get the wild type conformation, if possible
		var conf = Arrays.stream(confSpace.positions)
			.mapToInt(pos ->
				Arrays.stream(pos.confs)
					.filter(c -> c.id.startsWith("wt-"))
					.mapToInt(c -> c.index)
					.findFirst()
					.orElse(0) // otherwise, pick an arbitrary conformation
			)
			.toArray();

		var precision = Structs.Precision.Float64;
		var ecalcFactories = new ArrayList<Supplier<ConfEnergyCalculator>>();
		ecalcFactories.add(() -> new CPUConfEnergyCalculator(confSpace));
		ecalcFactories.add(() -> new NativeConfEnergyCalculator(confSpace, precision));
		if (CudaConfEnergyCalculator.isSupported()) {
			ecalcFactories.add(() -> new CudaConfEnergyCalculator(confSpace, precision, Parallelism.make(1, 1)));
		}

		// just compute the static-static interactions, since they're the simplest
		var inters = PosInterDist.staticStatic();

		for (var ecalcFactory : ecalcFactories) {
			log("making energy calculator ...");
			var ecalc = ecalcFactory.get();
			log("Ecalc: %s", ecalc.getClass().getSimpleName());
			log("\tcalc:        %f", ecalc.calc(conf, inters).energy);
			log("\tcalc energy: %f", ecalc.calcEnergy(conf, inters));
			log("\tmin:         %f", ecalc.minimize(conf, inters).energy);
			log("\tmin energy:  %f", ecalc.minimizeEnergy(conf, inters));
		}

		log("done");
	}
}
