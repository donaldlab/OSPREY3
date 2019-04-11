package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.Benchmark;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Random;

import static edu.duke.cs.osprey.tools.Log.log;

public class MathLab {

	private static double randDouble(Random r) {

		double d = Double.NaN;

		while (Double.isNaN(d)) {

			// Random.nextDouble() only returns numbers in [0,1]
			int lo = r.nextInt();
			int hi = r.nextInt();
			long bits = ((long)hi << 32) | (long)lo;

			d = Double.longBitsToDouble(bits);
		}

		return d;
	}

	public static void main(String[] args) {

		// get some random doubles
		Random rand = new Random(12345);

		int nums = 1024;
		double[] as = new double[nums];
		double[] bs = new double[nums];
		for (int i=0; i<nums; i++) {
			as[i] = randDouble(rand);
			bs[i] = randDouble(rand);
		}

		double[] lnas = new double[nums];
		double[] lnbs = new double[nums];
		for (int i=0; i<nums; i++) {
			lnas[i] = Math.log(as[i]);
			lnbs[i] = Math.log(bs[i]);
		}

		int pairs = 1024;
		int[] indices = new int[pairs*2];
		for (int i=0; i<pairs*2; i++) {
			indices[i] = rand.nextInt(nums);
		}

		double[] cs = new double[pairs];

		Benchmark dplus = new Benchmark(1, 10000, 5000000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i = 0; i < n; ) {
				cs[o++] = as[indices[i++]] + bs[indices[i++]];
			}
		});
		log("double +: %s", dplus);

		Benchmark dmult = new Benchmark(1, 10000, 5000000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				cs[o++] = as[indices[i++]] * bs[indices[i++]];
			}
		});
		log("double *: %s", dmult);


		class BoxDouble {

			final double val;

			public BoxDouble(double val) {
				this.val = val;
			}

			public BoxDouble add(BoxDouble other) {
				return new BoxDouble(this.val + other.val);
			}

			public BoxDouble mult(BoxDouble other) {
				return new BoxDouble(this.val * other.val);
			}
		}

		BoxDouble[] boxas = new BoxDouble[nums];
		BoxDouble[] boxbs = new BoxDouble[nums];
		for (int i=0; i<nums; i++) {
			boxas[i] = new BoxDouble(as[i]);
			boxbs[i] = new BoxDouble(bs[i]);
		}

		BoxDouble[] boxcs = new BoxDouble[pairs];

		Benchmark boxdplus = new Benchmark(1, 10000, 1000000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i = 0; i < n; ) {
				boxcs[o++] = boxas[indices[i++]].add(boxbs[indices[i++]]);
			}
		});
		log("BoxDbl +: %s", boxdplus);

		Benchmark boxdmult = new Benchmark(1, 10000, 1000000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				boxcs[o++] = boxas[indices[i++]].mult(boxbs[indices[i++]]);
			}
		});
		log("BoxDbl *: %s", boxdmult);

		BigExp[] beas = new BigExp[nums];
		BigExp[] bebs = new BigExp[nums];
		BigExp[] becs = new BigExp[pairs];
		for (int i=0; i<nums; i++) {
			beas[i] = new BigExp(as[i]);
			bebs[i] = new BigExp(bs[i]);
			becs[i] = new BigExp(0);
		}

		Benchmark bemult = new Benchmark(1, 10000, 1000000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				BigExp out = becs[o++];
				out.set(beas[indices[i++]]);
				out.mult(bebs[indices[i++]]);
			}
		});
		log("BigExp *: %s", bemult);

		Benchmark bemax = new Benchmark(1, 1000, 100000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				BigExp out = becs[o++];
				out.set(beas[indices[i++]]);
				out.max(bebs[indices[i++]]);
			}
		});
		log("BigExp max: %s", bemax);


		BigDecimal[] bigas = new BigDecimal[nums];
		BigDecimal[] bigbs = new BigDecimal[nums];
		for (int i=0; i<nums; i++) {
			bigas[i] = MathTools.biggen(as[i]);
			bigbs[i] = MathTools.biggen(bs[i]);
		}

		BigDecimal[] bigcs = new BigDecimal[pairs];


		MathContext mc32 = new MathContext(32, RoundingMode.HALF_UP);
		log("fixed num bytes: %d", new BigDecimalIO.Fixed(mc32).numBytes);
		BoltzmannCalculator bcalc32 = new BoltzmannCalculator(mc32);

		Benchmark bdplus = new Benchmark(1, 100, 50000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				bigcs[o++] = bigas[indices[i++]].add(bigbs[indices[i++]], mc32);
			}
		});
		log("BigD32 +: %s", bdplus);

		Benchmark bdmult = new Benchmark(1, 100, 50000, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				bigcs[o++] = bigas[indices[i++]].multiply(bigbs[indices[i++]], mc32);
			}
		});
		log("BigD32 *: %s", bdmult);

		Benchmark bdexp = new Benchmark(1, 10, 500, () -> {
			int n = indices.length;
			int o = 0;
			for (int i=0; i<n; ) {
				bigcs[o++] = bcalc32.exp(lnas[indices[i++]] + lnbs[indices[i++]]);
			}
		});
		log("BigD32 +exp: %s", bdexp);

		// compute ops for a mix of stuff
		//log("double mult, exp: %.2f", 1.0/(1000/dmult.opsPerSecond + 1/bdexp.opsPerSecond));
		//log("BigDecimal mult: %.2f", 1.0/(1000/bdmult.opsPerSecond));
	}
}
