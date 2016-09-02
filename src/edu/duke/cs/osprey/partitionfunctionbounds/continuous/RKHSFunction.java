package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.function.ToDoubleFunction;
import Jama.*;

public class RKHSFunction {

    private Kernel k;
    private double[] domainLB;
    private double[] domainUB;

    private FeatureMap[] featureMaps;
    private double[] coeffs;

    /**
     * Constructor - takes as parameters feature maps and linear coefficients
     *
     * @param featureMaps
     * @param coeffs
     */
    public RKHSFunction(
	    FeatureMap[] featureMaps,
	    double[] coeffs) {

	// check featuremap list for kernel consistency
	for (FeatureMap fmap : featureMaps) {
	    Kernel fK = fmap.getKernel();
	    if (!k.equals(fK)) {
		throw new RuntimeException("RKHS Function kernels don't match.");
	    }
	}

	this.featureMaps = featureMaps;
	this.coeffs = coeffs;

	// set up kernel, check for consistency 
	k = featureMaps[0].getKernel();
	domainLB = new Matrix(k.bounds).transpose().getArray()[0];
	domainUB = new Matrix(k.bounds).transpose().getArray()[1];
    }

    /**
     * Constructor - given a kernel, domain, and function -- fits the RKHS function
     *
     * @param k
     * @param domainLB
     * @param domainUB
     * @param f
     */
    public RKHSFunction(Kernel k, double[] domainLB, double[] domainUB, ToDoubleFunction<double[]> f) {
	// gridSample from dimensions
	// construct feature maps
	// fit coefficients
	// test fit --> increase samples if necessary
	gridSampleFit(k, f, domainLB, domainUB);
    }

    /**
     * Fits coefficients for a linear combination of feature maps to approximate a given function
     *
     * @param fMaps
     * @param func
     * @param k
     * @return
     */
    public double[] fitCoeffs(FeatureMap[] fMaps, ToDoubleFunction<double[]> func, Kernel k) {
	Matrix gMat = new Matrix(fMaps.length, fMaps.length);
	for (int i = 0; i < fMaps.length; i++) {
	    for (int j = 0; j < fMaps.length; j++) {
		gMat.set(i, j, k.eval(fMaps[i].getLoc(), fMaps[j].getLoc()));
	    }
	}

	Matrix fMat = new Matrix(fMaps.length, 1);
	for (int i = 0; i < fMaps.length; i++) {
	    fMat.set(i, 0, func.applyAsDouble(fMaps[i].getLoc()));
	}

	Matrix cMat = gMat.solve(fMat);
	return cMat.transpose().getArray()[0];
    }

    /**
     * Evaluates the RKHS function at a given point
     *
     * @param point
     * @return
     */
    public double eval(double[] point) {
	double result = 0.0;
	for (int i = 0; i < featureMaps.length; i++) {
	    result += featureMaps[i].eval(point) * coeffs[i];
	}
	return result;
    }

    /**
     * Grid samples and fits the RKHS function to the DoubleFunction provided Increases number of
     * samples as necessary
     *
     * @param k
     * @param f
     * @param domainLB
     * @param domainUB
     */
    private void gridSampleFit(Kernel k, ToDoubleFunction<double[]> f, double[] domainLB, double[] domainUB) {
	System.out.println("Beginning fit...");

	int numSamples = 90;
	double rmsd = Double.POSITIVE_INFINITY;

	FeatureMap[] fMaps = new FeatureMap[numSamples];
	double[] fitCoeffs = new double[numSamples];

	while (rmsd > 0.05) {
	    numSamples += 10;
	    double[][] samples = gridSample(numSamples, domainLB, domainUB);
	    fMaps = new FeatureMap[numSamples];
	    for (int i = 0; i < numSamples; i++) {
		fMaps[i] = new FeatureMap(k, samples[i]);
	    }
	    fitCoeffs = fitCoeffs(fMaps, f, k);
	    rmsd = computeFitRMSD(fMaps, fitCoeffs, domainLB, domainUB, f);
	    System.out.println("With " + numSamples + " functions RMSD is: " + rmsd);
	}

	System.out.println("Required " + numSamples + " functions, total RMSD: " + rmsd + ".");

	this.featureMaps = fMaps;
	this.coeffs = fitCoeffs;
    }

    /**
     * Grid samples from n dimensions -- row is sample, column is dimension Note that numSamples
     * needs to be a cube of some integer
     *
     * @param numSamples
     * @param lb
     * @param ub
     * @return
     */
    public double[][] gridSample(int numSamples, double[] lb, double[] ub) {
	int numDimensions = lb.length;
	double[][] samples = new double[numSamples][numDimensions];
	for (int i = 0; i < numDimensions; i++) {
	    double[] dimSamples = gridSample1D(numSamples, lb[i], ub[i]);
	    for (int j = 0; j < numSamples; j++) {
		samples[j][i] = dimSamples[j];
	    }
	}
	return samples;
    }

    /**
     * Grid samples from one dimension
     *
     * @param nSamples
     * @param lowerBound
     * @param upperBound
     * @return
     */
    public double[] gridSample1D(int nSamples, double lowerBound, double upperBound) {
	double[] gridSamples = new double[nSamples];
	int ind = 0;
	double increment = ((upperBound - lowerBound) / nSamples) * 0.95;
	double i = lowerBound;
	while (ind < nSamples) {
	    gridSamples[ind] = i;
	    i += increment;
	    ind++;
	}
	return gridSamples;
    }

    /**
     * Computes the RMSD of a fitted linear combination
     *
     * @param fmaps
     * @param coeffs
     * @param lBound
     * @param uBound
     * @param func
     * @return
     */
    public double computeFitRMSD(
	    FeatureMap[] fmaps,
	    double[] coeffs,
	    double[] lBound,
	    double[] uBound,
	    ToDoubleFunction<double[]> func) {
	double rmsd = 0.0;
	int numSamples = coeffs.length;
	double[][] testSamples = gridSample(numSamples * 50, lBound, uBound);
	for (double[] sample : testSamples) {
	    double funcValue = func.applyAsDouble(sample);
	    double fitValue = evalLC(fmaps, coeffs, sample);
	    double delta = funcValue - fitValue;
	    rmsd += Math.pow(delta, 2);
	}
	return Math.sqrt(rmsd / testSamples.length);
    }

    /**
     * Evaluates a linear combination of feature maps at a given point
     *
     * @param fMaps
     * @param coeffs
     * @param x
     * @return
     */
    public double evalLC(FeatureMap[] fMaps, double[] coeffs, double[] x) {
	double result = 0.0;
	for (int i = 0; i < fMaps.length; i++) {
	    result += fMaps[i].eval(x) * coeffs[i];
	}
	return result;
    }

    /**
     * Computes the empirical mean map for this function
     * @return RKHS function representing the mean map of the function
     */
    public RKHSFunction computeMeanMap() {
	// sample a point u.a.r from n dimensions
	// binary search for that point
	// repeat k times
	int numFeatureMaps = this.featureMaps.length;

	// compute value of the integral of the mass function over the whole domain
	int numDimensions = this.domainLB.length;
	double domainVolume = 0;
	for (int i=0; i<numDimensions; i++) { 
	    domainVolume = domainVolume * (domainUB[i] - domainLB[i]);
	}
	RKHSFunction lebesgueMeasure = getLebesgueMeasure(domainLB, domainUB);
	double functionIntegral = domainVolume * this.innerProduct(lebesgueMeasure);
	
	FeatureMap[] newFeatureMaps = new FeatureMap[numFeatureMaps];
	
	for (int fMapIndex = 0; fMapIndex < numFeatureMaps; fMapIndex++) {
	    double cdfValue = 1;
	    for (int i = 0; i < numDimensions; i++) {
		cdfValue = cdfValue * ((Math.random() * (domainUB[i] - domainLB[i])) + domainLB[i]);
	    }

	    double[] cdfPoint = binarySearchCDF(cdfValue, domainLB, domainUB, functionIntegral);
	    
	    newFeatureMaps[fMapIndex] = new FeatureMap(this.k, cdfPoint);
	}

	ToDoubleFunction<double[]> thisFunction = (point)->(this.eval(point));
	double[] newCoeffs = fitCoeffs(newFeatureMaps, thisFunction, this.k);
	return new RKHSFunction(newFeatureMaps, newCoeffs);
    }

    /**
     * Binary searches the domain to find a point whose CDF matches the value provided
     * @param cdfValue
     * @param LBs
     * @param UBs
     * @param totalIntegral
     * @return 
     */
    public double[] binarySearchCDF(
	    double cdfValue,
	    double[] LBs,
	    double[] UBs,
	    double totalIntegral) {

	double domainVolume = 1;
	int numDimensions = UBs.length;
	double epsilon = 0.005;
	double[] point = new double[UBs.length];
	for (int i=0; i<numDimensions; i++) {
	    point[i] = 0.5*(UBs[i] + UBs[i]);
	}
	
	// compute CDF of the point 
	RKHSFunction subspaceMeasure = getLebesgueMeasure(domainLB, point);
	double subspaceVolume = 1;
	for (int i=0; i<numDimensions; i++) { 
	    subspaceVolume = subspaceVolume * (point[i] - domainLB[i]);
	}
	double pointCDF = (this.innerProduct(subspaceMeasure) * subspaceVolume)/totalIntegral;
	
	if (Math.abs(pointCDF - cdfValue) < epsilon) { 
	    return point;
	} else if (pointCDF < cdfValue) { 
	    // change the bounds
	    int coordinate = (int) (Math.random() * numDimensions);
	    double[] newLBs = new double[LBs.length];
	    for (int i=0; i<newLBs.length; i++) { 
		newLBs[i] = (i==coordinate) ? point[i] : LBs[i];
	    }
	    return binarySearchCDF(cdfValue, newLBs, UBs, totalIntegral);
	} else if (pointCDF > cdfValue) { 
	    // change the bounds
	    int coordinate = (int) (Math.random() * numDimensions);
	    double[] newUBs = new double[UBs.length];
	    for (int i=0; i<newUBs.length; i++) { 
		newUBs[i] = (i==coordinate) ? point[i] : UBs[i];
	    }
	    return binarySearchCDF(cdfValue, LBs, newUBs, totalIntegral);
	} else {
	    throw new RuntimeException("Point CDF is neither equal to, greater than, or less than "
		    + "the required CDF value. What on earth did you do?");
	}
    }

    /**
     * Computes the inner product of this RKHSFunction with another RKHSFunction Assumes that the
     * kernel for both functions is the same
     *
     * @param f
     * @return
     */
    public double innerProduct(RKHSFunction f) {

	// check if both kernels are the same
	if (!this.k.equals(f.k)) {
	    throw new RuntimeException("Inner product between RKHS functions with different kernels.");
	}

	double innerProduct = 0;
	for (int i = 0; i < f.featureMaps.length; i++) {
	    for (int j = 0; j < this.featureMaps.length; j++) {
		innerProduct
			+= f.coeffs[i] * this.coeffs[j] * this.k.eval(f.featureMaps[i].loc, this.featureMaps[j].loc);
	    }
	}
	return innerProduct;
    }

    /**
     * Computes set of samples from the space designed to permit optimal approximations over
     * arbitrary RKHS functions
     *
     * @param nSamples number of samples
     * @param lb lower bound on each dimension
     * @param ub upper bound on each dimension
     * @return
     */
    public double[][] getOptimalSamples(int nSamples, double[] lb, double[] ub) {
	double[][] currentSamples = gridSample(nSamples, lb, ub);
	double epsilonPct = 0.01;
	boolean converged = false;

	while (!converged) {
	    double currentVal = getSamplingObjFcnValue(currentSamples);
	    // pick a sample 
	    int sampleIndex = (int) (Math.random() * (currentSamples.length));
	    // pick a coordinate 
	    int sampleCoord = (int) (Math.random() * (currentSamples[sampleIndex].length));
	    // optimize that coordinate 
	    ToDoubleFunction<Double> pivotFunction
		    = (x) -> getSamplingObjFcnValue(changeSample(currentSamples, sampleIndex, sampleCoord, x));
	    currentSamples[sampleIndex][sampleCoord]
		    = findOptimalCoord(currentSamples[sampleIndex], pivotFunction, lb[sampleCoord], ub[sampleCoord]);
	    double newVal = getSamplingObjFcnValue(currentSamples);
	    if (((newVal - currentVal) / currentVal) < epsilonPct) {
		converged = true;
	    }
	}
	return currentSamples;
    }

    /**
     * Changes the sample by altering the value at the given array coordinates to match the value
     *
     * @param arr
     * @param coord1
     * @param coord2
     * @param val
     * @return
     */
    public double[][] changeSample(double[][] arr, int coord1, int coord2, double val) {
	double[][] newSamples = new double[arr.length][arr[0].length];
	for (int i = 0; i < arr.length; i++) {
	    for (int j = 0; j < arr.length; j++) {
		if (i == coord1 && j == coord2) {
		    newSamples[i][j] = val;
		} else {
		    newSamples[i][j] = arr[i][j];
		}
	    }
	}
	return newSamples;
    }

    /**
     * Performs a line search to find the value minimizing the function
     *
     * @param point
     * @param fcn
     * @param lBound
     * @param uBound
     * @return
     */
    public double findOptimalCoord(
	    double[] point,
	    ToDoubleFunction<Double> fcn,
	    double lBound,
	    double uBound) {
	double maxValue = Double.NEGATIVE_INFINITY;
	double maxCoord = lBound;
	for (double candidate = lBound; candidate <= uBound; candidate += (lBound - uBound) / 1000) {
	    double value = fcn.applyAsDouble(candidate);
	    if (value > maxValue) {
		maxValue = value;
		maxCoord = candidate;
	    }
	}
	return maxCoord;
    }

    /**
     * Computes the value of the function to be maximized to get optimal sampling set
     *
     * @param samples
     * @return
     */
    public double getSamplingObjFcnValue(double[][] samples) {
	// return tr(K^(1/2) * (k[x])^-1 * K^(1/2))

	int numSamples = samples.length;

	// compute gramMatrix
	double[][] gMat = new double[numSamples][numSamples];
	for (int i = 0; i < gMat.length; i++) {
	    for (int j = 0; j < gMat[0].length; j++) {
		gMat[i][j] = this.k.eval(samples[i], samples[j]);
	    }
	}
	Matrix gramMatrix = new Matrix(gMat);
	Matrix gramMatrixInverse = gramMatrix.inverse();

	// get RKHS representation of the Lebesgue measure for computing k matrix
	RKHSFunction lebesgueMeasure = getLebesgueMeasure(domainLB, domainUB);

	// fill in the k matrix
	double[][] kDoubleMat = new double[numSamples][numSamples];
	for (int i = 0; i < kDoubleMat.length; i++) {
	    for (int j = 0; j < kDoubleMat[i].length; j++) {
		// get phi_{xi} phi_{xj}
		FeatureMap phixi = new FeatureMap(this.k, samples[i]);
		FeatureMap phixj = new FeatureMap(this.k, samples[j]);
		FeatureMap[] phis = new FeatureMap[]{phixi, phixj};
		ToDoubleFunction<double[]> prodFunction
			= (point) -> (phixi.eval(point) * phixj.eval(point));
		double[] lcs = fitCoeffs(phis, prodFunction, this.k);
		RKHSFunction productFunction = new RKHSFunction(phis, lcs);

		// take inner product of productFunction and lebesgueMeasure
		kDoubleMat[i][j] = productFunction.innerProduct(lebesgueMeasure);
	    }
	}

	// trace k^1/2 K^-1 k^1/2
	Matrix kMatrix = new Matrix(kDoubleMat);
	// sqrt k = v * d^1/2 v^-1
	double[][] eigenVecMatrix = kMatrix.eig().getV().getArray();
	double[][] eigenValMatrix = kMatrix.eig().getD().getArray();

	for (int i = 0; i < eigenValMatrix.length; i++) {
	    if (eigenValMatrix[i][i] < 0) {
		throw new RuntimeException("Negative eigenvalue.");
	    }
	}

	double[][] sqrtEigenvalues = new double[eigenValMatrix.length][eigenValMatrix[0].length];
	for (int i = 0; i < sqrtEigenvalues.length; i++) {
	    sqrtEigenvalues[i][i] = Math.sqrt(eigenValMatrix[i][i]);
	}
	Matrix eVec = new Matrix(eigenVecMatrix);
	Matrix sqrtEVal = new Matrix(sqrtEigenvalues);
	Matrix sqrtK = eVec.times(sqrtEVal).times(eVec.inverse());

	Matrix resultMatrix = sqrtK.times(gramMatrixInverse).times(sqrtK);
	double tr = resultMatrix.trace();

	return tr;
    }

    /**
     * Generates the RKHSFunction approximating the uniform measure over the specified domain. 
     * This amounts to uniform sampling within the bounds. 
     * @param lBounds
     * @param uBounds
     * @return 
     */
    public RKHSFunction getLebesgueMeasure(double[] lBounds, double[] uBounds) { 
	int numDimensions= lBounds.length;
	
	int numUniformSamples = 1000;
	double[][] uniformSamples = new double[numUniformSamples][numDimensions];
	for (double[] uniformSample : uniformSamples) {
	    for (int dim = 0; dim < uniformSamples.length; dim++) {
		uniformSample[dim] = Math.random() * (uBounds[dim] - lBounds[dim]) + lBounds[dim];
	    }
	}
	FeatureMap[] measureFMs = new FeatureMap[numUniformSamples];
	for (int fm = 0; fm < measureFMs.length; fm++) {
	    measureFMs[fm] = new FeatureMap(this.k, uniformSamples[fm]);
	}
	double[] measureCoeffs = new double[numUniformSamples];
	for (int coeff = 0; coeff < measureCoeffs.length; coeff++) {
	    measureCoeffs[coeff] = 1 / numUniformSamples;
	}
	RKHSFunction lebesgueMeasure = new RKHSFunction(measureFMs, measureCoeffs);
	return lebesgueMeasure;
    }
}
