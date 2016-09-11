package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.function.ToDoubleFunction;
import Jama.*;

public class RKHSFunction {

    public Kernel k;
    public double[] domainLB;
    public double[] domainUB;

    public FeatureMap[] featureMaps;
    public double[] coeffs;
    
    public ToDoubleFunction<double[]> referenceFunction;

    public int numSamplesPerDimension = 10; // we can change this as the function changes
    
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
	Kernel firstK = featureMaps[0].getKernel();
	for (FeatureMap fmap : featureMaps) {
	    Kernel fK = fmap.getKernel();
	    if (!firstK.equals(fK)) {
		throw new RuntimeException("RKHS Function kernels don't match.");
	    }
	}

	this.featureMaps = featureMaps;
	this.coeffs = coeffs;

	// set up kernel, check for consistency 
	this.k = firstK;
	this.domainLB = new Matrix(firstK.bounds).transpose().getArray()[0];
	this.domainUB = new Matrix(firstK.bounds).transpose().getArray()[1];
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
	this.gridSampleFit(k, f, domainLB, domainUB);
	this.referenceFunction = f;
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
	Matrix gramMatrix = new Matrix(fMaps.length, fMaps.length);
	for (int i=0; i<gramMatrix.getRowDimension(); i++) { 
	    for (int j=0; j<gramMatrix.getColumnDimension(); j++) { 
		double kernelValue = k.eval(fMaps[i].getLoc(), fMaps[j].getLoc());
		double perturbation = Math.random() * (kernelValue*0.001);
		gramMatrix.set(i, j, kernelValue+perturbation);
	    }
	}
	
	double[] fV = new double[fMaps.length];
	for (int i=0; i<fV.length; i++) { 
	    fV[i] = func.applyAsDouble(fMaps[i].getLoc());
	}
	Matrix fVals = new Matrix(fV, fMaps.length);
	Matrix ans;
	try {
	    ans = gramMatrix.solve(fVals);  
	} catch(Exception e) { 
	    throw new RuntimeException("Coefficient fit failed.");
	}
	return  ans.transpose().getArray()[0];	    
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
     * Grid samples the domain and fits the linear combination of feature maps to the provided
     * function. If fixSampleNumber is false, then the number of feature maps is increased until 
     * a specific accuracy level is reached. 
     * @param k
     * @param f
     * @param domainLB
     * @param domainUB
     * @param fixSampleNumber 
     */
    public void gridSampleFit(
	    Kernel k, 
	    ToDoubleFunction<double[]> f, 
	    double[] domainLB, 
	    double[] domainUB, 
	    boolean fixSampleNumber) {
	System.out.println("Beginning fit...");

	
	int numSamples = this.numSamplesPerDimension; // samples per dimension
	
	double rmsd = Double.POSITIVE_INFINITY;

	FeatureMap[] fMaps = new FeatureMap[numSamples];
	double[] fitCoeffs = new double[numSamples];

	while (rmsd > 0.05 && !fixSampleNumber) {
	    numSamples += 1;
	    double[][] samples = gridSample(numSamples, domainLB, domainUB);
	    fMaps = new FeatureMap[samples.length];
	    for (int i = 0; i < samples.length; i++) {
		fMaps[i] = new FeatureMap(k, samples[i]);
	    }
	    fitCoeffs = this.fitCoeffs(fMaps, f, k);
	    
	    rmsd = computeFitRMSDRatio(fMaps, fitCoeffs, domainLB, domainUB, f);
	    
	    if (Double.isNaN(rmsd)) { 
		throw new RuntimeException("RMSD is NaN.");
	    }
	    
	    System.out.println("With " + fMaps.length + " functions RMSD is: " + rmsd);
	}

	System.out.println("Used " + fMaps.length + " functions, total RMSD: " + rmsd + ".");

	this.featureMaps = fMaps;
	this.coeffs = fitCoeffs;
	this.domainLB = domainLB;
	this.domainUB = domainUB;
	this.k = k;
	this.referenceFunction = f;
	this.numSamplesPerDimension = numSamples;
    }

    /** 
     * Grid samples the domain and fits the linear combination of feature maps to the provided 
     * function. Increases the number of samples until a sufficient accuracy condition is reached. 
     * @param k
     * @param f
     * @param domainLB
     * @param domainUB 
     */
    private void gridSampleFit(
	    Kernel k, 
	    ToDoubleFunction<double[]> f,
	    double[] domainLB,
	    double[] domainUB) { 
	this.gridSampleFit(k, f, domainLB, domainUB, false);
    }
    
    /**
     * Grid samples from n dimensions -- row is sample, column is dimension
     *
     * @param numSamplesPerDimension
     * @param lb
     * @param ub
     * @return
     */
    public double[][] gridSample(int numSamplesPerDimension, double[] lb, double[] ub) {
	int numDimensions = lb.length;
	int samplesPerDim = numSamplesPerDimension;
	
	int[][] coords = getCoordinates(new int[0][0], samplesPerDim, numDimensions);
		
	double[][] samples = new double[coords.length][numDimensions];
	
	for (int i=0; i<samples.length; i++) { 
	    for (int j=0; j<samples[0].length; j++) { 
		samples[i][j] = lb[j] + (ub[j]-lb[j])*(((double) coords[i][j])/(samplesPerDim-1));
	    }
	}
	
	return samples;
    }
    
    /**
     * Returns a list of coordinates of cells in a matrix of length n with numDimensions dimensions
     * That is, this returns the list of (numDimensions)-tuples in {0...n-1}^(numDimensions)
     * 
     * @param coords
     * @param n
     * @param numDimensions
     * @return 
     */
    public int[][] getCoordinates(int[][] coords, int n, int numDimensions) { 
	if (numDimensions == 0) { 
	    return coords;
	}
	if (coords.length == 0) { // blank array, so just seed it
	    coords = new int[n][1];
	    for (int i=0; i<coords.length; i++) { 
		coords[i][0] = i;
	    }
	    return getCoordinates(coords, n, numDimensions-1);
	}
	
	int[] concats = new int[n];
	for (int i=0; i<concats.length; i++) { concats[i] = i; }
		
	int[][] newCoords = new int[coords.length * n][coords[0].length+1];
	for (int i=0; i<newCoords.length; i++) { 
	    int oldCoordIndex = (i / n);
	    int newCoordIndex = (i % n);
	    System.arraycopy(coords[oldCoordIndex], 0, newCoords[i], 0, newCoords[i].length-1);
	    newCoords[i][newCoords[i].length-1] = concats[newCoordIndex];
	}
	return getCoordinates(newCoords, n, numDimensions-1);
    }

    /**
     * Computes the accuracy of this RKHS function as compared to a reference function passed in
     * @param func
     * @return 
     */
    public double computeRKHSAccuracy(ToDoubleFunction<double[]> func) { 
	return this.computeFitRMSDRatio(
		this.featureMaps, 
		this.coeffs, 
		this.domainLB, 
		this.domainUB, 
		func);
    }
        
    /**
     * Computes the RMSD of a fitted linear combination by sampling more densely than the function
     * gridSamples used to construct the RKHS representation. It's not a great way to check
     * RMSD, but I guess it'll do for now. 
     * 
     * Ten bucks says it's never getting better... 
     *
     * @param fmaps
     * @param coeffs
     * @param lBound
     * @param uBound
     * @param func
     * @return
     */
    public double computeFitRMSDRatio(
	    FeatureMap[] fmaps,
	    double[] coeffs,
	    double[] lBound,
	    double[] uBound,
	    ToDoubleFunction<double[]> func) {
	double rmsd = 0.0;
	double funcSum = 0.0;
	
	// total samples should be 10 times the number of samples the function uses for the RKHS
	// representation
	int numSamples =  
		(int) Math.exp((1.0/lBound.length)*
			Math.log(10*Math.pow(this.numSamplesPerDimension, lBound.length)));
	double[][] testSamples = gridSample(numSamples, lBound, uBound);
	for (double[] sample : testSamples) {
	    double funcValue = func.applyAsDouble(sample);
	    double fitValue = evalLC(fmaps, coeffs, sample);
	    
	    double delta = (funcValue - fitValue);
	    
	    funcSum += Math.pow(funcValue, 2);
	    rmsd += Math.pow(delta, 2);
	}
	
	double rmsdRatio = 
		Math.sqrt(rmsd / testSamples.length)/Math.sqrt(funcSum / testSamples.length);
	
	return Math.sqrt(rmsdRatio);
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
	// we need to sample i.i.d. from this function -- so we use inverse transform sampling
	// and an RKHS-based approximation to the integral of the function to binary search 
	// the domain for the function with CDF equal to a value sampled u.a.r. from the unit
	// cube 
	
	int numFeatureMaps = this.featureMaps.length;

	// compute value of the integral of the mass function over the whole domain
	int numDimensions = this.domainLB.length;
	double domainVolume = 1.0;
	for (int i=0; i<numDimensions; i++) { 
	    domainVolume = domainVolume * (domainUB[i] - domainLB[i]);
	}
	RKHSFunction lebesgueMeasure = getGridMeasure(domainLB, domainUB);
	double functionIntegral = domainVolume * this.innerProduct(lebesgueMeasure);
	
	// set up array of new feature maps for the mean map 
	FeatureMap[] newFeatureMaps = new FeatureMap[numFeatureMaps];
	
	double[] cdfPoint;
	for (int fMapIndex = 0; fMapIndex < numFeatureMaps; fMapIndex++) {
	    cdfPoint = new double[numDimensions];
	    for (int i = 0; i < numDimensions; i++) {
		cdfPoint[i] = Math.random(); // u.a.r. CDF sampling
	    }
	    
	    double cdfValue = 1.0;
	    for (double coord : cdfPoint) { cdfValue = cdfValue * coord; }

//    	    System.out.println("Getting feature map "+(fMapIndex+1)+"/"+numFeatureMaps+" ("+cdfValue+")");
	    cdfPoint = binarySearchCDF(cdfPoint, domainLB, domainUB, functionIntegral);
	    newFeatureMaps[fMapIndex] = new FeatureMap(this.k, cdfPoint);
//	    System.out.println("Got feature map "+(fMapIndex+1)+"/"+numFeatureMaps);
	}

	ToDoubleFunction<double[]> thisFunction = (point)->(this.eval(point));
	double[] newCoeffs = fitCoeffs(newFeatureMaps, thisFunction, this.k);
	return new RKHSFunction(newFeatureMaps, newCoeffs);
    }

    /**
     * Returns the RKHS function with feature maps at "optimal" locations
     * @return RKHSFunction 
     */
    public RKHSFunction convertToOptimalSamples() { 
	int samplesPerDim = this.numSamplesPerDimension;
	double[][] optimalSamples = this.getOptimalSamples(samplesPerDim, domainLB, domainUB);
	FeatureMap[] fmaps = new FeatureMap[this.featureMaps.length];
	for (int i=0; i<fmaps.length; i++) { 
	    fmaps[i] = new FeatureMap(this.k, optimalSamples[i]);
	}
	double[] newCoeffs = this.fitCoeffs(fmaps, this.referenceFunction, k);
	return new RKHSFunction(fmaps, newCoeffs);
    }
    
    /**
     * Binary searches the domain to find a point whose CDF matches the value provided
     * @param cdfPoint
     * @param LBs
     * @param UBs
     * @param totalIntegral
     * @return 
     */
    public double[] binarySearchCDF(
	    double[] cdfPoint,
	    double[] LBs,
	    double[] UBs,
	    double totalIntegral) {

	int numDimensions = UBs.length;
	double epsilon = 0.005;
	double[] point = new double[UBs.length];
	for (int i=0; i<numDimensions; i++) { point[i] = 0.5*(UBs[i] + LBs[i]); }
	
	// compute CDF of the point 
	RKHSFunction subspaceMeasure = getGridMeasure(domainLB, point);
	double subspaceVolume = 1.0;
	for (int i=0; i<numDimensions; i++) { 
	    subspaceVolume = subspaceVolume * (point[i] - domainLB[i]);
	}
	double pointCDF = this.innerProduct(subspaceMeasure) * subspaceVolume / totalIntegral;
	double cdfVal = 1.0; for (double c : cdfPoint) { cdfVal = cdfVal * c; }
		
	if (Math.abs((pointCDF - cdfVal)/cdfVal) <= epsilon) { // if the point is close enough 
	    return point;
	} else if (pointCDF < cdfVal) {
	    return binarySearchCDF(cdfPoint, point, UBs, totalIntegral);
	} else if (pointCDF > cdfVal) { 
	    return binarySearchCDF(cdfPoint, LBs, point, totalIntegral);
	} else {
	    throw new RuntimeException("The point has CDF neither greater than, equal to, or "
		    + "less than the required CDF value. Congraulations, you broke math. ");
	}
    }

    /**
     * Computes the inner product of this RKHSFunction with another RKHSFunction. Assumes that the 
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
		    = findOptimalCoord(pivotFunction, lb[sampleCoord], ub[sampleCoord]);
	    double newVal = getSamplingObjFcnValue(currentSamples); 
	    double error = newVal - currentVal;
	    if (error < epsilonPct) {
		converged = true;
	    }
	}
	return currentSamples;
    }

    /**
     * Changes the sample by altering the value at the given array coordinates to match the value. 
     * I know this is the kludgiest thing ever. 
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
	    for (int j = 0; j < arr[0].length; j++) {
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
     * @param fcn
     * @param lBound
     * @param uBound
     * @return
     */
    public double findOptimalCoord(ToDoubleFunction<Double> fcn, double lBound, double uBound) {
	double maxValue = Double.NEGATIVE_INFINITY;
	double maxCoord = lBound;
	for (double candidate = lBound; candidate < uBound; candidate += (uBound - lBound) / 100) {
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
	RKHSFunction lebesgueMeasure = getGridMeasure(domainLB, domainUB);

	// fill in the k matrix
	double[][] kDoubleMat = new double[numSamples][numSamples];
	for (int i = 0; i < kDoubleMat.length; i++) {
	    for (int j = 0; j < kDoubleMat[i].length; j++) {
		// get phi_{xi} phi_{xj}
		FeatureMap phixi = new FeatureMap(this.k, samples[i]);
		FeatureMap phixj = new FeatureMap(this.k, samples[j]);
		FeatureMap[] phis = new FeatureMap[2];
		phis[0] = phixi; phis[1] = phixj;
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
     * Computes integral of this function over its domain 
     * @return totalIntegral
     */
    public double computeIntegral() { 
	System.out.println("Computing integral...");
	RKHSFunction measure = this.getGridMeasure(this.domainLB, this.domainUB);
	double domainVolume = 1.0;
		
	for (int i=0; i<domainLB.length; i++) { 
	    domainVolume = domainVolume * (domainUB[i] - domainLB[i]);
	}
	double totalIntegral = domainVolume * this.innerProduct(measure);
	
	return totalIntegral;
    }

    /**
     * Generates the RKHSFunction approximating the uniform measure over the specified domain. 
     * This amounts to uniform sampling within the bounds. 
     * This really doesn't work too well - gridMeasure seems to be a better approximation 
     * 
     * @param lBounds
     * @param uBounds
     * @return 
     */
    public RKHSFunction getLebesgueMeasure(double[] lBounds, double[] uBounds) { 
	int numDimensions= lBounds.length;
	
	int numUniformSamples = 1000;
	double measureCoeff = 1.0/numUniformSamples;
	double[][] uniformSamples = new double[numUniformSamples][numDimensions];
	for (double[] uniformSample : uniformSamples) {
	    for (int dim = 0; dim < uniformSample.length; dim++) {
		uniformSample[dim] = Math.random() * (uBounds[dim] - lBounds[dim]) + lBounds[dim];
	    }
	}
	FeatureMap[] measureFMs = new FeatureMap[numUniformSamples];
	for (int fm = 0; fm < measureFMs.length; fm++) {
	    measureFMs[fm] = new FeatureMap(this.k, uniformSamples[fm]);
	}
	double[] measureCoeffs = new double[numUniformSamples];
	for (int coeff = 0; coeff < measureCoeffs.length; coeff++) {
	    measureCoeffs[coeff] = measureCoeff;
	}
	
	RKHSFunction lebesgueMeasure = new RKHSFunction(measureFMs, measureCoeffs);
	return lebesgueMeasure;
    }
    
    /**
     * Approximates the Lebesgue measure over the domain using feature maps that uniformly grid
     * the domain. This isn't as theoretically pretty as getLebesgueMeasure, but it seems to provide 
     * superior approximations. 
     * @param lBounds
     * @param uBounds
     * @return 
     */
    public RKHSFunction getGridMeasure(double[] lBounds, double[] uBounds) {  
	int numDimensions = lBounds.length;
	double[][] samples = this.gridSample(this.numSamplesPerDimension, lBounds, uBounds);
	FeatureMap[] measureFMs = new FeatureMap[samples.length];
	for (int fm = 0; fm < measureFMs.length; fm++) {
	    measureFMs[fm] = new FeatureMap(this.k, samples[fm]);
	}
	
	double domainVolume = 1.0;
	for (int i=0; i<numDimensions; i++) { 
	    domainVolume = domainVolume * (uBounds[i] - lBounds[i]);
	}
	final double dVol = domainVolume;
	
	ToDoubleFunction<double[]> measureFunc = (point)->(1.0/dVol);
	
	double[] measureCoeffs = fitCoeffs(measureFMs, measureFunc, this.k);
	
	RKHSFunction lebesgueMeasure = new RKHSFunction(measureFMs, measureCoeffs);
	return lebesgueMeasure;
	
    }
}
