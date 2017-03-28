package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.function.ToDoubleFunction;
import Jama.*;
import cern.colt.Arrays;
import java.util.ArrayList;
import java.util.Random;

/**
 * A reproducing kernel Hilbert space (RKHS) representation of a function that represents arbitrary
 * ToDoubleFunction<double[]>'s as linear combinations of FeatureMap's.
 * @author am439
 */
public class RKHSFunction {
    
    // a word of warning: this code here is a complete mess, and iff I have time I'll get to it, but in all honesty
    // the best bet for future code sanitation is probably a total rewrite of this steaming pile
    
    // kernel and domain bounds
    public Kernel k;
    public double[] domainLB;
    public double[] domainUB;  
    
    // list of feature maps and coefficients
    public FeatureMap[] featureMaps;
    public double[] coeffs;
    
    // the function we care about
    public ToDoubleFunction<double[]> referenceFunction;
    
    // number of samples per dimension -- we can change the later, and in any case we probably want to use a
    // hard cap to avoid having ridiculously long LCs
    public int numSamplesPerDimension = 10; // we can change this as the function changes
    
    // max number of samples we want to allow
    public int maxSamples = 50;
    
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
    
    public RKHSFunction(Kernel k, double[] domainLB, double[] domainUB, ToDoubleFunction<double[]> f, int x) { 
        int numDims = domainLB.length;
        this.numSamplesPerDimension = x;
        this.domainLB = domainLB;
        this.domainUB = domainUB;
        this.k = k;
        this.referenceFunction = f;
        this.gridSampleFit(k, f, domainLB, domainUB, true);
    }
    
    /**
     * Constructor -- given a kernel, domain, and function -- fits the RKHS function
     * This is the constructor you'll want to use the most.
     *
     * @param k
     * @param domainLB
     * @param domainUB
     * @param f
     */
    public RKHSFunction(
            Kernel k,
            double[] domainLB,
            double[] domainUB,
            ToDoubleFunction<double[]> f) {
        // constructs an RKHS approximation to a function
        // points are uniformly sampled with a fixed number of samples
        
        // we want to roughly end up taking 1000 samples
        //   take the floor to err on the side of fewer samples to prevent the
        //   computer from kind of blowing up
        int numDims = domainLB.length;
        this.numSamplesPerDimension =
                (int) Math.round(Math.ceil(Math.exp(Math.log(maxSamples)/numDims)));
        
        this.domainLB = domainLB;
        this.domainUB = domainUB;
        this.k = k;
        this.referenceFunction = f;
        this.gridSampleFit(k, f, domainLB, domainUB, true);
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
//        System.out.println("    Fitting function:");
//        System.out.println("      lower bounds: "+Arrays.toString(domainLB));
//        System.out.println("      upper bounds: "+Arrays.toString(domainUB));
//        System.out.println("      samples/dim: "+this.numSamplesPerDimension);
        
        
        int numSamples = this.numSamplesPerDimension; // samples per dimension
        
        double rmsd = Double.POSITIVE_INFINITY;
        boolean computedRMSDOnce = false;
        
        FeatureMap[] fMaps = new FeatureMap[numSamples];
        double[] fitCoeffs = new double[numSamples];
        
        while ((rmsd > 0.05 && !fixSampleNumber) || !computedRMSDOnce) {
            computedRMSDOnce = true;
            numSamples += 1;
            double[][] samples = gridSample(numSamples, domainLB, domainUB);
            fMaps = new FeatureMap[samples.length];
            for (int i = 0; i < samples.length; i++) {
                fMaps[i] = new FeatureMap(k, samples[i]);
            }
            fitCoeffs = this.fitCoeffs(fMaps, f, k);
        }
                
        this.featureMaps = fMaps;
        this.coeffs = fitCoeffs;
        this.domainLB = domainLB;
        this.domainUB = domainUB;
        this.k = k;
        this.referenceFunction = f;
        this.numSamplesPerDimension = numSamples;
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
     * Computes the RMSD of a fitted linear combination by sampling more densely than the function
     * gridSamples used to construct the RKHS representation. It's not a great way to check
     * RMSD, but I guess it'll do for now.
     *
     * Ten bucks says it's never getting better...
     *
     * Update (2017/01/17) -- for the love of all that is holy don't take these values seriously, they're wildly wrong
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
        
        // total samples should be 10 times the number of samples the function uses for the RKHS
        // representation
        int numSamples =
                (int) Math.exp((1.0/lBound.length)*
                        Math.log(10*Math.pow(this.numSamplesPerDimension, lBound.length)));
        double[][] testSamples = gridSample(numSamples, lBound, uBound);
        for (double[] sample : testSamples) {
            double funcValue = func.applyAsDouble(sample);
            if (funcValue == 0) { funcValue = 0.000000005; }
            double fitValue = evalLC(fmaps, coeffs, sample);
            
            double delta = (funcValue - fitValue)/funcValue;
            
            rmsd += Math.pow(delta, 2);
        }
        
        return Math.sqrt(rmsd*(1.0/numSamples));
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
        RKHSFunction lebesgueMeasure = getLebesgueMeasure(domainLB, domainUB);
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
        RKHSFunction subspaceMeasure = getLebesgueMeasure(domainLB, point);
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
     * Computes integral of this function over its domain
     * @return totalIntegral
     */
    public double computeIntegral() {

    	RKHSFunction measure = this.getLebesgueMeasure(this.domainLB, this.domainUB);

    	double domainVolume = 1.0;
    	for (int i=0; i<domainLB.length; i++) {
    		domainVolume = domainVolume * (domainUB[i] - domainLB[i]);
    	}
    	final double vol = domainVolume;

    	double totalIntegral = domainVolume * this.innerProduct(measure);

    	return totalIntegral;
    }
    
    public double computeAreaUnderCurve() { 
    	RKHSFunction measure = this.getLebesgueMeasure(this.domainLB, this.domainUB);

    	double domainVolume = 1.0;
    	for (int i=0; i<domainLB.length; i++) {
    		domainVolume = domainVolume * (domainUB[i] - domainLB[i]);
    	}
    	final double vol = domainVolume;
    	
    	double result = 0.0;
    	for (int i=0; i<featureMaps.length; i++) { 
    		for (int j=0; j<measure.featureMaps.length; j++) {
    			double val = this.coeffs[i] * measure.coeffs[j] * this.k.eval(featureMaps[i].loc, measure.featureMaps[j].loc); 
    			if (val > 0) { result += val; }
    		}
    	}
    	    	
    	return result;    	
    }

    /**
     * Computes the average of this function over its domain
     * @return expectation
     */
    public double computeExpectation() { 
        RKHSFunction measure = this.getLebesgueMeasure(this.domainLB, this.domainUB);
	return this.innerProduct(measure);
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
    	Random r = new Random(0);
    	int numUniformSamples = 10000;

    	double measureCoeff = 1.0/numUniformSamples;
    	double[][] uniformSamples = new double[numUniformSamples][numDimensions];
    	for (double[] uniformSample : uniformSamples) {
    		for (int dim = 0; dim < uniformSample.length; dim++) {
    			uniformSample[dim] = r.nextDouble() * (uBounds[dim] - lBounds[dim]) + lBounds[dim];
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
     * Given two RKHSFunctions, returns a new function over the Cartesian product of their domains where
     * the first set of coordinates are passed to the first function and the second set to the second function
     *
     * @param func1
     * @param func2
     * @param k
     * @return
     */
    public static RKHSFunction getCartesianProductFunction(
            RKHSFunction func1,
            RKHSFunction func2,
            Kernel k) {
        double[] newLB = CMRF.concatArrays(func1.domainLB, func2.domainLB);
        double[] newUB = CMRF.concatArrays(func1.domainUB, func2.domainUB);
        return new RKHSFunction(
                k,
                newLB,
                newUB,
                (point) -> (evaluateFunctionSplit(point, func1.domainLB.length, func1, func2))
        );
    }
    
    /**
     * Given a point concatenating coordinates from two domains with RKHS functions over those domains, returns the
     * product of those functions evaluated at the de-concatenated points
     * @param concatPoint
     * @param firstDimSize
     * @param func1
     * @param func2
     * @return
     */
    public static double evaluateFunctionSplit(
            double[] concatPoint,
            int firstDimSize,
            RKHSFunction func1,
            RKHSFunction func2) {
        ArrayList<double[]> points =  CMRF.splitArray(concatPoint, firstDimSize);
        double[] point1 = points.get(0);
        double[] point2 = points.get(1);
        return func1.eval(point1)*func2.eval(point2);
    }
    
    public Matrix dumpPoints() { 
        double[][] gridS = this.gridSample(numSamplesPerDimension, domainLB, domainUB);
        Matrix m = new Matrix(gridS.length, gridS[0].length+1);
        for (int sample = 0; sample<gridS.length; sample++) { 
            for (int coord=0; coord<gridS[sample].length; coord++) { 
                m.set(sample, coord, gridS[sample][coord]);
            }
            m.set(sample, gridS[sample].length, this.eval(gridS[sample]));
        }
        
        return m;        
    }
}
