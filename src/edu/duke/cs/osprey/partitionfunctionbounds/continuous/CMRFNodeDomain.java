/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.util.function.ToDoubleFunction;

/**
 * Represents a domain for a node in the cMRF -- this corresponds to a single region of continuous conformational 
 * flexibility, for example a single voxel. 
 * 
 * @author aditya
 */
public class CMRFNodeDomain {
    
    // domain has lower bounds, upper bounds, and corresponding volume 
    public double[] domainLB;
    public double[] domainUB;
    double volume;
    
    // kernel over the domain -- in general, this should always be the same for every domain across the MRF 
    Kernel k;
    
    // energy function and resulting PDF with their corresponding RKHS representations 
    ToDoubleFunction<double[]> energyFunction;
    RKHSFunction energyRKHS;
    ToDoubleFunction<double[]> pdf;
    RKHSFunction probabilityRKHS;

    // thermodynamic constants
    double constRT = PoissonBoltzmannEnergy.constRT;

    /**
     * Constructor accepts a set of lower bounds, a set of upper bounds, a kernel, and the energy function in the form
     * of a ToDoubleFunction<double[]> that takes a point within the bounded region and returns a double corresponding
     * to the intra-rotamer energy at that point 
     * @param lBound
     * @param uBound
     * @param k
     * @param eFunc 
     */
    public CMRFNodeDomain(
    		double[] lBound, 
    		double[] uBound, 
    		Kernel k, 
    		ToDoubleFunction<double[]> eFunc) { 
    	this.domainLB = lBound;
    	this.domainUB = uBound;

    	double v = 1.0;
    	for (int i=0; i<lBound.length; i++) { v *= uBound[i] - lBound[i]; }
    	volume = v;

    	this.k = k;
    	this.energyFunction = eFunc;
    	this.energyRKHS = new RKHSFunction(k, lBound, uBound, eFunc);
    	this.pdf = (x)->(Math.exp(-energyFunction.applyAsDouble(x)/constRT));
    	this.probabilityRKHS = new RKHSFunction(k, lBound, uBound, this.pdf);
    	double partFn = this.probabilityRKHS.computeIntegral();
    	this.probabilityRKHS = new RKHSFunction(k, lBound, uBound, (point)->(this.probabilityRKHS.eval(point)/partFn));

    }

    /**
     * Gets the pdf value at a particular point
     * @param point
     * @return probability
     */
    public double getProbAtPoint(double[] point) { 
	return this.pdf.applyAsDouble(point);
    }
    
    /**
     * Gets the value of the energy function at a particular point
     * @param point
     * @return energy
     */
    public double getEnergyAtPoint(double[] point) { 
	return this.energyFunction.applyAsDouble(point);
    }
    
    /**
     * Determines whether or not a function is in the node domain
     * @param point
     * @return true/false 
     */
    public boolean isPointInDomain(double[] point) { 
	return this.k.validInput(point, point);
    }
}
