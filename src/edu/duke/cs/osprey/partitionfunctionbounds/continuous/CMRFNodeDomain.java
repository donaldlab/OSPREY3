/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.util.function.ToDoubleFunction;

/**
 *
 * @author aditya
 */
public class CMRFNodeDomain {
    public double[] domainLB;
    public double[] domainUB;
    double volume;
    
    Kernel k;
    
    ToDoubleFunction<double[]> energyFunction;
    RKHSFunction energyRKHS;
    ToDoubleFunction<double[]> pdf;
    RKHSFunction probabilityRKHS;
    double constRT = PoissonBoltzmannEnergy.constRT;
    
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
	double energyInt = energyRKHS.computeIntegral();
	this.pdf = (x)->(Math.exp(-energyFunction.applyAsDouble(x)/constRT)/energyInt);
        this.probabilityRKHS = new RKHSFunction(k, lBound, uBound, this.pdf);
    }
    
    public double getEnergyIntegral() { 
	return energyRKHS.computeIntegral();
    }
    
    public double getProbAtPoint(double[] point) { 
	return this.pdf.applyAsDouble(point);
    }
    
    public double getEnergyAtPoint(double[] point) { 
	return this.energyFunction.applyAsDouble(point);
    }
    
    public boolean pointInDomain(double[] point) { 
	return this.k.validInput(point, point);
    }
}
