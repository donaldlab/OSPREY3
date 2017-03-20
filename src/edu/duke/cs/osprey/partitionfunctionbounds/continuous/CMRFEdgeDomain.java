/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import Jama.Matrix;
import cern.colt.Arrays;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.util.function.ToDoubleFunction;

/**
 * Represents a domain for an edge of the MRF, corresponding to a continuous region of pairwise conformational space 
 * 
 * @author aditya
 */
public class CMRFEdgeDomain {
    
    // domains for the two residues
    CMRFNodeDomain resOneDomain;
    CMRFNodeDomain resTwoDomain;
    
    // bounds for the two residues 
    double[] resOneLB;
    double[] resOneUB;
    double[] resTwoLB;
    double[] resTwoUB;
    
    double[] domainLB;
    double[] domainUB;
    
    // volme for both residues and total volume 
    double resOneVolume;
    double resTwoVolume;
    double volume;
    
    // kernel for both residues, and a kernel for the pairwise domain 
    // probably want to make all of these the same 
    Kernel resOneK;
    Kernel resTwoK;
    Kernel resAllK;
    
    // energy and probability functions 
    ToDoubleFunction<double[]> eFunc; //eFunc and pFunc take concatenations of the coordinates
    ToDoubleFunction<double[]> pFunc;
    RKHSFunction eFuncRKHS;
    RKHSFunction pFuncRKHS;
    
    RKHSFunction pseudomarginal;
    
    double constRT = PoissonBoltzmannEnergy.constRT;
    
    /**
     * Constructor analogous to the one in CMRFNodeDomain. 
     * @param rOneLB
     * @param rOneUB
     * @param rTwoLB
     * @param rTwoUB
     * @param rOneK
     * @param rTwoK
     * @param prodK
     * @param energyFunc 
     */
    public CMRFEdgeDomain(
	    double[] rOneLB, double[] rOneUB,
	    double[] rTwoLB, double[] rTwoUB,
	    Kernel rOneK, Kernel rTwoK, Kernel prodK,
	    CMRFNodeDomain dOne, CMRFNodeDomain dTwo,
	    ToDoubleFunction<double[]> energyFunc) { 
	
        // we store the upper and lower bounds for both domains involved in the edge
	this.resOneLB = rOneLB;
	this.resOneUB = rOneUB;
	this.resTwoLB = rTwoLB;
	this.resTwoUB = rTwoUB;
	
	this.domainLB = concatArrays(resOneLB, resTwoLB);
	this.domainUB = concatArrays(resOneUB, resTwoUB);
	
        // we also store individual/collective volumes 
	this.resOneVolume = CMRFEdgeDomain.getVolume(rOneLB, rOneUB);
	this.resTwoVolume = CMRFEdgeDomain.getVolume(rTwoLB, rTwoUB);
	this.volume = resOneVolume * resTwoVolume;
	
        // we store individual and pairwise kernels -- for simplicity, we may want to keep these the same for now 
	this.resOneK = rOneK;
	this.resTwoK = rTwoK;
	this.resAllK = prodK;
	
	// screw it let's just store the bloody domains
	this.resOneDomain = dOne;
	this.resTwoDomain = dTwo;
	
        // energy function, pdf, and associated RKHS representations 
	this.eFunc = energyFunc;
	this.eFuncRKHS = new RKHSFunction(
		resAllK, 
		CMRFEdgeDomain.concatArrays(resOneLB, resTwoLB),
		CMRFEdgeDomain.concatArrays(resOneUB, resTwoUB),
		eFunc);
	double totalEnergy = eFuncRKHS.computeIntegral();
        
        
	this.pFunc = (point)->(Math.exp(-eFunc.applyAsDouble(point)/constRT)/totalEnergy);
        this.pFuncRKHS = new RKHSFunction(
                resAllK,
                CMRFEdgeDomain.concatArrays(resOneLB, resTwoLB),
                CMRFEdgeDomain.concatArrays(resOneUB, resTwoUB),
                pFunc);
        double edgeZ = pFuncRKHS.computeIntegral();
        pFuncRKHS = new RKHSFunction(
                resAllK,
                CMRFEdgeDomain.concatArrays(resOneLB, resTwoLB),
                CMRFEdgeDomain.concatArrays(resOneUB, resTwoUB),
                (point)->Math.max(pFuncRKHS.eval(point)/edgeZ, 0));
    }
    
    /**
     * Concatenates two arrays
     * @param arrOne
     * @param arrTwo
     * @return concatenated array
     */
    public static double[] concatArrays(double[] arrOne, double[] arrTwo) { 
	double[] arr = new double[arrOne.length+arrTwo.length];
        System.arraycopy(arrOne, 0, arr, 0, arrOne.length);
        System.arraycopy(arrTwo, 0, arr, arrOne.length, arrTwo.length);
	return arr;
    }
    
    /**
     * Returns the volume of the region with given bounds
     * @param LB
     * @param UB
     * @return volume
     */
    public static double getVolume(double[] LB, double[] UB) {
	if (LB.length != UB.length) {
	    throw new RuntimeException("Uneven dimensions.");
	}
	double vol = 1.0;
	for (int i=0; i<LB.length; i++) { 
	    vol *= UB[i] - LB[i];
	}
	return vol;
    }
    
    /**
     * Returns the energy at a particular point 
     * @param point
     * @return 
     */
    public double getEnergyAtPoint(double[] point) { 
	return this.eFunc.applyAsDouble(point);
    }
    
    /**
     * Checks if the coordinates are valid/within the confspace 
     * @param point1
     * @param point2
     * @return 
     */
    public boolean isValidPoint(double[] point1, double[] point2) { 
	boolean valid =  resOneDomain.isPointInDomain(point1) && resTwoDomain.isPointInDomain(point2);
	if (!valid) { System.out.println("\ninvalid point: \n" +
		"\tpoints:  " + Arrays.toString(point1) + ", " + Arrays.toString(point2) + "\n" + 
		"\tlbounds: " + Arrays.toString(this.resOneLB) + Arrays.toString(this.resTwoLB) + "\n" + 
		"\tubounds: " + Arrays.toString(this.resOneUB) + Arrays.toString(this.resTwoUB) + "\n"); } 
	return valid;
    }
    
    
}
