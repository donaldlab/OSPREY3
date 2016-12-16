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
public class CMRFEdgeDomain {
    double[] resOneLB;
    double[] resOneUB;
    double[] resTwoLB;
    double[] resTwoUB;
    
    double resOneVolume;
    double resTwoVolume;
    double volume;
    
    Kernel resOneK;
    Kernel resTwoK;
    Kernel resAllK;
    
    ToDoubleFunction<double[]> eFunc; //eFunc and pFunc take concatenations of the coordinates
    ToDoubleFunction<double[]> pFunc;
    RKHSFunction eFuncRKHS;
    
    double constRT = PoissonBoltzmannEnergy.constRT;
    
    public CMRFEdgeDomain(
	    double[] rOneLB, double[] rOneUB,
	    double[] rTwoLB, double[] rTwoUB,
	    Kernel rOneK, Kernel rTwoK, Kernel prodK,
	    ToDoubleFunction<double[]> energyFunc) { 
	
	this.resOneLB = rOneLB;
	this.resOneUB = rOneUB;
	this.resTwoLB = rTwoLB;
	this.resTwoUB = rTwoUB;
	
	this.resOneVolume = this.getVolume(rOneLB, rOneUB);
	this.resTwoVolume = this.getVolume(rTwoLB, rTwoUB);
	this.volume = resOneVolume * resTwoVolume;
	
	this.resOneK = rOneK;
	this.resTwoK = rTwoK;
	this.resAllK = prodK;
	
	this.eFunc = energyFunc;
	this.eFuncRKHS = new RKHSFunction(
		resAllK, 
		this.concatArrays(resOneLB, resTwoLB),
		this.concatArrays(resOneUB, resTwoUB),
		eFunc);
	double totalEnergy = eFuncRKHS.computeIntegral();
	this.pFunc = (point)->(Math.exp((-1*eFunc.applyAsDouble(point))/constRT)/totalEnergy);
    }
    
    public static double[] concatArrays(double[] arrOne, double[] arrTwo) { 
	double[] arr = new double[arrOne.length+arrTwo.length];
	for (int i=0; i<arrOne.length; i++) { 
	    arr[i] = arrOne[i];
	}
	for (int i=0; i<arrTwo.length; i++) { 
	    arr[i+arrOne.length] = arrTwo[i];
	}
	return arr;
    }
    
    public double getVolume(double[] LB, double[] UB) {
	if (LB.length != UB.length) {
	    throw new RuntimeException("Uneven dimensions.");
	}
	double vol = 1.0;
	for (int i=0; i<LB.length; i++) { 
	    vol *= UB[i] - LB[i];
	}
	return vol;
    }
    
    public double getEnergyAtPoint(double[] point) { 
	return this.eFunc.applyAsDouble(point);
    }
    
    public boolean isValidPoint(double[] point1, double[] point2) { 
	for (int i=0; i<point1.length; i++) { 
	    if ((resOneLB[i] > point1[i]) || (resOneUB[i] < point1[i])) { 
		return false;
	    }
	}
	
	for (int i=0; i<point2.length; i++) { 
	    if ((resTwoLB[i] > point2[i]) || (resTwoUB[i] < point2[i])) { 
		return false;
	    }
	}
	
	return true;
    }
    
}
