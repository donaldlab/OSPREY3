/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.ToDoubleFunction;

/**
 * Represents an edge connecting two domains/nodes in the cMRF 
 * 
 * @author aditya
 */
public class CMRFEdge {
    CMRFNode node1;
    CMRFNode node2;
    
    CMRFEdgeDomain[] domainLinks;
    HashMap<CMRFEdgeDomain, RKHSFunction> pseudomarginals;
    
    /**
     * Constructs an edge given two nodes and a map of pairwise energy functions 
     * Energy function map takes a domain from node 1 and a domain from node 2 and outputs the corresponding function
     * @param n1
     * @param n2 
     * @param eFuncs 
     */
    public CMRFEdge(
            CMRFNode n1, 
            CMRFNode n2, 
            HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> eFuncs) { 
	this.node1 = n1;
	this.node2 = n2;
	
	ArrayList<CMRFEdgeDomain> doms = new ArrayList<>();
	
	for (CMRFNodeDomain d1 : node1.domains) { 
	    for (CMRFNodeDomain d2 : node2.domains) { 
		double[] lBounds = CMRFEdgeDomain.concatArrays(d1.domainLB, d2.domainLB);
		double[] uBounds = CMRFEdgeDomain.concatArrays(d1.domainUB, d2.domainUB);
		double[][] bounds = new double[lBounds.length][2];
		for (int i=0; i<lBounds.length; i++) { 
		    bounds[i][0] = lBounds[i];
		    bounds[i][1] = uBounds[i];
		}
		
		double diff = 0.0;
		for (int i=0; i<lBounds.length; i++) { diff += uBounds[i] - lBounds[i]; }
		diff = diff/lBounds.length;
		KernelGaussian prodK = new KernelGaussian(bounds, diff);
		
		doms.add(new CMRFEdgeDomain(
			d1.domainLB, d1.domainUB,
			d2.domainLB, d2.domainUB,
			d1.k, d2.k, prodK,
			eFuncs.get(d1).get(d2)));
	    }
	}
	
	domainLinks = new CMRFEdgeDomain[doms.size()];
	for (int i=0; i<domainLinks.length; i++) { 
	    domainLinks[i] = doms.get(i);
	}
	
	pseudomarginals = new HashMap<>();
    }
    
    /**
     * Gets the energy at a point 
     * @param point
     * @return 
     */
    public double getEnergyAtPoint(double[] point) { 
	ArrayList<double[]> coords = CMRF.splitArray(point, node1.domains[0].domainLB.length);
	double[] coord1 = coords.get(0);
	double[] coord2 = coords.get(1);
	for (CMRFEdgeDomain d : domainLinks) { 
	    if (d.isValidPoint(coord1, coord2)) { 
		return d.getEnergyAtPoint(point);
	    }
	}
	throw new RuntimeException("Point is not valid input to pairwise energy function.");
    }
    
    /**
     * Gets the edge domain corresponding two the two CMRFNodeDomains
     * @param d1
     * @param d2
     * @return 
     */
    public CMRFEdgeDomain getEdgeDomain(CMRFNodeDomain d1, CMRFNodeDomain d2) { 
        for (CMRFEdgeDomain edgeDomain : domainLinks) { 
            if (((Arrays.equals(d1.domainLB, edgeDomain.resOneLB)) &&
                    (Arrays.equals(d1.domainUB, edgeDomain.resOneUB)) &&
                    (Arrays.equals(d2.domainLB, edgeDomain.resTwoLB)) && 
                    (Arrays.equals(d2.domainUB, edgeDomain.resTwoUB))) ||                     
                ((Arrays.equals(d1.domainLB, edgeDomain.resTwoLB)) &&
                    (Arrays.equals(d1.domainUB, edgeDomain.resTwoUB)) &&
                    (Arrays.equals(d2.domainLB, edgeDomain.resOneLB)) && 
                    (Arrays.equals(d2.domainUB, edgeDomain.resOneUB)))) { 
                return edgeDomain;
            }
        }
        throw new RuntimeException("Couldn't find the required edge domain.");
    }
}
