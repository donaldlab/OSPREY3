/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.HashMap;

/**
 * Represents a cMRF Node
 * 
 * @author aditya
 */
public class CMRFNode {
    // each node can have multiple domains (i.e. multiple disjoint domains of continuous flexibility) 
    public CMRFNodeDomain[] domains;
    public HashMap<CMRFNode, HashMap<CMRFNodeDomain, RKHSFunction>> outMessages; // messages for TRBP 
    public HashMap<CMRFNodeDomain, RKHSFunction> pseudomarginals; // pseudomarginals for TRBP
    public HashMap<CMRFNodeDomain, RKHSFunction> marginals; // marginals for SCMF 
    
    /**
     * Constructor -- a node is built off of its domains -- AND they all should have the same kernel 
     * @param doms 
     */
    public CMRFNode(CMRFNodeDomain[] doms) { 
	domains = doms;
        marginals = new HashMap<>();
        outMessages = new HashMap<>();
	pseudomarginals = new HashMap<>();
    }
    
    /**
     * Gets the domain corresponding to a particular point 
     * @param point
     * @return CMRFNodeDomain domain 
     */
    public CMRFNodeDomain getDomainForPoint(double[] point) { 
	for (CMRFNodeDomain d : domains ) {
	    if (d.isPointInDomain(point)) { 
		return d;
	    }
	}
	throw new RuntimeException("Point invalid.");
    }
    
    /**
     * Gets kernel for the CMRFNode
     * @return 
     */
    public Kernel getKernel() { 
	Kernel k = domains[0].k; 
	for (int i=1; i<domains.length; i++) { 
	    if (!(k.equals(domains[i].k))) { 
		throw new RuntimeException("Kernels aren't the same -- what are you doing?!");
	    }
	}
	return k;
    }
}
