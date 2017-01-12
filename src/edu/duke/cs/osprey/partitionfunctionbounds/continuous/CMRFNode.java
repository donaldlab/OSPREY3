/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import java.util.ArrayList;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
 *
 * @author aditya
 */
public class CMRFNode {
    public CMRFNodeDomain[] domains;
    public Map<CMRFNode, Map<CMRFNodeDomain, RKHSFunction>> outMessages;
    
    public CMRFNode(CMRFNodeDomain[] doms) { 
	domains = doms;
    }
    
    public CMRFNodeDomain getDomainForPoint(double[] point) { 
	for (CMRFNodeDomain d : domains ) {
	    if (d.pointInDomain(point)) { 
		return d;
	    }
	}
	throw new RuntimeException("Point invalid.");
    }
    
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
