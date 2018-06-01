/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.voxq;

/**
 *
 * Evaluates \int_a^b exp(-f(x)/RT) dx by Simpson's rule
 * 
 * @author mhall44
 */
public abstract class BoltzmannIntegrator1D {
    
    double a,b;
    int numSlices=50;//preferably even.  
    
    public BoltzmannIntegrator1D(double a, double b){
        this.a = a;
        this.b = b;
    }
    
    public abstract double f(double x);
    
    
    private double evalBoltz(double x){
        //value of Boltzmann factor as function of integration variable
        return Math.exp( - f(x) / IntraVoxelSampler.RT );
    }
     
    
    public double doIntegral(){
        
        double sliceWidth = (b-a)/numSlices;
        double num = evalBoltz(a) + evalBoltz(b);
        double denom = 2;
        
        for(int slice=1; slice<numSlices; slice++){
            double weight = 4;
            if(slice%2==0)
                weight = 2;
            
            num += weight * evalBoltz(a + slice*sliceWidth);
            denom += weight;
        }
        
        double integ = num * (b-a) / denom;
        return integ;
    }

}
