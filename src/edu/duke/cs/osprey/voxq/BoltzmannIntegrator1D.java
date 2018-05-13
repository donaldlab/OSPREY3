/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

