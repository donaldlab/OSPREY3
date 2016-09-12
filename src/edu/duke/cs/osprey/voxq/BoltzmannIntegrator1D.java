/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
