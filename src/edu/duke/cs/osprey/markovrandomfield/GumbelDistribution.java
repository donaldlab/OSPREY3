/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.markovrandomfield;

/**
 *
 * @author hmn5
 */
public class GumbelDistribution {
    //CDF: F(x;mu,beta) = exp(- exp( -(x-mu)/beta)))
    //E[x] = mu + gamma*beta, gamme = 0.5772 (Euler-Mascheroni constant)
    
    double mu;
    double beta;
    final static public double gamma = 0.5772156649;
    
    public GumbelDistribution(double mu, double beta){
        this.mu = mu;
        this.beta = beta;
    }
    
    //returns one random sample from the distribution
    public double sample(){
        return mu - beta*Math.log(- Math.log(Math.random()));
    }
    
    //returns an array of N random samples
    public double[] sampleN(int N){
        double[] samples = new double[N];
        for (int i=0; i<N; i++){
            samples[i] = this.sample();
        }
        return samples;
    }
    
}
