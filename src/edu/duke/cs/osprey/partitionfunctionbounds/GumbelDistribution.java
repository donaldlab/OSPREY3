/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

/**
 *
 * @author hmn5
 */
public class GumbelDistribution {
    //CDF: F(x;mu,beta) = exp(- exp( -(x-mu)/beta)))
    //E[x] = mu + gamma*beta, gamme = 0.5772 (Euler-Mascheroni constant)
    
    final static public double gamma = 0.5772156649;
    

    //returns one random sample from the distribution
    public static double sample(double mu, double beta){
        return mu - beta*Math.log(- Math.log(Math.random()));
    }

    //For now, beta is fixed at 1
    public static double sampleTruncated(double mu, double maximum){
        return mu - Math.log(Math.exp(-maximum + mu) - Math.log(Math.random()));
    }
    
    //returns an array of N random samples
    public static double[] sampleN(int N, double mu, double beta){
        double[] samples = new double[N];
        for (int i=0; i<N; i++){
            samples[i] = GumbelDistribution.sample(mu, beta);
        }
        return samples;
    }
    
}
