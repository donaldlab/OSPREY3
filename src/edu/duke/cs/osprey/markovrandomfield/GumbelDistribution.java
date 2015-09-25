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
    
    final static public double gamma = 0.5772156649;
    //returns one random sample from the distribution defined by mu and beta
    public static double sample(double mu, double beta){
        return mu - beta*Math.log(- Math.log(Math.random()));
    }
}