/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import java.io.Serializable;
import java.util.ArrayList;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.ConjugateGradient;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * Regularized tuple fitter
 * Using Lagrange multipliers, it minimizes the sum of squares of the tuple coefficients
 * subject to them satisfying the normal equations
 * 
 * @author mhall44
 */
public class RegTupleFitter implements Serializable {
    
    
    //Start with conj grad on normal eqs...
    
    
    RealLinearOperator AtA;
    RealLinearOperator MtM;//left hand-side operator for regularized fitting
    RealVector Atb;
    RealVector Mtb0;//right-hand side of regularized fitting
    
    //fit A * params = b
    //So solve A^T A params = A^T b 
    //(A is matrix defined by samp, b is true energies)
    //NOW SOLVING MtM[x;lambda] = Mt[x;lambda]
    
    int numSamp, numTup;
    ArrayList<ArrayList<Integer>> sampTuples;
    
    public RegTupleFitter (ArrayList<ArrayList<Integer>> sampleTuples, int numTuples, double[] trueVals){
        //sampleTuples.get(s) lists tuple numbers for sample s
        
        numSamp = sampleTuples.size();
        numTup = numTuples;
        sampTuples = sampleTuples;
        
        
        AtA = new RealLinearOperator(){

            @Override
            public int getRowDimension() {
                return numTup;
            }

            @Override
            public int getColumnDimension() {
                return numTup;
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                //first apply A
                double Arv[] = new double[numSamp];
                for(int s=0; s<numSamp; s++){
                    ArrayList<Integer> sampTup = sampTuples.get(s);
                    for(int t : sampTup)
                        Arv[s] += rv.getEntry(t);
                }
                
                //then apply A^T to Arv
                double ans[] = new double[numTup];
                for(int s=0; s<numSamp; s++){
                    ArrayList<Integer> sampTup = sampTuples.get(s);
                    for(int t : sampTup)
                        ans[t] += Arv[s];
                }
                
                return new ArrayRealVector(ans,false);//make RealVector without copying ans
            }
            
        };
        
        
        double atb[] = new double[numTup];
        //apply A^T to true vals
        for(int s=0; s<numSamp; s++){
            ArrayList<Integer> sampTup = sampTuples.get(s);
            for(int t : sampTup)
                atb[t] += trueVals[s];
        }
        
        Atb = new ArrayRealVector(atb);
        
        
        //Build "big" operators
        MtM = new RealLinearOperator() {

            @Override
            public int getRowDimension() {
                return 2*numTup;
            }

            @Override
            public int getColumnDimension() {
                return 2*numTup;
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                
                //rv is x followed by lambda
                RealVector x = rv.getSubVector(0, numTup);
                RealVector lambda = rv.getSubVector(numTup, numTup);
                
                RealVector opx = AtA.operate(AtA.operate(x));
                opx = opx.add( x.mapMultiply(4) );
                opx = opx.add( AtA.operate(lambda).mapMultiply(-2) );
                
                RealVector opl = AtA.operate(x).mapMultiply(-2);
                opl = opl.add( AtA.operate(AtA.operate(lambda)) );
                
                return opx.append(opl);
            }
            
        };
        
        
        //Mtb0 is [A^TAA^Tb; 0]
        Mtb0 = AtA.operate(Atb).append( new ArrayRealVector(new double[numTup],false) );
    }
    
    
    double[] doFit(){
        //return fit tuple coefficients
        
        ConjugateGradient cg = new ConjugateGradient(100000,1e-6,false);//max_iter; delta; whether to check pos def
        //delta is target ratio of residual norm to true vals norm
        
        long startTime = System.currentTimeMillis();
        
        RealVector ans = cg.solve(MtM, Mtb0).getSubVector(0, numTup);
        
        System.out.println( "Regularized conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
        
        double norm2 = ans.getNorm();
        double residM = AtA.operate(ans).getDistance(Atb);
        
        System.out.println( "Regularized fitting norm: " + norm2 + " Normal eq residual: " + residM );
        
        
        //Now let's refine it some...
        startTime = System.currentTimeMillis();
        
        RealVector ansA = cg.solve(AtA, Atb, ans).getSubVector(0, numTup);//unregularized fit
        System.out.println( "Second-round conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
        
        
        double norm2A = ansA.getNorm();
        double residA = AtA.operate(ansA).getDistance(Atb);
        
        System.out.println( "Second-round fitting norm: " + norm2A + " Normal eq residual: " + residA );
        
        
        //DEBUG!!!!  See if init val helped
        RealVector ansNoReg = cg.solve(AtA, Atb).getSubVector(0, numTup);//unregularized fit
        System.out.println( "No-reg conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
        
        
        double norm2NoReg = ansNoReg.getNorm();
        double residNoReg = AtA.operate(ansNoReg).getDistance(Atb);
        
        System.out.println( "No-reg fitting norm: " + norm2NoReg + " Normal eq residual: " + residNoReg );
        
        
        return ansA.toArray();
    }
}
