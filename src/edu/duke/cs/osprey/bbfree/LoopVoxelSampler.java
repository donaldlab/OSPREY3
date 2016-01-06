/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import edu.duke.cs.osprey.dof.deeper.perts.TripeptideClosure;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * Sample conformations of a loop within a voxel, CA to CA flexible)
 * 
 * @author mhall44
 */
public class LoopVoxelSampler {
    
    /*
    int numRes;//number of residues (including anchor residues) in the loop
    double[][] NCoord, CACoord, CCoord;//their main BB atom coords
    TripeptideClosure tc;//closure for first tripeptide (has bondlen, etc. data from orig conf)
    
    
    double voxelConstr[][];//voxel constraints on N and CA coords
    //indices: 0/1, coord number (over all freely moving N's and CA's)
    
    //Data (for each res) on circles where we'll find N's and CA's relative to other coords
    double NCircleRad[], NCircleHeight[], CACircleRad[], CACircleHeight[];
    
    PepPlaneLinModel[] pepPlanes;
    //CONSTRUCTOR TO BE CALLED AT ORIG CONF SO CAN INIT ALL THESE THINGS
    
    
    public LoopVoxelSampler(List<Residue> residues, double constr[][], PepPlaneLinModel[] planes){
        numRes = residues.size();
        pepPlanes = planes;
        voxelConstr = constr;
        
        NCoord = new double[numRes][];
        CACoord = new double[numRes][];
        CCoord = new double[numRes][];
        
        //record anchor coords (the rest will be sampled)...
        NCoord[0] = residues.get(0).getCoordsByAtomName("N");
        CACoord[0] = residues.get(0).getCoordsByAtomName("CA");
        CACoord[numRes-1] = residues.get(numRes-1).getCoordsByAtomName("CA");
        CCoord[numRes-1] = residues.get(numRes-1).getCoordsByAtomName("C");
        
        
        tc = new TripeptideClosure( residues.subList(0,3) );
        
        NCircleRad = new double[numRes];
        NCircleHeight = new double[numRes];
        CACircleRad = new double[numRes];
        CACircleHeight = new double[numRes];
        
        for(int resNum=0; resNum<numRes; resNum++){
            
            double N[] = residues.get(resNum).getCoordsByAtomName("N");
            double CA[] = residues.get(resNum).getCoordsByAtomName("CA");
            double C[] = residues.get(resNum).getCoordsByAtomName("C");
            
            NCircleRad[resNum] = measureCircleRad(N,CA,C);
            NCircleHeight[resNum] = measureCircleHeight(N,CA,C);
            
            if(resNum<numRes-1){//there's a next residue, to measure CA circle in terms of
                double nextN[] = residues.get(resNum+1).getCoordsByAtomName("N");
                double nextCA[] = residues.get(resNum+1).getCoordsByAtomName("CA");
                
                CACircleRad[resNum] = measureCircleRad(CA,nextN,nextCA);
                CACircleHeight[resNum] = measureCircleHeight(CA,nextN,nextCA);
            }
        }
    }
    
    
    
    DoubleMatrix2D sampleFullDOFs(int numSamples){
        //Sample valid confs, in full-DOF format, in voxel (for N's and CA's)
        //Return indices: DOF, sample
        
        
        double[][] ans = new double[6*numRes-9][numSamples];
        
        
        int sampleCount = 0;
        
        
        //OK so we're trying to fill in a CA-to-CA loop
        //we'll sample peptide planes backwards from the end,
        //then close the initial tripeptide
        
        int TCStartRes = 0;
        int TCEndRes = 2;
        
        //going to sample backwards
        
        while(sampleCount<numSamples){//still need to draw samples...
            

            
            /*for(int resNum=0; resNum<=TCStartRes; resNum++)
                sampleResBBForwards(resNum);*/ // sampling from end...
     /*       
            boolean sampledPep = true;
            
            for(int resNum=numRes-2; resNum>=TCEndRes; resNum--){
                if( ! samplePepPlaneBackwards(resNum) ){//got stuck sampling...start over
                    //Note: because the orig conf is in the middle of the voxel,
                    //there are always samples we can draw that won't get stuck
                    //if we keep trying
                    sampledPep = false;
                    break;
                }
            }
            
            if(!sampledPep)
                continue;
            
            //OK now pep planes have been sampled except for the first two
            //close them
            
            double r_soln_n[][][] = new double[16][3][3];
            double r_soln_a[][][] = new double[16][3][3];
            double r_soln_c[][][] = new double[16][3][3];
            
            int numSoln = tc.solve_3pep_poly( NCoord[TCStartRes], CACoord[TCStartRes], 
                    CACoord[TCEndRes], CCoord[TCEndRes], r_soln_n, r_soln_a, r_soln_c);
            
            
            boolean foundGoodSoln = false;
            
            //see if there's a soln satisfying the voxel constraints...
            //there are constraints on the middle N and C, and the last N
            for(int soln=0; soln<numSoln; soln++){
                
                
                if( checkNVoxelConstr(r_soln_n[soln][1], TCStartRes+1)
                        && checkCAVoxelConstr(r_soln_a[soln][1], TCStartRes+1)
                        && checkNVoxelConstr(r_soln_n[soln][2], TCEndRes) ){
                    
                    foundGoodSoln = true;
                    
                    //record the N and CA coord values...
                    NCoord[TCStartRes+1] = r_soln_n[soln][1];
                    CACoord[TCStartRes+1] = r_soln_a[soln][1];
                    NCoord[TCEndRes] = r_soln_n[soln][2];
                    
                    break;
                }
            }
            
            if(foundGoodSoln){
                //record the solution
                for(int resNum=1; resNum<numRes; resNum++){
                    ans[6*resNum-6][sampleCount] = NCoord[resNum][0];
                    ans[6*resNum-5][sampleCount] = NCoord[resNum][1];
                    ans[6*resNum-4][sampleCount] = NCoord[resNum][2];
                    
                    if(resNum<numRes-1){
                        ans[6*resNum-3][sampleCount] = CACoord[resNum][0];
                        ans[6*resNum-2][sampleCount] = CACoord[resNum][1];
                        ans[6*resNum-1][sampleCount] = CACoord[resNum][2];
                    }
                }
                
                sampleCount++;
            }
        }
        
        return DoubleFactory2D.dense.make(ans);
    }
    
    
    
    private boolean samplePepPlaneBackwards(int resNum){
        //Given the coordinates of the main BB chain atoms
        //from CA of resNum+1 onward,
        //sample coords for resNum's CA and C' and resNum+1's N
        //making sure that they're in our voxel (for N and CA)
        //Rejection sampling to stay in voxel
        //Return false if try a bunch of time and keep failing to get good ones
        //(We'll start from the beginning then)
       
               
        int maxNumTries = 50;//going to try this many samples
        //if it fails then we'll go back and sample all DOFs,
        //for we could be shoved out of the voxel by the previous res' samples        
        
        for (int tryNum=0; tryNum<maxNumTries; tryNum++) {
            //start by sampling N, constrained to voxel and to circle (only psi can rotate)...
            NCoord[resNum+1] = sampleOnCircle(NCircleRad[resNum+1], NCircleHeight[resNum+1],
                    CACoord[resNum+1], CCoord[resNum+1]);

            if( ! checkNVoxelConstr(NCoord[resNum+1], resNum+1) )
                continue;//failed...redraw (will overwrite NCoord[resNum+1])
            
            CACoord[resNum] = sampleOnCircle(CACircleRad[resNum], CACircleHeight[resNum],
                    NCoord[resNum+1], CACoord[resNum+1]);
            
            if( checkCAVoxelConstr(CACoord[resNum], resNum) ){
                CCoord[resNum] = pepPlanes[resNum].calcCCoords(CACoord[resNum],NCoord[resNum+1],CACoord[resNum+1]);
                return true;//looks good!  And our coords are in place
            }
        }
        
        //If we get here, we have failed maxNumTries times
        return false;
    }
    
    
    //Measure rad and height such that coord is on the circle with radius rad
    //around the point b+height(b-a)
    private double measureCircleHeight(double coord[], double b[], double a[]){
        double axis[] = VectorAlgebra.subtract(b,a);
        double relCoord[] = VectorAlgebra.subtract(coord, b);
        double upVector[] = VectorAlgebra.parallelComponent(relCoord, axis);
        return VectorAlgebra.norm(upVector) / VectorAlgebra.norm(axis);
    }
    
    private double measureCircleRad(double coord[], double b[], double a[]){
        double axis[] = VectorAlgebra.subtract(b,a);
        double relCoord[] = VectorAlgebra.subtract(coord, b);
        double radVector[] = VectorAlgebra.perpendicularComponent(relCoord, axis);
        return VectorAlgebra.norm(radVector);
    }
    
    
    private double[] sampleOnCircle(double rad, double height, double[] b, double[] a){
        //Sample a point uniformly a point on the circle with radius rad
        //around the point b+height(b-a)
        
        double axis[] = VectorAlgebra.subtract(b,a);
        double circleCenter[] = VectorAlgebra.add(b, VectorAlgebra.scale(axis,height) );
        
        double[] ax2 = VectorAlgebra.getPerpendicular(axis);
        double[] ax3 = VectorAlgebra.cross(ax2, axis);
        
        double angle = 2 * Math.random() * Math.PI;

        
        ax2 = VectorAlgebra.scale( ax2, Math.cos(angle) * rad / VectorAlgebra.norm(ax2) );
        ax3 = VectorAlgebra.scale( ax3, Math.sin(angle) * rad / VectorAlgebra.norm(ax3) );
        
        
        double[] ans = VectorAlgebra.add( circleCenter, VectorAlgebra.add(ax2,ax3) );
        
        //DEBUG!!!
        //double checkr = measureCircleRad(ans,b,a);
        //double checkh = measureCircleHeight(ans,b,a);
        
        
        return ans;
    }
    
    
    private boolean checkNVoxelConstr(double coord[], int resNum){
        //Are these coords in range to be coordinate of this residue's nitrogen?
        for(int dim=0; dim<3; dim++){
            if( coord[dim]>voxelConstr[1][6*(resNum-1)+dim]
                    || coord[dim]<voxelConstr[0][6*(resNum-1)+dim] ){
                return false;//out of range!
            }
        }
        
        return true;
    }
    
    
    private boolean checkCAVoxelConstr(double coord[], int resNum){
        //Are these coords in range to be coordinate of this residue's alpha carbon?
        for(int dim=0; dim<3; dim++){
            if( coord[dim]>voxelConstr[1][6*(resNum-1)+3+dim]
                    || coord[dim]<voxelConstr[0][6*(resNum-1)+3+dim] ){
                return false;//out of range!
            }
        }
        
        return true;
    }*/
}
