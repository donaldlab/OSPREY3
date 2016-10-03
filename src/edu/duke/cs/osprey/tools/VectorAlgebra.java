/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tools;

import java.nio.DoubleBuffer;

/**
 *
 * @author mhall44
 */
public class VectorAlgebra {

    public static final double[] UnitX = { 1, 0, 0 };
    public static final double[] UnitY = { 0, 1, 0 };
    public static final double[] UnitZ = { 0, 0, 1 };

    public static double[] make() {
        return new double[3];
    }
    
    public static double[] make(double[] src) {
        double[] out = make();
        copy(src, out);
        return out;
    }
    
    public static void set(double[] out, double x, double y, double z) {
        out[0] = x;
        out[1] = y;
        out[2] = z;
    }
    
    public static void copy(double[] src, double[] dest) {
        System.arraycopy(src, 0, dest, 0, 3);
    }
   
    public static void copy(double[] src, int srcOffset, double[] dest) {
        System.arraycopy(src, srcOffset, dest, 0, 3);
    }
    
    public static void copy(double[] src, double[] dest, int destOffset) {
        System.arraycopy(src, 0, dest, destOffset, 3);
    }
    
    public static void copy(double[] src, DoubleBuffer dest, int destOffset) {
    	copy(src, 0, dest, destOffset);
    }
    
    public static void copy(double[] src, int srcOffset, DoubleBuffer dest, int destOffset) {
    	dest.position(destOffset);
    	dest.put(src, srcOffset, 3);
    }
    
    public static void copy(DoubleBuffer src, int srcOffset, double[] dest) {
    	copy(src, srcOffset, dest, 0);
    }
    
    public static void copy(DoubleBuffer src, int srcOffset, double[] dest, int destOffset) {
    	src.position(srcOffset);
    	dest[destOffset + 0] = src.get();
    	dest[destOffset + 1] = src.get();
    	dest[destOffset + 2] = src.get();
    }
    
    public static double distance(double[] coords1, int atNum1, double[] coords2, int atNum2){
        //we're given two arrays of coordinates (concatenated 3-D coordinates of atoms)
        //get the distance between the atoms numbered atNum1 in coords1 and atNum2 in coords2
        double dx = coords1[3*atNum1] - coords2[3*atNum2];//difference in x coordinates between the two atoms
        double dy = coords1[3*atNum1+1] - coords2[3*atNum2+1];
        double dz = coords1[3*atNum1+2] - coords2[3*atNum2+2];
        double dist = Math.sqrt( dx*dx + dy*dy + dz*dz );
        
        return dist;
    }
    
    
    
    
    
    //A bunch of 3-D vector operations
    
    //special case of distance for two 3-D vectors
    public static double distance(double[] vec1, double[] vec2){
        return distance(vec1,0,vec2,0);
    }

    public static double dot(double[] vec1, double[] vec2){
        return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    }

    public static double norm(double[] vec){
        return (double)Math.sqrt(dot(vec,vec));
    }

    public static double normsq(double[] vec){
        return dot(vec, vec);
    }

    public static double[] add(double[] vec1,double[] vec2){
        double ans[]=new double[3];
        for(int a=0;a<3;a++){
            ans[a]=vec1[a]+vec2[a];
        }
        return ans;
    }
    
    public static void addInPlace(double[] a, double[] b) {
        a[0] += b[0];
        a[1] += b[1];
        a[2] += b[2];
    }

    public static double[] subtract(double[] x, double[] y){
        double ans[] = new double[x.length];
        for(int a=0; a<x.length; a++)
            ans[a] = x[a]-y[a];

        return ans;
    }
    
    public static void subtractInPlace(double[] a, double[] b) {
       a[0] -= b[0];
       a[1] -= b[1];
       a[2] -= b[2];
    }


    public static double[] scale(double[] vec, double scalar){
        double ans[]=new double[3];
        for(int a=0;a<3;a++){
            ans[a]=vec[a]*scalar;
        }
        return ans;
    }

    public static double[] cross(double[] vec1,double[] vec2){
        double ans[]=new double[3];
        ans[0]=vec1[1]*vec2[2]-vec2[1]*vec1[2];
        ans[1]=vec1[2]*vec2[0]-vec2[2]*vec1[0];
        ans[2]=vec1[0]*vec2[1]-vec2[0]*vec1[1];
        return ans;
    }


    public static double[] parallelComponent(double[] vec1, double[] vec2){//Component of vec1 parallel to vec2
        return scale( vec2, dot(vec1,vec2) / dot(vec2,vec2) );
    }

    public static double[] perpendicularComponent(double[] vec1, double[] vec2){//Component of vec1 perpendicular to vec2
        return subtract(vec1, parallelComponent(vec1,vec2) );
    }


    public static double[] getPerpendicular(double vec[]){//Generate a vector perpendicular to the argument
        if( vec[1]==0 && vec[2] == 0 ){
            double yhat[] = {0,1,0};
            return cross(vec, yhat);
        }
        else{
            double xhat[] = {1,0,0};
            return cross(vec, xhat);
        }
    }

     public static double[] average(double v1[], double v2[]){//Average two vectors
         return scale( add(v1,v2), 0.5f );
     }


    public static double[] normalize(double vec[]){
        return scale( vec, 1.0/norm(vec) );
    }
    
    public static void normalizeInPlace(double v[]) {
        double length = norm(v);
        v[0] /= length;
        v[1] /= length;
        v[2] /= length;
    }


        /**
         * From KiNG's driftwood.r3.Builder:
    * Given three points A, B, and C,
    * construct a line segment from C to D
    * of length len
    * at angle ang to BC (in degrees, 0 to 180)
    * and with a dihedral angle dihe to ABC (in degrees)
    * return D
    * Used in sidechain idealization
    */

    public static double[] get4thPoint(double[] a, double[] b, double[] c, double len, double ang, double dihe)
    {
        double d[] = subtract(b,c);
        d = scale(d, len/norm(d) );

        // Not robust to a/b/c colinear
        // Doesn't matter since that makes dihe undef.
        double x1[] = subtract(a,b);
        double x2[] = subtract(c,b);
        x1 = cross(x1,x2);

        RotationMatrix rot1 = new RotationMatrix(x1[0], x1[1], x1[2], ang, false);
        d = rot1.rotateVector(d);

        rot1 = new RotationMatrix(x2[0], x2[1], x2[2], dihe, false);
        d = rot1.rotateVector(d);

        return add(d,c);
    }

    

    //Given three points of the form (x,y,z),
    //return the coefficients for z in terms of x and y
    //nonsingularity assumed
    //i.e. coord[i][2] = coeffs[0] + coeffs[1]*coord[i][0] + coeffs[2]*coord[i][1] for i=0,1,2
    public static double[] get2DLinCoeffs(double coord[][]){

        double[] coeffs = new double[3];

        double det = (coord[0][1]-coord[1][1])*(coord[1][0]-coord[2][0]) - (coord[1][1]-coord[2][1])*(coord[0][0]-coord[1][0]);
        //determinant of coefficients from eliminating coeffs[0].  Must be nonzero for nonsingular problems


        coeffs[2] = ( (coord[0][2]-coord[1][2])*(coord[1][0]-coord[2][0])
                - (coord[1][2]-coord[2][2])*(coord[0][0]-coord[1][0]) ) / det;

        coeffs[1] = ( (coord[1][2]-coord[2][2])*(coord[0][1]-coord[1][1]) 
                - (coord[0][2]-coord[1][2])*(coord[1][1]-coord[2][1]) ) / det;
        
        coeffs[0] = coord[0][2] - coeffs[1]*coord[0][0] - coeffs[2]*coord[0][1];


        return coeffs;
    }
    
    public static void rotate(double[] inout, double[] x, double[] y, double[] z) {
        double out0 = inout[0]*x[0] + inout[1]*y[0] + inout[2]*z[0];
        double out1 = inout[0]*x[1] + inout[1]*y[1] + inout[2]*z[1];
        double out2 = inout[0]*x[2] + inout[1]*y[2] + inout[2]*z[2];
        inout[0] = out0;
        inout[1] = out1;
        inout[2] = out2;
    }
    
    public static void rotateInverse(double[] inout, double[] x, double[] y, double[] z) {
        double out0 = dot(inout, x);
        double out1 = dot(inout, y);
        double out2 = dot(inout, z);
        inout[0] = out0;
        inout[1] = out1;
        inout[2] = out2;
    }
}
