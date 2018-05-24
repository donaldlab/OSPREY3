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

package edu.duke.cs.osprey.tools;


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.math.Functions;
import com.joptimizer.functions.BarrierFunction;
import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.LogarithmicBarrier;
import com.joptimizer.optimizers.BarrierMethod;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import java.io.Serializable;
import java.util.ArrayList;


//@author mhall44
//This is largely a Java translation of MinVolEllipse.m by Nima Mostagh
public class MinVolEllipse implements Serializable {
    	
    DoubleMatrix2D A;
    DoubleMatrix1D c;//center of ellipse
    DoubleMatrix1D nc;//-c
            
    public MinVolEllipse(DoubleMatrix2D P, double tol, boolean speedup){
        //get an ellipse around the points P
        //MAY WANT TO CHANGE TOLERANCE TO BE BASED ON FURTHEST YOU CAN GO OUTSIDE
        

        //long startTime = System.currentTimeMillis();


        /*%%%%%%%%%%%%%%%%%%%%% Solving the Dual problem%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % ---------------------------------
        % data points 
        % -----------------------------------*/
        int d = P.rows();
        int N = P.columns();


        DoubleMatrix2D Q = DoubleFactory2D.dense.make(d+1,N);
        
        for(int j=0; j<N; j++){
            for(int i=0; i<d; i++)
                Q.set(i,j,P.get(i,j));
            
            Q.set(d,j,1);
        }

        //long qtime = System.currentTimeMillis();
        //System.out.println("Q built at "+(qtime-startTime));


        // initializations
        // -----------------------------------
        int count = 1;
        double err = 1;
        DoubleMatrix1D u = DoubleFactory1D.dense.make(N, 1.0/N);
        //u = (1/N) * ones(N,1);          % 1st iteration


        //long iterTime = System.currentTimeMillis();
        //System.out.println("Iterations starting at "+(iterTime-qtime));

        // Khachiyan Algorithm
        // -----------------------------------
        while(err > tol){

            DoubleMatrix1D M = null;

            if(speedup){
                DoubleMatrix2D X = QuQt(Q,u);

                DoubleMatrix2D invX = null;
                try{
                    invX = Algebra.DEFAULT.inverse(X);
                }
                catch(Exception e){
                    System.err.println("ERROR: Singular matrix in MinVolEllipse calculation");
                }
                
                M = diagMult( Q.zMult(invX,null,1,0,true,false), Q );
            }
            else{
                DoubleMatrix2D udiag = DoubleFactory2D.dense.diagonal(u);
                DoubleMatrix2D Qt = Algebra.DEFAULT.transpose(Q);
                DoubleMatrix2D X = Algebra.DEFAULT.mult(Q, Algebra.DEFAULT.mult(udiag,Qt));
                //X = Q * diag(u) * Q';       % X = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix

                DoubleMatrix2D invX = Algebra.DEFAULT.inverse(X);
                DoubleMatrix2D Mfull = Algebra.DEFAULT.mult(Qt, Algebra.DEFAULT.mult(invX,Q));
                M = DoubleFactory2D.dense.diagonal(Mfull);
            }
            //M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix
            
            int j = -1;
            double maximum = Double.NEGATIVE_INFINITY;
            for(int k=0; k<N; k++){
                if(M.get(k)>maximum){
                    maximum = M.get(k);
                    j = k;
                }
            }
            //[maximum j] = max(M);

            double step_size = (maximum - d -1)/((d+1)*(maximum-1));
            //step_size = (maximum - d -1)/((d+1)*(maximum-1));
            
            DoubleMatrix1D new_u = u.copy().assign(Functions.mult(1-step_size));
            //new_u = (1 - step_size)*u ;
            new_u.set( j, new_u.get(j)+step_size );
            //new_u(j) = new_u(j) + step_size;
            count = count + 1;

            u.assign(Functions.mult(-1));
            u.assign(new_u,Functions.plus);
            err = Math.sqrt( u.zDotProduct(u) ); 
            //err = norm(new_u - u);

            u = new_u;

            //long newIterTime = System.currentTimeMillis();
            //System.out.println("Iteration done at "+(newIterTime-iterTime));
            //iterTime = newIterTime;
        }



        /* %%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
        % Finds the ellipse equation in the 'center form': 
        % (x-c)' * A * (x-c) = 1
        % It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
        % of the ellipse. */

        //DoubleMatrix2D U = DoubleFactory2D.dense.diagonal(u);
        //U = diag(u);

        DoubleMatrix2D Pt = Algebra.DEFAULT.transpose(P);
        c = Algebra.DEFAULT.mult(P,u);     
        A = Algebra.DEFAULT.multOuter(c,c,null);
        A.assign(Functions.mult(-1));
        
        //A.assign( Algebra.DEFAULT.mult(P,Algebra.DEFAULT.mult(U,Pt)), Functions.plus );
        A.assign( QuQt(P,u), Functions.plus );
         
        A = Algebra.DEFAULT.inverse(A);
        A.assign(Functions.mult(1./d));    
        
        /* 
        % the A matrix for the ellipse
        % --------------------------------------------
        A = (1/d) * inv(P * U * P' - (P * u)*(P*u)' );
 
        % center of the ellipse 
        % --------------------------------------------
        c = P * u;*/
        
        nc = c.copy();
        nc.assign(Functions.mult(-1));


        //long doneTime = System.currentTimeMillis();
        //System.out.println("Ellipse built at "+(doneTime-iterTime));
    }
    
    
    
    public MinVolEllipse(DoubleMatrix2D P, double tol, MinVolEllipse[] littleEllipses, MinVolEllipse initEllipse){
        //This version supports enclosing other ellipses as well as points
        //can also consider enclosing voxel level-sets likewise...
        //we take initEllipse (if not null) as an initial value and refine from there
        
        

        //prepare Q as usual
        int d = P.rows();//must also be the dimension of the little ellipses
        int N = P.columns();
        
        if(d==1){//analytically solvable case
            double maxPt = Double.NEGATIVE_INFINITY;
            double minPt = Double.POSITIVE_INFINITY;
            for(int p=0; p<N; p++){
                maxPt = Math.max(maxPt,P.get(0,p));
                minPt = Math.min(minPt,P.get(0,p));
            }
            //ellipse simply stretches from minPt to maxPt
            A = DoubleFactory2D.dense.make( 1, 1, 4./((maxPt-minPt)*(maxPt-minPt)) );
            //(x-c)' * A * (x-c) = 1
            c = DoubleFactory1D.dense.make(1,(maxPt+minPt)/2);
            nc = DoubleFactory1D.dense.make(1,-(maxPt+minPt)/2);
                        
            return;
        }


        DoubleMatrix2D Q = DoubleFactory2D.dense.make(d+1,N);
        
        for(int j=0; j<N; j++){
            for(int i=0; i<d; i++)
                Q.set(i,j,P.get(i,j));
            
            Q.set(d,j,1);
        }
        
        double err = 1;
        DoubleMatrix2D X = null;
        //if(littleEllipses.length==1){
        if(initEllipse!=null){
            //this case is encountered when wrapping bigger ellipses around similar smaller ones...
            //let's initialize with that
            //this is robust to having too few points
            //construct X for the previous ellipse from blocks...
            DoubleMatrix2D blocks[][] = new DoubleMatrix2D[2][2];
            //blocks[0][0] = littleEllipses[0].A.copy();
            blocks[0][0] = initEllipse.A.copy();
            blocks[0][0].assign(Functions.mult(d));
            //DoubleMatrix2D nccol = DoubleFactory2D.dense.make( littleEllipses[0].nc.toArray(), d );//nc as column vector
            DoubleMatrix2D nccol = DoubleFactory2D.dense.make( initEllipse.nc.toArray(), d );//nc as column vector
            blocks[0][1] = blocks[0][0].zMult(nccol, null);
            blocks[1][0] = Algebra.DEFAULT.transpose(blocks[0][1]);
            //double lastBlock = Algebra.DEFAULT.mult( blocks[0][0], littleEllipses[0].c ).zDotProduct(littleEllipses[0].c) + 1;
            double lastBlock = Algebra.DEFAULT.mult( blocks[0][0], initEllipse.c ).zDotProduct(initEllipse.c) + 1;
            blocks[1][1] = DoubleFactory2D.dense.make(1,1,lastBlock);
            
            DoubleMatrix2D invX = DoubleFactory2D.dense.compose(blocks);
            X = invertX(invX);
            
            
            
            //DEBUG!!!
            /*DoubleMatrix2D checkA = invX.viewPart(0,0,d,d).copy();
            checkA.assign(Functions.mult(1./d));
            DoubleMatrix1D checkc = X.viewPart(0,d,d,1).viewColumn(0);
            
            DoubleMatrix2D checkA2 = invX.viewPart(0,0,d,d);
            DoubleMatrix2D invA = invertX(checkA2);
            DoubleMatrix1D lc = invX.viewPart(0,d,d,1).viewColumn(0); 
            DoubleMatrix1D checknc = Algebra.DEFAULT.mult(invA, lc);
            double fac = d+1+lc.zDotProduct(checknc)-invX.get(d,d);
            checkA2.assign(Functions.mult(1./fac));*
            int gg = 0;*/
        }
        else{
            //initialize with just points
            DoubleMatrix1D u = DoubleFactory1D.dense.make(N, 1.0/N);
            X = QuQt(Q,u);
        }
        
        
        while(err > tol){
            //using speedup by default
            //first compute M (matrix of how much points stick out)
            //for points in P
            DoubleMatrix2D invX = invertX(X);

            DoubleMatrix1D M = diagMult( Q.zMult(invX,null,1,0,true,false), Q );
            
            
            DoubleMatrix1D q = null;//direction of maximum sticking out
            double maximum = Double.NEGATIVE_INFINITY;
            for(int k=0; k<N; k++){
                if(M.get(k)>maximum){
                    maximum = M.get(k);
                    q = Q.viewColumn(k);
                }
            }
            
            //now see if any littleEllipses stick out more than the points in P do
            for( MinVolEllipse ell : littleEllipses ){
                DoubleMatrix1D amax = DoubleFactory1D.dense.make(d+1);
                double maxEllNorm = ell.maximizeQuadratic2(invX,amax);
                amax.set(d,1);//homogeneous coordinates as usual
                if(maxEllNorm>maximum){
                    maximum = maxEllNorm;
                    q = amax;
                }
            }
            
            //Looks like the point sticking out the most is q
            double step_size = (maximum - d -1)/((d+1)*(maximum-1));
            //step_size = (maximum - d -1)/((d+1)*(maximum-1));
            
            DoubleMatrix2D newX = X.copy().assign(Functions.mult(1-step_size));
            DoubleMatrix2D dirTerm = Algebra.DEFAULT.multOuter(q, q, null);
            dirTerm.assign(Functions.mult(step_size));
            newX.assign(dirTerm,Functions.plus);
            

            /*X.assign(Functions.mult(-1));
            X.assign(newX,Functions.plus);
            err = Math.sqrt( Algebra.DEFAULT.normF(X) ); */
            //redefining err to be epsilon achieved!!!
            err = maximum/(d+1) - 1;
            //can also use maximum to establish error limit!!

            X = newX;
        }



        /* %%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
        % Finds the ellipse equation in the 'center form': 
        % (x-c)' * A * (x-c) = 1
        % It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
        % of the ellipse. */

        
        
        DoubleMatrix2D invX = invertX(X);
        
        
        //DEBUG!!!!
        /*DoubleMatrix2D B = invX.viewPart(0,0,d,d);
        DoubleMatrix1D bvec = invX.viewPart(0,d,d,1).viewColumn(0);
        double b = invX.get(d,d);
        double checkm[] = new double[N];
        for(int nn=0; nn<N; nn++)
            checkm[nn] = Algebra.DEFAULT.mult(B,P.viewColumn(nn)).zDotProduct(P.viewColumn(nn))
                    + 2*bvec.zDotProduct(P.viewColumn(nn)) + b;
                    * /
        
        
        
        /*A = invX.viewPart(0,0,d,d);
        c = X.viewPart(0,d,d,1).viewColumn(0);        
        A.assign(Functions.mult(1./d));
        nc = c.copy();
        nc.assign(Functions.mult(-1));
        */
        //direct conversion from homogeneous to centered form
        A = invX.viewPart(0,0,d,d);
        DoubleMatrix2D invA = invertX(A);
        DoubleMatrix1D lc = invX.viewPart(0,d,d,1).viewColumn(0); 
        nc = Algebra.DEFAULT.mult(invA, lc);
        double fac = d+1+lc.zDotProduct(nc)-invX.get(d,d);
        A.assign(Functions.mult(1./fac));
        
        c = nc.copy();
        c.assign(Functions.mult(-1));


        //long doneTime = System.currentTimeMillis();
        //System.out.println("Ellipse built at "+(doneTime-iterTime));
    }
    
    
    
 
    
    static DoubleMatrix2D invertX(DoubleMatrix2D X){
            
            if(X.rows()!=X.columns())
                throw new RuntimeException("ERROR: Can't invert square matrix");
            
            if(X.rows()==1){//easy case of 1-D
                if(Math.abs(X.get(0,0))<1e-8)
                    return DoubleFactory2D.dense.make(1,1,0);
                else
                    return DoubleFactory2D.dense.make(1,1,1./X.get(0,0));
            }
        
            /*DoubleMatrix2D invX = null;
        
            try{
                invX = Algebra.DEFAULT.inverse(X);
            }
            catch(Exception e){
                System.err.println("ERROR: Singular matrix in MinVolEllipse calculation");
            }*/
            
            //Invert X, but if singular, this probably indicates a lack of points
            //but we can still make a meaningful lower-dimensional ellipsoid for the dimensions we have data
            //so we construct a pseudoinverse by setting those eigencomponents to 0
            //X is assumed to be symmetric (if not would need SVD instead...)
        
            EigenvalueDecomposition edec = new EigenvalueDecomposition(X);
            
            DoubleMatrix2D newD = edec.getD().copy();
            for(int m=0; m<X.rows(); m++){
                if( Math.abs(newD.get(m,m)) < 1e-8 ){
                    //System.out.println("Warning: singular X in MinVolEllipse calculation");
                    //this may come up somewhat routinely, especially with discrete DOFs in place
                    newD.set(m,m,0);
                }
                else
                    newD.set(m,m,1./newD.get(m,m));
            }

            DoubleMatrix2D invX = Algebra.DEFAULT.mult(edec.getV(),newD).zMult(edec.getV(), null, 1, 0, false, true);
            
            return invX;
    }
    
    
    
    
    double getScaledDistFromCenter(DoubleMatrix1D x){
        //How far are we from the center of the ellipse, relative to the boundary?
        
        if(x==null)
            System.out.println("woah...");
        
        DoubleMatrix1D xrel = x.copy().assign(nc,Functions.plus);
        return xrel.zDotProduct( Algebra.DEFAULT.mult(A, xrel) );
    }
    
    
    double getScaling(DoubleMatrix1D x){
        //give the scaling factor needed to get x on the ellipse
        //(we are scaling x)
        double co1 = x.zDotProduct( Algebra.DEFAULT.mult(A, x) );
        double co2 = -2*x.zDotProduct( Algebra.DEFAULT.mult(A, c) );
        double co3 = c.zDotProduct( Algebra.DEFAULT.mult(A, c) ) - 1;
        
        if(co1==0 && co2==0 ){
            //this means basically scaling x has no effect on whether we're in the ellipse
            
            
            //because we use getScaling as a measure of how "far out" we are relative to the ellipse,
            //if we're outside (or basically on) the ellipse, we'll consider this to be a scaling of 0.  Otherwise infinity.
            
            if(co3>=-1e-10)//outside
                return 0;
            else
                return Double.POSITIVE_INFINITY;
            
            /*if(Math.abs(co3)<1e-10)//we're already basically on the ellipse: scaling x to 0 is sufficient
                return 0;
            
            System.out.print("Bad quadratic in getScaling! co's: "+co1+" "+co2+" "+co3+" A:");
            System.out.print(A);
            System.out.print("x: ");
            System.out.print(x);
            System.out.print("c: ");
            System.out.println(c);*/
        }
        
        double soln1 = quadraticFormula(co1,co2,co3,true,false);
        double soln2 = quadraticFormula(co1,co2,co3,false,false);
        
        //DEBUG!!
        if(Double.isNaN(co1)||Double.isNaN(co2)||Double.isNaN(co3)){
            System.out.print("NaN in getScaling! co's: "+co1+" "+co2+" "+co3+" A:");
            System.out.print(A);
            System.out.print("x: ");
            System.out.print(x);
            System.out.print("c: ");
            System.out.println(c);
        }
        
        
        if(soln1>0 && soln2>0){
            //System.err.println("ERROR: Ellipse doesn't contain origin!");//This may be fairly routine if origin is at edge...
            return Math.max(soln1,soln2);//when we finally get out of ellipse
        }
            
        if( (!(soln1>0)) && (!(soln2>0)) ){
            //System.err.println("ERROR: Ellipse doesn't contain origin!");//and our ray misses it!//May be a small error for small ellipses...
            return 0;
        }
        
        if(soln1>0)
            return soln1;
        else//(soln2>0)
            return soln2;
    }
    
    
    boolean isPointInside(DoubleMatrix1D x){
        //is x inside the ellipse?
        if(getScaledDistFromCenter(x)<1)
            return true;
        
        return false;
    }



    //Sped-up linear algebra operations for the iterations
    //in the constructor
    //These are based on DoubleMatrix2D.zMult (which is used to multiply matrices by DoubleMatrix2D)

    static DoubleMatrix2D QuQt(DoubleMatrix2D Q, DoubleMatrix1D u){
        //Return matrix product of Q * diagonal version of u * Q^T
        //answer = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix

        int m = Q.rows();
        int n = Q.columns();
        if(u.size()!=n){
            throw new Error("Size mismatch in QuQt: "+n+" columns in Q, u length="+u.size());
        }

        DoubleMatrix2D ans = DoubleFactory2D.dense.make(m,m);

        for(int i=0; i<m; i++){
            for(int j=0; j<m; j++){
                double s = 0;
                for(int k=0; k<n; k++)
                    s += Q.getQuick(i,k)*Q.getQuick(j,k)*u.get(k);
                ans.setQuick(i, j, s);
            }
        }

        return ans;
    }


    static DoubleMatrix1D diagMult(DoubleMatrix2D A, DoubleMatrix2D B){
        //Return the diagonal of the product of two matrices
        //A^T and B should have the same dimensions
        int m = A.rows();
        int n = A.columns();
        if(B.rows()!=n || B.columns()!=m){
            throw new Error("Size mismatch in diagMult: A is "+m+"x"+n+
                    ", B is "+B.rows()+"x"+B.columns());
        }

        DoubleMatrix1D ans = DoubleFactory1D.dense.make(m);

        for(int i=0; i<m; i++){
            double s = 0;
            for(int k=0; k<n; k++)
                s += A.getQuick(i, k)*B.getQuick(k, i);
            ans.setQuick(i,s);
        }

        return ans;
    }
    
            /*
        % [A , c] = MinVolEllipse(P, tolerance)
        % Finds the minimum volume enclsing ellipsoid (MVEE) of a set of data
        % points stored in matrix P. The following optimization problem is solved: 
        %
        % minimize       log(det(A))
        % subject to     (P_i - c)' * A * (P_i - c) <= 1
        %                
        % in variables A and c, where P_i is the i-th column of the matrix P. 
        % The solver is based on Khachiyan Algorithm, and the final solution 
        % is different from the optimal value by the pre-spesified amount of 'tolerance'.
        %
        % inputs:
        %---------
        % P : (d x N) dimnesional matrix containing N points in R^d.
        % tolerance : error in the solution with respect to the optimal value.
        %
        % outputs:
        %---------
        % A : (d x d) matrix of the ellipse equation in the 'center form': 
        % (x-c)' * A * (x-c) = 1 
        % c : 'd' dimensional vector as the center of the ellipse. 
        % 
        % example:
        % --------
        %      P = rand(5,100);
        %      [A, c] = MinVolEllipse(P, .01)
        %
        %      To reduce the computation time, work with the boundary points only:
        %      
        %      K = convhulln(P');  
        %      K = unique(K(:));  
        %      Q = P(:,K);
        %      [A, c] = MinVolEllipse(Q, .01)
        %
        %
        % Nima Moshtagh (nima@seas.upenn.edu)
        % University of Pennsylvania
        %
        % December 2005
        % UPDATE: Jan 2009
         * 
         */

   /* 
    double maximizeQuadratic(DoubleMatrix2D augMatrix,DoubleMatrix1D amax){
        //maximize [x^T 1]augMatrix[x;1] for x in this ellipse
        //return highest value, attained at amax (which is in [x;1] form)
        
        //use Joptimizer and Yildirim method
        //we will optimize over the quadratic part of a series, as in ConvexityConstraint etc.
        
        return maximizeQuadratic2(augMatrix,amax);
        
        OptimizationRequest or = new OptimizationRequest();
        
        
        //to get a feasible initial point,
        //we'll need to minimize the second inequality constraint
        int n = augMatrix.rows();
        int dim = SeriesFitter.getNumParamsForOrder(n,2);
        
        
        //inequality constraints
        DoubleMatrix2D hm = homMatrix();
        double ellConstrCoeffs[] = flattenMatrix( hm, true, true );
        LinearMultivariateRealFunction ellConstr = new LinearMultivariateRealFunction(ellConstrCoeffs,hm.get(n-1,n-1));
        
        ConvexityConstraint convConstr = new ConvexityConstraint(n,false,dim,true);
        //no linear terms, just quadratic!
        convConstr.firstQuadCoeff = 0;        
        or.setFi( new ConvexMultivariateRealFunction[] {convConstr} );
        
                
        
        //optimization problem
        or.setToleranceFeas(1.E-12);//as online
        or.setTolerance(1.E-12);
        or.setInteriorPointMethod(JOptimizer.BARRIER_METHOD);
        
        
        or.setF0(ellConstr);
        //with just equality and convexity constraints, can start with identity matrix
        double initPoint1[] = flattenMatrix(DoubleFactory2D.dense.identity(n),true,false);
        or.setInitialPoint( initPoint1 );
        
        JOptimizer opt = new JOptimizer();
        opt.setOptimizationRequest(or);
        
        try{
            int returnCode = opt.optimize();
            System.out.println("JOptimizer return code: "+returnCode);
        } catch(Exception ex){
            throw new Error("maximizeQuadratic JOptimizer error", ex);
        }
        
        double initPoint[] = opt.getOptimizationResponse().getSolution();
        double iv = ellConstr.value(initPoint);
        if( iv >=0 ){
            throw new Error("maximizeQuadratic: no feasible initial point!");
        }
        
        //This gives us an initial point at the edge of the feasible region
        //but since the feasible region is convex, we can step back towards initPoint1
        //to get a strictly feasible point
        //let's use the value iv/2 for ellConstr there (unless iv1 is better than that already)
        //and use the linearity of ellConstr...
        double iv1 = ellConstr.value(initPoint1);
        if(iv/2 < iv1){
            double comp1 = (-iv/2)/(iv1-iv);//that is, (iv/2-iv)/(iv1-iv)
            for(int a=0; a<n; a++)
                initPoint[a] = comp1*initPoint1[a] + (1-comp1)*initPoint[a];
            
            //DEBUG!!
            double ivCheck = ellConstr.value(initPoint);//should be iv/2
            double qqq = 0;
        }
        else
            initPoint = initPoint1;
        
        double convCheck = convConstr.value(initPoint);
        double ellosaurus = ellConstr.value(initPoint);
        or.setInitialPoint(initPoint);
        
        //or.setFi( new ConvexMultivariateRealFunction[] {convConstr,ellConstr} );
        or.setFi(null);
        
        
        //Joptimizer minimizes by default, so negate augMatrix for this purpose
        double objCoeffs[] = flattenMatrix( augMatrix.copy().assign(Functions.mult(-1)), true, true );
        LinearMultivariateRealFunction objf = new LinearMultivariateRealFunction(objCoeffs,-augMatrix.get(n-1,n-1));
        or.setF0(objf);
        
        
        BarrierFunction bf = new LogarithmicBarrier( new 
                    ConvexMultivariateRealFunction[] {convConstr,ellConstr}, dim-1 );
            BarrierMethod bopt = new BarrierMethod(bf);
            bopt.setOptimizationRequest(or);
        //optimization

        try{
            int returnCode = bopt.optimize();
            System.out.println("JOptimizer return code: "+returnCode);
            
            //int returnCode = opt.optimize();
            //System.out.println("JOptimizer return code: "+returnCode);
        } catch(Exception ex){
            throw new Error("JOpitimizer error", ex);
        }
        
        
        double v[] = bopt.getOptimizationResponse().getSolution();//was just opt
        DoubleMatrix2D X = unflattenMatrix(v,n,true,false);
        
        
        //DEBUG!!!
        double detCheck = Algebra.DEFAULT.det(X);
        double ellConstrCheck = ellConstr.value(v);
        
        //eigenize X
        EigenvalueDecomposition edec = new EigenvalueDecomposition(X);
        DoubleMatrix1D evals = edec.getRealEigenvalues();
        DoubleMatrix2D V = edec.getV();
        
        int colTaken = 0;//we take our solution from an arbitrary column
        //but need our vector to have a a nonzero last element...
        //eval<0 is going to be basically 0 because of the convexity constraint
        while( Math.abs(V.viewColumn(colTaken).get(n-1))<1e-6 || evals.get(colTaken)<1e-6 )
            colTaken++;
        DoubleMatrix1D p = V.viewColumn(colTaken).copy();
        p.assign(Functions.mult(Math.sqrt(evals.get(colTaken))));//convert to rank-one decomposition form
        
        double alpha = p.get(n-1);
        
        //finally, copy over solution
        for(int a=0; a<n-1; a++)
            amax.set( a, p.get(a)/alpha );
        
        
        //DEBUG!!
        double objCheck = objf.value(v);
        double objRecheck = -augMatrix.copy().assign(X,Functions.mult).zSum();
        
        //CHECKING...
        for(int col=0; col<V.columns(); col++){
            if( ! ( Math.abs(V.viewColumn(col).get(n-1))<1e-6 || evals.get(col)<1e-6 )  ){
                DoubleMatrix1D p2 = V.viewColumn(col).copy();
                p2.assign(Functions.mult(1./p2.get(n-1)));
                
                DoubleMatrix2D X2 = Algebra.DEFAULT.multOuter(p2,p2,null);
                
                double detCheck2 = Algebra.DEFAULT.det(X2);
                double[] v2 = flattenMatrix(X2,true,false);
                double ellCheck2 = ellConstr.value(v2);
                double objCheck2 = objf.value(v2);
                
                double checkIn2 = hm.zMult(p2,null).zDotProduct(p2);
                double ans2 = augMatrix.zMult(p2, null).zDotProduct(p2);
                
                double qq=0;
            }
        }
        
        
        
        double ans = augMatrix.zMult(p, null).zDotProduct(p) / (alpha*alpha);
        if(Double.isNaN(ans))
            System.out.println("bonster...");
        double checkInside = hm.zMult(p,null).zDotProduct(p);
        
        for(int a=0; a<n-1; a++){
            for(double dx=-0.01; dx<=0.01; dx+=0.02){
                DoubleMatrix1D qp = p.copy();
                qp.set( a, qp.get(a)+dx );
                DoubleMatrix1D qpshort = qp.viewPart(0, qp.size()-1);
                qpshort.assign( Functions.mult( getScaling(qpshort) ) );//scale qpshort to be on ell1
                //ok now qp.get(n-1) should still be 1...

                double checkInp = getScaledDistFromCenter(qpshort);
                double pans = augMatrix.zMult(qp, null).zDotProduct(qp) / (qp.get(n-1)*qp.get(n-1));
                DoubleMatrix2D X2 = Algebra.DEFAULT.multOuter(qp,qp,null);
                X2.assign( DoubleFactory2D.dense.identity(n).assign(Functions.mult(1e-12)), Functions.plus );
                double pansCheck = augMatrix.copy().assign(X2,Functions.mult).zSum() / (qp.get(n-1)*qp.get(n-1));
                double[] v2 = flattenMatrix(X2,true,false);
                if(pans>ans){
                    double ocheck = objf.value(v2);
                    double ellCheck = ellConstr.value(v2);
                    double convvCheck = convConstr.value(v2);
                    
                    //v2 seems to be better than v ... hmm
                    double oocheck = objf.value(v);
                    or.setInitialPoint(v2);
                    
                    try{
                        int returnCode = opt.optimize();
                        System.out.println("JOptimizer return code: "+returnCode);
                    }
                    catch(Exception ex){
                        throw new Error("JOptimizer error", ex);
                    }


                    double v3[] = opt.getOptimizationResponse().getSolution();
                    double ocheck3 = objf.value(v3);
                    
                    double ggg = 0;
                }
            }
        }
        
        
        for(int a=0; a<n-1; a++){
            DoubleMatrix1D pertp = p.copy();
            pertp.set(a, pertp.get(a)/1.01);
            if( hm.zMult(pertp,null).zDotProduct(pertp) < 0 ){//pertp in ellipse
                //if it's exceeding ans we have a problem...
                double check  = augMatrix.zMult(pertp, null).zDotProduct(pertp) / (alpha*alpha);
                if(check>ans)
                    System.out.println("umm");
            }
        }
        
        //objective function at x...
        return augMatrix.zMult(p, null).zDotProduct(p) / (alpha*alpha);
    }
    */
    
    
    DoubleMatrix2D homMatrix(){
        //Return M, where this ellipse is represented as [x^T 1]M[x;1]=0
        int n = c.size();
        DoubleMatrix2D ans = DoubleFactory2D.dense.make(n+1,n+1);
        DoubleMatrix1D Qnc = A.zMult(nc, null);
        for(int a=0; a<n; a++){
            for(int b=0; b<n; b++)
                ans.setQuick( a, b, A.get(a,b) );
            
            ans.setQuick(a, n, Qnc.get(a));
            ans.setQuick(n, a, Qnc.get(a));
        }
        
        ans.setQuick(n, n, Qnc.zDotProduct(nc)-1 );
        
        return ans;
    }
    
    
    
    static double[] flattenMatrix(DoubleMatrix2D M, boolean removeLast, boolean coeffForm){
        //Given a symmetric matrix M
        //represent it in series form as in SeriesFitter, ConvexityConstraint, etc.
        //(and usable for maximization in maximizeQuadratic)
        //if removeLast, the lower right corner of M is assumed to be 1 and removed from the flattened form
        //if coeffForm, we are representing a linear function of a symmetric matrix, 
        //else just the matrix, so to make the functions line up, !coeffForm means double coefficients
        //as necessary so that flattenMatrix(A,removeLast,true) dot flattenMatrix(B,removeLast,false)
        // + (1 if removeLast) = Frobenius product(A,B)
        int count=0;

        int n = M.rows();
        int paramCount = SeriesFitter.getNumParamsForOrder(n, 2);
        if(removeLast)
            paramCount--;
        double ans[] = new double[paramCount];

        for(int a=0; a<n; a++){
            for(int b=0; b<=a; b++){
                if(b==a){
                    ans[count] = M.get(a,a);
                }
                else if(!coeffForm)//M_ab and M_ba collapsed into one series coefficient (using symmetry)
                    ans[count] = M.get(a,b)+M.get(b,a);           
                else
                    ans[count] = M.get(a,b);

                count++;
                if(count==paramCount)//relevant for removeLast
                    break;
            }
        }

        return ans;
    }
    
    
    
    //Reverse conversion (very similar to ConvexityConstraint.makeSeriesHessian)
    static DoubleMatrix2D unflattenMatrix(double serCoeffs[], int nd, boolean removeLast, boolean coeffForm){
        //nd = number of dimensions
        
        int count=0;

        DoubleMatrix2D hess = DoubleFactory2D.dense.make(nd,nd);

        for(int a=0; a<nd; a++){
            for(int b=0; b<=a; b++){
                if(b==a){
                    if(removeLast && a==nd-1)
                        hess.set(a,a,1);
                    else
                        hess.set(a,a,serCoeffs[count]);
                }
                else{
                    double co = serCoeffs[count];
                    if(!coeffForm)
                        co /= 2;
                    hess.set(a,b,co);//Hessian is symmetric
                    hess.set(b,a,co);//but computing by Hessian double-counts off-diagonal series coefficients
                }

                count++;
            }
        }

        return hess;
    }
    
    
    
    
    
    
    double maximizeQuadratic2(DoubleMatrix2D augMatrix,DoubleMatrix1D amax){
        //We are maximizing maximize [x^T 1]augMatrix[x;1] for x in this ellipse
        // max xT M x + gx w/ (x-c)T A (x-c) <= 1
        // If A = C CT, let y = CT(x-c) so x = C^-Ty + c --> maximing y^T (C^-1 M C^-T) y
        // + ( 2 * c^T M C^-T + g C^-T ) y + constant with y^T y <= 1
        //(Note C CT doesn't have to be Cholesky, this is just an efficient way to get it)
        //For this problem see http://www.springerreference.com/docs/html/chapterdbid/72617.html:
        //To minimize 0.5 x^T Q x + cT x over ball of radius r, for non-positive-semidefinite Q,
        //there is a unique sol'n (x,u) w/ necessary & sufficient conditions (1) x is on the ball, 
        //(2) Q+uI is pos. semidef., (3) (Q+uI)x = -c
        //If lambda is the (negative) lowest eigenvalue of Q, u lies between -lambda and -lambda + norm(c)/r
        //so we can bisect on that interval 
        
        
        /*DoubleMatrix2D C = new CholeskyDecomposition(A).getL();
        new CholeskyDecomposition(A).isSymmetricPositiveDefinite();
        
        //DoubleMatrix2D invC = Algebra.DEFAULT.inverse(C);
        //pseudoinverse if C is singular...this will keep us operating on the lower-dimensional subspace of interest
        //calling invertX on A because it's symmetric
        DoubleMatrix2D Ainv = invertX(A);
        DoubleMatrix2D invC = Algebra.DEFAULT.transpose( new CholeskyDecomposition(Ainv).getL() );
        //A = C CT --> A^-1 = C^-T C^-1
        
        if( Double.isNaN(invC.zSum()) )//Error getting C...can result from singular A
            invC = Algebra.DEFAULT.transpose( altCholesky(Ainv) );
        
        
        DoubleMatrix2D checkA = Algebra.DEFAULT.mult(C, Algebra.DEFAULT.transpose(C));
        DoubleMatrix2D checkAinv = Algebra.DEFAULT.mult( Algebra.DEFAULT.transpose(invC), invC );
        DoubleMatrix2D checkiden = Algebra.DEFAULT.mult(C, invC);*/
        
        //Make an A=C C^T decomposition
        //if A is nonsingular just use Cholesky
        CholeskyDecomposition chol = new CholeskyDecomposition(A);
        DoubleMatrix2D C=null, invC=null;
        if(chol.isSymmetricPositiveDefinite()){
            C = chol.getL();
            invC = Algebra.DEFAULT.inverse(C);
        }
        else {//singular...use an eigenvalue-based decomposition
            EigenvalueDecomposition edec = new EigenvalueDecomposition(A);
        
            DoubleMatrix2D sqrtD = edec.getD().copy();
            DoubleMatrix2D sqrtDinv = sqrtD.copy();
            for(int m=0; m<A.rows(); m++){
                double eigVal = sqrtD.get(m,m);
                if(eigVal<0){
                    if(eigVal<-1e-8){
                        System.err.println("ERROR: Non-positive-semidefinite A matrix in ellipse!");
                    }
                    eigVal = 0;
                }
                sqrtD.set( m, m, Math.sqrt(eigVal) );
                if(eigVal>1e-8)
                    sqrtDinv.set(m, m, 1./Math.sqrt(eigVal));
                else
                    sqrtDinv.set(m, m, 0);
            }

            DoubleMatrix2D V = edec.getV();
            C = Algebra.DEFAULT.mult(V,sqrtD);
            invC = Algebra.DEFAULT.mult(sqrtDinv,Algebra.DEFAULT.transpose(V));
        }
        
        
        int n = c.size();
        DoubleMatrix2D M = augMatrix.viewPart(0, 0, n, n);
        
        DoubleMatrix1D g = augMatrix.viewPart(0, n, n, 1).viewColumn(0).copy();
        g.assign(Functions.mult(2));
        
        
        //Set up problem on unit ball
        //note negation from the maximization form above (and doubling of Hessian)
        DoubleMatrix2D Qball = Algebra.DEFAULT.mult( invC, M.zMult(invC,null,1,0,false,true) );
        Qball.assign(Functions.mult(-2));
        
        DoubleMatrix1D cball =  Algebra.DEFAULT.mult( Algebra.DEFAULT.mult(invC,M), c );
        cball.assign(Functions.mult(-2));
        cball.assign( Algebra.DEFAULT.mult(invC,g), Functions.minus );
        
        double lambda = 0;//lowest eigenvalue of Qball
        EigenvalueDecomposition edec = new EigenvalueDecomposition(Qball);
        DoubleMatrix1D eigVals = edec.getRealEigenvalues();
        DoubleMatrix1D lvec = null;
        for(int m=0; m<n; m++){
            if(eigVals.get(m)<lambda){
                lambda = eigVals.get(m);
                lvec = edec.getV().viewColumn(m);
            }
        }
        
        if(lambda >= 0){
            //this is a convex problem, not what this class is intended for
            System.err.println("ERROR: maximizeQuadratic2 used with convex obj function...");
        }
        
        
        DoubleMatrix1D y = null;
        
        
        DoubleMatrix2D ncball = DoubleFactory2D.dense.make( cball.toArray(), n );//nc as column vector
        ncball.assign(Functions.mult(-1));

        boolean useBisection = true;
        
        DoubleMatrix2D qi = null;
        DoubleMatrix1D x0 = null;//FOR DEBUGGING: WANT TO SEE THESE
        double discr = -42;//special value to show lack of tampering

        
        
        if( Algebra.DEFAULT.norm2(cball) < 1e-8 ){
            //singularity: in this case to satisfy the conditions
            //we must have y in the nullspace of Q+uI where u=-lambda
            //i.e. y is on the ball and is an eigenvector of Q with the lowest eigenvalue
            y = lvec;
            useBisection = false;
        }
        /*else if( Math.abs( cball.zDotProduct(lvec) - Math.sqrt(Algebra.DEFAULT.norm2(cball)) ) < 1e-6 ){//cball in lvec direction
            //this would give us a singularity in our bisection method
            //but since lambda is the lowest eigenvalue and its eigenvector is aligned with cball,
            //then we know our optimum is in ncball direction
            y = lvec.copy();
            if(cball.zDotProduct(lvec)>0)
                lvec.assign(Functions.mult(-1));
        }*/        
        else if ( Math.abs(cball.zDotProduct(lvec)) < 1e-8 ){
            //(Qball-lambda*I)*x0 = -cball
            //Get inverse of Qball-lambda*I with lvec component removed
            DoubleMatrix2D newD = edec.getD().copy();
            for(int m=0; m<n; m++){
                if( Math.abs(newD.get(m,m)-lambda) < 1e-8 )
                    newD.set(m,m,0);
                else
                    newD.set(m,m,1./(newD.get(m,m)-lambda));
            }
            qi = Algebra.DEFAULT.mult(edec.getV(),newD).zMult(edec.getV(), null, 1, 0, false, true);
            x0 = Algebra.DEFAULT.mult(qi,ncball).viewColumn(0);
            
            //now we can add a multiple of lvec to x0 to get y
            //we'll still have (Qball-lambda*I)*y = -cball
            //but now we'll get y on the ball.  So our necessary and sufficient conditions will be satisfied.  
            //if there is no real solution to the quadratic,
            //then we have no solution on the ball, so we need to use bisection (move to other Qball+uI)
            //to find the correct solution.
            //Either way we find the unique solution (u,y) to the conditions

            double qb = 2*lvec.zDotProduct(x0);
            double qc = x0.zDotProduct(x0) - 1;
            double a1 = quadraticFormula( 1, qb, qc, true, false);
            
            if(!Double.isNaN(a1)){
                double a2 = quadraticFormula( 1, qb, qc, false, true);
                double a = a1;
                if(a2*a2>a1*a1)
                    a = a2;
                //since cball is perpendicular to lvec, adding multiples of lvec won't affect objective function --> maximize a^2

                y = lvec.copy();
                y.assign(Functions.mult(a));
                y.assign(x0,Functions.plus);

                useBisection = false;
            }
            else{
                //DEBUG!!!
                /*discr = qb*qb - 4*qc;
                
                DoubleMatrix2D Qshifted = Qball.copy().assign( 
                        DoubleFactory2D.dense.diagonal(DoubleFactory1D.dense.make(n,1e-10-lambda)), Functions.plus );
                y = Algebra.DEFAULT.solve( Qshifted, ncball ).viewColumn(0);//y should be basically x0!
                double check = y.zDotProduct(y);
                DoubleMatrix1D checkncball = Algebra.DEFAULT.mult(Qshifted, y);
                DoubleMatrix1D checkncball0 = Algebra.DEFAULT.mult(Qshifted, x0);
                
                
                double gggg = 0;*/
                
            }

            //DEBUG!!
            /*DoubleMatrix2D Qshifted = Qball.copy().assign( 
                        DoubleFactory2D.dense.diagonal(DoubleFactory1D.dense.make(n,-lambda)), Functions.plus );
            DoubleMatrix1D nccheck = Qshifted.zMult(y, null);
            double normcheck = y.zDotProduct(y);
            double qqq = 0;*/
        }
        
        if(useBisection) {
        
            //bisect

            //wait probably don't need to bisect if we have the eigendecomposition...
            //(1) y is on the ball, 
            //(2) Q+uI is pos. semidef., (3) (Q+uI)y = -c
            //in orthonormal eigenbasis of Q, (Q+uI)^-1 * c is a diagonal transformation
            // Q = S lam ST --> Q + uI = S (lam + uI) St --> (Q+uI)^-1 = S(lam+uI)^-1 St
            // y = - S (lam+uI)^-1 Stc, so y dot y = cT S (lam+uI)^-2 St c = 1
            // let q = St c, then we have 1 = sum_i q_i^2 / (lam_i + u)^2
            //hmm ok maybe do bisect haha

            double top = Math.sqrt(Algebra.DEFAULT.norm2(cball))-lambda;//should give solution inside sphere
            double bottom = -lambda;//should be outside...


            double bisecPrecision = 1e-10;


            //DEBUG!!!
            top *= 1+1e-10;//some tolerance...
            for(double mid : new double[] {top,bottom*(1+1e-10)} ){//shift bottom up to avoid singularity
                //make sure top gives y inside sphere and bottom outside!
                DoubleMatrix2D Qshifted = Qball.copy().assign( 
                        DoubleFactory2D.dense.diagonal(DoubleFactory1D.dense.make(n,mid)), Functions.plus );
                y = Algebra.DEFAULT.solve( Qshifted, ncball ).viewColumn(0);//FIX!!!
                double check = y.zDotProduct(y);
                if( (check>1) == (mid==top) ){
                    System.err.println("ERROR IN maximizeQuadratic2 bisection: root not bracketed!");
                    
                    DoubleMatrix1D checkncball = Algebra.DEFAULT.mult(Qshifted, y);
                    DoubleMatrix1D checkncball0 = Algebra.DEFAULT.mult(Qshifted, x0);
                    double ggg = 0;
                    
                }
            }

            double mid = 0;

            while( top-bottom > bisecPrecision*Math.max(Math.abs(top),Math.abs(bottom)) ){
                mid = (top+bottom)/2;

                DoubleMatrix2D Qshifted = Qball.copy().assign( 
                        DoubleFactory2D.dense.diagonal(DoubleFactory1D.dense.make(n,mid)), Functions.plus );
                y = Algebra.DEFAULT.solve( Qshifted, ncball ).viewColumn(0);

                if( y.zDotProduct(y) > 1 )
                    bottom = mid;
                else
                    top = mid;
            }
        }

        
        DoubleMatrix1D x = c.copy();
        x = invC.zMult(y, x, 1, 1, true);
        
        for(int a=0; a<n; a++)
            amax.set(a, x.get(a));
        
        amax.set(n, 1);
        return augMatrix.zMult(amax, null).zDotProduct(amax);
        //return highest value, attained at amax (which is in [x;1] form)
    }
    
    
    
    
    static double quadraticFormula (double a, double b, double c, boolean firstRoot, 
            boolean errorOnComplex) {
        //firstRoot indicates which root to use
        //for negative discriminants, we call an error if errorOnComplex; otherwise we flag
        //the complex value by returning nan
        double discr = b*b-4*a*c;

        if(a==0){//probably need to deal with linear case explicitly...
            if(b==0){//c==0
                throw new RuntimeException("ERROR: Can't solve equation "+c+"=0 by quadratic formula");
            }
            else//bx+c=0
                return -c/b;
        }

        if(discr<0){
            if(errorOnComplex){
                throw new RuntimeException("ERROR: Negative discrimant in quadratic formula");
            }
            else
                return Double.NaN;
        }

        if(b>0){
            if(firstRoot)
                return (-b-Math.sqrt(discr))/(2*a);
            else
                return 2*c/(-b-Math.sqrt(discr));
        }
        else{
            if(firstRoot)
                return (-b+Math.sqrt(discr))/(2*a);
            else
                return 2*c/(-b+Math.sqrt(discr));
        }
    }
    
    public DoubleMatrix2D getA() { return A.copy(); }
    public DoubleMatrix1D getC() { return c.copy(); }
    public DoubleMatrix1D getNC() { return nc.copy(); }
    
}
