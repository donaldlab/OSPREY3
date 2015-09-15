/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper.perts;


/* Modified version of the code downloaded
 * from http://www.acm.org/pubs/tog/GraphicsGems/gems/Sturm/
 * Modified by Chaok Seok 2003
     * (that version came with the paper
     * "A Kinematic View of Loop Closure",  E. A. Coutsias, C. Seok,
M. P. Jacobson, and K. A. Dill, Journal of Computational Chemistry,
Volume 25, Issue 4, Pages 510 - 528 (2004)
     * )
 * Then translated from C to Java for use in OSPREY by MH
 *
 * Using Sturm Sequences to Bracket Real Roots of Polynomial Equations
 * by D.G. Hook and P.R. McAree
 * from "Graphics Gems", Academic Press, 1990
*/

import java.io.Serializable;

public class SturmSolver implements Serializable{


    final int PRINT_LEVEL=0;
    final int MAX_ORDER=16;
    final int MAXPOW=32;
    final double SMALL_ENOUGH=1.0e-18;

    double RELERROR;
    int MAXIT, MAX_ITER_SECANT;


    private class poly {
        int	ord;
        double[] coef = new double[MAX_ORDER+1];
    }

    /* set termination criteria for polynomial solver */
    //In C this was a function called initialize_sturm
    SturmSolver(double tol_secant, int max_iter_sturm, int max_iter_secant)
    {
      RELERROR = tol_secant;
      MAXIT = max_iter_sturm;
      MAX_ITER_SECANT = max_iter_secant;
    }

    
    //void solve_sturm(int *p_order, int *n_root, double poly_coeffs[], double roots[])
    public int solve_sturm(int order, double poly_coeffs[], double roots[])
    //Returns the number of roots
    {
      poly sseq[] = new poly[MAX_ORDER*2];
      for(int a=0; a<MAX_ORDER*2; a++)
          sseq[a] = new poly();

      double min, max;
      int i, j, nroots, nchanges, np, atmin = 0, atmax = 0;

      for (i = order; i >= 0; i--)
        {
          sseq[0].coef[i] = poly_coeffs[i];
        }

      if (PRINT_LEVEL > 0)
        {
          for (i = order; i >= 0; i--)
            {
              System.out.printf("coefficients in Sturm solver\n");
              System.out.printf("%d %lf\n", i, sseq[0].coef[i]);
            }
        }

      /*
       * build the Sturm sequence
       */
      np = buildsturm(order, sseq);

      if (PRINT_LEVEL > 0)
        {
          System.out.printf("Sturm sequence for:\n");
          for (i = order; i >= 0; i--)
            System.out.printf("%lf ", sseq[0].coef[i]);
          System.out.printf("\n\n");
          for (i = 0; i <= np; i++)
            {
              for (j = sseq[i].ord; j >= 0; j--)
                System.out.printf("%lf ", sseq[i].coef[j]);
              System.out.printf("\n");
            }
          System.out.printf("\n");
        }

      /*
       * get the number of real roots
       */

      int at[] = {atmin, atmax};
      nroots = numroots(np, sseq, at);
      atmin = at[0];
      atmax = at[1];


      if (nroots == 0)
          return 0;
      // printf("solve: no real roots\n");


      if (PRINT_LEVEL > 0)
        System.out.printf("Number of real roots: %d\n", nroots);

      /*
       * calculate the bracket that the roots live in
       */
      min = -1;
      nchanges = numchanges(np, sseq, min);

      for (i = 0; nchanges != atmin && i != MAXPOW; i++) {
        min *= 10.0;
        nchanges = numchanges(np, sseq, min);
      }

      if (nchanges != atmin) {
        System.out.printf("Sturm solver: unable to bracket all negative roots\n");
        atmin = nchanges;
      }

      max = 1;
      nchanges = numchanges(np, sseq, max);
      for (i = 0; nchanges != atmax && i != MAXPOW; i++) {
        max *= 10.0;
        nchanges = numchanges(np, sseq, max);
      }

      if (nchanges != atmax) {
        System.out.printf("Sturm solver: unable to bracket all positive roots\n");
        atmax = nchanges;
      }

      nroots = atmin - atmax;

      /*
       * perform the bisection.
       */

      sbisect(np, sseq, min, max, atmin, atmax, roots);

      /*
       * write out the roots...
       */
      if (PRINT_LEVEL > 0)
        {
          if (nroots == 1) {
            System.out.printf("\n1 distinct real root at x = %f\n", roots[0]);
          } else {
            System.out.printf("\n%d distinct real roots for x: \n", nroots);

            for (i = 0; i != nroots; i++)
              {
                System.out.printf("%f\n", roots[i]);
              }
          }
        }

      return nroots;
    }

   

    double hyper_tan(double a, double x)
    {
      double exp_x1, exp_x2, ax;

      ax = a*x;
      exp_x1 = Math.exp(ax);
      exp_x2 = Math.exp(-ax);
      
      if( exp_x1 == Double.POSITIVE_INFINITY )
          return 1;
      else if( exp_x2 == Double.POSITIVE_INFINITY )
          return -1;
      else
          return (exp_x1 - exp_x2)/(exp_x1 + exp_x2);
    }


     /*
     * modp
     *
     *	calculates the modulus of u(x) / v(x) leaving it in r, it
     *  returns 0 if r(x) is a constant.
     *  note: this function assumes the leading coefficient of v
     *	is 1 or -1
     */
    public int modp(poly u, poly v, poly r)
    {
            int		k, j;

            int nr = 0;
            int end = u.ord;
            int uc = 0;

            while( uc<=end )
                r.coef[nr++] = u.coef[uc++];

            
            if (v.coef[v.ord] < 0) {


                            for (k = u.ord - v.ord - 1; k >= 0; k -= 2)
                                    r.coef[k] = -r.coef[k];

                            for (k = u.ord - v.ord; k >= 0; k--)
                                    for (j = v.ord + k - 1; j >= k; j--)
                                            r.coef[j] = -r.coef[j] - r.coef[v.ord + k]
                                            * v.coef[j - k];
            } else {
                            for (k = u.ord - v.ord; k >= 0; k--)
                                    for (j = v.ord + k - 1; j >= k; j--)
                                        r.coef[j] -= r.coef[v.ord + k] * v.coef[j - k];
            }

            k = v.ord - 1;
            while (k >= 0 && Math.abs(r.coef[k]) < SMALL_ENOUGH) {
                    r.coef[k] = 0;
                    k--;
            }

            r.ord = (k < 0) ? 0 : k;

            return(r.ord);
    }

    /*
     * buildsturm
     *
     *	build up a sturm sequence for a polynomial in smat, returning
     * the number of polynomials in the sequence
     */
    int buildsturm(int ord, poly sseq[])
    {
            int		i;
            double	f;
            int sp;

            sseq[0].ord = ord;
            sseq[1].ord = ord - 1;

            /*
             * calculate the derivative and normalise the leading
             * coefficient.
             */
            f = Math.abs(sseq[0].coef[ord]*ord);


            int fp = 0;
            int fc = 1;

            for (i=1; i<=ord; i++)
                sseq[1].coef[fp++] = sseq[0].coef[fc++] * i / f;


            // construct the rest of the Sturm sequence

            for ( sp = 2; modp( sseq[sp-2], sseq[sp-1], sseq[sp] ) != 0 ; sp++ ){

                //reverse the sign and normalise
                f = -Math.abs(sseq[sp].coef[sseq[sp].ord]);
                for( fp=sseq[sp].ord; fp>=0; fp-- )
                    sseq[sp].coef[fp] /= f;

            }

            sseq[sp].coef[0] = - sseq[sp].coef[0];//reverse the sign

            return sp;
    }

    /*
     * numroots
     *
     *	return the number of distinct real roots of the polynomial
     * described in sseq.
     */
    int numroots(int np, poly sseq[], int at[])
    {
                    int	atposinf, atneginf, s;
                    double	f, lf;

                    atposinf = atneginf = 0;


            /*
             * changes at positive infinity
             */
            lf = sseq[0].coef[sseq[0].ord];

            for ( s = 1; s <= np; s++ ){
                            f = sseq[s].coef[sseq[s].ord];
                            if (lf == 0.0 || lf * f < 0)
                                    atposinf++;
                    lf = f;
            }

            /*
             * changes at negative infinity
             */
            if ( (sseq[0].ord & 1) != 0 )
                            lf = -sseq[0].coef[sseq[0].ord];
            else
                            lf = sseq[0].coef[sseq[0].ord];

            for (s = 1; s <= np; s++ ) {
                            if ( ( sseq[s].ord & 1 ) != 0 )
                                    f = -sseq[s].coef[sseq[s].ord];
                            else
                                    f = sseq[s].coef[sseq[s].ord];
                            if (lf == 0 || lf * f < 0)
                                    atneginf++;
                            lf = f;
            }

            at[0] = atneginf;
            at[1] = atposinf;
            //	printf("atneginf, atposinf = %d %d\n", atneginf, atposinf);
            return(atneginf - atposinf);
    }

    /*
     * numchanges
     *
     *	return the number of sign changes in the Sturm sequence in
     * sseq at the value a.
     */
    int numchanges(int np, poly sseq[], double a)
    //	int		np;
    //	poly	*sseq;
    //	double	a;

    {
            int	changes;
            double f, lf;
            int s;

            changes = 0;

            lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

            for (s = 1; s <= np; s++) {
                            f = evalpoly(sseq[s].ord, sseq[s].coef, a);
                            if (lf == 0.0 || lf * f < 0)
                                    changes++;
                            lf = f;
    //			printf("lf %lf %d \n", f, changes);
            }

            //	printf("%d \n", changes);
            return(changes);
    }

    /*
     * sbisect
     *
     *	uses a bisection based on the sturm sequence for the polynomial
     * described in sseq to isolate intervals in which roots occur,
     * the roots are returned in the roots array in order of magnitude.
     */
    void sbisect(int np, poly sseq[], double min, double max, int atmin, int atmax, double roots[])
    {
      double mid=0;
      int n1 = 0, n2 = 0, its, atmid, nroot = atmin - atmax;


      if (nroot == 1)
        {

          /*
           * first try a less expensive technique.
           */
          //      printf("min max %lf, %lf \n", min, max);
          if (modrf(sseq[0].ord, sseq[0].coef, min, max, roots))
            return;

          //      printf("try hard way\n");
          /*
           * if we get here we have to evaluate the root the hard
           * way by using the Sturm sequence.
           */
          for (its = 0; its < MAXIT; its++)
            {
              mid = (min + max) / 2;

              atmid = numchanges(np, sseq, mid);

              if ( Math.abs(mid) > RELERROR)
                {
                  if ( Math.abs((max - min) / mid) < RELERROR)
                    {
                      roots[0] = mid;
                      return;
                    }
                }
              else if (Math.abs(max - min) < RELERROR)
                {
                  roots[0] = mid;
                  return;
                }

              if ((atmin - atmid) == 0)
                min = mid;
              else
                max = mid;
            }

          if (its == MAXIT)
            {
              System.err.printf("sbisect: overflow min %f max %f" +
                                            "diff %e nroot %d n1 %d n2 %d\n",
                      min, max, max - min, nroot, n1, n2);
              roots[0] = mid;
            }
          return;
        }

      /*
       * more than one root in the interval, we have to bisect...
       */

      for (its = 0; its < MAXIT; its++)
        {

          mid = (min + max) / 2;

          atmid = numchanges(np, sseq, mid);

          n1 = atmin - atmid;
          n2 = atmid - atmax;

          if (n1 != 0 && n2 != 0)
            {
              sbisect(np, sseq, min, mid, atmin, atmid, roots);

              double rootsn1[] = new double[roots.length - n1];//Root #n1 through the last one
              sbisect(np, sseq, mid, max, atmid, atmax, rootsn1);
              System.arraycopy(rootsn1, 0, roots, n1, rootsn1.length);
              
              break;
            }

          if (n1 == 0)
            min = mid;
          else
            max = mid;
        }

      if (its == MAXIT)
        {
          /*
          fprintf(stderr, "sbisect: roots too close together\n");
          fprintf(stderr, "sbisect: overflow min %f max %f diff %e\
                                    nroot %d n1 %d n2 %d\n",
                  min, max, max - min, nroot, n1, n2);
          */
          for (n1 = atmax; n1 < atmin; n1++)
            roots[n1 - atmax] = mid;
        }
    }

    /*
     * evalpoly
     *
     *	evaluate polynomial defined in coef returning its value.
     */
    double evalpoly(int ord, double coef[], double x)
    {
            double f = coef[ord];

            for( int fp=ord-1; fp>=0; fp-- )
                f = x*f + coef[fp];

            return f;
    }


    /*
     * modrf
     *
     *	uses the modified regula-falsi method to evaluate the root
     * in interval [a,b] of the polynomial described in coef. The
     * root is returned is returned in *val. The routine returns zero
     * if it can't converge.
     */
    boolean modrf(int ord, double coef[], double a, double b, double val[])
    {
      int its;
      double fa, fb, x, fx, lfx;
      
      int fp;

      fb = coef[ord];
      fa = coef[ord];

      for (fp = ord-1; fp >= 0; fp--) {
        fa = a * fa + coef[fp];
        fb = b * fb + coef[fp];
      }

      /*
       * if there is no sign difference the method won't work
       */
      if (fa * fb > 0.0)
        return false;

      /*  commented out to avoid duplicate solutions when the bounds are close to the roots

      if (fabs(fa) < RELERROR)
        {
          *val = a;
          return(1);
        }

      if (fabs(fb) < RELERROR)
        {
          *val = b;
          return(1);
        }
      */

      lfx = fa;

      for (its = 0; its < MAX_ITER_SECANT; its++)
        {
          x = (fb * a - fa * b) / (fb - fa);

          // constrain that x stays in the bounds
          if (x < a || x > b)
            x = 0.5f * (a+b);

          fx = coef[ord];
          for (fp = ord - 1; fp >= 0; fp--)
            fx = x * fx + coef[fp];

          if (Math.abs(x) > RELERROR)
            {
              if (Math.abs(fx / x) < RELERROR)
                {
                  val[0] = x;
                  //	      printf(" x, fx %lf %lf\n", x, fx);
                  return true;
                }
            }
          else if (Math.abs(fx) < RELERROR)
            {
              val[0] = x;
              //	  printf(" x, fx %lf %lf\n", x, fx);
              return true;
            }

          if ((fa * fx) < 0)
            {
              b = x;
              fb = fx;
              if ((lfx * fx) > 0)
                fa /= 2;
            }
          else
            {
              a = x;
              fa = fx;
              if ((lfx * fx) > 0)
                fb /= 2;
            }

          lfx = fx;
        }

      //fprintf(stderr, "modrf overflow %f %f %f\n", a, b, fx);

      return false;
    }




}

