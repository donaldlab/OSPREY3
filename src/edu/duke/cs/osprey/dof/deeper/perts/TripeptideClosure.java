/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper.perts;



/**
 *
 *
 * This class is basically a translation of "tripep_closure.h" from the CSJD paper's C distribution
 * into Java
 * The paper is "A Kinematic View of Loop Closure",  E. A. Coutsias, C. Seok,
M. P. Jacobson, and K. A. Dill, Journal of Computational Chemistry,
Volume 25, Issue 4, Pages 510 - 528 (2004)
 *
 *
 */
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;


public class TripeptideClosure implements Serializable {

//!----------------------------------------------------------------------
//! Copyright (C) 2003
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and Univeristy of New Mexico
//! Witten by Chaok Seok 2003.
//!----------------------------------------------------------------------------
//!----------------------------------------------------------------------------
//!*******************************************************************
//! subroutine  initialize_loop_closure(b_len, b_ang, t_ang)
//!*******************************************************************
//! connectivity of atoms:
//!   N1-A1-C1-N2-A2-C2-N3-A3-C3
//!
//! input:
//!
//  * b_len(1:6): bond lengths (A1-C1, C1-N2, ..., N3-A3)
//!  * b_ang(1:7): bond angles (N1-A1-C1, A1-C1-N2, ..., N3-A3-C3)
//!  * t_ang(1:2): torsion angles (A1-C1-N2-A2, A2-C2-N3-A3)
//!*******************************************************************
//!
//!*******************************************************************
//! subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
//!     r_soln_n, r_soln_a, r_soln_c, n_soln)
//!*******************************************************************
//! input:
//!  * r_n1(3), r_a1(3), r_a3(3), r_c3(3):
//!       Cartesian coordinates of N and CA atoms of the first residue and
//!        CA and C atoms of the last (third) residue.
//! output:
//!  * n_soln: number of alternative loop closure solutions.
//!  * r_soln_n(3,3,8), r_soln_a(3,3,8), r_soln_c(3,3,8):
//!       Cartesian coordinates of loop closure solutions.
//!       first dimension: x, y, z component
//!       second dim: residue number
//!       third dim: solution number
//!*******************************************************************


  double pi = Math.PI;
  double two_pi = 2*pi;
  double deg2rad = pi/180;
  double rad2deg = 180/pi;

  int max_soln = 16;
  int deg_pol = 16;
  int print_level = 0;
//  ! parameters for tripeptide loop (including bond lengths & angles)
  double[] len0 = new double[6], b_ang0 = new double[7], t_ang0 = new double[2];
  double aa13_min_sqr, aa13_max_sqr;
  double[] delta = new double[4], xi = new double[3], eta = new double [3],
          alpha = new double[3], theta = new double[3];

  double[] cos_alpha = new double[3], sin_alpha = new double[3], cos_theta = new double[3],
          sin_theta = new double[3];

  double[] cos_delta = new double[4], sin_delta = new double[4];
  double[] cos_xi = new double[3], cos_eta = new double[3], sin_xi = new double[3], sin_eta = new double[3];
  double[] r_a1a3 = new double[3], r_a1n1 = new double[3], r_a3c3 = new double[3];
  double[] b_a1a3 = new double[3], b_a1n1 = new double[3], b_a3c3 = new double[3];
  double[] len_na = new double[3], len_ac = new double[3], len_aa = new double[3];
//  ! used for polynomial coefficients
  double[][] C0 = new double[3][3], C1 = new double[3][3], C2 = new double[3][3];
  double[][] Q = new double[5][17], R = new double[3][17];

  int n_soln;

  //RotMatrix rm = new RotMatrix();
  SturmSolver sturm;

  //This constructor gets the necessary bond lengths and and angles from the specified residues of the molecule
  //This was not part of the C version
  public TripeptideClosure (List<Residue> tripepRes){

      double b_len[] = new double[6];
      double b_ang[] = new double[7];
      double t_ang[] = new double[2];

      //Coordinates of CAs, Ns, Cs
      double NCoord[][] = new double[3][3];
      double CACoord[][] = new double[3][3];
      double CCoord[][] = new double[3][3];

      for(int a=0;a<3;a++){
          NCoord[a] = tripepRes.get(a).getCoordsByAtomName("N");
          CACoord[a] = tripepRes.get(a).getCoordsByAtomName("CA");
          CCoord[a] = tripepRes.get(a).getCoordsByAtomName("C");
      }

      double vec[][] = new double[8][];

      //Vectors between these mainchain atoms, in order
      vec[0] = VectorAlgebra.subtract(CACoord[0], NCoord[0]);
      vec[1] = VectorAlgebra.subtract(CCoord[0], CACoord[0]);
      vec[2] = VectorAlgebra.subtract(NCoord[1], CCoord[0]);
      vec[3] = VectorAlgebra.subtract(CACoord[1], NCoord[1]);
      vec[4] = VectorAlgebra.subtract(CCoord[1], CACoord[1]);
      vec[5] = VectorAlgebra.subtract(NCoord[2], CCoord[1]);
      vec[6] = VectorAlgebra.subtract(CACoord[2], NCoord[2]);
      vec[7] = VectorAlgebra.subtract(CCoord[2], CACoord[2]);

      for(int a=0; a<6; a++)
          b_len[a] = VectorAlgebra.norm(vec[a+1]);

      for(int a=0; a<7; a++)
          b_ang[a] = Protractor.getAngleRadians( VectorAlgebra.scale(vec[a], -1), vec[a+1]);

      t_ang[0] = calc_dih_ang( vec[1], vec[2], vec[3] );
      t_ang[1] = calc_dih_ang( vec[4], vec[5], vec[6] );

      initialize_loop_closure(b_len, b_ang, t_ang);
  }

  

  public void initialize_loop_closure(double b_len[], double b_ang[], double t_ang[])
 {
  double len1, len2, a_min, a_max;
  double[] rr_a1 = new double[3], rr_c1 = new double[3], rr_n2 = new double[3],
          rr_a2 = new double[3], rr_n2a2_ref = new double[3], rr_c1a1 = new double[3];
  double[] rr_a1a2 = new double[3], dr = new double[3], bb_c1a1 = new double[3], bb_a1a2 = new double[3],
          bb_a2n2 = new double[3];

  double p[] = new double[4];
  double[] mulpro = new double[3];
  double[] tmp_val = new double[3];

  double tol_secant = 1.0e-15;

  int max_iter_sturm = 100;
  int max_iter_secant = 20;

  int i, j;

  sturm = new SturmSolver( tol_secant, max_iter_sturm, max_iter_secant);

  for(i=0;i<6;i++)
    len0[i] = b_len[i];

  for(i=0;i<7;i++)
    b_ang0[i] = b_ang[i];

  for(i=0;i<2;i++)
    t_ang0[i] = t_ang[i];


  for(i=0;i<3;i++)
    rr_c1[i] = 0;

//  do i = 0, 1
  for(i=0;i<2;i++)
   {
     rr_a1[0] = Math.cos(b_ang0[3*i+1])*len0[3*i];
     rr_a1[1] = Math.sin(b_ang0[3*i+1])*len0[3*i];
     rr_a1[2] = 0;
     rr_n2[0] = len0[3*i+1];
     rr_n2[1] = 0;
     rr_n2[2] = 0;

     for(j=0;j<3;j++)
       rr_c1a1[j] = rr_a1[j] - rr_c1[j];

     rr_n2a2_ref[0] = -Math.cos(b_ang0[3*i+2])*len0[3*i+2];
     rr_n2a2_ref[1] = Math.sin(b_ang0[3*i+2])*len0[3*i+2];
     rr_n2a2_ref[2] = 0;
     
     RotationMatrix Us = new RotationMatrix( 1., 0., 0., t_ang0[i], true );
     mulpro = Us.rotateVector(rr_n2a2_ref);
     
     for(j=0;j<3;j++)
      {
       rr_a2[j] =  mulpro[j] + rr_n2[j];
       rr_a1a2[j] = rr_a2[j] - rr_a1[j];
       dr[j] = rr_a1a2[j];
      }

     len2 = VectorAlgebra.dot(dr, dr);
     len1 = Math.sqrt(len2);

     len_aa[i+1] = len1;

     for(j=0;j<3;j++)
      {
       bb_c1a1[j] = rr_c1a1[j]/len0[3*i];
       bb_a1a2[j] = rr_a1a2[j]/len1;
       bb_a2n2[j] = (rr_n2[j] - rr_a2[j])/len0[3*i+2];
      }

//     ! xi

     for(j=0;j<3;j++)
       tmp_val[j] = -bb_a1a2[j];

     xi[i+1] = Protractor.getAngleRadians(tmp_val, bb_a2n2);

//     ! eta

     for(j=0;j<3;j++)
       tmp_val[j] = -bb_c1a1[j];

     eta[i] = Protractor.getAngleRadians(bb_a1a2, tmp_val);
//     ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)

     delta[i+1] = calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2);

     delta[i+1] = pi - delta[i+1];
   }



  a_min = b_ang[3] - (xi[1] + eta[1]);
  a_max = Math.min((b_ang[3] + (xi[1] + eta[1])), pi);

  aa13_min_sqr = Math.pow(len_aa[1],2) + Math.pow(len_aa[2],2) - 2*len_aa[1]*len_aa[2]*Math.cos(a_min);

  aa13_max_sqr = Math.pow(len_aa[1],2) + Math.pow(len_aa[2],2) - 2*len_aa[1]*len_aa[2]*Math.cos(a_max);
 }





public int solve_3pep_poly( double r_n1[], double r_a1[], double r_a3[], double r_c3[], double r_soln_n[][][],
        double r_soln_a[][][], double r_soln_c[][][])
        //Returns the number of solutions
{
  double[] poly_coeff = new double[deg_pol+1], roots = new double[max_soln];

  get_input_angles(r_n1, r_a1, r_a3, r_c3);

  if (n_soln == 0)
    return 0;

  get_poly_coeff(poly_coeff);

  n_soln = sturm.solve_sturm(deg_pol, poly_coeff, roots);

  if (n_soln == 0)//This probably needs to be changed
    return 0;

  coord_from_poly_roots(roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);

  return n_soln;
 }


public void get_input_angles(double r_n1[], double r_a1[], double r_a3[], double r_c3[])
 {
  double dr_sqr;
  double[] tmp_val = new double[3];
  int i;
  char[] cone_type = new char[2];

  n_soln = max_soln;

  for(i=0;i<3;i++)
    r_a1a3[i] = r_a3[i] - r_a1[i];

  dr_sqr = VectorAlgebra.dot(r_a1a3,r_a1a3);
  len_aa[0] = Math.sqrt(dr_sqr);

  if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr))
   {
    n_soln = 0;
    return;
   }

//  ! bond lengths
  for(i=0;i<3;i++)
    r_a1n1[i] = r_n1[i] - r_a1[i];

  len_na[0] = Math.sqrt(VectorAlgebra.dot(r_a1n1,r_a1n1));
  len_na[1] = len0[2];
  len_na[2] = len0[5];

  for(i=0;i<3;i++)
    r_a3c3[i] = r_c3[i] - r_a3[i];

  len_ac[0] = len0[0];
  len_ac[1] = len0[3];
  len_ac[2] = Math.sqrt(VectorAlgebra.dot(r_a3c3,r_a3c3));

//  ! unit vectors
  for(i=0;i<3;i++)
   {
    b_a1n1[i] = r_a1n1[i]/len_na[0];
    b_a3c3[i] = r_a3c3[i]/len_ac[2];
    b_a1a3[i] = r_a1a3[i]/len_aa[0];
   }

//  ! delta(3): dih of N(1)CA(1)CA(3)C(3)
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1n1[i];

  delta[3] = calc_dih_ang(tmp_val, b_a1a3, b_a3c3);

  delta[0] = delta[3];

  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1a3[i];
  xi[0] = Protractor.getAngleRadians(tmp_val, b_a1n1);

  eta[2] = Protractor.getAngleRadians(b_a1a3, b_a3c3);

  for(i=0;i<3;i++)
   {
     cos_delta[i+1] = Math.cos(delta[i+1]);
     sin_delta[i+1] = Math.sin(delta[i+1]);

     cos_xi[i] = Math.cos(xi[i]);
     sin_xi[i] = Math.sin(xi[i]);
     cos_eta[i] = Math.cos(eta[i]);
     sin_eta[i] = Math.sin(eta[i]);
   }

  cos_delta[0] = cos_delta[3];
  sin_delta[0] = sin_delta[3];

//  ! theta (N, CA, C) bond angle
  theta[0] = b_ang0[0];
  theta[1] = b_ang0[3];
  theta[2] = b_ang0[6];

  for(i=0;i<3;i++)
    cos_theta[i] = Math.cos(theta[i]);

//  ! alpha
  cos_alpha[0] = -(Math.pow(len_aa[0],2) + Math.pow(len_aa[1],2) - Math.pow(len_aa[2],2))/(2*len_aa[0]*len_aa[1]);
  alpha[0] = Math.acos(cos_alpha[0]);
  sin_alpha[0] = Math.sin(alpha[0]);
  cos_alpha[1] = (Math.pow(len_aa[1],2) + Math.pow(len_aa[2],2) - Math.pow(len_aa[0],2))/(2*len_aa[1]*len_aa[2]);
  alpha[1] = Math.acos(cos_alpha[1]);
  sin_alpha[1] = Math.sin(alpha[1]);
  alpha[2] = pi - alpha[0] + alpha[1];
  cos_alpha[2] = Math.cos(alpha[2]);
  sin_alpha[2] = Math.sin(alpha[2]);

  if (print_level > 0)
   {
//     write(*,'(a,3f9.4)') 'xi = ', xi(1:3)*rad2deg
     System.out.printf("xi = %9.4f%9.4f%9.4f\n", xi[0]*rad2deg, xi[1]*rad2deg, xi[2]*rad2deg);
//     write(*,'(a,3f9.4)') 'eta = ', eta(1:3)*rad2deg
     System.out.printf("eta = %9.4f%9.4f%9.4f\n", eta[0]*rad2deg, eta[1]*rad2deg, eta[2]*rad2deg);
//     write(*,'(a,3f9.4)') 'delta = ', delta(1:3)*rad2deg
     System.out.printf("delta = %9.4f%9.4f%9.4f\n", delta[1]*rad2deg, delta[2]*rad2deg, delta[3]*rad2deg);
//     write(*,'(a,3f9.4)') 'theta = ', theta(1:3)*rad2deg
     System.out.printf("theta = %9.4f%9.4f%9.4f\n", theta[0]*rad2deg, theta[1]*rad2deg, theta[2]*rad2deg);
//     write(*,'(a,3f9.4)') 'alpha = ', alpha(1:3)*rad2deg
     System.out.printf("alpha = %9.4f%9.4f%9.4f\n", alpha[0]*rad2deg, alpha[1]*rad2deg, alpha[2]*rad2deg);
//  end if
   }

//  ! check for existence of soln
  for(i=0;i<3;i++)
   {
    test_two_cone_existence_soln(theta[i], xi[i], eta[i], alpha[i], cone_type);
    if (n_soln == 0)
      return;
   }

  return;
//end subroutine get_input_angles
 }



public void test_two_cone_existence_soln(double tt, double kx, double et, double ap, char cone_type[])
 {

  double at, ex, abs_at, ap1, kx1, et1;
  double cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2;
  int s1, s2, t1, t2;
  boolean complicated = false;

  n_soln = max_soln;


  ap1 = ap;
  kx1 = kx;
  et1 = et;

  at = ap1 - tt;
  ex = kx1 + et1;
  abs_at = Math.abs(at);

//  ! case of no soln
  if (abs_at > ex)
   {
    n_soln = 0;
    return;
   }

  
//     ! find type of intersection
  if (complicated)
   {
    cos_tx1 = Math.cos(tt+kx1);
    cos_tx2 = Math.cos(tt-kx1);
    cos_te1 = Math.cos(tt+et1);
    cos_te2 = Math.cos(tt-et1);
    cos_ea1 = Math.cos(et1+ap1);
    cos_ea2 = Math.cos(et1-ap1);
    cos_xa1 = Math.cos(kx1+ap1);
    cos_xa2 = Math.cos(kx1-ap1);

    s1 = 0;
    s2 = 0;
    t1 = 0;
    t2 = 0;

    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0e0)
      s1 = 0;

    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0e0)
      s2 = 0;

    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0e0)
      t1 = 0;

    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0e0)
      t2 = 0;
   }

  return;
 }



public void get_poly_coeff(double poly_coeff[])
 {
  int i, j;
  double A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42;

  double[] B0 = new double[3], B1 = new double[3], B2 = new double[3], B3 = new double[3],
          B4 = new double[3], B5 = new double[3], B6 = new double[3], B7 = new double[3], B8 = new double[3];

  double[][] u11 = new double[5][5], u12 = new double[5][5], u13 = new double[5][5], u31 = new double[5][5],
          u32 = new double[5][5], u33 = new double[5][5];

  double[][] um1 = new double[5][5], um2 = new double[5][5], um3 = new double[5][5], um4 = new double[5][5],
          um5 = new double[5][5], um6 = new double[5][5], q_tmp = new double[5][5];

  int[] p1 = new int[2], p3 = new int[2], p_um1 = new int[2], p_um2 = new int[2], p_um3 = new int[2],
          p_um4 = new int[2], p_um5 = new int[2], p_um6 = new int[2], p_Q = new int[2];

  int p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, p_f9;
  int p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, p_f18;
  int p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26;

  int p_final;

  double[] f1 = new double[17], f2 = new double[17], f3 = new double[17], f4 = new double[17], f5 = new double[17],
          f6 = new double[17], f7 = new double[17], f8 = new double[17], f9 = new double[17];
  
  double[] f10 = new double[17], f11 = new double[17], f12 = new double[17], f13 = new double[17], f14 = new double[17],
          f15 = new double[17], f16 = new double[17], f17 = new double[17], f18 = new double[17];

  double[] f19 = new double[17], f20 = new double[17], f21 = new double[17], f22 = new double[17], f23 = new double[17],
          f24 = new double[17], f25 = new double[17], f26 = new double[17];

  for(i=0;i<3;i++)
   {
    A0 = cos_alpha[i]*cos_xi[i]*cos_eta[i] - cos_theta[i];
    A1 = -sin_alpha[i]*cos_xi[i]*sin_eta[i];
    A2 = sin_alpha[i]*sin_xi[i]*cos_eta[i];
    A3 = sin_xi[i]*sin_eta[i];
    A4 = A3*cos_alpha[i];
    
    j = i;
    A21 = A2*cos_delta[j];
    A22 = A2*sin_delta[j];
    A31 = A3*cos_delta[j];
    A32 = A3*sin_delta[j];
    A41 = A4*cos_delta[j];
    A42 = A4*sin_delta[j];
    B0[i] = A0 + A22 + A31;
    B1[i] = 2*(A1 + A42);
    B2[i] = 2*(A32 - A21);
    B3[i] = -4*A41;
    B4[i] = A0 + A22 - A31;
    B5[i] = A0 - A22 - A31;
    B6[i] = -2*(A21 + A32);
    B7[i] = 2*(A1 - A42);
    B8[i] = A0 - A22 + A31;
   }


  i = 0;
  C0[i][0] = B0[i];
  C0[i][1] = B2[i];
  C0[i][2] = B5[i];
  C1[i][0] = B1[i];
  C1[i][1] = B3[i];
  C1[i][2] = B7[i];
  C2[i][0] = B4[i];
  C2[i][1] = B6[i];
  C2[i][2] = B8[i];

  for(i=1;i<3;i++)
   {
    C0[i][0] = B0[i];
    C0[i][1] = B1[i];
    C0[i][2] = B4[i];
    C1[i][0] = B2[i];
    C1[i][1] = B3[i];
    C1[i][2] = B6[i];
    C2[i][0] = B5[i];
    C2[i][1] = B7[i];
    C2[i][2] = B8[i];
   }

//  ! first determinant
  for(i=0;i<3;i++)
   {
     u11[0][i] = C0[0][i];
     u12[0][i] = C1[0][i];
     u13[0][i] = C2[0][i];
     u31[i][0] = C0[1][i];
     u32[i][0] = C1[1][i];
     u33[i][0] = C2[1][i];
   }

  p1[0] = 2;
  p1[1] = 0;
  p3[0] = 0;
  p3[1] = 2;

  poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
  poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
  poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
  poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
  poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
  poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
  poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      Q[i][j] = q_tmp[i][j];

//  ! second determinant
  for(i=0;i<3;i++)
    for(j=0;j<17;j++)
      R[i][j] = 0;

  for(i=0;i<3;i++)
   {
    R[0][i] = C0[2][i];
    R[1][i] = C1[2][i];
    R[2][i] = C2[2][i];
   }

  p2 = 2;
  p4 = 4;


  p_f1 = poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1);
  p_f2 = poly_mul1(R[1], R[2], p2, p2, f2);
  p_f3 = poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1, p2, p_f2, f3);
  p_f4 = poly_mul1(R[2], f1, p2, p_f1, f4);
  p_f5 = poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3, p2, p_f4, f5);

  p_f6 = poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6);
  p_f7 = poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1, p2, p_f6, f7);
  p_f8 = poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3, p2, p_f7, f8);
  p_f9 = poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5, p2, p_f8, f9);

  p_f10 = poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10);
  p_f11 = poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1, p2, p_f10, f11);
  p_f12 = poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3, p2, p_f11, f12);

  p_f13 = poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13);
  p_f14 = poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1, p2, p_f13, f14);
  p_f15 = poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15);
  p_f16 = poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1, p2, p_f15, f16);
  p_f17 = poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14, p4, p_f16, f17);

  p_f18 = poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18);
  p_f19 = poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19);
  p_f20 = poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19, p4, p_f18, f20);
  p_f21 = poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21);
  p_f22 = poly_mul1(Q[4], f21, p4, p_f21, f22);
  p_f23 = poly_sub1(f20, f22, p_f20, p_f22, f23);
  p_f24 = poly_mul1(R[0], f23, p2, p_f23, f24);
  p_f25 = poly_sub1(f17, f24, p_f17, p_f24, f25);
  p_f26 = poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12, p2, p_f25, f26);
  p_final = poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9, p2, p_f26, poly_coeff);


  if (p_final != 16)
   {
    throw new Error("Degree of polynomial is not 16!");
   }


  if (poly_coeff[16] < 0)
    for(i=0;i<17;i++)
     poly_coeff[i] *= -1.0;

  if (print_level > 0)
   {
     System.out.printf("poly_coeff\n");
     for(i=0;i<17;i++)
        System.out.printf("%5d%15.6f\n", i, poly_coeff[i]);
   }

  return;
 }



void poly_mul_sub2(double u1[][], double u2[][], double u3[][], double u4[][], int p1[], int p2[], int p3[],
        int p4[], double u5[][], int p5[])
 {

  double[][] d1 = new double[5][5], d2 = new double[5][5];
  int[] pd1 = new int[2], pd2 = new int[2];

  poly_mul2(u1, u2, p1, p2, d1, pd1);
  poly_mul2(u3, u4, p3, p4, d2, pd2);
  poly_sub2(d1, d2, pd1, pd2, u5, p5);


 }



void poly_mul2(double u1[][], double u2[][], int p1[], int p2[], double u3[][], int p3[])
 {
  int i1, j1, i2, j2, i3, j3, p11, p12, p21, p22;
  int i, j;

  double u1ij;

  for(i=0;i<2;i++)
    p3[i] = p1[i] + p2[i];
  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      u3[i][j] = 0;

  p11 = p1[0];
  p12 = p1[1];
  p21 = p2[0];
  p22 = p2[1];

  for(i1=0;i1<=p12;i1++)
   {
     for(j1=0;j1<=p11;j1++)
      {
        u1ij = u1[i1][j1];
	for(i2=0;i2<=p22;i2++)
	 {
           i3 = i1 + i2;
	   for (j2=0;j2<=p21;j2++)
	    {
              j3 = j1 + j2;
              u3[i3][j3] = u3[i3][j3] + u1ij*u2[i2][j2];
	    }
	 }
      }
   }

 }

void poly_sub2(double u1[][], double u2[][], int p1[], int p2[], double u3[][], int p3[])
 {
  int i, j, p11, p12, p21, p22;
  boolean i1_ok, i2_ok;

  p11 = p1[0];
  p12 = p1[1];

  p21 = p2[0];
  p22 = p2[1];

  p3[0] = Math.max(p11,p21);
  p3[1] = Math.max(p12,p22);

  for(i=0;i<=p3[1];i++)
   {
    i1_ok = (i > p12);
    i2_ok = (i > p22);
    for(j=0;j<=p3[0];j++)
     {
      if (i2_ok || (j > p21))
        u3[i][j] = u1[i][j];
      else if (i1_ok || (j > p11))
        u3[i][j] = -u2[i][j];
      else
        u3[i][j] = u1[i][j] - u2[i][j];
     }
   }

  return;
 }



public int poly_mul_sub1(double u1[], double u2[], double u3[], double u4[], int p1, int p2, int p3, int p4, double u5[])
 {
  double[] d1 = new double[17], d2 = new double[17];
  int pd1, pd2;


  pd1 = poly_mul1(u1, u2, p1, p2, d1);
  pd2 = poly_mul1(u3, u4, p3, p4, d2);
  int p5 = poly_sub1(d1, d2, pd1, pd2, u5);

  return p5;
 }



public int poly_mul1(double u1[], double u2[], int p1, int p2, double u3[])
 {
  int i, i1, i2, i3;
  double u1i;

  int p3 = p1 + p2;

  for(i=0;i<17;i++)
    u3[i] = 0;

  for(i1=0;i1<=p1;i1++)
   {
    u1i = u1[i1];

    for(i2=0;i2<=p2;i2++)
     {
        i3 = i1 + i2;
        u3[i3] = u3[i3] + u1i*u2[i2];
     }
   }

  return p3;
 }



public int poly_sub1(double u1[], double u2[], int p1, int p2, double u3[])
 {
  int i;


  int p3 = Math.max(p1, p2);

  for(i=0;i<=p3;i++)
   {
    if (i > p2)
      u3[i] = u1[i];

    else if (i > p1)
      u3[i] = -u2[i];

    else
      u3[i] = u1[i] - u2[i];
   }

  return p3;
 }



void coord_from_poly_roots(double roots[], double r_n1[], double r_a1[], double r_a3[], double r_c3[], double r_soln_n[][][], double r_soln_a[][][], double r_soln_c[][][])
 {
  double[] ex = new double[3], ey = new double[3], ez = new double[3], b_a1a2 = new double[3],
          b_a3a2 = new double[3], r_tmp = new double[3];

  double[][] p_s = new double[3][3], s1 = new double[3][3], s2 = new double[3][3], p_t = new double[3][3],
          t1 = new double[3][3], t2 = new double[3][3];

  double[][] p_s_c = new double[3][3], s1_s = new double[3][3], s2_s = new double[3][3], p_t_c = new double[3][3],
          t1_s = new double[3][3], t2_s = new double[3][3];

  double angle, sig1_init, half_tan[] = new double[3];

  double[] cos_tau = new double[4], sin_tau = new double[4], cos_sig = new double[3], sin_sig = new double[3];
  double ht, tmp, sig1;

  double[] r_s = new double[3], r_t = new double[3], r0 = new double[3];
  double[][] r_n = new double[3][3], r_a = new double[3][3], r_c = new double[3][3];
  double[] p = new double[4];

  int i_soln, i, j;
  double a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3;
  double[] rr_a1c1 = new double[3], rr_c1n2 = new double[3], rr_n2a2 = new double[3], rr_a2c2 = new double[3],
          rr_c2n3 = new double[3], rr_n3a3 = new double[3], rr_a1a2 = new double[3], rr_a2a3 = new double[3];

  double a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3;
  double tmp_value, ex_tmp[] = new double[3];
  double[] tmp_array = new double[3], tmp_array1 = new double[3], tmp_array2 = new double[3], tmp_array3 = new double[3];
  double[] mat1 = new double[3], mat2 = new double[3], mat3 = new double[3], mat4 = new double[3], mat5 = new double[3];
  double[] mat11 = new double[3], mat22 = new double[3], mat33 = new double[3], mat44 = new double[3], mat55 = new double[3];

  if ( n_soln == 0 )
    return;

//  ! Define body frame (ex, ey, ez)
  for(i=0;i<3;i++)
    ex[i] = b_a1a3[i];
//  call cross(r_a1n1, ex, ez)
  ez = VectorAlgebra.cross(r_a1n1, ex);

  tmp_value = Math.sqrt(VectorAlgebra.dot(ez,ez));
  for(i=0;i<3;i++)
    ez[i] = ez[i]/tmp_value;

  ey = VectorAlgebra.cross(ez, ex);
//  ! vertual bond vectors in the reference plane
  for(i=0;i<3;i++)
   {
    b_a1a2[i] = -cos_alpha[0]*ex[i] + sin_alpha[0]*ey[i];
    b_a3a2[i] = cos_alpha[2]*ex[i] + sin_alpha[2]*ey[i];
   }


//  !! Define cone coordinates for each angle joint.
//  ! (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems

//  ! residue 1
  for(i=0;i<3;i++)
   {
    p_s[0][i] = -ex[i];
    s1[0][i]  = ez[i];
    s2[0][i]  = ey[i];
    p_t[0][i] = b_a1a2[i];
    t1[0][i]  = ez[i];
    t2[0][i]  = sin_alpha[0]*ex[i] + cos_alpha[0]*ey[i];
   }

//  ! residue 2
  for(i=0;i<3;i++)
   {
    p_s[1][i] = -b_a1a2[i];
    s1[1][i]  = -ez[i];
    s2[1][i]  = t2[0][i];
    p_t[1][i] = -b_a3a2[i];
    t1[1][i]  = -ez[i];
    t2[1][i]  = sin_alpha[2]*ex[i] - cos_alpha[2]*ey[i];
   }

//  ! residue 3
  for(i=0;i<3;i++)
   {
    p_s[2][i] = b_a3a2[i];
    s2[2][i]  = t2[1][i];
    s1[2][i]  = ez[i];
    p_t[2][i] = ex[i];
    t1[2][i] =  ez[i];
    t2[2][i] = -ey[i];
   }

//  ! scale vectors
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
     {
      p_s_c[i][j] = p_s[i][j]*cos_xi[i];
      s1_s[i][j]  = s1[i][j]*sin_xi[i];
      s2_s[i][j]  = s2[i][j]*sin_xi[i];
      p_t_c[i][j] = p_t[i][j]*cos_eta[i];
      t1_s[i][j]  = t1[i][j]*sin_eta[i];
      t2_s[i][j]  = t2[i][j]*sin_eta[i];
     }

//  ! initial sig(1)
  for(i=0;i<3;i++)
    r_tmp[i] = (r_a1n1[i]/len_na[0] - p_s_c[0][i])/sin_xi[0];
  
  angle = Protractor.getAngleRadians(s1[0], r_tmp);
  sig1_init = Math.copySign(angle, VectorAlgebra.dot(r_tmp,s2[0]));

//  ! CA
  for(i=0;i<3;i++)
   {
    r_a[0][i] = r_a1[i];
    r_a[1][i] = r_a1[i] + len_aa[1]*b_a1a2[i];
    r_a[2][i] = r_a3[i];
    r0[i] = r_a1[i];
   }


  for(i_soln=0;i_soln<n_soln;i_soln++)
   {
     half_tan[2] = roots[i_soln];
     half_tan[1] = calc_t2(half_tan[2]);
     half_tan[0] = calc_t1(half_tan[2], half_tan[1]);

     for(i=1;i<=3;i++)
      {
       ht = half_tan[i-1];
       tmp = 1 + ht*ht;
       cos_tau[i] = (1 - ht*ht)/tmp;
       sin_tau[i] = 2*ht/tmp;
      }

     cos_tau[0] = cos_tau[3];
     sin_tau[0] = sin_tau[3];

     for(i=0;i<3;i++)
      {
       cos_sig[i] = cos_delta[i]*cos_tau[i] + sin_delta[i]*sin_tau[i];
       sin_sig[i] = sin_delta[i]*cos_tau[i] - cos_delta[i]*sin_tau[i];
      }

     for(i=0;i<3;i++)
      for(j=0;j<3;j++)
       {
        r_s[j] = p_s_c[i][j] + cos_sig[i]*s1_s[i][j] + sin_sig[i]*s2_s[i][j];
        r_t[j] = p_t_c[i][j] + cos_tau[i+1]*t1_s[i][j] + sin_tau[i+1]*t2_s[i][j];
        r_n[i][j] = r_s[j]*len_na[i] + r_a[i][j];
        r_c[i][j] = r_t[j]*len_ac[i] + r_a[i][j];
       }

//     ! rotate back atoms by -(sig(1) - sig1_init) around -ex
     sig1 = Math.atan2(sin_sig[0], cos_sig[0]);

     ex_tmp[0] = -ex[0];
     ex_tmp[1] = -ex[1];
     ex_tmp[2] = -ex[2];
     tmp_value = -(sig1-sig1_init);
     
     RotationMatrix Us = new RotationMatrix(-ex[0], -ex[1], -ex[2], tmp_value, true);

     for(i=0;i<3;i++)
      {
       mat11[i] = r_c[0][i]-r0[i];
       mat22[i] = r_n[1][i]-r0[i];
       mat33[i] = r_a[1][i]-r0[i];
       mat44[i] = r_c[1][i]-r0[i];
       mat55[i] = r_n[2][i]-r0[i];
      }
     
     
     mat1 = Us.rotateVector(mat11);
     mat2 = Us.rotateVector(mat22);
     mat3 = Us.rotateVector(mat33);
     mat4 = Us.rotateVector(mat44);
     mat5 = Us.rotateVector(mat55);
     
     
     for(i=0;i<3;i++)
      {
       r_soln_n[i_soln][0][i] = r_n1[i];
       r_soln_a[i_soln][0][i] = r_a1[i];
       r_soln_c[i_soln][0][i] = mat1[i] + r0[i];
       r_soln_n[i_soln][1][i] = mat2[i] + r0[i];
       r_soln_a[i_soln][1][i] = mat3[i] + r0[i];
       r_soln_c[i_soln][1][i] = mat4[i] + r0[i];
       r_soln_n[i_soln][2][i] = mat5[i] + r0[i];
       r_soln_a[i_soln][2][i] = r_a3[i];
       r_soln_c[i_soln][2][i] = r_c3[i];
      }


     if (print_level > 0)
      {
        System.out.printf("roots: t0, t2, t1 %d\n", i_soln);
        System.out.printf("%15.6f %15.6f %15.6f\n", half_tan[2], half_tan[1], half_tan[0]);

	for(i=0;i<3;i++)
	 {
          rr_a1c1[i] = r_soln_c[i_soln][0][i] - r_soln_a[i_soln][0][i];
          rr_c1n2[i] = r_soln_n[i_soln][1][i] - r_soln_c[i_soln][0][i];
          rr_n2a2[i] = r_soln_a[i_soln][1][i] - r_soln_n[i_soln][1][i];
          rr_a2c2[i] = r_soln_c[i_soln][1][i] - r_soln_a[i_soln][1][i];
          rr_c2n3[i] = r_soln_n[i_soln][2][i] - r_soln_c[i_soln][1][i];
          rr_n3a3[i] = r_soln_a[i_soln][2][i] - r_soln_n[i_soln][2][i];
          rr_a1a2[i] = r_soln_a[i_soln][1][i] - r_soln_a[i_soln][0][i];
          rr_a2a3[i] = r_soln_a[i_soln][2][i] - r_soln_a[i_soln][1][i];
	 }


        a1c1 = Math.sqrt(VectorAlgebra.dot(rr_a1c1, rr_a1c1));
        c1n2 = Math.sqrt(VectorAlgebra.dot(rr_c1n2, rr_c1n2));
        n2a2 = Math.sqrt(VectorAlgebra.dot(rr_n2a2, rr_n2a2));
        a2c2 = Math.sqrt(VectorAlgebra.dot(rr_a2c2, rr_a2c2));
        c2n3 = Math.sqrt(VectorAlgebra.dot(rr_c2n3, rr_c2n3));
        n3a3 = Math.sqrt(VectorAlgebra.dot(rr_n3a3, rr_n3a3));
        a1a2 = Math.sqrt(VectorAlgebra.dot(rr_a1a2, rr_a1a2));
        a2a3 = Math.sqrt(VectorAlgebra.dot(rr_a2a3, rr_a2a3));

        System.out.printf("na: n2a2, n3a3 = %9.3f%9.3f%9.3f%9.3f\n", len0[2], n2a2, len0[5], n3a3);
        System.out.printf("ac: a1c1, a2c2 = %9.3f%9.3f%9.3f%9.3f\n", len0[0], a1c1, len0[3], a2c2);
        System.out.printf("cn: c1n2, c2n3 = %9.3f%9.3f%9.3f%9.3f\n", len0[1], c1n2, len0[4], c2n3);
        System.out.printf("aa: a1a2, a2a3 = %9.3f%9.3f%9.3f%9.3f\n", len_aa[1], a1a2, len_aa[2], a2a3);

	for(i=0;i<3;i++)
	  tmp_array[i] = rr_a1a2[i]/a1a2;
        a3a1a2 = Protractor.getAngleRadians(b_a1a3, tmp_array);

        for(i=0;i<3;i++)
	  tmp_array[i] = rr_a2a3[i]/a2a3;

        a2a3a1 = Protractor.getAngleRadians(tmp_array, b_a1a3);

        System.out.printf("alpha1, alpha3 = %9.3f%9.3f\n", (pi-a3a1a2)*rad2deg, (pi-a2a3a1)*rad2deg);

        for(i=0;i<3;i++)
	  tmp_array[i] = rr_a1c1[i]/a1c1;
        n1a1c1 = Protractor.getAngleRadians(b_a1n1, tmp_array);

        for(i=0;i<3;i++)
	  tmp_array[i] = -rr_n3a3[i]/n3a3;
        n3a3c3 = Protractor.getAngleRadians(b_a3c3, tmp_array);

        for(i=0;i<3;i++)
	 {
	  tmp_array1[i] = rr_a2c2[i]/a2c2;
	  tmp_array2[i] = -rr_n2a2[i]/n2a2;
	 }
        n2a2c2 = Protractor.getAngleRadians(tmp_array1, tmp_array2);

        System.out.printf("ang_nac = %9.3f%9.3f\n", b_ang0[0]*rad2deg, n1a1c1*rad2deg);
        System.out.printf("ang_nac = %9.3f%9.3f\n", b_ang0[3]*rad2deg, n2a2c2*rad2deg);
        System.out.printf("ang_nac = %9.3f%9.3f\n", b_ang0[6]*rad2deg, n3a3c3*rad2deg);

	for(i=0;i<3;i++)
	 {
	  tmp_array1[i] = rr_a1c1[i]/a1c1;
	  tmp_array2[i] = rr_c1n2[i]/c1n2;
	  tmp_array3[i] = rr_n2a2[i]/n2a2;
	 }
        a1c1n2a2 = calc_dih_ang(tmp_array1, tmp_array2, tmp_array3);

        for(i=0;i<3;i++)
	 {
	  tmp_array1[i] = rr_a2c2[i]/a2c2;
	  tmp_array2[i] = rr_c2n3[i]/c2n3;
	  tmp_array3[i] = rr_n3a3[i]/n3a3;
	 }
        a2c2n3a3 = calc_dih_ang(tmp_array1, tmp_array2, tmp_array3);

        System.out.printf("t_ang1 = %9.3f%9.3f\n", t_ang0[0]*rad2deg, a1c1n2a2*rad2deg);
        System.out.printf("t_ang2 = %9.3f%9.3f\n", t_ang0[1]*rad2deg, a2c2n3a3*rad2deg);
      }

//  end do
   }

  return;
 }



double calc_t2(double t0)
 {

  double tmp_value;
  double B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3;
  double K0, K1, K2, K3, t0_2, t0_3, t0_4;

  t0_2 = t0*t0;
  t0_3 = t0_2*t0;
  t0_4 = t0_3*t0;


  A0 = Q[0][0] + Q[0][1]*t0 + Q[0][2]*t0_2 + Q[0][3]*t0_3 + Q[0][4]*t0_4;
  A1 = Q[1][0] + Q[1][1]*t0 + Q[1][2]*t0_2 + Q[1][3]*t0_3 + Q[1][4]*t0_4;
  A2 = Q[2][0] + Q[2][1]*t0 + Q[2][2]*t0_2 + Q[2][3]*t0_3 + Q[2][4]*t0_4;
  A3 = Q[3][0] + Q[3][1]*t0 + Q[3][2]*t0_2 + Q[3][3]*t0_3 + Q[3][4]*t0_4;
  A4 = Q[4][0] + Q[4][1]*t0 + Q[4][2]*t0_2 + Q[4][3]*t0_3 + Q[4][4]*t0_4;


  B0 = R[0][0] + R[0][1]*t0 + R[0][2]*t0_2;
  B1 = R[1][0] + R[1][1]*t0 + R[1][2]*t0_2;
  B2 = R[2][0] + R[2][1]*t0 + R[2][2]*t0_2;


  B2_2 = B2*B2;
  B2_3 = B2_2*B2;


  K0 = A2*B2 - A4*B0;
  K1 = A3*B2 - A4*B1;
  K2 = A1*B2_2 - K1*B0;
  K3 = K0*B2 - K1*B1;

  tmp_value = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1);

  return tmp_value;
 }



double calc_t1(double t0, double t2)
 {
  double tmp_value;
  double U11, U12, U13, U31, U32, U33;
  double t0_2, t2_2;


  t0_2 = t0*t0;
  t2_2 = t2*t2;

  U11 = C0[0][0] + C0[0][1]*t0 + C0[0][2]*t0_2;
  U12 = C1[0][0] + C1[0][1]*t0 + C1[0][2]*t0_2;
  U13 = C2[0][0] + C2[0][1]*t0 + C2[0][2]*t0_2;

  U31 = C0[1][0] + C0[1][1]*t2 + C0[1][2]*t2_2;
  U32 = C1[1][0] + C1[1][1]*t2 + C1[1][2]*t2_2;
  U33 = C2[1][0] + C2[1][1]*t2 + C2[1][2]*t2_2;

  tmp_value = (U31*U13-U11*U33)/(U12*U33-U13*U32);

  return tmp_value;
 }


public double calc_dih_ang(double r1[], double r2[], double r3[])
 {
//!-----------------------------------------------------------------------
//! r1=Rab, r2=Rbc, r3=Rcd : angle between planes abc and bcd
//!-----------------------------------------------------------------------
  double[] p = new double[3], q = new double[3], s = new double[3];
  double arg;

  p = VectorAlgebra.cross(r1, r2);
  q = VectorAlgebra.cross(r2, r3);
  s = VectorAlgebra.cross(r3, r1);

  arg = VectorAlgebra.dot(p,q)/Math.sqrt(VectorAlgebra.dot(p,p)*VectorAlgebra.dot(q,q));

  arg = Math.copySign(Math.min(Math.abs(arg),1),arg);
  return Math.copySign(Math.acos(arg), VectorAlgebra.dot(s,r2));
 }


}



