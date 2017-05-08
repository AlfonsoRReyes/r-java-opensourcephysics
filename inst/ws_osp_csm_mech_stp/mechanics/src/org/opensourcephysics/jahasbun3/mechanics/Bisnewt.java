/*
@Author J E Hasbun 2007.
This finds the root of a function in a given interval. This function
is implemented through the OSP numerics.MultiVarFunction interface and
refers to the function fx() which is itself implemented through the OSP
numerics.Function interface.
@Copyright (c) 2007
 This software is to support Intermediate Classical Mechanics
 with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
 http://www.opensourcephysics.org under the terms of the GNU General Public
 License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.numerics.Function;

   public class Bisnewt {

   public double Bisnewt (double xleft, double xright, int icmax,
                           double tol, Function fx) {
     /*
     Uses Newton Raphson or bisection methods to find the root
     The root is between xleft and xright and xbest is the result.
     A better guess is accepted if it is within the known bounds,
     else the bisection step is taken.
     Ref: Computational Physics by P. L. Devries (J. Wiley, 1993)
     input: _x is an array of size 4 such that
     [xleft,xright] is the interval wherein fx has a root, icmax is the
     maximum iteration number, and tol is the tolerance level
     output: returns xbest as the value of the function
     Reasonable values of icmax and tol are 25, 5e-3
     */
     double xbest,fleft,fright,fbest,derfbest,delta,rtest=10*tol;
     int icount=0,iflag=0;//loop counter
     //variables
     fleft=fx.evaluate(xleft);
     fright=fx.evaluate(xright);
     if(fleft*fright >= 0) {iflag=1;}
     switch (iflag){
       case 1:
        System.out.println("No solution possible");
        break;
     }
     if(Math.abs(fleft) <= Math.abs(fright) ){
       xbest = xleft;
       fbest = fleft;
     }
     else {
       xbest = xright;
       fbest = fright;
     }
     derfbest = fxprime(xbest, fx);
     while (icount < icmax && rtest > tol){
       icount++;
       //decide Newton-Raphson or Bisection method to do:
       if((derfbest*(xbest-xleft)-fbest)*
         (derfbest*(xbest-xright)-fbest) <= 0){
         //Newton-Raphson step
         delta=-fbest/derfbest;
         xbest=xbest + delta;
         //System.out.println("Newton: count="+icount+", fx="+fbest);
       }
       else{
         //bisection step
         delta = (xright - xleft) / 2;
         xbest = (xleft + xright) / 2;
         //System.out.println("Bisection: count="+icount+", fx ="+fbest);
       }
       rtest=Math.abs(delta/xbest);
       //Compare the relative error to the tolerance
       if(rtest <= tol) {
         //if the error is small, the root has been found
         System.out.println("root found="+xbest);
       }
       else{
         //the error is still large, so loop
         fbest=fx.evaluate(xbest);
         derfbest=fxprime(xbest, fx);
         //adjust brackets
         if(fleft*fbest <= 0) {
           //root is in the xleft subinterval:
           xright = xbest;
           fright = fbest;
         }
         else{
           //root is in the xright subinterval:
           xleft = xbest;
           fleft = fbest;
         }
       }
     }
   //reach here if either the error is too large or icount reached icmax
   if (icount > icmax || rtest > tol){
     System.out.println(icmax + " trials made - no convergence achieved");
   }
   return xbest;
   }

   public double fxprime(double x, Function fx){
   //numerical derivative evaluation
     double del=(1.0/3.0)*1.e-3;
     double der=(fx.evaluate(x+del)-fx.evaluate(x-del))/del/2.0;
     //System.out.println("fxprime="+der);
     return der;
   }
 }