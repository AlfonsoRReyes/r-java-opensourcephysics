/**
 @Author J E Hasbun 2007
 Uses Simpson's rule to find the area of an array representing a
 function that's been evaluated at N intervals of size h, where N is
 an odd integer.

 Example of usage:
     int N=27;
     x=new double[N]; f=new double[N];
     double a=0, b=5, h=(b-a)/(N-1);
     for (int i=0; i< N;i++){
       x[i]=a+i*h;
       f[i]=x[i]*Math.exp(x[i]);
     }
      double sum=Simp.Simp(f,h);

 Results: sum=594.6615858178942

 @Copyright (c) 2007
 This software is to support the Open Source Physics library
 http://www.opensourcephysics.org under the terms of the GNU General Public
 License (GPL) as published by the Free Software Foundation.
**/

package org.opensourcephysics.jahasbun3.mechanics;
public class Simp{

   public static double Simp(double f[], double h) {
     //Simpson's rule for numerical integration
     //f is an odd array of evaluated functions in steps h
     int ip = f.length; //must be an odd number
     double sumOdd = 0.0, sumEven = 0.0;
     for(int i = 1; i < ip-1;i+=2) {
       sumOdd +=f[i];
       sumEven +=f[i+1];
     }
     return (4.0*sumOdd+2.0*sumEven+f[0]-f[ip-1])*h/3.0;
   }
}















