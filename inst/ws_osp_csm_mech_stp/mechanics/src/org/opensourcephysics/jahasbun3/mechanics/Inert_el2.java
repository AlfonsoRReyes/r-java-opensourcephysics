/*
@Author J E Hasbun 2007
//Inertia function integrands in cartesian coords
@Copyright (c) 2007
This software is to support the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;

  class Inert_el2  {
    public static double inert_el2(int m, double x, double y, double z) {
      //if m=6 it does the volume integrand
      double ff=0.;
      if (m == 0){ff= y*y + z*z;}
      if (m == 1){ff=-x*y;}
      if (m == 2){ff=-x*z;}
      if (m == 3){ff= x*x+z*z;}
      if (m == 4){ff =-y*z;}
      if (m == 5){ff = x*x+y*y;}
      if (m == 6){ff = 1.0;}
      if ( (m < 0) || (m > 6)) {
        System.out.println(" only 7 integrands are needed ");
      }
      return ff;
    }
  }

