/*
@Author J E Hasbun 2007
Performs printing to the Control window
@Copyright (c) 2007
This software is to support the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import java.text.DecimalFormat;

  public class MatPrint extends MatrixOps {
    Control myControl; //give myControl its value in the setControl method
                       //of the calling program
    DecimalFormat df= new DecimalFormat("+0.000E00");

    public void ContPrintMat1D(double A[], int dim, String name) {

      myControl.println(name);
      for(int i=0; i<dim;i++){
        myControl.print(df.format(A[i])+"   ");
      }
      myControl.println(" ");
    }

    public void ContPrintMat2D(double A[][], int rows, int cols, String name) {
      myControl.println(name);
      for(int i=0; i<rows;i++){
        for(int j = 0; j < cols; j++ ) {
          //myControl.print((float)Math.round(1000*A[i][j])/1000+"   ");
          myControl.print(df.format(A[i][j])+"   ");
        }
      myControl.println(" ");
      }
    }
  }
