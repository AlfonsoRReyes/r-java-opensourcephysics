/*
@Author J E Hasbun 2007
Applies the Jacobi method to find the eigenvalues and eigenvectors of a
symmetric matrix
 @Copyright (c) 2007
 This software is to support Intermediate Classical Mechanics
 with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
 http://www.opensourcephysics.org under the terms of the GNU General Public
 License (GPL) as published by the Free Software Foundation.
 */

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
//import java.text.NumberFormat;

public class eigenJacobiApp  implements Calculation{
  EigenJacobi eigenjacobi=new EigenJacobi();
  MatrixOps matrixops=new MatrixOps();
  MatPrint mprint=new MatPrint();
  private Control myControl;
  int NumColRow;
  int icmax;
  double tol;
  double [][] A,DiagD,V;

  public eigenJacobiApp (){

    //some examples of symmetric matrices
    //double[][]A ={{8,-1,3,-1},{-1,6,2,0},{3,2,9,1},{-1,0,1,7}};
    //double[][]A ={{1,-1}, {-1,1}};
    //double[][]A ={{4,3,2,1},{3,4,3,2},{2,3,4,3},{1,2,3,4}};
    //double[][]A={{8,-1,3,2,21,9},{-1,10,7,3,-5,8},{3,7,-6,17,3,2},{2,3,17,1,11,12},
    //            {21,-5,3,11,15,6},{9,8,2,12,6,7}};

  }

  public void calculate() {
    myControl.clearMessages();
    tol=myControl.getDouble("tolerance");
    icmax=myControl.getInt("max_iter");
    A=(double [][])myControl.getObject("A");
    if(A[0].length!=A.length){
      System.out.println("array is not square");
    }
    NumColRow=A.length;
    V=new double [NumColRow][NumColRow];
    DiagD=new double[NumColRow][NumColRow];
    eigenjacobi.EigenJacobi(A,DiagD,V,NumColRow,tol,icmax);//input: A, output: DiagD,V
    //matrixops.printMat2D(A,NumColRow,"Input Matrix");
    //matrixops.printMat2D(DiagD,NumColRow,"Eigenvalues");
    //matrixops.printMat2D(V,NumColRow,"Eigenvectors");
    mprint.ContPrintMat2D(A,NumColRow,NumColRow,"Input Matrix");
    mprint.ContPrintMat2D(DiagD,NumColRow,NumColRow,"Eigenvalues");
    mprint.ContPrintMat2D(V,NumColRow,NumColRow,"Eigenvectors");
  }

  public void resetCalculation() {
    clear();
    myControl.clearMessages();
    myControl.println (" ");
    icmax=100;
    tol=1.e-5;
    myControl.setValue("A",new double[][]{{1,-1,0},{-1,1,0},
                       {0,0,1}});
    myControl.setValue("tolerance",tol);
    myControl.setValue("max_iter",icmax);
  }

  public void clear  () {
  A=null;
  }

  public void setControl(Control control) {
    myControl = control;
    mprint.myControl = control;
    resetCalculation();
  }

   public static void main(String[] args) {
     Calculation model = new eigenJacobiApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(50, 10);
     myControl.setSize(475,475);
     myControl.setDividerLocation(135);
     model.setControl(myControl);
  }
}
