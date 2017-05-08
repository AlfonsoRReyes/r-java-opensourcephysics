/*
 @Author J E Hasbun 2007
 Given the particles masses, initial positions, velocities, and accelerations,
 the linear and angular momenta, energies, forces, and torques are calculated.
 @Copyright (c) 2007
 This software is to support Intermediate Classical Mechanics
 with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
 http://www.opensourcephysics.org under the terms of the GNU General Public
 License (GPL) as published by the Free Software Foundation.
 */

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
//import java.text.NumberFormat;

public class particle2App implements Calculation {
  MatrixOps  mops;
  MatPrint   mprint=new MatPrint();
  VectorMath vm;
  private Control myControl;
  double M,E;
  double [] m,rcm,vcm,acm,e,P1,P2,L,F1,F2,Tau;
  int NumPart, Ndim;
  double [][] r,v,a,p,l,f,tau;

  public particle2App (){

  }

  public void calculate() {
    myControl.clearMessages();
    m=(double [])myControl.getObject("m");
    r=(double [][])myControl.getObject("r");
    v=(double [][])myControl.getObject("v");
    a=(double [][])myControl.getObject("a");
    NumPart=r.length;       //number of rows => number of particles
    Ndim=r[0].length;       //number of columns => dimensions
    //myControl.println("NumPart="+NumPart+", Ndim="+Ndim);
    mprint.ContPrintMat1D(m,NumPart,"particles masses: m");
    mprint.ContPrintMat2D(r,NumPart,Ndim,"particles positions: r");
    mprint.ContPrintMat2D(v,NumPart,Ndim,"particles velocities: v");
    mprint.ContPrintMat2D(a,NumPart,Ndim,"particles accelerations: a");
    //total quantities with components
    rcm=new double[Ndim];
    vcm=new double[Ndim];
    acm=new double[Ndim];
    P1=new double[Ndim];
    P2=new double[Ndim];
    L=new double[Ndim];
    F1=new double[Ndim];
    F2=new double[Ndim];
    Tau=new double[Ndim];
    e=new double[NumPart];
    //individual particles quantities components
    p=new double[NumPart][Ndim];
    l=new double[NumPart][Ndim];
    f=new double[NumPart][Ndim];
    tau=new double[NumPart][Ndim];
    //linear momentum calculation
    M=0;
    for (int i=0; i<NumPart; i++){
      M=M+m[i];                      //total mass
      for (int j=0; j<Ndim; j++){
        p[i][j]=m[i]*v[i][j];        //jth momentum component of the ith particle
        f[i][j]=m[i]*a[i][j];        //jth force component of the ith particle
      }
    }
    //angular momentum calculation
    for (int i=0; i<NumPart; i++){
      l[i]=vm.cross3D(r[i],p[i]);    //angular momentum of each particle
      e[i]=vm.dot(p[i],p[i])/2/m[i]; //energy of each particle
      tau[i]=vm.cross3D(r[i],f[i]);  //torque on each particle
    }
    mprint.ContPrintMat2D(p,NumPart,Ndim,"particles momenta: p");
    mprint.ContPrintMat2D(f,NumPart,Ndim,"particles forces: f");
    mprint.ContPrintMat2D(l,NumPart,Ndim,"particles ang. momenta: l");
    mprint.ContPrintMat2D(tau,NumPart,Ndim,"particles torques: tau");
    mprint.ContPrintMat1D(e,NumPart,"particles energies: e");
    //center of mass values, total quantities, the sum way
    E=ArrayLib.sum(e);//total energy
    for (int i=0; i<Ndim; i++){
      rcm[i]=0;
      vcm[i]=0;
      acm[i]=0;
      P2[i]=0.;
      F2[i]=0.;
      L[i]=0.;
      Tau[i]=0.;
          for(int j=0; j<NumPart; j++){
            rcm[i]=rcm[i]+m[j]*r[j][i]/M;
            vcm[i]=vcm[i]+m[j]*v[j][i]/M;
            acm[i]=acm[i]+m[j]*a[j][i]/M;
            P2[i]=P2[i]+p[j][i];     //total ith Momentum component
            F2[i]=F2[i]+f[j][i];     //total ith Momentum component
            L[i]=L[i]+l[j][i];       //total ith Ang. Moment. component
            Tau[i]=Tau[i]+tau[j][i]; //total ith torque component
          }
    }
    //the center of mass way for the total momentum and the force
    for (int i=0; i<Ndim; i++){
      P1[i]=M*vcm[i];
      F1[i]=M*acm[i];
    }
    mprint.ContPrintMat1D(P1,Ndim,"first way total momentum: P1");
    mprint.ContPrintMat1D(P2,Ndim,"second way total momentum: P2");
    mprint.ContPrintMat1D(F1,Ndim,"first way total force: F1");
    mprint.ContPrintMat1D(F2,Ndim,"second way total force: F2");
    mprint.ContPrintMat1D(L,Ndim,"total ang. momentum: L");
    mprint.ContPrintMat1D(Tau,Ndim,"total torque: Tau");
    myControl.println("total energy: E="+E);
  }

  public void resetCalculation() {
    clear();
    myControl.clearMessages();
    myControl.println ("Given the particles masses, initial positions, velocities,");
    myControl.println ("and accelerations, the linear and angular momenta, energies,");
    myControl.println ("forces, and torques are calculated. This is set up to do");
    myControl.println ("five particles presently. Refer to the text for formulas used.");

    myControl.setValue("m",new double[]{1.,2.,3.,1.e-10,1.e-10});
    myControl.setValue("r",new double[][]{{5.,4.,0.},{-1.,2.,0.},
                       {1.,-3.,0.},{0.,0.,0.},{0.,0.,0.}});
    myControl.setValue("v",new double[][]{{4.,0.,0.},{-1.,2.,0.},
                       {0.,-6.,0.},{0.,0.,0.},{0.,0.,0.}});
    myControl.setValue("a",new double[][]{{4.,0.,0.},{2.,0.,0.},
                       {0.,-6.,0.},{0.,0.,0.},{0.,0.,0.}});
  }

  public void clear  () {
  r=null;
  v=null;
  a=null;
  }

  public void setControl(Control control) {
    myControl = control;
    mprint.myControl=control;
    resetCalculation();
  }

   public static void main(String[] args) {
     Calculation model = new particle2App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(50, 10);
     myControl.setSize(475,475);
     myControl.setDividerLocation(135);
     model.setControl(myControl);
  }
}
