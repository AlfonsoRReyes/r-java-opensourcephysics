/*
@Author J E Hasbun 2007.
Finds the angular momentum of a rigid body about an axis of rotation given the
angular speed.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import java.text.*;

public class r_energyApp implements Calculation {
  EigenJacobi eigenjacobi=new EigenJacobi();
  MatPrint mprint=new MatPrint();
  private Control myControl;
  int NumColRow=3;
  double [] L, T, w, w1, w2, w3, L1, L2, L3;
  double [][] I,Ip,V;
  double mw,h,wh1,wh2,wh3,wha,whb,whc,Ttotal, T1, T2, T3;
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public r_energyApp(){
    h=Math.PI/180.;
    nf.setMinimumFractionDigits(4);
  }

   public void calculate() {
     clear();
     I=(double [][])myControl.getObject("I-tensor");
     V=(double [][])myControl.getObject("Princ. Axes Dir. Cosines");
     Ip=(double [][])myControl.getObject("Princ. Axes Moments");
     mw=myControl.getDouble("rotational frquency");
     wh1=myControl.getDouble("frequency angle 1");
     wh2=myControl.getDouble("frequency angle 2");
     wh3=myControl.getDouble("frequency angle 3");
     w=new double[NumColRow];
     L=new double[NumColRow];
     w1=new double[NumColRow];
     w2=new double[NumColRow];
     w3=new double[NumColRow];
     L1=new double[NumColRow];
     L2=new double[NumColRow];
     L3=new double[NumColRow];
     //general rotation direction
     wha=wh1*h; whb=wh2*h; whc=wh3*h;
     w[0]=mw*Math.cos(wha);
     w[1]=mw*Math.cos(whb);
     w[2]=mw*Math.cos(whc);
     double wmag=VectorMath.magnitude(w);
     mprint.ContPrintMat1D(w,NumColRow,"rotation frequency");
     myControl.println("frequency magnitude="+nf.format(wmag));
     mprint.ContPrintMat2D(I,NumColRow,NumColRow,"Input inertia tensor");
     //angular momentum
     for (int i=0; i < NumColRow; i++){
       L[i]=0.0;
       for(int j=0; j < NumColRow; j++){
         L[i]=I[i][j]*w[j]+L[i];
       }
     }
     //total energy
     Ttotal=0.0;
     for (int i=0; i < NumColRow; i++){
         Ttotal=w[i]*L[i]/2.+Ttotal;
     }
     mprint.ContPrintMat1D(L,NumColRow,"Angular Momentum");
     myControl.println("total energy="+nf.format(Ttotal));
     //principal axes rotations way
     mprint.ContPrintMat2D(Ip,NumColRow,NumColRow,"Input Principal Axes Moments");
     w1[0]=V[0][0]*mw; w1[1]=V[1][0]*mw; w1[2]=V[2][0]*mw;
     w2[0]=V[0][1]*mw; w2[1]=V[1][1]*mw; w2[2]=V[2][1]*mw;
     w3[0]=V[0][2]*mw; w3[1]=V[1][2]*mw; w3[2]=V[2][2]*mw;
     //L1
     for (int i=0; i < NumColRow; i++){
       L1[i]=0.0;
       for(int j=0; j < NumColRow; j++){
         L1[i]=Ip[i][j]*w1[j]+L1[i];
       }
     }
     //T1
     T1=0.0;
     for (int i=0; i < NumColRow; i++){
         T1=w1[i]*L1[i]/2.+T1;
     }
     mprint.ContPrintMat1D(L1,NumColRow,"First Princ. Axis L1");
     mprint.ContPrintMat1D(w1,NumColRow,"w1 - Rotation");
     myControl.println("T1="+nf.format(T1));
     //L2
     for (int i=0; i < NumColRow; i++){
       L2[i]=0.0;
       for(int j=0; j < NumColRow; j++){
         L2[i]=Ip[i][j]*w2[j]+L2[i];
       }
     }
     //T2
     T2=0.0;
     for (int i=0; i < NumColRow; i++){
         T2=w2[i]*L2[i]/2.+T2;
     }
     mprint.ContPrintMat1D(L2,NumColRow,"Second Princ. Axis L2");
     mprint.ContPrintMat1D(w2,NumColRow,"w2 - Rotation");
     myControl.println("T2="+nf.format(T2));
     //L3
     for (int i=0; i < NumColRow; i++){
       L3[i]=0.0;
       for(int j=0; j < NumColRow; j++){
         L3[i]=Ip[i][j]*w3[j]+L3[i];
       }
     }
     //T3
     T3=0.0;
     for (int i=0; i < NumColRow; i++){
         T3=w3[i]*L3[i]/2.+T3;
     }
     mprint.ContPrintMat1D(L3,NumColRow,"Third Princ. Axis L3");
     mprint.ContPrintMat1D(w3,NumColRow,"w3 - Rotation");
     myControl.println("T3="+nf.format(T3));
   }

   public void clear  () {
     myControl.clearMessages();
   }

   public void resetCalculation() {
     clear();
     myControl.println ("Finds the angular momentum of a rigid body about an");
     myControl.println ("axis of rotation given the angular speed. Can also find");
     myControl.println ("L, and Kinetic energies for rotations about the principal axes,");
     myControl.println ("given their inertia moments and directions. For general rotations,");
     myControl.println ("we use L={I}w, and Trot=w'{I}w/2. For rotations about the principal");
     myControl.println ("axes, the inertia tensor is diagonal, and there are three");
     myControl.println ("possible rotations according to the direction cosines. With");
     myControl.println ("this understanding, the above equations are repeated for L,");
     myControl.println ("and Trot; thus obtaining three possible L's, and T's.");
     mw=Math.sqrt(2);
     wh1=45.; wh2=45.; wh3=90.;
     //inertia tensor about the origin
     double Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz;
     Ixx=1./3.; Ixy=-1./4.; Ixz=0.0; Iyy=1./3.; Iyz=0.0; Izz=2./3.;
     Iyx=Ixy; Izx=Ixz; Izy=Iyz;//use symmetry for the rest of the tensor elements
     double V11,V12,V13,V21,V22,V23,V31,V32,V33;
     V11=Math.cos(135*h); V12=Math.cos(135*h); V13=Math.cos(90*h);
     V21=Math.cos(135*h); V22=Math.cos(45*h); V23=Math.cos(90*h);
     V31=Math.cos(90*h); V32=Math.cos(90*h); V33=Math.cos(00*h);
     double p11=1./12, p22=7./12., p33=2./3.;
     //Inertia tensor
     myControl.setValue("I-tensor",new double[][]{{Ixx,Ixy,Ixz},{Iyx,Iyy,Iyz},
                        {Izx,Izy,Izz}});
     myControl.setValue("Princ. Axes Dir. Cosines",new double[][]{{V11,V12,V13},
                        {V21,V22,V23},{V31,V32,V33}});
     myControl.setValue("Princ. Axes Moments",new double[][]{{p11,0.0,0.0},
                        {0.0,p22,0.0},{0.0,0.0,p33}});
     myControl.setValue("rotational frquency",nf.format(mw));
     myControl.setValue("frequency angle 1",wh1);
     myControl.setValue("frequency angle 2",wh2);
     myControl.setValue("frequency angle 3",wh3);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
     mprint.myControl=control;
   }

   public static void main(String[] args) {
     Calculation model = new r_energyApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(50, 5);
     myControl.setSize(500,590);
     myControl.setDividerLocation(225);
     model.setControl(myControl);
   }
}