/*
@Author J E Hasbun 2007.
Shows the relationship between A1 and the initial angle for the non-linear
approximation of the pendulum's sin(theta) term.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
//import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class pend0App implements Calculation {
   PlotFrame plot=new PlotFrame
  ("Theta_{0} (degrees)","A_{1} (degrees)","Relation between A_1 & Theta_0");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double A1[], th[];
   double tol,w0,a3,th0,thmax,dth,x_initial,x_final,thofi;
   int NPTS,maxiter;
   Function fx=new fofx();
   Bisnewt bisnewt=new Bisnewt();

   public pend0App(){
     maxiter=10;
     tol=1.e-5;
      //plot.setVisible(true);
      plot.setConnected(1,true);// set all plot's connected
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     w0=myControl.getDouble("w0");
     th0=myControl.getDouble("th0");
     thmax=myControl.getDouble("thmax");
     maxiter=myControl.getInt("maxiter");
     NPTS=myControl.getInt("NPTS");
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth((tmax-th0)/NPTS)=",dth);
     A1     = new double[NPTS];
     th     = new double[NPTS];
     for (int i=0; i<=NPTS-1;i++){
       th[i]=th0+i*dth;
       x_initial=th[i]-4; x_final=th[i]+4; //root interval
       thofi=th[i]; // variable used in fx=fofx
       A1[i]=bisnewt.Bisnewt(x_initial, x_final, maxiter, tol,fx);
     }
     myControl.println("Red line is A1=Theta_0 approximation");
     myControl.println("Blue dots is A1 as the solution to:");
     myControl.println("f = A1+a3*A1^2/(27*a3*A1^2-32*w0^2)-theta0=0");
     plot.setLineColor(0,Color.blue);
     plot.setMarkerSize(0,1);
     plot.append(0, th, A1); //A1 self-consistent solution
     plot.setLineColor(1,Color.red);
     plot.setMarkerSize(1,0);
     plot.append(1, th, th); //A1=theta line
     }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Shows the relationship between A1 and the initial angle");
     myControl.println ("for the non-linear approximation of the pendulum's ");
     myControl.println ("sin(theta) term. Here A1 is the solution to the equation:");
     myControl.println ("f = A1+a3*A1^2/(27*a3*A1^2-32*w0^2)-theta0=0, where");
     myControl.println ("theta0 is the intial angle");

     A1    = null;
     th    = null;
     w0=1;
     a3=w0*w0/6;
     maxiter=10;
     th0=1.e-3; thmax=90;
     NPTS=25;
     myControl.setValue("w0",w0);
     myControl.setValue("th0",th0);
     myControl.setValue("thmax",thmax);
     myControl.setValue("maxiter",maxiter);
     myControl.setValue("NPTS",NPTS);
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth((tmax-th0)/NPTS)=",dth);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public class fofx implements Function {
     public double evaluate(double x) {
       //defines f=0 here x=A1 is sought, th0fi=initial angle defined before
       double f;
       f = x+a3*x*x*x/(27*a3*x*x-32*w0*w0)-thofi;
       return f;
     }
   }

   public static void main(String[] args) {
     Calculation model = new pend0App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
