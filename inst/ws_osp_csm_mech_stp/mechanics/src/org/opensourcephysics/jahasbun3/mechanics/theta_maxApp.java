/*
@Author J E Hasbun 2007.
Plots the maximum scattering angle theta_1 versus the m2/m1 ratio.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.text.*;

public class theta_maxApp implements Calculation {
   PlotFrame plot= new PlotFrame("m_{2}/m_{1}","theta_{max}(pi)",
                                 "case of m2 < m1 theta versus m2/m1");
   NumberFormat nf = NumberFormat.getInstance();
   //DecimalFormat df= new DecimalFormat("+0.00000E00");
   private Control myControl;
   double [] x,y;
   double dx,pi,x0,xmax;
   int NPTS;

   public theta_maxApp(){
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
      pi=Math.PI;
   }

   public void calculate() {
     myControl.clearMessages ();
     x0=myControl.getDouble  ("x0");
     xmax=myControl.getDouble("xmax");
     if(x0<0.0 || xmax>1.0){
       myControl.println("Improper Range");
       resetCalculation();
     }
     NPTS=myControl.getInt   ("NPTS");
     dx=(xmax-x0)/NPTS;
     myControl.setValue      ("dx=((xmax-x0)/NPTS)=",dx);
     x   = new double[NPTS];
     y   = new double[NPTS];
     for (int i = 0; i < NPTS; i++) {
     x[i]=x0+i*dx;
     y[i]=Math.acos(Math.sqrt(1-x[i]*x[i]))/pi;
     }
     plot.setPreferredMinMax(0,1,0,0.5);
     plot.setConnected(0,true);
     plot.setMarkerSize(0, 0);
     plot.append(0, x, y);
     plot.setLineColor(0,java.awt.Color.blue);
     myControl.println("red: atom1 \nblack: atom2 \ndots:average positions");
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("Plots the maximum scattering angle theta_1 versus the m2/m1");
     myControl.println("ratio as given by y=acos(sqrt(1-x*x))/pi in units of pi.");
     x    = null;
     y    = null;
     x0=0;
     xmax=1;
     NPTS=200;
     myControl.setValue("x0",x0);
     myControl.setValue("xmax",xmax);
     myControl.setValue("NPTS",NPTS);
     dx=(xmax-x0)/NPTS;
     myControl.setValue("dx=((xmax-x0)/NPTS)=",dx);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new theta_maxApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,350);
     myControl.setDividerLocation(150);
     model.setControl(myControl);
   }
}
