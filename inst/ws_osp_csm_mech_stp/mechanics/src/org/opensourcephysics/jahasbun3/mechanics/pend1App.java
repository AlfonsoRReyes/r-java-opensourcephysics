/*
@Author J E Hasbun 2007.
This program compares the period of a pendulum with the next approximations
beyond the linear simple pendulum formula.
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

public class pend1App implements Calculation {
   PlotFrame plot=new PlotFrame
  ("Period (tau_0)","Theta_{0} (degrees)","Comparison of the Periods");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double y0[], y1[], y2[], th[];
   double thmax, thmin, dth, a, b, th0;
   int NPTS,NINT;
   Function fint=new fintegrand();

   public pend1App(){
      //plot.setVisible(true);
      plot.setConnected(1,true);
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     thmin=myControl.getDouble("thmin");
     thmax=myControl.getDouble("thmax");
     NPTS=myControl.getInt("NPTS");
     dth=(thmax-thmin)/NPTS;
     myControl.setValue("dth((thmax-thmin)/NPTS)=",dth);
     NINT=myControl.getInt("N_elliptic_integral_even_points");
     if( NINT%2 > 0){       //make sure NINT is even - use modulus operator
       NINT = NINT + 1;
       myControl.setValue("N_elliptic_integral_even_points",NINT);
     }
     y0     = new double[NPTS];
     y1     = new double[NPTS];
     y2     = new double[NPTS];
     th     = new double[NPTS];
     a=0.0; b=Math.PI/2; //elliptic integral limits
     for (int i=0; i<=NPTS-1;i++){
       th[i] = thmin + i * dth;
       th0 = th[i]*Math.PI/180; // variable used in fx=fofx in radians
       y0[i]=1.0; // the small angle approximation
       y1[i]=1./(Math.sqrt(1-Math.pow(th[i]*Math.PI/180,2)/8));
       y2[i] =2*Integral.simpson(fint,a, b,NINT)/Math.PI;
     }
     myControl.println("Pendulum period calculated in different ways as follows.");
     myControl.println("Red line -  elliptic integral full solution");
     myControl.println("Blue dots - approximation: 1/(sqrt(1-theta0^2/8))");
     myControl.println("Black dots: small angle approximation here it is unity.");
     plot.setLineColor(0,Color.blue);
     plot.setMarkerSize(0,1);
     plot.append(0, th, y1); //non-linear approximation
     plot.setLineColor(1,Color.red);
     plot.setMarkerSize(1,0);
     plot.append(1, th, y2); //The elliptic integral full solution
     plot.setLineColor(2,Color.black);
     plot.setMarkerSize(2,1);
     plot.setMarkerColor(2,Color.black);
     plot.append(2, th, y0); //The small angle solution
     }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("This program compares the period of a pendulum with the");
     myControl.println ("next approximations beyond the linear simple pendulum formula");
     myControl.println ("The small angle approximation is: tau0=2*pi*sqrt(L/g).");
     myControl.println ("The nonlinear approximation uses:1./(sqrt(1-theta_0^/8)), where");
     myControl.println ("theta_0 is the initial angle. The full solution is the elliptic integral:");
     myControl.println ("2*integral of(1/sqrt(1-(sin(theta_0/2)*sin(x))^2))/pi over x on [0,pi/2].");
     myControl.println ("The periods are plotted in units of tau0.");
     y0    = null;
     y1    = null;
     y2    = null;
     th    = null;
     thmin=0; thmax=90;
     NPTS=25;
     NINT=50;//elliptic integral integration points
     myControl.setValue("thmin",thmin);
     myControl.setValue("thmax",thmax);
     myControl.setValue("NPTS",NPTS);
     dth=(thmax-thmin)/NPTS;
     myControl.setValue("dth((thmax-thmin)/NPTS)=",dth);
     myControl.setValue("N_elliptic_integral_even_points",NINT);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public class fintegrand implements Function {
     public double evaluate(double x) {
       //defines the integrand associated with an elliptic integral
       //th0 is an initial angle defined before
       double f;
       f = 1/Math.sqrt(1-Math.pow(Math.sin(th0/2)*Math.sin(x),2));
       return f;
     }
   }

   public static void main(String[] args) {
     Calculation model = new pend1App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
