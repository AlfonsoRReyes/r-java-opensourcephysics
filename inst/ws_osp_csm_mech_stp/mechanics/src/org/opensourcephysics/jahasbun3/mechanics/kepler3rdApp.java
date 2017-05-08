/*
@Author J E Hasbun 2007.
This program looks at the period-semimajor axis relationshipof the planets in
the solar system, in particular, Kepler's 3rd law.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;


public class kepler3rdApp implements Calculation {
   PlotFrame plot0= new PlotFrame("a^{3}(AU^{3})","period^{2}(yr^{2})",
                                  "period^2 vs a^3 & line fit");
   PlotFrame plot1= new PlotFrame("a","period","period vs a^(3/2) semilog plot");
   Function curvefit;
   private Control myControl;
   String [] planet={"Mercury ","Venus   ","Earth   ","Mars    ","Jupiter ",
                     "Saturn  ","Uranus  ","Neptune ","Pluto   "};//
   double [] a={0.387,0.723,1.000,1.524,5.203,9.539,19.191,30.061,39.529};//semimajor axes
   double [] p={0.241,0.615,1.000,1.881,11.86,29.46,84.01,164.8,248.5};   //%period
   double [] x,y,xx,yy,x3,y3;
   int N=9,N2=100;

   public kepler3rdApp(){
      plot0.setLocation(5,5);
      plot0.setSize(400,295);
      plot1.setLocation(5,300);
      plot1.setSize(400,295);
   }

   public void calculate() {
     double xmin,xmax,ymin,ymax,dx,af=(1+0.1);
     myControl.clearMessages();
     x   = new double[N];
     y   = new double[N];
     xx  = new double [N2];
     yy  = new double [N2];
     x3  = new double [N2];
     y3  = new double [N2];
     for (int i = 0; i < N; i++) {
       x[i]=Math.pow(a[i],3);  //cube the planet's semimajor axis
       y[i]=Math.pow(p[i],2);  //square the planet's period
     }
     curvefit=CurveFitting.linearRegression(x,y);//find the linear fit function m*x+b
     xmin=ArrayLib.min(x);
     xmax=ArrayLib.max(x);
     dx=(xmax-xmin)/(N2-1);
     for (int i = 0; i < N2; i++) {
       xx[i]=xmin+i*dx;
       yy[i]=curvefit.evaluate(xx[i]);
     }
     ymin=ArrayLib.min(yy);
     ymax=ArrayLib.max(yy);
     plot0.setPreferredMinMax(xmin,xmax,ymin*af,ymax*af);
     plot0.setConnected(0,false);
     plot0.setMarkerShape(0,1);
     plot0.setMarkerSize(0,2);
     plot0.append(0,x,y);          //planet data p^2 vs a^3
     plot0.setConnected(1,true);
     plot0.setMarkerSize(1,0);
     plot0.setLineColor(1,java.awt.Color.red);
     plot0.append(1,xx,yy);        //line fit to planet data
     // next we do y3=x3^(3/2) and plot in semilog scale
     xmin=ArrayLib.min(a);
     xmax=ArrayLib.max(a);
     dx=(xmax-xmin)/(N2-1);
     for (int i = 0; i < N2; i++) {
       x3[i]=xmin+i*dx;
       y3[i]=Math.pow(x3[i],3.0/2.0);
     }
     ymin=ArrayLib.min(p);
     ymax=ArrayLib.max(p);
     plot1.setLogScale(false,true);//semilog scale
     plot1.setPreferredMinMax(xmin,xmax,ymin*af,ymax*af);
     plot1.setConnected(0,false);
     plot1.setMarkerSize(0, 1);
     plot1.append(0,a,p);          //actual planet data p versus a
     plot1.setConnected(1,true);
     plot1.setMarkerSize(1,0);
     plot1.setLineColor(1,java.awt.Color.red);
     plot1.append(1,x3,y3);        //line y3=x3^(3/2)
     myControl.println("Upper graph line fit (red) found:\n"+curvefit+"");
     myControl.println("Here the dots are planetary data as: y=period^2, x=a^3.\n");
     myControl.println("Planet Data (black dots) being used:");
     myControl.println("Planet\t a (AU)\t p (yr)");
     for (int i = 0; i < N; i++) {
       myControl.println(planet[i]+"\t"+a[i]+"\t"+p[i]);
     }
     myControl.println("\nThe Red line on the lower graph is the period versus a^(3/2) plot.");
   }

   public void clear  () {
     plot0.clearData();
     plot1.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("This program looks at the period-semimajor axis relationship");
     myControl.println("of the planets in the solar system, in particular, Kepler's");
     myControl.println("3rd law.");
     x    = null;
     y    = null;
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new kepler3rdApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,590);
     myControl.setDividerLocation(175);
     model.setControl(myControl);
   }
}
