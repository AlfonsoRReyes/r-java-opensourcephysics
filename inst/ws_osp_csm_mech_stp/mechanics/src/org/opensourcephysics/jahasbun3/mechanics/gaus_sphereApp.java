/*
@Author J E Hasbun 2007.
plots the gravitational field for a sphere of mass M and inner radius a and
outer radius b.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;
import java.awt.Color;

public class gaus_sphereApp implements Calculation {
   PlotFrame plot=new PlotFrame
       ("r(m)","-g/(GM)","Gravitational Field From sphere with radii a, b");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double g[], r[];
   double a,b;
   double rmax,r0,dr;
   int NPTS;
   Arrow [] arrows=new Arrow[2];

   public gaus_sphereApp(){
      plot.setConnected(true);// set all plot's connected
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
      arrows[0] = new Arrow(0.,0.,0.,0.);
      arrows[1] = new Arrow(0.,0.,0.,0.);
   }

   public void calculate() {
     myControl.clearMessages();
     a=myControl.getDouble   ("a");
     b=myControl.getDouble   ("b");
     r0=myControl.getDouble  ("r0");
     rmax=myControl.getDouble("rmax");
     NPTS=myControl.getInt   ("NPTS");
     dr=(rmax-r0)/NPTS;
     myControl.setValue("dr((rmax-r0)/NPTS)=",dr);
     r     = new double[NPTS];
     g     = new double[NPTS];
     for (int i=0; i<=NPTS-1;i++){
       r[i]=r0+i*dr;
       g[i]=0;
       if ((r[i]>=a) & (r[i]<b)){
         g[i]=(r[i]-a)/(b-a)/r[i]/r[i];
       } else if (r[i]>=b){
         g[i] = 1./r[i]/r[i];
       }
     }
     double gmax=ArrayLib.max(g);
     arrows[0].setXY(a,0);
     arrows[0].setYlength(gmax/2.);
     arrows[0].setHeadSize(0.1f);
     arrows[0].setColor(Color.blue);
     arrows[1].setXY(b,0);
     arrows[1].setYlength(gmax/2.);
     arrows[1].setHeadSize(0.1f);
     arrows[1].setColor(Color.green);
     plot.setLineColor(0,Color.red);
     plot.setMarkerSize(0,0);
     plot.append(0, r, g);
     plot.addDrawable(arrows[0]);
     plot.addDrawable(arrows[1]);
     plot.render();
     myControl.println("Region I is to the left of the blue line.");
     myControl.println("Region II is between the blue and green lines.");
     myControl.println("Region III is to the right of the green line");
     }

   public void clear  () {
     plot.clearDrawables();
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("plots the gravitational field for a sphere of mass M");
     myControl.println ("and inner radius a and outer radius b. We have that,");
     myControl.println ("g=0 for r<a, g=-GM(r-a)/[(b-a)r^2] for a<=r<=b, and");
     myControl.println ("g=-GM/r^2 for r>b. The plots are in units of -g/(GM)");
     r    = null;
     g    = null;
     a=3.0; b=6.0;     //main parameters
     NPTS=250; r0=0.0; rmax=2*b;
     myControl.setValue("a",a);
     myControl.setValue("b",b);
     myControl.setValue("r0",r0);
     myControl.setValue("rmax",rmax);
     myControl.setValue("NPTS",NPTS);
     dr=(rmax-r0)/NPTS;
     myControl.setValue("dr((rmax-r0)/NPTS)=",dr);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new gaus_sphereApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
