/*
@Author J E Hasbun 2007.
Plots the coordinates of the atoms of a free falling molecule.
The molecular potential model used is 1/x^4-1/x^3.
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
import java.text.*;

public class molec_muApp implements Calculation {
   PlotFrame plot= new PlotFrame("x","y(a_{b})","atomic positions as molecule free falls");
   NumberFormat nf = NumberFormat.getInstance();
   //DecimalFormat df= new DecimalFormat("+0.00000E00");
   private Control myControl;
   double [] x1,x2,ycm,x,x1a,x2a;
   double v0y,y0,g,tmax,xb,xi,A1,xcm,pi,dt,t0,t;
   int NPTS;


   public molec_muApp(){
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
      pi=Math.PI;
   }

   public void calculate() {
     myControl.clearMessages ();
     y0=myControl.getDouble  ("y0");
     v0y=myControl.getDouble ("v0y");
     xi=myControl.getDouble  ("initial_position");
     xb=myControl.getDouble  ("bond length");
     g=myControl.getDouble   ("gravity");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt   ("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue      ("dt=((tmax-t0)/Nx)=",dt);
     x1   = new double[NPTS];
     x2   = new double[NPTS];
     x1a  = new double[NPTS];
     x2a  = new double[NPTS];
     x    = new double[NPTS];
     ycm  = new double[NPTS];
     A1=9*(-1+Math.sqrt(1+32*xi/9))/16;           //A1, at t=0 x-xb=xi
     for (int i = 0; i < NPTS; i++) {
     t=t0+i*dt;
     x[i]=xb+A1*Math.cos(2*pi*t)+
         4*A1*A1*(3-Math.cos(4*pi*t))/9; //relative coordinate
     ycm[i]=y0+v0y*t-0.5*g*t*t;          //y-center of mass position
     x1[i]=xcm-0.5*x[i];                 //coordinates of the atoms
     x2[i]=xcm+0.5*x[i];
     x1a[i]=-xb/2;  x2a[i]=xb/2;         //average atoms positions
     }
     plot.setPreferredMinMax(-A1-xb/2,xb/2+A1,ArrayLib.min(ycm),y0);
     plot.setConnected(0,true);
     plot.setMarkerSize(0, 0);//0-dataset index, use size 0 for symbol
     plot.append(0, x1, ycm); //colors based on DatasetManager.getLineColor(i)
     plot.setLineColor(0,java.awt.Color.red);
     plot.setConnected(1,true);
     plot.setMarkerSize(1, 0);//1-dataset index, use size 0 for symbol
     plot.append(1, x2, ycm); //colors based on DatasetManager.getLineColor(i)
     plot.setLineColor(1,java.awt.Color.black);
     plot.append(2, x1a, ycm); //colors based on DatasetManager.getLineColor(i)
     plot.setMarkerColor(2,java.awt.Color.blue);
     plot.setMarkerSize(2,1);
     plot.append(2, x2a, ycm); //colors based on DatasetManager.getLineColor(i)
     myControl.println("red: atom1 \nblack: atom2 \ndots:average positions");
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("Plots the coordinates of the atoms of a free falling molecule.");
     myControl.println("The molecular potential model used is 1/x^4-1/x^3. The method");
     myControl.println("of successive approximations solution used is as follows:");
     myControl.println("x=xb+A1*cos(2*pi*t)+4*A1^2*(3-cos(4*pi*t))/9, for the relative");
     myControl.println("coordinate. Each atomic coordinate is given by x1=xcm-0.5*x");
     myControl.println("and x2=xcm+0.5*x. Here A1=9*(-1+sqrt(1+32*xi/9))/16 and xb=bond");
     myControl.println("length. xi=initial position from equilibrium. Dimensionles units");
     myControl.println("are used (see text).");
     x1    = null;
     x2    = null;
     x1a   = null;
     x2a   = null;
     x     = null;
     ycm   = null;
     xcm=0.0;                  //x center of mass position
     v0y=0.0; y0=20.0; g=9.8;  //initial y-velocity, and y-position, g=gravity
     t0=0;                     //initial time
     tmax=v0y/g+Math.sqrt((v0y/g)*(v0y/g)+2*y0/g); //time to reach ground used for tmax
     xb=3/2;         //equilibrium relative coordinate (molecule bond length)
     xi=0.2;         //initial position measured from equilibrium
     A1=9*(-1+Math.sqrt(1+32*xi/9))/16;           //A1, at t=0 x-xb=xi
     NPTS=100;
     myControl.setValue("y0",y0);
     myControl.setValue("v0y",v0y);
     myControl.setValue("initial_position",xi);
     myControl.setValue("bond length",xb);
     myControl.setValue("gravity",g);
     myControl.setValue("tmax",nf.format(tmax));
     myControl.setValue("NPTS",NPTS);
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt=((tmax-t0)/Nx)=",nf.format(dt));
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new molec_muApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
