/*
@Author J E Hasbun 2007.
Plots the coordinate solutions for the full bimodal coupled spring mass
system.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
//import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.numerics.*;
import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class inter_spr2App implements Calculation {
   PlotFrame plot=new PlotFrame
       ("time","Position","Coupled Spring-Mass Bimodal System ");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double x1[], x2[], tt[];
   double m,k0,k,x10,x20,xs,xd,w1,w2,wmin;
   double tau,t0,tmax,dt;
   int NPTS;

   public inter_spr2App(){
      //plot.setVisible(true);
      plot.setConnected(true);// set all plot's connected
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     m=myControl.getDouble("m");
     k0=myControl.getDouble("k0");
     k=myControl.getDouble("k");
     x10=myControl.getDouble("x10");
     x20=myControl.getDouble("x20");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     w1=Math.sqrt(k/m); w2=Math.sqrt((k+2*k0)/m);
     x1     = new double[NPTS];
     x2     = new double[NPTS];
     tt     = new double[NPTS];
     xs=(x10+x20)/2;
     xd=(x10-x20)/2;
     for (int i=0; i<=NPTS-1;i++){
       tt[i]=t0+i*dt;
       x1[i]=xs*Math.cos(w1*tt[i])+xd*Math.cos(w2*tt[i]); //mass 1 position
       x2[i]=xs*Math.cos(w1*tt[i])-xd*Math.cos(w2*tt[i]); //mass 2 position
       /*can also write x1 and x2 as below
       x1[i]=x10*Math.cos((w1+w2)*tt[i]/2)*Math.cos((w2-w1)*tt[i]/2)
       +x20*Math.sin((w1+w2)*tt[i]/2)*Math.sin((w2-w1)*tt[i]/2);
       x2[i]=x10*Math.sin((w1+w2)*tt[i]/2)*Math.sin((w2-w1)*tt[i]/2)
       +x20*Math.cos((w1+w2)*tt[i]/2)*Math.cos((w2-w1)*tt[i]/2);
       */
     }
     myControl.println("blue: mass 1, red: mass 2");
     myControl.println("w1="+nf.format(w1)+", w2="+nf.format(w2));
     plot.setLineColor(0,Color.blue);
     plot.setMarkerSize(0,0);
     plot.append(0, tt, x1);
     plot.setLineColor(1,Color.red);
     plot.setMarkerSize(1,0);
     plot.append(1, tt, x2);
     }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the coordinate solutions for the full bimodal");
     myControl.println ("coupled spring mass system.");
     myControl.println ("k0=middle spring constant, k=ends spring constant,");
     myControl.println ("also m=mass=m1=m2 and the frequencies are w1=sqrt(k/m),");
     myControl.println ("w2=sqrt((k+2*k0)/m). Here x1=xs*cos(w1*t)+xd*cos(w2*t),");
     myControl.println ("x2=xs*cos(w1*t)-xd*cos(w2*t), where xs=(x10+x20)/2 and");
     myControl.println ("xd=(x10-x20)/2.");
     myControl.println ("Initial mass positions are x10, x20");
     x1    = null;
     x2    = null;
     tt    = null;
     m=1; k0=1.0; k=10; x10=1; x20=0;
     w1=Math.sqrt(k/m); w2=Math.sqrt((k+2*k0)/m);
     wmin=Math.min(w1,w2); tau=2*Math.PI/wmin;
     NPTS=300; t0=0; tmax=20*tau;
     myControl.setValue("m",m);
     myControl.setValue("k0",k0);
     myControl.setValue("k",k);
     myControl.setValue("x10",x10);
     myControl.setValue("x20",x20);
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("NPTS",NPTS);
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new inter_spr2App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
