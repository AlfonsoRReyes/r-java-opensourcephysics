/*
@Author J E Hasbun 2007.
Plots the amplitude of the driven HO solution.
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


public class drive_solnApp implements Calculation {
   PlotFrame plot1= new PlotFrame("t","Xf, F","Forced Solution(black) and Force(blue)");
   PlotFrame plot2= new PlotFrame("t","Xt=Xh+Xf","Total Xh+Xf solution");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double xf[],FF[],tt[],xh[],xt[];
   double x0,v0,t0,m,k,w0,theta,dt,tau,tmax,gamma,desc1;
   double delta,wd,F0,c,desc2,den,A,phase,B,th,wu,xe;
   int NPTS;

   public drive_solnApp(){
      //plot1 props
      plot1.setConnected (0,true);
      plot1.setMarkerSize(0, 0);
      plot1.setLineColor (0,Color.black);
      plot1.setConnected (1,true);
      plot1.setLineColor (1,Color.blue);
      plot1.setMarkerSize(1, 0);
      plot1.setLocation(5, 5);
      plot1.setSize(400, 290);
      //plot2 props
      plot2.setConnected(0,true);
      plot2.setLineColor(0,Color.red);
      plot2.setMarkerSize(0, 0);
      plot2.setLocation(5, 295);
      plot2.setSize(400, 290);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     x0=myControl.getDouble("x0");
     v0=myControl.getDouble("v0");
     k=myControl.getDouble("k");
     m=myControl.getDouble("m");
     F0=myControl.getDouble("F0");
     wd=myControl.getDouble("wd");
     c=myControl.getDouble("c");
     theta=myControl.getDouble("theta");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     gamma=c/m/2;
     w0=Math.sqrt(k/m);
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     myControl.setValue("gamma(c/m/2)=",gamma);
     myControl.setValue("w0",w0);
     desc2 = (2 * gamma * wd) * (2 * gamma * wd) +
           (w0 * w0 - wd * wd) * (w0 * w0 - wd * wd);
     A = F0 / m / Math.sqrt(desc2); //the driven HO amplitude
     den = w0 * w0 - wd * wd;
     if (wd < w0) {
       phase = Math.atan(2 * gamma * wd / den); //The phase difference
     }
     else {
      //shift by pi needed if w >w0;
      phase = Math.PI + Math.atan(2 * gamma * wd / den);
     }
     myControl.println("phase="+phase);
     delta = theta - phase;
     desc1=w0*w0-gamma*gamma;
     myControl.println("desc1=w0*w0-gamma*gamma="+desc1);
     if (desc1 <= 0){
        thereisaproblem();
        } else{
          tt     = new double[NPTS+1];
          xf     = new double[NPTS+1];
          FF     = new double[NPTS+1];
          xh     = new double[NPTS+1];
          xt     = new double[NPTS+1];
          for (int i = 0; i <= NPTS; i++) {
             tt[i] = t0 + i * dt;
             xf[i] = A * Math.cos(wd * tt[i] + delta); //forced or particular solution
             FF[i] = F0 * Math.cos(wd * tt[i] + theta); //driving force
             //the underdamped homogeneous solution
             wu = Math.sqrt(desc1); //underdamped hom. soln. freq.
             th = Math.atan(wu * x0 / (v0 + gamma * x0));
             B = Math.sqrt(x0 * x0 + (v0 + gamma * x0) * (v0 + gamma * x0) / wu / wu);
             xe = B * Math.exp( -gamma * tt[i]);
             xh[i] = xe * Math.sin(wu * tt[i] + th);
             xt[i] = xf[i] + xh[i];
             //System.out.println("i="+i+" xf="+xf[i]);
          }
          plot1.append(0, tt, xf);
          plot1.append(1, tt, FF);
          plot2.append(0, tt, xt);
        }
   }

   public void clear  () {
     plot1.clearData();
     plot2.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Driven HO solution and applied force. The particular solution is");
     myControl.println ("xf(t)=A*cos(wd*t+delta), wehere The amplitude is A=F0/m/sqrt(desc2)");
     myControl.println ("where desc2=(2*gamma*wd)^2+(w0*w0- wd*wd)^2, gamma=c/2/m, and");
     myControl.println ("delta = theta - phase, with theta the driving force phase. The");
     myControl.println ("quantity 'phase' has to be shifted beyond w=w0. If wd < w0 the");
     myControl.println ("phase=atan(2*gamma*wd/(w0^2-wd^2) else we have");
     myControl.println ("phase=pi+atan(2*gamma*wd/(w0^2-wd^2)). The driving force ");
     myControl.println ("is F=F0*cos(wd*t+theta), where F0=force amplitude,");
     myControl.println ("wd=driving frequency. We have taken the underdamped");
     myControl.println ("HO solution, xh(t)=B*exp(-gamma*t)*sin(w*t+th) as the ");
     myControl.println ("homogeneous solution for which k/m > gamma^2, m=mass and");
     myControl.println ("we need desc1=w0^2-gamma^2>0. Here w0=sqrt(k/m), and");
     myControl.println ("w=sqrt(desc1); B=sqrt(x0^2+(v0+gamma*x0)^2/w^2)");
     myControl.println ("th=arctan(w*x0/(v0+gamma*x0), x0=initial position, v0=initial");
     myControl.println ("speed, k=spring const., c=damping coefficient, t0=initial");
     myControl.println ("time tmax=upper time limit, NPTS=numberof points plotted");
     myControl.println ("dt is the time step. The combined solution is Xt(t)=Xh(t)+Xf(t).");
     tt    = null;
     FF    = null;
     xf    = null;
     x0=1; v0=5.0; m=0.5; k=0.5; F0=0.5; theta=0.0;
     wd=0.1; c=0.1;
     tau=2*Math.PI/wd;
     t0=0; tmax=3*tau;
     NPTS=500;
     F0=0.5;
     myControl.setValue("x0",x0);
     myControl.setValue("v0",v0);
     myControl.setValue("k",k);
     myControl.setValue("m",m);
     myControl.setValue("F0",F0);
     myControl.setValue("wd",wd);
     myControl.setValue("c",c);
     myControl.setValue("theta",theta);
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("NPTS",NPTS);
     dt=(tmax-t0)/NPTS;
     gamma=c/m/2;
     w0=Math.sqrt(k/m);
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     myControl.setValue("gamma(c/m/2)=",gamma);
     myControl.setValue("w0",w0);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public void thereisaproblem(){
      myControl.clearValues();
      myControl.clearMessages();
      myControl.println("*** c needs to be smaller since desc1 is <= zero ***");
      myControl.println("desc1=w0^2-(c/2/m)^2="+desc1);
      myControl.println("Please reset or fix and re-calculate");
      myControl.setValue("x0",x0);
      myControl.setValue("v0",v0);
      myControl.setValue("k",k);
      myControl.setValue("m",m);
      myControl.setValue("F0",F0);
      myControl.setValue("wd",wd);
      myControl.setValue("theta",theta);
      myControl.setValue("c",c);
      myControl.setValue("t0",t0);
      myControl.setValue("tmax",tmax);
      myControl.setValue("NPTS",NPTS);
      dt=(tmax-t0)/NPTS;
      gamma=c/m/2;
      w0=Math.sqrt(k/m);
      myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
      myControl.setValue("gamma(c/m/2)=",gamma);
      myControl.setValue("w0",w0);
    }

   public static void main(String[] args) {
     Calculation model = new drive_solnApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,580);
     myControl.setDividerLocation(295);
     model.setControl(myControl);
   }
}
