/*
@Author J E Hasbun 2007.
Plots the solutions of the two atom molecular potential model using the simple,
the nonlinear approximation, and the full solution.
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

public class molecApp implements Calculation {
   PlotFrame plot=new PlotFrame
   ("Time (tau_{0})","Position (a_{0})","Comparison of Solutions");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double x0[], x1[], x2[], tt[];
   double tmax, t0, dt, xi, xb, A1, pi;
   int NPTS;
   double[] state=new double [3];

   public molecApp(){
      //plot.setConnected(1,true);
      plot.setConnected(true);//set all connected
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
      pi=Math.PI;
   }

   public void calculate() {
     ODE ode = new Molecule();
     ODESolver odesolver = new RK4(ode);
     myControl.clearMessages();
     xb=3.0/2.0; //equilibrium position
     xi=myControl.getDouble("xi");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     x0     = new double[NPTS];
     x1     = new double[NPTS];
     x2     = new double[NPTS];
     tt     = new double[NPTS];
     state[0]=xi+xb;
     state[1]=0;
     state[2]=t0;
     x0[0]=xb+xi;
     x1[0]=xb+xi;
     x2[0]=xb+xi;
     tt[0]=t0;
     A1=9*(-1+Math.sqrt(1+32*xi/9))/16;
     odesolver.initialize(dt); // step size
     for (int i=1; i<=NPTS-1;i++){
       tt[i] = t0 + i * dt;
       x0[i]=xb+xi*Math.cos(2*pi*tt[i]); // the small angle approximation
       x1[i]=xb+A1*Math.cos(2*pi*tt[i])+4*A1*A1*(3-Math.cos(4*pi*tt[i]))/9;
                                        //nonlinear approx
       odesolver.step();
       x2[i]=ode.getState()[0];         //differential equation solution
     }
     myControl.println("Equilibrium position used: xb="+xb);
     myControl.println("Pendulum solution calculated in different ways as follows:");
     myControl.println("Small angle approximation - black line");
     myControl.println("Non-linear approximation - blue dotted line");
     myControl.println("Full differential equation solution - red line");
     plot.setLineColor(0,Color.black);
     plot.setMarkerSize(0,0);
     plot.setMarkerColor(0,Color.black);
     plot.append(0, tt, x0); //The small angle solution
     plot.setLineColor(1,Color.blue);
     plot.setMarkerSize(1,1);
     plot.append(1, tt, x1); //non-linear approximation
     plot.setLineColor(2,Color.red);
     plot.setMarkerSize(2,0);
     plot.append(2, tt, x2); //The elliptic integral full solution
     }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the solutions of the two atom molecular potential");
     myControl.println ("model using the simple, the nonlinear approximation,");
     myControl.println ("and the full solution. The small angle approximation");
     myControl.println ("is: xb+xi*cos(2*pi*t). The nonlinear approximation");
     myControl.println ("is: x1=xb+A1*cos(2*pi*t)+4*A1^2*(3-cos(4*pi*t))/9,");
     myControl.println ("where A1=9*(-1+sqrt(1+32*xi/9))/16");
     myControl.println ("The full case is the numerical solution to the differential");
     myControl.println ("equation: d^2(x)/dt^2=81*pi^2*(3/x^4-2/x^3)/8, with initial");
     myControl.println ("conditions such that at t=0, x=xb+xi, v=0.");
     x0    = null;
     x1    = null;
     x2    = null;
     tt    = null;
     t0=0; tmax=2;
     NPTS=100;
     xi=0.2; //initial position measured from equilibrium
     myControl.setValue("xi",xi);
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

    class Molecule implements ODE {
     //double[] state = new double[] {xb+xi, 0, t0};//theta,d(theta)/dt,t

     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2]=x,dx/dt,t
       rate[0]=state[1];                  //dx/dt=v
       rate[1]=81*pi*pi*(3./Math.pow(state[0],4)-2./Math.pow(state[0],3))/8;
                                          //dv/dt=81*pi^2*(3/x^4-2/x^3)/8
       rate[2]=1;                         //dt/dt=1
     }
   }

   public static void main(String[] args) {
     Calculation model = new molecApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }

}

