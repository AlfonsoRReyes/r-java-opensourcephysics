/*
@Author J E Hasbun 2007.
Plots the solutions of the pendulum using the simple, the nonlinear approximation,
and the full differential equation involving the sin(theta) term.
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

public class pend2App implements Calculation {
   PlotFrame plot=new PlotFrame
   ("Time (sec)","Amplitude (degrees)","Comparison of Solutions");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double y0[], y1[], y2[], tt[];
   double tmax, t0, dt, theta_0, w0, om, cf, tau0, A1, A3, a3;
   int NPTS;
   double[] state=new double [3];

   public pend2App(){
      //plot.setConnected(1,true);
      plot.setConnected(true);//set all connected
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
      cf=Math.PI/180; //convertion factor between degrees and radians
   }

   public void calculate() {
     ODE ode = new Pendulum();
     ODESolver odesolver = new RK4(ode);
     myControl.clearMessages();
     theta_0=myControl.getDouble("theta_0");
     w0=myControl.getDouble("w0");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     y0     = new double[NPTS];
     y1     = new double[NPTS];
     y2     = new double[NPTS];
     tt     = new double[NPTS];
     A1=theta_0*cf;
     a3=w0*w0/6;
     A3=a3*Math.pow(A1,3)/(27*a3*A1*A1-32*w0*w0);
     om=w0*Math.sqrt(1-A1*A1/8);
     state[0]=theta_0*cf;
     state[1]=0;
     state[2]=t0;
     y0[0]=theta_0;
     y1[0]=theta_0;
     y2[0]=theta_0;
     tt[0]=t0;
     odesolver.initialize(dt); // step size
     for (int i=1; i<=NPTS-1;i++){
       tt[i] = t0 + i * dt;
       y0[i]=theta_0*Math.cos(w0*tt[i]); // the small angle approximation
       y1[i]=theta_0*Math.cos(om*tt[i])+A3*Math.cos(3*om*tt[i]);//nonlinear approx
       odesolver.step();
       y2[i]=ode.getState()[0]/cf; //differential equation solution
     }
     myControl.println("Pendulum solution calculated in different ways as follows:");
     myControl.println("Small angle approximation - black line");
     myControl.println("Non-linear approximation - blue dotted line");
     myControl.println("Full differential equation solution - red line");
     plot.setLineColor(0,Color.black);
     plot.setMarkerSize(0,0);
     plot.setMarkerColor(0,Color.black);
     plot.append(0, tt, y0); //The small angle solution
     plot.setLineColor(1,Color.blue);
     plot.setMarkerSize(1,1);
     plot.append(1, tt, y1); //non-linear approximation
     plot.setLineColor(2,Color.red);
     plot.setMarkerSize(2,0);
     plot.append(2, tt, y2); //The elliptic integral full solution
     }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the solutions of the pendulum using the simple, ");
     myControl.println ("the nonlinear approximation, and the full differential");
     myControl.println ("equation involving the sin(theta) term.");
     myControl.println ("The small angle approximation is: theta_0*cos(w0*t).");
     myControl.println ("The nonlinear approximation is: theta_0*cos(om*t)+A3*cos(3*om*t),");
     myControl.println ("where om=w0*sqrt(1-A1*A1/8), A1=theta_0, and");
     myControl.println ("A3=a3*A1^3/(27*a3*A1*A1-32*w0*w0), a3=w0*w0/6.");
     myControl.println ("The full case is the numerical solution to the differential");
     myControl.println ("equation: d^2(theta)/dt^2=-w0^2*sin(theta) with initial conditions");
     myControl.println ("such that at t=0 theta=theta_0, and d(theta)/dt=0.");
     y0    = null;
     y1    = null;
     y2    = null;
     tt    = null;
     w0=1;   //frequency of the HO
     tau0=2*Math.PI/w0;
     t0=0; tmax=4*tau0;
     NPTS=100;
     theta_0=90;
     myControl.setValue("theta_0",theta_0);
     myControl.setValue("w0",w0);
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

    class Pendulum implements ODE {
     //double[] state = new double[] {thr, 0, t0};//theta,d(theta)/dt,t

     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2]=theta,d(theta)/dt,t
       rate[0]=state[1];                  //d(theta)/dt=w
       rate[1]=-w0*w0*Math.sin(state[0]); //d(w)/dt=-w0^2*sin(theta)
       rate[2]=1;                         //dt/dt=1
     }
   }

   public static void main(String[] args) {
     Calculation model = new pend2App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }

}

