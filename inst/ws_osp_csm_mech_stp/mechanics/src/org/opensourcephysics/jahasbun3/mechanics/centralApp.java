/*
@Author J E Hasbun 2007.
Plots the solution of a body of mass m under the action of a central force
of the form -a*r^p.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.display.axes.PolarType1;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.awt.Color;

public class centralApp implements Calculation {
  PlotFrame plot0= new PlotFrame("t","r,v","Central Force=a*r^p, p=1");
  PlotFrame plot1= new PlotFrame("theta(rad)","r(theta)","");
  PlotFrame plot2= new PlotFrame("t","theta","Central Force=a*r^p, p=1");
  PlottingPanel polarPanel = new PlottingPanel("Theta","r","Polar Plot");
  DrawingFrame  polarFrame = new DrawingFrame(polarPanel);
  PolarType1 axes = new PolarType1(polarPanel);
  Trail trail=new Trail();
  private Control myControl;
  double r[], v[], th[], t[];
  double[] state=new double [4];
  double r0,v0,th0,a,p,m,L,t0,tmax,dt;
  int NPTS;

  public centralApp(){
    plot0.setConnected(true);//set all connected
    plot0.setLocation(5,5);
    plot0.setSize(250,250);
    plot1.setConnected(true);//set all connected
    plot1.setLocation(5,260);
    plot1.setSize(250,250);
    plot2.setConnected(true);//set all connected
    plot2.setLocation(260,5);
    plot2.setSize(250,250);
    polarFrame.setLocation(260,260);
    polarFrame.setSize(250,250);
    polarFrame.setTitle("Polar Plot - r(theta)");
  }

   public void calculate() {
     clear();
     double xmin,xmax,ymin,ymax, af=(1+0.1);
     ODE ode = new central_der();
     ODESolver odesolver = new RK4(ode);
     myControl.clearMessages();
     m=myControl.getDouble("m");
     p=myControl.getDouble("p");
     a=myControl.getDouble("a");
     r0=myControl.getDouble("r0");
     v0=myControl.getDouble("v0");
     th0=myControl.getDouble("th0");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     r     = new double[NPTS];
     v     = new double[NPTS];
     th    = new double[NPTS];
     t     = new double[NPTS];
     double [] xy_coord = new double[2];
     //initial conditions
     L=m*v0*r0; //angular momentum
     state[0]=r0;
     state[1]=v0;
     state[2]=th0;
     state[3]=t0;
     xmin=0; xmax=0; ymin=0; ymax=0;
     odesolver.initialize(dt); // step size
     for (int i=0; i< NPTS;i++){
       t[i]=ode.getState()[3];
       odesolver.step();
       //differential equation solutions for the r,v,theta positions versus time
       r[i]=ode.getState()[0];
       v[i]=ode.getState()[1];
       th[i]=ode.getState()[2];
       xy_coord[0]=r[i]*Math.cos(th[i]);
       xy_coord[1]=r[i]*Math.sin(th[i]);
       trail.addPoint(xy_coord[0],xy_coord[1]);
       if(xy_coord[0]<xmin){xmin=xy_coord[0];}
       if(xy_coord[0]>xmax){xmax=xy_coord[0];}
       if(xy_coord[1]<ymin){ymin=xy_coord[1];}
       if(xy_coord[1]<ymax){ymin=xy_coord[1];}
     }
     plot0.setLineColor(0,Color.black);
     plot0.setMarkerSize(0,0);
     plot0.setMarkerColor(0,Color.black);
     plot0.append(0,t,r); //r(t)
     plot0.setLineColor(1,Color.blue);
     plot0.setMarkerSize(1,0);
     plot0.append(1,t,v); //v(t)
     plot1.setLineColor(0,Color.black);
     plot1.setMarkerSize(0,0);
     plot1.append(0,th,r); //r(theta)
     plot2.setLineColor(0,Color.blue);
     plot2.setMarkerSize(0,0);
     plot2.append(0,t,th); //theta(t)
     polarPanel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     //polarPanel.setSquareAspect(true);
     trail.color=java.awt.Color.red;
     polarPanel.addDrawable(trail);
     axes.setDeltaR(Math.sqrt(xmax*xmax+ymax*ymax)/5);
     axes.setDeltaTheta(Math.PI/6);
     polarFrame.render();
     myControl.println("Regarding the plots");
     myControl.println("Upper left: r(t) - black, v(t) - blue");
     myControl.println("Upper right: theta(t)");
     myControl.println("Lower left: r versus theta");
     myControl.println("Lower right: Polar plot of x(t),y(t)- orbit shape");
     }

   public void clear () {
     plot0.clearData();
     plot1.clearData();
     plot2.clearData();
     trail.clear();
     polarPanel.clear();
     polarFrame.render();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the solution of a body of mass m under");
     myControl.println ("the action of a central force of the form -a*r^p.");
     myControl.println ("The following differential equations are solved:");
     myControl.println ("dr/dt=v, dv/dt=-(a/m)*r^p+L^2/(m^2*r^3) and");
     myControl.println ("dtheta/dt=L/(m*r^2) with r0,v0, and th0 as");
     myControl.println ("init. conditions at t=t0. Here L=ang. momentum");
     myControl.println ("m*r*v0, & v0 the intial tangential speed.");
     r    = null;
     v    = null;
     th    = null;
     t    = null;
     m=1;
     a=108; p=1;
     t0=0; tmax=1.5;
     r0=1.0;          //initial position
     v0=6.0;          //initial tangential speed
     th0=0.0;         //initial angle in radians
     NPTS=500;
     myControl.setValue("m",m);
     myControl.setValue("p",p);
     myControl.setValue("a",a);
     myControl.setValue("r0",r0);
     myControl.setValue("v0",v0);
     myControl.setValue("th0",th0);
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

    class central_der implements ODE {

     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2,3]=r,v,theta,t
       //rates are the derivatives of the state, for example
       //rate[0]=dr/dt=v->state(1), and
       //dv/dt=-(a/m)*r^p+L^2/m^2/r^3, and dtheta/dt=L/m/r^2
       rate[0]=state[1];
       rate[1]=-(a/m)*Math.pow(state[0],p)+L*L/m/m/Math.pow(state[0],3);
       rate[2]=L/m/state[0]/state[0];
       rate[3]=1; //dt/dt=1
     }
   }

   public static void main(String[] args) {
     Calculation model = new centralApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(520, 5);
     myControl.setSize(275,505);
     myControl.setDividerLocation(250);
     model.setControl(myControl);
   }
}

