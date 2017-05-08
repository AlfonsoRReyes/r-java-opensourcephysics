/*
@Author J E Hasbun 2007.
Plots the solution of a body of mass m under the action of a central force of
the form -a*r^p, where r=1/u. Two different methods are used to get the orbit
one solves r(t), theta(t); the other solves u(fi). The two solutions are compared
and demonstrated to agree with each other.
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

public class centraluApp implements Calculation {
  PlotFrame plot1= new PlotFrame("theta(rad)","r(theta)","Cartesian plot");
  PlottingPanel polarPanel = new PlottingPanel("Theta","r","Polar Plot");
  DrawingFrame  polarFrame = new DrawingFrame(polarPanel);
  PolarType1 axes = new PolarType1(polarPanel);
  Trail trail0=new Trail();
  Trail trail1=new Trail();
  private Control myControl;
  double r[], v[], th[], t[], u[], r2[], fi[];
  double[] state=new double [7];
  double r0,v0,th0,u0,uv0,a,p,m,L,t0,tmax,dt, fimin,fimax,pi,dfi;
  int NPTS, system;

  public centraluApp(){
    //plot1.setConnected(true);//set all connected
    plot1.setLocation(5,5);
    plot1.setSize(400,295);
    polarFrame.setLocation(5,300);
    polarFrame.setSize(400,295);
    polarFrame.setTitle("Polar Plot - r(theta)");
    pi=Math.PI;
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
     fimin=myControl.getDouble("hmin");
     fimax=myControl.getDouble("hmax");
     dfi=(fimax-fimin)/NPTS;
     myControl.setValue("dh=(hmax-hmin)/NPTS=",dfi);
     r     = new double[NPTS];
     v     = new double[NPTS];
     th    = new double[NPTS];
     t     = new double[NPTS];
     u     = new double[NPTS];
     r2    = new double[NPTS];
     fi    = new double[NPTS];
     double [] xy_coord = new double[2];
     xmin=0; xmax=0; ymin=0; ymax=0;
     L=m*v0*r0; //angular momentum - used in both systems of equations being solved
     //----- system 1 -------------------------
     system=1;
     //initial conditions for the r(t), theta(t0) equations - system 1
     state[0]=r0;
     state[1]=v0;
     state[2]=th0;
     state[3]=t0;
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
       trail0.addPoint(xy_coord[0],xy_coord[1]);
     }
     //----- system 2 -------------------------
     system=2;
     //initial conditions for the u(theta) equation - system 2
     u0=1/r0; uv0=-m*v0/L;
     state[4]=u0;
     state[5]=uv0;
     state[6]=fimin;
     odesolver.initialize(dfi); // step size
     for (int i=0; i< NPTS;i++){
       fi[i]=ode.getState()[6];
       odesolver.step();
       //differential equation solutions for the r,v,theta positions versus time
       u[i]=ode.getState()[4];  //u
       r2[i]=1/u[i]; //r(fi)
       fi[i]=ode.getState()[6];
       xy_coord[0]=r2[i]*Math.cos(fi[i]);
       xy_coord[1]=(1/u[i])*Math.sin(fi[i]);
       trail1.addPoint(xy_coord[0],xy_coord[1]);
       if(xy_coord[0]<xmin){xmin=xy_coord[0];}
       if(xy_coord[0]>xmax){xmax=xy_coord[0];}
       if(xy_coord[1]<ymin){ymin=xy_coord[1];}
       if(xy_coord[1]<ymax){ymin=xy_coord[1];}
     }
     plot1.setPreferredMinMaxX(0,fimax);
     plot1.setConnected(0,true);
     plot1.setLineColor(0,Color.black);
     plot1.setMarkerSize(0,0);
     plot1.append(0,th,r); //r(theta)
     plot1.setConnected(1,false);
     //plot1.setLineColor(1,Color.red);
     plot1.setMarkerSize(1,1);
     plot1.setMarkerShape(1,1);
     plot1.setMarkerColor(1,Color.red);
     plot1.append(1,fi,r2); //should be equivalent to r(theta)
     polarPanel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     //polarPanel.setSquareAspect(true);
     trail0.color=java.awt.Color.blue;
     trail1.color=java.awt.Color.red;
     polarPanel.addDrawable(trail0);
     polarPanel.addDrawable(trail1);
     axes.setDeltaR(Math.sqrt(xmax*xmax+ymax*ymax)/5);
     axes.setDeltaTheta(Math.PI/6);
     polarFrame.render();
     myControl.println("Regarding the plots");
     myControl.println("Upper: black-r(theta) and red-1/u(fi)");
     myControl.println("Lower: polar plot of blue-r(theta) and red-1/u(fi).");
     myControl.println("The two systems are indistinguishable as expected.");
     }

   public void clear () {
     plot1.clearData();
     trail0.clear();
     trail1.clear();
     polarPanel.clear();
     polarFrame.render();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the solution of a body of mass m under");
     myControl.println ("the action of a central force of the form -a*r^p,");
     myControl.println ("where r=1/u. Two equivalent systems are solved and");
     myControl.println ("compared. System 1 uses differential equations:");
     myControl.println ("dr/dt=v, dv/dt=-(a/m)*r^p+L^2/(m^2*r^3) and");
     myControl.println ("dtheta/dt=L/(m*r^2) with r0,v0, and th0 as");
     myControl.println ("init. conditions at t=t0. Here L=ang. momentum");
     myControl.println ("m*r*v0, & v0 the intial tangential speed.");
     myControl.println ("System 2 solves the differential equation:");
     myControl.println ("d^2u/dfi^2=-u-((m/L^2)/u^2)*(-a./u^p)) with initial");
     myControl.println ("conditions that at fi=0, u=1/r0, and u'=-m*v0/L.");
     r    = null;
     r2   = null;
     u    = null;
     v    = null;
     th   = null;
     fi   = null;
     t    = null;
     m=1;
     a=108; p=1;
     t0=0; tmax=1.5;
     r0=1.0;          //initial position
     v0=6.0;          //initial tangential speed
     th0=0.0;         //initial angle in radians
     fimin=0; fimax=2*pi;
     NPTS=200;
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
     myControl.setValue("hmin",fimin);
     myControl.setValue("hmax",fimax);
     dfi=(fimax-fimin)/NPTS;
     myControl.setValue("dh=(hmax-hmin)/NPTS=",dfi);
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
       if(system==1){
         //state[0,1,2,3]=r,v,theta,t
         //rates are the derivatives of the state, for example
         //rate[0]=dr/dt=v->state(1), and
         //dv/dt=-(a/m)*r^p+L^2/m^2/r^3, and dtheta/dt=L/m/r^2
         rate[0] = state[1];
         rate[1] = - (a / m) * Math.pow(state[0], p) +
             L * L / m / m / Math.pow(state[0], 3);
         rate[2] = L / m / state[0] / state[0];
         rate[3] = 1; //dt/dt=1
       }
       if(system==2){
         //to solve the u(fi) where fi is now the angle
         //let state[4,5]=u,u'
         //rate[4]=du/dfi=state[5], and
         //rate[5]=du'/dfi=-u-((m/L^2)/u^2)*(-a./u^p)
         rate[4] = state[5];
         rate[5] = -state[4] -
             ( (m / L / L) / state[4] / state[4]) * ( -a / Math.pow(state[4], p));
         rate[6] = 1; //dfi/dfi
       }
     }
   }

   public static void main(String[] args) {
     Calculation model = new centraluApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(420, 5);
     myControl.setSize(300,590);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}

