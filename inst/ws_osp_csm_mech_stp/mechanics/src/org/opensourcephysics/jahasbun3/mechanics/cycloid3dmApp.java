/*
@Author J E Hasbun 2007.
Solves and plots the motion of a charged particle in an electromagnetic field in 3D.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import java.awt.Color;
import org.opensourcephysics.frames.*;

public class cycloid3dmApp implements Calculation {
  DrawingPanel3D panel = new DrawingPanel3D();
  DrawingFrame3D frame = new DrawingFrame3D(panel);
  PlotFrame plot=new PlotFrame
  ("t","x,y,z","black-x(t), blue-y(t), red-z(t)");
  ElementTrail trace = new ElementTrail();
  ElementCircle point=new ElementCircle();
  private Control myControl;
  double x[], y[], z[], t[];
  double[] state=new double [7];
  double [] E,B;
  double tmax, t0, dt, q, m, qm;
  double x0,vx0,y0,vy0,z0,vz0;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int NPTS;

  public cycloid3dmApp(){
    //panel.getVisualizationHints().setZFormat(null);
    //panel.setAlignmentX(panel.LEFT_ALIGNMENT);
    panel.getVisualizationHints().setDecorationType(VisualizationHints.DECORATION_AXES);
    //frame.setVisible(true);
    frame.setLocation(5,5);
    frame.setSize(400,295);
    //plot.setConnected(1,true);
    plot.setConnected(true);//set all connected
    plot.setLocation(5,300);
    plot.setSize(400,295);
    trace.setConnected(false);
    trace.getStyle().setLineColor(java.awt.Color.red);
    point.getStyle().setFillColor(java.awt.Color.GREEN);
    point.getInteractionTarget(Element.TARGET_POSITION).setEnabled(true);
  }

   public void calculate() {
     clear();
     ODE ode = new Cycloid3d_der();
     ODESolver odesolver = new RK4(ode);
     myControl.clearMessages();
     m=myControl.getDouble("mass");
     q=myControl.getDouble("charge");
     x0=myControl.getDouble("x0");
     vx0=myControl.getDouble("vx0");
     y0=myControl.getDouble("y0");
     vy0=myControl.getDouble("vy0");
     z0=myControl.getDouble("z0");
     vz0=myControl.getDouble("vz0");
     B=(double [])myControl.getObject("B");
     E=(double [])myControl.getObject("E");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     x     = new double[NPTS];
     y     = new double[NPTS];
     z     = new double[NPTS];
     t     = new double[NPTS];
     //initial conditions
     qm=q/m;
     state[0]=x0; state[1]=vx0;
     state[2]=y0; state[3]=vy0;
     state[4]=z0; state[5]=vz0;
     state[6]=t0;
     xmin=0; xmax=0; ymin=0; ymax=0; zmin=0; zmax=0;
     odesolver.initialize(dt); // step size
     for (int i=0; i< NPTS;i++){
       t[i]=ode.getState()[6];
       odesolver.step();
       //differential equation solutions for the x,y,z positions versus time
       x[i]=ode.getState()[0];
       y[i]=ode.getState()[2];
       z[i]=ode.getState()[4];
       trace.addPoint(x[i], y[i], z[i]);
       if(x[i]<xmin){xmin=x[i];}
       if(y[i]<ymin){ymin=y[i];}
       if(z[i]<zmin){zmin=z[i];}
       if(x[i]>xmax){xmax=x[i];}
       if(y[i]>ymax){ymax=y[i];}
       if(z[i]>xmax){zmax=z[i];}
     }
     point.setXYZ(x[NPTS-1], y[NPTS-1], z[NPTS-1]);//place a point at the end of C
     point.setSizeXYZ(0.1, 0.1, 0.1);     //point size
     panel.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
     panel.addElement(point);//Interactive point at the end position
     panel.addElement(trace);
     frame.setTitle("Charge in general E&B fields - 3D");
     frame.getJFrame().setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
     plot.setLineColor(0,Color.black);
     plot.setMarkerSize(0,0);
     plot.setMarkerColor(0,Color.black);
     plot.append(0, t, x); //x
     plot.setLineColor(1,Color.blue);
     plot.setMarkerSize(1,0);
     plot.append(1, t, y); //y
     plot.setLineColor(2,Color.red);
     plot.setMarkerSize(2,0);
     plot.append(2, t, z); //z
     myControl.println("The upper left 3D plot is z versus x and y. The green interactive");
     myControl.println("point at the end of the path taken can be moved. The lower left");
     myControl.println("graph is a plot of x(t)-black, y(t)-blue, and z(t)-red.");
     }

   public void clear  () {
     plot.clearData();
     trace.clear();
     point.setXYZ(0,0,0);
     point.setSizeXYZ(0, 0, 0);
     point.getStyle().setFillColor(java.awt.Color.white);
     panel.addElement(trace);
     panel.addElement(point);
     panel.setPreferredMinMax(0,0,0,0,0,0);
     frame.render();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Solves the motion of a charged particle in an electromagnetic");
     myControl.println ("field in 3D. The differential equations used are as follows.");
     myControl.println ("x part: dx/dt=vx, dvx/dt=qm*(vy*Bz-vz*By+Ex),");
     myControl.println ("y part: dy/dt=vy, dvy/dt=qm*(vz*Bx-vx*Bz+Ey),");
     myControl.println ("z part: dz/dt=vz, dvz/dt=qm*(vx*By-vy*Bx+Ez).");
     myControl.println ("Initial conditions at t=t0, x=x0, vx=vx0, y=y0, vy=vy0.");
     myControl.println ("The fields are inputs: B=(Bx,By,Bz), E=(Ex,Ey,Ez).");
     x    = null;
     y    = null;
     z    = null;
     t    = null;
     q=1.6e-19; m=1.67e-27;
     x0=0.0; vx0=5.0;
     y0=0.0; vy0=0.0;
     z0=0.0; vz0=0.5;
     t0=0; tmax=10.0;
     NPTS=500;
     myControl.setValue("mass",m);
     myControl.setValue("charge",q);
     myControl.setValue("x0",x0);
     myControl.setValue("vx0",vx0);
     myControl.setValue("y0",y0);
     myControl.setValue("vy0",vy0);
     myControl.setValue("z0",z0);
     myControl.setValue("vz0",vz0);
     myControl.setValue("B",new double[]{1e-8,-1e-9,5.13e-8});
     myControl.setValue("E",new double[]{0.5e-8,1e-9,-3e-9});
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

    class Cycloid3d_der implements ODE {

     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2,3,4,5]=x,dx/dt,y,dy/dt,z,dz/dt
       //rates are the derivatives of the state, for example
       //rate[0]=dx/dt=vx->state(1), and
       //dvx/dt=q*(vy*Bz-vz*By+Ex)/m->dstate(1)/dt=rate[1], etc
       rate[0]=state[1];
       rate[1]=qm*(state[3]*B[2]-state[5]*B[1]+E[0]);
       rate[2]=state[3];
       rate[3]=qm*(state[5]*B[0]-state[1]*B[2]+E[1]);
       rate[4]=state[5];
       rate[5]=qm*(state[1]*B[1]-state[3]*B[0]+E[2]);
       rate[6]=1;                         //dt/dt=1
     }
   }

   public static void main(String[] args) {
     Calculation model = new cycloid3dmApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(375,550);
     myControl.setDividerLocation(300);
     model.setControl(myControl);
   }

}

