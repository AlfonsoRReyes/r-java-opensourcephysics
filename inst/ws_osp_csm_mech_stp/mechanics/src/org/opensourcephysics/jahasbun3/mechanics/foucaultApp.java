/*
@Author J E Hasbun 2007.
Plots the solution of the Foucault pendulum formulas.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class foucaultApp implements Calculation {
  PlotFrame plot0= new PlotFrame("t","x,y","Foucault Pendulum - x,y versus time");
  PlotFrame plot1= new PlotFrame("x","y","y versus x precession plot");
  NumberFormat nf = NumberFormat.getInstance();
  private Control myControl;
  double [] t, x, y;
  double x0, y0, t0, tmax, dt, w0, wf, L, lth, th, tf, w1, w2, tp, g;
  double xmin,xmax,ymin,ymax, pi, cf;
  int N;

  public foucaultApp(){
    plot0.setLocation(5,5);
    plot0.setSize(400,250);
    plot1.setLocation(5,255);
    plot1.setSize(340,340);
    nf.setMaximumFractionDigits(5);
    pi=Math.PI; cf=pi/180;
 }

   public void calculate() {
     clear();
     lth=myControl.getDouble("Latitude_angle(degrees)");
     L=myControl.getDouble("Pendulum_length(m)");
     g=myControl.getDouble("g_acceleration");
     x0=myControl.getDouble("x0(m)");
     y0=myControl.getDouble("y0(m)");
     t0=myControl.getDouble("t0(sec)");
     tmax=myControl.getDouble("tmax(sec)");
     N=myControl.getInt("N");
     dt=(tmax-t0)/N;
     myControl.setValue("dt=((tmax-t0)/N)=",dt);
     th=lth*cf;                //convert to radians
     tf=2*pi*Math.sqrt(L/g);   //Foucault pendulum period
     wf=2*pi/tf;               //Foucault pendulum frequency, same as sqrt(g/L)
     w0=2*pi/24/3600;          //Earth's rotational frequency
     w1=w0*Math.sin(th);       //precession frequency
     w2=Math.sqrt(w1*w1+wf*wf);//roughly pendulum frequency
     tp=Math.abs(2*pi/w1);     //precession period
     x   = new double[N];      //x coord
     y   = new double[N];      //y coord
     t   = new double[N];      //time
     xmin=0; xmax=0; ymin=0; ymax=0;
     for (int i = 0; i < N; i++) {
       t[i]= t0+i*dt;
       x[i]=x0*Math.cos(w1*t[i])*Math.cos(w2*t[i])+
           y0*Math.cos(w1*t[i])*Math.sin(w2*t[i]);
       y[i]=-x0*Math.sin(w1*t[i])*Math.cos(w2*t[i])+
           y0*Math.cos(w1*t[i])*Math.cos(w2*t[i]);
       if(x[i]<xmin){xmin=x[i];}
       if(y[i]<ymin){ymin=y[i];}
       if(x[i]>xmax){xmax=x[i];}
       if(y[i]>ymax){ymax=y[i];}
     }
     //System.out.println("i="+i+" x="+x[i]+" y="+y[i]);
     plot0.setPreferredMinMax(t0,tmax,Math.min(xmin,ymin),Math.max(xmax,ymax));
     plot0.setConnected(0,true);
     plot0.setMarkerSize(0, 0);//0-dataset index, use size 0 for symbol
     plot0.append(0, t, x); //colors based on DatasetManager.getLineColor(i)
     plot0.setLineColor(0,java.awt.Color.black);
     plot0.setConnected(1,true);
     plot0.setMarkerSize(1, 0);//1-dataset index, use size 0 for symbol
     plot0.append(1, t, y); //colors based on DatasetManager.getLineColor(i)
     plot0.setLineColor(1,java.awt.Color.blue);
     plot1.setConnected(0,true);
     plot1.setMarkerSize(0, 0);//1-dataset index, use size 0 for symbol
     plot1.append(0, x, y); //colors based on DatasetManager.getLineColor(i)
     plot1.setLineColor(0,java.awt.Color.red);
     plot1.setSquareAspect(true);
     myControl.println("Upper left graph - black: x(t); blue: y(t). Lower is y(x).");
     myControl.println("A Foucault pendulum traces out a path given by y(x).");
     myControl.println("Precession frequency, Period -> "+nf.format(w1*3600)+
                " rad/hr\t"+nf.format(tp/3600)+" hrs\n"+"Pendulum frequency, period -> "+
                nf.format(wf*3600)+" rad/hr\t"+nf.format(tf/3600)+" hrs");
   }

   public void clear() {
     myControl.clearMessages();
     plot0.clearData();
     plot1.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("Plots the solution of the Foucault pendulum formulas");
     myControl.println("x=x0*cos(w1*t)*cos(w2*t)+y0*cos(w1*t)*sin(w2*t) and");
     myControl.println("y=-x0*sin(w1*t)*cos(w2*t)+y0*cos(w1*t)*cos(w2*t).");
     myControl.println("The example run uses a large pedulum, whose period is");
     myControl.println("one hour so that at the north pole each swing will equal");
     myControl.println("1/24 part of a full pendulum plane precession period. Here");
     myControl.println("w1=w0*sin(latitude_angle), w0=Earth's rotational frequency");
     myControl.println("and w2=sqrt(w1^2+wf^2) with wf=sqrt(g/L).");
     t    = null;
     x    = null;
     y    = null;
     lth=34;                   //latitude angle in degrees
     tf=3600;                  //pendulum period - 1 hour for convenience
     g=9.8;                    //acceleration due to gravity
     L=Math.round((tf/2/pi)*(tf/2/pi)*g*1000)/1000;//pendulum length
     x0=1.0; y0=0.0;           //initial conditions
     t0=0; tmax=24*3600;       //time range
     N=500;                    //number of t points calculated
     myControl.setValue("Latitude_angle(degrees)",lth);
     myControl.setValue("Pendulum_length(m)",L);
     myControl.setValue("g_acceleration",g);
     myControl.setValue("x0(m)",x0);
     myControl.setValue("y0(m)",y0);
     myControl.setValue("t0(sec)",t0);
     myControl.setValue("tmax(sec)",tmax);
     myControl.setValue("N",N);
     dt=(tmax-t0)/N;
     myControl.setValue("dt=((tmax-t0)/N)=",dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new foucaultApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(375,520);
     myControl.setDividerLocation(250);
     model.setControl(myControl);
   }
}
