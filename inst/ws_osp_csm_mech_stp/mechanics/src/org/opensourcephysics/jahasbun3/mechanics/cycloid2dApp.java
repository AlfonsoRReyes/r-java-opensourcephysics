/*
@Author J E Hasbun 2007.
Plots the motion of a charged particle in the presence of an electric
field (+x direction) and a magnetic field (+z direction).
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;


public class cycloid2dApp implements Calculation {
   PlotFrame plot0= new PlotFrame("t","x,y","x,y versus time");
   PlotFrame plot1= new PlotFrame("x","y","y versus x");
   private Control myControl;
   double [] t, x, y;
   double q, m, B0, E0, v0x, v0y, f, wc, A, a, b, x0, y0, t0, tmax, dt;
   double xmin,xmax,ymin,ymax;
   int N;

   public cycloid2dApp(){
      plot0.setLocation(5,5);
      plot0.setSize(400,295);
      plot1.setLocation(5,300);
      plot1.setSize(400,295);
   }

   public void calculate() {
     myControl.clearMessages();
     m=myControl.getDouble("mass");
     q=myControl.getDouble("charge");
     x0=myControl.getDouble("x0");
     y0=myControl.getDouble("y0");
     v0x=myControl.getDouble("v0x");
     v0y=myControl.getDouble("v0y");
     B0=myControl.getDouble("B0");
     E0=myControl.getDouble("E0");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     N=myControl.getInt("N");
     dt=(tmax-t0)/N;
     myControl.setValue("dt=((tmax-t0)/N)=",dt);
     f=Math.atan(v0x/(v0y+E0/B0)); // the angle fi
     wc=q*B0/m;         //cyclotron frequency
     myControl.setValue("wc=q*B0/m",wc);
     A=-Math.sqrt(v0x*v0x+
       (v0y+E0/B0)*(v0y+E0/B0))/wc; //(pick sign for convenience) - constant A
     //a, b are constants based on initial conditions
     a=x0-A*Math.cos(f);
     b=y0+A*Math.sin(f);
     x   = new double[N];       //x coord
     y   = new double[N];       //y coord
     t   = new double[N];       //time
     xmin=0; xmax=0; ymin=0; ymax=0;
     for (int i = 0; i < N; i++) {
       t[i]= t0+i*dt;
       x[i]=A*Math.cos(wc*t[i]+f)+a;
       y[i]=-A*Math.sin(wc*t[i]+f)-E0*t[i]/B0+b;
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
     myControl.println("Upper left graph - black: x(t); blue: y(t). Lower is y(x).");
     myControl.println("E field is in the +x direction, the B field is in the +z direction.");
   }

   public void clear  () {
     plot0.clearData();
     plot1.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("Plots the motion of a charged particle in the presence of an");
     myControl.println("electric field (+x direction) and a magnetic field (+z direction).");
     myControl.println("We have: x=A*cos(wc*t+f)+a, y=-A*sin(wc*t+f)-E0*t/B0+b, where");
     myControl.println("f=atan(v0x/(v0y+E0/B0)), A=sqrt(v0x*v0x+(v0y+E0/B0)^2)/wc");
     myControl.println("a=x0-A*cos(f), b=y0+A*sin(f).Here, wc=q*B0/m is the cyclotron");
     myControl.println("frequency.");
     t    = null;
     x    = null;
     y    = null;
     m=1.67e-27;        //particle mass
     q=1.6e-19;         //particle charge
     B0=3.13e-8;        //B field
     E0=5.50e-8;        //E field
     v0x=5.0; v0y=0.0;  //initial speeds in (m/s)
     x0=0.0; y0=0.0;    //initial positions (m)
     t0=0; tmax=10;     //time range
     f=Math.atan(v0x/(v0y+E0/B0)); // the angle fi
     wc=q*B0/m;         //cyclotron frequency
     A=-Math.sqrt(v0x*v0x+
       (v0y+E0/B0)*(v0y+E0/B0))/wc; //(pick sign for convenience) - constant A
     //a, b are constants based on initial conditions
     a=x0-A*Math.cos(f);
     b=y0+A*Math.sin(f);
     N=500;             //number of t points calculated
     myControl.setValue("mass",m);
     myControl.setValue("charge",q);
     myControl.setValue("x0",x0);
     myControl.setValue("y0",y0);
     myControl.setValue("v0x",v0x);
     myControl.setValue("v0y",v0y);
     myControl.setValue("B0",B0);
     myControl.setValue("E0",E0);
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("N",N);
     dt=(tmax-t0)/N;
     myControl.setValue("dt=((tmax-t0)/N)=",dt);
     myControl.setValue("wc=q*B0/m",wc);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new cycloid2dApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(375,500);
     myControl.setDividerLocation(290);
     model.setControl(myControl);
   }
}
