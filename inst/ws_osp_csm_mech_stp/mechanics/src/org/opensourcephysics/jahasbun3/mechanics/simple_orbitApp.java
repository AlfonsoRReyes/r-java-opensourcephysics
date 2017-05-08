/*
@Author J E Hasbun 2007.
Plots the zero force case orbit u=C*sin(theta)=1/r.
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
import org.opensourcephysics.frames.*;
import java.awt.Color;

public class simple_orbitApp implements Calculation {
  PlotFrame plot1= new PlotFrame("theta(rad)","r(theta)","r(theta) - absence of a force");
  PlottingPanel polarPanel = new PlottingPanel("Theta","r","Zero force straigt line motion");
  DrawingFrame  polarFrame = new DrawingFrame(polarPanel);
  PolarType1 axes = new PolarType1(polarPanel);
  Trail trail=new Trail();
  private Control myControl;
  double [] th,r;
  double r0,th0,thmax,dth,pi,C;
  int NPTS, system;

  public simple_orbitApp(){
    //plot1.setConnected(true);//set all connected
    plot1.setLocation(5,5);
    plot1.setSize(400,295);
    polarFrame.setLocation(5,300);
    polarFrame.setSize(400,295);
    polarFrame.setTitle("Zero force straight line motion - polar plot");
    pi=Math.PI;
  }

   public void calculate() {
     clear();
     double xmin,xmax,ymin,ymax, af=(1+0.1);
     myControl.clearMessages();
     C=myControl.getDouble("C");
     th0=myControl.getDouble("th0");
     thmax=myControl.getDouble("thmax");
     NPTS=myControl.getInt("NPTS");
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth((thmax-th0)/NPTS)=",dth);
     r     = new double[NPTS];
     th    = new double[NPTS];
     double [] xy_coord = new double[2];
     xmin=0; xmax=0; ymin=0; ymax=0;
     for (int i=0; i< NPTS;i++){
       th[i]=th0+i*dth;
       r[i]=C/Math.sin(th[i]);  //r=C*cosecant(theta)
       xy_coord[0]=r[i]*Math.cos(th[i]);
       xy_coord[1]=r[i]*Math.sin(th[i]);
       trail.addPoint(xy_coord[0],xy_coord[1]);
       if(xy_coord[0]<xmin){xmin=xy_coord[0];}
       if(xy_coord[0]>xmax){xmax=xy_coord[0];}
       if(xy_coord[1]<ymin){ymin=xy_coord[1];}
       if(xy_coord[1]<ymax){ymin=xy_coord[1];}
     }
     plot1.setPreferredMinMaxX(th0,thmax);
     plot1.setConnected(0,true);
     plot1.setLineColor(0,Color.black);
     plot1.setMarkerSize(0,0);
     plot1.append(0,th,r); //r(theta)
     polarPanel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     polarPanel.setSquareAspect(true);
     trail.color=java.awt.Color.red;
     polarPanel.addDrawable(trail);
     axes.setDeltaR(Math.sqrt(xmax*xmax+ymax*ymax)/5);
     axes.setDeltaTheta(Math.PI/6);
     polarFrame.render();
     myControl.println("Regarding the zero-force plots");
     myControl.println("Upper: black-r(theta) in Cartesian coordinates");
     myControl.println("Lower: red-r(theta) in polar coordinates");
     myControl.println("The straight line motion is due to the absence of forces.");
     }

   public void clear () {
     plot1.clearData();
     trail.clear();
     polarPanel.clear();
     polarFrame.render();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the zero force case orbit u=C*sin(theta)=1/r.");
     r    = null;
     th   = null;
     C=1.5;
     th0=0.3;
     thmax=Math.ceil((pi-th0)*100)/100;
     r0=1.0;          //initial position
     NPTS=200;
     myControl.setValue("C",C);
     myControl.setValue("th0",th0);
     myControl.setValue("thmax",thmax);
     myControl.setValue("NPTS",NPTS);
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth((thmax-th0)/NPTS)=",dth);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new simple_orbitApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(420, 5);
     myControl.setSize(300,590);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}

