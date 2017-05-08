/*
@Author J E Hasbun 2007.
Draws an ellipse of minimum radius rmin and eccentricity e.
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
//import org.opensourcephysics.display3d.simple3d.*;
//import java.awt.Color;
import java.text.NumberFormat;


public class ellipseApp implements Calculation {
  PlottingPanel polarPanel = new PlottingPanel("Theta","r","Ellipse");
  DrawingFrame  polarFrame = new DrawingFrame(polarPanel);
  PolarType1 axes = new PolarType1(polarPanel);
  Trail trail=new Trail();
  Circle [] point=new Circle[2];
  private Control myControl;
  NumberFormat nf = NumberFormat.getInstance();
  double [] th,r;
  double rmin,a,rmax,e,thmax,dth,pi,th0,ff,f1,f2;
  int NPTS, system;

  public ellipseApp(){
    nf.setMaximumFractionDigits(3);
    polarFrame.setLocation(5,5);
    polarFrame.setSize(400,590);
    polarFrame.setTitle("Ellipse of eccentricity e");
    pi=Math.PI;
    point[0]=new Circle();
    point[1]=new Circle();
  }

   public void calculate() {
     clear();
     double xmin,xmax,ymin,ymax, af=(1+0.1);
     myControl.clearMessages();
     e=myControl.getDouble("eccentricity");
     rmin=myControl.getDouble("rmin");
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
       r[i]=rmin*(1+e)/(1+e*Math.cos(th[i]));  //ellipse in AU units
       xy_coord[0]=r[i]*Math.cos(th[i]);
       xy_coord[1]=r[i]*Math.sin(th[i]);
       trail.addPoint(xy_coord[0],xy_coord[1]);
       if(xy_coord[0]<xmin){xmin=xy_coord[0];}
       if(xy_coord[0]>xmax){xmax=xy_coord[0];}
       if(xy_coord[1]<ymin){ymin=xy_coord[1];}
       if(xy_coord[1]<ymax){ymin=xy_coord[1];}
     }
     rmax=rmin*(1+e)/(1+e*Math.cos(pi));// max distance
     a=(rmin+rmax)/2;                   // semimajor axis
     ff=2*a*e;                          //distance between focal points
     f1=0; f2=f1-ff;                    //focal points
     point[0].setXY(f1,0);               //place the first focal point
     point[0].color=java.awt.Color.blue;
     point[0].pixRadius=2;
     point[1].setXY(f2, 0);             //place the 2nd focal point
     point[1].color=java.awt.Color.blue;
     point[1].pixRadius=2;
     polarPanel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     polarPanel.setSquareAspect(true);
     trail.color=java.awt.Color.red;
     polarPanel.addDrawable(trail);
     polarPanel.addDrawable(point[0]);  //1st focal point
     polarPanel.addDrawable(point[1]);  //2nd focal point
     axes.setDeltaR(Math.sqrt(xmax*xmax+ymax*ymax)/5);
     axes.setDeltaTheta(Math.PI/6);
     polarFrame.render();
     myControl.println("Ellipse Properties");
     myControl.println("Maximum distance: rmax="+nf.format(rmax));
     myControl.println("Semimajor axis: a="+nf.format(a));
     myControl.println("Blue Focal points at: "+f1+", and "+f2);
     }

   public void clear () {
     trail.clear();
     polarPanel.clear();
     polarFrame.render();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Draws an ellipse of minimum radius rmin and eccentricity");
     myControl.println ("e according to r=rmin*(1+e)/(1+e*cos(theta)) with ");
     myControl.println ("rmax=rmin*(1+e)/(1-e), a=(rmin+rmax)/2, and distance");
     myControl.println ("between focal points=2*a*e. Here a=semimajor axis, and");
     myControl.println ("a maximum radius.");
     r    = null;
     th   = null;
     e=0.6;
     rmin=2.5;
     th0=0.0;
     thmax=2*pi;
     NPTS=200;
     myControl.setValue("eccentricity",e);
     myControl.setValue("rmin",rmin);
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
     Calculation model = new ellipseApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(420, 5);
     myControl.setSize(350,590);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}

