/*
@Author J E Hasbun 2007.
Draws Earth's elliptical orbit around the sun, given its mass, its eccentricity,
its semimajor axis, the sun's mass, and G.
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
import java.text.NumberFormat;

public class earthorbApp implements Calculation {
  PlottingPanel polarPanel = new PlottingPanel
                          ("Theta","r","Polar Plot - Earth Orbit");
  DrawingFrame  polarFrame = new DrawingFrame(polarPanel);
  PolarType1 axes = new PolarType1(polarPanel);
  Trail trail=new Trail();
  Circle [] point=new Circle[2];
  private Control myControl;
  double r[],th[];
  double m,G,M,a,e,K,vc,v0,r0,rmax,v_rmax,th0,thmax,dth;
  double pi;
  int NPTS;
  NumberFormat nf = NumberFormat.getInstance();

  public earthorbApp (){
    polarFrame.setLocation(5,5);
    polarFrame.setSize(400,590);
    polarFrame.setTitle("Earth Orbit Polar Plot - r(theta)");
    nf.setMaximumFractionDigits(3);
    pi=Math.PI;
    point[0]=new Circle();
    point[1]=new Circle();
  }

   public void calculate () {
     clear();
     double xmin,xmax,ymin,ymax, af=(1+0.1);
     myControl.clearMessages();
     m=myControl.getDouble("Earth_m");
     M=myControl.getDouble("Sun_M");
     a=myControl.getDouble("semimajor_axis");
     e=myControl.getDouble("eccentricity");
     G=myControl.getDouble("G");
     th0=myControl.getDouble("th0");
     thmax=myControl.getDouble("thmax");
     NPTS=myControl.getInt("NPTS");
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth((thmax-th0)/NPTS)=",dth);
     r     = new double[NPTS];
     th    = new double[NPTS];
     double [] xy_coord = new double[2];
     //initial conditions
     r0=a*(1-e);                    //same as rmin
     K=-G*M*m/Math.pow(1000,3);     //in km^3.kg/s^2
     vc=Math.sqrt(-K/r0/m);         //speed in km/s (circular orbit)
     v0=vc*Math.sqrt(1+e);          //speed at r0
     rmax=r0*(v0/vc)*(v0/vc)/(2-(v0/vc)*(v0/vc)); //same as r0*(1+e)/(1-e)
     v_rmax=r0*v0/rmax;             //using angular momentum conservation
     xmin=0; xmax=0; ymin=0; ymax=0;
     //orbit calculation
     for (int i=0; i< NPTS;i++){
       th[i]=th0+i*dth;             //angle variable
       r[i]=r0*(1+e)/(1+e*Math.cos(th[i]));
       xy_coord[0]=r[i]*Math.cos(th[i]);
       xy_coord[1]=r[i]*Math.sin(th[i]);
       trail.addPoint(xy_coord[0],xy_coord[1]);
       if(xy_coord[0]<xmin){xmin=xy_coord[0];}
       if(xy_coord[0]>xmax){xmax=xy_coord[0];}
       if(xy_coord[1]<ymin){ymin=xy_coord[1];}
       if(xy_coord[1]<ymax){ymin=xy_coord[1];}
     }
     point[0].setXY(xmin,0);            //place a point at rmin
     point[0].color=java.awt.Color.blue;
     point[0].pixRadius=2;
     point[1].setXY(xmax,0);            //place apoint at rmax
     point[1].color=java.awt.Color.blue;
     point[1].pixRadius=2;
     polarPanel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     polarPanel.setSquareAspect(true);
     trail.color=java.awt.Color.red;
     polarPanel.addDrawable(trail);
     polarPanel.addDrawable(point[0]);
     polarPanel.addDrawable(point[1]);
     axes.setDeltaR(Math.sqrt(xmax*xmax+ymax*ymax)/5);
     axes.setDeltaTheta(Math.PI/6);
     polarFrame.render();
     myControl.println("Red line: Earth's orbit");
     myControl.println("Blue dots: (rmin;rmax)=("+nf.format(r0)+"; "+nf.format(rmax)+") km");
     myControl.println("(v0;v_rmax)=("+nf.format(v0)+"; "+nf.format(v_rmax)+") km/s");
     }

     public void clear () {
       trail.clear();
       polarPanel.clear();
       polarFrame.render();
     }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Draws Earth's elliptical orbit around the sun, given its");
     myControl.println ("mass, its eccentricity, its semimajor axis, the sun's mass,");
     myControl.println ("and G. Other planets can similarly be done.");
     myControl.println ("The formulas used are: for the orbit r=r0*(1+e)/(1+e*cos(th))");
     myControl.println ("with r0=a*(1-e),  K=-G*M*m, vc=sqrt(-K/r0/m), v0=vc*sqrt(1+e),");
     myControl.println ("rmax=r0*(v0/vc)*(v0/vc)/(2-(v0/vc)*(v0/vc)), v_rmax=r0*v0/rmax.");
     r    = null;
     th    = null;
     m=5.98e24;       //Earth's mass in kg
     G=6.67e-11;      //Universal gravitational constant
     M=1.991e30;      //Sun's mass in kg
     a=1.496e8;       //average Earth-Sun distance in km, or semimajor axis
     e=0.017;         //Earth's eccentricity
     th0=0; thmax=2*pi;
     NPTS=500;
     myControl.setValue("Earth_m",m);
     myControl.setValue("Sun_M",M);
     myControl.setValue("semimajor_axis",a);
     myControl.setValue("eccentricity",e);
     myControl.setValue("G",G);
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
     Calculation model = new earthorbApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(350,590);
     myControl.setDividerLocation(250);
     model.setControl(myControl);
   }
}