/*
@Author J E Hasbun 2007.
Binary star system given the eccentricity. We use the center of mass - relative
coordinate method.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;
import java.text.NumberFormat;

  public class binary1App extends AbstractAnimation {
  PlottingPanel panel = new PlottingPanel
                          ("x(AU)","y(AU)","Binary System");
  DrawingFrame  frame = new DrawingFrame(panel);
  Trail [] trail=new Trail[3];
  Circle [] point=new Circle[4];
  private Control myControl;
  double [] r, r1, r2,th;
  double [] x,y;
  double m,m1,m2,e,rcm,rmin,rmax,a,tau,dth,th0,thmax,af;
  double xmin,xmax,ymin,ymax;
  double pi;
  int NPTS,istep;
  NumberFormat nf = NumberFormat.getInstance();

  public binary1App (){
    frame.setLocation(5,10);
    frame.setSize(400,400);
    frame.setTitle("Binary System y vs x");
    nf.setMaximumFractionDigits(3);
    pi=Math.PI;
    trail[0]=new Trail();
    trail[1]=new Trail();
    trail[2]=new Trail();
  }

   protected void doStep() {
     if(istep+1<NPTS){
       istep++;
       calculate();
       //reduced mass
       x[0] = r[istep] * Math.cos(th[istep]);
       y[0] = r[istep] * Math.sin(th[istep]);
       trail[0].addPoint(x[0], y[0]);
       point[0].setXY(x[0],y[0]);
       //mass 1
       x[1] = r1[istep] * Math.cos(th[istep]);
       y[1] = r1[istep] * Math.sin(th[istep]);
       trail[1].addPoint(x[1],y[1]);
       point[1].setXY(x[1],y[1]);
       //mass 2
       x[2] = r2[istep] * Math.cos(th[istep]);
       y[2] = r2[istep] * Math.sin(th[istep]);
       trail[2].addPoint(x[2],y[2]);
       point[2].setXY(x[2],y[2]);
       //draw
       panel.addDrawable(trail[0]);
       panel.addDrawable(point[0]);
       panel.addDrawable(trail[1]);
       panel.addDrawable(point[1]);
       panel.addDrawable(trail[2]);
       panel.addDrawable(point[2]);
       checkPanel();
     }else{
       //place the center of mass point also
       point[3].setXY(rcm, 0);
       panel.addDrawable(point[3]);
       frame.render();
       stopAnimation();
       myControl.println("black:reduced mass, blue:m1, red:m2, grey:cm");
       myControl.println("(rmin,rmax)=\t(" + nf.format(rmin) + "," +
                         nf.format(rmax) + ") AU");
       myControl.println("tau=\t" + nf.format(tau) + " tau0, a=\t" + nf.format(a) +
                         " AU");
       myControl.println("(m,m1,m2)=\t(" + m + "," + m1 + "," + m2 + ")Ms");
       myControl.println("reduced mass=\t" + nf.format(m1 * m2 / m) + " Ms");
     }
   }

   public void calculate () {
     //orbit calculation
     th[istep]=th0+istep*dth;             //angle variable
     r[istep]=rmin*(1+e)/(1+e*Math.cos(th[istep]));
     r1[istep]=rcm-m2*r[istep]/m;
     r2[istep]=rcm+m1*r[istep]/m;
   }

   public void checkPanel() {
     if(x[0]<xmin||x[1]<xmin||x[2]<xmin){xmin=ArrayLib.min(x);}
     if(x[0]>xmax||x[1]>xmax||x[2]>xmax){xmax=ArrayLib.max(x);}
     if(y[0]<ymin||y[1]<ymin||y[2]<ymin){ymin=ArrayLib.min(y);}
     if(y[0]>ymax||y[1]>ymax||y[2]>ymax){ymax=ArrayLib.max(y);}
     panel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     frame.render();
    }

   public void clear () {
       trail[0].clear();
       trail[1].clear();
       trail[2].clear();
       panel.clear();
       frame.render();
   }

   public void resetAnimation () {
     clear();
     myControl.clearMessages();
     myControl.println ("Binary star system given the eccentricity. We use the");
     myControl.println (" center of mass-relative coordinate method. The relative");
     myControl.println ("coordinate is: r=rmin*(1+e)/(1+e*Math.cos(theta)) for");
     myControl.println ("the reduced mass. m1's coordinate is r1=rcm-m2*r/m, and");
     myControl.println ("m2's coordinate is r2=rcm+m1*r/m. Reduced mass=m1*m2/m.");
     r    = null;
     r    = null;
     r    = null;
     th    = null;
     m=5; m1=3; m2=m-m1;   //initial total and individual masses in Ms
     rcm=0; rmin=2.5;      //center of mass, rmin for relative coordinate in AU
     e=0.6;                //eccentricity
     th0=0; thmax=2*pi;
     NPTS=150;
     delayTime=25;         //time in between animation steps (Abstract Animation)
     myControl.setValue("m",m);
     myControl.setValue("m1",m1);
     myControl.setValue("eccentricity",e);
     myControl.setValue("rmin",rmin);
     myControl.setValue("rcm",rcm);
     myControl.setValue("th0",th0);
     myControl.setValue("thmax",thmax);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     dth=(thmax-th0)/(NPTS-1);
     myControl.setValue("dth((thmax-th0)/(NPTS-1))=",dth);
     m2=m-m1;
     myControl.setValue("m2(m-m1)=",m2);
   }

   public void initializeAnimation() {
     //initial conditions
     istep=-1;  //so that istep+1 in doStep =0 the first time
     clear();
     myControl.clearMessages();
     m=myControl.getDouble("m");
     m1=myControl.getDouble("m1");
     e=myControl.getDouble("eccentricity");
     rmin=myControl.getDouble("rmin");
     rcm=myControl.getDouble("rcm");
     th0=myControl.getDouble("th0");
     thmax=myControl.getDouble("thmax");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dth=(thmax-th0)/(NPTS-1);
     myControl.setValue("dth((thmax-th0)/(NPTS-1))=",dth);
     m2=m-m1;
     myControl.setValue("m2(m-m1)=",m2);
     r     = new double[NPTS];
     r1    = new double[NPTS];
     r2    = new double[NPTS];
     th    = new double[NPTS];
     x = new double[3];
     y = new double[3];
     rmax=rmin*(1+e)/(1-e);
     a=(rmin+rmax)/2;
     tau=Math.sqrt(Math.pow(a,3.0)/m); //reduced mass, m1, & m2 have the same period
     af=(1+0.1);
     xmin=-rmax/2; xmax=rmax/2;
     ymin=xmin; ymax=xmax;
     panel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     panel.setSquareAspect(true);
     trail[0].color=java.awt.Color.black;
     trail[1].color=java.awt.Color.blue;
     trail[2].color=java.awt.Color.red;
     point[0]=new Circle();
     point[1]=new Circle();
     point[2]=new Circle();
     point[3]=new Circle();
     point[0].color=java.awt.Color.black;
     point[0].pixRadius=3;
     point[1].color=java.awt.Color.blue;
     point[1].pixRadius=3;
     point[2].color=java.awt.Color.red;
     point[2].pixRadius=3;
     point[3].color=java.awt.Color.gray;
     point[3].pixRadius=3;
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new binary1App();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(410, 5);
     myControl.setSize(360,450);
     myControl.setDividerLocation(250);
     app.setControl(myControl);
   }
}