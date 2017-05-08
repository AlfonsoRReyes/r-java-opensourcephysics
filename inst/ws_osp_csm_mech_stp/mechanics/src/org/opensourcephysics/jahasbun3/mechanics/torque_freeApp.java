/*
@Author J E Hasbun 2007.
Plots the frequency and angular momentum for torque free motion of a top versus
time in the body (S') frame.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.controls.*;
import java.text.*;

public class torque_freeApp extends AbstractAnimation {
  Display3DFrame panel = new Display3DFrame("Torque Free Motion of a Top");
  private Control myControl;
  double w1[], w2[], t[], L1[], L2[], L, w3, L3;
  double tmax, dt, pi, cf, fi_b, fi_L, fi_s, c, omb, fi;
  double I1, I2, I3, gam, omL, w;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int NPTS, istep, Iarrows=3;
  ElementTrail [] trace=new ElementTrail[2];
  ElementArrow [] arrow=new ElementArrow [Iarrows];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public torque_freeApp(){
    panel.setDecorationType(VisualizationHints.DECORATION_AXES);
    panel.setLocation(5,5);
    panel.setSize(400,550);
    pi=Math.PI; cf=180./pi;
    nf.setMaximumFractionDigits(4);
    trace[0]=new ElementTrail();
    trace[1]=new ElementTrail();
    trace[0].setConnected(false);
    trace[0].getStyle().setLineColor(java.awt.Color.black);
    trace[1].setConnected(false);
    trace[1].getStyle().setLineColor(java.awt.Color.red);
    for(int i=0; i< Iarrows; i++){
      arrow[i]=new ElementArrow();
    }
  }

   protected void doStep() {
     if (istep==0){
       myControl.clearMessages();
       myControl.println("Values of important variables:");
       myControl.println("magnitudes: w=" + nf.format(w) + " rad/s, L=" +
       nf.format(L) +" kgm^2/s\nangles (degrees): fi_b=" +nf.format(fi_b * cf) +
       ", fi_L=" +nf.format(fi_L * cf) + ", fi_s=" + nf.format(fi_s * cf) +
       "\nfrequencies (rad/s) omb=" + nf.format(omb) +", omL=" + nf.format(omL));
     }
     if (istep < NPTS){
       calculate();
       trace[0].addPoint(w1[istep], w2[istep], w3);      //w cone (body)
       trace[1].addPoint(L1[istep], L2[istep], L3);      //L cone (space)
       //arrows start at origin set by cleararrows()
       arrow[0].setSizeXYZ(w1[istep],w2[istep],w3);      //w arrow
       arrow[1].setSizeXYZ(L1[istep],L2[istep],L3);      //L arrow
       panel.addElement(trace[0]);
       panel.addElement(trace[1]);
       panel.addElement(arrow[0]);
       panel.addElement(arrow[1]);
       panel.render();
       istep++;
     } else {
         stopAnimation();
         myControl.println("black trace and arrow - w vector");
         myControl.println("red trace and arrow - L vector");
         myControl.println("blue arrow - w3");
         myControl.println("Press 'Stop' before beginning a new animation");
     }
   }

   public void calculate() {
       t[istep]=istep*dt;
       w1[istep]=c*Math.cos(omb*t[istep]+fi/cf); //w components vs time
       w2[istep]=c*Math.sin(omb*t[istep]+fi/cf);
       L1[istep]=I1*w1[istep];                   //Ang. Mom. components
       L2[istep]=I2*w2[istep];
     }

   public void clear  () {
     trace[0].clear();
     trace[1].clear();
     panel.addElement(trace[0]);
     panel.addElement(trace[1]);
     for(int i=0; i< Iarrows; i++){
       arrow[i].setXYZ(0, 0, 0); //x,y,z arrow position
       //arrow length corresponds to the vector function x,y,z components
       arrow[i].setSizeXYZ(0, 0, 0);
       //clear the arrows from the frame
       panel.addElement(arrow[i]);
     }
     panel.setPreferredMinMax(0,0,0,0,0,0);
     panel.render();
   }

   public void resetAnimation() {
     istep=0;                     //doStep index counter
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the frequency and angular momentum for torque");
     myControl.println ("free motion of a top versus time in the body (S') frame.");
     myControl.println ("w=(w1,w2,w3), where w1=c*cos(omb*t+fi), w2=c*sin(omb*t+fi),");
     myControl.println ("L=(L1,L2,L3), where L1=I1*w1, L2=I2*w2, L3=I3*w3,");
     myControl.println ("gam=(I3-I1)/I1, omb=gam*w3, angle between w and w3:");
     myControl.println ("fi_b=atan(c/w3); angle between w and L: fi_L - see text,");
     myControl.println ("angle between L and w3: fi_s=atan(I1*tan(fi_b)/I3);");
     myControl.println ("prec. freq. of w about L: omL=omb*sin(fi_b)/sin(fi_L);");
     myControl.println ("magnitude of w: sqrt(c*c+w3*w3); ");
     myControl.println ("magnitude of L: sqrt(I1^2*c^2+I3^2*w3^3).");
     myControl.println("Press 'Initialize' before beginning the animation");
     w1    = null;
     w2    = null;
     L1    = null;
     L2    = null;
     t     = null;
     I1=1.3; I2=I1; I3=1.5;  //principal moments of inertia
     w3=1.0; fi=0.0; c=1.0;  //initial values
     gam=(I3-I1)/I1;         //gamma
     omb=gam*w3;
     tmax=2.*pi/omb;
     NPTS=200;
     delayTime=40;        //time in between animation steps (Abstract Animation)
     myControl.setValue("I1",I1);
     myControl.setValue("I2",I2);
     myControl.setValue("I3",I3);
     myControl.setValue("w3",w3);
     myControl.setValue("c",c);
     myControl.setValue("fi_angle(degrees)",fi);
     myControl.setValue("tmax(s)",nf.format(tmax));
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",nf.format(dt));
   }

   public void initializeAnimation() {
     clear();
     istep=0;             //increased in doStep()
     myControl.clearMessages();
     I1=myControl.getDouble("I1");
     I2=myControl.getDouble("I2");
     I3=myControl.getDouble("I3");
     w3=myControl.getDouble("w3");
     c=myControl.getDouble("c");
     fi=myControl.getDouble("fi_angle(degrees)");
     tmax=myControl.getDouble("tmax(s)");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",nf.format(dt));
     w1    = new double[NPTS];
     w2    = new double[NPTS];
     L1    = new double[NPTS];
     L2    = new double[NPTS];
     t     = new double[NPTS];
     //initial conditions
     gam=(I3-I1)/I1;                        //gamma
     omb=gam*w3;                            //precessional frequency
     L3=I3*w3;                              //z-Ang. Mom. component
     fi_b=Math.atan(c/w3);                  //angle between w and w3
     //angle between w and L
     fi_L=Math.acos((I1*c*c+I3*w3*w3)/(Math.sqrt((c*c+w3*w3)*(I1*I2*c*c+
          I3*I3*w3*w3))));
     fi_s=Math.atan(I1*Math.tan(fi_b)/I3);     //angle between L and w3
     omL=omb*Math.sin(fi_b)/Math.sin(fi_L);    //prec. freq. of w about L
     //omL=Math.sqrt(I1*I1*c*c+I3*I3*w3*w3)/I1 //another formula for omL
     w=Math.sqrt(c*c+w3*w3);                   //magnitude of w
     L=Math.sqrt(I1*I1*c*c+I3*I3*w3*w3);       //magnitude of L
     xmin=-L;
     xmax=-xmin;
     ymin=xmin;
     ymax=-ymin;
     zmin=0;
     zmax=L3;
     panel.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
     panel.setAltitude(0.25);
     panel.setAzimuth(0.5);
     arrow[0].getStyle().setLineColor(java.awt.Color.black);//w components
     arrow[1].getStyle().setLineColor(java.awt.Color.red);  //Ang. Mom.
     arrow[2].getStyle().setLineColor(java.awt.Color.blue); //w3
     arrow[2].setSizeXYZ(0.0,0.0,w3);                       //w3 arrow
     panel.addElement(arrow[2]);
     panel.render();
     myControl.println("Press 'start' to begin animation");
     myControl.println("Feel free to mouse-rotate the 3D display frame for better view.");
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation model = new torque_freeApp();
     AnimationControl myControl;
     myControl = new AnimationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(375,550);
     myControl.setDividerLocation(235);
     model.setControl(myControl);
   }
}

