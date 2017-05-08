/*
@Author J E Hasbun 2007.
Plots the frequency and angular momentum for torque free motion of a top versus
time in the body (S') frame as well as in the space frame (S).
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

public class torque_free_sApp extends AbstractAnimation {
  Display3DFrame panel = new Display3DFrame("Torque Free Motion of a Top");
  private Control myControl;
  double w1[], w2[], t[], L1, L2, L3, w3, L;
  double tmax, dt, pi, cf, fi_b, fi_L, fi_s, c, omb, fi;
  double I1, I2, I3, gam, omL, w;
  double e31[], e32[], e33;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int NPTS, istep, Iarrows=3;
  ElementTrail trace=new ElementTrail();
  ElementArrow [] arrow=new ElementArrow [Iarrows];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public torque_free_sApp(){
    panel.setDecorationType(VisualizationHints.DECORATION_AXES);
    panel.setLocation(5,5);
    panel.setSize(400,550);
    pi=Math.PI; cf=180./pi;
    nf.setMaximumFractionDigits(4);
    trace=new ElementTrail();
    trace.setConnected(false);
    trace.getStyle().setLineColor(java.awt.Color.black);
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
       trace.addPoint(w1[istep], w2[istep], w3);         //S frame w
       //arrows start at origin set by cleararrows()
       arrow[0].setSizeXYZ(w1[istep],w2[istep],w3);      //S frame w arrow
       arrow[1].setSizeXYZ(e31[istep],e32[istep],e33);                    //e3 arrow
       panel.addElement(trace);
       panel.addElement(arrow[0]);
       panel.render();
       istep++;
     } else {
         stopAnimation();
         myControl.println("black trace and arrow - w vector S-frame");
         myControl.println("red arrow  - L vector ");
         myControl.println("blue arrow - w3 ");
         myControl.println("Press 'Stop' before beginning a new animation");
     }
   }

   public void calculate() {
       t[istep]=istep*dt;
       //before L1, L2 were arrays and w1, w2 did not involve fi_s
       //w1[istep]=c*Math.cos(omb*t[istep]+fi/cf); //w components vs time
       //w2[istep]=c*Math.sin(omb*t[istep]+fi/cf);
       //L1[istep]=I1*w1[istep];                //Ang. Mom. components
       //L2[istep]=I2*w2[istep];
       //=== shift old e3 by fi_s to swap places with old L, preserve magnitude
       e31[istep]=1.*Math.sin(fi_s)*Math.cos(omL*t[istep]+fi);//must move with omL
       e32[istep]=1.*Math.sin(fi_s)*Math.sin(omL*t[istep]+fi);//must move with omL
       //new w and components, shift by fi_s, preserve magnitude, change rate
       w1[istep]=c*Math.cos(omL*t[istep]+fi);   //use omL precession rate
       w2[istep]=c*Math.sin(omL*t[istep]+fi);
     }

   public void clear  () {
     trace.clear();
     panel.addElement(trace);
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
     myControl.println ("free motion of a top versus time in the body (S') frame");
     myControl.println ("as well as in the space frame (S). See text for details ");
     myControl.println ("regarding the shift of angle made and rate change made");
     myControl.println ("compared to torque_freeApp.");
     myControl.println("Press 'Initialize' before beginning the animation");
     w1    = null;
     w2    = null;
     e31   = null;
     e32   = null;
     t     = null;
     I1=1.3; I2=I1; I3=1.8;  //principal moments of inertia
     w3=1.0; fi=0.0; c=1.0;  //initial values
     gam=(I3-I1)/I1;         //gamma
     omb=gam*w3;
     tmax=0.9*2.*pi/omb;
     NPTS=100;
     delayTime=40;           //time in between animation steps (Abstract Animation)
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
     e31   = new double[NPTS];
     e32   = new double[NPTS];
     t     = new double[NPTS];
     //initial conditions
     gam=(I3-I1)/I1;                        //gamma
     omb=gam*w3;                            //precessional frequency
     //angle between w and L
     fi_b=Math.atan(c/w3);                     //angle between w and w3
     fi_L=Math.acos((I1*c*c+I3*w3*w3)/(Math.sqrt((c*c+w3*w3)*(I1*I2*c*c+
          I3*I3*w3*w3))));
     fi_s=Math.atan(I1*Math.tan(fi_b)/I3);     //angle between L and w3
     omL=omb*Math.sin(fi_b)/Math.sin(fi_L);    //prec. freq. of w about L
     w=Math.sqrt(c*c+w3*w3);                   //magnitude of w
     L=Math.sqrt(I1*I1*c*c+I3*I3*w3*w3);       //magnitude of L
     //shift old e3 by fi_s to swap places with old L, preserve magnitude
     e33=1.*Math.cos(fi_s);
     //new w and components, shift by fi_s, preserve magnitude, change rate
     c=w*Math.sin(fi_b-fi_s);
     w3=w*Math.cos(fi_b-fi_s);
     //new L now in old w3 direction and fixed, preserve magnitude
     L1=0; L2=0; L3=L;
     L3=I3*w3;                                 //z-Ang. Mom. component
     xmin=-Math.sin(fi_s);                     //amplitude of e31
     xmax=-xmin;
     ymin=xmin;
     ymax=-ymin;
     zmin=0;
     zmax=L3;
     panel.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
     panel.setAltitude(0.25);
     panel.setAzimuth(0.5);
     arrow[0].getStyle().setLineColor(java.awt.Color.black); //S frame w
     arrow[1].getStyle().setLineColor(java.awt.Color.blue);  //e3 arrow
     arrow[2].getStyle().setLineColor(java.awt.Color.red);   //L arrow
     arrow[2].setSizeXYZ(L1,L2,L3);                          //L arrow
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
     Animation model = new torque_free_sApp();
     AnimationControl myControl;
     myControl = new AnimationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(375,550);
     myControl.setDividerLocation(235);
     model.setControl(myControl);
   }
}