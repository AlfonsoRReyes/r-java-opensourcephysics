/*
@Author J E Hasbun 2007.
Animates the position of a rod-mass system, its angular momentum, and its
torque for rotation about a fixed axis.
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
import org.opensourcephysics.numerics.*;
import java.text.*;

public class fixed_axisApp extends AbstractAnimation {
  Display3DFrame panel = new Display3DFrame
       ("Fixed-axis massless rod-point mass system rotation");
  private Control myControl;
  double x[], y[], t[], Lx[], Ly[], Lz, Tx[], Ty[], Tz;
  double tmax, t0, dt, w, tau, pi, th, a, thr, s, c;
  double R, I, m, A, ct, st, z, co;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int NPTS, istep, Iarrows=4;
  ElementTrail trace    = new ElementTrail();
  ElementCircle point   = new ElementCircle();
  ElementArrow [] arrow = new ElementArrow[Iarrows];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public fixed_axisApp(){
    panel.setDecorationType(VisualizationHints.DECORATION_AXES);
    panel.setLocation(5,5);
    panel.setSize(400,550);
    trace.setConnected(false);
    trace.getStyle().setLineColor(java.awt.Color.blue);
    for(int i=0; i< Iarrows; i++){
      arrow[i]=new ElementArrow();
    }
    pi=Math.PI;
    nf.setMaximumFractionDigits(3);
  }

   protected void doStep() {
     if (istep < NPTS){
       calculate();
       trace.addPoint(x[istep], y[istep], z);      //particle position
       //arrows start at origin from the start by clear()
       arrow[0].setSizeXYZ(x[istep],y[istep],z);   //particle position vector
       arrow[1].setSizeXYZ(Lx[istep],Ly[istep],Lz);//angular momentum components
       arrow[2].setSizeXYZ(Tx[istep],Ty[istep],Tz);//torque components
       arrow[3].setSizeXYZ(0.,0.,z);               //z-axis
       point.setXYZ(x[istep],y[istep],z);          //particle
       panel.addElement(point);
       panel.addElement(trace);
       panel.addElement(arrow[0]);
       panel.addElement(arrow[1]);
       panel.addElement(arrow[2]);
       panel.addElement(arrow[3]);
       panel.render();
       istep++;
     } else {
         stopAnimation();
         myControl.clearMessages();
         myControl.println("black line: particle (black dot) position with blue trace");
         myControl.println("red line:   Angular momentum ");
         myControl.println("green line: Torque");
         myControl.println("cyan line:  z axis");
         myControl.println("Some values calculated: \nR(m)="+nf.format(R)+", I(kgm^2)="+
                       nf.format(I)+", A(kgm^2/s)="+nf.format(A)+"\nz(m)="+nf.format(z)+
                       ", theta="+nf.format(th)+" degrees or "+nf.format(thr)+
                       " radians");
         myControl.println("Press 'Stop' before beginning a new animation");
     }
   }

   public void calculate() {
       t[istep]=t0+istep*dt;
       ct=Math.cos(w*t[istep]);
       st=Math.sin(w*t[istep]);
       x[istep]=R*ct;            //time dependent positions
       y[istep]=R*st;
       Lx[istep]=-A*co*ct;       //time dependent angular momentum components
       Ly[istep]=-A*co*st;
       Tx[istep]=A*w*co*st;      //time dependent torque components
       Ty[istep]=-A*w*co*ct;
     }

   public void clear () {
     trace.clear();
     point.setXYZ(0,0,0);
     point.setSizeXYZ(0, 0, 0);
     point.getStyle().setFillColor(java.awt.Color.white);
     panel.addElement(trace);
     panel.addElement(point);
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
     myControl.println ("Animates the position of a rod-mass system, its angular");
     myControl.println ("momentum, and its torque for rotation about a fixed axis.");
     myControl.println ("Here, the particle position is:");
     myControl.println ("x=R*cos(w*t), y=R*cos(w*t), z=a*cos(theta).");
     myControl.println ("The angular momentum is Lx=A*cot(theta)*cos(w*t),");
     myControl.println ("Ly=A*cot(theta)*sin(w*t), Lz=A=I*w, I=m*R^2, and the");
     myControl.println ("torque is Tx=dLx/dt, Ty=dLy/dt, Tz=dLz/dt. Also thta is");
     myControl.println ("the rod tilt angle and R=a*sin(theta), with a=rod length.");
     myControl.println("Press 'Initialize' before beginning the animation");
     x    = null;
     y    = null;
     t    = null;
     m=1.0;               //mass
     a=1.0;               //rod length
     //w=pi;              //angular frequency previously used
     w=1.0; tau=2*pi/w;   //angular velocity and period
     th=25;               //angle of tilt
     t0=0.0; tmax=tau;    //time range
     NPTS=100;
     delayTime=50;        //time in between animation steps (Abstract Animation)
     myControl.setValue("m-mass(kg)",m);
     myControl.setValue("a-rod_length(m)",a);
     myControl.setValue("w-frequency(rad/s)",nf.format(w));
     myControl.setValue("theta-rod_tilt_angle(degrees)",th);
     myControl.setValue("t0(s)",t0);
     myControl.setValue("tmax(s)",nf.format(tmax));
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",nf.format(dt));
   }

   public void initializeAnimation() {
     clear();
     istep=0;             //increased in doStep()
     myControl.clearMessages();
     m=myControl.getDouble("m-mass(kg)");
     a=myControl.getDouble("a-rod_length(m)");
     w=myControl.getDouble("w-frequency(rad/s)");
     th=myControl.getDouble("theta-rod_tilt_angle(degrees)");
     t0=myControl.getDouble("t0(s)");
     tmax=myControl.getDouble("tmax(s)");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",nf.format(dt));
     x     = new double[NPTS];
     y     = new double[NPTS];
     Lx    = new double[NPTS];
     Ly    = new double[NPTS];
     Tx    = new double[NPTS];
     Ty    = new double[NPTS];
     t     = new double[NPTS];
     double [][] v=new double [2][3];
     //initial conditions
     thr=Math.toRadians(th);
     s=Math.sin(thr);
     c=Math.cos(thr);
     R=a*s;
     z=a*c;               //z-position
     co=1./Math.tan(thr);
     I=m*R*R;
     A=I*w;
     Lz=A;                //angular momentum z-component
     Tz=0.0;              //torque z-component
     v[0][0]=R; v[0][1]=z; v[0][2]=A*w*co;
     v[1][0]=R; v[1][1]=z; v[1][2]=A;
     xmin=-ArrayLib.max(v[0]);
     xmax=-xmin;
     ymin=xmin;
     ymax=-ymin;
     zmin=0;
     zmax=ArrayLib.max(v[1]);
     panel.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
     panel.setAltitude(0.8);
     panel.setAzimuth(0.1);
     arrow[0].getStyle().setLineColor(java.awt.Color.black);//position x,y,z
     arrow[1].getStyle().setLineColor(java.awt.Color.red);  //Ang. Mom.
     arrow[2].getStyle().setLineColor(java.awt.Color.green);//Torque
     arrow[3].getStyle().setLineColor(java.awt.Color.cyan); //z axiz
     point.getStyle().setFillColor(java.awt.Color.black);//the particle
     double scale_factor=Math.sqrt(xmax*xmax+ymax*ymax+zmax*zmax);
     double p_size=0.015*scale_factor;
     point.setSizeXYZ(p_size,p_size,p_size);//scale the point according to frame
     myControl.println("Press 'start' to begin animation");
     myControl.println("Feel free to mouse-rotate the 3D display frame for better view.");
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation model = new fixed_axisApp();
     AnimationControl myControl;
     myControl = new AnimationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(375,550);
     myControl.setDividerLocation(300);
     model.setControl(myControl);
   }
}

