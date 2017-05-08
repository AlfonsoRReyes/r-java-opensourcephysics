/*
@Author J E Hasbun 2007.
Program to solve Euler's equations, plot them and their respective angular
speeds, as well as angular momentum properties.
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
import org.opensourcephysics.frames.*;
import java.text.*;

public class topApp extends AbstractAnimation {
  Display3DFrame frame3d = new Display3DFrame("Ellipsoid System w, L animation");
  PlotFrame plot0=new PlotFrame ("t","phi,theta,psi"," Angles versus time");
  PlotFrame plot1=new PlotFrame ("t","phi_{d},theta_{d},psi_{d}",
                                 " Angle rates versus time");
  PlotFrame plot2=new PlotFrame ("theta","V,E"," Veff and E' vs theta");
  private Control myControl;
  double w1[], w2[], w3[], w4[], w5[], w6[], t[], Lx[], Ly[], Lz[];
  double Lzppx[], Lzppy[], Lzppz[], Vef[], Ekp[], Ep[];
  double tmax, dt, pi, ph0, th0, ps0, ph0d, th0d;
  double I, Is, g, M, a, tau0, ws;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int NPTS, istep, Iarrows=4;
  double[] state= new double [6];
  ODE ode = new top_der();
  ODESolver odesolver = new RK4(ode);
  //ODESolver odesolver = new RK45MultiStep(ode);
  ElementTrail [] trace=new ElementTrail[2];
  ElementArrow [] arrow=new ElementArrow [Iarrows];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public topApp(){
    //the animation frame
    frame3d.setDecorationType(VisualizationHints.DECORATION_CUBE);
    frame3d.setLocation(5,5);
    frame3d.setSize(260,295);
    //The angles frame
    plot0.setConnected(true);//set all connected
    plot0.setLocation(270,5);
    plot0.setSize(260,295);
    plot0.setLineColor(0,java.awt.Color.black);//for phi
    plot0.setMarkerSize(0,0);
    plot0.setLineColor(1,java.awt.Color.blue);  //for theta
    plot0.setMarkerSize(1,0);
    plot0.setLineColor(2,java.awt.Color.red); //for psi
    plot0.setMarkerSize(2,0);
    //The rates frame
    plot1.setConnected(true);//set all connected
    plot1.setLocation(5,305);
    plot1.setSize(260,295);
    plot1.setLineColor(0,java.awt.Color.black);//for phi_d
    plot1.setMarkerSize(0,0);
    plot1.setLineColor(1,java.awt.Color.blue);  //for theta_d
    plot1.setMarkerSize(1,0);
    plot1.setLineColor(2,java.awt.Color.red); //for psi_d
    plot1.setMarkerSize(2,0);
    //The energies frame
    plot2.setConnected(true);//set all connected
    plot2.setLocation(270,305);
    plot2.setSize(260,295);
    plot2.setLineColor(0,java.awt.Color.green);//for Veff
    plot2.setMarkerSize(0,0);
    plot2.setLineColor(1,java.awt.Color.blue); //for Ep
    plot2.setMarkerSize(1,0);
    pi=Math.PI;
    nf.setMaximumFractionDigits(4);
    trace[0]=new ElementTrail();
    trace[1]=new ElementTrail();
    trace[0].setConnected(false);                           //L
    trace[0].getStyle().setLineColor(java.awt.Color.black);
    trace[1].setConnected(false);                           //Lzpp
    trace[1].getStyle().setLineColor(java.awt.Color.blue);
    for(int i=0; i< Iarrows; i++){
      arrow[i]=new ElementArrow();
    }
    arrow[0].getStyle().setLineColor(java.awt.Color.red);     //Lz
    arrow[1].getStyle().setLineColor(java.awt.Color.black);   //L
    arrow[2].getStyle().setLineColor(java.awt.Color.magenta); //Lzppz
    arrow[3].getStyle().setLineColor(java.awt.Color.blue);    //Lzpp
  }

   protected void doStep() {
     if (istep < NPTS){
       arrow[0].setSizeXYZ(0.0,0.0,Lz[istep]);                      //Lz
       arrow[1].setSizeXYZ(Lx[istep],Ly[istep],Lz[istep]);          //L arrow
       arrow[2].setSizeXYZ(0.0,0.0,Lzppz[istep]);                   //Lzppz
       arrow[3].setSizeXYZ(Lzppx[istep],Lzppy[istep],Lzppz[istep]); //Lzpp arrow
       trace[0].addPoint(Lx[istep], Ly[istep], Lz[istep]);          //L trace
       trace[1].addPoint(Lzppx[istep], Lzppy[istep], Lzppz[istep]); //Lz trace
       frame3d.addElement(arrow[0]);
       frame3d.addElement(arrow[1]);
       frame3d.addElement(arrow[2]);
       frame3d.addElement(arrow[3]);
       frame3d.addElement(trace[0]);
       frame3d.addElement(trace[1]);
       frame3d.render();
       istep++;
     } else {
         stopAnimation();
         plot0.setPreferredMinMax(0.0,tmax,0.0,1.5*ArrayLib.max(w3));
         plot0.append(0,t,w1); //phi
         plot0.append(1,t,w3); //theta
         plot0.append(2,t,w5); //psi
         plot0.render();
         plot1.setPreferredMinMax(0.0,tmax,ArrayLib.min(w4),10.*pi/tmax);
         plot1.append(0,t,w2); //phi_d
         plot1.append(1,t,w4); //theta_d
         plot1.append(2,t,w6); //psi_d
         plot1.render();
         plot2.setPreferredMinMax(ArrayLib.min(w3),ArrayLib.max(w3),
                     (1.-0.005)*ArrayLib.min(Vef),(1.+0.005)*ArrayLib.max(Vef));
         plot2.append(0,w3,Vef); //Potential
         plot2.append(1,w3,Ep);  //E'
         plot2.render();
         myControl.println("animation: black trace and arrow - L vector");
         myControl.println("animation: blue trace and arrow - Lz'' vector ");
         myControl.println("animation: red arrow - Lz vector");
         myControl.println("animation: magenta arrow - Lz''_z vector");
         myControl.println("top-right: phi-black, theta-blue, psi-red");
         myControl.println("bottom-left: phi_d-black, theta_d-blue, psi_d-red");
         myControl.println("bottom-right: E'-blue, Veff-green");
         myControl.println("E'="+nf.format(Ep[0])+", theta_1="+nf.format(w3[0])
                                           +", theta_2="+nf.format(w3[NPTS-1]));
     }
   }

   public void calculate() {
     //w1:phi, w2:phi_dot, w3:theta, w4:theta_dot, w5:psi, w6:psi_dot,t=time
     for(int i=1; i<NPTS; i++){
       odesolver.step();
       w1[i]=state[0];//phi
       w2[i]=state[1];//phi_dot
       w3[i]=state[2];//theta
       w4[i]=state[3];//theta_dot
       w5[i]=state[4];//psi
       w6[i]=ws-w2[i]*Math.cos(w3[i]);   //psi_dot
       t[i]=state[5];
       //Angular momentum components Lx, Ly, Lz, L''z, L''y, L''z
       Lx[i]=I*w4[i]*Math.cos(w1[i])-(I*w2[i]*Math.cos(w3[i])-Is*ws)*
           Math.sin(w3[i])*Math.sin(w1[i]);
       Ly[i]=I*w4[i]*Math.sin(w1[i])+(I*w2[i]*Math.cos(w3[i])-Is*ws)*
           Math.sin(w3[i])*Math.cos(w1[i]);
       Lz[i]=I*w2[i]*Math.sin(w3[i])*Math.sin(w3[i])+Is*ws*Math.cos(w3[i]);
       //Lz" along x, y, z
       Lzppx[i]=Is*ws*Math.sin(w3[i])*Math.sin(w1[i]);
       Lzppy[i]=-Is*ws*Math.sin(w3[i])*Math.cos(w1[i]);
       Lzppz[i]=Is*ws*Math.cos(w3[i]);
       //Effective Potential
       Vef[i]=M*g*a*Math.cos(w3[i])+(Lz[i]-Lzppz[i])*(Lz[i]-Lzppz[i])/
             (2*I*Math.sin(w3[i])*Math.sin(w3[i]));
       Ekp[i]=0.5*I*w4[i]*w4[i];
       Ep[i]=Ekp[i]+Vef[i];                //kinetic, and prime energy
     }
   }

  class top_der implements ODE {
    public double[] getState() {
    return state;
    }

    public void getRate(double[] state, double[] rate) {
      //state(0):w1=phi, state(1):w2=phi_dot, state(2):w3=theta
      //state[3]:w4=theta_dot, state(4):w5=psi, state(5):t=time
      //main program produces w6:psi_dot, and time
      rate[0]=state[1];
      rate[1]=(Is*ws-
               2*I*state[1]*Math.cos(state[2]))*state[3]/(I*Math.sin(state[2]));
      rate[2]=state[3];
      rate[3]=(tau0-(Is*ws-
           I*state[1]*Math.cos(state[2]))*state[1])*Math.sin(state[2])/I;
      rate[4]=ws-state[1]*Math.cos(state[2]);
      rate[5]=1.0;    //dt/dt
    }
  }

   public void clear  () {
     trace[0].clear();
     trace[1].clear();
     frame3d.addElement(trace[0]);
     frame3d.addElement(trace[1]);
     //panelb.addElement(trace[2]);
     for(int i=0; i< Iarrows; i++){
       arrow[i].setXYZ(0, 0, 0); //x,y,z arrow position
       //arrow length corresponds to the vector function x,y,z components
       arrow[i].setSizeXYZ(0, 0, 0);
       //clear the arrows from the frame
       frame3d.addElement(arrow[i]);
     }
     frame3d.setPreferredMinMax(0,0,0,0,0,0);
     frame3d.render();
     //panelb.setPreferredMinMax(0,0,0,0,0,0);
     //panelb.render();
     //plot.clearData();
  }

   public void resetAnimation() {
     istep=0;                     //doStep index counter
     clear();
     myControl.clearMessages();
     myControl.println ("program to solve Euler's equations, plot");
     myControl.println ("them and their respective angular speeds,");
     myControl.println ("as well as angular momentum properties.");
     myControl.println ("See the text for details. Press 'Initialize'");
     myControl.println ("to begin, or 'reset' if needed.");
     w1    = null;
     w2    = null;
     w3    = null;
     w4    = null;
     w5    = null;
     w6    = null;
     Lx    = null;
     Ly    = null;
     Lz    = null;
     Lzppx = null;
     Lzppy = null;
     Lzppz = null;
     Vef   = null;
     Ekp   = null;
     Ep    = null;
     t     = null;
     I=1.0; Is=1.5;            //Inertia
     g=9.8; M=1.; a=0.1;       //gravity, mass, lever arm
     tau0=M*g*a;               //torque value at theta=0
     ph0=0.; th0=0.6*pi/2.+0.01; ps0=0.;  //initial angles in rad
     ph0d=0.0; th0d=0.0;       //initial angular speeds in rad/s
     ws=1.5+Math.sqrt(4.*I*tau0*Math.cos(th0))/Is;
     tmax=10;
     NPTS=200;
     delayTime=40;           //time in between animation steps (Abstract Animation)
     myControl.setValue("I",I);
     myControl.setValue("M",M);
     myControl.setValue("a",a);
     myControl.setValue("Is",Is);
     myControl.setValue("g",g);
     myControl.setValue("phi_0 (rad)",ph0);
     myControl.setValue("theta_0 (rad)",nf.format(th0));
     myControl.setValue("psi_0 (rad)",ps0);
     myControl.setValue("phi_0d (rad/s)",ph0d);
     myControl.setValue("theta_0d (rad/s)",th0d);
     myControl.setValue("ws (rad/s)",nf.format(ws));
     myControl.setValue("tmax(s)",tmax);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",nf.format(dt));
   }

   public void initializeAnimation() {
     clear();
     istep=0;             //increased in doStep()
     myControl.clearMessages();
     I=myControl.getDouble     ("I");
     M=myControl.getDouble     ("M");
     a=myControl.getDouble     ("a");
     Is=myControl.getDouble    ("Is");
     g=myControl.getDouble     ("g");
     ph0=myControl.getDouble   ("phi_0 (rad)");
     th0=myControl.getDouble   ("theta_0 (rad)");
     ps0=myControl.getDouble   ("psi_0 (rad)");
     ph0d=myControl.getDouble  ("phi_0d (rad/s)");
     th0d=myControl.getDouble  ("theta_0d (rad/s)");
     ws=myControl.getDouble    ("ws (rad/s)");
     tmax=myControl.getDouble  ("tmax(s)");
     NPTS=myControl.getInt     ("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",nf.format(dt));
     w1    = new double[NPTS];
     w2    = new double[NPTS];
     w3    = new double[NPTS];
     w4    = new double[NPTS];
     w5    = new double[NPTS];
     w6    = new double[NPTS];
     Lx    = new double[NPTS];
     Ly    = new double[NPTS];
     Lz    = new double[NPTS];
     Lzppx = new double[NPTS];
     Lzppy = new double[NPTS];
     Lzppz = new double[NPTS];
     Vef   = new double[NPTS];
     Ekp   = new double[NPTS];
     Ep    = new double[NPTS];
     t     = new double[NPTS];
     //initial conditions
     //w1:phi, w2:phi_dot, w3:theta, w4:theta_dot, w5:psi, w6:psi_dot,t=time
     tau0=M*g*a;               //torque value at theta=0
     state[0]=ph0;
     state[1]=ph0d;
     state[2]=th0;
     state[3]=th0d;
     state[4]=ps0;
     state[5]=0.0;     //t0=0.
     w1[istep]=state[0];
     w2[istep]=state[1];
     w3[istep]=state[2];
     w4[istep]=state[3];
     w5[istep]=state[4];
     w6[istep]=ws-w2[istep]*Math.cos(w3[istep]);
     t[istep]=state[5];
     odesolver.initialize(dt); // step size
     //Angular momentum components Lx, Ly, Lz, L''z, L''y, L''z
     Lx[istep]=I*w4[istep]*Math.cos(w1[istep])-(I*w2[istep]*Math.cos(w3[istep])-
               Is*ws)*Math.sin(w3[istep])*Math.sin(w1[istep]);
     Ly[istep]=I*w4[istep]*Math.sin(w1[istep])+(I*w2[istep]*Math.cos(w3[istep])-
               Is*ws)*Math.sin(w3[istep])*Math.cos(w1[istep]);
     Lz[istep]=I*w2[istep]*Math.sin(w3[istep])*Math.sin(w3[istep])+
               Is*ws*Math.cos(w3[istep]);
     //Lz" along x, y, z
     Lzppx[istep]=Is*ws*Math.sin(w3[istep])*Math.sin(w1[istep]);
     Lzppy[istep]=-Is*ws*Math.sin(w3[istep])*Math.cos(w1[istep]);
     Lzppz[istep]=Is*ws*Math.cos(w3[istep]);
     //Effective Potential
     Vef[istep]=M*g*a*Math.cos(w3[istep])+
                (Lz[istep]-Lzppz[istep])*(Lz[istep]-Lzppz[istep])/
               (2*I*Math.sin(w3[istep])*Math.sin(w3[istep]));
     Ekp[istep]=0.5*I*w4[istep]*w4[istep];
     Ep[istep]=Ekp[istep]+Vef[istep];       //kinetic, and prime energy
     calculate();
     double va1x, va1y, va1z, va2x, va2y, va2z;
     double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z;
     double v1x,v1y, v1z, v2x, v2y, v2z;
     va1x=ArrayLib.min(Lx);va1y=ArrayLib.min(Ly);va1z=ArrayLib.min(Lz); //get largest # for window view
     va2x=ArrayLib.max(Lx);va2y=ArrayLib.max(Ly);va2z=ArrayLib.max(Lz);
     vb1x=ArrayLib.min(Lzppx);vb1y=ArrayLib.min(Lzppy);vb1z=ArrayLib.min(Lzppz);
     vb2x=ArrayLib.max(Lzppx);vb2y=ArrayLib.max(Lzppy);vb2z=ArrayLib.max(Lzppz);
     v1x=Math.min(va1x,vb1x);v1y=Math.min(va1y,vb1y);v1z=Math.min(va1z,vb1z);
     v2x=Math.max(va2x,vb2x);v2y=Math.max(va2y,vb2y);v2z=Math.max(va2z,vb2z);
     double [] vv=new double[] {Math.abs(v1x),v2x,Math.abs(v1y),v2y,Math.abs(v1z),v2z};
     double v=ArrayLib.max(vv);
     xmin=-v; ymin=-v; zmin=0.0; xmax=v; ymax=v; zmax=v;
     frame3d.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
     frame3d.getCamera().setXYZ(1.25*(xmax-xmin),1.25*(ymax-ymin),1.5*(zmax-zmin));
     frame3d.setAltitude(0.4);
     frame3d.setAzimuth(-0.7);
     myControl.println("Press 'start' to begin animation");
     myControl.println("Feel free to mouse-rotate the 3D display.");
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation model = new topApp();
     AnimationControl myControl;
     myControl = new AnimationControl(model);
     myControl.setLocation(540, 5);
     myControl.setSize(250,595);
     myControl.setDividerLocation(320);
     model.setControl(myControl);
   }
}