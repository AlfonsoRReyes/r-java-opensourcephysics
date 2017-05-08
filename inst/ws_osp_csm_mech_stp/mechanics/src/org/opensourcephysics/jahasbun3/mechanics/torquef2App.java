/*
@Author J E Hasbun 2007.
Solves Euler's equations for an ellipsoid without torques.
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
import org.opensourcephysics.frames.*;
import java.text.*;

public class torquef2App extends AbstractAnimation {
  Display3DFrame panela = new Display3DFrame("Ellipsoid System w, L animation");
  Display3DFrame panelb = new Display3DFrame("Ellipsoid System w1, w2, w3 space");
  PlotFrame plot=new PlotFrame
  ("t","w1,w2,w3","black-w1(t), blue-w2(t), red-w3(t)");
  private Control myControl;
  double w1[], w2[], w3[], t[], L1[], L2[], L3[];
  double tmax, dt, pi, cf, w10, w20, w30, Iv;
  double I1, I2, I3, gam1, gam2, gam3;
  double xmina,xmaxa,ymina,ymaxa,zmina,zmaxa;
  double xminb,xmaxb,yminb,ymaxb,zminb,zmaxb;
  int NPTS, istep, Iarrows=2;
  double[] state= new double [4];
  ODE ode = new torquef2_der();
  ODESolver odesolver = new RK4(ode);
  //ODESolver odesolver = new RK45MultiStep(ode);
  ElementTrail [] trace=new ElementTrail[3];
  ElementArrow [] arrow=new ElementArrow [Iarrows];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public torquef2App(){
    //animation panel of w, L
    panela.setDecorationType(VisualizationHints.DECORATION_AXES);
    panela.setLocation(35,5);
    panela.setSize(450,250);
    //w - space results panel
    panelb.setDecorationType(VisualizationHints.DECORATION_CUBE);
    panelb.setLocation(275,260);
    panelb.setSize(265,295);
    //individual w plots
    plot.setConnected(true);//set all connected
    plot.setLocation(5,260);
    plot.setSize(265,295);
    pi=Math.PI; cf=180./pi;
    nf.setMaximumFractionDigits(4);
    trace[0]=new ElementTrail();
    trace[1]=new ElementTrail();
    trace[2]=new ElementTrail();
    trace[0].setConnected(false);                           //w
    trace[0].getStyle().setLineColor(java.awt.Color.black);
    trace[1].setConnected(false);                           //L
    trace[1].getStyle().setLineColor(java.awt.Color.blue);
    trace[2].getStyle().setLineColor(java.awt.Color.black); //w-space by it self
    trace[2].getStyle().setLineWidth(1.5f);
    for(int i=0; i< Iarrows; i++){
      arrow[i]=new ElementArrow();
    }
    arrow[0].getStyle().setLineColor(java.awt.Color.black); //w
    arrow[1].getStyle().setLineColor(java.awt.Color.blue);  //L
  }

   protected void doStep() {
     if (istep < NPTS){
       trace[0].addPoint(w1[istep], w2[istep], w3[istep]);    //w trace
       trace[1].addPoint(L1[istep], L2[istep], L3[istep]);    //L trace
       trace[2].addPoint(w1[istep], w2[istep], w3[istep]);    //w trace again
       arrow[0].setSizeXYZ(w1[istep],w2[istep],w3[istep]);    //w arrow
       arrow[1].setSizeXYZ(L1[istep],L2[istep],L3[istep]);    //L arrow
       panela.addElement(trace[0]);
       panela.addElement(trace[1]);
       panela.addElement(arrow[0]);
       panela.addElement(arrow[1]);
       panela.render();
       istep++;
     } else {
         stopAnimation();
         panelb.addElement(trace[2]);
         panelb.render();
         plot.setLineColor(0,java.awt.Color.black);
         plot.setMarkerSize(0,0);
         plot.append(0, t, w1); //w1
         plot.setLineColor(1,java.awt.Color.blue);
         plot.setMarkerSize(1,0);
         plot.append(1, t, w2); //w2
         plot.setLineColor(2,java.awt.Color.red);
         plot.setMarkerSize(2,0);
         plot.append(2, t, w3); //w3
         plot.render();
         myControl.println("animation: black trace and arrow - w vector S-frame");
         myControl.println("animation: blue trace and arrow - L vector ");
         myControl.println("Lower left: w1(black), w2(blue), w3(red) vs t");
         myControl.println("Lower right: w1, w2, w3 phase plot");
         myControl.println("Press 'Stop' before beginning a new animation");
         myControl.println("gamma1="+nf.format(gam1)+", gamma2="+nf.format(gam2)
                           +", gamm3="+nf.format(gam3));
     }
   }

   public void calculate() {
     //w, L calculation
     for(int i=1; i<NPTS; i++){
       odesolver.step();
       w1[i]=state[0];
       w2[i]=state[1];
       w3[i]=state[2];
       t[i]=state[3];
       L1[i]=I1*w1[i]/Iv;
       L2[i]=I2*w2[i]/Iv;
       L3[i]=I3*w3[i]/Iv;
     }
  }

     class torquef2_der implements ODE {
       public double[] getState() {
         return state;
       }

       public void getRate(double[] state, double[] rate) {
         // state(0):w1, state(1):w2, state(2):w3, state[3]:time
         rate[0]=-state[1]*state[2]*gam1;
         rate[1]=-state[2]*state[0]*gam2;
         rate[2]=-state[0]*state[1]*gam3;
         rate[3]=1;    //dt/dt
       }
     }

   public void clear  () {
     trace[0].clear();
     trace[1].clear();
     trace[2].clear();
     panela.addElement(trace[0]);
     panela.addElement(trace[1]);
     panelb.addElement(trace[2]);
     for(int i=0; i< Iarrows; i++){
       arrow[i].setXYZ(0, 0, 0); //x,y,z arrow position
       //arrow length corresponds to the vector function x,y,z components
       arrow[i].setSizeXYZ(0, 0, 0);
       //clear the arrows from the frame
       panela.addElement(arrow[i]);
     }
     panela.setPreferredMinMax(0,0,0,0,0,0);
     panela.render();
     panelb.setPreferredMinMax(0,0,0,0,0,0);
     panelb.render();
     plot.clearData();
  }

   public void resetAnimation() {
     istep=0;                     //doStep index counter
     clear();
     myControl.clearMessages();
     myControl.println ("Solves Euler's equations for an ellipsoid");
     myControl.println ("without torques. The differential equations");
     myControl.println ("are: dw1/dt=-gamma1*w2*w3, dw2/dt=-gamma2*w1*w3,");
     myControl.println ("and dw3/dt=-gamma3*w1*w2, where gamma1=(I3-I2)/I1,");
     myControl.println ("gamma2=(I1-I3)/I2, gamma3=(I2-I1)/I3, and the");
     myControl.println ("angular momentum is and Li=Ii*wi.");
     myControl.println("Press 'Initialize' before beginning the animation");
     w1    = null;
     w2    = null;
     w3    = null;
     L1    = null;
     L2    = null;
     L3    = null;
     t     = null;
     I1=145./46.; I2=1589./251.; I3=1242./152.;  //principal moments of inertia
     gam1=(I3-I2)/I1;         //gamma1
     gam2=(I1-I3)/I2;         //gamma2
     gam3=(I2-I1)/I3;         //gamma2
     w10=1.; w20=1.; w30=1.;  //initial conditions
     tmax=12;
     NPTS=150;
     delayTime=40;           //time in between animation steps (Abstract Animation)
     myControl.setValue("I1",nf.format(I1));
     myControl.setValue("I2",nf.format(I2));
     myControl.setValue("I3",nf.format(I3));
     myControl.setValue("w10",w10);
     myControl.setValue("w20",w20);
     myControl.setValue("w30",w30);
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
     I1=myControl.getDouble("I1");
     I2=myControl.getDouble("I2");
     I3=myControl.getDouble("I3");
     w10=myControl.getDouble("w10");
     w20=myControl.getDouble("w20");
     w30=myControl.getDouble("w30");
     tmax=myControl.getDouble("tmax(s)");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",nf.format(dt));
     w1    = new double[NPTS];
     w2    = new double[NPTS];
     w3    = new double[NPTS];
     L1    = new double[NPTS];
     L2    = new double[NPTS];
     L3    = new double[NPTS];
     t     = new double[NPTS];
     //initial conditions
     gam1=(I3-I2)/I1;         //gamma1
     gam2=(I1-I3)/I2;         //gamma2
     gam3=(I2-I1)/I3;         //gamma2
     Iv=(I1+I2+I3)/3.;        //Average Moment
     state[0]=w10;
     state[1]=w20;
     state[2]=w30;
     state[3]=0.0;            //t0=0.
     w1[istep]=state[0];
     w2[istep]=state[1];
     w3[istep]=state[2];
     t[istep]=state[3];
     odesolver.initialize(dt); // step size
     calculate();
     double v1a, v2a, v3a, v1b, v2b, v3b;
     v1a=ArrayLib.max(L1); v2a=ArrayLib.max(L1); v3a=ArrayLib.max(L1);
     v1b=ArrayLib.max(w1); v2b=ArrayLib.max(w2); v3b=ArrayLib.max(w3);
     xmina=-Math.max(v1a,v1b);
     xmaxa=-xmina;
     ymina=-Math.max(v2a,v2b);
     ymaxa=-ymina;
     zmina=-Math.max(v3a,v3b);
     zmaxa=-zmina;
     panela.setPreferredMinMax(xmina, xmaxa, ymina, ymaxa, zmina, zmaxa);
     panela.getCamera().setXYZ(0.85*(xmaxa-xmina),0.85*(ymaxa-ymina),1.6*(zmaxa-zmina));
     panela.setAltitude(0.3);
     panela.setAzimuth(0.5);
     xminb=ArrayLib.min(w1);
     xmaxb=ArrayLib.max(w1);
     yminb=ArrayLib.min(w2);;
     ymaxb=ArrayLib.max(w2);
     zminb=ArrayLib.min(w3);;
     zmaxb=ArrayLib.max(w3);
     panelb.setPreferredMinMax(xminb, xmaxb, yminb, ymaxb, zminb, zmaxb);
     panelb.getCamera().setXYZ(1.3*(xmaxb-xminb),1.3*(ymaxb-yminb),1.3*(zmaxb-zminb));
     panelb.setAltitude(0.3);
     panelb.setAzimuth(0.6);
     myControl.println("Press 'start' to begin animation");
     myControl.println("Feel free to mouse-rotate the 3D display ");
     myControl.println("frame for better view.");
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation model = new torquef2App();
     AnimationControl myControl;
     myControl = new AnimationControl(model);
     myControl.setLocation(540, 5);
     myControl.setSize(250,555);
     myControl.setDividerLocation(235);
     model.setControl(myControl);
   }
}