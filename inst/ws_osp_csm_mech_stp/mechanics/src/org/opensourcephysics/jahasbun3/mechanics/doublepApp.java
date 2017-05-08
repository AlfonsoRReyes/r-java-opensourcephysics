/*
@Author J E Hasbun 2007.
Program to solve the double pendulum equations of motion numerically and plot
their solutions.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.text.*;

public class doublepApp extends AbstractAnimation {
  PlottingPanel panel = new PlottingPanel  ("x1, x2","y1, y2","Double Pendulum");
  DrawingFrame frame = new DrawingFrame(panel);
  PlotFrame plot0=new PlotFrame ("t","theta_{1}, theta_{2}"," Angles versus time");
  PlotFrame plot1=new PlotFrame ("theta_{1}","theta_{2}"," theta2 vs theta1");
  PlotFrame plot2=new PlotFrame ("t","theta_{d1}, theta_{d2}"," Angles versus time");
  private Control myControl;
  double w1[], w2[], w3[], w4[], t[];
  double x1[], x2[], y1[], y2[], r1[], r2[];
  double tmax, dt, pi,av1,av2;
  double L1, L2, m1, m2, g, tau, th10, th20, th10d, th20d, v;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int NPTS, istep, Iarrows=2;
  double[] state= new double [5];
  ODE ode = new doublep_der();
  ODESolver odesolver = new RK4(ode);
  //ODESolver odesolver = new RK45MultiStep(ode);
  Trail [] trail=new Trail[2];
  Arrow [] arrow=new Arrow [Iarrows];
  Circle [] circle=new Circle[2];
  DrawableShape rectangle;
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public doublepApp(){
    //the animation frame
    //frame3d.setDecorationType(VisualizationHints.DECORATION_CUBE);
    frame.setLocation(5,5);
    frame.setSize(260,295);
    //The angles frame
    plot0.setConnected(true);//set all connected
    plot0.setLocation(270,5);
    plot0.setSize(260,295);
    plot0.setLineColor(0,java.awt.Color.black);//for theta1
    plot0.setMarkerSize(0,0);
    plot0.setLineColor(1,java.awt.Color.blue);  //for theta2
    plot0.setMarkerSize(1,0);
    //phase plot
    plot1.setConnected(true);//set all connected
    plot1.setLocation(5,305);
    plot1.setSize(260,295);
    plot1.setLineColor(0,java.awt.Color.black);//for theta2 vs theta_1
    plot1.setMarkerSize(0,0);
    //The rates frame
    plot2.setConnected(true);//set all connected
    plot2.setLocation(270,305);
    plot2.setSize(260,295);
    plot2.setLineColor(0,java.awt.Color.black);//for theta_1d
    plot2.setMarkerSize(0,0);
    plot2.setLineColor(1,java.awt.Color.blue);  //for theta_2d
    plot2.setMarkerSize(1,0);
    pi=Math.PI;
    nf.setMaximumFractionDigits(4);
    trail[0]=new Trail();
    trail[1]=new Trail();
    trail[0].color=java.awt.Color.black;
    trail[0].setDashedStroke(1,4);//dash thickness=1, & length=4
    trail[0].setConnected(true);
    trail[1].color=java.awt.Color.blue;
    trail[1].setDashedStroke(1,4);
    trail[1].setConnected(true);
    //instantiate arrows
    arrow[0]=new Arrow(0.0,0.0,0.0,0.0);
    arrow[1]=new Arrow(0.0,0.0,0.0,0.0);
    //instantiate circles
    circle[0]=new Circle();
    circle[1]=new Circle();
  }

   protected void doStep() {
     if (istep < NPTS){
       //first mass arrow and circle
       arrow[0].setXY(0.,v);
       arrow[0].setXlength(x1[istep]);
       arrow[0].setYlength(y1[istep]-v);
       circle[0].setXY(x1[istep],y1[istep]);
       trail[0].addPoint(x1[istep],y1[istep]);
       //2nd mass arrow and circle
       arrow[1].setXY(x1[istep],y1[istep]);
       arrow[1].setXlength(x2[istep]-x1[istep]);
       arrow[1].setYlength(y2[istep]-y1[istep]);
       circle[1].setXY(x2[istep],y2[istep]);
       trail[1].addPoint(x2[istep],y2[istep]);
       panel.addDrawable(arrow[0]);
       panel.addDrawable(circle[0]);
       panel.addDrawable(trail[0]);
       panel.addDrawable(arrow[1]);
       panel.addDrawable(circle[1]);
       panel.addDrawable(trail[1]);
       panel.addDrawable(rectangle);
       frame.render();
       istep++;
     } else {
         stopAnimation();
         //thetas vs t
         double wmax=Math.max(ArrayLib.max(w1),ArrayLib.max(w3));
         double wmin=Math.min(ArrayLib.min(w1),ArrayLib.min(w3));
         plot0.setPreferredMinMax(0.0,tmax,wmin,wmax);
         plot0.append(0,t,w1); //theta_1
         plot0.append(1,t,w3); //theta_2
         plot0.render();
         plot1.setPreferredMinMax(ArrayLib.min(w1),ArrayLib.max(w1),
                                  ArrayLib.min(w3),ArrayLib.max(w3));
         //phase plot
         plot1.append(0,w1,w3); //theta_2 versus theta_1
         plot1.render();
         //theta_dots
         wmax=Math.max(ArrayLib.max(w2),ArrayLib.max(w4));
         wmin=Math.min(ArrayLib.min(w2),ArrayLib.min(w4));
         plot2.setPreferredMinMax(0.0,tmax,wmin,wmax);
         plot2.append(0,t,w2); //theta_1d
         plot2.append(1,t,w4); //theta_2d
         plot2.render();
         myControl.println("animation: black trace and arrow - m1");
         myControl.println("animation: blue trace and arrow - m2");
         myControl.println("top-right: theta_1-black, theta_2-blue");
         myControl.println("bottom-left: phase plot theta_2 vs theta_1");
         myControl.println("bottom-right: theta_1d-black, theta_d2-blue");
         myControl.println("vs t.");
         myControl.println("Time average of the lengths:\nr1_av="+nf.format(av1)+
                           ", r2_av="+nf.format(av2));
     }
   }

   public void calculate() {
     //start at 1, case i=0, done in initialize
     for(int i=1; i<NPTS; i++){
       odesolver.step();
       //state[0]=w1:theta1, state[1]=w2:theta1_dot,
       //state[2]=w3:theta2, state[3]=w4:theta2_dot, state[4]=t
       w1[i]=state[0];
       w2[i]=state[1];
       w3[i]=state[2];
       w4[i]=state[3];
       t[i]=state[4];
       //Coordinates versus time
       x1[i]=L1*Math.sin(w1[i]);
       y1[i]=L1*Math.cos(w1[i]);
       x2[i]=x1[i]+L2*Math.sin(w3[i]);
       y2[i]=y1[i]+L2*Math.cos(w3[i]);
       //length conservation check if needed
       r1[i]=Math.sqrt(x1[i]*x1[i]+y1[i]*y1[i]);
       r2[i]=Math.sqrt(x2[i]*x2[i]+y2[i]*y2[i]);
       //support is at L1+L2, where y's are measured from-shift the y's for simulation purpose
       y1[i]=v-y1[i];
       y2[i]=v-y2[i];
     }
     av1=Util.computeAverage(r1,0,NPTS);//average length mass 1 check
     av2=Util.computeAverage(r2,0,NPTS);//average length mass 2 check
   }

  class doublep_der implements ODE {
    public double[] getState() {
    return state;
    }

    public void getRate(double[] state, double[] rate) {
     //state[0]=w1:theta1, state[1]=w2:theta1_dot,
     //state[2]=w3:theta2, state[3]=w4:theta2_dot, state[4]=t
     double tp=state[0]-state[2];
     double cs=Math.sin(tp), cc=Math.cos(tp);
     double ta=m2*L2*cc/(m1+m2)/L1;
     double tb=(m2*L2*state[3]*state[3]*cs+
                (m1+m2)*g*Math.sin(state[0]))/(m1+m2)/L1;
     double tc=(L1*tb*cc+L1*state[1]*state[1]*cs-
                g*Math.sin(state[2]))/(L2-ta*L1*cc);
     rate[0]=state[1];
     rate[1]=-ta*tc-tb;
     rate[2]=state[3];
     rate[3]=tc;
     rate[4]=1.0;    //dt/dt
     }
  }

   public void clear  () {
     trail[0].clear();
     trail[1].clear();

     for(int i=0; i< Iarrows; i++){
       arrow[i].setXY(0.,0.);
       arrow[i].setXlength(0.0);
       arrow[i].setYlength(0.0);
       arrow[i].setColor(java.awt.Color.white);
       arrow[i].setHeadSize(0.0f);
       circle[i].pixRadius=0;
     }
     frame.clearDrawables();
     frame.render();
  }

   public void resetAnimation() {
     istep=0;                     //doStep index counter
     clear();
     myControl.clearMessages();
     myControl.println ("program to solve the double pendulum");
     myControl.println ("equations of motion numerically and plot");
     myControl.println ("their solutions. Refer to the text for the");
     myControl.println ("formulas used. theta_1 and theta_2 are ");
     myControl.println ("masses 1 & 2 initial angles, theta_1d and ");
     myControl.println ("theta_2d are respective initial angular");
     myControl.println ("speeds. Press 'initialize' to begin");
     myControl.println ("or 'reset' if needed.");
     w1    = null;
     w2    = null;
     w3    = null;
     w4    = null;
     x1    = null;
     x2    = null;
     y1    = null;
     y2    = null;
     r1    = null;
     r2    = null;
     t     = null;
     L1=1.0; L2=2.0;           //lengths
     m1=1.0; m2=2.0;           //masses
     g=9.8;                    //gravity
     tau=Math.sqrt(g/(L1+L2)); //a time unit
     tmax=3.0*tau;             //max time
     NPTS=125;
     delayTime=40;           //time in between animation steps (Abstract Animation)
     th10=pi/4; th20=pi/3;   //init angles in rad
     th10d=0.0; th20d=0.0;   //init angle speeds in rad/s
     myControl.setValue("L1",L1);
     myControl.setValue("L2",L2);
     myControl.setValue("m1",m1);
     myControl.setValue("m2",m2);
     myControl.setValue("g",g);
     myControl.setValue("theta_1(rad)",nf.format(th10));
     myControl.setValue("theta_2(rad)",nf.format(th20));
     myControl.setValue("theta_1d(rad)",th10d);
     myControl.setValue("theta_2d (rad/s)",th20d);
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
     L1   =myControl.getDouble("L1");
     L2   =myControl.getDouble("L2");
     m1   =myControl.getDouble("m1");
     m2   =myControl.getDouble("m2");
     g    =myControl.getDouble("g");
     th10 =myControl.getDouble("theta_1(rad)");
     th20 =myControl.getDouble("theta_2(rad)");
     th10d=myControl.getDouble("theta_1d(rad)");
     th20d=myControl.getDouble("theta_2d (rad/s)");
     tmax =myControl.getDouble("tmax(s)");
     NPTS =myControl.getInt   ("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",nf.format(dt));
     w1    = new double[NPTS];
     w2    = new double[NPTS];
     w3    = new double[NPTS];
     w4    = new double[NPTS];
     x1    = new double[NPTS];
     x2    = new double[NPTS];
     y1    = new double[NPTS];
     y2    = new double[NPTS];
     r1    = new double[NPTS];
     r2    = new double[NPTS];
     t     = new double[NPTS];
     //initial conditions
     //state[0]=w1:theta1, state[1]=w2:theta1_dot,
     //state[2]=w3:theta2, state[3]=w4:theta2_dot, state[4]=t
     state[0]=th10;
     state[1]=th10d;
     state[2]=th20;
     state[3]=th20d;
     state[4]=0.0;     //t0=0.
     w1[istep]=state[0];
     w2[istep]=state[1];
     w3[istep]=state[2];
     w4[istep]=state[3];
     t[istep]=state[4];
     //Coordinates versus time
     x1[istep]=L1*Math.sin(w1[istep]);
     y1[istep]=L1*Math.cos(w1[istep]);
     x2[istep]=x1[istep]+L2*Math.sin(w3[istep]);
     y2[istep]=y1[istep]+L2*Math.cos(w3[istep]);
     //length conservation check
     r1[istep]=Math.sqrt(x1[istep]*x1[istep]+y1[istep]*y1[istep]);
     r2[istep]=Math.sqrt(x2[istep]*x2[istep]+y2[istep]*y2[istep]);
     //support is at L1+L2, where y's are measured from-shift the y's for simulation purpose
     v=L1+L2;
     y1[istep]=v-y1[istep]; //shift the y's for simulation purpose
     y2[istep]=v-y2[istep];
     odesolver.initialize(dt); // step size
     calculate();
     //arrows properties
     arrow[0].setColor(java.awt.Color.black);
     arrow[0].setHeadSize(2.0f);
     arrow[1].setColor(java.awt.Color.blue);
     arrow[1].setHeadSize(2.0f);
     //circle properties
     circle[0].color=java.awt.Color.black;
     circle[0].pixRadius=3;
     circle[1].color=java.awt.Color.blue;
     circle[1].pixRadius=3;
     //the pendulum support
     rectangle=DrawableShape.createRectangle(0.0,v,0.2,0.2);//size is 0.2
     rectangle.color=new java.awt.Color(255, 0, 0, 50);//transparent red
     //panel frame properies
     double x1min=ArrayLib.min(x1),y1min=ArrayLib.min(y1);
     double x1max=ArrayLib.max(x1),y1max=ArrayLib.max(y1);
     double x2min=ArrayLib.min(x2),y2min=ArrayLib.min(y2);
     double x2max=ArrayLib.max(x2),y2max=ArrayLib.max(y2);
     double [] vx=new double[] {x1min,x2min,x1max,x2max};
     double [] vy=new double[] {v,y1min,y2min,y1max,y2max};
     xmin=ArrayLib.min(vx); xmax=ArrayLib.max(vx);
     ymin=ArrayLib.min(vy); ymax=ArrayLib.max(vy);
     frame.setPreferredMinMax(xmin,xmax,ymin,ymax);
     myControl.println("Press 'start' to begin animation");
     myControl.println("Feel free to expand or move the");
     myControl.println("windows");
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation model = new doublepApp();
     AnimationControl myControl;
     myControl = new AnimationControl(model);
     myControl.setLocation(540, 5);
     myControl.setSize(250,595);
     myControl.setDividerLocation(290);
     model.setControl(myControl);
   }
}