/*
@Author J E Hasbun 2007.
Binary star system solved numerically.
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
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

  public class binary2App extends AbstractAnimation {
  PlottingPanel panel = new PlottingPanel
                          ("x(AU)","y(AU)","Binary System - numerical");
  DrawingFrame  frame = new DrawingFrame(panel);
  PlotFrame plot1= new PlotFrame("t(yr)","r1,r2(AU)","r vs t - numeric binary");
  PlotFrame plot2= new PlotFrame("t(yr)","v1,v2(Au/yr)","v vs t - numeric binary");
  Trail [] trail=new Trail[2];
  Circle [] point=new Circle[3];
  ODE ode = new binary2_der();
  ODESolver odesolver = new RK4(ode);
  private Control myControl;
  double [] x1,y1,x2,y2,r1,r2,v1,v2,t;
  double[] state= new double [9];
  double m,m1,m2,e,rcm,tau,dt,t0,tmax,d;
  double x10,y10,x20,y20,x00;
  double v1x0,v1y0,v2x0,v2y0,rmin1,rmin2,rmax1,rmax2,a1,a2,a;
  double pi,c;
  double xmin,xmax,ymin,ymax,af;
  int NPTS,istep;
  NumberFormat nf = NumberFormat.getInstance();

  public binary2App (){
    frame.setLocation(5,10);
    frame.setSize(430,290);
    frame.setTitle("Binary System y vs x numerical soln");
    plot1.setLocation(5,300);
    plot1.setSize(215,290);
    plot2.setLocation(220,300);
    plot2.setSize(215,290);
    nf.setMaximumFractionDigits(3);
    pi=Math.PI; c=4.*pi*pi;
    trail[0]=new Trail();
    trail[1]=new Trail();
  }

   protected void doStep() {
     if(istep<NPTS){
       //mass 1
       trail[0].addPoint(x1[istep],y1[istep]);
       point[0].setXY(x1[istep],y1[istep]);
       panel.addDrawable(trail[0]);
       panel.addDrawable(point[0]);
       //mass 2
       trail[1].addPoint(x2[istep],y2[istep]);
       point[1].setXY(x2[istep],y2[istep]);
       panel.addDrawable(trail[1]);
       panel.addDrawable(point[1]);
       checkPanel();
       istep++;   //increase the count by one
       calculate();//calculate the next step
     }else{
     //place the center of mass point also
     point[2].setXY(rcm,0);
     panel.addDrawable(point[2]);
     frame.render();
     doOtherPlots();
     stopAnimation();
     rmin1=ArrayLib.min(r1);  rmax1=ArrayLib.max(r1);
     rmin2=ArrayLib.min(r2);  rmax2=ArrayLib.max(r2);
     a1=(rmin1+rmax1)/2;//semimajor axes
     a2=(rmin2+rmax2)/2;
     a=a1+a2;           //the reduced mass semimajor axis
     tau=Math.sqrt(Math.pow(a,3.0)/m); //the system's period
     myControl.println("blue:m1, red:m2, grey:cm");
     myControl.println("(rmin1,rmax1)=\t("+nf.format(rmin1)+","+nf.format(rmax1)+") AU");
     myControl.println("(rmin2,rmax2)=\t("+nf.format(rmin2)+","+nf.format(rmax2)+") AU");
     myControl.println("(a1,a2,a)=\t("+nf.format(a1)+","+nf.format(a2)+
                       ","+nf.format(a)+") AU");
     myControl.println("tau=\t"+nf.format(tau)+" tau0");
     myControl.println("(m,m1,m2)=\t("+m+","+m1+","+m2+")Ms");
     }
   }
   public void calculate () {
     //orbits calculation
     //state[0,1,2,3,4,5,6,7,8]=x1,v1x,y1,v1y,x2,v2x,y2,v2y,t
     t[istep]=ode.getState()[8];
     odesolver.step();
     x1[istep]=ode.getState()[0];
     y1[istep]=ode.getState()[2];
     x2[istep]=ode.getState()[4];
     y2[istep]=ode.getState()[6];
     r1[istep]=Math.sqrt(x1[istep]*x1[istep]+y1[istep]*y1[istep]);
     r2[istep]=Math.sqrt(x2[istep]*x2[istep]+y2[istep]*y2[istep]);
     double v1x=ode.getState()[1], v1y=ode.getState()[3];
     double v2x=ode.getState()[5], v2y=ode.getState()[7];
     v1[istep]=Math.sqrt(v1x*v1x+v1y*v1y);
     v2[istep]=Math.sqrt(v2x*v2x+v2y*v2y);
   }

   public void doOtherPlots () {
     plot1.setConnected(0,true);
     plot1.setLineColor(0,java.awt.Color.blue);
     plot1.setMarkerSize(0,0);
     plot1.append(0,t,r1);
     plot1.setConnected(1,false);
     plot1.setMarkerSize(1,1);
     plot1.setMarkerShape(1,1);//circle
     plot1.setMarkerColor(1,java.awt.Color.red);
     plot1.append(1,t,r2);
     plot1.render();
     plot2.setConnected(0,true);
     plot2.setLineColor(0,java.awt.Color.blue);
     plot2.setMarkerSize(0,0);
     plot2.append(0,t,v1);
     plot2.setConnected(1,false);
     plot2.setMarkerSize(1,1);
     plot2.setMarkerShape(1,1);//circle
     plot2.setMarkerColor(1,java.awt.Color.red);
     plot2.append(1,t,v2);
     plot2.render();
     //for (int i=0; i<NPTS;i++){
     //  myControl.println("i,r1,r2,v1,v2="+i+": "+r1[i]+", "+r2[i]+", "+v1[i]+", "+v2[i]);
     //}
   }

   public void checkPanel() {
     if(x1[istep]<xmin||x2[istep]<xmin){xmin=Math.min(x1[istep],x2[istep]);}
     if(x1[istep]>xmax||x2[istep]>xmax){xmax=Math.max(x1[istep],x2[istep]);}
     if(y1[istep]<ymin||y2[istep]<ymin){ymin=Math.min(y1[istep],y2[istep]);}
     if(y1[istep]>ymax||y2[istep]>ymax){ymax=Math.max(y1[istep],y2[istep]);}
     panel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     frame.render();
    }

   class binary2_der implements ODE {
     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2,3,4,5,6,7,8]=x1,v1x,y1,v1y,x2,v2x,y2,v2y,t
       //rates are the derivatives of the state, for example
       //rate[0]=dr/dt=v->state(1), and dv/dt=c*m2*
       //(state(4)-state(0))/((state(4)-state(0))^2+(state(6)-state(2))^2)^(3/2),
       //etc
       double wr=Math.pow((Math.pow(state[4]-state[0],2.0)+
                           Math.pow(state[6]-state[2],2.0)),3.0/2.0);
       rate[0]=state[1];
       rate[1]=c*m2*(state[4]-state[0])/wr;
       rate[2]=state[3];
       rate[3]=c*m2*(state[6]-state[2])/wr;
       rate[4]=state[5];
       rate[5]=-c*m1*(state[4]-state[0])/wr;
       rate[6]=state[7];
       rate[7]=-c*m1*(state[6]-state[2])/wr;
       rate[8]=1; //dt/dt
     }
   }

   public void clear () {
       trail[0].clear();
       trail[1].clear();
       panel.clear();
       frame.render();
   }

   public void resetAnimation () {
     clear();
     myControl.clearMessages();
     myControl.println ("Binary star system solved numerically. The system to be");
     myControl.println ("solved is:dx1/dt=v1x, dv1x/dt=4pi^2*m2*(x2-x1)/r^3,");
     myControl.println ("where r=((x2-x1)^2+(y2-y1)^2)^(1/2), etc. The binary system");
     myControl.println ("turns out to involve four such pairs of DE's. The center");
     myControl.println ("of mass is set at zero. Other ic's are shown in the input");
     myControl.println ("values. The eccentricity only plays a role to set other");
     myControl.println ("initial values to compare with the results of binary1App.");
     x1    = null;
     x2    = null;
     y1    = null;
     y2    = null;
     r1    = null;
     r2    = null;
     v1    = null;
     v2    = null;
     t     = null;
     m=5; m1=3; m2=m-m1;   //initial total and individual masses in Ms
     e=0.6;                //eccentricity
     d=Math.sqrt(1+e);     //e>0 =>deviation from circular orbit
     x10=-1;  y10=0;       //initial position of first mass
     x20=-m1*x10/m2; y20=0;//initial position of second mass
     x00=x20-x10;          //relative coordinate circular orbit parameter
     //initial speed value components for each mass follow
     v1x0=0; v1y0=d*(-m2*2*pi*Math.sqrt(m/x00)/m);
     v2x0=0; v2y0=-m1*v1y0/m2;
     t0=0; tmax=6.988;     //in units of tau0 (tau0=1 year)
     NPTS=200;
     delayTime=25;         //time in between animation steps (Abstract Animation)
     myControl.setValue("x10",x10);
     myControl.setValue("x20",x20);
     myControl.setValue("y10",y10);
     myControl.setValue("y20",y20);
     myControl.setValue("v1x0",v1x0);
     myControl.setValue("v1y0",v1y0);
     myControl.setValue("v2x0",v2x0);
     myControl.setValue("v2y0",v2y0);
     myControl.setValue("m",m);
     myControl.setValue("m1",m1);
     myControl.setValue("eccentricity",e);
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     dt=(tmax-t0)/(NPTS-1);
     myControl.setValue("dt((tmax-t0)/(NPTS-1))=",dt);
     m2=m-m1;
     myControl.setValue("m2(m-m1)=",m2);
   }

   public void initializeAnimation() {
     //initial conditions
     istep=0;  //so that istep+1 in doStep =0 the first time
     clear();
     myControl.clearMessages();
     x10=myControl.getDouble("x10");
     x20=myControl.getDouble("x20");
     y10=myControl.getDouble("y10");
     y20=myControl.getDouble("y20");
     v1x0=myControl.getDouble("v1x0");
     v1y0=myControl.getDouble("v1y0");
     v2x0=myControl.getDouble("v2x0");
     v2y0=myControl.getDouble("v2y0");
     m=myControl.getDouble("m");
     m1=myControl.getDouble("m1");
     e=myControl.getDouble("eccentricity");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=(tmax-t0)/(NPTS-1);
     myControl.setValue("dt((tmax-t0)/(NPTS-1))=",dt);
     m2=m-m1;
     myControl.setValue("m2(m-m1)=",m2);
     x1     = new double[NPTS+1];
     x2     = new double[NPTS+1];
     y1     = new double[NPTS+1];
     y2     = new double[NPTS+1];
     r1     = new double[NPTS+1];
     r2     = new double[NPTS+1];
     v1     = new double[NPTS+1];
     v2     = new double[NPTS+1];
     t      = new double[NPTS+1];
     rcm=0; //the center of mass
     af=(1+0.1);
     xmin=-Math.abs(x10);
     xmax=-xmin;
     ymin=-Math.abs(y10);
     ymax=-ymin;
     panel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     //panel.setPreferredMinMax(-6.,4.,-3.,3.);
     panel.setSquareAspect(true);
     trail[0].color=java.awt.Color.blue;
     trail[1].color=java.awt.Color.red;
     point[0]=new Circle();
     point[1]=new Circle();
     point[2]=new Circle();
     point[0].color=java.awt.Color.blue;
     point[0].pixRadius=3;
     point[1].color=java.awt.Color.red;
     point[1].pixRadius=3;
     point[2].color=java.awt.Color.gray;
     point[2].pixRadius=3;
     state[0]=x10;
     state[1]=v1x0;
     state[2]=y10;
     state[3]=v1y0;
     state[4]=x20;
     state[5]=v2x0;
     state[6]=y10;
     state[7]=v2y0;
     state[8]=t0;
     x1[istep]=state[0];
     y1[istep]=state[2];
     x2[istep]=state[4];
     y2[istep]=state[6];
     r1[istep]=Math.sqrt(x1[istep]*x1[istep]+y1[istep]*y1[istep]);
     r2[istep]=Math.sqrt(x2[istep]*x2[istep]+y2[istep]*y2[istep]);
     v1[istep]=Math.sqrt(v1x0*v1x0+v1y0*v1y0);
     v2[istep]=Math.sqrt(v2x0*v2x0+v2y0*v2y0);
     odesolver.initialize(dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new binary2App();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(445, 10);
     myControl.setSize(350,580);
     myControl.setDividerLocation(350);
     app.setControl(myControl);
   }
}