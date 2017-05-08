/*
@Author J E Hasbun 2007.
Solves the variable mass rocket equation in the precense of a gravitational
force.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;
import java.text.*;

  public class rocketApp extends AbstractAnimation {
    PlottingPanel panel = new PlottingPanel
        ("x(a_{b})","y(a_{b})","Rocket Simulation");
    DrawingFrame frame = new DrawingFrame(panel);
    PlotFrame plot0= new PlotFrame("t(tau)","x,y,r (R)","Positions versus time");
    PlotFrame plot1= new PlotFrame("t(tau)","vx,vy,v (v_{0})","Speeds versus time");
    Trail trail=new Trail();
    Circle point=new Circle();  //used for the rocket
    DrawableShape earth;
    ODE ode = new rocket_der();
    ODESolver odesolver = new RK4(ode);
    //ODESolver odesolver = new RK45MultiStep(ode);
    private Control myControl;
    double[] state= new double [6];
    double [] t,x,y,vx,vy,r,v,m;
    double G,R,m0,M,tau,v0,F0,mi,mi0,ff,Thrust,u,an,th,ux,uy,mp;
    double alpha,mf,tfmax,t0,tmax,x0,y0,vx0,vy0,dt,Range;
    double mL,vmin,vmax,xmin,xmax,ymin,ymax,aL,rmin,rmax;
    double x10,y10,v1x0,v1y0;
    double pi,cf;
    int NPTS,istep,istep_max,iflag;
    NumberFormat nf = NumberFormat.getInstance();
    DecimalFormat df = new DecimalFormat("+0.00000E00");

  public rocketApp(){
    frame.setLocation(50,5);
    frame.setSize(300,250);
    frame.setTitle("Rocket Simulation");
    plot0.setLocation(5,260);
    plot0.setSize(235,300);
    plot0.setConnected(true);
    plot1.setLocation(240,260);
    plot1.setSize(235,300);
    plot1.setConnected(true);
    nf.setMaximumFractionDigits(3);
    df.setMaximumFractionDigits(3);
    pi=Math.PI; cf=pi/180.;
    trail=new Trail();
    trail.setDashedStroke(1,3); //trail: thickness=2, dash length=4
  }

   protected void doStep() {
     if(istep < istep_max){  //only plot up to when the rocket touches ground
       trail.addPoint(x[istep], y[istep]); //draw line every even point
       point.setXY(x[istep], y[istep]);
       //draw
       checkPanel();
       panel.addDrawable(trail);
       panel.addDrawable(point);
       frame.render();
       panel.setMessage("x="+df.format(x[istep]),0);
       panel.setMessage("y="+df.format(y[istep]),1);
       panel.setMessage("t="+Util.f2(t[istep]),3);
       istep++;
     } else{
       stopAnimation();
       plot0.setMarkerSize(0,0);
       plot0.setMarkerSize(1,0);
       plot0.setMarkerSize(2,0);
       plot0.setLineColor(0,java.awt.Color.black);
       plot0.setLineColor(1,java.awt.Color.blue);
       plot0.setLineColor(2,java.awt.Color.red);
       plot0.append(0,t,x);
       plot0.append(1,t,y);
       plot0.append(2,t,r);
       plot0.render();
       plot1.setMarkerSize(0,0);
       plot1.setMarkerSize(1,0);
       plot1.setMarkerSize(2,0);
       plot1.setLineColor(0,java.awt.Color.black);
       plot1.setLineColor(1,java.awt.Color.blue);
       plot1.setLineColor(2,java.awt.Color.red);
       plot1.append(0,t,vx);
       plot1.append(1,t,vy);
       plot1.append(2,t,v);
       plot1.render();
     }
   }

   public void calculate () {
     //state[0,1,2,3,4,5]=x,vx,y,vy,m,t
     double [] tmp=new double[3];
     for(int i=1; i<=NPTS; i++){
       odesolver.step();
       x[i] = ode.getState()[0];
       y[i] = ode.getState()[2];
       r[i] = Math.sqrt(x[i]*x[i]+y[i]*y[i]);
       vx[i] = ode.getState()[1];
       vy[i] = ode.getState()[3];
       v[i] = Math.sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
       t[i] = ode.getState()[5];
       if(iflag==0){
         tmp[0] = x[i]; tmp[1] = y[i]; tmp[2] = r[i];
         if (ArrayLib.min(tmp) < rmin) { rmin = ArrayLib.min(tmp); }
         if (ArrayLib.max(tmp) > rmax) { rmax = ArrayLib.max(tmp); }
         tmp[0] = vx[i]; tmp[1] = vy[i];  tmp[2] = v[i];
         if (ArrayLib.min(tmp) < vmin) { vmin = ArrayLib.min(tmp); }
         if (ArrayLib.max(tmp) > vmax) { vmax = ArrayLib.max(tmp); }
         if (r[i] < 0.99) { //used in DoStep to plot above ground point only
                            //ground is 1.0 R, include points slightly below
           istep_max = i;
           iflag = 1;
         }
       }
     //myControl.println ("i,x,y,t="+i+", "+x[i]+", "+y[i]+", "+t[i]+", istep_max="+istep_max);
     }
     //landing point => range (as measured from the initial point
     Range=R*Math.abs(Math.atan(y[istep_max]/x[istep_max])-Math.atan(y[0]/x[0]));
   }

   class rocket_der implements ODE {
     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2,3,4,5]=x,vx,y,vy,m,t
       //rates are the derivatives of the state
       double wr=Math.sqrt(state[0]*state[0]+state[2]*state[2]);
       double tmp=alpha*stepf(state[5],tfmax,1);
       rate[0]=state[1];
       rate[1]=ux*tmp/state[4]-M*state[0]/Math.pow(wr,3.);
       rate[2]=state[3];
       rate[3]=uy*tmp/state[4]-M*state[2]/Math.pow(wr,3.);
       rate[4]=-tmp;
       rate[5]=1;    //dt/dt
     }
   }

   public double stepf (double x, double s, int ss){
   // This is a step function stepf(x,s,ss)
   // x is the input array, s the value where the step occurs
   // if ss=1 the step occurs at 's' from 1 to 0 for all x
   // if ss=2 the step occurs at 's' from 0 to 1 for all x
   // if ss is neither 1 nor 2 the function returns the value of zero
   double step;
   if (ss==1){
     step=1./(1+Math.exp(100*(x-s)));
     }else{
       if (ss == 2){
         step = 1 - 1. / (1 + Math.exp(100 * (x - s)));
       } else {
         step = 0.0;
       }
     }
     return step;
   }

   public void checkPanel() {
     if(x[istep]<xmin){xmin=xmin*aL;}
     if(x[istep]>xmax){xmax=xmax*aL;}
     if(y[istep]<ymin){ymin=ymin*aL;}
     if(y[istep]>ymax){ymax=ymax*aL;}
     panel.setPreferredMinMax(xmin,xmax,ymin,ymax);
     frame.render();
    }

   public void clear () {
     myControl.clearMessages();
     trail.clear();
     panel.clear();
     frame.render();
     plot0.clearData();
     plot1.clearData();
   }

   public void resetAnimation () {
     clear();
     myControl.println ("Solves the variable mass rocket equation in the");
     myControl.println ("precense of a gravitational force. The equations");
     myControl.println ("are: dx/dt=vx, dy/dt=vy, dm/dt=-alpha(t).");
     myControl.println ("dvx/dt=ux*alpha(t)/m(t)-M*x/(x^2+y^2)^(3/2)");
     myControl.println ("dvx/dt=uy*alpha(t)/m(t)-M*y/(x^2+y^2)^(3/2),");
     myControl.println ("Refer to the text for further details.");
     myControl.println ("Decrease dt by increasing NPTS or decreasing tmax.");
     myControl.println ("dt=solver step size. Press 'Initialize' to proceed.");
     x     = null;
     y     = null;
     r     = null;
     vx    = null;
     vy    = null;
     t     = null;
     G=6.67e-11;                    //universal gravitational constant (Nm^2/kg^2)
     R=6.37e6;                      //unit of distance (m) - massive body radius
     m0=5.98e24;                    //unit of mass (kg) - could be the massive body
     M=5.98e24/m0;                  //massive body mass in units of m0
     tau=Math.sqrt(R*R*R/(G*m0));   //unit of time in seconds
     v0=R/tau; F0=m0*R/tau/tau;     //unit of speed and unit of force
     mi0=2.8e6;                     //initial payload+fuel mass in units of m0
     ff=0.96;                       //fuel fraction
     mi=mi0/m0;                     //initial payload+fuel mass in units of m0
     Thrust=1.5*(G*mi*M/R/R)*m0*m0; //Let the Thrust be # times initial mass weight
     u=0.35;                        //gas exhaust velocity in units of v0
     an=55;                         //angle of burn in degrees (determines launch angle)
     mp=mi*ff;                      //fuel mass
     mf=mi-mp;                      //final mass is payload mass (after fuel burnout)
     alpha=(Thrust/F0/u);           //alpha in units of m0/tau
     tfmax=(mi-mf)/alpha;           //fuel burnout time in units of tau
     t0=0.;
     tmax=Math.round(50*tfmax);     //simulation run time
     x0=0;y0=1;vx0=0;vy0=0;         //initial positions (R), speeds (R/tau)
     NPTS=300;
     delayTime=50;  //time in between animation steps (Abstract Animation)
     myControl.setValue("x0 (R) units",x0);
     myControl.setValue("y0 (R) units",y0);
     myControl.setValue("vx0 (R/tau) units",vx0);
     myControl.setValue("vy0 (R/tau) units",vy0);
     myControl.setValue("Earth mass(m0 units)",M);
     myControl.setValue("payload+fuel (kg)",df.format(mi0));
     myControl.setValue("fuel mass fraction",ff);
     myControl.setValue("Thrust(Newtons)",df.format(Thrust));
     myControl.setValue("exhaust speed (v0 units)",u);
     myControl.setValue("angle of burn (degrees)",an);
     myControl.setValue("tmax",tmax);
     myControl.setValue("t0",t0);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     dt=Math.round(1000.*(tmax-t0)/(NPTS-1))/1000.;
     myControl.setValue("dt((tmax-t0)/(NPTS-1))=",dt);
   }

   public void initializeAnimation() {
     clear();
     //initial conditions
     istep=0;  //step counter
     istep_max=NPTS; //calculate() can change this to animate above ground orbit points only
     iflag=0;        //changes value in calculate() if istep_max also changes value
     clear();
     myControl.clearMessages();
     x0=myControl.getDouble    ("x0 (R) units");
     y0=myControl.getDouble    ("y0 (R) units");
     vx0=myControl.getDouble   ("vx0 (R/tau) units");
     vy0=myControl.getDouble   ("vy0 (R/tau) units");
     M=myControl.getDouble     ("Earth mass(m0 units)");
     mi0=myControl.getDouble   ("payload+fuel (kg)");
     ff=myControl.getDouble    ("fuel mass fraction");
     Thrust=myControl.getDouble("Thrust(Newtons)");
     u=myControl.getDouble     ("exhaust speed (v0 units)");
     an=myControl.getDouble    ("angle of burn (degrees)");
     tmax=myControl.getDouble  ("tmax");
     t0=myControl.getDouble    ("t0");
     NPTS=myControl.getInt     ("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     dt=Math.round(1000.*(tmax-t0)/(NPTS-1))/1000.;
     myControl.setValue        ("dt((tmax-t0)/(NPTS-1))=",dt);
     Thrust=Thrust/F0;       //Thrust in units of force
     mi=mi0/m0;              //initial payload+fuel mass in units of m0
     mp=mi*ff;               //fuel mass
     th=an*cf;               //angle of burn in rad determines launch angle
     ux=u*Math.cos(th);      //exhaust velocity components in units of v0
     uy=u*Math.sin(th);
     alpha=(Thrust/u);       //alpha in units of m0/tau
     mf=mi-mp;               //final mass is payload mass (after fuel burnout)
     tfmax=(mi-mf)/alpha;    //fuel burnout time in units of tau
     x    = new double[NPTS+1];
     y    = new double[NPTS+1];
     vx   = new double[NPTS+1];
     vy   = new double[NPTS+1];
     v    = new double[NPTS+1];
     t    = new double[NPTS+1];
     m    = new double[NPTS+1];
     r    = new double[NPTS+1];
     v    = new double[NPTS+1];
     //initial conditions
     //state[0..6]:x,vx,y,vy,m,t
     state[0]=x0;
     state[1]=vx0;
     state[2]=y0;
     state[3]=vy0;
     state[4]=mi;
     state[5]=t0;
     x[istep]=state[0];
     vx[istep]=state[1];
     y[istep]=state[2];
     vy[istep]=state[3];
     m[istep]=state[4];
     t[istep]=state[5];
     r[istep]=Math.sqrt(x[istep]*x[istep]+y[istep]*y[istep]);
     v[istep]=Math.sqrt(vx[istep]*vx[istep]+vy[istep]*vy[istep]);
     aL=1+0.1;                             //used in checkPanel
     mL=2.0*Math.sqrt(x0*x0+y0*y0);        //initial windows view range
     xmin=-mL; xmax=mL; ymin=-mL; ymax=mL; //these can be modified by checkPanel()
     vmin=0; vmax=0; rmin=0; rmax=0;       //other values are determined by calculate()
     odesolver.initialize(dt);
     calculate();                          //get the whole orbit
     plot0.setPreferredMinMax(t0,t[istep_max],rmin,rmax);
     plot1.setPreferredMinMax(t0,t[istep_max],vmin,vmax);//vmin, vmax obtained by calculate()
     panel.setPreferredMinMax(xmin,xmax,ymin,ymax);
     panel.setSquareAspect(true);
     panel.setMessage("x="+df.format(x[istep]),0);
     panel.setMessage("y="+df.format(y[istep]),1);
     panel.setMessage("t="+Util.f2(t[istep]),3);
     trail.color=java.awt.Color.red;
     point.color=java.awt.Color.blue; //the rockect
     point.pixRadius=2;
     earth=DrawableShape.createCircle(0,0,2.0);//circle at [x,y]=(0,0), diameter=2
     earth.color=java.awt.Color.green;
     earth.edgeColor=java.awt.Color.darkGray;
     panel.addDrawable(earth);
     myControl.println ("Press 'start' to get the animation going.");
     myControl.println("Lower plots: black, blue, red: (left) x,y,r; (right) vx,vy,v");
     myControl.println("Units: v0="+nf.format(v0)+" m/s, R="+df.format(R)+" m");
     myControl.println("tau="+nf.format(tau)+" s, m0="+df.format(m0)+" kg");
     myControl.println("exhaust speed="+nf.format(u*v0)+" m/s, alpha="+
                       nf.format(alpha*m0/tau)+" kg/s");
     myControl.print("fuel burnout time="+nf.format(tfmax)+" tau, flight time="+
                       nf.format(t[istep_max])+" tau"+",\nfinal mass="+df.format(mf*m0)+
                       " kg, ");
     if(iflag==1){
       myControl.println("Range="+df.format(Range)+" m"+"\nRocket (blue), Earth (green)");
      }else{
       myControl.println("Rocket does not land. \nRocket (blue), Earth (green)");
      }
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new rocketApp();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(480, 5);
     myControl.setSize(315,555);
     myControl.setDividerLocation(315);
     app.setControl(myControl);
   }
}