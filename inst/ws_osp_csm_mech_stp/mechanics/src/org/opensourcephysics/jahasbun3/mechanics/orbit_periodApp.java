/*
@Author J E Hasbun 2007.
This finds the time it takes to go from rmin to rmax in an orbit
due to a force of the form F(r)=-a*r^p.
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
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.awt.Color;
import java.text.NumberFormat;

public class orbit_periodApp implements Calculation {
  PlotFrame plot0= new PlotFrame("r","Integrand",
                                 "Central Force=a*r^p - period integrand");
  PlottingPanel polarPanel = new PlottingPanel
                          ("Theta","r","Polar Plot - Central Force=a*r^p");
  DrawingFrame  polarFrame = new DrawingFrame(polarPanel);
  PolarType1 axes = new PolarType1(polarPanel);
  Trail trail=new Trail();
  private Control myControl;
  double r[], v[], th[], t[];
  double[] state=new double [4];
  double r0,v0,th0,a,p,m,L,t0,tmax,dt,ra,vr,E;
  double rmin,rmax;
  int NPTS;
  NumberFormat nf = NumberFormat.getInstance();
  Dataset dataset = new Dataset();

  public orbit_periodApp (){
    polarFrame.setLocation(5,5);
    polarFrame.setSize(400,295);
    polarFrame.setTitle("Polar Plot - r(theta)");
    plot0.setConnected(true);//set all connected
    plot0.setLocation(5,300);
    plot0.setSize(400,295);
    nf.setMaximumFractionDigits(3);
  }

   public void calculate () {
     clear();
     double xmin,xmax,ymin,ymax, af=(1+0.1);
     ODE ode = new central_der();
     ODESolver odesolver = new RK4(ode);
     myControl.clearMessages();
     m=myControl.getDouble("m");
     p=myControl.getDouble("p");
     a=myControl.getDouble("a");
     r0=myControl.getDouble("r0");
     v0=myControl.getDouble("v0");
     th0=myControl.getDouble("th0");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     r     = new double[NPTS];
     v     = new double[NPTS];
     th    = new double[NPTS];
     t     = new double[NPTS];
     double [] xy_coord = new double[2];
     //initial conditions
     L=m*v0*r0; //angular momentum
     state[0]=r0;
     state[1]=v0;
     state[2]=th0;
     state[3]=t0;
     //differential equation solutions for the r,v,theta positions versus time
     odesolver.initialize(dt); // step size
     xmin=0; xmax=0; ymin=0; ymax=0;
     for (int i=0; i< NPTS;i++){
       t[i]=ode.getState()[3];
       odesolver.step();
       r[i]=ode.getState()[0];
       v[i]=ode.getState()[1];
       th[i]=ode.getState()[2];
       xy_coord[0]=r[i]*Math.cos(th[i]);
       xy_coord[1]=r[i]*Math.sin(th[i]);
       trail.addPoint(xy_coord[0],xy_coord[1]);
       if(xy_coord[0]<xmin){xmin=xy_coord[0];}
       if(xy_coord[0]>xmax){xmax=xy_coord[0];}
       if(xy_coord[1]<ymin){ymin=xy_coord[1];}
       if(xy_coord[1]<ymax){ymin=xy_coord[1];}
     }
     //ra is needed to estimate vr in the next line
     ra=ArrayLib.min(r);
     vr=v0*Math.sqrt((r0/ra)*(r0/ra)-1.0);//the radial velocity component
     if (p==-1){
       E=0.5*m*vr*vr+L*L/m/r0/r0/2+a*Math.log(r0);
     }
     else{
       E=0.5*m*vr*vr+L*L/m/r0/r0/2+a*Math.pow(r0,(p+1))/(p+1);
     }
    //search for the lower and upper limits where the square
    //of the orbit integrand is positive
    double st=0.00015, rlim=10*r0;
    double x=r0;
    int ilim=Math.round((float)(rlim/st));
    for (int i=0; i<ilim; i++){
      if(Ix(x)<=0){rmin=x+st; break;}
      x-=st;
    }
    x=r0;
    for (int i=0; i<ilim; i++){
      if(Ix(x)<=0){rmax=x-st; break;}
      x+=st;
    }
     int N=5001;
     if(N%2==0){N=N+1;}//make sure N is odd
     double [] r2=new double[N];
     double [] f=new double[N]; //double [] f2=new double[N];
     double a=rmin, b=rmax, h=(b-a)/(N-1);
     //the period is 4 * the integral of (Math.sqrt(Ix)) on {rmin,rmax}
     for (int i=0; i< N;i++){
       r2[i]=a+i*h;
       f[i]=Math.sqrt(Ix(r2[i]));
       //myControl.println("i="+i+", r2="+r2[i]+", f="+f[i]);
     }
     double period=4.*Simp.Simp(f,h);
     myControl.println("Upper plot is the body's orbit, the lower plot is");
     myControl.println("the integrand value plotted in the range [rmin,rmax].");
     myControl.println("The period is obtained by multiplying the shaded area");
     myControl.println("under this curve by 4.");
     myControl.println("(rmin,rmax)=("+nf.format(rmin)+","+nf.format(rmax)+
                       "), period="+nf.format(period));
     polarPanel.setPreferredMinMax(xmin*af,xmax*af,ymin*af,ymax*af);
     //polarPanel.setSquareAspect(true);
     trail.color=java.awt.Color.red;
     polarPanel.addDrawable(trail);
     axes.setDeltaR(Math.sqrt(xmax*xmax+ymax*ymax)/5);
     axes.setDeltaTheta(Math.PI/6);
     polarFrame.render();
     //The integrand plot
     plot0.setPreferredMinMax(rmin,rmax,0.0,0.6);
     plot0.setLineColor(0,Color.black);
     plot0.setMarkerSize(0,0);
     plot0.append(0,r2,f); //Integrand used in the period versus r
     dataset.append(r2,f); //used to shade the area under the integrand
     //dataset.setMarkerShape(dataset.AREA);
     //dataset.setMarkerColor(Color.blue);
     dataset.setMarkerShape(dataset.BAR);
     dataset.setMarkerColor(new java.awt.Color(113, 157, 170, 1));
     plot0.addDrawable(dataset);//shade the area under the integrand
     }

     public void clear () {
       dataset.clear();
       plot0.clearData();
       trail.clear();
       polarPanel.clear();
       polarFrame.render();
     }


   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("This finds the time it takes to go from rmin to rmax");
     myControl.println ("in an orbit due to a force of the form F(r)=-a*r^p.");
     myControl.println ("The orbit differential equations are solved: dr/dt=v");
     myControl.println ("dv/dt=-(a/m)*r^p+L^2/(m^2*r^3) and dtheta/dt=L/(m*r^2)");
     myControl.println ("with r0,v0, and th0 as init. conditions at t=t0. Here");
     myControl.println ("L=ang. momentum m*r*v0, & v0 the intial tangential speed.");
     myControl.println ("The period is obtained from 4*Integral of");
     myControl.println ("{1./sqrt(2*(E-V(r)-L^2/m/r^2/2)/m)} on {rmin,rmax}.");
     myControl.println ("where E=m*vr^2/2+m*vt^2/2+V(r), where dV/dr=-F(r), and");
     myControl.println ("vr=radial velocity, vt=tangential velocity=r*theta_dot.");
     myControl.println ("The angular momentum L=m*vt*r.");
     r    = null;
     v    = null;
     th    = null;
     t    = null;
     m=1;
     a=108; p=1;
     t0=0; tmax=1.5;
     r0=1.0;          //initial position
     v0=6.0;          //initial tangential speed
     th0=0.0;         //initial angle in radians
     NPTS=500;
     myControl.setValue("m",m);
     myControl.setValue("p",p);
     myControl.setValue("a",a);
     myControl.setValue("r0",r0);
     myControl.setValue("v0",v0);
     myControl.setValue("th0",th0);
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("NPTS",NPTS);
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
   }

   double Ix (double x) {
     //this function defines the square of the integrand for the
     //orbit period calculation
     double f=0;
       if(p==-1){
         f=1./(2.*(E-a*Math.log(x)-L*L/m/x/x/2)/m);
       } else {
         f=1./(2*(E-a*Math.pow(x,(p+1))/(p+1)-L*L/m/x/x/2)/m);
       }
     return f;
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   class central_der implements ODE {
     public double[] getState() {
     return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2,3]=r,v,theta,t
       //rates are the derivatives of the state, for example
       //rate[0]=dr/dt=v->state(1), and
       //dv/dt=-(a/m)*r^p+L^2/m^2/r^3, and dtheta/dt=L/m/r^2
       rate[0]=state[1];
       rate[1]=-(a/m)*Math.pow(state[0],p)+L*L/m/m/Math.pow(state[0],3);
       rate[2]=L/m/state[0]/state[0];
       rate[3]=1; //dt/dt=1
     }
   }

   public static void main(String[] args) {
     Calculation model = new orbit_periodApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(350,590);
     myControl.setDividerLocation(250);
     model.setControl(myControl);
   }
}