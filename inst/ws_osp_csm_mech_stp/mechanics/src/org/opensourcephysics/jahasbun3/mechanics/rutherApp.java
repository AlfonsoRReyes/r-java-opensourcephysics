/*
@Author J E Hasbun 2007.
Simulates the Rutherford scattering alpha particle path numerically, given the
energy and the impact parameter.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
//import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;
import java.text.*;

  public class rutherApp extends AbstractAnimation {
    PlottingPanel panel = new PlottingPanel
        ("x(a_{b})","y(a_{b})","Apha-particle projectile, gold target target");
    DrawingFrame  frame = new DrawingFrame(panel);
    Trail trail=new Trail();
    Trail traila=new Trail();
    Trail trailb=new Trail();
    Circle point=new Circle();
    Circle point_t=new Circle();
    ODE ode = new ruther_der();
    ODESolver odesolver = new RK4(ode);
    //ODESolver odesolver = new RK45MultiStep(ode);
    private Control myControl;
    double[] state= new double [5];
    double [] x1,y1,r,xa,ya1,ya2;
    double [] t;
    double rmin,ye,za,zt,t0,tmax,dt;
    double m,vb,ma,Ene,v0,b,K,sK,sV,sa,dxa,THETA;
    double mLx,pLx,mLy,pLy;
    double x10,y10,v1x0,v1y0;
    double pi;
    int NPTS,istep,istart_plot,iend_plot;
    NumberFormat nf = NumberFormat.getInstance();
    DecimalFormat df= new DecimalFormat("+0.000E00");

  public rutherApp (){
    frame.setLocation(5,5);
    frame.setSize(400,400);
    frame.setTitle("Apha-particle projectile, gold target target");
    nf.setMaximumFractionDigits(4);
    pi=Math.PI;
    trail=new Trail();
    traila=new Trail();
    trailb=new Trail();
    trail.setDashedStroke(2,4); //trail: thickness=2, dash length=4
    traila.setDashedStroke(2,4);
    trailb.setDashedStroke(2,4);
  }

   protected void doStep() {
     if(istep<iend_plot){
     //create a dashed trail
     trail.addPoint(x1[istep], y1[istep]);   //draw line every even point
     point.setXY(x1[istep],y1[istep]);
     //draw
     panel.addDrawable(trail);
     panel.addDrawable(point);
     frame.render();
     //myControl.println("istep,x1,y1="+istep+", "+x1[istep]+", "+y1[istep]);
     //myControl.println("t,x1,y1="+t[istep]+", "+x1[istep]+", "+y1[istep]);
     istep++;
     }else{
       asymptotes();  //finds the asymptotes and scattering angle THETA
       panel.addDrawable(traila);//trails a,b come from asymptote()
       panel.addDrawable(trailb);
       frame.render();
       myControl.println("rmin=\t"+nf.format(rmin)+" a_b");//rmin is foun in calculate();
       myControl.println("THETA="+nf.format(THETA)+" rad, or "+nf.format(THETA*180/pi)+
                         " degrees");
       stopAnimation();
     }
   }

   public void calculate () {
     //orbits calculation
     //state[0,1,2,3,4]=x1,v1x,y1,v1y,t
     int flag1=0, flag2=0;
     rmin=Double.POSITIVE_INFINITY;
     for (istep=1; istep<NPTS; istep++){
       t[istep] = ode.getState()[4];
       odesolver.step();
       x1[istep] = ode.getState()[0];
       y1[istep] = ode.getState()[2];
       r[istep] = Math.sqrt(x1[istep]*x1[istep]+y1[istep]*y1[istep]);
       if(r[istep]<rmin){rmin=r[istep];}
       //istart_plot & iend_plot will be used in doStep() to plot within the viewing range
       if(x1[istep]<pLx*(1+0.1) && flag1==0){istart_plot=istep; flag1=1;}
       if(y1[istep]>pLy*(1+0.1) && flag2==0){iend_plot=istep; flag2=1;}
       //myControl.println("istep,t,r="+istep+", "+t[istep]+", "+r[istep]);
     }
       //myControl.println("istart_plot, iend_plot="+istart_plot+", "+iend_plot);
   }

   class ruther_der implements ODE {
     public double[] getState() {
       return state;
     }

     public void getRate(double[] state, double[] rate) {
       //state[0,1,2,3,4]=x1,v1x,y1,v1y,t
       //rates are the derivatives of the state, for example
       //rate[0]=dr/dt=v->state(1), and dv/dt=K*state[0]/wr/m;
       //etc
       double wr=Math.pow(state[0]*state[0]+state[2]*state[2],3.0/2.0);
       rate[0]=state[1];
       rate[1]=K*state[0]/wr/m;
       rate[2]=state[3];
       rate[3]=K*state[2]/wr/m;
       rate[4]=1; //dt/dt
     }
   }

   public void clear () {
       trail.clear();
       traila.clear();
       trailb.clear();
       panel.clear();
       frame.render();
   }

   public void resetAnimation () {
     //initial conditions
     istep=-1;  //so that istep+1 in doStep =0 the first time
     clear();
     myControl.clearMessages();
     myControl.println ("Press 'Initialize' then press 'start.' This Simulates");
     myControl.println ("the Rutherford scattering alpha particle path numerically,");
     myControl.println ("given the energy and the impact parameter. The equations");
     myControl.println ("are: dx/dt=v, dvx/dt=K*x/r^3/m with r=(x^2+y^2)^(1/2);");
     myControl.println ("also dy/dt=vx, dvy/dt=K*y/r^3/m. The scattering angle");
     myControl.println ("THETA is obtained from the projectile's entrance and exit");
     myControl.println ("trajectories. Other quantities (see text): K=za*zt*k*e^2,");
     myControl.println ("v0=sqrt(2*E/m_alpha), e=-sign(K)*sqrt(1+m^2*v0^4*b^2/K^2)");
     myControl.println ("and rmin=minimum value of (r).");
     x1    = null;
     y1    = null;
     r     = null;
     t     = null;
     xa    = null;
     ya1   = null;
     ya2   = null;
     za=2; zt=79; //za=projectile, zt=target charges (2=alpha, 70=gold)
     m=1;         //projectile mass in units of alpha particle mass
     vb=0.01965;  //velocity unit (a_b/tau_b) in units of light speed
     ma=3730.0e6; //alpha particle mass energy in eV
     Ene=5.0e6;   //initial projectile energy in eV
     b=20;        //impact parameter in units of ab=1.e-15 meters
     NPTS=600;
     tmax=2500.0; t0=0.0; //time range in units of tau_b=1.7e-22 sec
     delayTime=50;  //time in between animation steps (Abstract Animation)
     mLx=-50; pLx=100; mLy=-50; pLy=100; //windows view range
     myControl.setValue("projectile_energy(eV)",df.format(Ene));
     myControl.setValue("impact_parameter(a_b)",b);
     myControl.setValue("projectile_mass(m_alpha)",m);
     myControl.setValue("projectile_z",za);
     myControl.setValue("target_z",zt);
     myControl.setValue("tmax",tmax);
     myControl.setValue("t0",t0);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     myControl.setValue("max_view_-x",mLx);
     myControl.setValue("max_view_+x",pLx);
     myControl.setValue("max_view_-y",mLy);
     myControl.setValue("max_view_+y",pLy);
     dt=(tmax-t0)/(NPTS-1);
     myControl.setValue("dt((tmax-t0)/(NPTS-1))=",nf.format(dt));
   }

   public void initializeAnimation() {
     //initial conditions
     istep=0;  //so that istep+1 in doStep =0 the first time
     clear();
     myControl.clearMessages();
     Ene=myControl.getDouble("projectile_energy(eV)");
     b=myControl.getDouble("impact_parameter(a_b)");
     m=myControl.getDouble("projectile_mass(m_alpha)");
     za=myControl.getDouble("projectile_z");
     zt=myControl.getDouble("target_z");
     tmax=myControl.getDouble("tmax");
     t0=myControl.getDouble("t0");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     mLx=myControl.getDouble("max_view_-x");
     pLx=myControl.getDouble("max_view_+x");
     mLy=myControl.getDouble("max_view_-y");
     pLy=myControl.getDouble("max_view_+y");
     dt=(tmax-t0)/(NPTS-1);
     myControl.setValue("dt((tmax-t0)/(NPTS-1))=",dt);
     istart_plot=istep;//lower limit points to plot - changed in calculate()
     iend_plot=NPTS;   //upper limit points to plot - changed in calculate()
     v0=Math.sqrt(2*Ene/m/ma)/vb; //initial projectile speed in units of vb
     K=za*zt;          //dimensionless force constant
     sK=K/Math.abs(K); //sign of K
     ye=-sK*Math.sqrt(1+m*m*Math.pow(v0,4)*b*b/K/K); //eccentricity
     rmin=-m*v0*v0*b*b/K/(1+ye); //closest approach distance
     x10=60*100;  //initial x position - far away (a_b=1 fermi unit)
     y10=b;       //initial y position is the impact parameter
     v1y0=0.0;    //zero y-speed - projectile coming in horizontally
     //x10 by energy conservation is as follows
     if(x10!=0){
       sV=x10/Math.abs(x10); //get the sign of x10
       v1x0=-sV*Math.sqrt(v0*v0-K/x10/m/2);
     }else{v1x0=0;}
     x1    = new double[NPTS+1];
     y1    = new double[NPTS+1];
     r     = new double[NPTS+1];
     t     = new double[NPTS+1];
     xa    = new double[NPTS+1];
     ya1   = new double[NPTS+1];
     ya2   = new double[NPTS+1];
     //initial conditions
     state[0]=x10;
     state[1]=v1x0;
     state[2]=y10;
     state[3]=v1y0;
     state[4]=t0;
     x1[istep]=state[0];
     y1[istep]=state[2];
     r[istep]=Math.sqrt(x1[istep]*x1[istep]+y1[istep]*y1[istep]);
     odesolver.initialize(dt);
     calculate();       //get the whole orbit
     istep=istart_plot; //so that in doStep() we only plot in the viewing range
     panel.setPreferredMinMax(mLx,pLx,mLy,pLy);
     trail.color=java.awt.Color.red;
     traila.color=java.awt.Color.blue;
     trailb.color=java.awt.Color.blue;
     point=new Circle();           //used for the projectile
     point.color=java.awt.Color.red;
     point.pixRadius=3;
     point_t=new Circle();         //the target
     point_t.color=java.awt.Color.black;
     point_t.pixRadius=3;
     point_t.setXY(0.0,0.0);
     panel.addDrawable(point_t);
     myControl.println ("Press 'start' to get the animation going.");
     myControl.println("The simulated motion of the + projectile (red) interacting ");
     myControl.println("(electrically) with a + heavy target (black). Blue lines are");
     myControl.println("asymptotes. Units follow.");
     myControl.println("speed="+nf.format(vb)+"c, m_alpha="+
                       df.format(ma)+"eV;\na_b=1.e-15m, time=1.695e-22s");
     myControl.println("eccentricity e=\t"+nf.format(ye));
   }

   public void asymptotes(){
     //finds the asymptotes abd scattering angle
     sa=(y1[NPTS]-y1[NPTS-1])/(x1[NPTS]-x1[NPTS-1]);//outgoing asymptote slope - use last point
     dxa=(pLx-mLx)/NPTS;
     for(int i=0; i<=NPTS; i++){
       xa[i]=mLx+i*dxa;   //the variable
       ya1[i]=y1[NPTS]+sa*(xa[i]-x1[NPTS]);//outgoing asymtote
       ya2[i]=y1[0];                       //incoming asymptote
       //create a dashed trail
       traila.addPoint(xa[i], ya1[i]);     //asymtotes trails
       trailb.addPoint(xa[i], ya2[i]);     //lines on even points
     }
     if(sa >= 0){
       THETA=pi-Math.atan(sa);                     //first quadrant case
     }else{
       THETA=pi-(Math.atan(Math.abs(1./sa))+pi/2); //2nd quadrant case
     }
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new rutherApp();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(410, 5);
     myControl.setSize(360,555);
     myControl.setDividerLocation(300);
     app.setControl(myControl);
   }
}