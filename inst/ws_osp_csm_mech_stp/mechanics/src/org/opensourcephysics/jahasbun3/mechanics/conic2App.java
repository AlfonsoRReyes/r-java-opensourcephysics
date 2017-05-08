/*
@Author J E Hasbun 2007.
Obtains the hyperbolic projectile orbit of a positive projectile on a positive
and fixed target. The conic section formula with negative eccentricity is used.
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
import java.text.NumberFormat;

  public class conic2App extends AbstractAnimation {
  PlottingPanel panel = new PlottingPanel
                          ("x","y","Projectile orbit, fixed target");
  DrawingFrame  frame = new DrawingFrame(panel);
  PlotFrame plot= new PlotFrame("theta","cos(theta), -1/e","Asymptotes' angles");
  Trail trail=new Trail();
  Circle point=new Circle();
  Circle point_t=new Circle();
  Circle pointa=new Circle();
  Circle pointb=new Circle();
  private Control myControl;
  double [] x,y;
  double [] th;
  double rmin,ye,thmin,thmax,r,THETA,dth,thmin_s,thmax_s;
  double mLx,pLx,mLy,pLy;
  double pi;
  int NPTS,istep;
  int maxiter=25;
  double tolerance=5e-3;
  Bisnewt bisnewt=new Bisnewt();
  Function fc;
  NumberFormat nf = NumberFormat.getInstance();

  public conic2App (){
    plot.setLocation(5,5);
    plot.setSize(400,295);
    plot.setConnected(true);
    plot.setName("Asymptotes' angles");
    frame.setLocation(5,300);
    frame.setSize(400,295);
    frame.setTitle("Projectile orbit, fixed target");
    nf.setMaximumFractionDigits(4);
    pi=Math.PI;
    trail=new Trail();
  }

   protected void doStep() {
     if(istep+1<NPTS){
       istep++;
       calculate();
       //create a dashed trail as follows
       if(istep%2!=0){
         trail.moveToPoint(x[istep], y[istep]);//skip line every odd point
       }
       if(istep%2==0){
         trail.addPoint(x[istep], y[istep]);   //draw line every even point
       }
       point.setXY(x[istep],y[istep]);
       //draw
       panel.addDrawable(trail);
       panel.addDrawable(point);
       frame.render();
       //myControl.println("istep,x,y="+istep+", "+x[istep]+", "+y[istep]);
     }else{
       stopAnimation();
     }
   }

   public void calculate () {
     //orbit calculation
     th[istep]=thmin+istep*dth;             //angle variable
     r=rmin*(1+ye)/(1+ye*Math.cos(th[istep]));
     x[istep]=r*Math.cos(th[istep]);
     y[istep]=r*Math.sin(th[istep]);
     //myControl.println("istep,r="+istep+", "+th[istep]+", "+r);
   }

   public void clear () {
       trail.clear();
       panel.clear();
       frame.render();
       plot.clearDrawables();
       plot.render();
   }

   public void resetAnimation () {
     //initial conditions
     istep=-1;  //so that istep+1 in doStep =0 the first time
     clear();
     myControl.clearMessages();
     myControl.println ("Initialize: initializes; start: starts the simulation.");
     myControl.println ("Obtains the hyperbolic projectile(+) orbit with negative");
     myControl.println ("eccentricity incident on a positive and fixed target.");
     myControl.println ("We use the conic section r=rmin*(1+e)/(1+e*cos(theta), with");
     myControl.println ("y=r*sin(theta), x=r*cos(theta). The path is simulated");
     myControl.println ("on the interval [thmin,thmax] which should lie within");
     myControl.println ("the asymptotes angles obtained by solving 0=1+e*cos(x).");
     myControl.println ("The scattering angle THETA is the angle between the");
     myControl.println ("projectile's entrance and exit trajectories.");
     x    = null;
     y    = null;
     th   = null;
     rmin=1.0;      //minimum distance
     ye=-1.2;       //the eccetricity
     NPTS=200;
     delayTime=25;  //time in between animation steps (Abstract Animation)
     mLx=-1; pLx=5; mLy=-5; pLy=5; //windows view range
     //seek thmin_s on [-pi/2,0] and thmax_s on [0,pi/2]
     asymptotes();  //finds the assymptotes angles thmin_s, thmax_s for a given ye
     thmin=thmin_s; thmax=thmax_s;
     myControl.setValue("eccentricity",ye);
     myControl.setValue("rmin",rmin);
     myControl.setValue("theta_min",thmin);
     myControl.setValue("theta_max",thmax);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("delayTime(ms)",delayTime);
     myControl.setValue("max_view_-x",mLx);
     myControl.setValue("max_view_+x",pLx);
     myControl.setValue("max_view_-y",mLy);
     myControl.setValue("max_view_+y",pLy);
     dth=(thmax-thmin)/(NPTS-1);
     myControl.setValue("dth((thmax-thmin)/(NPTS-1))=",dth);
   }

   public void initializeAnimation() {
     //initial conditions
     istep=-1;  //so that istep+1 in doStep =0 the first time
     clear();
     myControl.clearMessages();
     ye=myControl.getDouble("eccentricity");
     asymptotes();  //finds the assymtotes angles thmin_s, thmax_s for a given ye
     rmin=myControl.getDouble("rmin");
     NPTS=myControl.getInt("NPTS");
     thmin=myControl.getDouble("theta_min");
     thmax=myControl.getDouble("theta_max");
     if(thmin<thmin_s){
       thmin=thmin_s;
       myControl.setValue("theta_min",thmin);
     }
     if(thmax>thmax_s){
       thmax=thmax_s;
       myControl.setValue("theta_max",thmax);
     }
     delayTime=myControl.getInt("delayTime(ms)");
     mLx=myControl.getDouble("max_view_-x");
     pLx=myControl.getDouble("max_view_+x");
     mLy=myControl.getDouble("max_view_-y");
     pLy=myControl.getDouble("max_view_+y");
     dth=(thmax-thmin)/(NPTS-1);
     myControl.setValue("dth((thmax-thmin)/(NPTS-1))=",dth);
     x    = new double[NPTS];
     y    = new double[NPTS];
     th   = new double[NPTS];
     panel.setPreferredMinMax(mLx,pLx,mLy,pLy);
     trail.color=java.awt.Color.red;
     point=new Circle();           //used for the projectile
     point.color=java.awt.Color.red;
     point.pixRadius=3;
     pointa=new Circle();          //used to indicate the lower asymptote angle
     pointb=new Circle();          //used to indicate the upper asymptote angle
     pointa.color=java.awt.Color.black; pointb.color=java.awt.Color.black;
     pointa.pixRadius=3; pointb.pixRadius=3;
     point_t=new Circle();         //the target
     point_t.color=java.awt.Color.black;
     point_t.pixRadius=3;
     point_t.setXY(0.0,0.0);
     panel.addDrawable(point_t);
     calc_cos_func();              //plots the cos(theta) and -1/e functions
     THETA=pi-(thmax_s-thmin_s);   //scattering angle
     myControl.println("Upper graph: the asymptotes angles as solutions of");
     myControl.println("0=1+e*cos(x). The black dots are the root positions");
     myControl.println("given below. The scattering angle THETA is given below.");
     myControl.println("Press 'start' to get the simulation going.");
     myControl.println("Lower graph: The simulated motion of the + projectile (red)");
     myControl.println("interacting (electrically) with a + heavy target (black).");
     myControl.println("Asymptote angle_1=\t"+nf.format(thmin_s)+" rad");
     myControl.println("Asymptote angle_2=\t"+nf.format(thmax_s)+" rad");
     myControl.println("THETA=pi-(thmax_s-thmin_s)\t"+nf.format(THETA)+" rad");
   }

   public void asymptotes(){
     //finds the asymptote angles by the newton raphson method
     fc=new fcos();
     thmin_s=bisnewt.Bisnewt(-pi/2,1.e-5, maxiter, tolerance,fc);
     thmax_s=bisnewt.Bisnewt(1.e-5,pi/2, maxiter, tolerance,fc);
   }

   public class fcos implements Function {
     public double evaluate(double x) {
       //defines f(x) whose value vanishes at x=root of f(x)
       double f;
       f=1+ye*Math.cos(x);
       return f;
     }
   }

   public void calc_cos_func() {
     //plots the cos(theta) and -1/e functions
     int N=200;
     double [] xx=new double [N];
     double [] yy=new double [N];
     double [] yo=new double [N];
     double t0=-2*pi, tf=2*pi, ts;
     ts=(tf-t0)/(N-1);
     for (int i=0; i<N; i++){
       xx[i]=t0+i*ts;
       yy[i]=Math.cos(xx[i]);
       yo[i]=-1.0/ye;
     }
     plot.append(0,xx,yy);
     plot.setMarkerSize(0,0);
     plot.append(1,xx,yo);
     plot.setMarkerSize(1,0);
     pointa.setXY(thmin_s,-1/ye);
     pointb.setXY(thmax_s,-1/ye);
     plot.addDrawable(pointa);
     plot.addDrawable(pointb);
     plot.render();
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new conic2App();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(410, 5);
     myControl.setSize(360,590);
     myControl.setDividerLocation(250);
     app.setControl(myControl);
   }
}