/*
@Author J E Hasbun 2007.
Uses the conic section formula with rmin and eccentricity to simulate
Rutherford scattering given the projectile energy and the impact parameter.
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

  public class conic3App extends AbstractAnimation {
  PlottingPanel panel = new PlottingPanel
  ("x(a_{b})","y(a_{b})","Apha-particle projectile, gold target target");
  DrawingFrame  frame = new DrawingFrame(panel);
  PlotFrame plot= new PlotFrame("theta","cos(theta), -1/e","Asymptotes' angles");
  Trail trail=new Trail();
  Circle point=new Circle();
  Circle point_t=new Circle();
  Circle pointa=new Circle();
  Circle pointb=new Circle();
  private Control myControl;
  double [] x,y;
  double [] ths;
  double rmin,ye,thmin,thmax,r,THETA,dth,za,zt,th;
  double m,vb,ma,Ene,v0,b,K,sK;
  double mLx,pLx,mLy,pLy;
  double pi;
  int NPTS,istep;
  int maxiter=25;
  double tolerance=5e-3;
  Bisnewt bisnewt=new Bisnewt();
  Function fc;
  NumberFormat nf = NumberFormat.getInstance();

  public conic3App (){
    plot.setLocation(5,5);
    plot.setSize(400,240);
    plot.setConnected(true);
    plot.setName("Asymptotes' angles");
    frame.setLocation(5,245);
    frame.setSize(400,350);
    frame.setTitle("Apha-particle projectile, gold target target");
    nf.setMaximumFractionDigits(4);
    pi=Math.PI;
    trail=new Trail();
    trail.setDashedStroke(2,4); //trail: thickness=2, dash length=4
  }

   protected void doStep() {
     if(istep<NPTS-1){
       //create a dashed trail
       trail.addPoint(x[istep], y[istep]);   //draw line every even point
       point.setXY(x[istep],y[istep]);

       istep++;
       calculate();
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
     th=thmin+istep*dth;    //angle variable
     r=rmin*(1+ye)/(1+ye*Math.cos(th));
     ths[istep]=th-thmin;   //rotate orbital path counter-clockwise by 'thmin'
                            //to align asymptote with +x axis
     x[istep]=r*Math.cos(ths[istep]);
     y[istep]=r*Math.sin(ths[istep]);
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
     istep=0;  //doStep index counter
     clear();
     myControl.clearMessages();
     myControl.println ("Initialize: initializes; start: starts the simulation.");
     myControl.println ("Uses the conic section formula with rmin and eccentricity");
     myControl.println ("to simulate Rutherford scattering given the projectile");
     myControl.println ("energy and the impact parameter.");
     myControl.println ("We use the conic section r=rmin*(1+e)/(1+e*cos(theta), with");
     myControl.println ("y=r*sin(theta), x=r*cos(theta). The path is simulated");
     myControl.println ("on the interval [thmin,thmax] which should lie within");
     myControl.println ("the asymptotes angles obtained by solving 0=1+e*cos(x).");
     myControl.println ("The scattering angle THETA is the angle between the");
     myControl.println ("projectile's entrance and exit trajectories.");
     myControl.println ("Other quantities are (see text):v0=sqrt(2*E/m_alpha),");
     myControl.println ("K=za*zt*k*e^2, e=-sign(K)*sqrt(1+m^2*v0^4*b^2/K^2), and");
     myControl.println ("rmin=-m*v0^2*b^2/K/(1+e).");
     x    = null;
     y    = null;
     ths  = null;
     za=2; zt=79; //za=projectile, zt=target charges (2=alpha, 70=gold)
     m=1;         //projectile mass in units of alpha particle mass
     vb=0.01965;  //velocity unit (a_b/tau_b) in units of light speed
     ma=3730.0e6; //alpha particle mass energy in eV
     Ene=5.0e6;   //initial projectile energy in eV
     b=20;        //impact parameter in units of ab=1.e-15 meters
     NPTS=200;
     thmin=-0.7; thmax=0.7; //for now to be safe, later changed by asymptotes()
     delayTime=25;  //time in between animation steps (Abstract Animation)
     mLx=-50; pLx=100; mLy=-50; pLy=100; //windows view range
     myControl.setValue("projectile_energy(eV)",Ene);
     myControl.setValue("impact_parameter(a_b)",b);
     myControl.setValue("projectile_mass(m_alpha)",m);
     myControl.setValue("projectile_z",za);
     myControl.setValue("target_z",zt);
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
     istep=0;  //doStep counter
     clear();
     myControl.clearMessages();
     Ene=myControl.getDouble("projectile_energy(eV)");
     b=myControl.getDouble("impact_parameter(a_b)");
     m=myControl.getDouble("projectile_mass(m_alpha)");
     za=myControl.getDouble("projectile_z");
     zt=myControl.getDouble("target_z");
     NPTS=myControl.getInt("NPTS");
     delayTime=myControl.getInt("delayTime(ms)");
     mLx=myControl.getDouble("max_view_-x");
     pLx=myControl.getDouble("max_view_+x");
     mLy=myControl.getDouble("max_view_-y");
     pLy=myControl.getDouble("max_view_+y");
     v0=Math.sqrt(2*Ene/m/ma)/vb; //initial projectile speed in units of vb
     K=za*zt;          //dimensionless force constant
     sK=K/Math.abs(K); //sign of K
     ye=-sK*Math.sqrt(1+m*m*Math.pow(v0,4)*b*b/K/K); //eccentricity
     rmin=-m*v0*v0*b*b/K/(1+ye); //closest approach distance
     asymptotes();  //finds the assymtotes angles thmin, thmax for a given ye
     THETA=pi-(thmax-thmin);     //scattering angle
     dth=(thmax-thmin)/(NPTS-1);
     myControl.setValue("dth((thmax-thmin)/(NPTS-1))=",dth);
     x    = new double[NPTS];
     y    = new double[NPTS];
     ths  = new double[NPTS];
     //initial conditions (istep=0)
     th=thmin;
     r=rmin*(1+ye)/(1+ye*Math.cos(th));
     ths[istep]=th-thmin;   //rotate orbital path counter-clockwise by 'thmin'
                            //to align asymptote with +x axis
     x[istep]=r*Math.cos(ths[istep]);
     y[istep]=r*Math.sin(ths[istep]);
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
     myControl.println("Upper graph: the asymptotes angles as solutions of");
     myControl.println("0=1+e*cos(x). The black dots are the roots from which");
     myControl.println("THETA is found and given below. Press 'start' for the simulation");
     myControl.println("to get going. Lower graph: The simulated motion of the + ");
     myControl.println("projectile (red) interacting (electrically) with a + heavy ");
     myControl.println("target (black). Units follow");
     myControl.println("speed="+nf.format(vb)+"c, m_alpha="+
                       nf.format(ma)+"eV;\na_b=1.e-15m, time=1.695e-22s");
     myControl.println("rmin=\t"+nf.format(rmin)+" a_b");
     myControl.println("eccentricity e=\t"+nf.format(ye));
     myControl.println("Asymptote angle_1=\t"+nf.format(thmin)+" rad");
     myControl.println("Asymptote angle_2=\t"+nf.format(thmax)+" rad");
     myControl.println("THETA=\t"+nf.format(THETA)+" rad"+
                       ", or "+nf.format(THETA*180/pi)+" degrees");
   }

   public void asymptotes(){
     //finds the asymptote angles by the newton raphson method
     fc=new fcos();
     thmin=bisnewt.Bisnewt(-pi/2,1.e-5, maxiter, tolerance,fc);
     thmax=bisnewt.Bisnewt(1.e-5,pi/2, maxiter, tolerance,fc);
     thmin=thmin+15.e-3; //displaced a little from exact
     thmax=thmax-15.e-3; //displaced a little from exact
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
     pointa.setXY(thmin,-1/ye);
     pointb.setXY(thmax,-1/ye);
     plot.addDrawable(pointa);
     plot.addDrawable(pointb);
     plot.render();
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new conic3App();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(410, 5);
     myControl.setSize(360,590);
     myControl.setDividerLocation(275);
     app.setControl(myControl);
   }
}