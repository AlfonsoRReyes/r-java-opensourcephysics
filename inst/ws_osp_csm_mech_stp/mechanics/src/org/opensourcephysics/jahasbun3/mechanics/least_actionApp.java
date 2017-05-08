/*
@Author J E Hasbun 2007.
Simulates Hamilton's Least Action principle for a particle under the action
of gravity.
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
import java.util.Random;


  public class least_actionApp extends AbstractAnimation {
  PlottingPanel panel = new PlottingPanel
                     ("t","y(t)","        Monte-Carlo Hamilton's principle");
  DrawingFrame  frame = new DrawingFrame(panel);
  PlotFrame plot= new PlotFrame("trial","Action (S)","Action versus trial number");

  private Control myControl;
  double [] t,y, yan, yg, yn, yp, T, V, L, S, per, N;
  double v0,v00,y0,y00,g,m,yf,t0,dy,tf,dt,S1,S2,dS;
  double pi, tol;
  int NPTS, NPTS2, NTRIALS,itrial,je,n,itrails=4;
  Trail [] trail=new Trail[itrails];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.000E00");
  Random random;

  public least_actionApp (){
    frame.setLocation(5,5);
    frame.setSize(400,295);
    frame.setTitle("Monte-Carlo Hamilton's principle");
    plot.setLocation(5,300);
    plot.setSize(400,295);
    plot.setConnected(true);
    plot.setName("Action versus trial number");
    nf.setMaximumFractionDigits(4);
    pi=Math.PI;
    for (int i=0; i<itrails; i++){
      trail[i] = new Trail();
      trail[i].setConnected(true);
    }
     trail[0].color=java.awt.Color.black; //analytic curve
     trail[0].setConnected(true);
     trail[1].color=java.awt.Color.blue;  //guess
     trail[1].setDashedStroke(1,8);
     trail[2].color=java.awt.Color.red;   //Least action
     trail[2].setDashedStroke(3,3);
     trail[3].color=java.awt.Color.blue;  //action
     trail[3].setConnected(true);
     trail[3].setDashedStroke(2,5);
  }

   protected void doStep() {
     //System's animation controller
     if (dS > tol && itrial < NTRIALS){
       calculate();
       panel.addDrawable(trail[2]);  //trail[2] is set up in calculate[]
       panel.setMessage("trials="+df.format(N[itrial]),0);
       panel.setMessage("Exchanges="+je,1);
       panel.setMessage("dS="+df.format(dS),3);
       frame.render();
       //plots S values versus trial number for each trial
       plot.setPreferredMinMax(0,itrial,S[itrial],S[0]);
       plot.addDrawable(trail[3]);
       plot.setMessage("S="+nf.format(S[itrial]),0);
       plot.render();
       itrial++;
     }else{
       stopAnimation();
       //2nd way - plot all points at once
       //ArrayLib.resize(N,itrial);//resize N array as needed
       //ArrayLib.resize(S,itrial);//resize S array as needed
       //plot.setPreferredMinMax(0,ArrayLib.max(N),S[itrial],S[0]);
       //plot.append(0,N,S); //works with resized arrays- not all points needed
       //3rd way - also plots all points at once
       //plot.setPreferredMinMax(0,ArrayLib.max(N),S[itrial-1],S[0]);
       //plot.addDrawable(trail[3]);
       //plot.render();
       myControl.println("Animation: Analytic result -black, guess - blue,");
       myControl.println("-----------Monte-Carlo, Hamilton's Principle - red");
       myControl.println("Lower graph: S versus trial number");
       myControl.println("Final Results");
       myControl.println("Total trials="+(itrial-1)+", Number of exchanges="+je);
       myControl.println("Final action change dS="+df.format(dS));
       myControl.println("Percent error="+df.format(per[itrial-1]));
     }
   }

   public void calc0 (){
     //guess, analytic solution
     for(int i=0; i<NPTS; i++){
       t[i]=t0+i*dt;                       //time array
       yan[i]=y0+v0*t[i]-g*t[i]*t[i]/2.;   //exact analytic trajectory
       yg[i]=y0+(yf-y0)*(t[i]-t0)/(tf-t0); //trajectory guess: interpolate end points
       yn[i]=yg[i];                        //initial trajectory is the guess
       trail[0].addPoint(t[i],yan[i]);
       trail[1].addPoint(t[i],yg[i]);
     }
     panel.setPreferredMinMax(t0,tf,ArrayLib.min(yan),ArrayLib.max(yan));
     panel.addDrawable(trail[0]);
     panel.addDrawable(trail[1]);

     //Kinetic and Potential energies, Lagrangian
     for(int i=0; i<NPTS2; i++){
       T[i]=m*Math.pow((yg[i+1]-yg[i])/dt,2.)/2.;
       V[i]=m*g*yg[i];
       L[i]=T[i]-V[i];
     }
     S1=dt*ArrayLib.sum(L);   //initial action, rectangular rule - use NPTS-1 terms
     //myControl.println("S1 calc="+df.format(S1));
     dS=Math.abs(S1);         //change in the action, initially it is equal to S1
   }


   public void calculate () {
       trail[2].clear(); //clear guess trail every time
       //start with trajectory in memory, use the previously saved trajectory
       for (int i = 0; i < NPTS; i++) {
         yp[i]=yn[i];
       }
       //begin new trajectory calculation process
       //we use random numbers - Monte Carlo method
       /*
       //one example of random numbers (double)
       double A=-5.0, B=5.0;
       for(int i=0; i<NPTS; i++){
         double number=A+(B-A)*random.nextDouble(); //random number between A and B
         myControl.println("number="+number);
       }
       */
       //y(0), y(NPTS-1) don't change so start at k=1, end at k=NPTS-2
       for (int k = 1; k < NPTS2; k++) {
         yp[k]=yp[k]+2.*dy*(random.nextDouble()-0.5); //modify kth step in trajectory
                                                      //memory, use Monte-Carlo step
         for (int j = 0; j < NPTS2; j++) {
           T[j]=m*Math.pow((yp[j+1]-yp[j])/dt,2.)/2.;//energies
           V[j]=m*g*yp[j];
           L[j]=T[j]-V[j];                           //Lagrangian
         }
         S2=dt*ArrayLib.sum(L);    //new action integral by rectangular rule - use NPTS-1 terms
         if (S2<S1){               //accept step if the action S decreases
           dS=Math.abs(S2-S1);     //variation of the action (must be small before stopping)
           S1=S2;                  //save the lower value of the action
           for (int p=1; p<NPTS2; p++) {
             yn[p]=yp[p];          //modify trajectory based on accepted step, and keep it
           }
           je=je+1;                //count exchanges made
         }
         else {
           for (int p=1; p<NPTS2; p++) {
             yp[p]=yn[p];           //if the step is not accepted, discard memory changes made,
                                    //and refresh the trajectory memory
           }
         }
       }
       //trial finished, save and create current trial plots
       S[itrial]=S1;                 //keep track of the action for each trial
       N[itrial]=itrial;             //trial number
       //percent error estimate
       double sum=0;
       for (int p=0; p<NPTS; p++) {
         sum=sum+(yn[p]-yan[p])*(yn[p]-yan[p])/(yan[p]*yan[p]+1.e-6);
         trail[2].addPoint(t[p],yn[p]);//for plotting purpose
       }
       trail[3].addPoint(N[itrial],S[itrial]); //for plotting purpose
       per[itrial]=100.*Math.sqrt(sum)/NPTS;   //percent error
       //myControl.println("itrial="+N[itrial]+", S="+nf.format(S[itrial])+
       //                  ", je="+je);
   }

   public void clear () {
     for (int i=0; i<itrails; i++){
       trail[i].clear();
     }
     panel.clear();
     frame.render();
     plot.clearDrawables();
     plot.render();
   }

   public void resetAnimation () {
     //initial conditions
     clear();
     myControl.clearMessages();
     myControl.println ("Simulates Hamilton's Least Action principle for a");
     myControl.println ("particle under the action of gravity. The trajectory");
     myControl.println ("is compared with what is expected analytically. The");
     myControl.println ("analytic result is y=y0+v0*t-g*t^/2. Hamilton's principle");
     myControl.println ("is based on the minimum of the action. The action is");
     myControl.println ("an integral of the Lagrangian over time, see text.");
     myControl.println ("Initialize: initializes; start: starts the simulation.");
     t    = null;
     y    = null;
     yan  = null;
     yg   = null;
     yn   = null;
     yp   = null;
     T    = null;
     V    = null;
     L    = null;
     S    = null;
     N    = null;
     tol=5.5e-8;         //tolerance reasonably small
     v0=5; g=9.8; m=1;   //initial speed, gravity and mass, given
     y0=0; yf=1;         //initial, final height, given
     t0=0; dy=.01;       //initial time, max change in y allowed
     v00=v0; y00=y0;     //backup values of v0 anf y0
     //final time estimate from analytic soln
     tf=v0/g+Math.sqrt((v0/g)*(v0/g)-2*(yf-y0)/g);
     NPTS=15;            //max number of points
     NTRIALS=1000;       //max number of trials
     delayTime=50;       //time in between animation steps (Abstract Animation)
     myControl.setValue("mass",m);
     myControl.setValue("v0",v0);
     myControl.setValue("g",g);
     myControl.setValue("y0",y0);
     myControl.setValue("yf",yf);
     myControl.setValue("t0",t0);
     myControl.setValue("tf",nf.format(tf));
     myControl.setValue("tolerance",tol);
     myControl.setValue("max change dy",dy);
     myControl.setValue("delayTime(ms)",delayTime);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("NTRIALS",NTRIALS);
     dt=(tf-t0)/(NPTS-1);
     myControl.setValue("dt[(tf-t0)/(NPTS-1)]=",nf.format(dt));
   }

   public void initializeAnimation() {
     //initial conditions
     random = new Random(19580427); //own seed - reproduce results
     //random = new Random(); //if using the system's own seed
     itrial=0;  //itrial number counter increased in doStep
     je=0;      //exchange counter increased in calculate()
     clear();
     myControl.clearMessages();
     m=myControl.getDouble("mass");
     v0=myControl.getDouble("v0");
     g=myControl.getDouble("g");
     y0=myControl.getDouble("y0");
     yf=myControl.getDouble("yf");
     t0=myControl.getDouble("t0");
     tf=myControl.getDouble("tf");
     tol=myControl.getDouble("tolerance");
     dy=myControl.getDouble("max change dy");
     delayTime=myControl.getInt("delayTime(ms)");
     NPTS=myControl.getInt("NPTS");
     NTRIALS=myControl.getInt("NTRIALS");
     dt=(tf-t0)/(NPTS-1);
     myControl.setValue("dt[(tf-t0)/(NPTS-1)]=",nf.format(dt));
     NPTS2=NPTS-1;
     if ((v0/g)*(v0/g)-2*(yf-y0)/g < 0){  //work with proper v0, y0, yf in this case
        myControl.println("v0 is too small, or yf is too large: values reset");
        myControl.setValue("v0",v00);
        myControl.setValue("y0",y00);
        v0=myControl.getDouble("v0");
        y0=myControl.getDouble("y0");
     }
     t    = new double[NPTS];
     y    = new double[NPTS];
     yan  = new double[NPTS];
     yg   = new double[NPTS];
     yn   = new double[NPTS];
     yp   = new double[NPTS];
     T    = new double[NPTS2];
     V    = new double[NPTS2];
     L    = new double[NPTS2];
     S    = new double[NTRIALS];
     per  = new double[NTRIALS];
     N    = new double[NTRIALS];
     calc0();
     myControl.println ("Press 'start' to get the animation going.");
     myControl.println("Upper left: Analytic result -black, guess - blue,");
     myControl.println("If running: Monte-Carlo, Hamilton's Principle - red");
     myControl.println("------------Lower graph: S versus trial number");
   }

   public void setControl(Control control) {
     myControl = control;
     resetAnimation();
   }

   public static void main(String[] args) {
     Animation app = new least_actionApp();
     AnimationControl myControl;
     myControl = new AnimationControl(app);
     myControl.setLocation(410, 5);
     myControl.setSize(360,590);
     myControl.setDividerLocation(285);
     app.setControl(myControl);
   }
}