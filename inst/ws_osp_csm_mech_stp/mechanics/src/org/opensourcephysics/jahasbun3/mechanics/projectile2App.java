/*
@Author J E Hasbun 2007.
Compares the free fall w/o drag projectile trajectories z(x) for a given angle.
The free fall case is z(x)=zi+viz*(x-xi)/vix-0.5*g*((x-xi)/vix)^2. The case with
drag is z(x)=(m*g/c+viz)*(x-xi)/vix+m^2*g*log(1-c*(x-xi)/m/vix)/c^2.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class projectile2App implements Calculation {
   PlotFrame plot= new PlotFrame("x","Z(x)","Trajectories for different angles");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double [][] z;
   double [] x;
   double vi, zi, g, xi, xmax, xs, th, R, m, tland;
   double vix, viz, pi, cf, zmax, xmin, c, Dtland, DR;
   int Nx;
   int maxiter=25;
   double tolerance=5e-3;
   Bisnewt bisnewt=new Bisnewt();
   Function ft=new fyoft();

   public projectile2App(){
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
      pi=Math.PI; cf=pi/180;
   }

   public void calculate() {
     myControl.clearMessages();
     xi=myControl.getDouble("xi");
     zi=myControl.getDouble("zi");
     vi=myControl.getDouble("init_vel");
     th=myControl.getDouble("init_angle: th");
     m=myControl.getDouble("mass");
     c=myControl.getDouble("Drag_coeff:");
     g=myControl.getDouble("gravity");
     xmin=myControl.getDouble("xmin");
     xmax=myControl.getDouble("xmax");
     if(xmax < xi){xmax=xmax+xi;}
     myControl.setValue("xmax",xmax);
     Nx=myControl.getInt("N_x_points");
     xs=(xmax-xmin)/Nx;
     myControl.setValue("dx=((xmax-xmin)/Nx)=",xs);
     x   = new double[Nx];       //x coord
     z   = new double[2][Nx];//z coords
     vix   = vi*Math.cos(th*cf); //initial x speed
     viz   = vi*Math.sin(th*cf); //initial z speed
     R  = 2*viz*vix/g+xi;        //free fall range
     zmax=viz*viz/(2*g)+zi;      //free fall max height
     tland=2*viz/g;              //free fall landing time
     double x_initial=0.001, x_final=tland;
     //solve for the landing time iteratively (see function fyoft)
     Dtland=bisnewt.Bisnewt(x_initial, x_final, maxiter, tolerance,ft);
     DR=m*vix*(1-Math.exp(-c*Dtland/m))/c; //range when there is drag
     for (int j = 0; j < Nx; j++) {
       x[j]= xmin+j*xs;
       z[0][j] = zi+viz*(x[j]-xi)/vix-
                 0.5*g*Math.pow((x[j]-xi)/vix,2);      //free fall
       z[1][j]=(m*g/c+viz)*(x[j]-xi)/vix+
               m*m*g*Math.log(1-c*(x[j]-xi)/m/vix)/c/c;//with drag
     }
     //System.out.println("i="+i+", j="+j+", th="+th[i]+", x="+x[j]+
     //                     " z="+z[i][j]);
     plot.setPreferredMinMax(xmin,xmax,0,zmax);
     plot.setConnected(0,true);
     plot.setMarkerSize(0, 0);//0-dataset index, use size 0 for symbol
     plot.append(0, x, z[0]); //colors based on DatasetManager.getLineColor(i)
     plot.setLineColor(0,java.awt.Color.blue);
     plot.setConnected(1,true);
     plot.setMarkerSize(1, 0);//1-dataset index, use size 0 for symbol
     plot.append(1, x, z[1]); //colors based on DatasetManager.getLineColor(i)
     plot.setLineColor(1,java.awt.Color.red);
     myControl.println("The blue curve is due to free fall and the red curve is the");
     myControl.println("projectile motion with drag present. The associated landing");
     myControl.println("times and ranges follow.");
     myControl.println("Free Fall: t_land="+nf.format(tland)+"\tRange="+nf.format(R)+"\n"+
                       "with Drag: t_land="+nf.format(Dtland)+"\tRange="+nf.format(Dtland));
   }

   public class fyoft implements Function {
     public double evaluate(double t) {
       //defines y(t=tmax)=0 on landing, here t=time
       double f;
       f=(m/c+viz/g)*(1-Math.exp(-c*t/m))-t;
       return f;
     }
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("Compares the free fall w/o drag projectile trajectories z(x)");
     myControl.println("for a given angle. The free fall case is");
     myControl.println("z(x)=zi+viz*(x-xi)/vix-0.5*g*((x-xi)/vix)^2. The case with drag is");
     myControl.println("z(x)=(m*g/c+viz)*(x-xi)/vix+m^2*g*log(1-c*(x-xi)/m/vix)/c^2");
     x    = null;
     z    = null;
     m=1; c=0.5;        //particle mass, and drag coefficient
     xi=0;              //initial position
     zi=0;              //initial height
     vi=25;             //initial velocity
     g=9.8;             //acceleration due to gravity
     th=45;             //initialangle
     vix=vi*Math.cos(45*cf); viz=vi*Math.sin(45*cf);
     xmax=2*viz*vix/g+xi;//maximum range to plot
     xmax=Math.ceil(xmax); xmin=xi;
     Nx=300;             //number of x points calculated
     myControl.setValue("xi",xi);
     myControl.setValue("zi",zi);
     myControl.setValue("init_vel",vi);
     myControl.setValue("init_angle: th",th);
     myControl.setValue("mass",m);
     myControl.setValue("Drag_coeff:",c);
     if(th==0){th=1.e-3;}
     myControl.setValue("gravity",g);
     myControl.setValue("xmin",xmin);
     myControl.setValue("xmax",xmax);
     myControl.setValue("N_x_points",Nx);
     xs=(xmax-xmin)/Nx;
     myControl.setValue("dx=((xmax-xmin)/Nx)=",xs);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new projectile2App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
