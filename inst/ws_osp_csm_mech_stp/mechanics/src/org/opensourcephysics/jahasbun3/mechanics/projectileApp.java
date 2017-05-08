/*
@Author J E Hasbun 2007.
Plots the free fall projectile trajectories z(x) for various angles,where
 z(x)=zi+viz*(x-xi)/vix-0.5*g*((x-xi)/vix)^2.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
//import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.numerics.*;
import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;


public class projectileApp implements Calculation {
   PlotFrame plot= new PlotFrame("x","Z(x)","Trajectories for different angles");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double [][] z;
   double [] x, th, H, R, xh;
   double vi, zi, g, xi, xmax, xs, thi, thm, ths;
   double vix, viz, pi, cf, zmax, xmin;
   int Nx, Nth;

   public projectileApp(){
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
     thi=myControl.getDouble("init_angle: thi");
     thm=myControl.getDouble("final angle: thm");
     g=myControl.getDouble("gravity");
     xmin=myControl.getDouble("xmin");
     xmax=myControl.getDouble("xmax");
     if(xmax < xi){xmax=xmax+xi;}
     myControl.setValue("xmax",xmax);
     Nx=myControl.getInt("N_x_points");
     xs=(xmax-xmin)/Nx;
     myControl.setValue("dx=((xmax-xmin)/Nx)=",xs);
     Nth=myControl.getInt("Nth");
     ths=(thm-thi)/Nth;
     myControl.setValue("dth=((thm-thi)/Nth)=",ths);
     x   = new double[Nx];       //x coord
     xh  = new double[Nth+1];    //half ranges
     th  = new double[Nth+1];    //eccentricities
     H   = new double[Nth+1];    //maximum heights
     R   = new double[Nth+1];    //Ranges
     z   = new double[Nth+1][Nx];//z coords
     zmax=zi;
     for (int i=0; i <= Nth;i++){
       th[i] = thi+i*ths;             //eccentricities
       vix   = vi*Math.cos(th[i]*cf); //initial x speed
       viz   = vi*Math.sin(th[i]*cf); //initial z speed
       R[i]  = 2*viz*vix/g+xi;        //range
       xh[i] = xi+(R[i]-xi)/2;        //half range
       H[i]  = viz*viz/(2*g)+zi;      //max height
       if(H[i]>zmax){zmax=H[i];}      //largest of all heights used for plotting
         for (int j = 0; j < Nx; j++) {
           if(i==0){x[j]= xmin+j*xs;} //do only once - use for plotting range
           z[i][j] = zi+viz*(x[j]-xi)/vix-0.5*g*Math.pow((x[j]-xi)/vix,2);
         }
         //System.out.println("i="+i+", j="+j+", th="+th[i]+", x="+x[j]+
         //                     " z="+z[i][j]);
     }

       myControl.println("Identify the initial angles of the various curves");
       myControl.println("according to their associated ranges and maximum heights.");
       myControl.println("The maximum heights are indicated by the green dots.");
       plot.setPreferredMinMax(xmin,xmax,0,zmax);
       for (int i=0; i <= Nth;i++){
         plot.setConnected(i,true);
         plot.setMarkerSize(i, 0);//ith=dataset index, use size 0 for symbol
         plot.append(i, x, z[i]); //colors based on DatasetManager.getLineColor(i)
         myControl.println("Theta="+nf.format(th[i])+"\tRange="+nf.format(R[i])+
                           "\tMax Height="+nf.format(H[i]));
       }
       //place dots at the maximum heights
       plot.setConnected(Nth+1,false);//do not connect
       plot.setMarkerSize(Nth+1,2);   //Nth+1=dataset index, use size 2 for symbol
       plot.setMarkerShape(Nth+1,1);  //circle
       plot.setMarkerColor(Nth+1,Color.green,Color.red);
       plot.append(Nth+1,xh,H);       //Max height H at positions xh
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println("Plots the free fall projectile trajectories z(x) for various angles,");
     myControl.println("where z(x)=zi+viz*(x-xi)/vix-0.5*g*((x-xi)/vix)^2.");
     x    = null;
     z    = null;
     xi=0;              //initial position
     zi=0;              //initial height
     vi=5;              //initial velocity
     g=9.8;             //acceleration due to gravity
     thi=15; thm=75;    //angle range calculation
     Nth=6;
     vix=vi*Math.cos(45*cf); viz=vi*Math.sin(45*cf);
     xmax=2*viz*vix/g+xi;//maximum range to plot
     xmax=Math.ceil(xmax); xmin=xi;
     Nx=300;             //number of x points calculated
     myControl.setValue("xi",xi);
     myControl.setValue("zi",zi);
     myControl.setValue("init_vel",vi);
     myControl.setValue("init_angle: thi",thi);
     myControl.setValue("final angle: thm",thm);
     if(thi==0){thi=1.e-3;}
     if(thm==0){thm=1.e-3;}
     myControl.setValue("gravity",g);
     myControl.setValue("xmin",xmin);
     myControl.setValue("xmax",xmax);
     myControl.setValue("N_x_points",Nx);
     xs=(xmax-xmin)/Nx;
     myControl.setValue("dx=((xmax-xmin)/Nx)=",xs);
     myControl.setValue("Nth",Nth);
     ths=(thm-thi)/Nth;
     myControl.setValue("dth=((thm-thi)/Nth)=",ths);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new projectileApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
