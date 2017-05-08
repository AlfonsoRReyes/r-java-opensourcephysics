/*
@Author J E Hasbun 2007.
Plots Plots the coordinate solutions for the single mode couple spring mass
system without walls.
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

public class inter_spr1App implements Calculation {
   PlotFrame plot=new PlotFrame
       ("time","Position","Single Mode Spring-Mass system without walls");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double x1[], x2[], xcm[], tt[];
   double m1,m2,k0,x10,x20,v10,v20,mu,xcm0,vcm;
   double xr,xr0,vr0,w0,tau,A,B,t0,tmax,dt;
   int NPTS;

   public inter_spr1App(){
      //plot.setVisible(true);
      plot.setConnected(true);// set all plot's connected
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     m1=myControl.getDouble("m1");
     m2=myControl.getDouble("m2");
     k0=myControl.getDouble("k0");
     x10=myControl.getDouble("x10");
     x20=myControl.getDouble("x20");
     v10=myControl.getDouble("v10");
     v20=myControl.getDouble("v20");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     mu=m1*m2/(m1+m2);            //reduced mass
     w0=Math.sqrt(k0/mu);          //frequency
     tau=2*Math.PI/w0;            //period
     xcm0=(m1*x10+m2*x20)/(m1+m2);//initial cm
     vcm=(m1*v10+m2*v20)/(m1+m2); //vcm speed
     xr0=(x20-x10);               //relative coordinate
     vr0=(v20-v10);               //relative velocity
     A=vr0/w0;                    //amplitudes
     B=xr0;
     x1     = new double[NPTS];
     x2     = new double[NPTS];
     xcm    = new double[NPTS];
     tt     = new double[NPTS];
     for (int i=0; i<=NPTS-1;i++){
       tt[i]=t0+i*dt;
       xr=A*Math.sin(w0*tt[i])+B*Math.cos(w0*tt[i]);
       xcm[i]=xcm0+vcm*tt[i];      //center of mass
       x1[i]=xcm[i]-m2*xr/(m1+m2); //mass 1 position
       x2[i]=xcm[i]+m1*xr/(m1+m2); //mass 2 position
     }
     myControl.println("red: center of mass, black: mass 1, blue: mass 2");
     myControl.println("reduced mass mu="+nf.format(mu)+"\n"+
                       "frequency w0=sqrt(k/mu)="+nf.format(w0)+"\n"+
                       "period tau=2*pi/w0="+nf.format(tau)+"\n"+
                       "initial cm position xcm0="+nf.format(xcm0)+"\n"+
                       "initial cm speed vcm=vcm0="+nf.format(vcm)+"\n"+
                       "initial relative position xr0="+nf.format(xr0)+"\n"+
                       "initial relative speed vr0="+nf.format(vr0));
     plot.setLineColor(0,Color.red);
     plot.setMarkerSize(0,0);
     plot.append(0, tt, xcm);
     plot.setLineColor(1,Color.black);
     plot.setMarkerSize(1,0);
     plot.append(1, tt, x1);
     plot.setLineColor(2,Color.blue);
     plot.setMarkerSize(2,0);
     plot.append(2, tt, x2);
     }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the coordinate solutions for the single mode");
     myControl.println ("coupled spring mass system without walls");
     myControl.println ("k0=spring constant, m1=mass 1, m2= mass 2, w0=sqrt(k/mu)");
     myControl.println ("is the frequency , mu=m1*m2/(m1+m2) is the reduced mass.");
     myControl.println ("The center of mass posiiton is xcm=xcm0+vcm*t, where");
     myControl.println ("xcm0=(m1*x10+m2*x20)/(m1+m2), vcm=(m1*v10+m2*v20)/(m1+m2).");
     myControl.println ("The mass' positions are x1=xcm-m2*xr/(m1+m2),");
     myControl.println ("x2=xcm+m1*xr/(m1+m2), where xr=A*sin(w0*t)+B*cos(w0*t)");
     myControl.println ("xr0=(x20-x10), vr0=(v20-v10), A=vr0/w0, B=xr0");
     myControl.println ("Initial mass positions & speeds are x10, x20, v10, v20");
     x1    = null;
     x2    = null;
     xcm   = null;
     tt    = null;
     m1=1; m2=2; k0=0.5; x10=1; x20=-1; v10=0.02; v20=0.04;
     mu=m1*m2/(m1+m2); w0=Math.sqrt(k0/mu); tau=2*Math.PI/w0;
     NPTS=75; t0=0; tmax=2*tau;
     myControl.setValue("m1",m1);
     myControl.setValue("m2",m2);
     myControl.setValue("k0",k0);
     myControl.setValue("x10",x10);
     myControl.setValue("x20",x20);
     myControl.setValue("v10",v10);
     myControl.setValue("v20",v20);
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("NPTS",NPTS);
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new inter_spr1App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
