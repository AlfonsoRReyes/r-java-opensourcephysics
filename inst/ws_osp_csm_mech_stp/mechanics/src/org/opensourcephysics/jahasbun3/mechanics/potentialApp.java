/*
@Author J E Hasbun 2007.
Plots an attractive potential, energy, etc., for a body under a central force.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class potentialApp implements Calculation {
   PlotFrame plot= new PlotFrame("r","Potential","Potential Plot");
   NumberFormat nf = NumberFormat.getInstance();
   Circle [] point=new Circle[3];
   private Control myControl;
   double [] r, V, CT, Veff,r2,E2;
   double m,ri,vi,thi,L,E,K,e,Ec,rc,rmin,rmax,dr,rlow,rhigh,pi,cf,vt,vr;
   int NPTS;

   public potentialApp(){
      plot.setConnected(true);
      plot.setLocation(10,5);
      plot.setSize(390,575);
      point[0]=new Circle();
      point[1]=new Circle();
      point[2]=new Circle();
      pi=Math.PI;
      cf=pi/180;
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     m=myControl.getDouble("m");
     ri=myControl.getDouble("ri");
     vi=myControl.getDouble("vi");
     thi=myControl.getDouble("thi");
     K=myControl.getDouble("K");
     NPTS=myControl.getInt("NPTS");
     rlow=myControl.getDouble("rlow");
     rhigh=myControl.getDouble("rhigh");
     dr=(rhigh-rlow)/NPTS;
     myControl.setValue("dr((rhigh-rlow)/NPTS)=",dr);
     r     = new double[NPTS];
     V     = new double[NPTS];
     Veff  = new double[NPTS];
     CT    = new double[NPTS];
     vt=vi*Math.sin(thi*cf);          //tangential speed
     vr=vi*Math.cos(thi*cf);          //radia speed
     L=m*ri*vt;                      //energy is constant
     E=m*vr*vr/2+L*L/m/ri/ri/2+K/ri; //angular momentum is constant
     e=Math.sqrt(1+2*E*L*L/m/K/K);   //eccentricity
     Ec=-m*K*K/L/L/2;                //cicular orbit (e=0 case) energy
     rc=-L*L/m/K;                    //cicular orbit (e=0 case) radius
     e=Math.sqrt(1+2*E*L*L/m/K/K);   //eccentricity
     rmin=1/(-m*K/L/L+Math.sqrt((m*K/L/L)*(m*K/L/L)+2*m*E/L/L));
     rmax=1/(-m*K/L/L-Math.sqrt((m*K/L/L)*(m*K/L/L)+2*m*E/L/L));
     for (int i=0; i < NPTS;i++){
       r[i]=rlow+i*dr;
       V[i]=K/r[i];                  //potential term
       CT[i]=L*L/r[i]/r[i]/2/m;      //centrifugal term
       Veff[i]=V[i]+CT[i];           //effective potential
     }
     r2=new double[2]; E2=new double[2];
     r2[0]=rmin; r2[1]=rmax; E2[0]=E; E2[1]=E; // straight line for the energy
     point[0].setXY(rmin,E);             //place a point at rmin
     point[0].color=java.awt.Color.black; point[0].pixRadius=2;
     point[1].setXY(rmax,E);             //place a point at rmax
     point[1].color=java.awt.Color.black; point[1].pixRadius=2;
     point[2].setXY(rc,Ec);              //place a point at the minimum of Veff
     point[2].color=java.awt.Color.blue;  point[2].pixRadius=2;
     plot.setPreferredMinMax(rlow,rhigh,-60,50);
     plot.setLineColor(0,Color.black);
     plot.setMarkerSize(0,0);
     plot.append(0,r,V);
     plot.setLineColor(1,Color.blue);
     plot.setMarkerSize(1,0);
     plot.append(1,r,CT);
     plot.setLineColor(2,Color.red);
     plot.setMarkerSize(2,0);
     plot.append(2,r,Veff);
     plot.setLineColor(3,Color.green);
     plot.setMarkerSize(3,0);
     plot.append(3,r2,E2);
     plot.addDrawable(point[0]);
     plot.addDrawable(point[1]);
     plot.addDrawable(point[2]);
     myControl.println("V(r): black, centrifugal term: blue, Veff(r): red");
     myControl.println("Ecentricity = \t"+nf.format(e));
     myControl.println("Angular Momentum L= \t"+nf.format(L));
     myControl.println("Energy E (green line)= \t"+nf.format(E));
     myControl.println("black dots: (rmin,rmax)= \t("+nf.format(rmin)+","+
                       nf.format(rmax)+")");
     myControl.println("Potential Min-blue dot: (rc,Ec)=\t("+nf.format(rc)+","+
                       nf.format(Ec)+")");
   }

   public void clear  () {
     plot.clearData();
     plot.clearDrawables();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots an attractive potential, energy, etc., for a body");
     myControl.println ("under a central force. In which case, we have the energy");
     myControl.println ("E=m*vr^2/2+L^2/m/ri^2/2+K/ri, where L=m*ri*vt, and where");
     myControl.println ("vt=vi*sin(thi), vr=vi*cos(thi). The eccentricity is");
     myControl.println ("e=sqrt(1+2*E*L^2/m/K^2), and V(r)=K/r, Veff(r)=V(r)+CT(r)");
     myControl.println ("CT(r)=L^2/r^2/2/m. We also have that the energy minimum is");
     myControl.println ("Ec=-m*K^2/L^2/2 and occurs ar r=rc=-L^2/m/K, finally,");
     myControl.println ("rmin=1/(-m*K/L^2+sqrt((m*K/L^2)*(m*K/L^2)+2*m*E/L^2))");
     myControl.println ("rmax=1/(-m*K/L^2-sqrt((m*K/L^2)*(m*K/L^2)+2*m*E/L^2)).");
     r     = null;
     V     = null;
     CT    = null;
     m=1.0; ri=1.0; K=-45.0; //mass, initial position, constant K
     vi=5; thi=75;           //initial velocity and angle
     NPTS=500;
     rlow=0.2; rhigh=4.0;
     myControl.setValue("m",m);
     myControl.setValue("ri",ri);
     myControl.setValue("vi",vi);
     myControl.setValue("thi",thi);
     myControl.setValue("K",K);
     myControl.setValue("NPTS",NPTS);
     myControl.setValue("rlow",rlow);
     myControl.setValue("rhigh",rhigh);
     dr=(rhigh-rlow)/NPTS;
     myControl.setValue("dr((rhigh-rlow)/NPTS)=",dr);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new potentialApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(350,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
