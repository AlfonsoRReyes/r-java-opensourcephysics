/*
@Author J E Hasbun 2007.
Plots possible conic section curves for various values of the eccentricity e.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class conic1App implements Calculation {
   PlotFrame plot= new PlotFrame(
   "x=r*cos(theta)","y=r*sin(theta)","r=rmin*(1+e)/(1+e*cos(theta)) versus theta for various e values");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double [][] x,y;
   double [] th, e;
   double rmin, dth, emin, emax, de, theta_min,theta_max,pi,r,L;
   int Nt, Ne;

   public conic1App(){
     plot.setLocation(10,5);
     plot.setSize(390,575);
     nf.setMaximumFractionDigits(3);
     pi=Math.PI;
   }

   public void calculate() {
     clear();
     myControl.clearMessages();
     rmin=myControl.getDouble("rmin");
     emin=myControl.getDouble("emin");
     emax=myControl.getDouble("emax");
     theta_min=myControl.getDouble("theta_min");
     theta_max=myControl.getDouble("theta_max");
     Ne=myControl.getInt("Ne");
     Nt=myControl.getInt("Nt");
     L=myControl.getDouble("view_range");
     de=(emax-emin)/(Ne-1);
     myControl.setValue("de=((emax-emin)/(Ne-1))=",de);
     dth=(theta_max-theta_min)/(Nt-1);
     myControl.setValue("dth=((theta_max-theta_min)/(Nt-1))=",dth);
     x  = new double[Ne][Nt];  //x coord
     y  = new double[Ne][Nt];  //y coord
     e  = new double[Ne];      //eccentricities
     th = new double[Nt];      //angle variable
     for (int i=0; i < Ne;i++){
       e[i] = emin+i*de;       //eccentricities
         for (int j = 0; j < Nt; j++) {
           if(i==0){th[j]=theta_min+j*dth;}
           r=rmin*(1+e[i])/(1+e[i]*Math.cos(th[j]));
           x[i][j]=r*Math.cos(th[j]);
           y[i][j]=r*Math.sin(th[j]);
         }
         //System.out.println("i="+i+", j="+j+", e="+e[i]+", " r="+r[i][j]);
     }
       myControl.println("Values of e for the curves are listed below. The curves for which");
       myControl.println("-1<e<0, are closed but smaller in radius than the closed ones");
       myControl.println("for which 0<e<1. When e<-1, the curves are hyperbolae and their");
       myControl.println("vertices remain close to each other. The curves are hyperbolae when");
       myControl.println("e>1 but their vertices get closer as e increased beyond e=1.");
       myControl.println("Also get a dot if e=-1, a circle if e=0, a parabola if e=1.");
       myControl.println("Finally notice that is a pair of curves for every |e|>1.");
       plot.setPreferredMinMax(-L,L,-L,L);
       plot.setSquareAspect(true);
       plot.setConnected(false);//applying to all the lines
       for (int i=0; i < Ne;i++){
         plot.setMarkerSize(i,1);//ith=dataset index, use size 0 for symbol
         plot.append(i, x[i], y[i]); //colors based on DatasetManager.getLineColor(i)
         myControl.println("e="+nf.format(e[i]));
       }
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots possible conic section curves for various values");
     myControl.println ("of the eccentricity e, where r=rmin*(1+e)/(1+e*cos(theta)),");
     myControl.println ("and x=r*cos(theta), y=r*sin(theta).");
     x    = null;
     y    = null;
     th   = null;
     e    = null;
     rmin=1;
     theta_min=0.0; theta_max=2*pi;
     emin=-1.8; emax=1.8;
     Nt=300; Ne=7;
     L=14;
     myControl.setValue("rmin",rmin);
     myControl.setValue("emin",emin);
     myControl.setValue("emax",emax);
     myControl.setValue("theta_min",theta_min);
     myControl.setValue("theta_max",theta_max);
     myControl.setValue("Ne",Ne);
     myControl.setValue("Nt",Nt);
     myControl.setValue("view_range",L);
     de=(emax-emin)/(Ne-1);
     myControl.setValue("de=((emax-emin)/(Ne-1))=",de);
     dth=(theta_max-theta_min)/(Nt-1);
     myControl.setValue("dth=((theta_max-theta_min)/(Nt-1))=",dth);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new conic1App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(240);
     model.setControl(myControl);
   }
}
