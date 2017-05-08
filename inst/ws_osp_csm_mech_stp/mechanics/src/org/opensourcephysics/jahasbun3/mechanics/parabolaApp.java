/*
@Author J E Hasbun 2007.
Plots parabolas: (x-x0)^2=4e(z-z0) whose curvatures=1/2e with vertex (x0,z0)
and focal points (x0,z0+e).
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class parabolaApp implements Calculation {
   PlotFrame plot= new PlotFrame("x","Z(x)","Parabolas: (x-x0)^2=4e(z-z0)");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double [][] z;
   double [] x, e, f, xf;
   double x0, z0, xmin, xmax, dx, emin, emax, de;
   int Nx, Ne;

   public parabolaApp(){
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     x0=myControl.getDouble("x0");
     z0=myControl.getDouble("z0");
     xmin=myControl.getDouble("xmin");
     xmax=myControl.getDouble("xmax");
     Nx=myControl.getInt("Nx");
     dx=(xmax-xmin)/Nx;
     myControl.setValue("dx=((xmax-xmin)/Nx)=",dx);
     emin=myControl.getDouble("emin");
     emax=myControl.getDouble("emax");
     Ne=myControl.getInt("Ne");
     de=(emax-emin)/Ne;
     myControl.setValue("de=((emax-emin)/Ne)=",de);
     x  = new double[Nx];    //x coord
     xf = new double[Ne+1];    //will store x0 positions for the focci
     e  = new double[Ne+1];    //eccentricities
     f  = new double[Ne+1];    //focci
     z  = new double[Ne+1][Nx];//z coords

     for (int i=0; i <= Ne;i++){
       e[i] = emin+i*de;//eccentricities
       f[i] = z0+e[i];  //focci
       xf[i]= x0;
         for (int j = 0; j < Nx; j++) {
           x[j]    = xmin+j*dx+x0;
           z[i][j] = z0+(x[j]-x0)*(x[j]-x0)/4/e[i];
         }
         //System.out.println("i="+i+", j="+j+", e="+e[i]+", x="+x[j]+
         //                     " z="+z[i][j]);
     }

       myControl.println("Values of e and respective focal points (green dots) corresponding");
       myControl.println("to each curve. To identify the curves, note that the farther the");
       myControl.println("focal point is from the vertex (x0,z0), the flatter is the parabola.");
       myControl.println("Also, low values of eccentricity (e) correspond to large curvatures.");
       for (int i=0; i <= Ne;i++){
         plot.setConnected(i,true);
         plot.setMarkerSize(i, 0);//ith=dataset index, use size 0 for symbol
         plot.append(i, x, z[i]); //colors based on DatasetManager.getLineColor(i)
         myControl.println("e="+nf.format(e[i])+",\t f(e="+
                            nf.format(e[i])+")="+nf.format(f[i]));
       }
       //place dots at the focal points (x0,z0+e) for each parabola
       plot.setConnected(Ne+1,false);//do not connect
       plot.setMarkerSize(Ne+1,3);   //Ne+1=dataset index, use size 3 for symbol
       plot.setMarkerShape(Ne+1,1);  //circle
       plot.setMarkerColor(Ne+1,Color.green,Color.red);
       plot.append(Ne+1,xf,f);       //focal point f at positions xf
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots parabolas: (x-x0)^2=4e(z-z0) whose curvatures=1/2e");
     myControl.println ("with vertex (x0,z0) and focal points (x0,z0+e)");
     x    = null;
     z    = null;
     x0=1;
     z0=1;
     xmin=-3; xmax=3;
     emin=-0.5; emax=-1.5;
     Nx=60; Ne=4;
     myControl.setValue("x0",x0);
     myControl.setValue("z0",z0);
     myControl.setValue("xmin",xmin);
     myControl.setValue("xmax",xmax);
     myControl.setValue("Nx",Nx);
     dx=(xmax-xmin)/Nx;
     myControl.setValue("dx=((xmax-xmin)/Nx)=",dx);
     myControl.setValue("emin",emin);
     myControl.setValue("emax",emax);
     myControl.setValue("Ne",Ne);
     de=(emax-emin)/Ne;
     myControl.setValue("de=((emax-emin)/Ne)=",de);

   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new parabolaApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
