/*
@Author J E Hasbun 2007.
Plots the underdamped HO solution.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.numerics.*;
import java.awt.Color;

public class under_dampApp implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("x axis ", "B*exp(t), x(t)",
                                "B*exp(t)-black, x(t)-red");
   DrawingFrame  drawingFramex = new DrawingFrame ("Plot", plotPanelx);
   DataTable dataTable = new DataTable();
   DataTableFrame framea= new DataTableFrame("Calculations", dataTable);
   DatasetManager datasets = new DatasetManager();
   Dataset dataseta = new Dataset ();
   Dataset datasetb = new Dataset ();
   private Control myControl;
   double [] xe, xu, tt;
   double dt,t0,x0,v0,tmax,gam,desc,m,k,c,B,w0,w,th;
   int NPTS;

   public under_dampApp(){
//-- plotting properties
     drawingFramex.setSize(350, 350);
     drawingFramex.setLocation(10, 5);
     plotPanelx.enableInspector(true);
 //-- datasetsa and its table properties
      dataseta.setSorted(true);//sorts the data
      dataseta.setConnected(true);
      dataseta.setMarkerShape(Dataset.PIXEL);
      dataseta.setMarkerColor(Color.black, Color.black, Color.black);
//-- datasetb properties
     datasetb.setSorted(true);//sorts the data
     datasetb.setConnected(true);
     datasetb.setMarkerShape(Dataset.PIXEL);
     datasetb.setLineColor(Color.RED);
//-- datasets and their table properties
     datasets.setSorted(true);
     datasets.setXYColumnNames(0, "t", "B*exp(t)");
     datasets.setSorted(true);
     datasets.setXYColumnNames(1, "", "x(t)");
     datasets.setXColumnVisible(1, false);
//-- dataTable properties
     dataTable.setRowNumberVisible(true);
     dataTable.setRowSelectionAllowed(true);
     dataTable.setMaximumFractionDigits(5);
     dataTable.add(datasets);
     dataTable.refreshTable();
     framea.setSize(350, 225);
     framea.setLocation(10, 355);
   }

   public void calculate() {
     myControl.clearMessages();
     dataseta.clear();
     datasets.clear();
     dataTable.clear();
     datasetb.clear();
     plotPanelx.clear();
     x0=myControl.getDouble("x0");
     v0=myControl.getDouble("v0");
     m=myControl.getDouble("m");
     k=myControl.getDouble("k");
     c=myControl.getDouble("c");
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     w0=Math.sqrt(k/m);
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     myControl.setValue("w0=",w0);
     tt = new double[NPTS+1];
     xe = new double[NPTS+1];
     xu = new double[NPTS+1];
     gam=c/2/m;
     desc=w0*w0-gam*gam;
     myControl.setValue("desc{w0^2-(c/2/m)^2}=",desc);
     if (desc <= 0){
        w=Double.NaN                           ;
        thereisaproblem();
     } else{
       w=Math.sqrt(desc);
       myControl.setValue("w=",w);
       B =Math.sqrt(x0*x0+(v0+gam*x0)*(v0+gam*x0)/w/w);
       th=Math.atan(w*x0/(v0+gam*x0));
       for (int i = 0; i <= NPTS; i++) {
         tt[i] = t0 + i * dt;          //time
         //xo=overdamped and xc=critically damped
         xe[i] = B*Math.exp(-gam*tt[i]);
         xu[i] =xe[i]*Math.sin(w*tt[i]+th);
         //System.out.println("t="+tt[i]+" xe="+xe[i]+" xu="+xu[i]);
       }
       dataseta.append(tt, xe);
       datasets.append(0, tt, xe);
       datasetb.append(tt, xu);
       datasets.append(1, tt, xu);
       dataTable.add(datasets);
       dataTable.refreshTable();
       plotPanelx.addDrawable(dataseta);
       plotPanelx.addDrawable(datasetb);
       plotPanelx.render();
     }
   }

   public void clear  () {
     plotPanelx.clear();
     plotPanelx.render();
     dataTable.clear();
     dataTable.refreshTable();
     dataseta.clear();
     datasets.clear();
     datasetb.clear();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Underdamped HO solution for which k/m > gam*gam.");
     myControl.println ("x0=initial position, v0=initial speed, m=mass");
     myControl.println ("k=spring constant, c=damping coefficient, t0=initial time");
     myControl.println ("tmax=upper time limit, NPTS=number of points plotted");
     myControl.println ("dt is the time step. Also, we have that");
     myControl.println ("x(t)=B*exp(-gam*t)*sin(w*t+th), with ");
     myControl.println ("B =sqrt(x0^2+(v0+gam*x0)^2/w^2), th=arctan(w*x0/(v0+gam*x0),");
     myControl.println ("gam=c/2/m, w0=sqrt(k/m), w=sqrt(desc), desc=w0^2-gam^2.");
     x0=1.0;
     v0=5.0;
     m=0.05;
     k=1;
     c=0.08;
     gam=c/2/m;
     t0=0.0;
     tmax=5.0;
     NPTS = 100;
     dt=(tmax-t0)/NPTS;
     w0=Math.sqrt(k/m);
     desc=w0*w0-gam*gam;
     w=Math.sqrt(desc);
     tt=null;
     xe=null;
     xu=null;
     myControl.setValue("x0", x0);
     myControl.setValue("v0", v0);
     myControl.setValue("m", m);
     myControl.setValue("k", k);
     myControl.setValue("c", c);
     myControl.setValue("t0", t0);
     myControl.setValue("tmax", tmax);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     myControl.setValue("w0=",w0);
     myControl.setValue("desc{w0^2-(c/2/m)^2}=",desc);
     myControl.setValue("w=",w);
   }

    public void thereisaproblem(){
       myControl.clearValues();
       myControl.clearMessages();
       myControl.println("*** c needs to be smaller since desc is <= zero ***");
       myControl.println("desc{w0^2-(c/2/m)^2}="+desc);
       myControl.println("Please reset or fix and re-calculate");
       try{Thread.sleep(100);//pause for a bit
       }catch(InterruptedException e){
       myControl.println(e.getMessage());}
       myControl.setValue("x0", x0);
       myControl.setValue("v0", v0);
       myControl.setValue("m", m);
       myControl.setValue("k", k);
       myControl.setValue("c", c);
       myControl.setValue("t0", t0);
       myControl.setValue("tmax", tmax);
       myControl.setValue("NPTS", NPTS);
       myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
       myControl.setValue("w0=",w0);
       myControl.setValue("desc{w0^2-(c/2/m)^2}=",desc);
       myControl.setValue("w=",w);
     }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new under_dampApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(370, 5);
     myControl.setSize(370,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
