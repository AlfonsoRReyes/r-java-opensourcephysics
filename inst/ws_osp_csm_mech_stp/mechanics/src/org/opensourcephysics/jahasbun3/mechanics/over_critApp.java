/*
@Author J E Hasbun 2007.
Plots the overdamped and critically damped HO solutions.
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

public class over_critApp implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("x axis ", "Xo(t), Xc(t)",
                                "Xo(t)-black, Xc(t)-red");
   DrawingFrame  drawingFramex = new DrawingFrame ("Plot", plotPanelx);
   DataTable dataTable = new DataTable();
   DataTableFrame framea= new DataTableFrame("Calculations", dataTable);
   DatasetManager datasets = new DatasetManager();
   Dataset dataseta = new Dataset ();
   Dataset datasetb = new Dataset ();
   private Control myControl;
   double [] xo, xc, tt;
   double dt,t0,x0,v0,tmax,gam,desc,Bo,Ao,Ac,Bc,gam1,gam2,m,k,c;
   int NPTS;

   public over_critApp(){
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
     datasets.setXYColumnNames(0, "t", "Xo(t)");
     datasets.setSorted(true);
     datasets.setXYColumnNames(1, "", "Xc(t)");
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
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     tt = new double[NPTS+1];
     xo = new double[NPTS+1];
     xc = new double[NPTS+1];
     gam=c/2/m;
     desc=gam*gam-k/m;
     myControl.setValue("desc{(c/2/m)^2-k/m}=",desc);
     if (desc <= 0){
        thereisaproblem();
     } else{
       gam1 = gam + Math.sqrt(desc);
       gam2 = gam - Math.sqrt(desc);
       Bo = (v0 + gam1 * x0) / (gam1 - gam2);
       Ao = x0 - Bo;
       Ac = x0;
       Bc = v0 + gam * x0;
       for (int i = 0; i <= NPTS; i++) {
         tt[i] = t0 + i * dt; //time
         //xo=overdamped and xc=critically damped
         xo[i] = Ao * Math.exp( -gam1 * tt[i]) + Bo * Math.exp( -gam2 * tt[i]);
         xc[i] = Ac * Math.exp( -gam * tt[i]) +
             Bc * tt[i] * Math.exp( -gam * tt[i]);
         //System.out.println("t="+tt[i]+" xo="+xo[i]+" xc="+xc[i]);
       }
       dataseta.append(tt, xo);
       datasets.append(0, tt, xo);
       datasetb.append(tt, xc);
       datasets.append(1, tt, xc);
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
     myControl.println ("Plots the overdamped and critically damped HO solutions");
     myControl.println ("for which gam*gam is > and = k/m, respectively.");
     myControl.println ("x0=initial position, v0=initial speed, m=mass");
     myControl.println ("k=spring constant, c=damping coefficient, t0=initial time");
     myControl.println ("tmax=upper time limit, NPTS=number of points plotted");
     myControl.println ("dt is the time step. Also, we have that");
     myControl.println ("xo=Ao*exp(-gam1*t)+Bo*exp(-gam2*t),");
     myControl.println ("xc= Ac*exp(-gam*t+Bc*t*exp(-gam*t), with ");
     myControl.println ("gam1=gam+sqrt(desc), gam2=gam-sqrt(desc), gam=c/2/m,");
     myControl.println ("desc=gam*gam-k/m, Ao=x0-Bo, Bo=(v0+gam1*x0)/(gam1-gam2),");
     myControl.println ("Ac=x0, and Bc=v0+gam*x0.");
     x0=1.0;
     v0=5.0;
     m=0.05;
     k=1;
     c=0.5;
     gam=c/2/m;
     desc=gam*gam-k/m;
     t0=0.0;
     tmax=2.0;
     NPTS = 100;
     dt=(tmax-t0)/NPTS;
     tt=null;
     xo=null;
     xc=null;
     myControl.setValue("x0", x0);
     myControl.setValue("v0", v0);
     myControl.setValue("m", m);
     myControl.setValue("k", k);
     myControl.setValue("c", c);
     myControl.setValue("t0", t0);
     myControl.setValue("tmax", tmax);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("dt((tmax-t0)/NPTS)=",dt);
     myControl.setValue("desc{(c/2/m)^2-k/m}=",desc);
   }

    public void thereisaproblem(){
       myControl.clearValues();
       myControl.clearMessages();
       myControl.println("*** c needs to be higher since desc is <= zero ***");
       myControl.println("desc{(c/2/m)^2-k/m}="+desc);
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
       myControl.setValue("desc{(c/2/m)^2-k/m}=",desc);
    }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new over_critApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(370, 5);
     myControl.setSize(370,575);
     model.setControl(myControl);
   }
}
