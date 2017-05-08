/*
@Author J E Hasbun 2007.
Plots the position, speed, and acceleration due to a force=F0*cos(w*t)
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

public class foftApp implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("t axis ", "x axis", "x(t)");
   DrawingFrame  drawingFramex = new DrawingFrame ("Plot", plotPanelx);
   PlottingPanel plotPanely    = new PlottingPanel ("t axis ", "v axis", "v(t)");
   DrawingFrame  drawingFramey = new DrawingFrame ("Plot", plotPanely);
   PlottingPanel plotPanelz    = new PlottingPanel ("t axis ", "a axis", "a(t)");
   DrawingFrame  drawingFramez = new DrawingFrame ("Plot", plotPanelz);
   DataTable dataTable = new DataTable();
   DataTableFrame framea= new DataTableFrame("Calculations", dataTable);
   DatasetManager datasets = new DatasetManager();
   Dataset dataseta = new Dataset ();
   Dataset datasetb = new Dataset ();
   Dataset datasetc = new Dataset ();
   private Control myControl;
   double [] xx, vv, aa, tt;
   double dt, t, t0, tmax, w, m, x0, v0, a0, F0;
   int NPTS;

   public foftApp(){
//-- plotting properties
     drawingFramex.setSize(250, 250);
     drawingFramex.setLocation(1, 30);
     drawingFramey.setSize(250, 250);
     drawingFramey.setLocation(260, 30);
     drawingFramez.setSize(250, 250);
     drawingFramez.setLocation(260, 285);
     plotPanelx.enableInspector(true);
     plotPanely.enableInspector(true);
     plotPanelz.enableInspector(true);
 //-- datasetsa and its table properties
      dataseta.setSorted(true);//sorts the data
      dataseta.setConnected(true);
      dataseta.setMarkerShape(Dataset.PIXEL);
      dataseta.setMarkerColor(Color.red, Color.red, Color.red);
//-- datasetb properties
     datasetb.setSorted(true);//sorts the data
     datasetb.setConnected(true);
     datasetb.setMarkerShape(Dataset.PIXEL);
     datasetb.setMarkerColor(Color.red, Color.red, Color.red);
//-- datasetc properties
     datasetc.setSorted(true);//sorts the data
     datasetc.setConnected(true);
     datasetc.setMarkerShape(Dataset.PIXEL);
     datasetc.setMarkerColor(Color.red, Color.red, Color.red);
//-- datasets and their table properties
     datasets.setSorted(true);
     datasets.setXYColumnNames(0, "t", "x");
     datasets.setSorted(true);
     datasets.setXYColumnNames(1, "", "v");
     datasets.setXColumnVisible(1, false);
     datasets.setXYColumnNames(2, "", "a");
     datasets.setXColumnVisible(2, false);
//-- dataTable properties
     dataTable.setRowNumberVisible(true);
     dataTable.setRowSelectionAllowed(true);
     dataTable.setMaximumFractionDigits(5);
     dataTable.add(datasets);
     dataTable.refreshTable();
     framea.setSize(250, 250);
     framea.setLocation(1, 285);
   }

   public void calculate() {
     dataseta.clear();
     datasets.clear();
     dataTable.clear();
     datasetb.clear();
     datasetc.clear();
     plotPanelx.clear();
     plotPanely.clear();
     plotPanelz.clear();
     x0=myControl.getDouble("x0");
     v0=myControl.getDouble("v0");
     m=myControl.getDouble("m");
     w=myControl.getDouble("w");
     F0=myControl.getDouble("F0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",dt);
     t0=0.0;
     tt = new double[NPTS+1];
     xx = new double[NPTS+1];
     vv = new double[NPTS+1];
     aa = new double[NPTS+1];
     for(int i = 0; i <= NPTS-1; i++) {
       tt[i]=t0+i*dt;                       //time
       aa[i]=F0*Math.cos(w*tt[i])/m;      //acceleration
       vv[i]=v0+F0*Math.sin(w*tt[i])/m/w; //velocity
       xx[i]=x0+v0*tt[i]-F0*(Math.cos(w*tt[i])-1)/m/w/w; //position
       //System.out.println("xx="+xx[i]+" vv="+vv[i]); //if needed
     }
     dataseta.append(tt, xx);
     datasets.append(0, tt, xx);
     datasetb.append(tt, vv);
     datasets.append(1, tt, vv);
     datasetc.append(tt, aa);
     datasets.append(2, tt, aa);
     dataTable.add(datasets);
     dataTable.refreshTable();
     plotPanelx.addDrawable(dataseta);
     plotPanely.addDrawable(datasetb);
     plotPanelz.addDrawable(datasetc);
     plotPanelx.render();
     plotPanelz.render();
     plotPanely.render();
   }

   public void clear  () {
     plotPanelx.clear();
     plotPanelx.render();
     plotPanely.clear();
     plotPanely.render();
     plotPanelz.clear();
     plotPanelz.render();
     dataTable.clear();
     dataTable.add(datasets);
     dataTable.refreshTable();
     dataseta.clear();
     datasets.clear();
     datasetb.clear();
     datasetc.clear();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Time dependent force (F0*cos(w*t)) plot for x, v, a");
     myControl.println ("a=F0*cos(w*t)/m,  v=v0+F0*sin(w*t)/m/w");
     myControl.println ("x=x0+v0*t-F0*(cos(w*t)-1)/m/w/w");
     myControl.println ("m is the particle's mass");
     myControl.println ("F0 is the force amplitude, w is the frequency");
     myControl.println ("x0 is the initial position");
     myControl.println ("v0 is the initial velocity");
     myControl.println ("tmax is the maximum time");
     myControl.println ("NPTS is the number of steps");
     myControl.println ("dt=tmax/NPTS is the time step used");
     t0=0;
     F0 = 1.0;
     m = 1.0;
     w = 3.0;
     x0 = 0.0;
     v0 = 0.05;
     tmax = 10.0;
     NPTS = 100;
     dt=tmax/NPTS;
     xx=null;
     vv=null;
     aa=null;
     tt=null;
     myControl.setValue("x0", x0);
     myControl.setValue("v0", v0);
     myControl.setValue("m", m);
     myControl.setValue("w", w);
     myControl.setValue("F0", F0);
     myControl.setValue("tmax", tmax);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("dt(tmax/NPTS)=",dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new foftApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(520, 30);
     myControl.setSize(250,505);
     model.setControl(myControl);
   }
}
