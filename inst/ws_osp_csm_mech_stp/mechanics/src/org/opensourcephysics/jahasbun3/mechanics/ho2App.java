/*
@Author J E Hasbun 2007.
Calculation of position, velocity, and acceleration for a body in
free fall with air resistance versus time. The equations of motion are used for small time
intervals.
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

public class ho2App implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("t axis ", "y axis", "y(t)");
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
   double [] yy, vv, aa, tt;
   double dt, t, t0, tmax, g, m, C, y0, v0, a0, F, Vt;
   int NPTS, model;

   public ho2App(){
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
     model=myControl.getInt("model(0 or 1)");
     g=myControl.getDouble("g");
     m = myControl.getDouble("m");
     C = myControl.getDouble("C");
     if(C <=0){
     C=1.e-9;
     myControl.setValue("C", C);}
     y0 = myControl.getDouble("y0");
     v0 = myControl.getDouble("v0");
     tmax = myControl.getDouble("tmax");
     NPTS = myControl.getInt("NPTS");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",dt);
     if(model ==0){
         F=-m*g-C*v0;
         Vt=Math.abs(m*g/C);}
       else if (model==1){
         F=-m*g-C*v0*Math.abs(v0);
         Vt=Math.sqrt(m*g/C);}
     myControl.setValue("Vt=",Vt);
     t0=0.0;
     tt = new double[NPTS+1];
     yy = new double[NPTS+1];
     vv = new double[NPTS+1];
     aa = new double[NPTS+1];
     a0=F/m;
     vv[0]=v0;
     yy[0]=y0;
     aa[0]=a0;
     tt[0]=t0;
     for(int i = 0; i <= NPTS-1; i++) {
       vv[i+1]=vv[i]+aa[i]*dt;                      //new velocity
       yy[i+1]=yy[i]+vv[i+1]*dt;                    //new position
       tt[i+1]=tt[i]+dt;                            //new time
       if(model==0){                               //new force
         F=-m*g-C*vv[i+1];}
       else if (model==1){
         F = -m * g - C * vv[i + 1] * Math.abs(vv[i + 1]);}
       else{
       resetCalculation();
       return;}
       aa[i+1]=F/m;                                 //new acceleration
     }
     dataseta.append(tt, yy);
     datasets.append(0, tt, yy);
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
     myControl.println ("Falling with gravity/air resistance - Euler-Cromer Method");
     myControl.println ("a=F/m, with F=-m*g-C*v (model=0), or F=-m*g-C*v^2 (model=1)");
     myControl.println ("v=v+a*dt, x=x+v*dt, t=t+dt");
     myControl.println ("g is the gravitational acceleration");
     myControl.println ("m is the particle's mass");
     myControl.println ("C is the drag coefficient");
     myControl.println ("y0 is the initial position");
     myControl.println ("v0 is the initial velocity");
     myControl.println ("tmax is the maximum time");
     myControl.println ("NPTS is the number of steps");
     myControl.println ("dt=tmax/NPTS is the time step used");
     myControl.println ("terminal speed Vt=(m*g/C) (model 0) or sqrt(m*g/C) (model=1)");
     t0=0;
     g = 9.8;
     m = 1.0;
     C = 0.05;
     y0 = 0.0;
     v0 = 110.0;
     tmax = 20.0;
     model=0;
     Vt=Math.abs(m*g/C);
     NPTS = 200;
     dt=tmax/NPTS;
     yy=null;
     vv=null;
     aa=null;
     tt=null;
     myControl.setValue("model(0 or 1)", model);
     myControl.setValue("g", g);
     myControl.setValue("m", m);
     myControl.setValue("C", C);
     myControl.setValue("y0", y0);
     myControl.setValue("v0", v0);
     myControl.setValue("tmax", tmax);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("dt(tmax/NPTS)=",dt);
     myControl.setValue("Vt=",Vt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new ho2App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(520, 30);
     myControl.setSize(250,505);
     model.setControl(myControl);
   }
}
