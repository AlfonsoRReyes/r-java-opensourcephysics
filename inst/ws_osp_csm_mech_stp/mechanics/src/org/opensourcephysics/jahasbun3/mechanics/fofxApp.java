/*
@Author J E Hasbun 2007.
Plots the kinetic, potential, total energy and force for a spring-mass system
assuming the mass is initially at the origin moving with initial speed v0.
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

public class fofxApp implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("x axis ", "KE, E-total", "KE, E-total");
   DrawingFrame  drawingFramex = new DrawingFrame ("Plot", plotPanelx);
   PlottingPanel plotPanely    = new PlottingPanel ("x axis ", "PE, E-total", "PE, E-total");
   DrawingFrame  drawingFramey = new DrawingFrame ("Plot", plotPanely);
   PlottingPanel plotPanelz    = new PlottingPanel ("x axis ", "F axis", "F(x)");
   DrawingFrame  drawingFramez = new DrawingFrame ("Plot", plotPanelz);
   DataTable dataTable = new DataTable();
   DataTableFrame framea= new DataTableFrame("Calculations", dataTable);
   DatasetManager datasets = new DatasetManager();
   Dataset dataseta = new Dataset ();
   Dataset datasetb = new Dataset ();
   Dataset datasetc = new Dataset ();
   Dataset datasetd = new Dataset ();
   private Control myControl;
   double [] KE, PE, EE, F, xx;
   double dt, t, t0, T, w, m, A, v0, k, E0, v;
   int NPTS;

   public fofxApp(){
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
     datasetc.setMarkerShape(Dataset.CIRCLE);
     datasetc.setMarkerColor(Color.red, Color.red, Color.red);
 //-- datasetd properties
      datasetd.setSorted(true);//sorts the data
      datasetd.setConnected(true);
      datasetd.setMarkerShape(Dataset.PIXEL);
      datasetd.setMarkerColor(Color.red, Color.red, Color.red);
//-- datasets and their table properties
     datasets.setSorted(true);
     datasets.setXYColumnNames(0, "x", "KE");
     datasets.setSorted(true);
     datasets.setXYColumnNames(1, "", "PE");
     datasets.setXColumnVisible(1, false);
     datasets.setXYColumnNames(2, "", "E");
     datasets.setXColumnVisible(2, false);
     datasets.setXYColumnNames(3, "", "F");
     datasets.setXColumnVisible(3, false);
//-- dataTable properties
     dataTable.setRowNumberVisible(false);
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
     datasetd.clear();
     plotPanelx.clear();
     plotPanely.clear();
     plotPanelz.clear();
     v0=myControl.getDouble("v0");
     m=myControl.getDouble("m");
     k=myControl.getDouble("k");
     NPTS=myControl.getInt("NPTS");
     w=Math.sqrt(k/m);
     T=2*Math.PI/w;
     E0=m*v0*v0/2;
     myControl.setValue("w",w);
     myControl.setValue("T",T);
     myControl.setValue("E0",E0);
     dt=T/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",dt);
     t0=0.0;
     A=Math.sqrt(2*E0/k);
     xx = new double[NPTS+1];
     KE = new double[NPTS+1];
     PE = new double[NPTS+1];
     EE = new double[NPTS+1];
     F = new double[NPTS+1];
     for(int i = 0; i <= NPTS; i++) {
       t=t0+i*dt;               //time
       xx[i]=A*Math.sin(w*t);   //position
       v=v0*Math.cos(w*t);      //velocity
       PE[i]=0.5*k*xx[i]*xx[i]; //Potential Energy
       KE[i]=0.5*m*v*v;         //Kinetic Energy
       EE[i]=PE[i]+KE[i];       //Total Energy
       F[i]=-k*xx[i];           //The force
       //System.out.println("x="+xx[i]+" v="+v+" PE="+PE[i]+" KE="+KE[i]);
     }
     dataseta.append(xx, KE);
     datasets.append(0, xx, KE);
     datasetb.append(xx, PE);
     datasets.append(1, xx, PE);
     datasetc.append(xx, EE);
     datasets.append(2, xx,EE);
     datasetd.append(xx,F);
     datasets.append(3,xx,F);
     dataTable.add(datasets);
     dataTable.refreshTable();
     plotPanelx.addDrawable(dataseta);
     plotPanelx.addDrawable(datasetc);
     plotPanely.addDrawable(datasetb);
     plotPanely.addDrawable(datasetc);
     plotPanelz.addDrawable(datasetd);
     plotPanelx.render();
     plotPanely.render();
     plotPanelz.render();
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
     datasetd.clear();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("spring-mass kinetic, Potential, total energies and force");
     myControl.println ("given that the mass is at x=0 initially, moving at v0");
     myControl.println ("x=sqrt(2*E0/k)*sin(w*t), v=v0*cos(w*t), total energy: E0=m*v0^2/2");
     myControl.println ("PE=0.5*x^2, KE=0.5*m*v^2, F=-k*x");
     myControl.println ("m is the particle's mass");
     myControl.println ("k is the spring constant, w=sqrt(k/m) is the frequency");
     myControl.println ("x0 is the initial position");
     myControl.println ("v0 is the initial velocity");
     myControl.println ("T=2*pi/w is the period");
     myControl.println ("NPTS is the number of steps");
     myControl.println ("dt=T/NPTS is the time step used");
     myControl.println ("E-total=PE+KE (red dots) and is constant (=E0) when conserved");

     t0=0;
     k = 0.01;
     m = 1.0;
     w=Math.sqrt(k/m);
     v0 = 0.5;
     NPTS = 100;
     E0=m*v0*v0/2;
     T=2*Math.PI/w;
     dt=T/NPTS;
     xx=null;
     KE=null;
     PE=null;
     F=null;
     myControl.setValue("v0", v0);
     myControl.setValue("m", m);
     myControl.setValue("k", k);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("w",w);
     myControl.setValue("T",T);
     myControl.setValue("E0",E0);
     myControl.setValue("dt(tmax/NPTS)=",dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new fofxApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(520, 30);
     myControl.setSize(250,505);
     model.setControl(myControl);
   }
}
