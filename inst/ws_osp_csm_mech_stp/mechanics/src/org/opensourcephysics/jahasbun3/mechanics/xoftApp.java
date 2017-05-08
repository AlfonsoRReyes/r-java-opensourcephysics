/*
@Author J E Hasbun 2007.
Plots Plots a function of time on interval [t0,tmax].
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import java.awt.Color;

public class xoftApp implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("t axis ", "f axis", "f(t)");
   DrawingFrame  drawingFramex = new DrawingFrame ("Plot", plotPanelx);
   DataTable dataTable = new DataTable();
   DataTableFrame framea= new DataTableFrame("Calculations", dataTable);
   DatasetManager datasets = new DatasetManager();
   Dataset dataseta = new Dataset ();
   private Control myControl;
   double [] xx,tt;
   double dt,t,t0,tau,tmax;
   int NPTS;

   public xoftApp(){
//-- plotting properties
     drawingFramex.setSize(250, 250);
     drawingFramex.setLocation(10, 30);
     plotPanelx.enableInspector(true);
 //-- dataseta and its table properties
      dataseta.setSorted(true);//sorts the data
      dataseta.setConnected(true);
      dataseta.setMarkerShape(Dataset.PIXEL);
      dataseta.setMarkerColor(Color.red, Color.red, Color.red);
//-- datasets and their table properties
     datasets.setSorted(true);
     datasets.setXYColumnNames(0, "t", "f");
     datasets.setSorted(true);
//-- dataTable properties
     dataTable.setRowNumberVisible(true);
     dataTable.setRowSelectionAllowed(true);
     dataTable.setMaximumFractionDigits(5);
     dataTable.add(datasets);
     dataTable.refreshTable();
     framea.setSize(250, 250);
     framea.setLocation(10, 285);
   }

   public void calculate() {
     Function f;
     dataseta.clear();
     datasets.clear();
     dataTable.clear();
     plotPanelx.clear();
     String fx = myControl.getString("f(x=time)");
     try { // read in function
     f = new ParsedFunction(fx);
     } catch(ParserException ex) {
       myControl.println(ex.getMessage());
       return;
     }
     t0=myControl.getDouble("t0");
     tmax=myControl.getDouble("tmax");
     NPTS=myControl.getInt("NPTS");
     dt=(tmax-t0)/NPTS;
     myControl.setValue("dt([tmax-t0]/NPTS)=",dt);
     xx = new double[NPTS+1];
     tt = new double[NPTS+1];
     for(int i = 0; i <= NPTS; i++) {
       tt[i]=t0+i*dt;               //time
       xx[i]=f.evaluate(tt[i]);   //position
       //System.out.println("t="+tt[i]+" x="+xx[i]);
     }
     dataseta.append(tt, xx);
     datasets.append(0, tt, xx);
     dataTable.add(datasets);
     dataTable.refreshTable();
     plotPanelx.addDrawable(dataseta);
   }

   public void clear  () {
     plotPanelx.clear();
     plotPanelx.render();
     dataTable.clear();
     dataTable.add(datasets);
     dataTable.refreshTable();
     dataseta.clear();
     datasets.clear();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots a function of time on interval [t0,tmax]");
     myControl.println ("where NPTS is the number of points and dt=(tmax-t0)/NPTS");
     t0=0;
     tmax=4;
     NPTS = 100;
     dt=(tmax-t0)/NPTS;
     myControl.setValue("f(x=time)", "2*cos(2*pi*x/2)");
     myControl.setValue("t0",t0);
     myControl.setValue("tmax",tmax);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("dt([tmax-t0]/NPTS)=",dt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new xoftApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     //myControl.setLocation(520, 30);
     //myControl.setSize(250,505);
     myControl.setLocation(270, 30);
     myControl.setSize(430,505);
     model.setControl(myControl);
   }
}
