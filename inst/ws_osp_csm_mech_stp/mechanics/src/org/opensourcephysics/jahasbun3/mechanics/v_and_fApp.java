/*
@Author J E Hasbun 2007.
Plots the molecule model potential V(x)=(1/x^3-1/x^2), and related force.
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

public class v_and_fApp implements Calculation {
   PlottingPanel plotPanelx    = new PlottingPanel ("x axis ", "V(x), F(x)",
                                "V(x)-black, F(x)-red");
   DrawingFrame  drawingFramex = new DrawingFrame ("Plot", plotPanelx);
   DataTable dataTable = new DataTable();
   DataTableFrame framea= new DataTableFrame("Calculations", dataTable);
   DatasetManager datasets = new DatasetManager();
   Dataset dataseta = new Dataset ();
   Dataset datasetb = new Dataset ();
   Dataset datasetc = new Dataset ();
   Dataset datasetd = new Dataset ();
   private Control myControl;
   double [] V, F, xx;
   double dx,x0, xmax,Vmin,Xb;
   int NPTS;

   public v_and_fApp(){
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
     datasets.setXYColumnNames(0, "x", "V(x)");
     datasets.setSorted(true);
     datasets.setXYColumnNames(1, "", "F(x)");
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
     Function Fx, Vx;
     dataseta.clear();
     datasets.clear();
     dataTable.clear();
     datasetb.clear();
     plotPanelx.clear();
     String vx = myControl.getString("V(x)");
     String fx = myControl.getString("F(x)=-dV(x)/dx");
     try { // read in function
     Vx = new ParsedFunction(vx);
     Fx = new ParsedFunction(fx);
   } catch(ParserException ex) {
     myControl.println(ex.getMessage());
     return;
   }
     x0=myControl.getDouble("x0");
     xmax=myControl.getDouble("xmax");
     NPTS=myControl.getInt("NPTS");
     dx=(xmax-x0)/NPTS;
     myControl.setValue("dx((xmax-x0)/NPTS)=",dx);
     xx = new double[NPTS+1];
     V = new double[NPTS+1];
     F = new double[NPTS+1];
     Vmin=1000;
     for(int i = 0; i <= NPTS; i++) {
       xx[i]=x0+i*dx;            //position
       V[i]=Vx.evaluate(xx[i]); //Potential Energy
       F[i]=Fx.evaluate(xx[i]); //The force
       //System.out.println("x="+xx[i]+" V="+V[i]+" F="+F[i]);
       if (V[i]<=Vmin){
         Vmin = V[i];
         Xb=xx[i];
       }
     }
     myControl.setValue("Vmin",Vmin);
     myControl.setValue("Xb",Xb);
     dataseta.append(xx, V);
     datasets.append(0, xx, V);
     datasetb.append(xx, F);
     datasets.append(1, xx, F);
     dataTable.add(datasets);
     dataTable.refreshTable();
     plotPanelx.addDrawable(dataseta);
     plotPanelx.addDrawable(datasetb);
     plotPanelx.render();
   }

   public void clear  () {
     plotPanelx.clear();
     plotPanelx.render();
     dataTable.clear();
     dataTable.add(datasets);
     dataTable.refreshTable();
     dataseta.clear();
     datasets.clear();
     datasetb.clear();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Molecule Potential V=1/x^3-1/x^2, and force F=3/x^4+2/x^3");
     myControl.println ("where F=-dV/dx, and the minimum of the potential Vmin");
     myControl.println ("corresponds to the minimum value of V(x) at x=Xb, where");
     myControl.println ("Xb is the bond length, at which point the force takes");
     myControl.println ("a value of zero. Vmin, and Xb are found at calculation time.");
     myControl.println ("If a different V(x) is entered, also enter the appropriate");
     myControl.println ("F(x). Increasing NPTS gives more accurate values");
     myControl.println ("of Xb and Vmin");
     Xb=0.0;
     Vmin=0.0;
     x0=1.0;
     xmax=4.0;
     NPTS = 100;
     dx=(xmax-x0)/NPTS;
     xx=null;
     myControl.setValue("V(x)", "1/x^3-1/x^2"); // the potential
     myControl.setValue("F(x)=-dV(x)/dx", "3/x^4-2/x^3"); // the force
     myControl.setValue("x0", x0);
     myControl.setValue("xmax", xmax);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("dx((xmax-x0)/NPTS)=",dx);
     myControl.setValue("Xb",Xb);
     myControl.setValue("Vmin",Vmin);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new v_and_fApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(370, 30);
     myControl.setSize(330,505);
     model.setControl(myControl);
   }
}
