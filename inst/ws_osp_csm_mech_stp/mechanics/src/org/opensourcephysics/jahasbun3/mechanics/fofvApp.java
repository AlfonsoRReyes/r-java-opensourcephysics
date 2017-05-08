/*
@Author J E Hasbun 2007.
Calculation of position, velocity, and acceleration for a body in
free fall with air resistance versus time. The air resistance term
used is proportional to velocity. The analytical solution is used.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.Function;
import java.awt.Color;

public class fofvApp implements Calculation {
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
   double dt, t0, tmax, g, m, C, y0, v0, Vt;
   int NPTS;
   int maxiter=25;
   double tolerance=5e-3;
   Bisnewt bisnewt=new Bisnewt();
   Function ft=new fyoft();

   public fofvApp(){
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
     g=myControl.getDouble("g");
     m = myControl.getDouble("m");
     C = myControl.getDouble("C");
     if(C <=0){
     C=1.e-9;
     myControl.setValue("C", C);}
     y0 = myControl.getDouble("y0");
     v0 = myControl.getDouble("v0");
     //tmax = myControl.getDouble("tmax");
     //let landing time interval between 0 and free fall round trip time
     double x_initial = 0.0;
     double x_final   = v0/g+Math.sqrt(Math.pow(v0/g,2)+2*y0/g);
     //solve for the landing time iteratively (see function fyoft)
     tmax=bisnewt.Bisnewt(x_initial, x_final, maxiter, tolerance,ft);
     myControl.setValue("tmax", tmax);
     NPTS = myControl.getInt("NPTS");
     dt=tmax/NPTS;
     myControl.setValue("dt(tmax/NPTS)=",dt);
     Vt=Math.abs(m*g/C);
     myControl.setValue("Vt",Vt);
     t0=0.0;
     tt = new double[NPTS+1];
     yy = new double[NPTS+1];
     vv = new double[NPTS+1];
     aa = new double[NPTS+1];
     for(int i = 0; i <= NPTS; i++) {
       tt[i]=t0+i*dt;               //time
       vv[i]=(m*g/C+v0)*Math.exp(-C*tt[i]/m)-m*g/C;           //velocity
       yy[i]=y0-m*(g*tt[i]+(m*g/C+v0)*(Math.exp(-C*tt[i]/m)-1))/C;//position
       aa[i]=-g-C*vv[i]/m;                                    //acceleration
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
     myControl.println ("Falling with gravity/air resistance - Analytic Method");
     myControl.println ("F=-m*g-C*v=mdv/dt, so that v(t)=(m*g/C+v0)*exp(-C*t/m)-m*g/C");
     myControl.println ("y(t)=y0-m*(g*t+(m*g/C+v0)*(exp(-C*t/m)-1))/C, a=F/m");
     myControl.println ("g is the gravitational acceleration");
     myControl.println ("m is the particle's mass");
     myControl.println ("C is the drag coefficient");
     myControl.println ("y0 is the initial position");
     myControl.println ("v0 is the initial velocity");
     myControl.println ("NPTS is the number of steps");
     myControl.println ("tmax is the time at which the object returns to the gound,");
     myControl.println ("obtained by finding the time when y(t=tmax)=0");
     myControl.println ("dt=tmax/NPTS is the time step used");
     myControl.println ("terminal speed Vt=(m*g/C)");
     t0=0;
     g = 9.8;
     m = 1.0;
     C = 0.05;
     y0 = 10.0;
     v0 = 20.0;
     tmax = 5.0;
     Vt=Math.abs(m*g/C);
     NPTS = 200;
     dt=tmax/NPTS;
     yy=null;
     vv=null;
     aa=null;
     tt=null;
     myControl.setValue("g", g);
     myControl.setValue("m", m);
     myControl.setValue("C", C);
     myControl.setValue("y0", y0);
     myControl.setValue("v0", v0);
     myControl.setValue("NPTS", NPTS);
     myControl.setValue("tmax", tmax);
     myControl.setValue("dt(tmax/NPTS)=",dt);
     myControl.setValue("Vt",Vt);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public class fyoft implements Function {
     public double evaluate(double t) {
       //defines y(t=tmax)=0 on landing, here x=time
       double f;
       f = y0 - m * (g * t + (m * g / C + v0) * (Math.exp( -C * t / m) - 1)) / C;
       return f;
     }
   }

   public static void main(String[] args) {
     Calculation model = new fofvApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(520, 30);
     myControl.setSize(250,505);
     model.setControl(myControl);
   }
}