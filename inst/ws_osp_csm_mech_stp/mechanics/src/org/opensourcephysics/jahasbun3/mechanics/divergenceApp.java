/*
@Author J E Hasbun 2007.
Plots the vector function (x,y,z) and does a surface plot of its divergence
evaluated at a particular value of z.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display2d.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.frames.Display3DFrame;

public class divergenceApp implements Calculation {
  Display3DFrame frame3D = new Display3DFrame("Vectors Application");
  DrawingPanel drawingPanel = new DrawingPanel();
  DrawingFrame frame = new DrawingFrame(drawingPanel);
  ElementArrow [][][] arrow;
  ElementCircle point=new ElementCircle();
  private Control myControl;
  double [][][][] F;
  double [][] Div;
  double dx,dy,dz,dv,vmax,xmax,ymax,z,zmax,vmin,xmin,ymin,zmin,dmin,dmax;
  int NPTS, NPTS_old;

  public divergenceApp(){
    frame.setCartesian("x","y","z");
    frame.setLocation(5,300);
    frame.setSize(350, 295);
    frame3D.setLocation(5,5);
    frame3D.setSize(350, 295);
    frame3D.setDecorationType(VisualizationHints.DECORATION_AXES);
    point = new ElementCircle(); //define an interactive point
  }

  public void calculate() {
    clear();
    cleararrows(NPTS_old); //clear old arrows before a new NPTS is defined
    SurfacePlot plot = new SurfacePlot();
    MultiVarFunction Fxyz_x,Fxyz_y,Fxyz_z, Div_xyz;
    String [] vars = {"x","y","z"};
    String Fx = myControl.getString("F_x"); // the function
    String Fy = myControl.getString("F_y"); // the function
    String Fz = myControl.getString("F_z"); // the function
    String div= myControl.getString("Div_at_zv"); // x-component of function gradient

    try { // read in function
    Fxyz_x = new ParsedMultiVarFunction(Fx,vars);
    Fxyz_y = new ParsedMultiVarFunction(Fy,vars);
    Fxyz_z = new ParsedMultiVarFunction(Fz,vars);
    Div_xyz = new ParsedMultiVarFunction(div,vars);
    } catch(ParserException ex) {
    myControl.println(ex.getMessage());
    return;
    }
    z=myControl.getDouble("zv");
    xmin=myControl.getDouble("xmin");
    xmax=myControl.getDouble("xmax");
    NPTS=myControl.getInt("NPTS");
    NPTS_old=NPTS;
    initialarrows(NPTS);//redefine new arrows
    dx=(xmax-xmin)/NPTS;
    myControl.setValue("dx((xmax-xmin)/NPTS)=",dx);
    ymin=myControl.getDouble("ymin");
    ymax=myControl.getDouble("ymax");
    dy=(ymax-ymin)/NPTS;
    myControl.setValue("dy((ymax-ymin)/NPTS)=",dy);
    zmin=myControl.getDouble("zmin");
    zmax=myControl.getDouble("zmax");
    dz=(zmax-zmin)/NPTS;
    myControl.setValue("dz((zmax-zmin)/NPTS)=",dz);
    //initialize Div
    dmin=0; dmax=0;
    Div=new double[NPTS][NPTS];// used for surface plot
    double [] xyz = new double[3];
    plot.setAll(Div,xmin, xmax, ymin, ymax);//2w-surface in plot, & sets the data + scale - needed
    xyz[2]=z;
    for(int i = 0; i < NPTS; i++) {
      double x0 = plot.indexToX(i); //x variable for surface - needed
      xyz[0]=x0;
        for(int j = 0; j< NPTS; j++) {
          double y0 = plot.indexToY(j); //y variable for surface - needed
          xyz[1]=y0;
          Div[i][j]=Div_xyz.evaluate(xyz);
          if(Div[i][j]<dmin){dmin=Div[i][j];}
          if(Div[i][j]>dmax){dmax=Div[i][j];}
        }
    }
    //The divergence surface plot
    plot.setAll(Div, xmin, xmax, ymin, ymax);//function data with coordinates
    plot.setAutoscaleZ(false,dmin*(1-0.2),dmax*(1+0.2));
    frame.setTitle("Vector Divergence="+div);
    drawingPanel.addDrawable(plot);
    SurfacePlotMouseController mouseController = new SurfacePlotMouseController(drawingPanel, plot);
    drawingPanel.addMouseListener(mouseController);
    drawingPanel.addMouseMotionListener(mouseController);
    frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    //next is the vector function plot
    F=new double [3][NPTS][NPTS][NPTS];// used for function vector plot
    for(int i = 0; i < NPTS; i++) {
      xyz[0]= xmin + i * dx;
      for(int j = 0; j < NPTS; j++) {
        xyz[1]= ymin + j * dy;
          for(int k = 0; k < NPTS; k++) {
            xyz[2] = zmin + k * dz;
            F[0][i][j][k] = Fxyz_x.evaluate(xyz);
            F[1][i][j][k] = Fxyz_y.evaluate(xyz);
            F[2][i][j][k] = Fxyz_z.evaluate(xyz);
            arrow[i][j][k].setXYZ(xyz[0], xyz[1], xyz[2]);//x,y,z arrow position
            //arrow length according vector function x,y,z components divided by NPTS
            arrow[i][j][k].setSizeXYZ(F[0][i][j][k] / NPTS,F[1][i][j][k] / NPTS,
                                      F[2][i][j][k] / NPTS);
          }
      }
    }
    for(int i = 0; i < NPTS; i++) {
      for(int j = 0; j < NPTS; j++) {
        for(int k = 0; k < NPTS; k++) {
          frame3D.addElement(arrow[i][j][k]);
        }
      }
    }
    point.setXYZ(0,0,0);                //place the interactive point at the origin
    point.setSizeXYZ(xmax*(1-0.95),ymax*(1-0.95),zmax*(1-0.95)); //point size
    point.getStyle().setFillColor(java.awt.Color.GREEN);
    point.getInteractionTarget(Element.TARGET_POSITION).setEnabled(true);
    frame3D.addElement(point);          //Interactive point at the origin
    frame3D.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
    frame3D.setAltitude(0.5);
    frame3D.setAzimuth(0.5);
    frame3D.setTitle("Vector Function F=("+Fx+","+Fy+","+Fz+")");
    myControl.println ("Upper left graph is a 3D representation (arrows) of the vector function.");
    myControl.println ("The Lower left is the x,y components of its divergence versus");
    myControl.println ("x,y at a value of z=zv. The inputs can be changed as needed.");
    myControl.println ("Rotate the 3d plots or expand the windows for better viewing.");
    myControl.println("The green interactive point can be moved. The Vector coordinates");
    myControl.println("have been divided by NPTS for better viewing.");
    myControl.println("Shift-drag the mouse within the surface plot panel to magnify the plot.");
  }

  public void initialarrows(int NPTS) {
    //initialize the arrow array
    arrow=new ElementArrow[NPTS][NPTS][NPTS];
    for(int i = 0; i < NPTS; i++) {
      for(int j= 0; j < NPTS; j++) {
        for(int k=0; k < NPTS; k++){
          arrow[i][j][k] = new ElementArrow();
        }
      }
    }
  }

  public void cleararrows(int NPTS) {
    //reset the arrows
    for(int i = 0; i < NPTS; i++) {
     for(int j = 0; j < NPTS; j++) {
         for(int k = 0; k < NPTS; k++) {
           arrow[i][j][k].setXYZ(0, 0, 0);//x,y,z arrow position
           //arrow length corresponds to the vector function x,y,z components
           arrow[i][j][k].setSizeXYZ(0,0,0);
         }
     }
  }
   //clear the arrows from the frame
   for(int i = 0; i < NPTS; i++) {
     for(int j = 0; j < NPTS; j++) {
       for(int k = 0; k < NPTS; k++) {
         frame3D.addElement(arrow[i][j][k]);
       }
     }
   }
   frame3D.repaint();
  }

  public void clear() {
    myControl.clearMessages();
    frame.clearDrawables();
    drawingPanel.clear();
    F=null;
    Div=null;
    point.setXYZ(0, 0, 0);
    point.setSizeXYZ(0, 0, 0);
    point.setXYZ(0, 0, 0);
    point.setSizeXYZ(0, 0, 0);
    frame3D.addElement(point);
  }

  public void resetCalculation() {
    clear();
    if(NPTS_old !=0){cleararrows(NPTS_old);}//NPTS_old aquires a value later
    myControl.println ("Plots the vector function (x,y,z) does a surface plot");
    myControl.println ("of its divergence evaluated at a particular value of z.");
    myControl.println ("The divergence of the (x,y,z) vector function is 3.");
    vmax=1.0; vmin=-vmax;
    xmax=vmax; ymax=vmax; zmax=vmax;
    xmin=vmin; ymin=vmin; zmin=vmin;
    NPTS = 5; NPTS_old=NPTS;
    initialarrows(NPTS); //initial arrows with the NPTS
    z=0.4;
    myControl.setValue("F_x", "x"); // x-component of vector function
    myControl.setValue("F_y", "y"); // x-component of vector function
    myControl.setValue("F_z", "z"); // x-component of vector function
    myControl.setValue("Div_at_zv", "3");//divergence of the vector
    myControl.setValue("zv",z);
    myControl.setValue("xmin", xmin);
    myControl.setValue("xmax", xmax);
    myControl.setValue("NPTS", NPTS);
    dx=(xmax-xmin)/NPTS;
    myControl.setValue("dx((xmax-xmin)/NPTS)=",dx);
    myControl.setValue("ymin", ymin);
    myControl.setValue("ymax", ymax);
    dy=(ymax-ymin)/NPTS;
    myControl.setValue("dy((ymax-ymin)/NPTS)=",dy);
    myControl.setValue("zmin", zmin);
    myControl.setValue("zmax", zmax);
    dz=(zmax-zmin)/NPTS;
    myControl.setValue("dz((zmax-zmin)/NPTS)=",dz);
  }

  public void setControl(Control control) {
    myControl = control;
    resetCalculation();
  }

  public static void main(String[] args) {
    Calculation model = new divergenceApp();
    CalculationControl myControl;
    myControl = new CalculationControl(model);
    myControl.setLocation(360, 5);
    myControl.setSize(420,590);
    myControl.setDividerLocation(350);
    model.setControl(myControl);
  }
}
