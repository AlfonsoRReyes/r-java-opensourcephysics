/*
@Author J E Hasbun 2007.
Plots the vector function (-y,x,0) and its curl.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.display2d.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.frames.Display3DFrame;


public class curlApp implements Calculation {
  Display3DFrame frame3Da = new Display3DFrame("Vectors Application");
  Display3DFrame frame3Db = new Display3DFrame("Vectors Application");
  ElementArrow [][][][] arrow;
  ElementCircle point=new ElementCircle();
  private Control myControl;
  double [][][][] V;
  double [][][][] Curl;
  double dx,dy,dz,dv,vmax,xmax,ymax,zmax,vmin,xmin,ymin,zmin,dmin,dmax;
  int NPTS, NPTS_old;

  public curlApp(){
    frame3Da.setLocation(5,5);
    frame3Da.setSize(350, 295);
    frame3Da.setDecorationType(VisualizationHints.DECORATION_AXES);
    frame3Db.setLocation(5,300);
    frame3Db.setSize(350, 295);
    frame3Db.setDecorationType(VisualizationHints.DECORATION_AXES);
    point = new ElementCircle(); //define an interactive point
  }

  public void calculate() {
    clear();
    cleararrows(NPTS_old); //clear old arrows before a new NPTS is defined
    //SurfacePlot plot = new SurfacePlot();
    MultiVarFunction Vxyz_x,Vxyz_y,Vxyz_z,Cxyz_x,Cxyz_y,Cxyz_z;
    String [] vars = {"x","y","z"};
    String Vx = myControl.getString("V_x"); // the function
    String Vy = myControl.getString("V_y"); // the function
    String Vz = myControl.getString("V_z"); // the function
    String Cx = myControl.getString("Curl_x"); // the function
    String Cy = myControl.getString("Curl_y"); // the function
    String Cz = myControl.getString("Curl_z"); // the function
    try { // read in functions
    Vxyz_x = new ParsedMultiVarFunction(Vx,vars);
    Vxyz_y = new ParsedMultiVarFunction(Vy,vars);
    Vxyz_z = new ParsedMultiVarFunction(Vz,vars);
    Cxyz_x = new ParsedMultiVarFunction(Cx,vars);
    Cxyz_y = new ParsedMultiVarFunction(Cy,vars);
    Cxyz_z = new ParsedMultiVarFunction(Cz,vars);
    } catch(ParserException ex) {
    myControl.println(ex.getMessage());
    return;
    }
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
    double [] xyz = new double[3];
    //vector functions plots
    Curl=new double[3][NPTS][NPTS][NPTS];// used for surface plot
    V=new double [3][NPTS][NPTS][NPTS];// used for function vector plot
    for(int i = 0; i < NPTS; i++) {
      xyz[0]= xmin + i * dx;
      for(int j = 0; j < NPTS; j++) {
        xyz[1]= ymin + j * dy;
          for(int k = 0; k < NPTS; k++) {
            xyz[2] = zmin + k * dz;
            V[0][i][j][k] = Vxyz_x.evaluate(xyz);
            V[1][i][j][k] = Vxyz_y.evaluate(xyz);
            V[2][i][j][k] = Vxyz_z.evaluate(xyz);
            Curl[0][i][j][k] = Cxyz_x.evaluate(xyz);
            Curl[1][i][j][k] = Cxyz_y.evaluate(xyz);
            Curl[2][i][j][k] = Cxyz_z.evaluate(xyz);
            arrow[0][i][j][k].setXYZ(xyz[0], xyz[1], xyz[2]);//x,y,z arrows positions
            arrow[1][i][j][k].setXYZ(xyz[0], xyz[1], xyz[2]);
            //arrows lengths according vector functions x,y,z components divided by NPTS
            arrow[0][i][j][k].setSizeXYZ(V[0][i][j][k] / NPTS,V[1][i][j][k] / NPTS,
                                      V[2][i][j][k] / NPTS);
            arrow[1][i][j][k].setSizeXYZ(Curl[0][i][j][k] / NPTS,Curl[1][i][j][k] / NPTS,
                                      Curl[2][i][j][k] / NPTS);
          }
      }
    }
    for(int i = 0; i < NPTS; i++) {
      for(int j = 0; j < NPTS; j++) {
        for(int k = 0; k < NPTS; k++) {
          frame3Da.addElement(arrow[0][i][j][k]);
          frame3Db.addElement(arrow[1][i][j][k]);
        }
      }
    }
    point.setXYZ(0,0,0);                //place the interactive point at the origin
    point.setSizeXYZ(xmax*(1-0.95),ymax*(1-0.95),zmax*(1-0.95)); //point size
    point.getStyle().setFillColor(java.awt.Color.GREEN);
    point.getInteractionTarget(Element.TARGET_POSITION).setEnabled(true);
    frame3Da.addElement(point);         //Interactive point at the origin for V
    frame3Da.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
    frame3Da.setAltitude(0.5);
    frame3Da.setAzimuth(0.5);
    frame3Da.setTitle("Vector Function V=("+Vx+","+Vy+","+Vz+")");
    frame3Db.addElement(point);         //Interactive point at the origin for Curl
    frame3Db.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
    frame3Db.setAltitude(0.5);
    frame3Db.setAzimuth(0.5);
    frame3Db.setTitle("Curl Function Curl=("+Cx+","+Cy+","+Cz+")");
    myControl.println ("Upper left graph is a 3D representation (arrows) of the vector function.");
    myControl.println ("The Lower graph is its curl representation. Rotate the 3d plots or expand");
    myControl.println ("the windows for better viewing. The green interactive point can be moved.");
    myControl.println("The Vector coordinates have been divided by NPTS for better viewing.");
    myControl.println("Shift-drag the mouse within the surface plot panel to magnify the plot.");
  }

  public void initialarrows(int NPTS) {
    //initialize the arrow array
    arrow=new ElementArrow[2][NPTS][NPTS][NPTS];
    for (int m = 0; m < 2; m++){
      for (int i = 0; i < NPTS; i++) {
        for (int j = 0; j < NPTS; j++) {
          for (int k = 0; k < NPTS; k++) {
            arrow[m][i][j][k] = new ElementArrow();
          }
        }
      }
    }
  }

  public void cleararrows(int NPTS) {
    //reset the arrows
    for(int m = 0; m < 2; m++) {
      for (int i = 0; i < NPTS; i++) {
        for (int j = 0; j < NPTS; j++) {
          for (int k = 0; k < NPTS; k++) {
            arrow[m][i][j][k].setXYZ(0, 0, 0); //x,y,z arrow position
            //arrow length corresponds to the vector function x,y,z components
            arrow[m][i][j][k].setSizeXYZ(0, 0, 0);
          }
        }
      }
    }
   //clear the arrows from the frame
      for (int i = 0; i < NPTS; i++) {
        for (int j = 0; j < NPTS; j++) {
          for (int k = 0; k < NPTS; k++) {
            frame3Da.addElement(arrow[0][i][j][k]);
            frame3Da.addElement(arrow[1][i][j][k]);
          }
        }
      }
   frame3Da.repaint();
   frame3Db.repaint();
  }

  public void clear() {
    myControl.clearMessages();
    V=null;
    Curl=null;
    point.setXYZ(0, 0, 0);
    point.setSizeXYZ(0, 0, 0);
    point.setXYZ(0, 0, 0);
    point.setSizeXYZ(0, 0, 0);
    frame3Da.addElement(point);
    frame3Db.addElement(point);
  }

  public void resetCalculation() {
    clear();
    if(NPTS_old !=0){cleararrows(NPTS_old);}//NPTS_old aquires a value later
    myControl.println ("Plots the vector function (-y,x,0) and its curl.");
    myControl.println ("The inputs can be changed as needed.");
    vmax=1.0; vmin=-vmax;
    xmax=vmax; ymax=vmax; zmax=vmax;
    xmin=vmin; ymin=vmin; zmin=vmin;
    NPTS = 5; NPTS_old=NPTS;
    initialarrows(NPTS); //initial arrows with the NPTS
    myControl.setValue("V_x", "-y"); //x-component of vector function
    myControl.setValue("V_y", "x"); // x-component of vector function
    myControl.setValue("V_z", "0"); // x-component of vector function
    myControl.setValue("Curl_x", "0"); // x-component of the Curl
    myControl.setValue("Curl_y", "0"); // x-component of the Curl
    myControl.setValue("Curl_z", "2"); // x-component of the Curl
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
    Calculation model = new curlApp();
    CalculationControl myControl;
    myControl = new CalculationControl(model);
    myControl.setLocation(360, 5);
    myControl.setSize(420,590);
    myControl.setDividerLocation(350);
    model.setControl(myControl);
  }
}
