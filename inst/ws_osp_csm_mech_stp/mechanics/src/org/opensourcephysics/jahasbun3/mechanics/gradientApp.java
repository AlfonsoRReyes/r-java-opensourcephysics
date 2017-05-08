/*
@Author J E Hasbun 2007.
Plots the function x*exp(-x^2-y^2-z^2), its x, y gradient
components evaluated at a particular value of z. Also does a contour plot.
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
import org.opensourcephysics.frames.Vector2DFrame;

public class gradientApp implements Calculation {
  DrawingPanel drawingPanela = new DrawingPanel();
  DrawingPanel drawingPanelb = new DrawingPanel();
  DrawingFrame framea = new DrawingFrame(drawingPanela);
  DrawingFrame frameb = new DrawingFrame(drawingPanelb);
  Vector2DFrame framec = new Vector2DFrame("x", "y", "");
  ContourPlot plotb = new ContourPlot();
  private Control myControl;
  double [][] F;
  double [][][] Grad_xy;
  double dx,dy,dz,dv,vmax,xmax,ymax,z,zmax,vmin,xmin,ymin,zmin;
  int NPTS,m;

  public gradientApp(){
    framea.setLocation(5,5);
    framea.setSize(350, 290);
    frameb.setLocation(360,400);
    frameb.setSize(310,190);
    framec.setLocation(5,300);
    framec.setSize(350, 290);
  }

  public void calculate() {
    clear();
    SurfacePlot plota = new SurfacePlot();
    MultiVarFunction Fxyz, Gxyz_x, Gxyz_y;
    String [] vars = {"x","y","z"};
    String fxyz = myControl.getString("F(x,y,z)"); // the function
    String gxyz_x= myControl.getString("Grad_x"); // x-component of function gradient
    String gxyz_y= myControl.getString("Grad_y"); // y-component of function gradient
    try { // read in function
    Fxyz = new ParsedMultiVarFunction(fxyz,vars);
    Gxyz_x = new ParsedMultiVarFunction(gxyz_x,vars);
    Gxyz_y = new ParsedMultiVarFunction(gxyz_y,vars);
    } catch(ParserException ex) {
    myControl.println(ex.getMessage());
    return;
    }
    z=myControl.getDouble("z");
    xmin=myControl.getDouble("xmin");
    xmax=myControl.getDouble("xmax");
    NPTS=myControl.getInt("NPTS");
    dx=(xmax-xmin)/NPTS;
    myControl.setValue("dx((xmax-xmin)/NPTS)=",dx);
    ymin=myControl.getDouble("ymin");
    ymax=myControl.getDouble("ymax");
    dy=(ymax-ymin)/NPTS;
    myControl.setValue("dy((ymax-ymin)/NPTS)=",dy);
    F=  new double [NPTS][NPTS];       // used for surface plot
    Grad_xy=new double [2][NPTS][NPTS];// used for gradient plot
    //** two ways to get his to work: 1w and 2w -- see comment lines
    //double [] xyz = new double[3];      //1w
    //xyz[2]=z;                           //1w-z variable - fixed
    double [][] xyz = new double[2][3];   //2w
    xyz[0][2]=z;                          //2w-z variable - fixed
    xyz[1][2]=z;                          //2w-z variable - fixed
    plota.setAll(F,xmin, xmax, ymin, ymax);//2w-surface in plota, & sets the data + scale - needed
    framec.setAll(Grad_xy); //2w vector field in framec - needed
    for(int i = 0; i < NPTS; i++) {
      //xyz[0]=xmin+i*dx;            //1w-x variable - not used
      double x0 = plota.indexToX(i); //2w-x variable for surface - needed
      double x1 = framec.indexToX(i);//2w-x variable for the frame - needed
      xyz[0][0]=x0;                  //2w
      xyz[1][0]=x1;                  //2w
      for(int j = 0; j< NPTS; j++) {
        //xyz[1]=-(ymin+j*dy);  //1w-y variable - notice important sign change
        //F[i][j]=Fxyz.evaluate(xyz);           //1w
        //Grad_xy[0][i][j]=Gxyz_x.evaluate(xyz);//1w
        //Grad_xy[1][i][j]=Gxyz_y.evaluate(xyz);//1w
        double y0 = plota.indexToY(j); //2w-y variable for surface - needed
        double y1 = framec.indexToY(j);//2w-y variable for the frame - needed
        xyz[0][1]=y0;                             //2w
        xyz[1][1]=y1;                             //2w
        F[i][j]=Fxyz.evaluate(xyz[0]);            //2w
        Grad_xy[0][i][j]=Gxyz_x.evaluate(xyz[1]); //2w
        Grad_xy[1][i][j]=Gxyz_y.evaluate(xyz[1]); //2w
        //myControl.println("i,j->"+xyz[0]+", "+xyz[1]+" F="+F[i][j]);
        //myControl.println("x,y->"+xyz[0]+", "+xyz[1]+" G_x="+Grad_xy[0][i][j]+", G_y="+
        //                 Grad_xy[1][i][j]);
      }
    }
    myControl.println ("Upper left graph contains the 3D function plot versus x,y evaluated");
    myControl.println ("at a particular value of z. The Lower left is the x,y components");
    myControl.println ("of its gradient (arrows) versus x,y at a value of z. The contour plot is at the");
    myControl.println ("bottom. Rotate the 3d plot or expand the 2d windows for better viewing.");
    myControl.println ("The contour color plot legend is behind the contour plot.");
    myControl.println("Shift-drag the mouse within the surface plot panel to magnify the plot.");
    plota.setAll(F, xmin, xmax, ymin, ymax);//function data with coordinates
    framea.setTitle("Function F="+fxyz);
    drawingPanela.addDrawable(plota);
    SurfacePlotMouseController mouseController = new SurfacePlotMouseController(drawingPanela, plota);
    drawingPanela.addMouseListener(mouseController);
    drawingPanela.addMouseMotionListener(mouseController);
    framea.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    // next is the gradient x,y vector components
    framec.setPreferredMinMax(xmin, xmax, ymin, ymax);
    framec.setTitle("Grad_x(F)="+gxyz_x+", Grad_y(F)="+gxyz_y);
    framec.setAll(Grad_xy);
    framec.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    plotb.setAll(F, xmin, xmax, ymin, ymax);//contour data with coordinates
    frameb.setTitle("Function F="+fxyz+" contours");
    drawingPanelb.addDrawable(plotb);       //contour
    frameb.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    plotb.showLegend().setTitle("Contour Legend");
    plotb.showLegend().setLocation(475,400);
    plotb.showLegend().setResizable(true);
    //plotb.showLegend().setSize(240,150);
  }

  public void clear  () {
    myControl.clearMessages();
    framea.clearDrawables();
    frameb.clearDrawables();
    framec.clearDrawables();
    drawingPanela.clear();
    drawingPanelb.clear();
    F=null;
    Grad_xy=null;
  }

  public void resetCalculation() {
    clear();
    myControl.println ("Plots the function x*exp(-x^2-y^2-z^2) and its x, y gradient,");
    myControl.println ("components evaluated at a particular value of z.");
    myControl.println ("The x, y gradient components of this function are:");
    myControl.println ("Grad_x(f)=(1-2*x^2)*exp(-x^2-y^2-z^2)");
    myControl.println ("Grad_y(f)=-2*x*y*exp(-x^2-y^2-z^2)");
    myControl.println ("The inputs can be changed as needed.");
    vmax=2.0; vmin=-vmax;
    xmax=vmax; ymax=vmax; zmax=vmax;
    xmin=vmin; ymin=vmin; zmin=vmin;
    NPTS = 40;
    dx=(xmax-xmin)/NPTS;
    dy=(ymax-ymin)/NPTS;
    dz=dx;
    m=Math.round(NPTS/2+5);
    z=zmin+(m-1)*dz;
    myControl.setValue("F(x,y,z)", "x*exp(-x^2-y^2-z^2)"); // the function
    myControl.setValue("Grad_x", "(1-2*x^2)*exp(-x^2-y^2-z^2)"); // x-component of function gradient
    myControl.setValue("Grad_y", "-2*x*y*exp(-x^2-y^2-z^2)");  // y-component of function gradient
    myControl.setValue("z",z);
    myControl.setValue("xmin", xmin);
    myControl.setValue("xmax", xmax);
    myControl.setValue("NPTS", NPTS);
    myControl.setValue("dx((xmax-xmin)/NPTS)=",dx);
    myControl.setValue("ymin", ymin);
    myControl.setValue("ymax", ymax);
    myControl.setValue("dy((ymax-ymin)/NPTS)=",dy);
  }

  public void setControl(Control control) {
    myControl = control;
    resetCalculation();
  }

  public static void main(String[] args) {
    Calculation model = new gradientApp();
    CalculationControl myControl;
    myControl = new CalculationControl(model);
    myControl.setLocation(360, 5);
    myControl.setSize(420,400);
    myControl.setDividerLocation(250);
    model.setControl(myControl);
  }
}
