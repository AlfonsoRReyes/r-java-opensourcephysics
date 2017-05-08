/*
@Author J E Hasbun 2007.
Shows Euler angles: ph, th, ps; the planes and the line of nodes.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.ejs.control.*;
import org.opensourcephysics.display3d.simple3d.*;
//import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import java.awt.Color;
import java.text.*;
//import java.util.ArrayList; can this be used to list the ejs properties with getPropertyList()?


public class euler_angApp implements Calculation {
  Display3DFrame frame3d = new Display3DFrame("Eulerian Angles, phi, theta, psi");
  //The drawpanel can be used with the slider if needed
  //DrawingPanel drawpanel = new PlottingPanel("x", "y", "EJS Controls Example");
  //Dataset dataset=new Dataset(); //can use with drawpanel and slider
  private Control myControl;
  double r[][], rp[][], rpp[][], rppp[][];
  double rphi[][], rthe[][], rpsi[][], rphie[][], rthee[][];
  double pi, rc, ph, th, ps, phi, the, psi;
  double x=1, y=1, z=1, xmin=-x, xmax=x, ymin=-y, ymax=y, zmin=-z, zmax=z;
  int N=3, Iarrows=12, Iflag;
  ElementArrow [] arrow=new ElementArrow [Iarrows];
  ElementEllipsoid [] ellipso=new ElementEllipsoid[2];
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");
  Matrix3DTransformation m3dt;
  GroupControl control;

  public euler_angApp(){
    frame3d.setDecorationType(VisualizationHints.CURSOR_NONE);
    frame3d.setLocation(5,5);
    frame3d.setSize(480,500);
    frame3d.setVisible(true);
    pi=Math.PI; rc=pi/180.;
    nf.setMaximumFractionDigits(4);
    for(int i=0; i< Iarrows; i++){
      arrow[i]=new ElementArrow();
      arrow[i].getStyle().setLineWidth(1.5f);
    }
    arrow[0].getStyle().setLineColor(Color.blue);
    arrow[1].getStyle().setLineColor(Color.blue);
    arrow[2].getStyle().setLineColor(Color.blue);
    arrow[3].getStyle().setLineColor(Color.black);
    arrow[4].getStyle().setLineColor(Color.black);
    arrow[5].getStyle().setLineColor(Color.black);
    arrow[6].getStyle().setLineColor(Color.red);
    arrow[7].getStyle().setLineColor(Color.red);
    arrow[8].getStyle().setLineColor(Color.red);
    arrow[9].getStyle().setLineColor(Color.green);
    arrow[10].getStyle().setLineColor(Color.green);
    arrow[11].getStyle().setLineColor(Color.green);
    ellipso[0]=new ElementEllipsoid();
    ellipso[1]=new ElementEllipsoid();

    //sliders follow

    //there are two similar blocks below
    //The first large group that's commented works, but I needed to place the
    //sliders into frame 3d, so worked with the block that follows this
    //so a 2nd block was created.

    //1st large block - commented
/*
    //the commented line below can be used with the commented drawpanel above
    Iflag=0; //flag used in sliderMoved()
    control = new GroupControl(this);
    control.add("Frame","name=plottingFrame; title=phi, theta, psi slider controls;exit=true; size=475,50; location=5,505");
    //control.addObject(drawpanel, "Panel", "name=drawpanel; parent=plottingFrame; position=center");
    control.add("Panel", "name=controlPanel; parent=plottingFrame; position=south; layout=hbox");
    //one way to do sliders and labels
    //control.add("Label", "position=north; parent=controlPanel; text=phi; font=italic,8; foreground=blue; horizontalAlignment=center");
    //control.add("Slider","parent=controlPanel; minimum=0; maximum=90;variable=x; dragaction=sliderMoved()");
    //control.add("Label", "position=north; parent=controlPanel; text=theta; font=italic,8; foreground=blue; horizontalAlignment=center");
    //control.add("Slider", "parent=controlPanel; minimum=0; maximum=90;variable=y; dragaction=sliderMoved()");
    //control.add("Label", "position=north; parent=controlPanel; text=psi; font=italic,8; foreground=blue; horizontalAlignment=center");
    //control.add("Slider", "parent=controlPanel; minimum=0; maximum=90;variable=z; dragaction=sliderMoved()");
    //nicer way to do sliders and labels
    ControlElement myComponent1 = control.add("Slider","parent=controlPanel; minimum=0; maximum=90;variable=x; dragaction=sliderMoved()");
    ControlElement myComponent2 = control.add("Slider", "parent=controlPanel; minimum=0; maximum=90;variable=y; dragaction=sliderMoved()");
    ControlElement myComponent3 = control.add("Slider", "parent=controlPanel; minimum=0; maximum=90;variable=z; dragaction=sliderMoved()");
    myComponent1.setProperty("format","phi=0.00");
    myComponent2.setProperty("format","theta=0.00");
    myComponent3.setProperty("format","psi=0.00");
    //a way to list the properties associated with the slider follows:
    //java.util.ArrayList list = myComponent3.getPropertyList();
    //for (java.util.Iterator it = list.iterator(); it.hasNext(); ) {
    //  System.out.println ("Property = "+it.next().toString());
    //}
    //control.add("Button", "parent=controlPanel; text=Clear; action=clearPlot()");//not used here
*/

    //2nd large block - uncommented
    //this puts the sliders into the frame3d instead of on their own
    javax.swing.JPanel panel = new javax.swing.JPanel();
    //frame3d.getContentPane().add(panel,java.awt.BorderLayout.NORTH);
    frame3d.getContentPane().add(panel,java.awt.BorderLayout.SOUTH);
    Iflag=0; //flag used in sliderMoved()
    //GroupControl control = new GroupControl(this);
    control = new GroupControl(this);
    control.addObject(panel, "Panel", "name=controlPanel; layout=hbox");
    //This last line brings a Swing panel that you have previously created and
    //placed inside a drawingFrame3D into the control group.
    ControlElement myComponent1 = control.add("Slider","parent=controlPanel; minimum=0; maximum=90;variable=x; dragaction=sliderMoved()");
    ControlElement myComponent2 = control.add("Slider", "parent=controlPanel; minimum=0; maximum=90;variable=y; dragaction=sliderMoved()");
    ControlElement myComponent3 = control.add("Slider", "parent=controlPanel; minimum=0; maximum=90;variable=z; dragaction=sliderMoved()");
    myComponent1.setProperty("format","phi=0.00");
    myComponent2.setProperty("format","theta=0.00");
    myComponent3.setProperty("format","psi=0.00");
  }

  public void sliderMoved() {
    //double x = control.getDouble("x");
    ph=control.getDouble("x");
    th=control.getDouble("y");
    ps=control.getDouble("z");
    myControl.setValue("phi",ph);
    myControl.setValue("theta",th);
    myControl.setValue("psi",ps);
    if(Iflag!=0){calculate();}
    //dataset.append(x, x*x);//hint how to plot on drawpanel if used
    //drawpanel.repaint();
  }

   public void calculate() {
     myControl.clearMessages();
     if(Iflag==0){Iflag=1;} //flag used in sliderMoved()
     ph=myControl.getDouble("phi");
     th=myControl.getDouble("theta");
     ps=myControl.getDouble("psi");
     rp     = new double [N][N];
     rpp    = new double [N][N];
     rppp   = new double [N][N];
     phi=ph*rc;
     the=th*rc;
     psi=ps*rc;
     double[][] r={{x,0.,0.},{0.,y,0.},{0.,0.,z}}; //starting cartesian axis
     //starting Cartesian axes
     arrow[0].setSizeXYZ(r[0][0],r[0][1],r[0][2]);
     arrow[1].setSizeXYZ(r[1][0],r[1][1],r[1][2]);
     arrow[2].setSizeXYZ(r[2][0],r[2][1],r[2][2]);
     //phi rotation about z by ph degrees use -ph to get ccw rotation here
     //we also create a cc rotation matrix rphie for later use
     double[][] rphi={{Math.cos(-phi),-Math.sin(-phi),0.},
             {Math.sin(-phi),Math.cos(-phi),0.},{0.,0.,1.0}}; //%ccw z matrix rotation
     double[][] rphie={{Math.cos(phi),-Math.sin(phi),0.},
             {Math.sin(phi),Math.cos(phi),0.},{0.,0.,1.0}};   //%cc z matrix rotation
     for (int i=0; i < N; i++){
       for(int j=0; j < N; j++){
         rp[i][j]=0.0;
         for (int k=0; k < N; k++){
           rp[i][j]=rp[i][j]+rphi[i][k]*r[k][j];
         }
       }
     }
     arrow[3].setSizeXYZ(rp[0][0],rp[0][1],rp[0][2]);
     arrow[4].setSizeXYZ(rp[1][0],rp[1][1],rp[1][2]);
     arrow[5].setSizeXYZ(rp[2][0],rp[2][1],rp[2][2]);
     //theta rotation about x' by th degrees from y axis use -th to get ccw rotation here
     //we also create a cc rotation matrix rthee for later use
     double [][] rthe={{1.0,0.,0.},{0.,Math.cos(-the),-Math.sin(-the)},
             {0.,Math.sin(-the),Math.cos(-the)}};   //ccw x' matrix rotation
     double [][] rthee={{1.0,0.,0.},{0.,Math.cos(the),-Math.sin(the)},
             {0.,Math.sin(the),Math.cos(the)}};     //cc x' matrix rotation
     for (int i=0; i < N; i++){
       for(int j=0; j < N; j++){
         rpp[i][j]=0.0;
         for (int k=0; k < N; k++){
           rpp[i][j]=rpp[i][j]+rthe[i][k]*rp[k][j];
         }
       }
     }
     arrow[6].setSizeXYZ(rpp[0][0],rpp[0][1],rpp[0][2]);
     arrow[7].setSizeXYZ(rpp[1][0],rpp[1][1],rpp[1][2]);
     arrow[8].setSizeXYZ(rpp[2][0],rpp[2][1],rpp[2][2]);
     //psi rotation about z'' by ps degrees use -ps to get ccw rotation here
     double[][] rpsi={{Math.cos(-psi),-Math.sin(-psi),0.},
                {Math.sin(-psi),Math.cos(-psi),0.},{0.,0.,1.0}}; //%z'' matrix rotation
     for (int i=0; i < N; i++){
       for(int j=0; j < N; j++){
       rppp[i][j]=0.0;
          for (int k=0; k < N; k++){
            rppp[i][j]=rppp[i][j]+rpsi[i][k]*rpp[k][j];
          }
       }
     }
     arrow[9].setSizeXYZ(rppp[0][0],rppp[0][1],rppp[0][2]);
     arrow[10].setSizeXYZ(rppp[1][0],rppp[1][1],rppp[1][2]);
     arrow[11].setSizeXYZ(rppp[2][0],rppp[2][1],rppp[2][2]);
     //Draw and rotate planes
     //first original plane, later work with its copy
     ellipso[0].setSizeXYZ(2.*x,2.*y,0.); //flat unit ellipsoids are the desired planes
     ellipso[1].setSizeXYZ(2.*x,2.*y,0.); //copy of flat unit ellipsoid for later use
     int red=0, green=0, blue=255, transparency=25; //comes out light blue
     Color col=new Color(red,green,blue,transparency);
     ellipso[0].getStyle().setFillColor(col);
     ellipso[0].getStyle().setLineColor(col);
     ellipso[0].getStyle().setResolution(new Resolution(1, 32, 8));
     //second plane
     //recall rphie cc rotation about z, rthee cc rotation about x' here in reverse order
     double[][] A=new double[N][N];
     for (int i=0; i < N; i++){
       for(int j=0; j < N; j++){
       A[i][j]=0.0;
          for (int k=0; k < N; k++){
            A[i][j]=A[i][j]+rphie[i][k]*rthee[k][j]; //create combined rotation
          }
       }
     }
     m3dt=new Matrix3DTransformation(A);       //create the transformation matrix
     ellipso[1].setTransformation(m3dt);       //transform the copy of the original ellipse
     red=255; green=0; blue=0; transparency=25;//comes out light pink
     col=new Color(red,green,blue,transparency);
     ellipso[1].getStyle().setFillColor(col);
     ellipso[1].getStyle().setLineColor(col);
     ellipso[1].getStyle().setResolution(new Resolution(1, 32, 8));
     frame3d.addElement(arrow[0]);
     frame3d.addElement(arrow[1]);
     frame3d.addElement(arrow[2]);
     frame3d.addElement(arrow[3]);
     frame3d.addElement(arrow[4]);
     frame3d.addElement(arrow[5]);
     frame3d.addElement(arrow[6]);
     frame3d.addElement(arrow[7]);
     frame3d.addElement(arrow[8]);
     frame3d.addElement(arrow[9]);
     frame3d.addElement(arrow[10]);
     frame3d.addElement(arrow[11]);
     frame3d.addElement(ellipso[0]);
     frame3d.addElement(ellipso[1]);
     frame3d.render();
     myControl.println("Can also move the bottom left sliders to");
     myControl.println("change the angles. Left is phi, center is theta,");
     myControl.println("and right is psi. For specific angles, enter them");
     myControl.println("above and press 'calculate'. Press");
     myControl.println("'reset' to start over.");
     myControl.println("Feel free to mouse-rotate the 3D plots.");
     myControl.println("x,y,z axes -blue");
     myControl.println("x',y',z' axes - black: phi rotation");
     myControl.println("x'',y'',z'' axes - red: theta rotation ");
     myControl.println("x''',y''',z''' axes - green: psi rotation");
     myControl.println("the x'=x'' axis ends up red and is the line of nodes.");
   }

   public void clear  () {
     myControl.clearMessages();
     for(int i=0; i< Iarrows; i++){
       arrow[i].setXYZ(0, 0, 0); //x,y,z arrow position
       //arrow length corresponds to the vector function x,y,z components
       arrow[i].setSizeXYZ(0, 0, 0);
       //clear the arrows from the frame
       frame3d.addElement(arrow[i]);
     }
     frame3d.render();
   }

   public void resetCalculation() {
     clear();
     myControl.println ("Shows Euler angles: phi, theta, psi; the planes and");
     myControl.println ("the line of nodes. Here the rotation matrix is used:");
     myControl.println ("cos(angle) -sin(angle)");
     myControl.println ("sin(angle)  cos(angle)");
     myControl.println ("where the angle can be any of the inputed ones.");
     myControl.println ("The rotation about x, y, or z axes is handled");
     myControl.println ("by a 3X3 matrix with the element 1 on the");
     myControl.println ("corresponding axis of rotation.");
     myControl.println("press 'calculate' to begin with starting values.");
     //set frame as desired
     frame3d.setPreferredMinMax(xmin,xmax,ymin,ymax,zmin,zmax);
     frame3d.setAltitude(0.2);
     frame3d.setAzimuth(0.5);
     frame3d.getCamera().setXYZ(2.*x,2.*y,2.*z);
     r      = null;
     rp     = null;
     rpp    = null;
     rppp   = null;
     ph=15;
     th=25;
     ps=35;
     myControl.setValue("phi",ph);
     myControl.setValue("theta",th);
     myControl.setValue("psi",ps);
     //set the initial slider values as well [see sliderMoved()]
     control.setValue("x",ph);
     control.setValue("y",th);
     control.setValue("z",ps);
  }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }


   public static void main(String[] args) {
     Calculation model = new euler_angApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(490, 5);
     myControl.setSize(300,500);
     myControl.setDividerLocation(135);
     model.setControl(myControl);
   }
}