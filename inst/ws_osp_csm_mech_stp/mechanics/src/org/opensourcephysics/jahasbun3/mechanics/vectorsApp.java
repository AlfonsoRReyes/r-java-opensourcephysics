/*
@Author J E Hasbun 2007.
//Application to add two 3-D vectors to obtain a resultant
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.frames.Display3DFrame;
import java.text.NumberFormat;

public class vectorsApp implements Calculation {
  Display3DFrame frame3D = new Display3DFrame("Vectors Application");
  //DrawingPanel3D panel3D= new DrawingPanel3D();// could be used also
  ElementArrow[] arrow = new ElementArrow[3]; // create an Arrow
  ElementCircle [] point=new ElementCircle[2];
  String [] vec_name =new String [3];
  NumberFormat nf = NumberFormat.getInstance();
  double [][] a;
  double [] min, max;
  double angle1,angle2,angle3, pi, a1,a2,a3,a3t,a3f,del=1.e-12;
  private Control myControl;

  public vectorsApp (){
    frame3D.setDecorationType(VisualizationHints.DECORATION_AXES);
    frame3D.setLocation(10,10);
    frame3D.setSize(400,500);
    frame3D.setAutoclear(true);
    arrow[0]=new ElementArrow(); //define 1st vector
    arrow[1]=new ElementArrow(); //define 2nd vector
    arrow[2]=new ElementArrow(); //define resultant vector
    point[0] = new ElementCircle(); //define origin point
    point[1] = new ElementCircle(); //define an interactive point
    vec_name[0]="A: 1st input vector (black): ";
    vec_name[1]="B: 2nd input vector ( blue): ";
    vec_name[2]="C: Resultant vector (  red): ";
    nf.setMaximumFractionDigits(3);
    pi=Math.PI;
  }

   public void calculate() {
     myControl.clearMessages();
     a= new double[4][3];             //Vector array
     a[0][0]=0; a[0][1]=0; a[0][2]=0; //origin position
     a[1]=(double [])myControl.getObject("A");
     a[2]=(double [])myControl.getObject("B");
     arrow[0].setXYZ(a[0][0],a[0][1],a[0][2]);     //1st vector starts a origin
     arrow[0].setSizeXYZ(a[1][0],a[1][1],a[1][2]); //1st vector
     arrow[0].getStyle().setFillColor(java.awt.Color.BLACK);
     arrow[0].getStyle().setLineColor(java.awt.Color.BLACK);
     arrow[1].setXYZ(a[1][0],a[1][1],a[1][2]);     //2nd vector starts at head of 1st arrow
     arrow[1].setSizeXYZ(a[2][0],a[2][1],a[2][2]); //2nd vector
     arrow[1].getStyle().setFillColor(java.awt.Color.BLUE);
     arrow[1].getStyle().setLineColor(java.awt.Color.BLUE);
     // add the 1st and second vectors
     for(int i = 0; i <  3; i++) {
       a[3][i] = a[1][i] + a[2][i];
      }
     myControl.setValue("C",new double[]{a[3][0],a[3][1],a[3][2]});
     arrow[2].setXYZ(a[0][0],a[0][1],a[0][2]);     //resultant starts at origin
     arrow[2].getStyle().setFillColor(java.awt.Color.RED);
     arrow[2].getStyle().setLineColor(java.awt.Color.RED);
     arrow[2].setSizeXYZ(a[3][0],a[3][1],a[3][2]); //the resultant
     point[0].setXYZ(a[0][0], a[0][1], a[0][2]);//place a point at the end of C
     point[0].setSizeXYZ(0.2, 0.2, 0.2);        // point size
     point[0].getStyle().setFillColor(java.awt.Color.darkGray);
     point[1].setXYZ(a[3][0], a[3][1], a[3][2]);//place a point at the end of C
     point[1].setSizeXYZ(0.25, 0.25, 0.25);     //point size
     point[1].getStyle().setFillColor(java.awt.Color.GREEN);
     //minimum x, y, z component of the three vectors - for plot use
     double [] min={0,0,0};
     double [] max={0,0,0};
     for (int i=1; i<4; i++){
       for (int j=0; j<3; j++){
         if(a[i][j] < min[j]){ min[j]=a[i][j];}
         if(a[i][j] > max[j]){ max[j]=a[i][j];}
       }
     }
     a1=VectorMath.magnitude(a[1]);
     a2=VectorMath.magnitude(a[2]);
     a3=VectorMath.magnitude(a[3]);
     angle1=Math.acos(VectorMath.dot(a[1],a[2])/(a1*a2+del));
     a3t=Math.acos(a[3][2]/(a3+del));
     if(a[3][0]==0){a[3][0]=del;}
     a3f=Math.atan(a[3][1]/a[3][0]);
     frame3D.setPreferredMinMax(min[0], max[0], min[1], max[1], min[2], max[2]);
     frame3D.setAltitude(0.5);
     frame3D.setAzimuth(0.5);
     frame3D.addElement(arrow[0]);//Vector A
     frame3D.addElement(arrow[1]);//Vector B
     frame3D.addElement(arrow[2]);//The resultant
     frame3D.addElement(point[0]);//The origin point
     point[1].getInteractionTarget(Element.TARGET_POSITION).setEnabled(true);
     frame3D.addElement(point[1]);//Interactive point at the end of C
     //could instead do the commented lines below with the DrawingPanel3D
/*
     panel3D.setPreferredMinMax(min[0], max[0], min[1], max[1], min[2], max[2]);
     panel3D.addElement(arrow[0]);
     panel3D.addElement(arrow[1]);
     panel3D.addElement(arrow[2]);
     panel3D.addElement(point[0]);
     panel3D.addElement(point[1]);
     frame3D.setDrawingPanel3D(panel3D);
*/
     //printing
     myControl.println("** Mouse-Rotate the coordinate system to change the view **");
     myControl.println("** Feel Free to move the green interactive point at the end of C **");
     myControl.println("-- The dark grey point marks the origin {0,0,0} --");
     //print the input vectors
     for(int i=1; i<4;i++){
     myControl.print(vec_name[i-1]+" {");
       for(int j = 0; j < 3; j++ ) {
         if(j!=2){myControl.print(a[i][j]+",");}
         else{myControl.print(a[i][j]+"");}
       }
     myControl.println("}");
     }
     myControl.println("Magnitudes: A="+nf.format(a1)+" B="+
                       nf.format(a2)+" C="+nf.format(a3));
     myControl.println("Angle between A and B="+nf.format(angle1)+" rad"+" or "+
                       nf.format(angle1*180/Math.PI)+" degrees");
     myControl.println("Polar Angle (C with z axis)="+nf.format(a3t)+" rad or "+
                       nf.format(a3t*180/pi)+" degrees");
     myControl.println("Azimuth Angle (C with x axis)="+nf.format(a3f)+" rad or "+
                       nf.format(a3f*180/pi)+" degrees");
     angle2=Math.acos(VectorMath.dot(a[3],a[1])/(a3*a1+del));
     angle3=Math.acos(VectorMath.dot(a[3],a[2])/(a3*a2+del));
     myControl.println("Angle between C and A="+nf.format(angle2)+" rad"+" or "+
                       nf.format(angle2*180/Math.PI)+" degrees");
     myControl.println("Angle between C and B="+nf.format(angle3)+" rad"+" or "+
                       nf.format(angle3*180/Math.PI)+" degrees");
     myControl.println("------------------  Done ----------------------------");
   }

   public void clear  () {
     myControl.clearMessages();
     a=null;
     arrow[0].setXYZ(0,0,0);
     arrow[0].setSizeXYZ(0,0,0);
     arrow[1].setXYZ(0,0,0);
     arrow[1].setSizeXYZ(0,0,0);
     arrow[2].setXYZ(0, 0, 0);
     arrow[2].setSizeXYZ(0, 0, 0);
     point[0].setXYZ(0, 0, 0);
     point[0].setSizeXYZ(0, 0, 0);
     point[1].setXYZ(0, 0, 0);
     point[1].setSizeXYZ(0, 0, 0);
     frame3D.addElement(arrow[0]);
     frame3D.addElement(arrow[1]);
     frame3D.addElement(arrow[2]);
     frame3D.addElement(point[0]);
     frame3D.addElement(point[1]);
     frame3D.repaint();
   }

   public void resetCalculation() {
     clear();
     myControl.println("Application to add two 3-D vectors to obtain resultant C.");
     myControl.println("Double-click the A, B arrays to change their component values.");
     myControl.setValue("A",new double[]{-3,4,8,});
     myControl.setValue("B",new double[]{2,-2,6});
     myControl.setValue("C",new double[]{0,0,0});
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

  public static void main(String[] args) {
    Calculation model = new vectorsApp();
    CalculationControl myControl;
    myControl = new CalculationControl(model);
    myControl.setLocation(410, 10);
    myControl.setSize(380, 500);
    myControl.setDividerLocation(150);
    model.setControl(myControl);
  }
}
