/*
@Author J E Hasbun 2007.
Draws a cube with the principal axes based on the entered symmetric inertia tensor.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import java.awt.Color;

public class cube_princ_axApp implements Calculation {
  Display3DFrame frame3d = new Display3DFrame("Unit Cube and Principal Axes");
  EigenJacobi eigenjacobi=new EigenJacobi();
  MatPrint mprint=new MatPrint();
  ElementTrail [] trace;
  ElementPolygon poly = new ElementPolygon();
  private Control myControl;
  int NumColRow=3, iflag=0, Nt=6;
  int icmax;
  double tol;
  double [][] A,DiagD,V, data;
  double xmin,xmax,ymin,ymax,zmin,zmax;

  public cube_princ_axApp(){
    frame3d.setLocation(5,5);
    frame3d.setSize(400,500);
  }

   public void calculate() {
     clear();
     cleartraces(Nt);
     myControl.clearMessages();
     A=(double [][])myControl.getObject("A");
     if(A[0].length!=A.length){
       System.out.println("array is not square");
     }
     tol=myControl.getDouble("tolerance");
     icmax=myControl.getInt("max_iter");
     V=new double [NumColRow][NumColRow];
     DiagD=new double[NumColRow][NumColRow];
     eigenjacobi.EigenJacobi(A,DiagD,V,NumColRow,tol,icmax);//input: A, output: DiagD,V
     //sort eigenvalues from low to high, and corresponding eigenvectors
     double tmp;
     for (int i=0; i< NumColRow; i++){
       for (int j=i; j<NumColRow; j++){
         if(DiagD[j][j] < DiagD[i][i]){       //swap eigenvalues
           tmp=DiagD[i][i];
           DiagD[i][i]=DiagD[j][j];
           DiagD[j][j]=tmp;
           for (int k=0; k < NumColRow; k++){ //swap corresponding eigenvectors
             tmp=V[k][i];
             V[k][i]=V[k][j];
             V[k][j]=tmp;
           }
         }
       }
     }
     //find largest element of the eigenvector matrix
     double VV=-9999.0;
     for(int i=0; i< NumColRow; i++){
       for(int j=0; j< NumColRow; j++){
         //if( (V[i][0]>VV) & (V[i][0]!=0) ){VV=V[i][0];}
         if(Math.abs(V[j][i])>VV){VV=Math.abs(V[j][i]);}
       }
     }
     if(VV<-9998)VV=1.0;  //that is if vv did not change, then just use 1.0
     //normalize V according to eigenvector with lowest eigenvalue (high symmetry)
     for (int i=0; i < NumColRow; i++){
       for(int j=0; j < NumColRow; j++){
         V[i][j]=V[i][j]/VV;
       }
     }
     mprint.ContPrintMat2D(A,NumColRow,NumColRow,"Input inertia tensor");
     mprint.ContPrintMat2D(DiagD,NumColRow,NumColRow,"Eigenvalues (principal moments)");
     mprint.ContPrintMat2D(V,NumColRow,NumColRow,"Eigenvectors (axes directions)");
     for(int i=0; i< Nt; i++){
       trace[i].addPoint(0., 0., 0.); //all traces begin from the origin
     }
     trace[0].addPoint(V[0][0],V[1][0],V[2][0]); //the main axes from eigen vectors
     trace[2].addPoint(V[0][1],V[1][1],V[2][1]);
     trace[4].addPoint(V[0][2],V[1][2],V[2][2]);
     int imax=50;  //used to make the mirror images (odd traces - dashed lines)
     for(int i=1; i< imax; i++){ //since the line is not connected - need more points
       double ss=i/(imax-1.0);
       trace[1].addPoint(-ss*V[0][0],-ss*V[1][0],-ss*V[2][0]);
       trace[3].addPoint(-ss*V[0][1],-ss*V[1][1],-ss*V[2][1]);
       trace[5].addPoint(-ss*V[0][2],-ss*V[1][2],-ss*V[2][2]);
     }
     //Six sided polygon with known vertices from the matlab cube_princ_ax.m script
     //The vertices are such as to touch half way the six sides of the cube.
     //Here we need the double[nPoints][3] array; i.e., three columns (x,y,z) per point
     double [][] data={   //the matlab numbers were shifted back by 0.5 here
         {0.5,-0.5,0.0},
         {0.0,-0.5,0.5},
         {-0.5,0.0,0.5},
         {-0.5,0.5,0.0},
         {0.0,0.5,-0.5},
         {0.5,0.0,-0.5}};
     poly.setClosed(true);
     poly.setSizeXYZ(2*VV,2*VV,2*VV);
     //poly.getStyle().setFillColor(java.awt.Color.green);
     //0,0,0,0=black transparent; 255,255,255,255=solid white
     int red=0, green=255, blue=0, transparency=50; //comes out light green
     //int red=255, green=0, blue=0, transparency=127; //comes out light pink
     Color col=new Color(red,green,blue,transparency);
     poly.getStyle().setFillColor(col);
     //poly.getStyle().setDrawingFill(false); //no color or hollow
     poly.setData(data);
     xmin=-VV; xmax=VV; ymin=-VV; ymax=VV; zmin=-VV; zmax=VV;
     frame3d.setPreferredMinMax(xmin, xmax, ymin, ymax, zmin, zmax);
     for(int i=0; i<Nt; i++){
       frame3d.addElement(trace[i]);
     }
     frame3d.getCamera().setXYZ(2.5*VV,2.5*VV,2.5*VV);
     frame3d.addElement(poly);
     frame3d.setAltitude(0.5);
     frame3d.setAzimuth(0.5);
     frame3d.render();
     myControl.println("Principal axes: black-high symmetry, & the two perpendicular");
     myControl.println("blue and red lie on the same plane as the hexagon.");
     myControl.println("The eigenvalues are the principal moments of inertia and");
     myControl.println("the eigenvectors are their corresponding direction cosines");
     }

   public void clear  () {
     poly.setSizeXYZ(0.0,0.0,0.0);
     frame3d.render();
   }

   public void resetCalculation() {
     myControl.clearMessages();
     if(iflag !=0){cleartraces(Nt);}
     clear();
     initialtraces(Nt);
     iflag=1;
     myControl.println ("Draws a cube with the principal axes based on the entered");
     myControl.println ("symmetric inertia tensor. The parameters tolerance and");
     myControl.println ("maxiter are used in the Jacobi iteration eigenvector finder.");
     icmax=100;
     tol=1.e-5;
     double a11=2./3., a12=-1./4.,a13=-1./4.;
     double a21=-1./4., a22=2./3., a23=-1./4.;
     double a31=-1./4., a32=-1./4., a33=2./3.;
     myControl.setValue("A",new double[][]{{a11,a12,a13},{a21,a22,a23},
                        {a31,a32,a33}});
     myControl.setValue("tolerance",tol);
     myControl.setValue("max_iter",icmax);
   }

   public void cleartraces(int N) {
    //clear the traces from the frame
       for (int i = 0; i < N; i++) {
         trace[i].clear();
         frame3d.addElement(trace[i]);
       }
   }

   public void initialtraces(int N) {
     //initialize the traces
     trace = new ElementTrail[N];
     for (int i = 0; i < N; i++) {
       trace[i] = new ElementTrail();
       trace[i].getStyle().setLineWidth((float)1.5);//all lines same thickness
     }
     for (int i = 1; i < N; i += 2) {
       trace[i].setConnected(false); //the odd ones are not connected
     }
     if ((N == Nt) && (Nt==6)){
        trace[0].getStyle().setLineColor(java.awt.Color.black);
        trace[1].getStyle().setLineColor(java.awt.Color.black);
        trace[2].getStyle().setLineColor(java.awt.Color.blue);
        trace[3].getStyle().setLineColor(java.awt.Color.blue);
        trace[4].getStyle().setLineColor(java.awt.Color.red);
        trace[5].getStyle().setLineColor(java.awt.Color.red);
     }else{
       myControl.println("N or Nt not equal to 6 in initialtraces()");
     }
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
     mprint.myControl=control;
   }

   public static void main(String[] args) {
     Calculation model = new cube_princ_axApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(375,500);
     myControl.setDividerLocation(150);
     model.setControl(myControl);
   }
}