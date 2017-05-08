/*
@Author J E Hasbun 2007.
Two dimensional collisions. It finds the approximate final velocities of two
particles and their scattering angles by a root finding mathod. It uses the
newton-Raphson method for a system of equations.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;
import java.text.*;

public class ecoll_2dApp implements Calculation {
   private Control myControl;
   PlottingPanel panel = new PlottingPanel("Px","Py","Elastic Collisions");
   DrawingFrame frame = new DrawingFrame(panel);
   NumberFormat nf = NumberFormat.getInstance();
   DecimalFormat df= new DecimalFormat("+0.000E00");
   MatPrint mprint=new MatPrint();
   int Iarrows=5;
   Trail [] trail=new Trail[2];
   Arrow [] arrows=new Arrow[Iarrows];
   double th10,th20,th1,m1,m2,v1i,v2i,tp10,tp20,tp1;
   double pi,cf,v1f,v2f,tp2;
   double p1ixs,p1ixe,p1iys,p1iye,p2ixs,p2ixe,p2iys,p2iye;
   double p1fxs,p1fxe,p1fys,p1fye,p2fxs,p2fxe,p2fys,p2fye;
   double Pxi,Pxf,Pyi,Pyf,Ei,Ef;
   double minx,maxx,miny,maxy;
   int maxiter=25, Ndim=3;
   double tol=5e-3;

   public ecoll_2dApp(){
      frame.setLocation(10,5);
      frame.setSize(390,500);
      nf.setMaximumFractionDigits(4);
      pi=Math.PI; cf=pi/180;
      for(int i=0; i < Iarrows; i++){
        arrows[i] = new Arrow(0.,0.,0.,0.);
      }
      trail[0]=new Trail();
      trail[0].setConnected(true);
      trail[1]=new Trail();
      trail[1].setConnected(true);
   }

   public void calculate() {
     double [] xx=new double[Ndim];
     double [] xxn=new double[Ndim];
     double [] F=new double[Ndim];
     //double [][] J    = new double [Ndim][Ndim];
     double [][] Jinv = new double [Ndim][Ndim];
     //double sum=10.0, a1=0.25, a2=0.25, a3=0.25;
     int Iterations=0;
     double err, relerr;
     clear();
     th10=myControl.getDouble("theta_10 (degrees)");
     th20=myControl.getDouble("theta_20 (degrees)");
     th1=myControl.getDouble ("theta_1 (degrees)");
     v1i=myControl.getDouble ("v1i (m/s)");
     v2i=myControl.getDouble ("v2i (m/s)");
     m1=myControl.getDouble  ("m1");
     m2=myControl.getDouble  ("m2");
     tp10=th10*cf;    //convert angles to radians
     tp20=th20*cf;
     tp1=th1*cf;
     //-- initial momenta plot setup
     p1ixs=-v1i*Math.cos(tp10);            //initial p m1 - line start point
     p1iys=v1i*Math.sin(tp10);
     p1ixe=0.0;                            //initial p m1 - line end point
     p1iye=0.0;
     p2ixs=-v2i*Math.cos(tp20);            //initial p m2 - line start point
     p2iys=-v2i*Math.sin(tp20);
     p2ixe=0.0;                            //initial p m2 - line end point
     p2iye=0.0;
     //myControl.println("Initial");
     //myControl.println("p1ixs="+p1ixs+", p1iys="+p1iys+", p1ixe="+p1ixe+", p1iye="+p1iye);
     //myControl.println("p2ixs="+p2ixs+", p2iys="+p2iys+", p2ixe="+p2ixe+", p2iye="+p2iye);
     //
     //use a grid method to search for a guess for the uknown final momenta and angle
     //xx[0]=v1i/2; xx[1]=v1i*Math.sqrt(0.75*m1/m2); xx[2]=pi/2-tp1; //old guesses
     //upper limits, steps
     double v1L, v2L, thL, v1s, v2s, ths;
     v1L=10.+v1i/2.;                      //upper search range for v1f
     v2L=10.+v1i*Math.sqrt(0.75*m1/m2);   //upper search range for v1f
     thL=pi/2.;                           //upper search range for thL
     v1s=0.5;                            //search steps
     v2s=0.5;
     ths=0.5;
     gridder(xx,v1L,v2L,thL,v1s,v2s,ths); //search approx roots - grid method
     //mprint.ContPrintMat1D(xx,Ndim,"xx -- guess");
     myControl.println("initial: p1(black), p2 (magenta); final: p1 (blue), p2 (red)");
     myControl.println("green arrow indicates the collision direction.");
     myControl.println("Guesses: v1f="+nf.format(xx[0])+", v2f="+nf.format(xx[1])+
                       ", tp2="+nf.format(xx[2]/cf)+" degrees");
     err=9999.;
     relerr=9999.;
     //Use the Newton-Raphson method for systems of equations - employs the Jacobian
     //Needs a good guess - use one found by the grid method above
     while (err > tol*1.e-6 && relerr > tol*1.e-6 && Iterations < maxiter){
       Iterations++;
       //J=Jacobian(Ndim,xx,tol/100.);  //Jacobian
       //double det=J_inverse(J,Jinv);  //inverse and determinant if needed
       F=feqs(xx);                    //the functions
       for (int i = 0; i < Ndim; i++) {
         double sumjF = 0;
         for (int j = 0; j < Ndim; j++) {
           sumjF = sumjF + Jinv[i][j] * F[j];
         }
         xxn[i] = xx[i] - sumjF;
       }
       err=(xx[0]-xxn[0])*(xx[0]-xxn[0]);
       relerr=xx[0]*xx[0];
       xx[0]=xxn[0];
       for (int i=1; i<Ndim; i++){
         err=err+(xx[i]-xxn[i])*(xx[i]-xxn[i]);
         relerr=relerr+xx[i]*xx[i];
         xx[i] = xxn[i];
       }
       err=Math.sqrt(err);
       relerr=err/(relerr+tol);
       //mprint.ContPrintMat1D(xx,Ndim,"xx-in while");
       //myControl.println("err="+err+", relerr="+relerr);
     }
     myControl.println("Iterations="+Iterations+", err="+df.format(err)+
                       "\nrelerr="+df.format(relerr));
     v1f=xx[0];
     v2f=xx[1];
     tp2=xx[2];
     myControl.println("Final speeds and exit angle");
     myControl.println("v1f="+nf.format(v1f)+" m/s, v2f="+nf.format(v2f)+
                       " m/s\ntheta_2="+nf.format(tp2)+" rad or "+
                       nf.format(tp2/cf)+" degrees");
     //-- final momenta plot setup
     p1fxs=0.0;                            //final p m1 - line start point
     p1fys=0.0;
     p1fxe=v1f*Math.cos(tp1);              //final p m1 - line end point
     p1fye=v1f*Math.sin(tp1);
     p2fxs=0.0;                            //final p m2 - start line point
     p2fys=0.0;
     p2fxe=v2f*Math.cos(tp2);              //final p m2 - end line point
     p2fye=-v2f*Math.sin(tp2);
     //plots
     double [] tmpx=new double[]{p1ixs,p1ixe,p2ixs,p2ixe,p1fxs,p1fxe,p2fxs,p2fxe};
     double [] tmpy=new double[]{p1iys,p1iye,p2iys,p2iye,p1fys,p1fye,p2fys,p2fye};
     minx=ArrayLib.min(tmpx);
     miny=ArrayLib.min(tmpy);
     maxx=ArrayLib.max(tmpx);
     maxy=ArrayLib.max(tmpy);
     //x line
     trail[0].addPoint(minx,0);
     trail[0].addPoint(maxx,0);
     trail[0].color=java.awt.Color.black;
     trail[0].setDashedStroke(1,3);
     //y line
     trail[1].addPoint(0,miny);
     trail[1].addPoint(0,maxy);
     trail[1].color=java.awt.Color.black;
     trail[1].setDashedStroke(1,3);
     //initial p m1  - arrow begin, end point
     arrows[0].setXY(p1ixs,p1iys);
     arrows[0].setXlength(p1ixe-p1ixs);
     arrows[0].setYlength(p1iye-p1iys);
     arrows[0].setColor(java.awt.Color.black);
     arrows[0].setHeadSize(6.0f);
     //initial p m2  - arrow begin, end point
     arrows[1].setXY(p2ixs,p2iys);
     arrows[1].setXlength(p2ixe-p2ixs);
     arrows[1].setYlength(p2iye-p2iys);
     arrows[1].setColor(java.awt.Color.magenta);
     arrows[1].setHeadSize(6.0f);
     //final p m1  - arrow begin, end point
     arrows[2].setXY(p1fxs,p1fys);
     arrows[2].setXlength(p1fxe-p1fxs);
     arrows[2].setYlength(p1fye-p1fys);
     arrows[2].setColor(java.awt.Color.blue);
     arrows[2].setHeadSize(6.0f);
     //final p m2  - line begin, end point
     arrows[3].setXY(p2fxs,p2fys);
     arrows[3].setXlength(p2fxe-p2fxs);
     arrows[3].setYlength(p2fye-p2fys);
     arrows[3].setColor(java.awt.Color.red);
     arrows[3].setHeadSize(6.0f);
     //arrows[4]=Arrow(0,maxy*(1-0.1)+0.1,0.25*maxx,0);//arrow indicated collision direction
     arrows[4].setXY(0.0,maxy*(1-0.1)+0.1);
     arrows[4].setXlength(0.25*maxx);
     arrows[4].setYlength(0.0);
     arrows[4].setColor(java.awt.Color.green);
     arrows[4].setHeadSize(6.0f);
     frame.setPreferredMinMax(minx,maxx,miny,maxy);
     frame.setSquareAspect(true);
     for(int i=0; i < Iarrows; i++){
        frame.addDrawable(arrows[i]);
     }
     frame.addDrawable(trail[0]);
     frame.addDrawable(trail[1]);
     frame.render();
     //myControl.println("Final");
     //myControl.println("p1fxs="+p1fxs+", p1fys="+p1fys+", p1fxe="+p1fxe+", p1fye="+p1fye);
     //myControl.println("p2fxs="+p2fxs+", p2fys="+p2fys+", p2fxe="+p2fxe+", p2fye="+p2fye);
     //total energy and momentum conservation check
     Ei=(m1*v1i*v1i+m2*v2i*v2i)/2;
     Ef=(m1*v1f*v1f+m2*v2f*v2f)/2;
     Pxi=m1*v1i*Math.cos(tp10)+m2*v2i*Math.cos(tp20);
     Pxf=m1*v1f*Math.cos(tp1)+m2*v2f*Math.cos(tp2);
     Pyi=-m1*v1i*Math.sin(tp10)+m2*v2i*Math.sin(tp20);
     Pyf=m1*v1f*Math.sin(tp1)-m2*v2f*Math.sin(tp2);
     myControl.println("Momentum, Energy check: results are best if final = initial");
     myControl.println("Ei="+nf.format(Ei)+", Ef="+nf.format(Ef));
     myControl.println("Pxi="+nf.format(Pxi)+", Pxf="+nf.format(Pxf));
     myControl.println("Pyi="+nf.format(Pyi)+", Pyf="+nf.format(Pyf));
   }

   public double[][] Jacobian (int n, double xx[],double tol) {
     //builds the jacobian
     //xxp, xxm contain the varied parameter values to evaluate the equations on
     //fp, fm comtain the transpose of the function evaluations needed in the Jacobian
     //The Jacobian is calculated by the finite derivative method
     double [][]   J=new double [n][n];
     double [][] xxp=new double [n][n];
     double [][] xxm=new double [n][n];
     double [][]  fp=new double [n][n];
     double [][]  fm=new double [n][n];
     //build the coordinates for the derivatives
     for (int i=0; i<n; i++){
       for (int j=0; j<n; j++){
         xxp[i][j]=xx[j];
         xxm[i][j]=xx[j];
       }
       xxp[i][i]=xxp[i][i]+tol;
       xxm[i][i]=xxm[i][i]-tol;
     }
     for (int i=0; i<n; i++){ //f's here are built in transpose form
       fp[i]=feqs(xxp[i]);    //i=1: f's at x+tol; i=2 f's at y+tol, etc
       fm[i]=feqs(xxm[i]);    //i=1: f's at x-tol; i=2 f's at y-tol, etc
     }
     //Build the Jacobian by the differences methods
     //becasue the f's above are in transpose form, we swap i, j in the derivative
     for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
           J[i][j]=(fp[j][i]-fm[j][i])/tol/2.; //the derivatives df_i/dx_j
      }
     }
     //mprint.ContPrintMat2D(xxp,n,n,"xxp");
     return J;
   }

   public double [] feqs (double [] X){
     //defines the equations whose roots we seek
     double [] equation=new double [Ndim];
     //from the energy equation
     equation[0]=m1*v1i*v1i+m2*v2i*v2i-m1*X[0]*X[0]-m2*X[1]*X[1];
     //from the x-momentum equation
     equation[1]=+m1*v1i*Math.cos(tp10)+m2*v2i*Math.cos(tp20)
                -m1*X[0]*Math.cos(tp1)-m2*X[1]*Math.cos(X[2]);
     //from the y-momentum equation
     equation[2]=-m1*v1i*Math.sin(tp10)+m2*v2i*Math.sin(tp20)
                -m1*X[0]*Math.sin(tp1)+m2*X[1]*Math.sin(X[2]);
     return equation;
   }

   public double J_inverse(double A [][], double B [][]){
     //input:   A   - a 2D matrix
     //output:  B   - the inverse of A
     //returns: det - the determinant of A
     double det;
     double [][] inv;
     LUPDecomposition lu=new LUPDecomposition(A);
     det=lu.determinant();
     inv=lu.inverseMatrixComponents();
     int n=inv.length;
     mprint.copy2d(inv,B,n);
     //myControl.println(" det="+det);
     //mprint.ContPrintMat2D(inv,3,3,"inv");
     return det;
   }

   public void gridder (double [] c, double L1, double L2, double L3,
                        double v1s, double v2s, double ths){
     // use a grid method to obtain initial guesses to a function
     // of three variable whose zero we seek
     // Input: upper limits L1,L2,L3
     // Input: step sizes v1s, v2s, ths;
     // Output: c[0],c[1],c[3]
     // the array c contains the approximate minima positions when finished
     double sr0=9999, sr;
     int N1=Math.round((float)(L1/v1s));
     int N2=Math.round((float)(L2/v2s));
     int N3=Math.round((float)(L3/ths));
     double [] x=new double [3];
     for (int i1=0; i1<N1; i1++){
       x[0] = i1 * v1s;
       for (int i2 = 0; i2 < N2; i2++) {
         x[1] = i2 * v2s;
         for (int i3 = 0; i3 < N3; i3++) {
           x[2] = i3 * ths;
           sr=totalAbsFs(x);   //evaluate the function whose lowest value we seek
           if (sr < sr0){
             c[0]=x[0];
             c[1]=x[1];
             c[2]=x[2];
             sr0=sr;
           }
         }
       }
     }
   }

   public double totalAbsFs (double[] X){
     //return the sum of the squares of the functions whose roots we seek;
     double [] M=new double [Ndim];
     M=feqs(X);
     for(int i=0; i<Ndim; i++){
       M[i]=M[i]*M[i];
     }
     return ArrayLib.sum(M);
   }

   public void clear  () {
     myControl.clearMessages();
     trail[0].clear();
     trail[1].clear();
     for(int i=0; i < Iarrows; i++){
       arrows[i].setXY(0.,0.);
       arrows[i].setXlength(0.0);
       arrows[i].setYlength(0.0);
       arrows[i].setColor(java.awt.Color.white);
       arrows[i].setHeadSize(0.0f);
     }
     frame.clearData();
     frame.clearDrawables();
   }

   public void resetCalculation() {
     clear();
     myControl.println("Two dimensional collisions. It finds approximately the");
     myControl.println("final velocities of two particles and their scattering angles");
     myControl.println("by the Newton-Raphson method for systems of equations.");
     myControl.println("The energy and Momentum equations are:");
     myControl.println("Ei=(m1*v1i^2+m2*v2i^2)/2, Ef=(m1*v1f^2+m2*v2f^2)/2");
     myControl.println("Pxi=m1*v1i*cos(theta_10)+m2*v2i*cos(theta_20),");
     myControl.println("Pxf=m1*v1f*cos(theta_1)+m2*v2f*cos(theta_2),");
     myControl.println("Pyi=-m1*v1i*sin(theta_10)+m2*v2i*sin(theta_20),");
     myControl.println("Pyf=m1*v1f*sin(theta_1)-m2*v2f*sin(theta_2).");
     //masses, init angles,speeds
     th10 = 40;
     th20 = 20;
     th1  = 30;
     v1i  = 1.5;
     v2i  = 0.75;
     m1   = 0.5;
     m2   = 1.5;
     myControl.setValue("theta_10 (degrees)",th10);
     myControl.setValue("theta_20 (degrees)",th20);
     myControl.setValue("theta_1 (degrees)",th1);
     myControl.setValue("v1i (m/s)",v1i);
     myControl.setValue("v2i (m/s)",v2i);
     myControl.setValue("m1",m1);
     myControl.setValue("m2",m2);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
     mprint.myControl=control;
   }

   public static void main(String[] args) {
     Calculation model = new ecoll_2dApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,500);
     myControl.setDividerLocation(190);
     model.setControl(myControl);
   }
}
