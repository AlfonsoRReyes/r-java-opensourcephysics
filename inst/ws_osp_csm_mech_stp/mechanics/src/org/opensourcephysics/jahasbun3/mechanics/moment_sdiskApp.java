/*
@Author J E Hasbun 2007.
Finds the integral of 4*f(x)/pi, where f(x) is associated with the moment of
a disk.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import java.text.*;

public class moment_sdiskApp implements Calculation {
   private Control myControl;
   double [] xx,ff;
   double dx,x,a,b,pi;
   int NPTS;
   NumberFormat nf = NumberFormat.getInstance();
   DecimalFormat df = new DecimalFormat("+0.00000E00");

   public moment_sdiskApp(){

     pi=Math.PI;
     nf.setMaximumFractionDigits(4);

   }

   public void calculate() {
     clear();
     Function f;
     String fx = myControl.getString("f(x)");
     try { // read in function
     f = new ParsedFunction(fx);
     } catch(ParserException ex) {
       myControl.println(ex.getMessage());
       return;
     }
     xx=new double [NPTS];
     ff=new double [NPTS];
     a=myControl.getDouble("a");
     b=myControl.getDouble("b");
     NPTS=myControl.getInt("NPTS");
     if(NPTS%2.==0){NPTS=NPTS+1;} //make sure NPTS is an odd integer
     dx=(b-a)/(NPTS-1);
     myControl.setValue("dx[(b-a)/(NPTS-1)]=",nf.format(dx));
     xx = new double[NPTS+1];
     ff = new double[NPTS+1];
     for(int i = 0; i < NPTS; i++) {
       xx[i]=a+i*dx;               //time
       ff[i]=f.evaluate(xx[i]);   //position
       //System.out.println("t="+tt[i]+" x="+xx[i]);
     }
     double sum=Simp.Simp(ff,dx);
     sum=sum*4.0/pi;
     myControl.println("the Integral of 4*"+fx+"/pi="+nf.format(sum));

   }

   public void clear  () {
     myControl.clearMessages();
   }

   public void resetCalculation() {
     clear();
     myControl.println ("Finds the integral of 4*f(x)/pi, where f(x) is");
     myControl.println ("associated with the moment of a disk:");
     myControl.println ("Iyy=Ms*R^2*J, where J=4*Integral(f(x))/pi on [0,1].");
     myControl.println ("Simpson's rule is used here, Ms is the disk's mass and");
     myControl.println ("R is the disk's radius.");
     xx=null;
     ff=null;
     a=0.0;
     b=1.0;
     NPTS=301;
     myControl.setValue("f(x)", "x^2*(1-x^2)^(0.5)");
     myControl.setValue("a",a);
     myControl.setValue("b",b);
     myControl.setValue("NPTS", NPTS);
     dx=(b-a)/(NPTS-1);
     myControl.setValue("dx[(b-a)/(NPTS-1)]=",nf.format(dx));
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new moment_sdiskApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     //myControl.setLocation(520, 30);
     //myControl.setSize(250,505);
     myControl.setLocation(270, 30);
     myControl.setSize(430,505);
     model.setControl(myControl);
   }
}
