/*
@Author J E Hasbun 2007.
Plots the phase difference between the driving force and
the solution of the driven HO.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
//import org.opensourcephysics.display.*;
import org.opensourcephysics.controls.*;
//import org.opensourcephysics.numerics.*;
//import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;


public class drive_phaseApp implements Calculation {
   PlotFrame plot= new PlotFrame("w_D","Phase(w)","Phase(w_D)");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double phase[][], ww[], cc[], gamma[];
   double den,m,k,cmin,cmax,w0,wmin,wmax,dw,dc;
   int Nc,Nw;

   public drive_phaseApp(){
      //plot.setConnected(0,true);
      //plot.setMarkerSize(0,1);
      //plot.setVisible(true);
      plot.setLocation(10,5);
      plot.setSize(390,575);
      nf.setMaximumFractionDigits(3);
   }

   public void calculate() {
     myControl.clearMessages();
     m=myControl.getDouble("m");
     k=myControl.getDouble("k");
     w0=Math.sqrt(k/m);
     wmin=myControl.getDouble("wmin");
     wmax=myControl.getDouble("wmax");
     Nw=myControl.getInt("Nw");
     dw=(wmax-wmin)/Nw;
     cmin=myControl.getDouble("cmin");
     cmax=myControl.getDouble("cmax");
     Nc=myControl.getInt("Nc");
     dc=(cmax-cmin)/Nc;
     myControl.setValue("w0=",w0);
     myControl.setValue("dc((cmax-cmin)/Nc)=",dc);
     myControl.setValue("dw((wmax-wmin)/Nw)=",dw);
     ww     = new double[Nw+1];
     cc     = new double[Nc+1];
     gamma  = new double[Nc+1];
     phase  = new double[Nc+1][Nw+1];
     for (int i=0; i<=Nc;i++){
       cc[i]=cmin+i*dc;
       gamma[i] = cc[i] / 2 / m;
         for (int j = 0; j <= Nw; j++) {
           if (i==0){ww[j]=wmin+j*dw;} // build ww array once
           den=w0*w0-ww[j]*ww[j];
              if (ww[j] < w0 ){
                phase[i][j]=Math.atan(2*gamma[i]*ww[j]/den); //The phase difference
              }else{
                //shift by pi needed if w >w0;
                phase[i][j]=Math.PI+Math.atan(2*gamma[i]*ww[j]/den);
              }
           //System.out.println("i="+i+", j="+j+", w="+ww[j]+", c="+cc[i]+
           //                     " ph="+phase[i][j]);
         }
      }
       myControl.println("Values of c, gamma and the phase(w) for a given");
       myControl.println("value of w is given below to identify the curves.");
       for (int i=0; i<=Nc;i++){
         plot.setConnected(i,true);
         plot.setMarkerSize(i, 1);//ith=dataset index, use size 1 for symbol
         plot.append(i, ww, phase[i]);//colors based on DatasetManager.getLineColor(i)
         myControl.println("c="+nf.format(cc[i])+", gamma="+nf.format(gamma[i])+
                           ", phase(w="+nf.format(ww[Nw])+")="+nf.format(phase[i][Nw]));
       }
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the phase difference between the driving force");
     myControl.println ("and the solution of the driven HO");
     myControl.println ("k=spring constant, c=damping coefficient, m=mass");
     myControl.println ("gamma=c/2/m, w0=sqrt(k/m), w=driving frequency");
     myControl.println ("The phase difference has to be shifted beyond w=w0.");
     myControl.println ("if w < w0 the phase=atan(2*gamma*w/(w0^2-w^2);");
     myControl.println ("else the phase=pi+atan(2*gamma*w/(w0^2-w^2))");
     myControl.println ("w range: [wmin,wmax], Nw: w points used, dw=step size");
     myControl.println ("c range: [cmin,cmax], Nc: c points used, dc=step size");
     ww    = null;
     cc    = null;
     gamma = null;
     phase = null;
     m=0.5;
     k=0.5;
     w0=Math.sqrt(k/m);
     wmin=0.0;
     double n=2.5;
     wmax=n*w0;
     Nw=(int)(double)(Math.floor(33*n+1));
     cmin=0.01;
     cmax=1.0; //maximum c for real resonant frequency
     Nc=6;
     myControl.setValue("m",m);
     myControl.setValue("k",k);
     myControl.setValue("wmin",wmin);
     myControl.setValue("wmax",wmax);
     myControl.setValue("Nw",Nw);
     dw=(wmax-wmin)/Nw;
     myControl.setValue("cmin",cmin);
     myControl.setValue("cmax",cmax);
     myControl.setValue("Nc",Nc);
     dc=(cmax-cmin)/Nc;
     myControl.setValue("w0=",w0);
     myControl.setValue("dc((cmax-cmin)/Nc)=",dc);
     myControl.setValue("dw((wmax-wmin)/Nw)=",dw);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new drive_phaseApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
