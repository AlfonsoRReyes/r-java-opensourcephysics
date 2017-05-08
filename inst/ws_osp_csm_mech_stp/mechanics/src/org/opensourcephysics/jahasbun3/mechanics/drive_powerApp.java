/*
@Author J E Hasbun 2007.
Plots the power supplied by the driving force vs frequency for the driven HO.
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


public class drive_powerApp implements Calculation {
   PlotFrame plot= new PlotFrame("w_D","Power(w)","Average Power vs w_D");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double pp[][], ww[], cc[], gamma[],pmax[],wpmax[];
   double den,m,k,cmin,cmax,w0,wmin,wmax,dw,dc,F0,desc,A,ph;
   int Nc,Nw;

   public drive_powerApp(){
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
     F0=myControl.getDouble("F0");
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
     pmax   = new double[Nw+1];
     wpmax  = new double[Nw+1];
     pp     = new double[Nc+1][Nw+1];
     for (int i=0; i<=Nc;i++){
     pmax[i]=0;
     cc[i]=cmin+i*dc;
     gamma[i] = cc[i] / 2 / m;
       for (int j = 0; j <= Nw; j++) {
         if (i==0){ww[j]=wmin+j*dw;} // build ww array once
         den=w0*w0-ww[j]*ww[j];
         if (ww[j] < w0 ){
            ph=Math.atan(2*gamma[i]*ww[j]/den); //The phase difference
         }else{
            //shift by pi needed if w >w0;
            ph=Math.PI+Math.atan(2*gamma[i]*ww[j]/den);
         }
         desc =(2*gamma[i]*ww[j])*(2*gamma[i]*ww[j])+
               (w0*w0-ww[j]*ww[j])*(w0*w0-ww[j]*ww[j]);
         A=F0/m/Math.sqrt(desc);
         pp[i][j]=0.5*F0*A*ww[j]*Math.sin(ph); //the power
         if(pp[i][j]>pmax[i]){
           pmax[i]=pp[i][j];     //store maximum power values
           wpmax[i]=ww[j];       //values of ww at maximum pp
         }
         //System.out.println("i="+i+", j="+j+", w="+ww[j]+", c="+cc[i]+
         //                     " power="+pp[i][j]);
       }
     }
     myControl.println("Values of c, gamma, power maximum, and the value of");
     myControl.println("w at which pmax occurs for easy curve identification.");
     for (int i=0; i<=Nc;i++){
       plot.setConnected(i,true);
       //plot.setMarkerSize(i, 1);//ith=dataset index, use size 1 for symbol
       plot.setMarkerSize(i,0);   // use size zero to have no marker
       plot.append(i, ww, pp[i]);//colors based on DatasetManager.getLineColor(i)
       myControl.println("c="+nf.format(cc[i])+", gamma="+nf.format(gamma[i])+
                         ", pmax(w="+nf.format(wpmax[i])+")="+nf.format(pmax[i]));
       }
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots the power supplied by the driving force vs frequency");
     myControl.println ("k=spring constant, c=damping coefficient, m=mass");
     myControl.println ("gamma=c/2/m, w0=sqrt(k/m), w=driving frequency");
     myControl.println ("p(w)=0.5*F0*A*w*sin(ph), where A=F0/m/sqrt(desc)");
     myControl.println ("with desc=(2*gamma*w)^2+(w0^2-w^2)^2), and where");
     myControl.println ("if w < w0 the ph=atan(2*gamma*w/(w0^2-w^2)");
     myControl.println ("else ph=pi+atan(2*gamma*w/(w0^2-w^2)). Here F0 is the");
     myControl.println ("driving force amplitude.");
     myControl.println ("w range: [wmin,wmax], Nw: w points used, dw=step size");
     myControl.println ("c range: [cmin,cmax], Nc: c points used, dc=step size");
     ww    = null;
     cc    = null;
     gamma = null;
     pp    = null;
     pmax  = null;
     wpmax = null;
     m=0.5;
     k=0.5;
     F0=0.5;
     w0=Math.sqrt(k/m);
     wmin=0.01;
     //double n=2.5;
     wmax=3;
     Nw=200;
     cmin=0.2;
     cmax=1.0; //maximum c for real resonant frequency
     Nc=3;
     myControl.setValue("m",m);
     myControl.setValue("k",k);
     myControl.setValue("F0",F0);
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
     Calculation model = new drive_powerApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
