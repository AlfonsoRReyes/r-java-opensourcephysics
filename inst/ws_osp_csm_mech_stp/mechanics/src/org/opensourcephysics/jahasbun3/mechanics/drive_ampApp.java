/*
@Author J E Hasbun 2007.
Plots the amplitude of the driven HO solution.
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
import java.awt.Color;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;


public class drive_ampApp implements Calculation {
   PlotFrame plot= new PlotFrame("w_D","A(w)","A versus w_D");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double AA[][], ww[], cc[], Amax[], om_res[], gamma[];
   double desc,m,k,cmin,cmax,F0,w0,wmin,wmax,dw,dc;
   int Nc,Nw;

   public drive_ampApp(){
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
     F0=myControl.getDouble("F0");
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
     om_res = new double[Nc+1];
     cc     = new double[Nc+1];
     Amax   = new double[Nc+1];
     gamma    = new double[Nc+1];
     AA     = new double[Nc+1][Nw+1];
     for (int i=0; i<=Nc;i++){
       cc[i]=cmin+i*dc;
       gamma[i] = cc[i] / 2 / m;
       om_res[i]=Math.sqrt(w0*w0-2*gamma[i]*gamma[i]); //resonant frequency
       Amax[i]=F0/2/m/gamma[i]/Math.sqrt(w0*w0-gamma[i]*gamma[i]);
         for (int j = 0; j <= Nw; j++) {
           if (i==0){ww[j]=wmin+j*dw;} // build ww array once
           desc =(2*gamma[i]*ww[j])*(2*gamma[i]*ww[j])+
                 (w0*w0-ww[j]*ww[j])*(w0*w0-ww[j]*ww[j]);
           AA[i][j] =F0/m/Math.sqrt(desc); // the amplitude vs ww
           //System.out.println("i="+i+", j="+j+", w="+ww[j]+", c="+cc[i]+
           //                     " A="+AA[i][j]);
         }
      }
       myControl.println("Below are: c, gamma, resonant frequencies, maximum amplitudes.");
       myControl.println("The graph's green circles are the resonant frequency positions");
       for (int i=0; i<=Nc;i++){
         plot.setMarkerSize(i, 1);
         plot.append(i, ww, AA[i]);
         myControl.println("c="+nf.format(cc[i])+", gamma="+nf.format(gamma[i])+
                           ", w_res="+nf.format(om_res[i])+", Amax="+
                           nf.format(Amax[i]));
       }
       //frequencies at which max amplitudes occur
        plot.setMarkerSize(Nc+1,3); //Nc+1=dataset index, use size 3 for symbol
        plot.setMarkerShape(Nc+1,1);//circle
        plot.setMarkerColor(Nc+1,Color.green,Color.red);
        plot.append(Nc+1,om_res,Amax);//Max amplitude positions om_res
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Amplitude of the driven HO solution versus driving frequency");
     myControl.println ("k=spring constant, c=damping coefficient, m=mass");
     myControl.println ("gamma=c/2/m, w0=sqrt(k/m), F0=driving force amplitude");
     myControl.println ("Amplitude: A=F0/2/m/gamma/sqrt(desc), w=driving frequency");
     myControl.println ("desc =(2*gamma*w)^2+(w0^2-w^2)*(w0^2-w^2)");
     myControl.println ("Resonant driving frequency: om_res=sqrt(w0^2-2*gamma^2)");
     myControl.println ("Max amplitude: Amax=F0/2/m/gamma/sqrt(w0^2-gamma^2)");
     myControl.println ("w range: [wmin,wmax], Nw: w points used, dw=step size");
     myControl.println ("c range: [cmin,cmax], Nc: c points used, dc=step size");
     ww    = null;
     AA    = null;
     cc    = null;
     om_res= null;
     Amax  = null;
     gamma = null;
     m=0.5;
     k=0.5;
     w0=Math.sqrt(k/m);
     F0=0.5;
     wmin=0.1;
     wmax=2;
     Nw=200;
     cmin=0.2;
     cmax=2*m*w0/Math.sqrt(2); //maximum c for real resonant frequency
     Nc=6;
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
     Calculation model = new drive_ampApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,575);
     myControl.setDividerLocation(280);
     model.setControl(myControl);
   }
}
