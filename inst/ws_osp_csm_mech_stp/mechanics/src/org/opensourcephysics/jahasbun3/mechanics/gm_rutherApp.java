/*
@Author J E Hasbun 2007.
Compares the experimental data of Geiger & Marsden (1913) with the Rutherford
Scattering formula for the number of particles.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.frames.*;
import java.text.NumberFormat;

public class gm_rutherApp implements Calculation {
   PlotFrame plot= new PlotFrame("THETA (degrees)","N(THETA)",
                                 "Comparison of Rutherford formula and 1913 Experiment");
   NumberFormat nf = NumberFormat.getInstance();
   private Control myControl;
   double [][] Nth;
   double [] th,thd;
   double v0,vb,m,ma,Ene,za,dth,th0,thmax,K,pi,cf,factor;
   int NPTS;
   double [] cr=new double[]{0.019,0.01};//fit coef for comparison with experiment
   double [] zt=new double[]{47,79};//target charges (47->silver, 79->gold)
   //Actual Data from H. Geiger & E. Marsden, Phil. Mag. Vol.25, 605 (1913)
   //for angle in degrees and number of scintillations obtained for Ag and Au
   double [] deh={15,22.5,30,37.5,45,60,75,105,120,135,150}; // Exp. angle in Degrees from paper
   double [] rth=new double [11];//to hold deh in radians
   double [] rsig1={105400,20300,5260,1760,989,320,136,47.3,33.0,27.4,22.2};//Silver
   double [] rsig2={132000,27300,7800,3300,1435,477,211,69.5,51.9,43.0,33.1};//Gold


   public gm_rutherApp(){
     plot.setLocation(10,5);
     plot.setSize(390,390);
     nf.setMaximumFractionDigits(3);
     pi=Math.PI; cf=pi/180;
   }

   public void calculate() {
     myControl.clearMessages();
     Ene=myControl.getDouble("projectile_energy(eV)");
     m=myControl.getDouble("projectile_mass(m_alpha)");
     za=myControl.getDouble("projectile_z");
     th0=myControl.getDouble("th0");
     thmax=myControl.getDouble("thmax");
     NPTS=myControl.getInt("NPTS");
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth=(thmax-th0)/NPTS)=",dth);
     th     = new double[NPTS];    //angle in radians
     thd    = new double[NPTS];    //angle in degrees
     Nth    = new double[2][NPTS]; //cross-section
     //other values
     v0=Math.sqrt(2*Ene/m/ma)/vb;  //initial speed in units of vb
     //theretical calculation first on silver then on gold
     for (int i=0; i<2; i++){
       K=za*zt[i];                 //dimensionless force constant
       factor=K*K/4/m/m/Math.pow(v0,4); //regular cross-section factor
       for (int j=0; j < NPTS;j++){
         if (i == 0) {
           th[j]=th0+j*dth;
           thd[j]=th[j]/cf;       //convert to degree for plotting
         }
         Nth[i][j]=cr[i]*factor/Math.pow(Math.sin(th[j]/2),4);
       }
       //myControl.println("i="+i+", cr="+cr[i]+" zt="+zt[i]+", factor="+factor);
     }
     myControl.println("Units used: speed="+nf.format(vb)+"c, m_alpha="+
                       nf.format(ma)+"eV;\na_b=1.e-15m, time=1.695e-22s");
     myControl.println("black: silver, blue: gold");
     myControl.println("dots: experimental, lines: Rutherford formula");
     plot.setLogScaleY(true);
     plot.setConnected(0,false);
     plot.setConnected(1,false);
     plot.setConnected(2,true);
     plot.setConnected(3,true);
     plot.setMarkerSize(0, 2);
     plot.setMarkerShape(0,1);
     plot.setMarkerColor(0,java.awt.Color.black);
     plot.append(0, deh, rsig1); //silver - experimental data
     plot.setMarkerSize(1, 2);
     plot.setMarkerShape(1,1);
     plot.setMarkerColor(1,java.awt.Color.blue);
     plot.append(1, deh, rsig2); //gold - experimental data
     plot.setLineColor(2,java.awt.Color.black);
     plot.setMarkerSize(2, 0);
     plot.append(2,thd,Nth[0]);  //silver theory line
     plot.setLineColor(3,java.awt.Color.blue);
     plot.setMarkerSize(3, 0);
     plot.append(3,thd,Nth[1]);  //gold theory line
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Compares the experimental data of Geiger & Marsden (1913)");
     myControl.println ("with the Rutherford Scattering formula for the number of");
     myControl.println ("particles:N(THETA)=c*(K^2/4/m^2/v0^4)/sin(THETA)/2)^4, where");
     myControl.println ("c is our fitting constant depending on the target.");
     th    = null;
     thd   = null;
     Nth   = null;
     m=1;                //projectile mass in units of alpha particle mass
     vb=0.01965;         //velocity units (a_b/tau_b) (in inits of c=light speed)
     ma=3730e6;          //alpha particle mass energy in eV
     Ene=1e6;            //initial projectile energy in eV
     za=2;               //za=alpha particle
     NPTS=50;
     for (int i=0; i< deh.length; i++){
       rth[i]=deh[i]*cf; //convert the experimental angle to radians
     }
     th0=Math.round(1000.*(ArrayLib.min(rth)-0.05))/1000.;
     thmax=Math.round(1000.*(ArrayLib.max(rth)+0.05))/1000.;
     myControl.setValue("projectile_energy(eV)",Ene);
     myControl.setValue("projectile_mass(m_alpha)",m);
     myControl.setValue("projectile_z",za);
     myControl.setValue("th0",th0);
     myControl.setValue("thmax",thmax);
     myControl.setValue("NPTS",NPTS);
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth=(thmax-th0)/NPTS)=",dth);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new gm_rutherApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,390);
     myControl.setDividerLocation(210);
     model.setControl(myControl);
   }
}
