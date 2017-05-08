/*
@Author J E Hasbun 2007.
Plots the Rutherford scattering cross-section versus scattering angle.
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
import java.text.*;

public class ruther_cross1App implements Calculation {
   PlotFrame plot= new PlotFrame("THETA (rad)","Sigma(a_{b}^{2})",
                                 "Rutherford Scattering Cross-Section");
   NumberFormat nf = NumberFormat.getInstance();
   DecimalFormat df= new DecimalFormat("+0.000E00");
   private Control myControl;
   double [] th, sigma;
   double v0,vb,m,ma,Ene,za,zt,dth,th0,thmax,K,pi;
   int NPTS;

   public ruther_cross1App(){
     plot.setLocation(10,5);
     plot.setSize(390,390);
     nf.setMaximumFractionDigits(3);
     pi=Math.PI;
   }

   public void calculate() {
     myControl.clearMessages();
     m=myControl.getDouble("projectile_mass(m_alpha)");
     Ene=myControl.getDouble("projectile_energy(eV)");
     za=myControl.getDouble("projectile_z");
     zt=myControl.getDouble("target_z");
     th0=myControl.getDouble("th0");
     thmax=myControl.getDouble("thmax");
     NPTS=myControl.getInt("NPTS");
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth=(thmax-th0)/NPTS)=",nf.format(dth));
     th     = new double[NPTS];    //angle
     sigma  = new double[NPTS];    //cross-section
     //other values
     v0=Math.sqrt(2*Ene/m/ma)/vb;  //initial speed in units of vb
     K=za*zt;                      //dimensionless force constant
     for (int i=0; i < NPTS;i++){
       th[i]=th0+i*dth;
       sigma[i]=K*K*2*pi*Math.sin(th[i])/Math.pow(v0*Math.sin(th[i]/2.),4)/4/m/m;
       //myControl.println("i="+i+", th="+th[i]+" sigma="+sigma[i]);
     }
     myControl.println("The units employed follow.");
     myControl.println("speed="+nf.format(vb)+"c, m_alpha="+
                       df.format(ma)+"eV;\na_b=1.e-15m, time=1.695e-22s");
     plot.setConnected(0,true);
     plot.setPreferredMinMax(0.0,thmax,0.0,ArrayLib.max(sigma));
     plot.setMarkerSize(0, 0);  //ith=dataset index, use size 0 for symbol
     plot.append(0, th, sigma); //colors based on DatasetManager.getLineColor(i)
   }

   public void clear  () {
     plot.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Plots scattering cross-section versus scattering angle.");
     myControl.println ("We use the Rutherford Scattering formula:");
     myControl.println ("sigma=K^2*2*pi*sin(theta)/sin(theta/2)^4/4/m^2/v0^4)");
     th    = null;
     sigma = null;
     m=1;                //projectile mass in units of alpha particle mass
     vb=0.01965;         //velocity units (a_b/tau_b) (in inits of c=light speed)
     ma=3730e6;          //alpha particle mass energy in eV
     Ene=5e6;            //initial projectile energy in eV
     za=2; zt=79;        //za=projectile, zt=target charges (2->alpha, 79->gold)
     NPTS=50;
     th0=0.4; thmax=pi;
     myControl.setValue("projectile_mass(m_alpha)",m);
     myControl.setValue("projectile_energy(eV)",df.format(Ene));
     myControl.setValue("projectile_z",za);
     myControl.setValue("target_z",zt);
     myControl.setValue("th0",th0);
     myControl.setValue("thmax",nf.format(thmax));
     myControl.setValue("NPTS",NPTS);
     dth=(thmax-th0)/NPTS;
     myControl.setValue("dth=(thmax-th0)/NPTS)=",nf.format(dth));
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new ruther_cross1App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(405, 5);
     myControl.setSize(390,390);
     myControl.setDividerLocation(210);
     model.setControl(myControl);
   }
}
