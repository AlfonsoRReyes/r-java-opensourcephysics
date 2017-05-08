/*
@Author J E Hasbun 2007.
Does Rutherford scattering cross-section versus atomic number.
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

public class ruther_cross2App implements Calculation {
  PlotFrame plot0= new PlotFrame("Zt","THETA (rad)",
                                "Alpha-particle scatt for various targets");
  PlotFrame plot1= new PlotFrame("Zt","Sigma(a_{b}^{2})",
                                "Cross-Section vs Zt");
  PlotFrame plot2= new PlotFrame("THETA (rad)","Sigma(a_{b}^{2})",
                                "Cross-Section vs THETA");
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df= new DecimalFormat("+0.000E00");
  private Control myControl;
  double [] th, sigma,zt;
  double v0,vb,m,ma,Ene,za,dth,K,pi,b,sK,ye;
  int NPTS,ztmin,ztmax,dzt;
  int maxiter=25;
  double tolerance=5e-3;
  Bisnewt bisnewt=new Bisnewt();
  Function fc;

   public ruther_cross2App(){
     plot0.setLocation(5,5);
     plot0.setSize(360, 230);
     plot1.setLocation(370,5);
     plot1.setSize(360, 230);
     plot2.setLocation(5,240);
     plot2.setSize(360,230);
     nf.setMaximumFractionDigits(3);
     pi=Math.PI;
   }

   public void calculate() {
     myControl.clearMessages();
     fc=new fcos();
     m=myControl.getDouble("projectile_mass(m_alpha)");
     Ene=myControl.getDouble("projectile_energy(eV)");
     b=myControl.getDouble("impact_parameter(a_b)");
     za=myControl.getDouble("projectile_z");
     ztmin=myControl.getInt("target_z_min");
     ztmax=myControl.getInt("target_z_max");
     dzt=myControl.getInt("dzt");
     NPTS=ztmax-ztmin+1;
     NPTS=myControl.getInt("NPTS");
     th     = new double[NPTS];    //angle
     sigma  = new double[NPTS];    //cross-section
     zt     = new double[NPTS];    //cross-section
     //other values
     v0=Math.sqrt(2*Ene/m/ma)/vb;  //initial speed in units of vb
     for (int i=0; i < NPTS;i++){
       zt[i]=i+ztmin;          //target array
       K=za*zt[i];             //dimensionless force constant
       if(K!=0){
         sK=K/Math.abs(K);
         ye=-sK*Math.sqrt(1+m*m*Math.pow(v0,4)*b*b/K/K); //eccentricity
       }else{ye=0;}
       th[i]=bisnewt.Bisnewt(1.e-5,pi/2, maxiter, tolerance,fc);//find alpha asymptote angle
       th[i]=pi-2*th[i];      //total scattering angle produced by target species
       sigma[i]=K*K*2*pi*Math.sin(th[i])/Math.pow(v0*Math.sin(th[i]/2.),4)/4/m/m;
       //myControl.println("i="+i+", zt="+zt[i]+", ye="+ye+" sigma="+sigma[i]);
     }
     myControl.println("The units employed follow.");
     myControl.println("speed="+nf.format(vb)+"c, m_alpha="+
                       df.format(ma)+"eV;\na_b=1.e-15m, time=1.695e-22s");
     plot0.setConnected(0,true);
     plot0.setPreferredMinMax(0,ztmax,0,ArrayLib.max(th));
     plot0.setMarkerSize(0, 0);  //ith=dataset index, use size 0 for symbol
     plot0.append(0, zt, th);
     plot1.setConnected(0,true);
     plot1.setPreferredMinMax(0,ztmax,0,ArrayLib.max(sigma));
     plot1.setMarkerSize(0, 0);  //ith=dataset index, use size 0 for symbol
     plot1.append(0, zt, sigma);
     plot2.setConnected(0,true);
     plot2.setPreferredMinMax(0.0,ArrayLib.max(th),0,ArrayLib.max(sigma));
     plot2.setMarkerSize(0, 0);  //ith=dataset index, use size 0 for symbol
     plot2.append(0, th, sigma); //colors based on DatasetManager.getLineColor(i)
   }

   public class fcos implements Function {
     public double evaluate(double x) {
       //defines f(x) whose value vanishes at x=root of f(x)
       double f;
       f=1+ye*Math.cos(x);
       return f;
     }
   }

   public void clear  () {
     plot0.clearData();
     plot1.clearData();
     plot2.clearData();
   }

   public void resetCalculation() {
     clear();
     myControl.clearMessages();
     myControl.println ("Does scattering cross-section versus atomic number");
     myControl.println ("We use the Rutherford Scattering formula:");
     myControl.println ("sigma=K^2*2*pi*sin(theta)/sin(theta/2).^4/4/m^2/v0^4)");
     th    = null;
     sigma = null;
     m=1;                //projectile mass in units of alpha particle mass
     vb=0.01965;         //velocity units (a_b/tau_b) (in inits of c=light speed)
     ma=3730e6;          //alpha particle mass energy in eV
     Ene=5e6;            //initial projectile energy in eV
     b=20;               //impact parameter
     za=2;               //za=projectile (2->alpha)
     ztmin=13; ztmax=79; //from aluminum to gold
     NPTS=ztmax-ztmin+1;
     dzt=1;
     myControl.setValue("projectile_mass(m_alpha)",m);
     myControl.setValue("projectile_energy(eV)",df.format(Ene));
     myControl.setValue("impact_parameter(a_b)",b);
     myControl.setValue("projectile_z",za);
     myControl.setValue("target_z_min",ztmin);
     myControl.setValue("target_z_max",ztmax);
     myControl.setValue("dzt",dzt);
     myControl.setValue("NPTS",NPTS);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
   }

   public static void main(String[] args) {
     Calculation model = new ruther_cross2App();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(370, 240);
     myControl.setSize(360,365);
     myControl.setDividerLocation(210);
     model.setControl(myControl);
   }
}
