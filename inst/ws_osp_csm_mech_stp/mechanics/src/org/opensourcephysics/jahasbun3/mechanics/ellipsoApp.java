/*
@Author J E Hasbun 2007.
Calculates an ellipsoid inertia tensor & mass numerically. Also plots it in 3D.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.display3d.simple3d.*;
import java.awt.Color;
import java.text.*;

public class ellipsoApp implements Calculation {
  Display3DFrame frame3d = new Display3DFrame("Ellipsoid");
  MatPrint mprint=new MatPrint();
  private Control myControl;
  int Nt=7, NumColRow=3;
  int Nx,Ny,Nz;
  double a,b,c,M,rho,ax,ay,az,bx,by,bz,dx,dy,dz,pi;
  double Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz,Mass;
  ElementEllipsoid ellipso=new ElementEllipsoid();
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public ellipsoApp(){
    frame3d.setDecorationType(VisualizationHints.DECORATION_AXES);
    frame3d.setLocation(5,5);
    frame3d.setSize(400,450);
    nf.setMinimumFractionDigits(4);
    pi=Math.PI;
  }

   public void calculate() {
     clear();
     a=myControl.getDouble("a");
     b=myControl.getDouble("b");
     c=myControl.getDouble("c");
     rho=myControl.getDouble("rho");
     Nx=myControl.getInt("x-integration points");
     Ny=myControl.getInt("y-integration points");
     Nz=myControl.getInt("z-integration points");
     if(Nx%2.==0){        //make sure Nx, Ny, Nz are odd integers
        Nx=Nx+1;
        myControl.setValue("x-integration points",Nx);
     }
     if(Ny%2.==0){
       Ny=Ny+1;
       myControl.setValue("y-integration points",Ny);
     }
     if(Nz%2.==0){
       Nz=Nz+1;
       myControl.setValue("z-integration points",Nz);
     }
     double [][] fx=new double [Nt][Nx];
     double [][] fy=new double [Nt][Ny];
     double [][] fz=new double [Nt][Nz];
     double [] x=new double [Nx];
     double [] y=new double [Ny];
     double [] z=new double [Nz];
     //create the inertia tensor integrands and perform the x-integration
     ax=-a; bx=-ax;                                 //x lower, upper limit
     dx=(bx-ax)/(Nx-1);                             //x spacing
     //x coordinates
     makegrid(ax,x,dx,Nx);
     double tmp;
     for(int k=0;k<Nx; k++){                        //x loop
       tmp=1.-x[k]*x[k]/a/a;
       if(tmp<0){tmp=0;}
       ay=-b*Math.sqrt(tmp);                        //y lower limit,
       by=-ay;                                      //y upper limit,
       dy=(by-ay)/(Ny-1);                           //y spacing
       makegrid(ay,y,dy,Ny);                        //y grid
       for(int j=0; j<Ny; j++){                     //y loop
          tmp=1.-x[k]*x[k]/a/a-y[j]*y[j]/b/b;
          if(tmp<0){tmp=0;}
          az=-c*Math.sqrt(tmp);                     //z lower limit
          bz=-az;                                   //z upper limit
          dz=(bz-az)/(Nz-1);                        //z spacing
          makegrid(az,z,dz,Nz);                     //z grid
          for(int i=0; i<Nz; i++){                  //z loop
             for(int m=0; m<Nt; m++){
               fz[m][i]=rho*Inert_el2.inert_el2(m,x[k],y[j],z[i]);
             }
          }                                         //end z loop
          for(int m=0; m<Nt; m++){
            fy[m][j]=Simp.Simp(fz[m],dz);           //Simpson rule
          }
       }                                            //end y loop
       for(int m=0; m<Nt; m++){
         fx[m][k]=Simp.Simp(fy[m],dy);              //Simpson rule
       }
     }                                              //end x loop
     //finally integrate over the x coord to get the moments
     double ff=0.;
     for (int m=0;m<Nt; m++){
         ff=Simp.Simp(fx[m],dx);                   //Simpson rule
         if(m==0){Ixx=ff;}
         if(m==1){Ixy=ff;}
         if(m==2){Ixz=ff;}
         if(m==3){Iyy=ff;}
         if(m==4){Iyz=ff;}
         if(m==5){Izz=ff;}
         if(m==6){Mass=ff;}
     }
     Iyx=Ixy; Izx=Ixz; Izy=Iyz;//use symmetry for the rest of the tensor elements
     double [][] A={{Ixx,Ixy,Ixz},{Iyx,Iyy,Iyz},{Izx,Izy,Izz}}; //Inertia tensor
     M=rho*4.*pi*a*b*c/3.;
     double p_e=(Mass-M)*100./M;                  //error on the mass
     mprint.ContPrintMat2D(A,NumColRow,NumColRow,"Inertia tensor");

     myControl.println("The solid's mass is Mass="+nf.format(Mass)+
                        ", % mass error="+nf.format(p_e));
     myControl.println("Feel free to mouse-rotate the ellipsoid image.");
     ellipso.setSizeXYZ(a,b,c);
     int red=255, green=0, blue=0, transparency=20; //comes out light pink
     Color col=new Color(red,green,blue,transparency);
     ellipso.getStyle().setFillColor(col);
     ellipso.getStyle().setResolution(new Resolution(3, 24, 24));
     frame3d.setAltitude(0.2);
     frame3d.setAzimuth(0.5);
     frame3d.addElement(ellipso);
     frame3d.setPreferredMinMax(-0.45*a,0.45*a,-0.45*b,0.45*b,-0.45*c,0.45*c);
     frame3d.getCamera().setXYZ(1.25*a,1.25*b,1.25*c);
     frame3d.render();
  }

  public void makegrid (double v0, double v[], double dv, int N) {
    for(int i=0; i<N; i++){
      v[i]=v0+i*dv;
    }
  }

   public void clear  () {
     myControl.clearMessages();
     ellipso.setSizeXYZ(0,0,0);
     frame3d.addElement(ellipso);
     frame3d.render();
   }

   public void resetCalculation() {
     clear();
     myControl.println ("Calculates an ellipsoid inertia tensor & mass numerically.");
     myControl.println ("Also plots it in 3D. The ellipsoid is given by: ");
     myControl.println ("(x-x0)^2/a^2+(y-y0)^2/b^2+(z-z)^2/c^2=1 with density rho.");
     myControl.println ("See text for the tensor formulas and integration limits.");
     a=3.0; b=2.0; c=1.0;            //semimajor axes
     rho=1./8.;                      //density
     Nx=35; Ny=35; Nz=35;            //x,y,z points
     myControl.setValue("a",a);
     myControl.setValue("b",b);
     myControl.setValue("c",c);
     myControl.setValue("rho",rho);
     myControl.setValue("x-integration points",Nx);
     myControl.setValue("y-integration points",Ny);
     myControl.setValue("z-integration points",Nz);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
     mprint.myControl=control;
   }

   public static void main(String[] args) {
     Calculation model = new ellipsoApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(410, 5);
     myControl.setSize(375,450);
     myControl.setDividerLocation(200);
     model.setControl(myControl);
   }
}