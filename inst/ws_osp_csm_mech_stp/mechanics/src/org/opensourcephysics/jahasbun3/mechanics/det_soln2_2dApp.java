/*
@Author J E Hasbun 2007.
Uses cartesian coordinates to finds a rectangle's inertia tensor numerically,
the inertia moments and the principal axes directions.
@Copyright (c) 2007
This software is to support Intermediate Classical Mechanics
with MATLAB Applications by J. E. Hasbun using the Open Source Physics library
http://www.opensourcephysics.org under the terms of the GNU General Public
License (GPL) as published by the Free Software Foundation.
*/

package org.opensourcephysics.jahasbun3.mechanics;
import org.opensourcephysics.controls.*;
import java.text.*;

public class det_soln2_2dApp implements Calculation {
  EigenJacobi eigenjacobi=new EigenJacobi();
  MatPrint mprint=new MatPrint();
  private Control myControl;
  int NumColRow=3, Nt=6;
  int icmax,Nx,Ny;
  double tol;
  double [][] A,DiagD,V, data;
  double a,b,M,rho,ax,ay,bx,by,dx,dy;
  double Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz;
  NumberFormat nf = NumberFormat.getInstance();
  DecimalFormat df = new DecimalFormat("+0.00000E00");

  public det_soln2_2dApp(){
    nf.setMinimumFractionDigits(3);
  }

   public void calculate() {
     clear();
     a=myControl.getDouble("a");
     b=myControl.getDouble("b");
     M=myControl.getDouble("M");
     ax=myControl.getDouble("ax");
     ay=myControl.getDouble("ay");
     tol=myControl.getDouble("tolerance");
     icmax=myControl.getInt("max_iterations");
     Nx=myControl.getInt("x-integration points");
     Ny=myControl.getInt("y-integration points");
     if(Nx%2.==0){Nx=Nx+1;} //make sure Nx, Ny are odd integers
     if(Ny%2.==0){Ny=Ny+1;}
     rho=M/(a*b);
     bx=a+ax; dx=(bx-ax)/(Nx-1.);
     by=b+ay; dy=(by-ay)/(Ny-1.);
     double [][] fx=new double [Nt][Nx];
     double [][] fy=new double [Nt][Ny];
     double [] x=new double [Nx];
     double [] y=new double [Ny];
     //x coordinates
     for(int i=0; i<Nx; i++){
       x[i]=ax+i*dx;
     }
     //y coordinates
     for(int i=0; i<Ny; i++){
       y[i]=ay+i*dy;
     }
     //create the inertia tensor integrands and perform the x-integration
     for (int j=0; j<Ny; j++){
       for (int i=0; i<Nx; i++){
             for (int k=0; k<Nt; k++){
                 fx[k][i]=rho*Inert_el.inert_el(k,x[i],y[j],0);//z=0 in 2d
             }
       }
         for (int k=0; k<Nt; k++){
             fy[k][j]=Simp.Simp(fx[k],dx); //Simpson's rule
         }
     }
     //finally integrate over the y coord to get the moments
     double ff=0.;
     for (int m=0; m<Nt; m++){
       ff=Simp.Simp(fy[m],dy); //Simpson's rule
         if (m==0){Ixx=ff;}
         if (m==1){Ixy=ff;}
         if (m==2){Ixz=ff;}
         if (m==3){Iyy=ff;}
         if (m==4){Iyz=ff;}
         if (m==5){Izz=ff;}
     }
     Iyx=Ixy; Izx=Ixz; Izy=Iyz;//use symmetry for the rest of the tensor elements
     double [][] A={{Ixx,Ixy,Ixz},{Iyx,Iyy,Iyz},{Izx,Izy,Izz}}; //Inertia tensor
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
     myControl.println("The eigenvalues are the principal moments of inertia and");
     myControl.println("the eigenvectors are their corresponding direction cosines");
     myControl.println ("The solid density is rho=M/(a*b)="+nf.format(rho));
     }

   public void clear  () {
     myControl.clearMessages();
   }

   public void resetCalculation() {
     clear();
     myControl.println ("Uses cartesian coordinates to finds a rectangle's");
     myControl.println ("inertia tensor numerically, the inertia moments and");
     myControl.println ("the principal axes directions. Here a, b are the plate's");
     myControl.println ("side lengths, M is its mass, and ax, ay are the origin of");
     myControl.println ("coordinates.");
     icmax=100;
     tol=1.e-5;
     a=1.0; b=1.0;             //solid plates sides
     M=1.0;                    //solid mass and density
     ax=0.; ay=0.;             //origin at corner as lower limits
     Nx=11; Ny=11;
     myControl.setValue("a",a);
     myControl.setValue("b",b);
     myControl.setValue("M",M);
     myControl.setValue("ax",ax);
     myControl.setValue("ay",ay);
     myControl.setValue("tolerance",tol);
     myControl.setValue("max_iterations",icmax);
     myControl.setValue("x-integration points",Nx);
     myControl.setValue("y-integration points",Ny);
   }

   public void setControl(Control control) {
     myControl = control;
     resetCalculation();
     mprint.myControl=control;
   }

   public static void main(String[] args) {
     Calculation model = new det_soln2_2dApp();
     CalculationControl myControl;
     myControl = new CalculationControl(model);
     myControl.setLocation(50, 5);
     myControl.setSize(500,590);
     myControl.setDividerLocation(225);
     model.setControl(myControl);
   }
}