//============================================================================
//  
// ODE sample :  ellipse
//  
//----------------------------------------------------------------------------
//  > gcc -g -o main00 main.c -lm
//  > ./main00
//============================================================================

#include <stdio.h>
#include <math.h> 
//#include "header.h"

//----------------------------------------------------
#define PI  3.1415926535897932384626433832795028841971
#define PI2 (2.0*PI)
#define PIH (0.5*PI)
#define EXP 2.718281828459
//----------------------------------------------------
#define FD double
//----------------------------------------------------
#define SOLVER (5)
//(0) Euler (1) Heun (2) Euler(2nd) (3) Runge-Kutta(3rd)
//(4) Runge-Kutta(4th) (5) Runge-Kutta-Gill
//----------------------------------------------------
void exact_ellipse(FD ph,FD c[],FD *u);
void diff_ellipse(FD x,FD y[],FD dydx[],FD c[],int n);
void diff(FD x,FD y[],FD dydx[],FD c[],int n);
//----------------------------------------------------
void euler(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	   void (*diff)(FD,FD [],FD [],FD [],int));
void heun(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	  void (*diff)(FD,FD [],FD [],FD [],int));
void euler2(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	    void (*diff)(FD,FD [],FD [],FD [],int));
void rungekutta3(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
		 void (*diff)(FD,FD [],FD [],FD [],int));
void rungekutta4(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
		 void (*diff)(FD,FD [],FD [],FD [],int));
void rkg(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	 void (*diff)(FD,FD [],FD [],FD [],int));
//----------------------------------------------------

int main(int argc, char *argv[])
{
  char outf[64];
  sprintf(outf, "result001.dat");
  FILE *fpo = fopen(outf,"w");
  
  FD dx=2.e-1;
  FD xmax=8.0;
  FD x=0.0;

  int i;

  FD e=0.4;
  FD a=10.0;
  FD l2 = 0.5*a*(1.0-e*e);
  
  FD c[2];
  c[0] = e; 
  c[1] = a;
    
  int n=2;
  FD y[n];
  
  y[0]=(1.0+e)*0.5/l2;
  y[1]=0.0;
  
  FD xnext;
  FD ynext[n];
  FD yexact;

  exact_ellipse(x,c,&yexact);
  
  FD xx = cos(x)/y[0];
  FD yy = sin(x)/y[0];
  FD xxexact = cos(x)/yexact;
  FD yyexact = sin(x)/yexact;
  
  fprintf(fpo,"%f %f %f %f %f %f %f %f \n",
	  x,yexact,y[0],y[1],xx,yy,xxexact,yyexact);
  
  while(x<=xmax){

    if(SOLVER==0) euler(x,y,dx,&xnext,ynext,c,n,diff_ellipse);
    if(SOLVER==1) heun(x,y,dx,&xnext,ynext,c,n,diff_ellipse);
    if(SOLVER==2) euler2(x,y,dx,&xnext,ynext,c,n,diff_ellipse);
    if(SOLVER==3) rungekutta3(x,y,dx,&xnext,ynext,c,n,diff_ellipse);
    if(SOLVER==4) rungekutta4(x,y,dx,&xnext,ynext,c,n,diff_ellipse);
    if(SOLVER==5) rkg(x,y,dx,&xnext,ynext,c,n,diff_ellipse);

    x=xnext;
    for(i=0;i<n;i++)  y[i] = ynext[i]; 

    exact_ellipse(x,c,&yexact);
    
    xx = cos(x)/y[0];
    yy = sin(x)/y[0];
    xxexact = cos(x)/yexact;
    yyexact = sin(x)/yexact;
  
    fprintf(fpo,"%f %f %f %f %f %f %f %f \n",
	    x,yexact,y[0],y[1],xx,yy,xxexact,yyexact);
  }
  
  fclose(fpo);
        
  return 0;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void exact_ellipse(FD ph,FD c[],FD *u){

  FD e=c[0];
  FD a=c[1];
  FD l2 = 0.5*a*(1.0-e*e);
  
  *u=(1.0+e*cos(ph))*0.5/l2;
}

void diff_ellipse(FD x,FD y[],FD dydx[],FD c[],int n){

  FD e=c[0];
  FD a=c[1];
  FD l2 = 0.5*a*(1.0-e*e);

  dydx[0] = y[1];
  dydx[1] = -y[0]+0.5/l2;
}

void diff(FD x,FD y[],FD dydx[],FD c[],int n){
  
  dydx[0] = y[1];
  dydx[1] = -4.0*y[0];
}
//----------------------------------------------------------------------------
void euler(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	   void (*diff)(FD,FD [],FD [],FD [],int)){

  int i;

  FD dydx[n];
  
  (*diff)(x,y,dydx,c,n);
  
  x += dx;

  for(i=0;i<n;i++)  y[i] += dydx[i]*dx;
  
  *xnext = x;

  for(i=0;i<n;i++)  ynext[i] = y[i]; 
}

void heun(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	  void (*diff)(FD,FD [],FD [],FD [],int)){

  int i;

  FD dydx1[n];
  
  (*diff)(x,y,dydx1,c,n);
  
  FD x2;
  FD y2[n];
  FD dydx2[n];

  x2 = x+dx;
  for(i=0;i<n;i++)  y2[i] = y[i]+dydx1[i]*dx;

  (*diff)(x2,y2,dydx2,c,n);
  
  x += dx;
  for(i=0;i<n;i++)  y[i] += (dydx1[i]+dydx2[i])*0.5*dx;
  
  *xnext = x;
  for(i=0;i<n;i++)  ynext[i] = y[i]; 
}

void euler2(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	    void (*diff)(FD,FD [],FD [],FD [],int)){

  int i;

  FD dydx1[n];
  
  (*diff)(x,y,dydx1,c,n);
  
  FD x2;
  FD y2[n];
  FD dydx2[n];
  FD dx2=0.5*dx;

  x2 = x+dx2;
  for(i=0;i<n;i++)  y2[i] = y[i]+dydx1[i]*dx2;

  (*diff)(x2,y2,dydx2,c,n);
  
  x += dx;
  for(i=0;i<n;i++)  y[i] += dydx2[i]*dx;
  
  *xnext = x;
  for(i=0;i<n;i++)  ynext[i] = y[i]; 
}

void rungekutta3(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
		 void (*diff)(FD,FD [],FD [],FD [],int)){

  int i;
  
  FD dydx1[n];
  
  (*diff)(x,y,dydx1,c,n);

  FD x2;
  FD y2[n];
  FD dydx2[n];
  FD dx2=0.5*dx;

  x2 = x+dx2;
  for(i=0;i<n;i++)  y2[i] = y[i]+dydx1[i]*dx2;

  (*diff)(x2,y2,dydx2,c,n);

  FD x3;
  FD y3[n];
  FD dydx3[n];

  x3 = x+dx;
  for(i=0;i<n;i++)  y3[i] = y[i]-dydx1[i]*dx+2.0*dydx2[i]*dx;

  (*diff)(x3,y3,dydx3,c,n);

  FD dx6=dx/6.0;
  
  x += dx;
  for(i=0;i<n;i++)  y[i] += (dydx1[i]+4.0*dydx2[i]+dydx3[i])*dx6;
  
  *xnext = x;
  for(i=0;i<n;i++)  ynext[i] = y[i]; 
}

void rungekutta4(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
		 void (*diff)(FD,FD [],FD [],FD [],int)){

  int i;
  
  FD dydx1[n];
  
  (*diff)(x,y,dydx1,c,n);

  FD x2;
  FD y2[n];
  FD dydx2[n];
  FD dx2=0.5*dx;

  x2 = x+dx2;
  for(i=0;i<n;i++)  y2[i] = y[i]+dydx1[i]*dx2;

  (*diff)(x2,y2,dydx2,c,n);

  FD x3;
  FD y3[n];
  FD dydx3[n];

  x3 = x+dx2;
  for(i=0;i<n;i++)  y3[i] = y[i]+dydx2[i]*dx2;

  (*diff)(x3,y3,dydx3,c,n);

  FD x4;
  FD y4[n];
  FD dydx4[n];

  x4 = x+dx;
  for(i=0;i<n;i++)  y4[i] = y[i]+dydx3[i]*dx;

  (*diff)(x4,y4,dydx4,c,n);
  
  FD e1 = 1.0/6.0;
  FD e2 = 2.0/6.0;
  FD e3 = 2.0/6.0;
  FD e4 = 1.0/6.0;
  
  x += dx;
  for(i=0;i<n;i++)
    y[i] += (e1*dydx1[i] +e2*dydx2[i] +e3*dydx3[i] +e4*dydx4[i]) * dx;
  
  *xnext = x;
  for(i=0;i<n;i++)  ynext[i] = y[i]; 
}

void rkg(FD x,FD y[],FD dx,FD *xnext,FD ynext[],FD c[],int n,
	 void (*diff)(FD,FD [],FD [],FD [],int)){

  int i;
  
  FD dydx1[n];
  
  (*diff)(x,y,dydx1,c,n);

  FD x2;
  FD y2[n];
  FD dydx2[n];
  FD dx2=0.5*dx;

  x2 = x+dx2;
  for(i=0;i<n;i++)  y2[i] = y[i]+dydx1[i]*dx2;

  (*diff)(x2,y2,dydx2,c,n);

  FD x3;
  FD y3[n];
  FD dydx3[n];
  FD c1 = (-1.0+sqrt(2.0))*0.5;
  FD c2 = 1.0-1.0/sqrt(2.0);

  x3 = x + dx2;
  for(i=0;i<n;i++)  y3[i] = y[i] +c1*dydx1[i]*dx +c2*dydx2[i]*dx;

  (*diff)(x3,y3,dydx3,c,n);

  FD x4;
  FD y4[n];
  FD dydx4[n];
  FD d2 = -1.0/sqrt(2.0);
  FD d3 = 1.0+1.0/sqrt(2.0);

  x4 = x + dx;
  for(i=0;i<n;i++)  y4[i] = y[i] +d2*dydx2[i]*dx +d3*dydx3[i]*dx;

  (*diff)(x4,y4,dydx4,c,n);

  FD e1 = 1.0/6.0;
  FD e2 = (2.0-sqrt(2.0))/6.0;
  FD e3 = (2.0+sqrt(2.0))/6.0;
  FD e4 = 1.0/6.0;
  
  x += dx;
  for(i=0;i<n;i++)
    y[i] += (e1*dydx1[i] +e2*dydx2[i] +e3*dydx3[i] +e4*dydx4[i]) * dx;
  
  *xnext = x;
  for(i=0;i<n;i++)  ynext[i] = y[i]; 
}

