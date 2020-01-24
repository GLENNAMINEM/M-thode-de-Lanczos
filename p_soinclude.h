#include <math.h>
#include <stdio.h>
#define ITERMAX 100000
#define NMAX 4
#include <mpi.h>


double norme_mat(double a[NMAX][NMAX],int n)
{
 double s,r;
 int i,j;
 r=0;
 for(i=0;i<n;i++)
 {
  s=0;
  for(j=0;j<n;j++)
  s += fabs(a[i][j]);
  if(s>r) r=s;
 }
 return(r);
}
double norme_vect(double x[NMAX],int n)
{
 double t;
 int i;
 t=0;
 for(i=0;i<n;i++)
 t += x[i]*x[i];
 t=sqrt(t);
 return(t);
}

void prod_mat(double a[NMAX][NMAX],double b[NMAX][NMAX],double c[NMAX][NMAX],int n)
{
 double s;
 int i,j,k;
 for(i=0;i<n;i++)
 for(j=0;j<n;j++)//pour le stop
 {
  s=0;
  for(k=0;k<n;k++)
  s += a[i][k]*b[k][j];
  c[i][j]=s;
 }
}

void prod_mat_vect(double a[NMAX][NMAX],int n,double x[NMAX],double y[NMAX])
{
 int i,j;
 for(i=0;i<n;i++)
 {
  y[i]=0;
  for(j=0;j<n;j++)
  y[i] += a[i][j]*x[j];
 }
}
double ps(double x[NMAX],double y[NMAX],int n)
{
 int i;
 double s;
 s=0.0;
 for(i=0;i<n;i++)
 s+=x[i]*y[i];
 return(s);
}


void somme_vect(double q1[NMAX],double q2[NMAX],double q3[NMAX],int n){


 for(int i=0;i<n;i++){

    q3[i] = q1[i] + q2[i];	
	}

}


void vect_sca(double q[NMAX], double q3[NMAX], int a,int k){

   for(int i = 0;i<k;i++){

	q3[i] = a * q[i];

 }

}

