#include <math.h>
#include <stdio.h>
#define NMAX 4
#define ITERMAX 2

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
 s=0;
 for(i=0;i<n;i++)
 s+=x[i]*y[i];
 return(s);
}

double ea[NMAX];double eb[NMAX]; int n;

double eq_f(double x)
{
 int k;
 double p0,p1,p2;
 p0=0.0;p1=0.1;
 for(k=0;k<=n-1;k++)
 {
  p2=(ea[k+1]-x)*p1-eb[k]*p0;
  p0=p1;p1=p2;
 }
 return(p2);
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


void transpose(double tab1[NMAX][NMAX],int n,double tab2[NMAX][NMAX]){

for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
		tab2[j][i]=tab1[i][j];

 }


int vp_lanczos(double A[NMAX][NMAX],double q[NMAX],
double T[NMAX][NMAX],double Q[NMAX][NMAX],int k)
{
 int i,j;
 double b,q1[k],q2[k], q3[k],q_t[k],a,w;


// on commence par stocker e vecteur dans la matrice 
 for(i=0;i<k;i++)
 {
  Q[i][0]=q[i];
  }
 
//calcul store alpha dans la grande matrice
 prod_mat_vect(A,k,q,q1); // A * q = q1
 a = ps(q1,q,k);
 T[0][0] = a;

// calcul de beta et store
vect_sca(q, q2, -a, k);    // q3 = a * q;
somme_vect(q1,q2, q3, k);
b = norme_vect(q3, k);

T[1][0] = b;
T[0][1] = b;

//calcul et stockage de q2
  for(i=0;i<k;i++)
 {
  Q[i][1]= q3[i] / b;
  }

// pour quelque soit i
 for(j=1;j<k-1;j++)
 {
 

//calcul store alpha dans la grande matrice
 prod_mat_vect(A,k,Q[j],q1); // A * q = q1
 vect_sca(Q[j-1], q_t, -b,k);
 somme_vect( q1, q_t, q1,k);
 a = ps(q1,Q[j],k);
 T[j][j] = a;

// calcul de beta et store
 vect_sca(Q[j], q2, -a, k);    // q3 = a * q;
 somme_vect(q1,q2, q3, k);
 b = norme_vect(q3, k);
 T[j][j+1] = b;
 T[j+1][j] = b;

 if(b==0) return(0);
 //calcul et stockage de q+


  for(i=0;i<k;i++)
 {
  Q[i][j+1]= q3[i] / b;
  }

    if(j == k-2){
 prod_mat_vect(A,k,Q[j],q1); // A * q = q1
 vect_sca(Q[j-1], q_t, -b,k);
 somme_vect( q1, q_t, q1,k);
 a = ps(q1,Q[j],k);
 T[j+1][j+1] = a;
 }

 }

 return 1;
}

//***************************
void arnoldi(double A[NMAX][NMAX],double q[NMAX],
double T[NMAX][NMAX],double Q[NMAX][NMAX],int k)
{
 int i,j, ep=0;
 double b,q1[k],q2[k], q3[k],q_t[k],a,w;

 for(i=0;i<k;i++)
	q_t[i]=0;

// on commence par stocker q vecteur dans la matrice 
 for(i=0;i<n;i++)
 {
  Q[i][0]=q[i];
  }

// while(ep<ITERMAX){

 for(j=0;j<k;j++)
 {
 

 for(i=0;i<=k;i++){
 prod_mat_vect(A,k,Q[j],q1); // A * q = q1
 T[i][j] = ps(q1,Q[i],k);
     }

 for (i = 0;i <= j-1;i++){
 vect_sca(Q[j], q2, T[i][j], k);    
 somme_vect(q2,q_t, q2, k);

 for(i=0;i<k;i++)
	q_t[i]=q2[i];

    }

 somme_vect(q1,q2, q3, k);
 T[j+1][j] = norme_vect(q3, k);

  for(i=0;i<k;i++)
 {
  Q[i][j+1]= q3[i] / T[j+1][j];
  }
 }


  ep++;
 //}

}
