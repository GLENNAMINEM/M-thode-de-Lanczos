#include <math.h>
#include "soinclude.h"
#include <stdlib.h>
#include <stdio.h>
 
void main()
{

 int i,j;
 double a[NMAX][NMAX],mu[NMAX][NMAX];
 double q[NMAX],v[NMAX];
 double T[NMAX][NMAX],tmp_T[NMAX][NMAX],Q[NMAX][NMAX];
 printf("Méthode de Lanczos\n");
 n=4; //la taille de la matrice
for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++){
	if(i==j){
	a[i][j] = 5*i+5;
	}
	a[i][j] = j * 10 +3*i+10 ; 
        a[j][i]=a[i][j];

  }
	q[i] = i*5 +10;
}

 printf("Matrice A\n");
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  printf("%17.9e",a[i][j]);
  printf("\n");
 }

 int p = norme_vect(q,n);
  for(j=0;j<n;j++)
	q[j] = q[j] /p;


for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		T[i][j]=0;




 int value = 0;


 value = vp_lanczos(a , q, T, Q, n);

 //arnoldi(a , q, T, Q, n);

 //transpose( Q, n, tmp_T);

 //prod_mat( Q, tmp_T, Q,n);
 //double c = ps( Q[1], Q[0],n);

 printf("Matrice Q\n");
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  printf("%17.9e",Q[i][j]);
  printf(" \n");
 }


 printf("Matrice T\n");
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  printf("%17.9e",T[i][j]);
  printf("\n");
 }

}
