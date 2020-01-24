#include <math.h>
#include "p_soinclude.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
 
void main(int argc, char** argv)
{

 int i,j,rang,nbP;
 double d,d1,d2[2];
 double Q[NMAX][NMAX];//store qi
 double q[NMAX];
 double T[NMAX][NMAX],tmp_T[NMAX][NMAX];//store hi
 double a0;
 int n=NMAX;
 int k = n;
 double b,q1[NMAX],q2[NMAX], q3[NMAX],q_t[NMAX];
 double a[NMAX][NMAX];



    MPI_Init(&argc, &argv);//***************init
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WORLD, &nbP);

//init array 0
 for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		T[i][j]=0;

//symetrique array
for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++){
	if(i==j){
	a[i][j] = 5*i+5;
	}
	a[i][j] = j * 10 +3*i + 10; 
        a[j][i]=a[i][j];

  }
	q[i] = i*5 +10;
}



if(rang == 0){ //juste rang 0

 printf("Méthode de Lanczos\n");

 printf("Matrice A\n");
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  printf("%17.9e",a[i][j]);
  printf("\n");
 }
}

//normalisée par tous les procs
 int p = norme_vect(q,n);
  for(j=0;j<n;j++)
	q[j] = q[j] /p;
 
// on commence par stocker q dans Q
 for(i=0;i<n;i++)
 {
  Q[i][0]=q[i];
  }

//calcul store alpha dans la grande matrice
//on remplace mat_vec par prodruit scalaire each proc

  d = ps(a[rang+rang],q,k);//each proc takes twno rows in array
  d1 = ps(a[rang+1+rang],q,k);
  d2[0]=d;d2[1]=d1;

//all proc sent it array to rang=0
  MPI_Gather(&d2, 2, MPI_DOUBLE, q1, 2 ,MPI_DOUBLE,0, MPI_COMM_WORLD);
 


  if(rang == 0){

   a0 = ps(q1,q,k); 
   T[0][0] = a0;  
	}


  if(rang == 0){
  vect_sca(q, q2, -a0, k);    
   }

  if(rang == 0){
  somme_vect(q1,q2, q3, k); 
  b = norme_vect(q3, k); 
  T[1][0] = b;
  T[0][1] = b;
  }

//rang=0 send b and q3 at other proc
MPI_Bcast(q3,k,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(&b,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

//calcul et stockage de q2
  for(i=0;i<k;i++)
 {
  Q[i][1]= q3[i] / b;
  }

 for(j=1;j<k-1;j++)
 {
 
  //calcul store alpha dans la grande matrice
  d = ps(a[rang+rang],Q[j],k);
  d1 = ps(a[rang+1+rang],Q[j],k);
  d2[0]=d;d2[1]=d1;

  MPI_Gather(&d2, 2, MPI_DOUBLE, q1, 2 ,MPI_DOUBLE,0, MPI_COMM_WORLD);

  if ( rang == 0){

      vect_sca(Q[j-1], q_t, -b,k);
      somme_vect( q1, q_t, q1,k);
      a0 = ps(q1,Q[j],k);
      T[j][j] = a0;
  }

   // calcul de beta et store
   if ( rang == 0){
       vect_sca(Q[j], q2, -a0, k);    
       somme_vect(q1,q2, q3, k);
       b = norme_vect(q3, k);
       T[j+1][j] = b;
       T[j][j+1] = b;
        }
  MPI_Bcast(q3,k,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&b,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if(b==0) break;

  //calcul et stockage de q+
  for(i=0;i<k;i++)
 {
  Q[i][j+1]= q3[i] / b;
  }

  if(j == k-2 ){

  d = ps(a[rang+rang],Q[j],k);
  d1 = ps(a[rang+1+rang],Q[j],k);
  d2[0]=d;d2[1]=d1;

  MPI_Gather(&d2, 2, MPI_DOUBLE, q1, 2 ,MPI_DOUBLE,0, MPI_COMM_WORLD);

  if(rang == 0){

      vect_sca(Q[j-1], q_t, -b,k);
      somme_vect( q1, q_t, q1,k);
      a0 = ps(q1,Q[j],k);
      T[j+1][j+1] = a0;
    }

   }

}


//rank=0 print all array
if(rang == 0){ //juste rang 0

 printf("Matrice Q\n");
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  printf("%17.9e",Q[i][j]);
  printf("\n");
 }

 printf("Matrice T\n");
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  printf("%17.9e",T[i][j]);
  printf("\n");
 }

}
 
     MPI_Finalize();

}
