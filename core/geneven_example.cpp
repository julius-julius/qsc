/* libraries to be included */  
#include <stdio.h>  
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <cln/cln.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cln/number.h>
#include <cln/integer.h>
#include <cln/real.h>
#include <cln/float.h>
#include <cln/io.h>
#include <cln/float_io.h>
#include <cln/malloc.h>

using namespace std;
using namespace cln;




/* Definitions of functions */



/**********************************************************************/



/* matrix inverse of a, e should be the unit martix at calling and after evaluation it stores the inverse, as well, n=dimension, 
 *value of a changes during evaluaton due to pivoting */
void inverzE(cl_N * a,cl_N * e, long int n)
{  
  
  long int d;
  
  long int i,j,k,i0,i1,xm;
  
  cl_N am,c;
  
   d = 1;
   for ( k = 0; k < n-1; k++ )
   {
  am=a[k*n+k]; xm=k;
  for ( i0 = k; i0 <= n-1; i0++ )
  {
    if (abs(am) < abs(a[i0*n+k]))
    {
      am=a[i0*n+k]; xm=i0;

  for ( i1 = 0; i1 <= n-1; i1++ )
  {
    c=a[xm*n+i1]; a[xm*n+i1]=a[k*n+i1]; a[k*n+i1]=c;
    c=e[xm*n+i1]; e[xm*n+i1]=e[k*n+i1]; e[k*n+i1]=c;
  }     
      
      d=d*(-1);
    }//if-ends
  }//i0
 
     
      for ( i = k+1; i <= n-1; i++ )
      {
	 a[i*n+k]=a[i*n+k]/a[k*n+k];
	 for ( j = k+1; j <= n-1; j++ )
	 {
	   a[i*n+j]=a[i*n+j]-a[i*n+k]*a[k*n+j];
	 }//j
      }//i
 }//k

 cl_N x[n];
 cl_N y[n];
 cl_N b[n];
  for ( j = 0; j <= n-1; j++ )
  {  
    for ( i = 0; i <= n-1; i++ )
    {
      b[i]=e[i*n+j];
    }//i
   
y[0]=b[0];
  for ( k = 1; k < n; k++)
  {
    y[k]=b[k];
    for ( i = 0; i <= k-1; i++)
    {
      y[k]=y[k]-a[k*n+i]*y[i];
    }//i  
 }//k

    x[n-1]=y[n-1]/a[(n-1)*n+n-1];
    for ( k = n-2; k >= 0; k--)
    {
      x[k]=complex(0,0);
      for ( i = k+1; i < n; i++ )
      {
	x[k]=x[k]-a[k*n+i]*x[i];
      }//i
      x[k]=(y[k]+x[k])/a[k*n+k];
    }//k
    
    for ( i = 0; i <= n-1; i++ )
    {
      e[i*n+j]=x[i];
    }//i
  }//j
 
}


/*************************************************************************************************/


/* Function to compute the inverse transpose of the nxn matrix a, and the result is written into e  */

void inversetransposeE(cl_N * a,cl_N * e, long int n){
  /* a-t transzponaljuk es e-t egyseggel toltjuk fel */
long int i,j,d;
d=1;
  cl_N * acp;
  acp=new cl_N[n*n];

for(i=0;i<=n-1;i++){
  for(j=0;j<=n-1;j++){
    acp[i*n+j]=a[j*n+i];
    if(i==j){e[i*n+j]=1;} else{e[i*n+j]=0;}
  }
}

  
  long int k,i0,i1,xm;
  
  cl_N am,c;

   d = 1;
   for ( k = 0; k < n-1; k++ )
   {

  am=acp[k*n+k]; xm=k;
  for ( i0 = k; i0 <= n-1; i0++ )
  {
    if (abs(am) < abs(acp[i0*n+k]))
    {
      am=acp[i0*n+k]; xm=i0;

  for ( i1 = 0; i1 <= n-1; i1++ )
  {
    c=acp[xm*n+i1]; acp[xm*n+i1]=acp[k*n+i1]; acp[k*n+i1]=c;
    c=e[xm*n+i1]; e[xm*n+i1]=e[k*n+i1]; e[k*n+i1]=c;
  }     
  
      d=d*(-1);
    }//if-ends
  }//i0

      for ( i = k+1; i <= n-1; i++ )
      {
	 acp[i*n+k]=acp[i*n+k]/acp[k*n+k];
	 for ( j = k+1; j <= n-1; j++ )
	 {
	   acp[i*n+j]=acp[i*n+j]-acp[i*n+k]*acp[k*n+j];
	 }//j
      }//i
 }//k

 cl_N x[n];
 cl_N y[n];
 cl_N b[n];
  for ( j = 0; j <= n-1; j++ )
  {  
    for ( i = 0; i <= n-1; i++ )
    {
      b[i]=e[i*n+j];
    }//i
    
y[0]=b[0];
  for ( k = 1; k < n; k++)
  {
    y[k]=b[k];
    for ( i = 0; i <= k-1; i++)
    {
      y[k]=y[k]-acp[k*n+i]*y[i];
    }//i  
 }//k
 
    x[n-1]=y[n-1]/acp[(n-1)*n+n-1];
    for ( k = n-2; k >= 0; k--)
    {
      x[k]=complex(0,0);
      for ( i = k+1; i < n; i++ )
      {
	x[k]=x[k]-acp[k*n+i]*x[i];
      }//i
      x[k]=(y[k]+x[k])/acp[k*n+k];
    }//k
   
    for ( i = 0; i <= n-1; i++ )
    {
      e[i*n+j]=x[i];
    }//i
  }//j
}


/* Function to solve the linear problem: a.x=b, soution is written into x, n=dimension of the problem */


void linsolvepE(cl_N * a,cl_N * b,cl_N x[],long int n)//decomposition with pivoting, then solving the equation
{
	long int d=1;
  	for(long int k=0;k<n-1;k++)
  	{
   		cl_N am = a[k*n+k]; //pivoting
	 	long int xm=k;
	 	for(long int i=k;i<n;i++)
	 	{
	   		if(abs(am)<abs(a[i*n+k]))
	  		{
	     		am=a[i*n+k];
	     		xm=i;
	   			cl_N c;
				for (long int i=0;i<n;i++)//switching lines
				{
					c=a[xm*n+i];
					a[xm*n+i]=a[k*n+i];
					a[k*n+i]=c;
	 			}
	 			c=b[xm];
	 			b[xm]=b[k];
	 			b[k]=c;
	   			d=-d;
	   		}
	 	}
     	for(long int i=k+1;i<n;i++)
     	{
	 		a[i*n+k]/=a[k*n+k];
	 		for(long int j=k+1;j<n;j++)
	 		{
	  			a[i*n+j]-=a[i*n+k]*a[k*n+j];
	 		}
     	}
  	}
	//inverting matrix
	cl_N y[n]; 
  	y[0]=b[0];
 	for(long int k=1;k<n;k++) //solving Ly=b
 	{
   		y[k]=b[k];
   		for(long int j=0;j<=k-1;j++)
   		{
     		y[k]-=a[k*n+j]*y[j];
   		}
   	}
   	x[n-1]=y[n-1]/a[(n-1)*n+n-1]; //solving Ux=y
   	for(long int k=n-2;k>=0;k--)
 	{
     	x[k]=complex(0,0);
     	for(long int j=k+1;j<n;j++)
     	{
			x[k]-=a[k*n+j]*x[j];
     	}
     	x[k]=(y[k]+x[k])/a[k*n+k];
   }
}



/* Real Binomial[z,n], with WorkingPrecision=prec and n=integer */
cl_R binomial(cl_R z,long int n,float_format_t prec){
cl_R z1;
cl_R zm=-z;
long int i;
long int n2=abs(n);
cl_R n1=cl_float((cl_I) abs(n),prec);
cl_R elojel;
if(n%2==0) {elojel=cl_float(1,prec);} else {elojel=cl_float(-1,prec);}
if(n==0){z1=cl_float(1,prec);}else{
if(z>0){z1=z/n1;
for(i=1;i<=n2-1;i++ ){
      z1=z1*((z-i)/( n1-i));
    }}else{ z1=(zm+n1-1)/n1;
for(i=1;i<=n2-1;i++ ){
      z1=z1*(zm+n1-1-i)/(n1- i);
    } z1=elojel*z1;
    }}
    return z1;
}


/* Complex Binomial[z,n], with WorkingPrecision=prec and n=integer */

cl_N cbinomial(cl_N z,long int n,float_format_t prec){
cl_N z1;
cl_N zm;
cl_N egy=complex(cl_float(1,prec),cl_float(0,prec));
cl_N megy=complex(cl_float(-1,prec),cl_float(0,prec));
zm=(-1)*z;
long int i;
long int n2=abs(n);
cl_N n1= abs(n);
cl_N elojel;
if(n%2==0) {elojel=egy;} else {elojel=megy;}
if(n==0){z1=egy;}else{
if(realpart(z)>0){z1=z/n1;
for(i=1;i<=n2-1;i++ ){
      z1=z1*((z-i)/(n1-i));
    }}else{ z1=(zm+n1-1)/n1;
for(i=1;i<=n2-1;i++ ){
      z1=z1*(zm+n1-1-i)/(n1-i);
    } z1=elojel*z1;
    }}
    return z1;
}




/* x[u] with short cut */
cl_N XS(cl_N u){
cl_N r1;
cl_N u2=u*u;
r1=u*(1+sqrt(1-4/u2))/2;
return r1;  }

/* 1/x[u] with short cut */
cl_N invXS(cl_N u){
cl_N r1;
cl_N u2=u*u;
r1=u*(1-sqrt(1-4/u2))/2;
return r1;  }




/* x[u] with long cut */
cl_N X(cl_N u){
cl_N r1;
cl_N u2=u*u;
cl_N Ifel=complex(0,1)/2;
r1=u/2-Ifel*sqrt(4-u2);
return r1;  }

/* 1/x[u] with long cut */
cl_N invX(cl_N u){
cl_N r1;
cl_N u2=u*u;
cl_N Ifel=complex(0,1/2);
r1=u/2+Ifel*sqrt(4-u2);
return r1;  }

/*  function: (-1)^m */
cl_N m1m(long int m){
cl_N r;
if(m % 2 == 0){ r=complex(1,0);}
else{ r=complex(-1,0);}
return r;
}

/*  function: (-1)^m with precison=prec*/
cl_N m1m_f(long int m,float_format_t prec){
cl_N r;
if(m % 2 == 0){ r=complex(cl_float(1,prec),cl_float(0,prec));}
else{ r=complex(cl_float(-1,prec),cl_float(0,prec));;}
return r;
}



/* An auxiliary function for Chebyshev-expansion */

cl_N Cffunc(long int k, long int i,long int lc,float_format_t prec){
  cl_N k1,i1,pi1,lcs1,r;
  k1=complex(k,0);
  i1=complex(i,0);
  pi1=pi(prec);
  lcs1=complex(lc,0);
  r=cos(pi1*(2*k1+1)*i1/(2*lcs1));
  return r;
};


/* An auxiliary function for the 1/u expansion of integer powers of x[u], n,s are integers */
cl_N kappa(long int n, long int s,float_format_t prec)
{cl_N r,n1,s1;
n1=complex(n,0);
s1=complex(s,0);;
if(s==0){r=complex(1,0); }
else{
     if(n==0){ r=complex(0,0); }
  else{ r=n1/s1*cbinomial(n1+2*s1-1,s-1,prec); } 
  
}
return r;
}

/* New auxiliary function for half integer cases: kappa(n,s)=kappabar(2*n,s). Useful in case of n=half-integer  */
cl_N kappabar(long int n, long int s,float_format_t prec)
{cl_N r,n1,s1,two;
two=cl_float(2,prec);
n1=complex(n,0);
s1=complex(s,0);
if(s==0){r=complex(1,0);}
else{
     if(n==0){r=complex(0,0);}
  else{ 
        r=n1/s1/two*cbinomial(n1/two+2*s1-1,s-1,prec);   
  }
}
return r;
}

/* Auxiliary function for the 1/u expansion of P_a[u] */

/*  twiceMta=2*Mt[a]=integer */
cl_N fsigma(long int twiceMta,long int n,long int r,cl_N g,float_format_t prec){
  long int s;
  cl_N R=complex(0,0);
  long int k,q0;
  ldiv_t kq=ldiv(n,2);
  k=kq.quot;
  q0=kq.rem;
  for(s=0;s<=k-r;s++){ R+=kappabar(twiceMta,s,prec)*kappa(2*r+q0,k-r-s,prec); }
  
  R*=expt(sqrt(g),twiceMta+2*n);
  
  return R;
  
}

/* Auxiliary function for the 1/u expansion of P^a[u] */

/*  twiceMta=2*Mt[a]=integer */
cl_N fsigmaf(long int twiceMta,long int n,long int r,cl_N g,float_format_t prec){
  long int s;
  cl_N R=complex(0,0);
  long int k,q0;
  ldiv_t kq=ldiv(n,2);
  k=kq.quot;
  q0=kq.rem;
  for(s=0;s<=k-r;s++){ R+=kappabar(2-twiceMta,s,prec)*kappa(2*r+q0,k-r-s,prec); }
  
  R*=expt(sqrt(g),2-twiceMta+2*n);
  
  return R;
  
}

/* Auxiliary function for the 1/u expansion of P_a[u] */

void sigmasubfunc2(cl_N * sigma/*(1+NQ)*(1+N0)*/,long int twiceMta,long int N00/*N0*/,long int NQ0/*NQ*/,cl_N g,float_format_t prec){
  long int NQp1=NQ0+1;
  long int n,r;
  
  for(n=0;n<=NQ0;n++){  
    for(r=0;r<=N00;r++){  sigma[r*NQp1+n]=complex(0,0);
      sigma[r*NQp1+n]+=fsigma(twiceMta,2*n,r,g,prec);
      
    }/*r*/
  }/*n*/
  
}

/* Auxiliary function for the 1/u expansion of P^a[u] */
void sigmasupfunc2(cl_N * sigma/*(1+NQ)*(1+N0)*/,long int twiceMta,long int N00/*N0*/,long int NQ0/*NQ*/,cl_N g,float_format_t prec){
  long int NQp1=NQ0+1;
  long int n,r;
  
  for(n=0;n<=NQ0;n++){  
    for(r=0;r<=N00;r++){  sigma[r*NQp1+n]=complex(0,0);
      sigma[r*NQp1+n]+=fsigmaf(twiceMta,2*n,r,g,prec);
      
    }/*r*/
  }/*n*/
  
}


/* An auxiliary function for the 1/u expansion of P_a[u]  */

void kanfunc2(cl_N * ksub[4]/*NQ+1*/,cl_N * c[4]/*N0+1*/,cl_N * sigmasub[4]/*(1+NQ)*(1+N0)*/,long int NQ0/*NQ*/,long int N00/*N0*/){
  long int n,k,q0,r,a;
  long int NQp1=NQ0+1;
  for(a=0;a<=3;a++){
  for(n=0;n<=NQ0;n++){   
    ksub[a][n]=complex(0,0);
for(r=0;r<=min(n,N00);r++){  ksub[a][n]+=c[a][r]*sigmasub[a][r*NQp1+n];  }/*r*/    
  }/*n*/
  }/*a*/
  }
  

/**************************************************************************************************************/


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void T41func(cl_N * T41/*lmax*(2*lmax-2)*/,cl_N * m4k/*lmax+2*/,long int lmax,float_format_t prec)
{ long int k,m,ref;
  for(m=0;m<=2*lmax-3;m++){ ref=ldiv(m,2L).quot;
    for(k=0;k<=lmax-1;k++){
     if(k<=ref){  T41[m*lmax+k]=binomial(cl_float(-2.0*(k+1) ,prec),m-2*k,prec)*m4k[k+1]; }else{ T41[m*lmax+k]=complex(0,0); }
    }
  }
}

/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */


void T2func(cl_N * T2/*(lmax+1)*(lmax-1) */,cl_N * m1p4k/*lmax+2*/,long int lmax,float_format_t prec){
  long int l,k,lmm1;
  lmm1=lmax-1;
    for(l=2;l<=lmax;l++){
  for(k=0;k<=lmax-2;k++){ 
    if(k<=l-2){ T2[l*lmm1+k]=cbinomial(cl_float((cl_I) -2*(k+1),prec),2*(l-k-1),prec)*m1p4k[l-k-1];  }else{ T2[l*lmm1+k]=complex(0,0); }
  }  
  }
  /* l=0,1 */
  for(l=0;l<=1;l++){
  for(k=0;k<=lmax-2;k++){ T2[l*lmm1+k]=complex(0,0); }
    }
}


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void T1func(cl_N * T1/*(lmax+1)*(lmax-1)*/,cl_N * m1p4k/*lmax+2*/,long int lmax,float_format_t prec){
  long int l,k,lmm1;
  lmm1=lmax-1;
    for(l=2;l<=lmax;l++){
  for(k=0;k<=lmax-2;k++){ 
    if(k<=l-2){ T1[l*lmm1+k]=cbinomial(cl_float((cl_I) -2*(k+1),prec),2*(l-k)-1,prec)*m1p4k[l-k-1];  }else{ T1[l*lmm1+k]=complex(0,0); }
  }  
  }
  /* l=0,1 */
  for(l=0;l<=1;l++){
  for(k=0;k<=lmax-2;k++){ T1[l*lmm1+k]=complex(0,0); }}
}

/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void S1func2(cl_N * S1/*(NQ-1)*(NQ+1)*/,cl_N * m1p4k/*NQ+2*/,long int NQ0,float_format_t prec){
  long int n,j,NQm1;
  NQm1=NQ0-1;
  for(n=0;n<=NQ0;n++){
   for(j=0;j<=NQ0-2;j++) {
     if(j<=n-2){ S1[n*NQm1+j]=cbinomial(cl_float((cl_I) -2*(j+1),prec),2*(n-j-1),prec)*m1p4k[n-j-1]; }else{ S1[n*NQm1+j]=complex(0,0); }
  }/*j*/
  }/*n*/
}

/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void S31func2(cl_N * S31/*(2*NQ-1)*(NQ-1)*/,cl_N * m4k/*NQ+2*/,long int NQ0,float_format_t prec){
  long int k,j,NQm1,ref;
  NQm1=NQ0-1;
  for(k=0;k<=2*NQ0-2;k++){ ref=(ldiv(k+1,2).quot)-1;
    for(j=0;j<=NQ0-2;j++){
      if(j<=ref){ S31[k*NQm1+j]=cbinomial(cl_float((cl_I) -2*(j+1),prec),k-2*j-1,prec)*m4k[j+1];  }else{ S31[k*NQm1+j]=complex(0,0); }
    }/*j*/
  }/*k*/
  
}


/* Auxiliary function to impose P_a P^a=0 at the level of the coefficients: c_{a,n} and c^{a}_n */

cl_N cfixalo(long int a1,long int a2,long int n,cl_N * c[4]/*N0+1*/,cl_N * cf[4]/*N0+1*/){
  long int m,a;
  cl_N r;
  r=complex(0,0);
  
  for(a=a1;a<=a2;a++){ r-=cf[a][n]*c[a][0]; }
  
 for(a=0;a<=3;a++){ for(m=0;m<=n-1;m++){ r-=c[a][n-m]*cf[a][m]; }/*m*/ }/*a*/
  
  return r;
  
}


/* A function giving the value of the A_{a0}A^{a0}=... constraint as a 4-dimensional vector */

void AAfunc(cl_N mt[4],cl_N mh[4],cl_N AAf[4]){
  
  cl_N II=complex(0,1);
  long int a,j;
  cl_N r;
  
  for(a=0;a<=3;a++){  r=II;
  
  for(j=0;j<=3;j++){ r*=(mt[a]-mh[j]);  }
  
  for(j=0;j<=3;j++){ if(j==a){ ; }else{ r/=(mt[a]-mt[j]); }  }
  
  AAf[a]=r;
    
    }/*a*/
  
}


/*Function computing the "unconstrained" variables of the problem out of  {Delta,c_{a,n},c^a_n} */

void CtoV2(cl_N Delta,cl_N * c[4],cl_N * cf[4],long int Mtint[4],cl_N * V/*dimV*/,long int Nch[4],long int N00){
  long int n,j,k,N1;
  cl_N * ctilde;
  
  N1=N00+Nch[1]+Nch[2]+Nch[3];
  
  V[0]=Delta;
  
  for(j=1;j<=N00;j++){ V[j]=c[0][j];  }/*j*/
  
  ctilde=new cl_N[N00];
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(2*n==(Mtint[0]-Mtint[1]-2) ){ k++; }else{ ctilde[n-k]=c[1][n+1]; }
  }/*n*/
  
  for(j=N00+1;j<=N00+Nch[1];j++){ V[j]=ctilde[j-N00-1];  }/*j*/
  
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(2*n==(Mtint[0]-Mtint[2]-2)  ||  2*n==(Mtint[1]-Mtint[2]-2)  ){ k++; }else{ ctilde[n-k]=c[2][n+1]; }
  }/*n*/
  
  for(j=N00+Nch[1]+1;j<=N00+Nch[1]+Nch[2];j++){ V[j]=ctilde[j-N00-Nch[1]-1];  }/*j*/  
    
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(2*n==(Mtint[0]-Mtint[3]-2) ||  2*n==(Mtint[1]-Mtint[3]-2) ||  2*n==(Mtint[2]-Mtint[3]-2)   ){ k++; }else{ ctilde[n-k]=c[3][n+1]; }
  }/*n*/
  
  for(j=N00+Nch[1]+Nch[2]+1;j<=N1;j++){ V[j]=ctilde[j-N00-Nch[1]-Nch[2]-1];  }/*j*/    
  
  for(j=N1+1;j<=N1+N00;j++){ V[j]=cf[3][j-N1];  }/*j*/
    
    for(j=N1+N00+1;j<=N1+2*N00;j++){ V[j]=cf[2][j-N1-N00];  }/*j*/
      
      for(j=N1+2*N00+1;j<=N1+3*N00;j++){ V[j]=cf[1][j-N1-2*N00];  }/*j*/
    
  
    
  delete [] ctilde;     
}


/********************************************************************/

/********************************************************************/

/* The inverse of the previous function, namely it  computes {Delta,c_{a,n},c^a_n} out of the the "unconstrained" variables of the problem  */

void VtoC2(cl_N * V/*7*N0-5*/,cl_N * ptrDelta/*&Delta*/,cl_N * c[4],cl_N * cf[4],long int Mtint[4],long int Nch[4],long int N00){
  long int n,j,k,N1;
  cl_N * ctilde;
  
  N1=N00+Nch[1]+Nch[2]+Nch[3];
  
  *ptrDelta=V[0];
  
   for(j=1;j<=N00;j++){ c[0][j]=V[j];  }/*j*/
    
  ctilde=new cl_N[N00];
 
  for(j=N00+1;j<=N00+Nch[1];j++){ ctilde[j-N00-1]=V[j];  }/*j*/
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(2*n==(Mtint[0]-Mtint[1]-2) ){ c[1][n+1]=complex(0,0);  k++; }else{ c[1][n+1]=ctilde[n-k]; }
  }/*n*/

  for(j=N00+Nch[1]+1;j<=N00+Nch[1]+Nch[2];j++){ ctilde[j-N00-Nch[1]-1]=V[j];  }/*j*/ 
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(2*n==(Mtint[0]-Mtint[2]-2)  ||  2*n==(Mtint[1]-Mtint[2]-2)  ){ c[2][n+1]=complex(0,0); k++; }else{ c[2][n+1]=ctilde[n-k]; }
  }/*n*/

  for(j=N00+Nch[1]+Nch[2]+1;j<=N1;j++){ ctilde[j-N00-Nch[1]-Nch[2]-1]=V[j];  }/*j*/ 
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(2*n==(Mtint[0]-Mtint[3]-2) ||  2*n==(Mtint[1]-Mtint[3]-2) ||  2*n==(Mtint[2]-Mtint[3]-2)   ){ c[3][n+1]=complex(0,0); k++; }
    else{ c[3][n+1]=ctilde[n-k]; }
  }/*n*/
  
  for(j=N1+1;j<=N1+N00;j++){ cf[3][j-N1]=V[j];  }/*j*/
    
    for(j=N1+N00+1;j<=N1+2*N00;j++){ cf[2][j-N1-N00]=V[j];  }/*j*/
      
      for(j=N1+2*N00+1;j<=N1+3*N00;j++){ cf[1][j-N1-2*N00]=V[j];  }/*j*/
    
 delete [] ctilde;   
    
} 




/*A function to compute A_a in the convention of Volin and Marboe: arxiv.org>hep.th: 1812.09238  */

cl_N VolinAfunc(long int a,cl_N mt[4],cl_N mh[4]){
  
  cl_N II=complex(0,1);
  long int j;
  cl_N r;
  
    r=(mh[0]-mt[a])*(mt[a]-mh[1]);
  
  for(j=a+1;j<=3;j++){ r/=(II*(mt[a]-mt[j]));  }
  
  return r;
     
} 

/*********************************************************************/

/* Function to select the proper equations from the whole set */

void UJCtoE2(long int Nch[4],cl_N * deltac,cl_N * deltacf,long int Mtint[4],cl_N * E,long int N00){
  long int n,j,k;
  cl_N * ctilde;
  
  E[0]=imagpart(deltac[0*4+0]); 
  
  for(j=1;j<=N00;j++){ E[j]=imagpart(deltac[j*4+0]);  }/*j*/
  
  ctilde=new cl_N[N00];
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if( (Mtint[0]-Mtint[1])%2==0 && n==((Mtint[0]-Mtint[1])/2-1)){ k++; }else{ ctilde[n-k]=realpart(deltac[(n+1)*4+1]);/*c[1][n+1];*/ }
  }/*n*/
  
  for(j=N00+1;j<=N00+1+Nch[1]-1;j++){ E[j]=ctilde[j-N00-1];  }/*j*/
  
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(((Mtint[0]-Mtint[2])%2==0  && n==((Mtint[0]-Mtint[2])/2-1) ) || ( ( Mtint[1]-Mtint[2])%2==0  && n==((Mtint[1]-Mtint[2])/2-1) ) ){ k++; }else{ ctilde[n-k]=imagpart(deltac[(n+1)*4+2]);/*c[2][n+1];*/ }
  }/*n*/
  
  for(j=N00+1+Nch[1];j<=N00+1+Nch[1]+Nch[2]-1;j++){ E[j]=ctilde[j-(N00+1+Nch[1])];  }/*j*/  
    
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(((Mtint[0]-Mtint[3])%2==0 &&  n==((Mtint[0]-Mtint[3])/2-1) ) || ((Mtint[1]-Mtint[3])%2==0 &&  n==((Mtint[1]-Mtint[3])/2-1) ) || ((Mtint[2]-Mtint[3])%2==0  && n==((Mtint[2]-Mtint[3])/2-1) )    ){ k++; }else{ 
      ctilde[n-k]=realpart(deltac[(n+1)*4+3]);  }
  }/*n*/
  
  for(j=N00+1+Nch[1]+Nch[2];j<=N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ E[j]=ctilde[j-(N00+1+Nch[1]+Nch[2])];  }/*j*/    

  
  for(j=N00+1+Nch[1]+Nch[2]+Nch[3];j<=2*N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ E[j]=imagpart(deltacf[(j+1-(N00+1+Nch[1]+Nch[2]+Nch[3]))*4+3]);  }/*j*/
    
  for(j=2*N00+1+Nch[1]+Nch[2]+Nch[3];j<=3*N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ E[j]=realpart(deltacf[(j+1-(2*N00+1+Nch[1]+Nch[2]+Nch[3]))*4+2]);  }/*j*/  
    
  for(j=3*N00+1+Nch[1]+Nch[2]+Nch[3];j<=4*N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ E[j]=imagpart(deltacf[(j+1-(3*N00+1+Nch[1]+Nch[2]+Nch[3]))*4+1]);  }/*j*/  
    
  delete [] ctilde;     
}

/* Function to select  only the proper equations out of the whole set for the derivative */

void UJCtoDE2(long int Nch[4],cl_N * deltac,cl_N * deltacf,long int Mtint[4],cl_N * DE,long int dimV,long int N00,long int n0){
  long int n,j,k;
  cl_N * ctilde;
  
  DE[0*dimV+n0]=imagpart(deltac[0*4+0]); 
  
  for(j=1;j<=N00;j++){ DE[j*dimV+n0]=imagpart(deltac[j*4+0]);  }/*j*/
  
  ctilde=new cl_N[N00];
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if( (Mtint[0]-Mtint[1])%2==0 && n==((Mtint[0]-Mtint[1])/2-1)){ k++; }else{ ctilde[n-k]=realpart(deltac[(n+1)*4+1]);/*c[1][n+1];*/ }
  }/*n*/
  
  for(j=N00+1;j<=N00+1+Nch[1]-1;j++){ DE[j*dimV+n0]=ctilde[j-N00-1];  }/*j*/
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(((Mtint[0]-Mtint[2])%2==0 &&  n==((Mtint[0]-Mtint[2])/2-1) ) || ( ( Mtint[1]-Mtint[2])%2==0 &&  n==((Mtint[1]-Mtint[2])/2-1) ) ){ k++; }else{ ctilde[n-k]=imagpart(deltac[(n+1)*4+2]); }
  }/*n*/
  
  for(j=N00+1+Nch[1];j<=N00+1+Nch[1]+Nch[2]-1;j++){ DE[j*dimV+n0]=ctilde[j-(N00+1+Nch[1])];  }/*j*/  
    
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(((Mtint[0]-Mtint[3])%2==0 &&  n==((Mtint[0]-Mtint[3])/2-1) ) || ((Mtint[1]-Mtint[3])%2==0 &&  n==((Mtint[1]-Mtint[3])/2-1) ) || ((Mtint[2]-Mtint[3])%2==0  && n==((Mtint[2]-Mtint[3])/2-1) )    ){ k++; }else{ 
      ctilde[n-k]=realpart(deltac[(n+1)*4+3]);  }
  }/*n*/
  
  for(j=N00+1+Nch[1]+Nch[2];j<=N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ DE[j*dimV+n0]=ctilde[j-(N00+1+Nch[1]+Nch[2])];  }/*j*/    
    
  for(j=N00+1+Nch[1]+Nch[2]+Nch[3];j<=2*N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ DE[j*dimV+n0]=imagpart(deltacf[(j+1-(N00+1+Nch[1]+Nch[2]+Nch[3]))*4+3]);  }/*j*/
    
  for(j=2*N00+1+Nch[1]+Nch[2]+Nch[3];j<=3*N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ DE[j*dimV+n0]=realpart(deltacf[(j+1-(2*N00+1+Nch[1]+Nch[2]+Nch[3]))*4+2]);  }/*j*/  
    
  for(j=3*N00+1+Nch[1]+Nch[2]+Nch[3];j<=4*N00+1+Nch[1]+Nch[2]+Nch[3]-1;j++){ DE[j*dimV+n0]=imagpart(deltacf[(j+1-(3*N00+1+Nch[1]+Nch[2]+Nch[3]))*4+1]);  }/*j*/  
    
  delete [] ctilde;     
}


/*A function to determine the set of equations obtained after the Fourier-transformation */ 

void QtoEtypeIInewton2(cl_N * deltac,cl_N * deltacup,cl_N * deltaP,cl_N * deltaPt,cl_N * deltaPup,cl_N * deltaPtup,cl_N * CT,cl_N * CU,long int Nas[4][2],cl_N * suA/*lc*/,long int lc,long int N00,cl_N g,int ncmax[4],int ncfmax[4],int kmax[4],int kfmax[4],int n0min[4],int nf0min[4],int kmin[4],int kfmin[4]){
  
  long int k,a,n,l;
  
  /*symmetric and antisymmetric combinations */
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  cl_N * fSup=new cl_N[4*lc];
  cl_N * fAup=new cl_N[4*lc];
 
  
  cl_N * cS[4];
  cl_N * cA[4];
  cl_N * cSup[4];
  cl_N * cAup[4];
  
  for(a=0;a<=3;a++){ cS[a]=new cl_N[ncmax[a]+1]; cA[a]=new cl_N[ncmax[a]+1]; 
                     cSup[a]=new cl_N[ncfmax[a]+1]; cAup[a]=new cl_N[ncfmax[a]+1];
  }//a
  
 
    for(a=0;a<=0;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaP[k*4+a]+deltaPt[k*4+a])/2; 
                          fA[k*4+a]=(deltaP[k*4+a]-deltaPt[k*4+a])/2/suA[k];
			  
			  fSup[k*4+a]=0;
			  fAup[k*4+a]=0;
    }/*k*/
  }/*a*/
  
  
  
  for(a=1;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaP[k*4+a]+deltaPt[k*4+a])/2; 
                          fA[k*4+a]=(deltaP[k*4+a]-deltaPt[k*4+a])/2/suA[k];
			  
			  fSup[k*4+a]=(deltaPup[k*4+a]+deltaPtup[k*4+a])/2;
			  fAup[k*4+a]=(deltaPup[k*4+a]-deltaPtup[k*4+a])/2/suA[k];
    }/*k*/
  }/*a*/

for(a=0;a<=0;a++){ 
  
 
  for(n=0;n<=ncmax[a];n++){  cS[a][n]=complex(0,0); }//n
    
  for(k=0;k<=kmax[a];k++){  
      for(l=0;l<=lc-1;l++){ cS[a][n0min[a]+2*k]+=fS[l*4+a]*CT[(lc-1-l)*lc+n0min[a]+2*k];  }//l
     }/*k*/
  
   for(n=0;n<=ncmax[a];n++){ cS[a][n]/=lc;   }/*n*/
  
  for(n=0;n<=ncfmax[a];n++){  cSup[a][n]=complex(0,0); }//n
  
}/*a*/


for(a=1;a<=3;a++){ 
  
  for(n=0;n<=ncmax[a];n++){  cS[a][n]=complex(0,0); }//n
    
  for(k=0;k<=kmax[a];k++){  
      for(l=0;l<=lc-1;l++){ cS[a][n0min[a]+2*k]+=fS[l*4+a]*CT[(lc-1-l)*lc+n0min[a]+2*k];  }//l
     }/*k*/
  
 for(n=0;n<=ncmax[a];n++){ cS[a][n]/=lc;   }/*n*/

  for(n=0;n<=ncfmax[a];n++){  cSup[a][n]=complex(0,0); }//n
    
  for(k=0;k<=kfmax[a];k++){  
      for(l=0;l<=lc-1;l++){ cSup[a][nf0min[a]+2*k]+=fSup[l*4+a]*CT[(lc-1-l)*lc+nf0min[a]+2*k];  }//l
     }/*k*/
  
 for(n=0;n<=ncfmax[a];n++){ cSup[a][n]/=lc;   }/*n*/
  
  
}/*a*/

for(a=0;a<=0;a++){ 
  
  cA[a][0]=complex(0,0); cAup[a][0]=complex(0,0);
  
  /*n runs from 1*/ 
  
  for(n=1;n<=ncmax[a];n++){  cA[a][n]=complex(0,0); }//n
    
  for(k=kmin[a];k<=kmax[a];k++){   
    
  for(l=0;l<=lc-1;l++){   cA[a][n0min[a]+2*k]+=fA[l*4+a]*CU[(lc-1-l)*lc+n0min[a]+2*k-1]; }/*l*/
    }//k
  for(n=1;n<=ncmax[a];n++){  cA[a][n]/=lc;  cA[a][n]*=2*complex(0,1)*g; }//n

  for(n=1;n<=ncfmax[a];n++){  cAup[a][n]=complex(0,0); }//n
  
}/*a*/


for(a=1;a<=3;a++){ 
  
  cA[a][0]=complex(0,0); cAup[a][0]=complex(0,0);
  
  
  /*n runs from 1 */ 

  for(n=1;n<=ncmax[a];n++){  cA[a][n]=complex(0,0); }//n
    
  for(k=kmin[a];k<=kmax[a];k++){   
    
  for(l=0;l<=lc-1;l++){   cA[a][n0min[a]+2*k]+=fA[l*4+a]*CU[(lc-1-l)*lc+n0min[a]+2*k-1]; }/*l*/
    }//k
  for(n=1;n<=ncmax[a];n++){  cA[a][n]/=lc;  cA[a][n]*=2*complex(0,1)*g; }//n
  
  for(n=1;n<=ncfmax[a];n++){  cAup[a][n]=complex(0,0); }//n
    
  for(k=kfmin[a];k<=kfmax[a];k++){   
    
  for(l=0;l<=lc-1;l++){   cAup[a][nf0min[a]+2*k]+=fAup[l*4+a]*CU[(lc-1-l)*lc+nf0min[a]+2*k-1]; }/*l*/
    }//k
  for(n=1;n<=ncfmax[a];n++){  cAup[a][n]/=lc;  cAup[a][n]*=2*complex(0,1)*g; }//n
  
}/*a*/

for(a=0;a<=0;a++){ 
  for(n=0;n<=N00;n++){  
    
    if(2*n>=Nas[a][0]){ deltac[n*4+a]=cS[a][(-Nas[a][0]+2*n)]+cA[a][(-Nas[a][0]+2*n)]; }else{ deltac[n*4+a]=cS[a][(Nas[a][0]-2*n)]-cA[a][(Nas[a][0]-2*n)]; } 
    
    deltacup[n*4+a]=0;
    
  }/*n*/
}/*a*/

  
for(a=1;a<=3;a++){ 
  for(n=0;n<=N00;n++){  
    
    if(2*n>=Nas[a][0]){ deltac[n*4+a]=cS[a][(-Nas[a][0]+2*n)]+cA[a][(-Nas[a][0]+2*n)]; }else{ deltac[n*4+a]=cS[a][(Nas[a][0]-2*n)]-cA[a][(Nas[a][0]-2*n)]; } 
   
   if(2*n>=Nas[a][1]){ deltacup[n*4+a]=cSup[a][(-Nas[a][1]+2*n)]+cAup[a][(-Nas[a][1]+2*n)]; }else{ deltacup[n*4+a]=cSup[a][(Nas[a][1]-2*n)]-cAup[a][(Nas[a][1]-2*n)]; }
    
  }/*n*/
}/*a*/


delete [] fS;
delete [] fA;
delete [] fSup;
delete [] fAup;

for(a=0;a<=3;a++){ delete [] cS[a]; 
                   delete [] cA[a];
                   delete [] cSup[a]; 
                   delete [] cAup[a];
}//a
  
}

/*Function to compute the norm squared of a vector of dimension N00*/

cl_N norma2(cl_N * E/*N00*/,long int N00){
  long int n;
  cl_N r;
  
  r=complex(0,0);
  
  for(n=0;n<=N00-1;n++){ r+=E[n]*conjugate(E[n]); }/*n*/
  
  return r;
}

/* A function giving the value of the A_{a0}A^{a0}=... constraint at a0=a */

cl_N AAfunc0(long int a,cl_N mt[4],cl_N mh[4]){
  
  cl_N II=complex(0,1);
  long int j;
  cl_N r;
  
    r=II;
  
  for(j=0;j<=3;j++){ r*=(mt[a]-mh[j]);  }
  
  for(j=0;j<=3;j++){ if(j==a){ ; }else{ r/=(mt[a]-mt[j]); }  }
  
  return r;
  
}

/*A function to compute some important Delta dependent quantities, needed to be updated at each iteration */

void Ca0func(cl_N Delta,cl_N * c[4],cl_N * cf[4], cl_N Mt[4],cl_N Mhat0[4],cl_N Mhat[4],cl_N A[4],cl_N Af[4],cl_N AA[4][4],
	   cl_N B[4],cl_N BB[4][4],cl_N alfa[4][4],cl_N g,float_format_t prec){
  
long int a,i;  
/* Mhat's update */
 
Mhat[0]=Mhat0[0]+Delta/cl_float(2,prec);
Mhat[1]=Mhat0[1]+Delta/cl_float(2,prec);
Mhat[2]=Mhat0[2]-Delta/cl_float(2,prec);
Mhat[3]=Mhat0[3]-Delta/cl_float(2,prec);  

/*A_a*/

for(a=0;a<=3;a++) { A[a]=VolinAfunc(a,Mt,Mhat); }

/*A^a*/


for(a=0;a<=3;a++){ Af[a]=AAfunc0(a,Mt,Mhat)/A[a]; }


/*other auxiliary 4x4 matrices e.g B_{a|i}, A_a A^i, etc.*/
/*B_i=1 convention is assumed*/

for(a=0;a<=3;a++){ 
  for(i=0;i<=3;i++){  alfa[a][i]=Mhat[i]-Mt[a];
                      
  BB[a][i]=complex(0,1)*A[a]*B[i]/(Mt[a]-Mhat[i]);
  AA[a][i]=A[a]*Af[i];
 }/*i*/
}/*a*/
 
/* new values of c_{a,0} and c^a_0  */ 
 
for(a=0;a<=3;a++){ c[a][0]=A[a]/expt(g,Mt[a]);
                  cf[a][0]=Af[a]/expt(g,1-Mt[a]); 
}/*a*/  
  
  
} /* Ca0func ends */

/***********************************************************/





/**********************************************************************************************************/
/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
/**********************************************************************************************************/
void T5funci(cl_N * T5,cl_N * alfaais,cl_N * m1p4k/*lmax+2*/,long int lmax){
  long int a,l,m,ref,lmp1;
  lmp1=lmax+1;
  
  for(a=0;a<=3;a++){
  
  for(l=2;l<=lmax;l++){ ref=2*l-3;
  for(m=0;m<=2*lmax-3;m++){ 
        if(m<=ref) {T5[((m*lmp1+l)*4+a)]=alfaais[((2*l-m-1)*4+a)]*m1p4k[l];} else{ T5[((m*lmp1+l)*4+a)]=complex(0,0); } 
    }/*for-m */  
  }/* for l */
/*l=0,1  */  
for(l=0;l<=1;l++){ 
  for(m=0;m<=2*lmax-3;m++){ 
         T5[((m*lmp1+l)*4+a)]=complex(0,0); 
    }/*for-m */  
  }/* for l */  
  
  //}/*for i */
}/* for a  */

}




/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void T3funci(cl_N * T3,cl_N * alfaais,cl_N * m1p4k/*lmax+2*/,long int lmax){
  long int a,l;
  for(a=0;a<=3;a++){
    
  for(l=0;l<=lmax;l++){ T3[(l*4+a)]=alfaais[((2*l+1)*4+a)]*m1p4k[l]; }
   
  }
} 

/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void S1nfunc2i(cl_N * S1n,cl_N * alfaais,cl_N * m1p4k/*NQ+2*/,long int NQ0){
  long int a,n;
  for(a=0;a<=3;a++){
  for(n=0;n<=NQ0;n++){
S1n[(n*4+a)]=alfaais[((2*n)*4+a)]*m1p4k[n];    
  }/*n*/
  }/*a*/
}                   
/**********/


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void S32func2i(cl_N * S32,cl_N * alfaais,cl_N * m1p4k/*NQ+2*/,long int NQ0){
  long int a,i,n,k,kNQm1;
  kNQm1=2*NQ0-1;
  for(a=0;a<=3;a++){
      for(n=0;n<=NQ0;n++){
	for(k=0;k<=2*NQ0-2;k++){
	  if(k<=2*n-2){ S32[((n*kNQm1+k)*4+a)]=alfaais[((2*n-k-1)*4+a)]*m1p4k[n]; }else{ S32[((n*kNQm1+k)*4+a)]=complex(0,0); }
	}/*k*/
      }/*n*/
  }/*a*/
}                

/*******************************************************************************/

/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void F2m2i(long int i,long int m,cl_N * F2m,cl_N AA[4][4],cl_N BB[4][4],cl_N * b,cl_N * q/*4*4*(NQ+1)*/,cl_N * S1n,cl_N * S1/*(NQ-1)*(NQ+1)*/,cl_N * S31/*(2*NQ-1)*(NQ-1)*/,cl_N * S32,long int NQ0,long int NQmax){
  long int a,b0,n,k,j,NQm1,kNQm1;
  NQm1=NQmax-1;
  kNQm1=2*NQ0-1;
  cl_N S0s,S1s,S2s,S3s,S3spart,Sq;
  
  for(a=0;a<=3;a++){  F2m[a]=complex(0,0);
    for(b0=0;b0<=3;b0++){
      
      S0s=complex(0,0);
    if(m>=2){ for(n=1;n<=m-1;n++){S0s+=b[(n*4+b0)]*q[((m-n)*4+a)*4+b0];    }/*n*/}else{}  
    
    Sq=complex(0,0);
    
    for(n=0;n<=m;n++){
      
      S1s=S1n[(n*4+b0)];
      
      S2s=complex(0,0);
      
      for(j=0;j<=n-2;j++){
	S2s+=b[((j+1)*4+b0)]*S1[n*NQm1+j];
      }/*j*/
     
     S3s=complex(0,0);
      
      for(k=0;k<=2*n-2;k++){ S3spart=complex(0,0);
	for(j=0;j<=(ldiv(k+1,2).quot)-1;j++){ 
	  S3spart+=b[((j+1)*4+b0)]*S31[k*NQm1+j];	  
	}
S3s+=S32[((n*kNQm1+k)*4+b0)]*S3spart;
}/*k*/
      
Sq+=q[((m-n)*4+a)*4+b0]*(S1s+S2s+S3s);	
      
}/*n*/


F2m[a]+=AA[a][b0]*BB[b0][i]*(S0s+Sq);      
           
    }/*b0*/
    
     }/*a*/
   
}    



/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void F1ps2i(long int i,long int l,cl_N * F12l,cl_N BB[4][4],cl_N  alfa[4][4],cl_N * b,cl_N * T1/*(lmax+1)*(lmax-1)*/,cl_N * T2/*(lmax+1)*(lmax-1) */,cl_N * T3,cl_N * T41/*lmax*(2*lmax-2)*/,cl_N * T5,long int NQ0,long int NQmax){
  long int a,k,m,ref;
  cl_N mI=complex(0,-1);
  cl_N T1s,T2s,T3s,T4s,T4spart;
  long int lmax=NQmax;
  long int lmp1=NQ0+1;
  long int lmm1=lmax-1;

    for(a=0;a<=3;a++){ 
      
      T1s=complex(0,0);
      for(k=0;k<=l-2;k++){ T1s+=b[((k+1)*4+a)]*T1[l*lmm1+k]; }/*k*/
	
T2s=complex(0,0);	
for(k=0;k<=l-2;k++){ T2s+=b[((k+1)*4+a)]*T2[l*lmm1+k]; }/*k*/
T2s*=alfa[a][i];	

T3s=b[(0*4+a)]*T3[(l*4+a)];

T4s=complex(0,0);

for(m=0;m<=2*l-3;m++){ ref=ldiv(m,2L).quot;  
T4spart=complex(0,0);
  for(k=0;k<=ref;k++){ T4spart+=b[((k+1)*4+a)]*T41[m*lmax+k]; }/*k*/
T4s+=T5[((m*lmp1+l)*4+a)]*T4spart;
}/*m*/


F12l[a]=mI*BB[a][i]*(T1s+T2s+T3s+T4s);	

    }/*a*/

}  

/* Function to compute important auxiliary arrays for solving the linear problems for b_{a|i,n} */

void AlfatoT2i(long int i,cl_N alfa[4][4],cl_N * alfaais,cl_N * T3,cl_N * T5,cl_N * S1n,cl_N * S32/*4*4*(NQ0+1)*(NQ0-1)*/
,cl_N * m1p4k,long int NQ0,float_format_t prec){

  long int a,m;

for(a=0;a<=3;a++){

    for(m=0;m<=2*NQ0+1;m++){ alfaais[(m*4+a)]=cbinomial(alfa[a][i],m,prec); 
      
    }/*m*/
 
}/*a*/


T5funci(T5,alfaais,m1p4k,NQ0);


T3funci(T3,alfaais,m1p4k,NQ0);


S1nfunc2i(S1n,alfaais,m1p4k,NQ0);


S32func2i(S32,alfaais,m1p4k,NQ0);


}




/* Function to compute the 4x4 matrices of the recursive linear problems for b_{a|i,n} */
void totalscTmaker2i(long int i,cl_N AA[4][4],cl_N B[4][4],cl_N alfa[4][4],long int NQ0,cl_N * scT){
long int m;
long int a,b0;
cl_N II=complex(0,1);
for(m=0;m<=NQ0;m++){
for(a=0;a<=3;a++){
for(b0=0;b0<=3;b0++){  if(m==0){scT[((m)*4+a)*4+b0]=complex(0,0);}
else{  if(a==b0){scT[((m)*4+a)*4+b0]=AA[a][b0]*B[b0][i]-II*B[a][i]*(2*m-alfa[a][i]);  }
       else{ scT[((m)*4+a)*4+b0]=AA[a][b0]*B[b0][i]; }     }
}/*b0*/
}/*a*/
} /*m*/ 
}  


/* The function computing the Q_i[u] functions and the difference of P_a[u] and P_a[u] computed back from the gluing equations at the sampling points of the cut */
void QconstructorUJ2i(cl_N * deltaP,cl_N * deltaPt,cl_N * deltaPup,cl_N * deltaPtup,cl_N * Qlower,cl_N * Qtlower,cl_N * Qupper,cl_N * Qtupper,cl_N AA[4][4],cl_N BB[4][4],cl_N alfa[4][4],cl_N * c[4],cl_N * cf[4],cl_N * sigmasub[4],cl_N * sigmasup[4],cl_N * PaT,cl_N * PtaT,cl_N * PfaT,cl_N * PftaT,cl_N * PujT,cl_N * PfujT,cl_N * T1,cl_N * T2,cl_N * T41,cl_N * T3[4],cl_N * T5[4],cl_N * S1n[4],cl_N * S1,cl_N * S31,cl_N * S32[4],cl_N * TuANIn[4],cl_N * TuMiMaNI,cl_N * m1p4k,cl_N g,long int N00,long int NQ0[4],long int NQmax,unsigned long int * NI,unsigned long int NImax,long int lc,cl_N * scTm[4],float_format_t prec){
  
 
  long int i,a,m,n,k,j,b0;
  

/* fixing c^1_n from the P.P=0 constraint*/


for(n=1;n<=N00;n++){  cf[0][n]=cfixalo(1,3,n,c,cf)/c[0][0];

}/*n*/

/*computation of the 1/u coefficients of P_a[u] and P^a[u] -> ksub and ksup */

cl_N * ksub[4];
cl_N * ksup[4];
for(a=0;a<=3;a++){ ksub[a]=new cl_N[NQmax+1];  
                   ksup[a]=new cl_N[NQmax+1]; }
	   
		   
kanfunc2(ksub,c,sigmasub,NQmax,N00);
kanfunc2(ksup,cf,sigmasup,NQmax,N00);


/*  Determination of q^ab_n ~ ksup[a]*ksup[b]  */

cl_N * q=new cl_N[4*4*(1+NQmax)];


for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=0;n<=NQmax;n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*ksup[b0][n-m];  }/*m*/
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   

for(a=0;a<=3;a++){ delete [] ksub[a];  
                   delete [] ksup[a]; }


/* solving the subsequent linear problems for b_{a|i,n}=b[i][n*4+a] */

cl_N * b[4];

for(i=0;i<=3;i++) { b[i]=new cl_N[4*(NQ0[i]+1)]; }/*i*/ 

 
  for(i=0;i<=3;i++){
for(m=0;m<=NQ0[i];m++){
  for(a=0;a<=3;a++){
 if(m==0){ b[i][(m*4+a)]=complex(1,0); } else{ b[i][(m*4+a)]=complex(0,0); }  
}/*a*/
}/*m*/
}/*i*/		   
		   
cl_N x4[4];
cl_N v4[4];
cl_N M4[4][4];


cl_N * f1=new cl_N[4];
cl_N * f2=new cl_N[4];

 long int l0;
 

  
for(i=0;i<=3;i++){  
  
for(m=1;m<=NQ0[i];m++){

F2m2i(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);  

F1ps2i(i,m,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);

  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector*/
    
    for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];    }/*b0*/ }/*a*/ /* matrix */
   
/* solving the 4x4 linear equations*/      

linsolvepE(&M4[0][0],&v4[0],x4,4); 


/* filling the up b[i] with the previously computed values */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];} 
  }/*m*/
}/*i*/


delete [] q; 
delete [] f1;
delete [] f2;
		   
/* Determination of Q_{a|i}[u] */	

cl_N * Q[4][4];
cl_N * Qup[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q[a][i]=new cl_N[lc];  Qup[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

/*Q_{a|i} at the large u samplig points */
for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  Q[a][i][k]=complex(0,0);
  for(n=0;n<=NQ0[i];n++){  Q[a][i][k]+=b[i][(n*4+a)]*TuANIn[i][n*lc+k]; }/*n*/
                       Q[a][i][k]*=(BB[a][i]*TuMiMaNI[(k*4+a)*4+i]);   
}/*k*/
}/*a*/
}/*i*/



for(i=0;i<=3;i++) { delete [] b[i]; }/*i*/


/* Construction of the tables of P_a(u_A+I*n)  and P^a(u_A+I*n) */

cl_N * Puj[4];
cl_N * Pfuj[4];

for(a=0;a<=3;a++){ Puj[a]=new cl_N[NImax*lc];
                   Pfuj[a]=new cl_N[NImax*lc];    }

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NImax-1;n++){ Puj[a][n*lc+k]=complex(0,0);
                          Pfuj[a][n*lc+k]=complex(0,0);
        for(m=0;m<=N00;m++){  Puj[a][n*lc+k]+=c[a][m]*PujT[((m*NImax+n)*lc+k)*4+a];
	                     Pfuj[a][n*lc+k]+=cf[a][m]*PfujT[((m*NImax+n)*lc+k)*4+a];
	  	}/*m*/			  
    }/*n*/
  }/*k*/
}/*a*/		   



cl_N vector2[4];/*auxiliary vector*/


/*coming down to the real cut */

for(k=0;k<lc;k++){
for(i=0;i<=3;i++){
for(n=(NI[i]-1);n>=0;n--){
    for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q[a][i][k]+=Pfuj[b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q[a][i][k]*=Puj[a][n*lc+k];
    Q[a][i][k]+=vector2[a];  
      }/*a*/
    }/*i*/
  }/*k*/
    }/*n*/
    
for(a=0;a<=3;a++){ delete [] Puj[a];
                   delete [] Pfuj[a];     }
/********************************************/

/*values of P_a[u], P^a[u] and Ptilde_a[u], Ptilde^a[u] at the sampling points of the cut */

cl_N * P[4];
cl_N * Pt[4];
cl_N * Pf[4];
cl_N * Pft[4];

for(a=0;a<=3;a++){ P[a]=new cl_N[lc];
                   Pt[a]=new cl_N[lc];
		   Pf[a]=new cl_N[lc];
		   Pft[a]=new cl_N[lc];   }
		   
for(a=0;a<=3;a++){ 
  for(k=0;k<=lc-1;k++){  P[a][k]=complex(0,0);  
                       Pt[a][k]=complex(0,0);
		       Pf[a][k]=complex(0,0);
		       Pft[a][k]=complex(0,0);
	for(n=0;n<=N00;n++){  P[a][k]+=c[a][n]*PaT[(n*lc+k)*4+a];
	                     Pt[a][k]+=c[a][n]*PtaT[(n*lc+k)*4+a];
			     Pf[a][k]+=cf[a][n]*PfaT[(n*lc+k)*4+a];
			     Pft[a][k]+=cf[a][n]*PftaT[(n*lc+k)*4+a];
                        }/*n*/	       
  }/*k*/
}/*a*/		   
		   
/*computing Q_i[u] at the cut */ 

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qlower[k*4+i]-=Pf[a][k]*Q[a][i][k];  }
}}


for(a=0;a<=3;a++){ delete [] Pf[a];   }
    
/*computing Qtilde_i[u] at the cut */    
    
for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qtlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qtlower[k*4+i]-=Pft[a][k]*Q[a][i][k];  }
}}
 
for(a=0;a<=3;a++){ delete [] Pft[a];   }


/*computing Q^i[u] and Qtilde^i[u] at the cut */ 


cl_N * Qe=new cl_N[4*4];
cl_N * Qa=new cl_N[4*4];


 
for(k=0;k<lc;k++){
  for(a=0;a<=3;a++){ 
    for(i=0;i<=3;i++){ Qa[a*4+i]=Q[a][i][k];    }/*i*/
  }/*a*/
  

inversetransposeE(Qa,Qe,4);

   
for(i=0;i<=3;i++){  Qupper[k*4+i]=complex(0,0);  Qtupper[k*4+i]=complex(0,0);
  
for(a=0;a<=3;a++){ Qupper[k*4+i]-=P[a][k]*Qe[a*4+i];  Qtupper[k*4+i]-=Pt[a][k]*Qe[a*4+i];  Qup[a][i][k]=complex(-1,0)*Qe[a*4+i]; }/*a*/
}/*i*/
}/*k*/

for(a=0;a<=3;a++){ delete [] P[a];   }

for(a=0;a<=3;a++){ delete [] Pt[a];   }

/* Computing the gluing constants */

cl_N alfa1;
cl_N alfa3;

alfa1=complex(0,0);
alfa3=complex(0,0);

for(k=0;k<=lc-1;k++){ alfa1+=Qlower[k*4+0]/conjugate(Qupper[k*4+1]);  alfa3+=Qlower[k*4+2]/conjugate(Qupper[k*4+3]); }/*k*/
for(k=0;k<=lc-1;k++){ alfa1+=Qtlower[k*4+0]/conjugate(Qtupper[k*4+1]);  alfa3+=Qtlower[k*4+2]/conjugate(Qtupper[k*4+3]); }/*k*/  
  
  alfa1/=(2*lc);
  alfa3/=(2*lc);
  
alfa1=complex(0,imagpart(alfa1));  
alfa3=complex(0,imagpart(alfa3)); 
  
/*computation of P_a[u]-PG_a[u] at the sampling points, with PG_a[u] being P_a[u] computed back from RHS of the gluing equations and its upper index version */

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){   
    
    deltaPup[k*4+a]=Qup[a][0][k]*(Qlower[k*4+0]-alfa1*conjugate(Qupper[k*4+1]))+Qup[a][1][k]*(Qlower[k*4+1]+alfa1*conjugate(Qupper[k*4+0]))
    +Qup[a][2][k]*(Qlower[k*4+2]-alfa3*conjugate(Qupper[k*4+3]))+Qup[a][3][k]*(Qlower[k*4+3]+alfa3*conjugate(Qupper[k*4+2]));
    
    deltaPtup[k*4+a]=Qup[a][0][k]*(Qtlower[k*4+0]-alfa1*conjugate(Qtupper[k*4+1]))+Qup[a][1][k]*(Qtlower[k*4+1]+alfa1*conjugate(Qtupper[k*4+0]))
    +Qup[a][2][k]*(Qtlower[k*4+2]-alfa3*conjugate(Qtupper[k*4+3]))+Qup[a][3][k]*(Qtlower[k*4+3]+alfa3*conjugate(Qtupper[k*4+2]));
    
                          
deltaP[k*4+a]=complex(0,0);    
    
deltaP[k*4+a]-= (Q[a][0][k]*(Qupper[k*4+0]+conjugate(Qlower[k*4+1]/alfa1))+Q[a][1][k]*(Qupper[k*4+1]-conjugate(Qlower[k*4+0]/alfa1))
+Q[a][2][k]*(Qupper[k*4+2]+conjugate(Qlower[k*4+3]/alfa3))+Q[a][3][k]*(Qupper[k*4+3]-conjugate(Qlower[k*4+2]/alfa3)));

deltaPt[k*4+a]=complex(0,0);    
    
deltaPt[k*4+a]-= (Q[a][0][k]*(Qtupper[k*4+0]+conjugate(Qtlower[k*4+1]/alfa1))+Q[a][1][k]*(Qtupper[k*4+1]-conjugate(Qtlower[k*4+0]/alfa1))
+Q[a][2][k]*(Qtupper[k*4+2]+conjugate(Qtlower[k*4+3]/alfa3))+Q[a][3][k]*(Qtupper[k*4+3]-conjugate(Qtlower[k*4+2]/alfa3)));
    
  }/*k*/
}/*a*/

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  delete [] Q[a][i]; delete [] Qup[a][i];
 }/*a*/
}/*i*/

delete [] Qe;
delete [] Qa;		   
		   
    
}/* QconstructorUJ2i ends */

/*The same as the previous function, but it does the computations two times with NI imaginary cutoff and with NI-1 imaginary cutoff */

void QconstructorUJ2icm(cl_N * deltaP,cl_N * deltaPt,cl_N * deltaPup,cl_N * deltaPtup,cl_N * deltaPv,cl_N * deltaPtv,cl_N * deltaPupv,cl_N * deltaPtupv,cl_N * Qlower,cl_N * Qtlower,cl_N * Qupper,cl_N * Qtupper,cl_N AA[4][4],cl_N BB[4][4],cl_N alfa[4][4],cl_N * c[4],cl_N * cf[4],cl_N * sigmasub[4],cl_N * sigmasup[4],cl_N * PaT,cl_N * PtaT,cl_N * PfaT,cl_N * PftaT,cl_N * PujT,cl_N * PfujT,cl_N * T1,cl_N * T2,cl_N * T41,cl_N * T3[4],cl_N * T5[4],cl_N * S1n[4],cl_N * S1,cl_N * S31,cl_N * S32[4],cl_N * TuANIn[4],cl_N * TuMiMaNI,cl_N * m1p4k,cl_N g,long int N00,long int NQ0[4],long int NQmax,unsigned long int * NI,unsigned long int NImax,long int lc,cl_N * scTm[4],float_format_t prec,cl_N * TuANInshift[4],cl_N * TuMiMaNIshift,cl_N * Qalfadev){
  
 
  long int i,a,m,n,k,j,b0;


for(n=1;n<=N00;n++){  cf[0][n]=cfixalo(1,3,n,c,cf)/c[0][0];

}/*n*/


cl_N * ksub[4];
cl_N * ksup[4];
for(a=0;a<=3;a++){ ksub[a]=new cl_N[NQmax+1];  
                   ksup[a]=new cl_N[NQmax+1]; }

kanfunc2(ksub,c,sigmasub,NQmax,N00);
kanfunc2(ksup,cf,sigmasup,NQmax,N00);

cl_N * q=new cl_N[4*4*(1+NQmax)];




for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=0;n<=NQmax;n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*ksup[b0][n-m];  }/*m*/
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   

/**/

for(a=0;a<=3;a++){ delete [] ksub[a];  
                   delete [] ksup[a]; }


cl_N * b[4];

for(i=0;i<=3;i++) { b[i]=new cl_N[4*(NQ0[i]+1)]; }/*i*/ 

  
 
  for(i=0;i<=3;i++){
for(m=0;m<=NQ0[i];m++){
  for(a=0;a<=3;a++){
 if(m==0){ b[i][(m*4+a)]=complex(1,0); } else{ b[i][(m*4+a)]=complex(0,0); }  
}/*a*/
}/*m*/
}/*i*/		   
		   
cl_N x4[4];
cl_N v4[4];
cl_N M4[4][4];

cl_N * f1=new cl_N[4];
cl_N * f2=new cl_N[4];


 long int l0;
 
  
for(i=0;i<=3;i++){  
  
for(m=1;m<=NQ0[i];m++){

F2m2i(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);  

F1ps2i(i,m,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);

  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
    for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];    }/*b0*/ }/*a*/ /* matrix */
    

linsolvepE(&M4[0][0],&v4[0],x4,4); 


/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
  
}/*m*/

}/*i*/

delete [] q; 
delete [] f1;
delete [] f2;
		   
/* Q is for Q_{a|i}[u] computed with NI imaginary cutoff*/

cl_N * Q[4][4];
cl_N * Qup[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q[a][i]=new cl_N[lc];  Qup[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

/* Q1 is for Q_{a|i}[u] computed with NI-1 imaginary cutoff*/

cl_N * Q1[4][4];
cl_N * Q1up[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q1[a][i]=new cl_N[lc];  Q1up[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

cl_N * Q2[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q2[a][i]=new cl_N[lc];  
 }/*a*/
}/*i*/



/* The values at large u sampling points */
for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  Q[a][i][k]=complex(0,0);
  for(n=0;n<=NQ0[i];n++){  Q[a][i][k]+=b[i][(n*4+a)]*TuANIn[i][n*lc+k]; }/*n*/
                       Q[a][i][k]*=(BB[a][i]*TuMiMaNI[(k*4+a)*4+i]);   
}/*k*/
}/*a*/
}/*i*/


for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  Q1[a][i][k]=complex(0,0);
  for(n=0;n<=NQ0[i];n++){  Q1[a][i][k]+=b[i][(n*4+a)]*TuANInshift[i][n*lc+k]; }/*n*/
                       Q1[a][i][k]*=(BB[a][i]*TuMiMaNIshift[(k*4+a)*4+i]);   
}/*k*/
}/*a*/
}/*i*/



for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  if(i==a){ Q2[a][i][k]=complex(1,0); }else{ Q2[a][i][k]=complex(0,0); }
  
}/*k*/
}/*a*/
}/*i*/



for(i=0;i<=3;i++) { delete [] b[i]; }/*i*/



cl_N * Puj[4];
cl_N * Pfuj[4];

for(a=0;a<=3;a++){ Puj[a]=new cl_N[NImax*lc];
                   Pfuj[a]=new cl_N[NImax*lc];    }

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NImax-1;n++){ Puj[a][n*lc+k]=complex(0,0);
                          Pfuj[a][n*lc+k]=complex(0,0);
        for(m=0;m<=N00;m++){  Puj[a][n*lc+k]+=c[a][m]*PujT[((m*NImax+n)*lc+k)*4+a];
	                     Pfuj[a][n*lc+k]+=cf[a][m]*PfujT[((m*NImax+n)*lc+k)*4+a];
	  	}/*m*/			  
    }/*n*/
  }/*k*/
}/*a*/		   



cl_N vector2[4];/*auxiliary vector*/


/**************************************************/


for(k=0;k<lc;k++){
for(i=0;i<=3;i++){
for(n=(NI[i]-1);n>=(NI[i]-1);n--){
    for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q[a][i][k]+=Pfuj[b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q[a][i][k]*=Puj[a][n*lc+k];
    Q[a][i][k]+=vector2[a];  
      }/*a*/
    }/*i*/
  }/*k*/
    }/*n*/
    
    
/*************************************/    
int s;    


for(k=0;k<lc;k++){
  
  
  
for(i=0;i<=3;i++){
for(n=(NI[i]-2);n>=0;n--){
    for(j=0;j<=3;j++){ vector2[j]=Q2[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q2[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q2[a][i][k]+=Pfuj[b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q2[a][i][k]*=Puj[a][n*lc+k];
    Q2[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
  
 /**** Q[a][i] and Q1[a][i] at u+I/2  ***/ 
  
 for(i=0;i<=3;i++){  for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; }
  for(a=0;a<=3;a++){ 
    
    Q[a][i][k]=complex(0,0);
  
  for(s=0;s<=3;s++){  Q[a][i][k]+=Q2[a][s][k]*vector2[s];  }//  
  }//i  
}//a


/**********/

for(i=0;i<=3;i++){  for(j=0;j<=3;j++){ vector2[j]=Q1[j][i][k]; }
  for(a=0;a<=3;a++){ 
    
    Q1[a][i][k]=complex(0,0);
  
  for(s=0;s<=3;s++){  Q1[a][i][k]+=Q2[a][s][k]*vector2[s];  }//  
  }//i  
}//a
  /**********************/
    }/*for-k ends*/    
 /***********************/   
    
    
 for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q2[a][i]; }/*i*/
}/*a*/    
    
    
for(a=0;a<=3;a++){ delete [] Puj[a];
                   delete [] Pfuj[a];     }
/********************************************/

cl_N * P[4];
cl_N * Pt[4];
cl_N * Pf[4];
cl_N * Pft[4];

for(a=0;a<=3;a++){ P[a]=new cl_N[lc];
                   Pt[a]=new cl_N[lc];
		   Pf[a]=new cl_N[lc];
		   Pft[a]=new cl_N[lc];   }
		   
for(a=0;a<=3;a++){ 
  for(k=0;k<=lc-1;k++){  P[a][k]=complex(0,0);  
                       Pt[a][k]=complex(0,0);
		       Pf[a][k]=complex(0,0);
		       Pft[a][k]=complex(0,0);
	for(n=0;n<=N00;n++){  P[a][k]+=c[a][n]*PaT[(n*lc+k)*4+a];
	                     Pt[a][k]+=c[a][n]*PtaT[(n*lc+k)*4+a];
			     Pf[a][k]+=cf[a][n]*PfaT[(n*lc+k)*4+a];
			     Pft[a][k]+=cf[a][n]*PftaT[(n*lc+k)*4+a];
                        }/*n*/	       
  }/*k*/
}/*a*/		   


/* Array for Q_i[u] from the decreased cutoff NI-1*/
cl_N * Qlowerm1=new cl_N[4*lc];     

/* Computing Q^i[u], Qtilde^i[u] with NI and with NI-1*/
for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qlower[k*4+i]=complex(0,0); Qlowerm1[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qlower[k*4+i]-=Pf[a][k]*Q[a][i][k]; Qlowerm1[k*4+i]-=Pf[a][k]*Q1[a][i][k]; }
}}


for(a=0;a<=3;a++){ delete [] Pf[a];   }
    

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qtlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qtlower[k*4+i]-=Pft[a][k]*Q[a][i][k];  }
}}
 
for(a=0;a<=3;a++){ delete [] Pft[a];   }


/*Array for Q^i[u] from the decreased cutoff NI-1*/
cl_N * Qupperm1=new cl_N[4*lc];

cl_N * Qe=new cl_N[4*4];
cl_N * Qa=new cl_N[4*4];


/* Computing Q^i[u], Qtilde^i[u] with NI */
for(k=0;k<lc;k++){
  for(a=0;a<=3;a++){ 
    for(i=0;i<=3;i++){ Qa[a*4+i]=Q[a][i][k];    }/*i*/
  }/*a*/
 
inversetransposeE(Qa,Qe,4);


   
for(i=0;i<=3;i++){  Qupper[k*4+i]=complex(0,0);  Qtupper[k*4+i]=complex(0,0);
  
for(a=0;a<=3;a++){ Qupper[k*4+i]-=P[a][k]*Qe[a*4+i];  Qtupper[k*4+i]-=Pt[a][k]*Qe[a*4+i];  Qup[a][i][k]=complex(-1,0)*Qe[a*4+i]; }/*a*/
}/*i*/
}/*k*/

/* Computing Q^i[u], Qtilde^i[u] with NI-1 */
for(k=0;k<lc;k++){
  for(a=0;a<=3;a++){ 
    for(i=0;i<=3;i++){ Qa[a*4+i]=Q1[a][i][k];    }/*i*/
  }/*a*/
  

inversetransposeE(Qa,Qe,4);

   
for(i=0;i<=3;i++){  Qupperm1[k*4+i]=complex(0,0); 
  
for(a=0;a<=3;a++){ Qupperm1[k*4+i]-=P[a][k]*Qe[a*4+i];  Q1up[a][i][k]=complex(-1,0)*Qe[a*4+i]; }/*a*/
}/*i*/
}/*k*/


for(a=0;a<=3;a++){ delete [] P[a];   }

for(a=0;a<=3;a++){ delete [] Pt[a];   }


/* Computating the norm of the differences of Q_i[u] with NI and Q_i[u] with NI-1 and the same for the upper case counterpart */

Qalfadev[0]=complex(0,0);

for(i=0;i<=3;i++){
  for(k=0;k<=lc-1;k++){ Qalfadev[0]+=(Qlower[k*4+i]-Qlowerm1[k*4+i])*conjugate(Qlower[k*4+i]-Qlowerm1[k*4+i]); 
                        Qalfadev[0]+=(Qupper[k*4+i]-Qupperm1[k*4+i])*conjugate(Qupper[k*4+i]-Qupperm1[k*4+i]); 
  }//k
}//i

Qalfadev[0]=sqrt(Qalfadev[0]);




/* Computing the gluing constants */

cl_N alfa1;
cl_N alfa3;

alfa1=complex(0,0);
alfa3=complex(0,0);

for(k=0;k<=lc-1;k++){ alfa1+=Qlower[k*4+0]/conjugate(Qupper[k*4+1]);  alfa3+=Qlower[k*4+2]/conjugate(Qupper[k*4+3]); }/*k*/
for(k=0;k<=lc-1;k++){ alfa1+=Qtlower[k*4+0]/conjugate(Qtupper[k*4+1]);  alfa3+=Qtlower[k*4+2]/conjugate(Qtupper[k*4+3]); }/*k*/  
  
  alfa1/=(2*lc);
  alfa3/=(2*lc);
  
alfa1=complex(0,imagpart(alfa1));  
alfa3=complex(0,imagpart(alfa3)); 
  


/*******************************************************/

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){   /*computing P_a[u]-P(from gluing)_a[u], and  P^a[u]-P(from gluing)^a[u], and their tilded versions*/
    
    deltaPup[k*4+a]=Qup[a][0][k]*(Qlower[k*4+0]-alfa1*conjugate(Qupper[k*4+1]))+Qup[a][1][k]*(Qlower[k*4+1]+alfa1*conjugate(Qupper[k*4+0]))
    +Qup[a][2][k]*(Qlower[k*4+2]-alfa3*conjugate(Qupper[k*4+3]))+Qup[a][3][k]*(Qlower[k*4+3]+alfa3*conjugate(Qupper[k*4+2]));
    
    deltaPtup[k*4+a]=Qup[a][0][k]*(Qtlower[k*4+0]-alfa1*conjugate(Qtupper[k*4+1]))+Qup[a][1][k]*(Qtlower[k*4+1]+alfa1*conjugate(Qtupper[k*4+0]))
    +Qup[a][2][k]*(Qtlower[k*4+2]-alfa3*conjugate(Qtupper[k*4+3]))+Qup[a][3][k]*(Qtlower[k*4+3]+alfa3*conjugate(Qtupper[k*4+2]));
  
	
deltaP[k*4+a]=complex(0,0);    
    
deltaP[k*4+a]-= (Q[a][0][k]*(Qupper[k*4+0]+conjugate(Qlower[k*4+1]/alfa1))+Q[a][1][k]*(Qupper[k*4+1]-conjugate(Qlower[k*4+0]/alfa1))
+Q[a][2][k]*(Qupper[k*4+2]+conjugate(Qlower[k*4+3]/alfa3))+Q[a][3][k]*(Qupper[k*4+3]-conjugate(Qlower[k*4+2]/alfa3)));


deltaPt[k*4+a]=complex(0,0);    
    
deltaPt[k*4+a]-= (Q[a][0][k]*(Qtupper[k*4+0]+conjugate(Qtlower[k*4+1]/alfa1))+Q[a][1][k]*(Qtupper[k*4+1]-conjugate(Qtlower[k*4+0]/alfa1))
+Q[a][2][k]*(Qtupper[k*4+2]+conjugate(Qtlower[k*4+3]/alfa3))+Q[a][3][k]*(Qtupper[k*4+3]-conjugate(Qtlower[k*4+2]/alfa3)));

/******************************************************************************************/
/*computing only P(from gluing)_a[u], and  P(from gluing)^a[u], and their tilded versions*/


 deltaPupv[k*4+a]=Qup[a][0][k]*((-1)*alfa1*conjugate(Qupper[k*4+1]))+Qup[a][1][k]*(alfa1*conjugate(Qupper[k*4+0]))
    +Qup[a][2][k]*((-1)*alfa3*conjugate(Qupper[k*4+3]))+Qup[a][3][k]*(alfa3*conjugate(Qupper[k*4+2]));
    
    deltaPtupv[k*4+a]=Qup[a][0][k]*((-1)*alfa1*conjugate(Qtupper[k*4+1]))+Qup[a][1][k]*(alfa1*conjugate(Qtupper[k*4+0]))
    +Qup[a][2][k]*((-1)*alfa3*conjugate(Qtupper[k*4+3]))+Qup[a][3][k]*(alfa3*conjugate(Qtupper[k*4+2]));
  
	
deltaPv[k*4+a]=complex(0,0);    
    
deltaPv[k*4+a]-= (Q[a][0][k]*(conjugate(Qlower[k*4+1]/alfa1))+Q[a][1][k]*((-1)*conjugate(Qlower[k*4+0]/alfa1))
+Q[a][2][k]*(conjugate(Qlower[k*4+3]/alfa3))+Q[a][3][k]*((-1)*conjugate(Qlower[k*4+2]/alfa3)));


deltaPtv[k*4+a]=complex(0,0);    
    
deltaPtv[k*4+a]-= (Q[a][0][k]*(conjugate(Qtlower[k*4+1]/alfa1))+Q[a][1][k]*((-1)*conjugate(Qtlower[k*4+0]/alfa1))
+Q[a][2][k]*(conjugate(Qtlower[k*4+3]/alfa3))+Q[a][3][k]*((-1)*conjugate(Qtlower[k*4+2]/alfa3)));
    
  }/*k*/
}/*a*/

/* computing the L2-norm of the deviations of the ratios "Q/Qtilde" from the guing constants*/

Qalfadev[1]=complex(0,0);

cl_N aux=complex(0,0);

for(k=0;k<=lc-1;k++){ 
  
Qalfadev[1]+=(alfa1-(Qlower[k*4+0]/conjugate(Qupper[k*4+1])))*conjugate(alfa1-(Qlower[k*4+0]/conjugate(Qupper[k*4+1]))); 
Qalfadev[1]+=(alfa1+(Qlower[k*4+1]/conjugate(Qupper[k*4+0])))*conjugate(alfa1+(Qlower[k*4+1]/conjugate(Qupper[k*4+0])));

    aux+=(alfa3-(Qlower[k*4+2]/conjugate(Qupper[k*4+3])))*conjugate(alfa3-(Qlower[k*4+2]/conjugate(Qupper[k*4+3])));
    
    aux+=(alfa3+(Qlower[k*4+3]/conjugate(Qupper[k*4+2])))*conjugate(alfa3+(Qlower[k*4+3]/conjugate(Qupper[k*4+2])));
    
}//k

Qalfadev[1]/=alfa1;

Qalfadev[1]+=aux/alfa3;

Qalfadev[1]=sqrt(Qalfadev[1]);
Qalfadev[1]/=(4*lc);
	   

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  delete [] Q[a][i]; delete [] Qup[a][i]; delete [] Q1[a][i]; delete [] Q1up[a][i];
 }/*a*/
}/*i*/

delete [] Qe;
delete [] Qa;	

delete [] Qlowerm1;
delete [] Qupperm1;
		   
    
}/* function ends */


/* Function computing c_{a,n0} from P and Ptilde computed back from gluing */ 

void QtoEtypeIIcm(long int n0/*-1,..,N0 */,cl_N * Ev,cl_N * deltaPv,cl_N * deltaPtv,cl_N * deltaPupv,cl_N * deltaPtupv,cl_N * CT,cl_N * CU,long int Nas[4][2],cl_N * suA/*lc*/,long int lc,cl_N g){
  
  long int k,a,n;
  
  
  /*auxiliary index vector*/
  long int na[4][2];
  
  for(a=0;a<=3;a++){ na[a][0]=labs(2*n0-Nas[a][0]); na[a][1]=labs(2*n0-Nas[a][1]); }//a
  
  /*symmetric and antisymmetric combinations */
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  cl_N * fSup=new cl_N[4*lc];
  cl_N * fAup=new cl_N[4*lc];
  
  cl_N * cS=new cl_N[4];
  cl_N * cSup=new cl_N[4];
  
  cl_N * cA=new cl_N[4];
  cl_N * cAup=new cl_N[4];
  
  for(a=0;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaPv[k*4+a]+deltaPtv[k*4+a])/2; 
                          fA[k*4+a]=(deltaPv[k*4+a]-deltaPtv[k*4+a])/2/suA[k];
			  
			  fSup[k*4+a]=(deltaPupv[k*4+a]+deltaPtupv[k*4+a])/2;
			  fAup[k*4+a]=(deltaPupv[k*4+a]-deltaPtupv[k*4+a])/2/suA[k];
    }/*k*/
  }/*a*/
  

for(a=0;a<=3;a++){ 
  
  cS[a]=complex(0,0); cSup[a]=complex(0,0);
    
  for(k=0;k<=lc-1;k++){   cS[a]+=fS[k*4+a]*CT[(lc-1-k)*lc+na[a][0]]; cSup[a]+=fSup[k*4+a]*CT[(lc-1-k)*lc+na[a][1]]; 
    
  }/*k*/
  
  cS[a]/=lc;
  cSup[a]/=lc;
    
}/*a*/

for(a=0;a<=3;a++){ 
    
    cA[a]=complex(0,0); cAup[a]=complex(0,0);
    
    if(na[a][0]==0){}else{
    
  for(k=0;k<=lc-1;k++){   cA[a]+=fA[k*4+a]*CU[(lc-1-k)*lc+na[a][0]-1];  
    
  }/*k*/
  
  cA[a]/=lc; 
  cA[a]*=2*complex(0,1)*g;
    }//else ends
    
    
/*******************************/
if(na[a][1]==0){}else{
    
  for(k=0;k<=lc-1;k++){  cAup[a]+=fAup[k*4+a]*CU[(lc-1-k)*lc+na[a][1]-1]; 
    
  }/*k*/
  
  
  cAup[a]/=lc;
  cAup[a]*=2*complex(0,1)*g;
    }//else ends

}/*a*/

for(a=0;a<=3;a++){ 
 
    if(2*n0>=Nas[a][0]){ Ev[a]=cS[a]+cA[a]; }else{ Ev[a]=cS[a]-cA[a]; } 
    
    if(2*n0>=Nas[a][1]){ Ev[4+a]=cSup[a]+cAup[a]; }else{ Ev[4+a]=cSup[a]-cAup[a]; } 
 
}/*a*/


delete [] fS;
delete [] fA;
delete [] fSup;
delete [] fAup;
delete [] cS;
delete [] cSup;
delete [] cA;
delete [] cAup;
  
}


/*************************************************************************/ 
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/ 
 
 
int main( int argc, char * argv[]) {
  
long int i,j,k,n,a,m,b0;  
  
 long int workingprecision=atoi(argv[1]); /*working precision*/

float_format_t precision=float_format(workingprecision); 

long int N0=atoi(argv[2]);/* cutoff parameter c_ {a, n}, c^a_n n = 0, ..., N0 */

long int NCb=atoi(argv[3]);
/*cutoff parameter, the last kept term in the series of Q_ {a | i} is u^(-Mhat[i])/u^(2*NCb) */

long int NIcutoff=atoi(argv[4]);/* pull down cutoff */

long int lc=atoi(argv[5]);/*number of Chebyshev sampling points on the cut */

long int DH=atoi(argv[6]);/* integer derivative parameter \[Epsilon] = 10^(-DH) */

long int precssf=atoi(argv[7]);/* integer quit condition parameter */

long int precDelta=atoi(argv[8]);/* integer quit condition parameter */

long int maxiter=atoi(argv[9]);/* iteration quit condition parameter, if the number of iterations is 
larger than maxiter, than the programm finishes running, and returns
the value, it obtained upto maxiter iterations */

/* The quit condition is as follows: either "number of iterations" > maxiter, 
                                      
 or the simultaneous fullfillment of two conditions is required:
                                         
|Delta^{(n)}-Delta^{(n-1)}|<10^(-precDelta) and 
"the sum of |c_{a,n}-(c_{a,n} from gluing)|^2"<10^(-precssf), 

where Delta^{(n)} is the value of the anomalous part dimension after the n-th iteration.*/



/* filling up the oscillator quantum numbers*/
long int nb[2];
long int na[2];
long int nf[4];

nb[0]=atoi(argv[10]);//n_{b_1}
nb[1]=atoi(argv[11]);//n_{b_2}

nf[0]=atoi(argv[12]);//n_{f_1}
nf[1]=atoi(argv[13]);//n_{f_2}
nf[2]=atoi(argv[14]);//n_{f_3}
nf[3]=atoi(argv[15]);//n_{f_4}

na[0]=atoi(argv[16]);//n_{a_1}
na[1]=atoi(argv[17]);//n_{a_2}

long int multiplicity=atoi(argv[18]); /*multiplicity*/

/*coupling constant*/

cl_N g;

stringstream stext;

stext<<argv[19];

string text=stext.str();

g=cl_float((cl_R) text.c_str(),precision);
/* clear stext stringstream*/
stext.str(std::string());


/*initial value of Delta=amomalous dimension in this code*/
stext<<argv[20];
cl_R Delta0=cl_float((cl_R) stext.str().c_str(),precision);
/* clear stext stringstream*/
stext.str(std::string());

/*The value of N0 for the input c-arrays */
long int N0be=(argc-29)/8;


/*initial values of c_{a,n}*/ 
cl_R * c0[4];
cl_R * cf0[4];

for(a=0;a<=3;a++){ c0[a]=new cl_R[N0be+1];  cf0[a]=new cl_R[N0be+1]; }/*a*/

for(j=0;j<=N0be;j++){    
                         stext<<argv[21+j]; c0[0][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         /* clear stext stringstream*/stext.str(std::string());
                       
			 stext<<argv[22+N0be+j]; c0[1][j]=cl_float((cl_R) stext.str().c_str(),precision);
                         /* clear stext stringstream*/stext.str(std::string());
			 
			 					 		      
		       stext<<argv[23+2*N0be+j]; c0[2][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         /* clear stext stringstream*/stext.str(std::string());
		       
		      
		       stext<<argv[24+3*N0be+j]; c0[3][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         /* clear stext stringstream*/stext.str(std::string());
		       
		       
		       stext<<argv[25+4*N0be+j]; cf0[0][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         /* clear stext stringstream*/stext.str(std::string());
			 
                     
		       stext<<argv[26+5*N0be+j]; cf0[1][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         /* clear stext stringstream*/stext.str(std::string());
			 
		       
		       stext<<argv[27+6*N0be+j]; cf0[2][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         /* clear stext stringstream*/stext.str(std::string());
			 
		      
		       stext<<argv[28+7*N0be+j]; cf0[3][j]=cl_float((cl_R) stext.str().c_str(),precision);
                         /* clear stext stringstream*/stext.str(std::string());
}/*j*/
  
  
  
  
/* some useful constants*/  

cl_N Pi=pi(precision);
cl_R two=cl_float(2,precision);
cl_R egy=cl_float(1,precision);
cl_R half=egy/two;
cl_N II=complex(0,1);

 
 /* "lattice constant" for differentiation */
 const cl_N h=expt(cl_float(10,precision),(-1)*DH);
 

/* Computing zero-mode levels */
long int ni[4]; 

ni[0]=nb[0]-nb[1]-1;
ni[1]=nb[1]-nb[0]+1;
ni[2]=na[1]-na[0]-1;
ni[3]=na[0]-na[1]+1; 

/* some auxiliary integers*/

long int nmin;
long int nmax;

long int imin;
long int imax;
/* auxiliary index */
long int index;

nmax=0;
imax=0;
index=0;

for(i=0;i<=3;i++) { if(ni[i]>0 && index==0){imin=i; nmin=ni[i]; index++; }else{}/*1. if */
                    if(ni[i]>0 && index==1){imax=i; nmax=ni[i];  }else{}/*2. if */}/*i*/

if(nmin>nmax){ index=nmax; nmax=nmin; nmin=index;
               index=imax; imax=imin; imin=index; /* interchange */
                
}else{}


/* Filling Mtilde vector with its integer values */

long int Mtint[4];

Mtint[0]=nf[0]+2;
Mtint[1]=nf[1]+1;
Mtint[2]=nf[2];
Mtint[3]=nf[3]-1;


cl_N Mt[4];

for(a=0;a<=3;a++){ Mt[a]=cl_float((cl_I) Mtint[a],precision);  }

/* The double of the Mtilde is always integer */

long int twiceMt[4];

for(a=0;a<=3;a++){ twiceMt[a]=2*Mtint[a];
  }/*a*/


/*The value of L */
cl_N L=cl_float((cl_I) nf[0]+nf[1]+nf[2]+nf[3]+na[0]+na[1]-nb[0]-nb[1],precision)/two;
long int Lint=(nf[0]+nf[1]+nf[2]+nf[3]+na[0]+na[1]-nb[0]-nb[1])/2;


/* initializing Delta*/
cl_N Delta=cl_float(Delta0,precision);

/* Mhat at g=0 */

cl_N Mhat0[4];

Mhat0[0]=L+cl_float((cl_I) nb[0]+1,precision);
Mhat0[1]=L+cl_float((cl_I) nb[1]+2,precision);
Mhat0[2]=cl_float((cl_I) -1-na[0],precision);
Mhat0[3]=cl_float((cl_I) -na[1],precision);

/* Mhat at the begining. I.e initial value*/
cl_N Mhat[4];

Mhat[0]=L+cl_float((cl_I) nb[0]+1,precision)+Delta/two;
Mhat[1]=L+cl_float((cl_I) nb[1]+2,precision)+Delta/two;
Mhat[2]=cl_float((cl_I) -1-na[0],precision)-Delta/two;
Mhat[3]=cl_float((cl_I) -na[1],precision)-Delta/two;



/* Filling up c_{a,n}'s with their initial values */
long int N0c=max(ldiv(nmax,2).quot+1,N0);

cl_N * c[4];
cl_N * cf[4];


for(a=0;a<=3;a++){ c[a]=new cl_N[N0c+1];
                   cf[a]=new cl_N[N0c+1];   
}



for(a=0;a<=3;a++){
  for(n=0;n<=N0;n++){  
    
    if(n<=N0be){/*if-n */
    
    if(a==0 || a==2){  c[a][n]=II*cl_float(c0[a][n],precision);
                      cf[a][n]=cl_float(cf0[a][n],precision);  }else{
                                                   c[a][n]=cl_float(c0[a][n],precision);
                                                   cf[a][n]=II*cl_float(cf0[a][n],precision);      }/*if-a ends*/
}else{ c[a][n]=complex(0,0);
      cf[a][n]=complex(0,0); }/*if -n ends*/
  }/*n*/
for(n=N0+1;n<=N0c;n++){
    c[a][n]=complex(0,0);
    cf[a][n]=complex(0,0);
  }/*n*/
}/*a*/

/* c0 can be deleted */
for(a=0;a<=3;a++){ delete [] c0[a]; delete [] cf0[a]; }/*a*/
  
  
  
N0=N0c;  
  


/* the A_a*A^a vector */

cl_N AAproduct[4];

AAfunc(Mt,Mhat,AAproduct);

cl_N A[4];
cl_N Af[4];

/* Giving the initial values for A_a, A^a and c_{a,0}, c^a_0 */
for(a=0;a<=3;a++){
A[a]=VolinAfunc(a,Mt,Mhat); 
Af[a]=AAproduct[a]/A[a];

c[a][0]=A[a]/expt(g,Mt[a]);
cf[a][0]=expt(g,Mt[a]-1)*Af[a];
  
}/*a*/


/* A_a*A^b matrix */

cl_N AA[4][4];
for(a=0;a<=3;a++){ 
  for(b0=0;b0<=3;b0++){  AA[a][b0]=A[a]*Af[b0]; 
  }/*b0*/
  }/*a*/


cl_N B[4]; /* B_i vector and its initialization */

for(i=0;i<=3;i++){ B[i]=cl_float(1,precision);  }


/* auxiliary array for computing the dimension of true unknown vector */
long int Nch[4];

/* a=1 case */
Nch[0]=N0;
/* a=2 case */
k=0;
if((Mtint[0]-Mtint[1])%2==0 && (Mtint[0]-Mtint[1])<=2*N0){ k++; }
Nch[1]=N0-k;
/* a=3 case */
k=0;
if((Mtint[0]-Mtint[2])%2==0  && (Mtint[0]-Mtint[2])<=2*N0){ k++; }
if((Mtint[1]-Mtint[2])%2==0  && (Mtint[1]-Mtint[2])<=2*N0){ k++; }
Nch[2]=N0-k;
/* a=4 case */
k=0;
if((Mtint[0]-Mtint[3])%2==0  && (Mtint[0]-Mtint[3])<=2*N0){ k++; }
if((Mtint[1]-Mtint[3])%2==0  && (Mtint[1]-Mtint[3])<=2*N0){ k++; }
if((Mtint[2]-Mtint[3])%2==0  && (Mtint[2]-Mtint[3])<=2*N0){ k++; }
Nch[3]=N0-k;

/* The dimension of the unknown vector */

long int dimV=1+4*N0+Nch[1]+Nch[2]+Nch[3];


/* Fixing gauge-fixed elements of c_{a,n}  */

if( (Mtint[0]-Mtint[1])%2==0 ) { c[1][(Mtint[0]-Mtint[1])/2]=0; }

if( (Mtint[0]-Mtint[2])%2==0 ) { c[2][(Mtint[0]-Mtint[2])/2]=0; }
if( (Mtint[1]-Mtint[2])%2==0 ) { c[2][(Mtint[1]-Mtint[2])/2]=0;  }

if( (Mtint[0]-Mtint[3])%2==0 ) { c[3][(Mtint[0]-Mtint[3])/2]=0;  }
if( (Mtint[1]-Mtint[3])%2==0 ) { c[3][(Mtint[1]-Mtint[3])/2]=0;  }
if( (Mtint[2]-Mtint[3])%2==0 ) { c[3][(Mtint[2]-Mtint[3])/2]=0;  }

/* Fixing c^1_n from P.P=0 */

for(n=1;n<=N0;n++){  cf[0][n]=cfixalo(1,3,n,c,cf)/c[0][0];
}/*n*/




/*  TheB_{a|i} matrix */
cl_N BB[4][4];

cl_N alfa[4][4];
cl_N beta[4][4];

/* and two other auxiliary matrices. See below their definitions. */

for(a=0;a<=3;a++){ 
  for(i=0;i<=3;i++){  alfa[a][i]=Mhat[i]-Mt[a];
                      beta[a][i]=Mt[a]-Mt[i];
  BB[a][i]=II*A[a]*B[i]/(Mt[a]-Mhat[i]);
 }/*i*/
}/*a*/

/*two important auxiliary vectors */

cl_N * gvmin;
gvmin=new cl_N[4];

for(b0=0;b0<=3;b0++){
  
gvmin[b0]=expt(g,1-Mt[b0]+nmin); }

cl_N * gvmax;
gvmax=new cl_N[4];

for(b0=0;b0<=3;b0++){
gvmax[b0]=expt(g,1-Mt[b0]+nmax); }



/* Important auxiliary array*/
long int Nas[4][2]; /*Na*/

for(a=0;a<=3;a++){ Nas[a][0]=(-1)*Mtint[a];  Nas[a][1]=Mtint[a]-1; }/*a*/
  
 

/* correcting lc, if it is too small for Fourier-transformation */
long int lct[4][4];

for(a=0;a<=3;a++){ 
   lct[a][0]=2*N0+1-Nas[a][0]; 
   lct[a][1]=3+Nas[a][0];
   lct[a][2]=2*N0+1-Nas[a][1]; 
   lct[a][3]=3+Nas[a][1];
}//a

for(a=0;a<=3;a++){ 
  for(j=0;j<=3;j++){ if(lc<lct[a][j]){ lc=lct[a][j]; }  }//j
}//a

/***/

/* auxiliary index arrays for FFT*/
int ncmax[4];
int ncfmax[4];

for(a=0;a<=3;a++){ ncmax[a]=abs(2*N0-Nas[a][0]); if(ncmax[a]<abs(2+Nas[a][0])){ ncmax[a]=abs(2+Nas[a][0]); }
                   ncfmax[a]=abs(2*N0-Nas[a][1]); if(ncfmax[a]<abs(2+Nas[a][1])){ ncfmax[a]=abs(2+Nas[a][1]);}  
}//a

int ncmin[4];
int ncfmin[4];

for(a=0;a<=3;a++){ ncmin[a]=abs(2*N0-Nas[a][0]); if(ncmin[a]>abs(2+Nas[a][0])){ ncmin[a]=abs(2+Nas[a][0]); }
                   ncfmin[a]=abs(2*N0-Nas[a][1]); if(ncfmin[a]>abs(2+Nas[a][1])){ ncfmin[a]=abs(2+Nas[a][1]);}  
}//a



int n0min[4];
int nf0min[4];

for(a=0;a<=3;a++){ if(-2<=Nas[a][0] && Nas[a][0]<=2*N0){ n0min[a]=abs(Nas[a][0]%2);  }
                   else{ n0min[a]=ncmin[a]; } 
                   
                   if(-2<=Nas[a][1] && Nas[a][1]<=2*N0){ nf0min[a]=abs(Nas[a][1]%2);  }
                   else{ nf0min[a]=ncfmin[a]; }
}//a

int kmax[4];
int kfmax[4];

for(a=0;a<=3;a++){  kmax[a]=(ncmax[a]-n0min[a])/2;  
                   kfmax[a]=(ncfmax[a]-nf0min[a])/2; 
}//a


/* auxiliary index arrays for FFT*/

int kmin[4];
int kfmin[4];

for(a=0;a<=3;a++){  if( Nas[a][0]%2==0 ){ kmin[a]=1-n0min[a]/2; if( kmin[a]<0 ){ kmin[a]=0; };  }
                    else{ kmin[a]=(1-n0min[a])/2; if( kmin[a]<0 ){ kmin[a]=0; };  }

                    if( Nas[a][1]%2==0 ){ kfmin[a]=1-nf0min[a]/2; if( kfmin[a]<0 ){ kfmin[a]=0; };  }
                    else{ kfmin[a]=(1-nf0min[a])/2; if( kfmin[a]<0 ){ kfmin[a]=0; };  }
}//a




/* determination of sampling points */
cl_R * uA=new cl_R[lc];

for(k=0;k<lc;k++){ uA[k]=(-2)*realpart(g)*cos(realpart(Pi)*(2*k + 1)/(2*lc));
}



/*Useful auxiliary matrix for Chebyshev-T expansion */
cl_N * CT=new cl_N[lc*lc];

for(n=0;n<=lc-1;n++){ 
  for(k=0;k<=lc-1;k++){ CT[n*lc+k]=Cffunc(n,k,lc,precision);  }/*k*/
}/*n*/

/*Useful auxiliary matrix for Chebyshev-U expansion in C++ index convention*/
cl_N * CU=new cl_N[lc*lc];

for(n=0;n<=lc-1;n++){ 
  for(k=0;k<=lc-1;k++){ if(k<=lc-3){ CU[n*lc+k]=(CT[n*lc+k]-CT[n*lc+k+2])/2; }else{ CU[n*lc+k]=CT[n*lc+k]/2; }  }/*k*/
}/*n*/

/*auxiliary array for Chebyshev expansion*/  
cl_N * suA=new cl_N[lc];

for(k=0;k<=lc-1;k++){ suA[k]=sqrt(4*g*g-uA[k]*uA[k]); }/*k*/
  

  
/*cutoff parameters of the large u-expansions of Q_{a|i}[u] */

long int NQT[4];

for(i=0;i<=3;i++){ NQT[i]=NCb;
}/*i*/

long int NQ;

NQ=NQT[0];
  
for(j=1;j<=3;j++){ if(NQT[j]>NQ){ NQ=NQT[j]; } }/*j*/
  
  
 /* useful constants */

cl_R fel=cl_float(1,precision)/cl_float(2,precision);
cl_R mnegyed=cl_float(-1,precision)/cl_float(4,precision);



/* important array size integer */

long int lmax=NQ;

/* (-4)^(k) vector k=0,...,lmax+1 */

cl_N * m4k=new cl_N[lmax+2];

for(j=0;j<=lmax+1;j++){ m4k[j]=expt(cl_float(-4,precision),j);
}

/* (-1/4)^(n) vector n=0,..,lmax+1 */

cl_N * m1p4k=new cl_N[lmax+2];

for(j=0;j<=lmax+1;j++){ m1p4k[j]=expt(mnegyed,j);
}


/************************************************/

/*auxiliary array fo computing b_{a|i,n}*/

cl_N * t41=new cl_N[lmax*(2*lmax-2)];

T41func(t41,m4k,lmax,precision);


/*auxiliary array fo computing b_{a|i,n}*/

cl_N * t2=new cl_N[(lmax+1)*(lmax-1)];

T2func(t2,m1p4k,lmax,precision);


/*auxiliary array fo computing b_{a|i,n}*/

cl_N * t1=new cl_N[(lmax+1)*(lmax-1)];

T1func(t1,m1p4k,lmax,precision);




/*auxiliary array fo computing b_{a|i,n}*/

cl_N * s1=new cl_N[(NQ-1)*(NQ+1)];


S1func2(s1,m1p4k,NQ,precision);




/*auxiliary array fo computing b_{a|i,n}*/

cl_N * s31=new cl_N[(2*NQ-1)*(NQ-1)];


S31func2(s31,m4k,NQ,precision);









/* auxiliary array size integer*/
long int LN0=ldiv(N0c,2).quot;

/* Auxiliary arrays for computing the large u-expansion of P_a[u] and  P^a[u] */

cl_N * SigmaSub[4];
for(a=0;a<=3;a++) { SigmaSub[a]=new cl_N[(1+NQ)*(1+N0c)]; }
cl_N * SigmaSup[4];
for(a=0;a<=3;a++) { SigmaSup[a]=new cl_N[(1+NQ)*(1+N0c)]; }


for(a=0;a<=3;a++){
sigmasubfunc2(SigmaSub[a],twiceMt[a],N0c,NQ,g,precision);
sigmasupfunc2(SigmaSup[a],twiceMt[a],N0c,NQ,g,precision);}



/* Integer cutoffs for imaginary "pull down process" to determine Q_i[u] at the sampling points of the cut*/

unsigned long int * NIT=new unsigned long int[4];

for(i=0;i<=3;i++){ NIT[i]=NIcutoff; }//a




/* integers for certain array sizes*/
long int lmaxT[4];


for(i=0;i<=3;i++) { lmaxT[i]=NQT[i]; }/*i*/
 
/*cutoff parameter used for all i for the large u-expansions of Q_{a|i}[u] */

NQ=NQT[0];
  
for(j=1;j<=3;j++){ if(NQT[j]>NQ){ NQ=NQT[j]; } }

/* Integer cutoff used for all i values as the imaginary pull down cutoff parameter for the determination of Q_i[u] at the sampling points of the cut*/

unsigned long int NI=NIT[0];
  
for(j=1;j<=3;j++){ if(NIT[j]>NI){ NI=NIT[j]; } }


lmax=NQ;






/**********************************************************************************************/
/* Memory allocatioin and initialization of some important Delta-dependent auxiliary arrays */
/*******************************************************************************************/
cl_N * alfaais[4];  
cl_N * T5[4]; 
cl_N * T3[4]; 
cl_N * S1n[4];  
cl_N * S32[4]; 

for(i=0;i<=3;i++){
 alfaais[i]=new cl_N[4*(2*NQT[i]+2)];
 T5[i]=new cl_N[4*(lmaxT[i]+1)*(2*lmaxT[i]-2)];
 T3[i]=new cl_N[4*(lmaxT[i]+1)];
 S1n[i]=new cl_N[4*(NQT[i]+1)];
 S32[i]=new cl_N[4*(NQT[i]+1)*(2*NQT[i]-1)];  }/*i*/




/* Auxiliary arrays for P_a[u] and Ptilde_a[u] and for their upper indexed counterpart at the sampling points of the cut */
cl_N * TXm=new cl_N[(N0+1)*lc];
cl_N * TXmi=new cl_N[(N0+1)*lc];



cl_N * Xrealaux=new cl_N[lc]; 

for(k=0;k<=lc-1;k++){  Xrealaux[k]=X(uA[k]/g); }


for(n=0;n<=N0;n++){
  for(k=0;k<=lc-1;k++){ TXm[n*lc+k]=expt(Xrealaux[k],2*n); 
                      TXmi[n*lc+k]=egy/TXm[n*lc+k];  
  }/*k*/
}/*n*/

cl_N * Txmsub=new cl_N[4*lc];
cl_N * Txmsubi=new cl_N[4*lc];
cl_N * Txmsup=new cl_N[4*lc];
cl_N * Txmsupi=new cl_N[4*lc];

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){  Txmsub[k*4+a]=expt(Xrealaux[k],Mt[a]); 
                         Txmsubi[k*4+a]=egy/Txmsub[k*4+a];
			 Txmsupi[k*4+a]=expt(Xrealaux[k],Mt[a]-1);
			 Txmsup[k*4+a]=egy/Txmsupi[k*4+a];
  }/*k*/
}/*a*/


delete [] Xrealaux;




/*Unified auxiliary array for P_a */
cl_N * PaT=new cl_N[4*lc*(N0+1)];

/*Unified auxiliary array for Ptilde_a */


cl_N * PtaT=new cl_N[4*lc*(N0+1)];

/*Unified auxiliary array for P^a */

cl_N * PfaT=new cl_N[4*lc*(N0+1)];

/*Unified auxiliary array for Ptilde^a */

cl_N * PftaT=new cl_N[4*lc*(N0+1)];

/* filling them up with values */

for(a=0;a<=3;a++){
 for(k=0;k<=lc-1;k++){
   for(n=0;n<=N0;n++){  PaT[(n*lc+k)*4+a]=Txmsub[k*4+a]*TXm[n*lc+k];
                        PtaT[(n*lc+k)*4+a]=Txmsubi[k*4+a]*TXmi[n*lc+k];
			PfaT[(n*lc+k)*4+a]=Txmsup[k*4+a]*TXm[n*lc+k];
			PftaT[(n*lc+k)*4+a]=Txmsupi[k*4+a]*TXmi[n*lc+k];     
  }/*n*/
}/*k*/ 
}/*a*/

delete [] Txmsub;
delete [] Txmsubi;
delete [] Txmsup;
delete [] Txmsupi;





/* Auxiliary arrays for complex shifted values  of P_a and Ptilde_a */
cl_N * TxAnm=new cl_N[lc*NI*(N0+1)];


cl_N * XAn=new cl_N[lc*NI];

for(k=0;k<=lc-1;k++){ 
  for(n=0;n<=NI-1;n++){  XAn[n*lc+k]=XS((uA[k]+II*(n+1))/g);  }/*n*/
}/*k*/

for(m=0;m<=N0;m++){
  for(k=0;k<=lc-1;k++){ 
    for(n=0;n<=NI-1;n++){  TxAnm[(m*NI+n)*lc+k]=egy/expt(XAn[n*lc+k],2*m);  }/*n*/
  }/*k*/
}/*m*/


cl_N * TxMt=new cl_N[4*NI*lc];
cl_N * TxMft=new cl_N[4*NI*lc];

for(a=0;a<=3;a++){
  for(n=0;n<=NI-1;n++){
    for(k=0;k<=lc-1;k++){ TxMt[(n*lc+k)*4+a]=egy/expt(XAn[n*lc+k],Mt[a]); 
                          TxMft[(n*lc+k)*4+a]=expt(XAn[n*lc+k],Mt[a]-1);  
    }/*k*/
  }/*n*/
}/*a*/

delete [] XAn;




/* Unified auxiliary array to compute P_a[u_A+I*k]*/

cl_N * PujT=new cl_N[4*lc*NI*(N0+1)];
cl_N * PfujT=new cl_N[4*lc*NI*(N0+1)];

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NI-1;n++){
      for(m=0;m<=N0;m++){  PujT[((m*NI+n)*lc+k)*4+a]=TxMt[(n*lc+k)*4+a]*TxAnm[(m*NI+n)*lc+k];
	                   PfujT[((m*NI+n)*lc+k)*4+a]=TxMft[(n*lc+k)*4+a]*TxAnm[(m*NI+n)*lc+k]; 
      }/*m*/
    }/*n*/
  }/*k*/
}/*a*/

delete [] TxMt;
delete [] TxMft;
delete [] TxAnm;


/* Auxiliary arrays for the computation of the large u vales of Q_{a|i}[u]*/
  
cl_N * TuANIn[4];

for(i=0;i<=3;i++){ TuANIn[i]=new cl_N[lc*(NQT[i]+1)];  }/*i*/

for(i=0;i<=3;i++){
for(k=0;k<=lc-1;k++){ 
  for(n=0;n<=NQT[i];n++){  TuANIn[i][(n*lc+k)]=egy/expt(uA[k]+II*(NIT[i]+half),2*n);
      }/*n*/ 
  }/*k*/
}/*i*/  
  

cl_N * TuMiMaNI=new cl_N[4*4*lc];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){ 
    TuMiMaNI[(k*4+a)*4+i]=expt(uA[k]+II*(NIT[i]+half),Mhat[i]-Mt[a]);   
}/*a*/
  }/*k*/  
}/*i*/

/* For error estimation: for the Qconstructor...cm() function the previous arrays corresponding to NI-1 cutoff*/

cl_N * TuANInshift[4];
for(i=0;i<=3;i++){ TuANInshift[i]=new cl_N[lc*(NQT[i]+1)]; }/*i*/

for(i=0;i<=3;i++){
for(k=0;k<=lc-1;k++){ 
  for(n=0;n<=NQT[i];n++){  TuANInshift[i][n*lc+k]=egy/expt(uA[k]+II*(((int) NIT[i])-1+half),2*n);
      }/*n*/ 
  }/*k*/
}/*i*/

  
cl_N * TuMiMaNIshift=new cl_N[4*4*lc];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){ 
    TuMiMaNIshift[(k*4+a)*4+i]=expt(uA[k]+II*(((int) NIT[i])-1+half),Mhat[i]-Mt[a]);   
}/*a*/
  }/*k*/  
}/*i*/

  
/*Arrays for Q_i[u], Q^i[u], Qtilde_i[u] and Qtilde^i[u] at the sampling points: Qlower[n*4+i-1]=Q_i[u_n] etc.*/
cl_N * Qlower=new cl_N[4*lc]; 
cl_N * Qtlower=new cl_N[4*lc];  
cl_N * Qupper=new cl_N[4*lc];  
cl_N * Qtupper=new cl_N[4*lc];  


cl_N * scT[4];

for(i=0;i<=3;i++){

scT[i]=new cl_N[4*4*(NQT[i]+1)];


/*computation of the 4x4 matrices for the linear problems for b_{a|i,n}*/
totalscTmaker2i(i,AA,BB,alfa,NQT[i],scT[i]);

/*computation of Delta dependent auxiliary arrays*/
AlfatoT2i(i,alfa,alfaais[i],T3[i],T5[i],S1n[i],S32[i],m1p4k,NQT[i],precision);


}/*i*/

/* P_a[u_n]-(P_a[u_n] from gluing)=deltaP[n*4+a-1] , and their tilded  and upper case versions, as well*/
cl_N * deltaP=new cl_N[4*lc];
cl_N * deltaPt=new cl_N[4*lc];
cl_N * deltaPup=new cl_N[4*lc];
cl_N * deltaPtup=new cl_N[4*lc];

/* array to store c_{a,-1} computed back from gluing condition */
cl_N * Ev=new cl_N[8];
cl_R Evnorma=888;


/* Arrays for storing some numbers for running diagnostics */
cl_N * Qalfadev=new cl_N[2] ;

cl_R * Qalfadevout=new cl_R[(maxiter+1)*2];


for(i=0;i<=(maxiter+1)*2-1;i++){ Qalfadevout[i]=0; }//i

/*Variables to store the values of P_a[u], Ptilde_a[u], P^a[u] and Ptilde^a[u], computed back from gluing condition*/
cl_N * deltaPv=new cl_N[4*lc];
cl_N * deltaPtv=new cl_N[4*lc];
cl_N * deltaPupv=new cl_N[4*lc];
cl_N * deltaPtupv=new cl_N[4*lc];

/* Computing the Q_i[u], P_a[u]-(P_a[u] from gluing) and the values of P_a[u] computed back from gluing condition and their tilded and upper case versions as well*/
QconstructorUJ2icm(deltaP,deltaPt,deltaPup,deltaPtup,deltaPv,deltaPtv,deltaPupv,deltaPtupv,Qlower,Qtlower,Qupper,Qtupper,AA,BB,alfa,c,cf,SigmaSub,SigmaSup,PaT,PtaT,PfaT,PftaT,PujT,PfujT,t1,t2,t41,T3,T5,S1n,s1,s31,S32,TuANIn,TuMiMaNI,
	    m1p4k,g,N0c,NQT,NQ,NIT,NI,lc,scT,precision,TuANInshift,TuMiMaNIshift,Qalfadev);


Qalfadevout[0*2+0]=realpart(Qalfadev[0]);
Qalfadevout[0*2+1]=realpart(Qalfadev[1]);


QtoEtypeIIcm(-1L,Ev,deltaPv,deltaPtv,deltaPupv,deltaPtupv,CT,CU,Nas,suA,lc,g);

Evnorma=realpart(norma2(Ev,8));


delete [] deltaPv;
delete [] deltaPtv;
delete [] deltaPupv;
delete [] deltaPtupv;



/* Putting the inital values of the equations into a vector denoted by E */
long int dimE=dimV;

cl_N * deltac=new cl_N[4*(N0+1)];
cl_N * deltacup=new cl_N[4*(N0+1)];

cl_N * E=new cl_N[dimE];

QtoEtypeIInewton2(deltac,deltacup,deltaP,deltaPt,deltaPup,deltaPtup,CT,CU,Nas,suA,lc,N0,g,ncmax,ncfmax,kmax,kfmax,n0min,nf0min,kmin,kfmin);

UJCtoE2(Nch,deltac,deltacup,Mtint,E,N0);







/* The proper set of unknowns  and their initialization */
cl_N * V=new cl_N[dimV];
/*Initialization of the independent variables of the problem */
CtoV2(Delta,c,cf,Mtint,V,Nch,N0);


/* A new copy of the true variables*/
cl_N * V0=new cl_N[dimV];

for(j=0;j<=dimV-1;j++){ V0[j]=V[j];  }/*j*/
  
 /* Making a copy of the initial equations */ 
  cl_N * E0=new cl_N[dimE];

for(j=0;j<=dimE-1;j++){ E0[j]=E[j];  }/*j*/
  
  

/* An auxiliary vector for bookkeeping the realness or imaginariness of the variables*/

long int N1=N0+Nch[1]+Nch[2]+Nch[3];

cl_N * PhiV=new cl_N[dimV]; 

 PhiV[0]=complex(1,0); 

for(k=1;k<=N0;k++){ PhiV[k]=complex(0,1); }

for(k=N0+1;k<=N0+Nch[1];k++){ PhiV[k]=complex(1,0); }

for(k=N0+Nch[1]+1;k<=N0+Nch[1]+Nch[2];k++){ PhiV[k]=complex(0,1); }

for(k=N0+Nch[1]+Nch[2]+1;k<=N1;k++){ PhiV[k]=complex(1,0); }

for(k=N1+1;k<=N1+N0;k++){ PhiV[k]=complex(0,1); }

for(k=N1+N0+1;k<=N1+2*N0;k++){ PhiV[k]=complex(1,0); }

for(k=N1+2*N0+1;k<=N1+3*N0;k++){ PhiV[k]=complex(0,1); }  
  

  
  
/*Norm square of "equations" */
cl_R Enormainf=realpart(norma2(E0,dimE));

/*Norm square of "equations" */
cl_R Enormainfold;


/* Actual value of Delta*/
cl_R Deltainf=realpart(V0[0]);


/*Derivative matrix*/ 
cl_N * DEplus=new cl_N[dimE*dimV];


long int l;

/* iteration number  */
long int * it=new long int[1];





long int * errorcode=new long int[1];
/*Setting inital value for the errocode */
errorcode[0]=0;

/* Variable to store |Delta^{(n)}-Delta^{(n-1)}| */
cl_R DDelta=cl_float(9,precision);

/* Upper bound for the norm squared of the equations in the quit condition*/
cl_R ssfmax=expt(cl_float(10,precision),(-1)*precssf);







/* The iteration cycle*/

it[0]=0L;
while(!( Enormainf<ssfmax && DDelta<expt(cl_float(10,precision),(-1)*precDelta) ) && it[0]<maxiter){// <- quit conditions

/* Making a copy of the norm squared of the equations obtained after the previous ietration */ 
Enormainfold=Enormainf;  


/*********************************************************/
/* Here starts the compuation of the derivative matrix */
/********************************************************/
for(n=0;n<=dimV-1;n++){
  
/* shift of variables with 1 lattice constant unit*/
V[n]+=PhiV[n]*h;


/*Delta and  the c's from the shifted variables*/
VtoC2(V,&Delta,c,cf,Mtint,Nch,N0);  



/* computing the Delta dependent auxiliary arrays  */
if(n==0 || n==1) {   Ca0func(Delta,c,cf,Mt,Mhat0,Mhat,A,Af,AA,B,BB,alfa,g,precision);  

for(i=0;i<=3;i++){

for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){ 
    TuMiMaNI[(k*4+a)*4+i]=expt(uA[k]+II*(NIT[i]+half),Mhat[i]-Mt[a]);   
}/*a*/
  }/*k*/   

totalscTmaker2i(i,AA,BB,alfa,NQT[i],scT[i]);


AlfatoT2i(i,alfa,alfaais[i],T3[i],T5[i],S1n[i],S32[i],m1p4k,NQT[i],precision);


}/*i*/

}/*if -n ends/



/*computing Q_i[u] and Qtilde_i[u] from the auxiliary arrays */


QconstructorUJ2i(deltaP,deltaPt,deltaPup,deltaPtup,Qlower,Qtlower,Qupper,Qtupper,AA,BB,alfa,c,cf,SigmaSub,SigmaSup,PaT,PtaT,PfaT,PftaT,PujT,PfujT,t1,t2,t41,T3,T5,S1n,s1,s31,S32,TuANIn,TuMiMaNI,m1p4k,g,N0c,NQT,NQ,NIT,NI,lc,scT,precision);

/*computing the "shifted" equations and puting them into a matrix */
QtoEtypeIInewton2(deltac,deltacup,deltaP,deltaPt,deltaPup,deltaPtup,CT,CU,Nas,suA,lc,N0,g,ncmax,ncfmax,kmax,kfmax,n0min,nf0min,kmin,kfmin);


UJCtoDE2(Nch,deltac,deltacup,Mtint,DEplus,dimV,N0,n);

/*shifting back the variables into their initial values*/
V[n]-=PhiV[n]*h;
  
}
/* cycle ends*/


/* computation of the derivative matrix out of the "shifted equations" */

for(j=0;j<=dimE-1;j++){
  for(k=0;k<=dimV-1;k++){  DEplus[j*dimV+k]-=E0[j];  DEplus[j*dimV+k]/=(h*PhiV[k]);
     }/**/
}/*j*/



/*********************************************************/
/* Here ends the compuation of the derivative matrix */
/********************************************************/


/*the unknown vector of the linear problem of the Newton-method*/ 
cl_N * dc=new cl_N[dimV];


/* Solving the linear problem to get the change of variables in the Newton-method*/
linsolvepE(DEplus,E0,dc,dimV);

/* Updating the variables with their new values */
for(j=0;j<=dimV-1;j++){      V[j]=V0[j]-dc[j];   }/*j*/

delete [] dc;   
  


/*The computation of the new values of Delta and the c's from the updated variables*/
VtoC2(V,&Delta,c,cf,Mtint,Nch,N0);   





/*Computing Delta-dependent auxiliary quantities*/

Ca0func(Delta,c,cf,Mt,Mhat0,Mhat,A,Af,AA,B,BB,alfa,g,precision);


for(i=0;i<=3;i++){

for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){ 
    TuMiMaNI[(k*4+a)*4+i]=expt(uA[k]+II*(NIT[i]+half),Mhat[i]-Mt[a]); 
    TuMiMaNIshift[(k*4+a)*4+i]=expt(uA[k]+II*(((int) NIT[i])-1+half),Mhat[i]-Mt[a]);
}/*a*/
  }/*k*/   

totalscTmaker2i(i,AA,BB,alfa,NQT[i],scT[i]);

AlfatoT2i(i,alfa,alfaais[i],T3[i],T5[i],S1n[i],S32[i],m1p4k,NQT[i],precision);


}/*i*/


/*variable to store the values of P_a[u] computed back from gluing condition, and thier upper case counterpart*/
deltaPv=new cl_N[4*lc];
deltaPtv=new cl_N[4*lc];
deltaPupv=new cl_N[4*lc];
deltaPtupv=new cl_N[4*lc];


/* Computing the Q_i[u], Qtilde_i[u], and the values of P_a[u] computed back from gluing condition, and their  upper case counterparts */
QconstructorUJ2icm(deltaP,deltaPt,deltaPup,deltaPtup,deltaPv,deltaPtv,deltaPupv,deltaPtupv,Qlower,Qtlower,Qupper,Qtupper,AA,BB,alfa,c,cf,SigmaSub,SigmaSup,PaT,PtaT,PfaT,PftaT,PujT,PfujT,t1,t2,t41,T3,T5,S1n,s1,s31,S32,TuANIn,TuMiMaNI,
	    m1p4k,g,N0c,NQT,NQ,NIT,NI,lc,scT,precision,TuANInshift,TuMiMaNIshift,Qalfadev);


Qalfadevout[(it[0]+1)*2+0]=realpart(Qalfadev[0]);
Qalfadevout[(it[0]+1)*2+1]=realpart(Qalfadev[1]);

/* Computing the equations out of the new data */

QtoEtypeIInewton2(deltac,deltacup,deltaP,deltaPt,deltaPup,deltaPtup,CT,CU,Nas,suA,lc,N0,g,ncmax,ncfmax,kmax,kfmax,n0min,nf0min,kmin,kfmin);


UJCtoE2(Nch,deltac,deltacup,Mtint,E,N0);

/*Computing c_{a,-1} and c^a_{-1} from gluing condition*/

QtoEtypeIIcm(-1L,Ev,deltaPv,deltaPtv,deltaPupv,deltaPtupv,CT,CU,Nas,suA,lc,g);


delete [] deltaPv;
delete [] deltaPtv;
delete [] deltaPupv;
delete [] deltaPtupv;

/*Computing |Deta^{(n)}-Delta_{(n-1)}|*/
DDelta=abs(realpart(Delta)-Deltainf);


//Updating the "initial" values of the equations
for(j=0;j<=dimE-1;j++){ E0[j]=E[j]; }/*j*/
  
 //Updating the "initial" values of the variables (unknowns)  
  for(j=0;j<=dimV-1;j++){ V0[j]=V[j]; }/*j*/  
    
    //Updating the norm square of the equations
    Enormainf=realpart(norma2(E0,dimE));
  
    
    //Updating the initial value of the  Delta 
    Deltainf=realpart(Delta);
   
    /* Computing the norm squared of c_{a,-1} and c^a_{-1} obtained after the iteration */ 
    Evnorma=realpart(norma2(Ev,8));
    
    
   
it[0]+=1;}/*The "big" iteration ends*/



/*If the required convergence and accuracy couldn't be reached within maxiter number of iterations, the value of errorcode variable changes from 0 to 1  */
if(it[0]==maxiter){ errorcode[0]=1; }




/********************************************************************/
/* Wrinting out the data in a way, so that Mathematica could use it*/
/******************************************************************/



stringstream strm; 

strm<<"{"<<errorcode[0]<<","<<realpart(g)<<"+I*"<<imagpart(g)<<",";

strm<<"{{"<<nb[0]<<","<<nb[1]<<",{"<<nf[0]<<","<<nf[1]<<","<<nf[2]<<","<<nf[3]<<"},"<<na[0]<<","<<na[1]<<"},"<<multiplicity<<"},";

strm<<realpart(Delta)<<"+I*"<<imagpart(Delta)<<",{";

/* c_{1,n} */

strm<<realpart(c[0][0])<<"+I*"<<imagpart(c[0][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(c[0][n])<<"+I*"<<imagpart(c[0][n]);     }/*n*/

/*kovetkezo elem zarojelezese */
  
  strm<<"},{";
  
/* c_{2,n} */

strm<<realpart(c[1][0])<<"+I*"<<imagpart(c[1][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(c[1][n])<<"+I*"<<imagpart(c[1][n]);     }/*n*/

  
  strm<<"},{"; 
  
/* c_{3,n} */

strm<<realpart(c[2][0])<<"+I*"<<imagpart(c[2][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(c[2][n])<<"+I*"<<imagpart(c[2][n]);     }/*n*/

  
  strm<<"},{";  
  
/* c_{4,n} */

strm<<realpart(c[3][0])<<"+I*"<<imagpart(c[3][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(c[3][n])<<"+I*"<<imagpart(c[3][n]);     }/*n*/

  
  strm<<"},{";  
  

/* c^{1,n} */

strm<<realpart(cf[0][0])<<"+I*"<<imagpart(cf[0][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(cf[0][n])<<"+I*"<<imagpart(cf[0][n]);     }/*n*/
  
  strm<<"},{";
  
/* c^{2,n} */

strm<<realpart(cf[1][0])<<"+I*"<<imagpart(cf[1][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(cf[1][n])<<"+I*"<<imagpart(cf[1][n]);     }/*n*/

  
  strm<<"},{"; 
  
/* c^{3,n} */

strm<<realpart(cf[2][0])<<"+I*"<<imagpart(cf[2][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(cf[2][n])<<"+I*"<<imagpart(cf[2][n]);     }/*n*/

  
  strm<<"},{";  
  
/* c^{4,n} */

strm<<realpart(cf[3][0])<<"+I*"<<imagpart(cf[3][0]);

for(n=1;n<=N0;n++){ strm<<","<<realpart(cf[3][n])<<"+I*"<<imagpart(cf[3][n]);     }/*n*/

  
  strm<<"},{";  
  
/* parameters of the running */  


strm<<it[0]<<","<<realpart(Enormainf)<<","<<realpart(DDelta)<<","<<N0<<","<<NCb<<","<<NIcutoff<<","<<lc<<","<<DH<<","<<workingprecision<<","
<<precDelta<<","<<precssf<<","<<maxiter<<","<<Evnorma<<"},{";



/*Q shift deviations */

strm<<Qalfadevout[0*2+0];

/*Deviations of Q/Qtilde ratios from the gluing constants */

for(i=1;i<=(it[0]+1)-1;i++){ 
  strm<<","<<Qalfadevout[i*2+0]; }//i

strm<<"},{";

strm<<Qalfadevout[0*2+1];

for(i=1;i<=(it[0]+1)-1;i++){ strm<<","<<Qalfadevout[i*2+1]; }//i

strm<<"}}"; 




/* The string to display */

string textout;

textout=strm.str(); 








long int strlen=textout.length();

/* Transforming the string to a form being interpretable by Mathematica */

for(i=0;i<=strlen;i++){
 if(textout[i]=='L'){    if(i>=3 && textout[i-1]=='0' && textout[i-2]=='.' && textout[i-3]=='0' ){ textout.replace(i,1,"0"); strlen=textout.length(); }  else{ textout.replace(i,1,"*10^"); strlen=textout.length();} } 
}/*i*/





cout<<textout;





/* Deleting previously allocated memory*/

delete [] errorcode; 

delete [] Qalfadev;
delete [] Qalfadevout;
delete [] Ev;
 

for(a=0;a<=3;a++){ delete [] c[a];
                   delete [] cf[a];
		   delete [] SigmaSub[a];
		   delete [] SigmaSup[a];
}/*a*/



delete [] uA;
delete [] m4k;
delete [] m1p4k;

delete [] TXm;
delete [] TXmi;

delete [] PaT;
delete [] PtaT;
delete [] PfaT;
delete [] PftaT;

delete [] PujT;
delete [] PfujT;



for(i=0;i<=3;i++){ delete [] TuANIn[i]; delete [] TuANInshift[i]; }



delete [] NIT;


delete [] TuMiMaNI;
delete [] TuMiMaNIshift;

delete [] Qlower;
delete [] Qupper;

delete [] Qtlower;
delete [] Qtupper;

delete [] V;
delete [] V0;

delete [] E;
delete [] E0;

delete [] PhiV;

delete [] DEplus;



delete [] t41;
delete [] t2;
delete [] t1;
delete [] s1;
delete [] s31;

for(i=0;i<=3;i++){
delete [] alfaais[i];
delete [] T3[i];
delete [] T5[i];
delete [] S1n[i];
delete [] S32[i];  }

delete [] gvmin;
delete [] gvmax;


delete [] deltaP;
delete [] deltaPt;
delete [] deltaPup;
delete [] deltaPtup;

delete [] CT;
delete [] CU;
delete [] suA;

delete [] deltac;
delete [] deltacup;


return 0;
}
/******************************************************************************************/
/******************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/******************************************************************************************/
