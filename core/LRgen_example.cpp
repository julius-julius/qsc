/*Libraries to be included*/
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



/* Definition of necessary functions */

/* Solver for linear equations */


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


/* Solver for real linear equations */

void linsolvepER(cl_R * a,cl_R * b,cl_R x[],long int n)//decomposition with pivoting, then solving the equation
{
	long int d=1;
  	for(long int k=0;k<n-1;k++)
  	{
   		cl_R am = a[k*n+k]; //pivoting
	 	long int xm=k;
	 	for(long int i=k;i<n;i++)
	 	{
	   		if(abs(am)<abs(a[i*n+k]))
	  		{
	     		am=a[i*n+k];
	     		xm=i;
	   			cl_R c;
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
	cl_R y[n]; 
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
     	x[k]=0;
     	for(long int j=k+1;j<n;j++)
     	{
			x[k]-=a[k*n+j]*x[j];
     	}
     	x[k]=(y[k]+x[k])/a[k*n+k];
   }
}





/* Real Binomial[z,n] with working precision=prec */
cl_R binomial(cl_R z,long int n,float_format_t prec){
cl_R z1;
cl_R zm=-z;
long int i;
long int n2=abs(n);
cl_R n1=cl_float((cl_I) n2,prec);
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


/* Complex Binomial[z,n] with working precision=prec */
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


/* (-1)^m  */
cl_N m1m(long int m){
cl_N r;
if(m % 2 == 0){ r=complex(1,0);}
else{ r=complex(-1,0);}
return r;
}


/* An auxiliary function for Chebyshev-expansion */

cl_N Cffunc(int k, int i,int lc,float_format_t prec){
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
  ldiv_t kq=ldiv(n,2L);
  k=kq.quot;
  q0=kq.rem;

  for(s=0;s<=k-r;s++){ R+=kappabar(twiceMta,s,prec)*kappa(2*r+q0,k-r-s,prec); }
  
  R*=expt(sqrt(g),twiceMta+2*n);
  
  return R;
  
}

/* Auxiliary function for the 1/u expansion of P^a[u] */

void sigmasubfunc(cl_N * sigma/*(1+NQ)*(1+[N0/2])*/,long int twiceMta,long int N00/*N0*/,long int NQ0/*NQ*/,cl_N g,float_format_t prec){
  long int intN0p2=ldiv(N00,2L).quot;
  long int NQp1=NQ0+1;
  long int n,r;
  
  for(n=0;n<=NQ0;n++){  
    for(r=0;r<=intN0p2;r++){  sigma[r*NQp1+n]=complex(0,0);
      sigma[r*NQp1+n]+=fsigma(twiceMta,n,r,g,prec);
      
    }/*r*/
  }/*n*/
  
}


/* An auxiliary function for the 1/u expansion of P_a[u]= Sum of{ksub[a][n]/u^n } in the range n=0,...,NQ0 */

void kanfunc(cl_N * ksub[4],cl_N * c[4],cl_N * sigmasub[4],long int NQ0,long int N00){
  long int n,k,q0,r,a,intN0p2,intN0p2m1;
  intN0p2=ldiv(N00,2L).quot;
  intN0p2m1=ldiv(N00-1,2L).quot;
  ldiv_t kq;
  long int NQp1=NQ0+1;
  for(a=0;a<=3;a++){
  for(n=0;n<=NQ0;n++){  kq=ldiv(n,2L); 
     q0=kq.rem; 
     if(q0==0){ k=min(kq.quot,intN0p2); }else{ k=min(kq.quot,intN0p2m1); } 
    ksub[a][n]=complex(0,0);
for(r=0;r<=k;r++){  ksub[a][n]+=c[a][2*r+q0]*sigmasub[a][r*NQp1+n];  }/*r*/    
  }/*n*/
  }/*a*/
  }



/* /* An auxiliary function for the 1/u expansion of P_a[u]= Sum of{ksub[a][n]/u^n } in the index range n being between nmin and nmax  */

void kanfuncpart(cl_N * ksub[4],cl_N * c[4],cl_N * sigmasub[4],long int NQ0,long int N00,long int nmin,long int nmax){
  long int n,k,q0,r,a,intN0p2,intN0p2m1;
  intN0p2=ldiv(N00,2L).quot;
  intN0p2m1=ldiv(N00-1L,2L).quot;
  ldiv_t kq;
  long int NQp1=NQ0+1;
  for(a=0;a<=3;a++){
  for(n=nmin;n<=nmax;n++){  kq=ldiv(n,2L); 
     q0=kq.rem; 
     if(q0==0){ k=min(kq.quot,intN0p2); }else{ k=min(kq.quot,intN0p2m1); } 
    ksub[a][n]=complex(0,0);
for(r=0;r<=k;r++){  ksub[a][n]+=c[a][2*r+q0]*sigmasub[a][r*NQp1+n];  }/*r*/    
  }/*n*/
  }/*a*/
  }

/* The functionally the previous function, but only 1 single element: ksub[a][n] is determined */

void dkanfunc(cl_N * dksub,cl_N * c[4],cl_N * sigmasub[4],long int NQ0,long int N00,long int n){
  long int k,q0,r,a,intN0p2,intN0p2m1;
  intN0p2=ldiv(N00,2L).quot;
  intN0p2m1=ldiv(N00-1L,2L).quot;
  ldiv_t kq;
  long int NQp1=NQ0+1L;
  for(a=0;a<=3;a++){ 
    kq=ldiv(n,2L); 
     q0=kq.rem; 
     if(q0==0){ k=min(kq.quot,intN0p2); }else{ k=min(kq.quot,intN0p2m1); } 
    dksub[a]=complex(0,0);
for(r=0;r<=k-1;r++){  dksub[a]+=c[a][2*r+q0]*sigmasub[a][r*NQp1+n];  }/*r*/    
  //}/*n*/
  }/*a*/
  }  



/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void T41func(cl_N * T41,cl_N * m4k,long int lmax,float_format_t prec)
{ long int k,m,ref;
  for(m=0;m<=2*lmax-3;m++){ ref=ldiv(m,2L).quot;
    for(k=0;k<=lmax-1;k++){
     if(k<=ref){  T41[m*lmax+k]=binomial(cl_float((cl_I) -2*(k+1) ,prec),m-2*k,prec)*m4k[k+1]; }else{ T41[m*lmax+k]=complex(0,0); }
    }
  }
}



/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void T2func(cl_N * T2,cl_N * m1p4k,long int lmax,float_format_t prec){
  long int l,k,lmm1;
  lmm1=lmax-1;
    for(l=2;l<=lmax;l++){
  for(k=0;k<=lmax-2;k++){ 
    if(k<=l-2){ T2[l*lmm1+k]=cbinomial(cl_float((cl_I) -2*(k+1),prec),2*(l-k-1),prec)*m1p4k[l-k-1];  }else{ T2[l*lmm1+k]=complex(0,0); }
  }  
  }
  /* l=0 , 1 cases */
  for(l=0;l<=1;l++){
  for(k=0;k<=lmax-2;k++){ T2[l*lmm1+k]=complex(0,0); }
    }
}

/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void T1func(cl_N * T1,cl_N * m1p4k,long int lmax,float_format_t prec){
  long int l,k,lmm1;
  lmm1=lmax-1;
    for(l=2;l<=lmax;l++){
  for(k=0;k<=lmax-2;k++){ 
    if(k<=l-2){ T1[l*lmm1+k]=cbinomial(cl_float((cl_I) -2*(k+1),prec),2*(l-k)-1,prec)*m1p4k[l-k-1];  }else{ T1[l*lmm1+k]=complex(0,0); }
  }  
  }
  /* l=0 , 1 cases */
  for(l=0;l<=1;l++){
  for(k=0;k<=lmax-2;k++){ T1[l*lmm1+k]=complex(0,0); }}
}






/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void Y1func(cl_N * Y1,cl_N * m1p4k,long int Lmax,float_format_t prec){
  long int l,k,Lmp1;
  Lmp1=Lmax+1;
    for(l=1;l<=Lmax;l++){
  for(k=0;k<=Lmax;k++){ 
    if(k<=l-1){ Y1[l*Lmp1+k]=cbinomial(cl_float((cl_I) -2*k-1,prec),2*(l-k)+1,prec)*m1p4k[l-k];  }else{ Y1[l*Lmp1+k]=complex(0,0); }
  }  
  }
  /* l=0 case */
  for(k=0;k<=Lmax;k++){ Y1[0*Lmp1+k]=complex(0,0); }
}




/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void Y2func(cl_N * Y2,cl_N * m1p4k,long int Lmax,float_format_t prec){
  long int l,k,Lmp1;
  Lmp1=Lmax+1;
    for(l=1;l<=Lmax;l++){
  for(k=0;k<=Lmax;k++){ 
    if(k<=l-1){ Y2[l*Lmp1+k]=cbinomial(cl_float((cl_I) -2*k-1,prec),2*(l-k),prec)*m1p4k[l-k];  }else{ Y2[l*Lmp1+k]=complex(0,0); }
  }  
  }
  /* l=0 case */
  for(k=0;k<=Lmax;k++){ Y2[0*Lmp1+k]=complex(0,0); }
}



/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void Y31func(cl_N * Y31,cl_N * m4k,long int Lmax,float_format_t prec){
  long int m,k,Lmp1,kLm,ref;
  Lmp1=Lmax+1;
  kLm=2*Lmax;
    for(m=0;m<=kLm;m++){ if(m==0){ ref=-1L;}else{ ref=ldiv(m-1L,2L).quot; }
  for(k=0;k<=Lmax;k++){ 
    if(k<=ref){ Y31[m*Lmp1+k]=cbinomial(cl_float((cl_I) -2*k-1,prec),m-2*k-1,prec)*m4k[k];  }else{ Y31[m*Lmp1+k]=complex(0,0); }
  }  
  }
  /* l=0 case */
  for(k=0;k<=Lmax;k++){ Y31[0*Lmp1+k]=complex(0,0); }
}




/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void S1func(cl_N * S1/*(NQ-1)*(NQ+1)*/,cl_N * ip2n/*NQ+1*/,long int NQ0,float_format_t prec){
  long int n,j,NQm1;
  NQm1=NQ0-1;
  for(n=0;n<=NQ0;n++){
   for(j=0;j<=NQ0-2;j++) {
     if(j<=n-2){ S1[n*NQm1+j]=cbinomial(cl_float((cl_I) -j-1,prec),n-j-1,prec)*ip2n[n-j-1]; }else{ S1[n*NQm1+j]=complex(0,0); }
  }/*j*/
  }/*n*/
}


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void S31func(cl_N * S31/*(NQ-1)*(NQ-1)*/,cl_N * m2ij/*NQ*/,long int NQ0,float_format_t prec){
  long int k,j,NQm1;
  NQm1=NQ0-1;
  for(k=0;k<=NQ0-2;k++){
    for(j=0;j<=NQ0-2;j++){
      if(j<=k){ S31[k*NQm1+j]=cbinomial(cl_float((cl_I) -j-1,prec),k-j,prec)*m2ij[j+1];  }else{ S31[k*NQm1+j]=complex(0,0); }
    }/*j*/
  }/*k*/
  
}





/************************/


/* Auxiliary function to fix c_{a,n} from the P_a P^a=0 condition */

cl_N cfixalo(long int a1,long int a2,long int n,cl_N * c[4]/*N0+1*/,cl_N * cf[4]/*N0+1*/){
  long int m,a;
  cl_N r;
  r=complex(0,0);
  
  for(a=a1;a<=a2;a++){ r-=cf[a][n]*c[a][0]; }
  
 for(a=0;a<=3;a++){ for(m=0;m<=n-1;m++){ r-=c[a][n-m]*cf[a][m]; }/*m*/ }/*a*/
  
  return r;
  
}


/* Auxiliary function for the deremination of zero mode constraints arising at the determination of b_{a|i,n} */

void JfRaconstructor(long int i,long int ni/*ni*/,cl_N * Jf/*[1]*/,cl_N R[4],cl_N A[4],cl_N BB[4][4],cl_N dqhat[4][4],cl_N dk[4],
		cl_N dF12[4]/*F1-dF2*/,cl_N trnull[4]/*Z*/,cl_N * gvminmax){
  
  int a,b0;
  
  cl_N ZA,ABB;
  
  ZA=complex(0,0);
  ABB=complex(0,0);
  
  for(a=0;a<=3;a++){ ZA+=trnull[a]*A[a]; ABB+=m1m(a+1)*A[3-a]*BB[a][i]; }//a
  
  /* R[a]  */
  
 for(a=0;a<=3;a++){ R[a]=gvminmax[a]*(trnull[a]*ABB+m1m(a)*BB[3-a][i]*ZA); }//a 
  
 /* Jf computation */ 
  /*auxiliary numbers*/
cl_N ZdF=complex(0,0); 
cl_N Zdk=complex(0,0); 
cl_N BBdk=complex(0,0);  
cl_N Zdqhat=complex(0,0);

for(a=0;a<=3;a++){ ZdF+=trnull[a]*dF12[a]; 
                   Zdk+=trnull[a]*dk[a];
		   BBdk+=m1m(a+1)*BB[a][i]*dk[3-a];
}//a 

/***/

for(a=0;a<=3;a++){ 
for(b0=0;b0<=3;b0++){  

    Zdqhat+=trnull[a]*A[a]*m1m(b0+1)*A[3-b0]*BB[b0][i]*dqhat[a][b0];
 
}//b0 
}//a 



Jf[0]=ZdF-Zdk*ABB-ZA*BBdk-Zdqhat;



}




/* Auxiliary function for computing the left nullvector of the  linear operator entering  at the determination of b_{a|i,n} */

void nullTRvectorfunc(long int i,cl_N A[4],cl_N alfa[4][4],cl_N beta[4][4],cl_N v0[4]){

switch(i){  
  

  case 0:  {
  
/* i=1 nullvector */

   
   
   v0[0]=alfa[0][0]*alfa[0][2]*alfa[0][3]*beta[1][3]*beta[2][1]*beta[2][3]/A[0];
   v0[1]=alfa[1][0]*alfa[1][2]*alfa[1][3]*beta[0][2]*beta[0][3]*beta[2][3]/A[1];
   v0[2]=alfa[2][0]*alfa[2][2]*alfa[2][3]*beta[0][3]*beta[1][0]*beta[1][3]/A[2];
   v0[3]=alfa[3][0]*alfa[3][2]*alfa[3][3]*beta[0][1]*beta[0][2]*beta[1][2]/A[3];   

v0[0]=v0[0]/v0[3];
v0[1]=v0[1]/v0[3];
v0[2]=v0[2]/v0[3];
v0[3]=v0[3]/v0[3];
 break; }

  case 1:{

/* i=2 nullvector */

   
   v0[0]=alfa[0][1]*alfa[0][2]*alfa[0][3]*beta[1][3]*beta[2][1]*beta[2][3]/A[0];
   v0[1]=alfa[1][1]*alfa[1][2]*alfa[1][3]*beta[0][2]*beta[0][3]*beta[2][3]/A[1];
   v0[2]=alfa[2][1]*alfa[2][2]*alfa[2][3]*beta[0][3]*beta[1][0]*beta[1][3]/A[2];
   v0[3]=alfa[3][1]*alfa[3][2]*alfa[3][3]*beta[0][1]*beta[0][2]*beta[1][2]/A[3];   

v0[0]=v0[0]/v0[3];
v0[1]=v0[1]/v0[3];
v0[2]=v0[2]/v0[3];
v0[3]=v0[3]/v0[3];
 break; }

  case 2:{

/* i=3 nullvector */
   
   
   v0[0]=alfa[0][0]*alfa[0][1]*alfa[0][2]*beta[1][2]*beta[1][3]*beta[2][3]/A[0];
   v0[1]=alfa[1][0]*alfa[1][1]*alfa[1][2]*beta[0][3]*beta[2][0]*beta[2][3]/A[1];
   v0[2]=alfa[2][0]*alfa[2][1]*alfa[2][2]*beta[0][1]*beta[0][3]*beta[1][3]/A[2];
   v0[3]=alfa[3][0]*alfa[3][1]*alfa[3][2]*beta[0][2]*beta[1][0]*beta[1][2]/A[3];   

v0[0]=v0[0]/v0[3];
v0[1]=v0[1]/v0[3];
v0[2]=v0[2]/v0[3];
v0[3]=v0[3]/v0[3];
 
break;    
  }
  
  case 3:{


/* i=4 nullvector */

   
   
   v0[0]=alfa[0][0]*alfa[0][1]*alfa[0][3]*beta[1][3]*beta[2][1]*beta[2][3]/A[0];
   v0[1]=alfa[1][0]*alfa[1][1]*alfa[1][3]*beta[0][2]*beta[0][3]*beta[2][3]/A[1];
   v0[2]=alfa[2][0]*alfa[2][1]*alfa[2][3]*beta[0][3]*beta[1][0]*beta[1][3]/A[2];
   v0[3]=alfa[3][0]*alfa[3][1]*alfa[3][3]*beta[0][1]*beta[0][2]*beta[1][2]/A[3];   
   

v0[0]=v0[0]/v0[3];
v0[1]=v0[1]/v0[3];
v0[2]=v0[2]/v0[3];
v0[3]=v0[3]/v0[3];

break;
    
  }

}/* switch ends */


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

/* A function giving the value of the B_{a0}B^{a0}=... constraint at a0=a */

cl_N BBfunc(long int a,cl_N mt[4],cl_N mh[4]){
  
  cl_N II=complex(0,1);
  long int j;
  cl_N r;
  
   r=II;
  
  for(j=0;j<=3;j++){ r*=(mh[a]-mt[j]);  }
  
  for(j=0;j<=3;j++){ if(j==a){ ; }else{ r/=(mh[a]-mh[j]); }  }
  
  return r;
  
}


/*A function to compute our choice for B_1, and B_2 */

cl_N Bnikafunc(long int a,cl_N mt[4],cl_N mh[4]){
  
  cl_N II=complex(0,1);
  long int j;
  cl_N r;
  
   r=complex(1,0);
  
  for(j=a+1;j<=3;j++){ r/=(II*(mh[j]-mh[a]));  }
  
  
  
  return r;
  
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




/*Function computing the "unconstrained" variables of the problem out of  {Delta,c_{a,n}} : the single case */

void singleCtoVLR(cl_N Delta,long int nmin,long int nmax,cl_N * c[4],long int Mtint[4],cl_N * V/*4*N0-5*/,long int N00){
  long int n,j,k;
  cl_N * ctilde;
  V[0]=Delta;
      
  ctilde=new cl_N[N00-1];
  
   k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(nmin-1) || n==(nmax-1) ){ k++; }else{ ctilde[n-k]=c[0][n+1]; }
  }/*n*/
  
  for(j=1;j<=N00-2;j++){ V[j]=ctilde[j-1];  }/*j*/
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[1]-1)){ k++; }else{ ctilde[n-k]=c[1][n+1]; }
  }/*n*/
  
  for(j=N00-1;j<=2*N00-3;j++){ V[j]=ctilde[j-N00+1];  }/*j*/
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[2]-1)  ||  n==(Mtint[1]-Mtint[2]-1)  ){ k++; }else{ ctilde[n-k]=c[2][n+1]; }
  }/*n*/
  
  for(j=2*N00-2;j<=3*N00-5;j++){ V[j]=ctilde[j-2*N00+2];  }/*j*/  
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[3]-1)   ){ k++; }else{ ctilde[n-k]=c[3][n+1]; }
  }/*n*/
  
  for(j=3*N00-4;j<=4*N00-6;j++){ V[j]=ctilde[j-3*N00+4];  }/*j*/    
    
  delete [] ctilde;  
    
}

/*Function computing the "unconstrained" variables of the problem out of  {Delta,c_{a,n}} : the double case */


void doubleCtoVLR(cl_N Delta,long int nmin,long int nmax,cl_N * c[4],long int Mtint[4],cl_N * V/*4*N0-5*/,long int N00){
  long int n,j,k;
  cl_N * ctilde;
  V[0]=Delta;
     
  ctilde=new cl_N[N00-1];
  
   k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(nmin-1) ){ k++; }else{ ctilde[n-k]=c[0][n+1]; }
  }/*n*/
  
  for(j=1;j<=N00-1;j++){ V[j]=ctilde[j-1];  }/*j*/
    
    
k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[1]-1)){ k++; }else{ ctilde[n-k]=c[1][n+1]; }
  }/*n*/
  
  for(j=N00;j<=2*N00-2;j++){ V[j]=ctilde[j-N00];  }/*j*/
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[2]-1)  ||  n==(Mtint[1]-Mtint[2]-1)  ){ k++; }else{ ctilde[n-k]=c[2][n+1]; }
  }/*n*/
  
  for(j=2*N00-1;j<=3*N00-4;j++){ V[j]=ctilde[j-2*N00+1];  }/*j*/  
  
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[3]-1)  ){ k++; }else{ ctilde[n-k]=c[3][n+1]; }
  }/*n*/
  
  for(j=3*N00-3;j<=4*N00-5;j++){ V[j]=ctilde[j-3*N00+3];  }/*j*/      
    
 
  delete [] ctilde;  
    
}


/* The inverse of the previous function, namely it  computes {Delta,c_{a,n}} out of the the "unconstrained" variables of the problem : single case */


void singleVtoCLR(cl_N *ptrDelta,long int nmin,long int nmax,cl_N * c[4],long int Mtint[4],cl_N * V/*4*N0-5*/,long int N00){
  long int n,j,k;
  cl_N * ctilde;
  *ptrDelta=V[0];
     
  ctilde=new cl_N[N00-1];
  
  
  for(j=1;j<=N00-2;j++){ ctilde[j-1]=V[j];  }/*j*/
  
   k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(nmin-1) || n==(nmax-1) ){ k++; }else{ c[0][n+1]=ctilde[n-k]; }/* ha nem fizikai erintetlenul hagyjuk, majd Qconstructor legyartja ahol kell*/
  }/*n*/
  
  for(j=N00-1;j<=2*N00-3;j++){ ctilde[j-N00+1]=V[j];  }/*j*/
  
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[1]-1)){ k++; }else{ c[1][n+1]=ctilde[n-k]; }
  }/*n*/
  

  for(j=2*N00-2;j<=3*N00-5;j++){ ctilde[j-2*N00+2]=V[j];  }/*j*/ 
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[2]-1)  ||  n==(Mtint[1]-Mtint[2]-1)  ){ k++; }else{ c[2][n+1]=ctilde[n-k]; }
  }/*n*/
   

  for(j=3*N00-4;j<=4*N00-6;j++){ ctilde[j-3*N00+4]=V[j];  }/*j*/ 
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[3]-1)   ){ k++; }else{ c[3][n+1]=ctilde[n-k]; }
  }/*n*/
    
    
  delete [] ctilde;  
    
}

/* The inverse of the previous function, namely it  computes {Delta,c_{a,n}} out of the the "unconstrained" variables of the problem  : double case*/


void doubleVtoCLR(cl_N * ptrDelta,long int nmin,long int nmax,cl_N * c[4],long int Mtint[4],cl_N * V/*4*N0-5*/,long int N00){
  long int n,j,k;
  cl_N * ctilde;
  *ptrDelta=V[0];
     
  ctilde=new cl_N[N00-1];
  
  for(j=1;j<=N00-1;j++){ ctilde[j-1]=V[j];  }/*j*/
  
   k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(nmin-1) ){ k++; }else{ c[0][n+1]=ctilde[n-k]; }
  }/*n*/
  
  for(j=N00;j<=2*N00-2;j++){ ctilde[j-N00]=V[j];  }/*j*/  
    
k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[1]-1)){ k++; }else{ c[1][n+1]=ctilde[n-k]; }
  }/*n*/
  

  for(j=2*N00-1;j<=3*N00-4;j++){ ctilde[j-2*N00+1]=V[j];  }/*j*/
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[2]-1)  ||  n==(Mtint[1]-Mtint[2]-1)  ){ k++; }else{ c[2][n+1]=ctilde[n-k]; }
  }/*n*/
    

  for(j=3*N00-3;j<=4*N00-5;j++){ ctilde[j-3*N00+3]=V[j];  }/*j*/  
  
  k=0;
  for(n=0;n<=N00-1;n++){ 
    if(n==(Mtint[0]-Mtint[3]-1)  ){ k++; }else{ c[3][n+1]=ctilde[n-k]; }
  }/*n*/
             
    
  delete [] ctilde;  
    
}


/*A function to compute some important Delta dependent quantities, needed to be updated at each iteration */

void Ca0funcLR(cl_N Delta,cl_N * c[4], cl_N Mt[4],cl_N Mhat0[4],cl_N Mhat[4],cl_N A[4],cl_N AA[4][4],
	   cl_N B[4],cl_N BB[4][4],cl_N alfa[4][4],cl_N g,float_format_t prec){
  
long int a,i;  
/* Mhat */
 
Mhat[0]=Mhat0[0]+Delta/cl_float(2,prec);
Mhat[1]=Mhat0[1]+Delta/cl_float(2,prec);
Mhat[2]=Mhat0[2]-Delta/cl_float(2,prec);
Mhat[3]=Mhat0[3]-Delta/cl_float(2,prec);  

/*A_1 and A_2 from Marboe & Volin*/

for(a=0;a<=1;a++) { A[a]=VolinAfunc(a,Mt,Mhat); }

/*A_3 and A_4 from the constraint A_a A^a=...*/

A[2]=AAfunc0(1,Mt,Mhat)/A[1];
A[3]=complex(-1,0)*AAfunc0(0,Mt,Mhat)/A[0];


/* B_1 and B_2 */

for(i=0;i<=1;i++){ B[i]=Bnikafunc(i,Mt,Mhat);  }

/*B_3 and B_4 from the constraint B_i B^i=...*/

B[2]=BBfunc(1,Mt,Mhat)/B[1];

B[3]=BBfunc(3,Mt,Mhat)/B[0];



/*updating other auxiliary matrices: A_a A_i, B_{a|i},   Mhat[i]-Mtilde[a] */

for(a=0;a<=3;a++){ 
  for(i=0;i<=3;i++){  alfa[a][i]=Mhat[i]-Mt[a];
                      
  BB[a][i]=complex(0,1)*A[a]*B[i]/(Mt[a]-Mhat[i]);
  AA[a][i]=A[a]*m1m(3-i)*A[3-i];
 }/*i*/
}/*a*/
 
/* Giving c_{a,0} from A_a */ 
 
for(a=0;a<=1;a++){ c[a][0]=A[a]/expt(g,Mt[a]);
                   
}/*a*/  
 
c[2][0]=AAfunc0(1,Mt,Mhat)/A[1]/expt(g,1-Mt[1]);  

c[3][0]=complex(-1,0)*AAfunc0(0,Mt,Mhat)/A[0]/expt(g,1-Mt[0]);  
 
  
} /* Ca0func ends */






/*Function to compute the norm squared of a vector of dimension N00*/

cl_N norma2(cl_N * E/*N00*/,long int N00){
  long int n;
  cl_N r;
  
  r=complex(0,0);
  
  for(n=0;n<=N00-1;n++){ r+=E[n]*conjugate(E[n]); }/*n*/
  
  return r;
}

/*********************/



/**********************************************************************************************************/
/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
/**********************************************************************************************************/
void T5funci(cl_N * T5,cl_N * alfaais,cl_N * m1p4k,long int lmax){
  long int a,l,m,ref,lmp1;
  lmp1=lmax+1;
  
  for(a=0;a<=3;a++){
  
  for(l=2;l<=lmax;l++){ ref=2*l-3;
  for(m=0;m<=2*lmax-3;m++){ 
        if(m<=ref) {T5[((m*lmp1+l)*4+a)]=alfaais[((2*l-m-1)*4+a)]*m1p4k[l];} else{ T5[((m*lmp1+l)*4+a)]=complex(0,0); } 
    }/*m*/  
  }/* l*/
/*l=0,1 cases */  
for(l=0;l<=1;l++){ 
  for(m=0;m<=2*lmax-3;m++){ 
         T5[((m*lmp1+l)*4+a)]=complex(0,0); 
    }/*m*/  
  }/* l*/  
  
}/* a */

}



/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void T3funci(cl_N * T3,cl_N * alfaais,cl_N * m1p4k,long int lmax){
  long int a,l;
  for(a=0;a<=3;a++){
   
  for(l=0;l<=lmax;l++){ T3[(l*4+a)]=alfaais[((2*l+1)*4+a)]*m1p4k[l]; }
   
  }
} 



/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void F1psiLR(long int i,long int l,cl_N * F12l/*4*/,cl_N BB[4][4],cl_N  alfa[4][4],cl_N * b,cl_N * T1,cl_N * T2,cl_N * T3,cl_N * T41,cl_N * T5,long int NQ0,long int NQmax){
  long int a,k,m,ref;
  cl_N mI=complex(0,-1);
  cl_N T1s,T2s,T3s,T4s,T4spart;
  long int lmax=ldiv(NQmax,2L).quot;
  
  long int lmp1=ldiv(NQ0,2L).quot+1L;
  long int lmm1=lmax-1;
  
    for(a=0;a<=3;a++){ 
      
      T1s=complex(0,0);
      for(k=0;k<=l-2;k++){ T1s+=b[(2*(k+1)*4+a)]*T1[l*lmm1+k]; }/*k*/
	
T2s=complex(0,0);	
for(k=0;k<=l-2;k++){ T2s+=b[(2*(k+1)*4+a)]*T2[l*lmm1+k]; }/*k*/
T2s*=alfa[a][i];	

T3s=b[(0*4+a)]*T3[(l*4+a)];

T4s=complex(0,0);

for(m=0;m<=2*l-3;m++){ ref=ldiv(m,2L).quot;      T4spart=complex(0,0);
  for(k=0;k<=ref;k++){ T4spart+=b[(2*(k+1)*4+a)]*T41[m*lmax+k]; }/*k*/
T4s+=T5[((m*lmp1+l)*4+a)]*T4spart;
}/*m*/


F12l[a]=mI*BB[a][i]*(T1s+T2s+T3s+T4s);	

    }/*a*/
  
} 




/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void Y32funci(cl_N * Y32,cl_N * alfaais,cl_N * m1p4k,long int Lmax){
  long int a,l,m,kLm,Lmp;
  Lmp=Lmax+1;
  kLm=2*Lmax;
  
  
 for(a=0;a<=3;a++){
   for(l=0;l<=Lmax;l++){
     for(m=0;m<=kLm;m++){
       if(m<=2*l){ Y32[((m*Lmp+l)*4+a)]=alfaais[((2*l+2-m)*4+a)]*m1p4k[l];   }else{ Y32[((m*Lmp+l)*4+a)]=complex(0,0); }
    }/*m*/
   }/*l*/ 
  }/*a*/
}               





/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */
void F1ptliLR(long int i,long int l,cl_N * F12lp1,cl_N BB[4][4],cl_N  alfa[4][4],cl_N * b,cl_N * Y1,cl_N * Y2,cl_N * Y31,cl_N * Y32,long int NQ0,long int NQmax){
  long int a,k,m,ref;
  cl_N mI=complex(0,-1);
  cl_N T1s,T2s,T3s,T3spart;
  
  long int Lmax=ldiv((NQmax-1L),2L).quot;
  long int Lmp1=Lmax+1;
  
  long int Lmp1i=ldiv((NQ0-1L),2L).quot+1L;
   
    for(a=0;a<=3;a++){ 
      
      T1s=complex(0,0);
      for(k=0;k<=l-1;k++){ T1s+=b[((2*k+1)*4+a)]*Y1[l*Lmp1+k]; }/*k*/
	
T2s=complex(0,0);	
for(k=0;k<=l-1;k++){ T2s+=b[((2*k+1)*4+a)]*Y2[l*Lmp1+k]; }/*k*/
T2s*=alfa[a][i];	

T3s=complex(0,0);

for(m=0;m<=2*l;m++){ if(m==0){ ref=-1L; }else{ ref=ldiv(m-1L,2L).quot; } 
T3spart=complex(0,0);
  for(k=0;k<=ref;k++){ T3spart+=b[((2*k+1)*4+a)]*Y31[m*Lmp1+k]; }/*k*/
T3s+=Y32[((m*Lmp1i+l)*4+a)]*T3spart;
}/*m*/



F12lp1[a]=mI*BB[a][i]*(T1s+T2s+T3s);	

    }/*a*/
}     




/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void S1nfunci(cl_N * S1n,cl_N * alfaais,cl_N * ip2n,long int NQ0){
  long int a,n;
  for(a=0;a<=3;a++){
    
  for(n=0;n<=NQ0;n++){
S1n[(n*4+a)]=alfaais[(n*4+a)]*ip2n[n];    
  }/*n*/
 
  }/*a*/
}  


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void S32funci(cl_N * S32,cl_N * alfaais,cl_N * ip2n,long int NQ0){
  long int a,n,k,NQm1;
  NQm1=NQ0-1;
  for(a=0;a<=3;a++){
    //for(i=0;i<=3;i++){
      for(n=0;n<=NQ0;n++){
	for(k=0;k<=NQ0-2;k++){
	  if(k<=n-2){ S32[((n*NQm1+k)*4+a)]=alfaais[((n-k-1)*4+a)]*ip2n[n]; }else{ S32[((n*NQm1+k)*4+a)]=complex(0,0); }
	}/*k*/
      }/*n*/
    //}/*i*/
  }/*a*/
}           


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void F2miLR(long int i,long int m,cl_N * F2m,cl_N AA[4][4],cl_N BB[4][4],cl_N * b,cl_N * q,cl_N * S1n,cl_N * S1,cl_N * S31,cl_N * S32,long int NQ0,long int NQmax){
  long int a,b0,n,k,j,NQm1,NQm1i;
  NQm1=NQmax-1;
  NQm1i=NQ0-1;
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
      
      for(k=0;k<=n-2;k++){ S3spart=complex(0,0);
	for(j=0;j<=k;j++){
	  S3spart+=b[((j+1)*4+b0)]*S31[k*NQm1+j];
	}
S3s+=S32[((n*NQm1i+k)*4+b0)]*S3spart;
}/*k*/
      
Sq+=q[((m-n)*4+a)*4+b0]*(S1s+S2s+S3s);	
      
}/*n*/



F2m[a]+=AA[a][b0]*BB[b0][i]*(S0s+Sq);      
      
      
    }/*b0*/
       
      }/*a*/
   
}  


/* Auxiliary function for sources of the linear problems arising  at the determination of b_{a|i,n} */

void dF2mi(long int m,long int i,cl_N * dF2m,cl_N AA[4][4],cl_N BB[4][4],cl_N * b,cl_N * q,cl_N * S1n,cl_N * S1,cl_N * S31,cl_N * S32,long int NQ0,long int NQmax){
  long int a,b0,n,k,j,NQm1,NQm1i;
  NQm1i=NQ0-1;
  NQm1=NQmax-1;
  cl_N S0s,S1s,S2s,S3s,S3spart,Sq;
  
  
  for(a=0;a<=3;a++){  dF2m[a]=complex(0,0);
    for(b0=0;b0<=3;b0++){
      
      S0s=complex(0,0);
    if(m>=2){ for(n=1;n<=m-1;n++){S0s+=b[(n*4+b0)]*q[((m-n)*4+a)*4+b0];    }/*n*/}else{}  
    
    Sq=complex(0,0);
    
    for(n=1;n<=m;n++){
      
      S1s=S1n[(n*4+b0)];
      
      S2s=complex(0,0);
      
      for(j=0;j<=n-2;j++){
	S2s+=b[((j+1)*4+b0)]*S1[n*NQm1+j];
      }/*j*/
      
     
     S3s=complex(0,0);
      
      for(k=0;k<=n-2;k++){ S3spart=complex(0,0);
	for(j=0;j<=k;j++){
	  S3spart+=b[((j+1)*4+b0)]*S31[k*NQm1+j];
	}
S3s+=S32[((n*NQm1i+k)*4+b0)]*S3spart;
}/*k*/
      
Sq+=q[((m-n)*4+a)*4+b0]*(S1s+S2s+S3s);	
      
}/*n*/



dF2m[a]+=AA[a][b0]*BB[b0][i]*(S0s+Sq);            
      
    }/*b0*/
    
    
      }/*a*/
   
  
}  





/* Function to compute the 4x4 matrices of the problems for b_{a|i,n} */
void totalscTmakerLRi(long int i,cl_N AA[4][4],cl_N B[4][4],cl_N alfa[4][4],long int NQ0,cl_N * scT){
long int m;
long int a,b0;
cl_N II=complex(0,1);
for(m=0;m<=NQ0;m++){
for(a=0;a<=3;a++){
for(b0=0;b0<=3;b0++){  if(m==0){scT[((m)*4+a)*4+b0]=complex(0,0);}
else{  if(a==b0){scT[((m)*4+a)*4+b0]=AA[a][b0]*B[b0][i]-II*B[a][i]*(m-alfa[a][i]);  }
       else{ scT[((m)*4+a)*4+b0]=AA[a][b0]*B[b0][i]; }     }
}/*b0*/
}/*a*/
//}/*i*/
} /*m*/ 
}    




/* Function to compute important auxiliary arrays for solving the linear problems for b_{a|i,n} */
void AlfatoTLRi(long int i,cl_N alfa[4][4],cl_N * alfaais,cl_N * T3,cl_N * T5,cl_N * Y32,cl_N * S1n,cl_N * S32,cl_N * ip2n,cl_N * m1p4k,long int lmax,long int Lmax,long int NQ0,float_format_t prec){

  long int a,m;

for(a=0;a<=3;a++){
  
    for(m=0;m<=NQ0+1;m++){ alfaais[(m*4+a)]=cbinomial(alfa[a][i],m,prec); 
      
    }/*m*/
  //}/*i*/
}/*a*/


T5funci(T5,alfaais,m1p4k,lmax);


T3funci(T3,alfaais,m1p4k,lmax);


Y32funci(Y32,alfaais,m1p4k,Lmax);


S1nfunci(S1n,alfaais,ip2n,NQ0);


S32funci(S32,alfaais,ip2n,NQ0);



}/*AlfatoTLRi function ends*/




/* The function computing the Q_i[u] functions and the difference of P_a[u] and P_a[u] computed back from the gluing equations at the sampling points of the cut: single case */

void QsingleLRuji(cl_N * Qlower,cl_N * Qtlower,cl_N * deltaP,cl_N * deltaPt,cl_N A[4],cl_N AA[4][4],cl_N BB[4][4],cl_N alfa[4][4],cl_N * c[4],cl_N * sigmasub[4],cl_N * PaT,cl_N * PtaT,cl_N * PujT,cl_N * T1,cl_N * T2,cl_N * T41,cl_N * T3[4],cl_N * T5[4],cl_N * Y1,cl_N * Y2,cl_N * Y31,cl_N * Y32[4],cl_N * S1n[4],cl_N * S1,cl_N * S31,cl_N * S32[4],cl_N * TuANIn[4],cl_N * TuMiMaNI,cl_N g,long int N00,long int NQ0[2],long int NQmax,unsigned long int * NI,unsigned long int NImax,long int lc,long int ni[4],cl_N lambda0[4],long int * niminmax/*4*/,cl_N * gvmin,cl_N * gvmax,cl_N trnullmin[4],cl_N trnullmax[4],cl_N * scTm[4],float_format_t prec){
  
 
  long int i,a,m,n,k,j,b0;
  
  
  long int imin,imax,n0,n1;
  
  /* zero-mode levels and their i-indeces */
  
  n0=niminmax[0];
  n1=niminmax[1];
  imin=niminmax[2];
  imax=niminmax[3];



cl_N * ksub[4];

for(a=0;a<=3;a++){ ksub[a]=new cl_N[NQmax+1];  }

/* Computing k_{a,n} upto n<=n0-1: P_a[u]Sum k_{a,n}/u^(n+Mtilde[a])*/		   
		   
kanfuncpart(ksub,c,sigmasub,NQmax,N00,0,min(n0-1,NQmax));


cl_N * q=new cl_N[4*4*(1+NQmax)];


/* Computing q^{a,b}_n=k_{a,n} k^{b}_n upto n<=n0-1: P^a[u]Sum k^{a}_n/u^(1-n-Mtilde[a])*/

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=0;n<=min(n0-1,NQmax);n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   


/* Initialize b_{a|i}_n */
cl_N * b[4];

for(i=0;i<=3;i++) { b[i]=new cl_N[4*(NQ0[i]+1)]; }/*i*/

for(i=0;i<=3;i++){
for(m=0;m<=NQ0[i];m++){
for(a=0;a<=3;a++){
 if(m==0){b[i][(m*4+a)]=complex(1,0); } else{b[i][(m*4+a)]=complex(0,0); } 
}/*a*/
}/*m*/
}/*i*/		   
		   
cl_N x4[4];
cl_N v4[4];
cl_N M4[4][4];

cl_N x3[3];
cl_N v3[3];
cl_N M3[3][3];


cl_N * f1=new cl_N[4];
cl_N * f2=new cl_N[4];


 long int l0;
 
/* Solving b_{a|i}_n upto n0-1 */
for(i=0;i<=3;i++){
  
  for(m=1;m<=n0-1;m++){

  
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);  
  
if(m%2==0){ l0=ldiv(m,2L).quot;  F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax); }
else{ l0=ldiv(m-1,2L).quot; F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  }  

  
  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
 /* matrix of the 4x4 linear problem*/
  for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 

    /*Solving M4.x4=v4 */

linsolvepE(&M4[0][0],&v4[0],x4,4); 



/* Updating b_{a|i}_n */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
      
}/*m*/
  
  
}/*for-i ends*/



/*********************************************************/
/* Procedure for fixing c_{1,n0}  */
/*********************************************************/
cl_N R[4];

cl_N * Jf=new cl_N[1];

cl_N * f1aux=new cl_N[4];

cl_N dk[4];
cl_N dqhat[4][4];
cl_N dF12[4];

/*initialization*/

for(a=0;a<=3;a++){ dk[a]=complex(0,0); dF12[a]=complex(0,0); R[a]=complex(0,0); 
  for(b0=0;b0<=3;b0++){ dqhat[a][b0]=complex(0,0); }//b0
}//a

/* dk  */

dkanfunc(&dk[0],c,sigmasub,NQmax,N00,n0);


/* dF12 */

dF2mi(n0,imin,&dF12[0],AA,BB,b[imin],q,S1n[imin],S1,S31,S32[imin],NQ0[imin],NQmax);  

if(m%2==0){ l0=ldiv(n0,2L).quot;  F1psiLR(imin,l0,f1aux,BB,alfa,b[imin],T1,T2,T3[imin],T41,T5[imin],NQ0[imin],NQmax); }
else{ l0=ldiv(n0-1,2L).quot; F1ptliLR(imin,l0,f1aux,BB,alfa,b[imin],Y1,Y2,Y31,Y32[imin],NQ0[imin],NQmax);  }  /*if ends*/

for(a=0;a<=3;a++){ dF12[a]=f1aux[a]-dF12[a]; }//a

/* dqhat */

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){   dqhat[a][b0]=complex(0,0);
    
    for(j=1;j<=n0-1;j++){ dqhat[a][b0]+=ksub[a][j]*m1m(b0+1)*ksub[3-b0][n0-j];  }//k
    
    dqhat[a][b0]/=AA[a][b0];
                         
  }//b0
}//a

/*  Jf and R_a belongigng to n0 */

JfRaconstructor(imin,n0,Jf,R,A,BB,dqhat,dk,dF12,trnullmin,gvmin);



/* Fixing c_{1,n0}  */

c[0][n0]=Jf[0];

for(a=1;a<=3;a++){ c[0][n0]-=c[a][n0]*R[a]; }//a

c[0][n0]/=R[0];


/*********************************************************/
/* Procedure for fixing c_{1,n0}  ends */
/*********************************************************/

/*****************************************************************************************************************/
/****************************************************************************************************************/
/* Computing k_{a,n} and q_{a,b}_n from  n=n0 to n1-1  */

kanfuncpart(ksub,c,sigmasub,NQmax,N00,n0,min(n1-1,NQmax));
 

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=n0;n<=min(n1-1,NQmax);n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   


/*Solving b-recursion from n0 to n1-1*/

/******************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n0;m<=min(n1-1,NQmax);m++){


F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  /*since f1 has been computed at n0 */

if(m==n0 && i==imin){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;  
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  /*if ends*/}/*m==n0 if ends*/


  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector  */
    
if(m==ni[i]){ /*the zero-mode case */
  
/* 3x3 minor*/
for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
/* 3-dimesional source */  
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/ 
/* Solving the 3x3 minor equation M3.x3=v3 */
linsolvepE(&M3[0][0],&v3[0],x3,3);
/*embedding the result into a 4-dimensional vector */
  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
 
/* Solution for b_{a|i}_m */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}/*the zero-mode case ends*/
else{    
/* regular case*/
for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* Updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} /* else ends*/
      
}/*m*/
  
  
}/*for-i ends */
/********************************************************************************************/

/********************************************/
/* Procedure for fixing c_{1,n1} */
/******************************************/


/*Re-initialization*/

for(a=0;a<=3;a++){ dk[a]=complex(0,0); dF12[a]=complex(0,0); R[a]=complex(0,0); 
  for(b0=0;b0<=3;b0++){ dqhat[a][b0]=complex(0,0); }//b0
}//a


/* dk at n1 */

dkanfunc(&dk[0],c,sigmasub,NQmax,N00,n1);



/* dF12 at n1 */

dF2mi(n1,imax,&dF12[0],AA,BB,b[imax],q,S1n[imax],S1,S31,S32[imax],NQ0[imax],NQmax);  

if(m%2==0){ l0=ldiv(n1,2L).quot;  F1psiLR(imax,l0,f1aux,BB,alfa,b[imax],T1,T2,T3[imax],T41,T5[imax],NQ0[imax],NQmax); }
else{ l0=ldiv(n1-1,2L).quot; F1ptliLR(imax,l0,f1aux,BB,alfa,b[imax],Y1,Y2,Y31,Y32[imax],NQ0[imax],NQmax);  }  /*if ends*/

for(a=0;a<=3;a++){ dF12[a]=f1aux[a]-dF12[a]; }//a

/* dqhat at n1 */

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){   dqhat[a][b0]=complex(0,0);
    
    for(j=1;j<=n1-1;j++){ dqhat[a][b0]+=ksub[a][j]*m1m(b0+1)*ksub[3-b0][n1-j];  }//k
    
    dqhat[a][b0]/=AA[a][b0];
                         
  }//b0
}//a

/*  Jf es R_a at n1 */

JfRaconstructor(imax,n1,Jf,R,A,BB,dqhat,dk,dF12,trnullmax,gvmax);



/* Fixing c_{1,n1}  */

c[0][n1]=Jf[0];

for(a=1;a<=3;a++){ c[0][n1]-=c[a][n1]*R[a]; }//a

c[0][n1]/=R[0];


/********************************************/
/* Procedure for fixing c_{1,n1} ends*/
/******************************************/

/* Computing k_{a,n} and q_{a,b}_n from n1 to the cutoff value*/

kanfuncpart(ksub,c,sigmasub,NQmax,N00,n1,NQmax);
 

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=n1;n<=NQmax;n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   



/*Solving the b-recursion*/

/******************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n1;m<=NQ0[i];m++){

 
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  /*since f1 is known at n1 */

if(m==n1 && i==imax){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;  
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  /*if ends*/}/*m==n1 if ends*/


  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
if(m==ni[i]){ /*zero-mode case  */
  
/* 3x3 minor*/
for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
/*3-dimsional source*/  
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  
/*solving the 3x3 minor equation */
linsolvepE(&M3[0][0],&v3[0],x3,3);
/*embedding the result into a 4-dimensional vector */
  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
 
/* solution for b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}/*zero-mode case ends */
else{    /*usual case  */

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} /* else ends*/
      
}/*m*/
  
  
}/*for-i ends */
/********************************************************************************************/

delete [] q; 		

delete [] f1;
delete [] f2;

delete [] f1aux;
delete [] Jf;

/**********************************/		   
/* Process to compute Q_{a|i}[u] */
/********************************/

for(a=0;a<=3;a++){ delete [] ksub[a];   }
                   
/*** Fixing the two ambiguous constants: alfa2 and alfa4  ****/                   
                   
cl_N A12,A14;
cl_N A22,A24;
cl_N A32,A34;

A12=complex(0,0);
A22=complex(0,0);
A32=complex(0,0);
A14=complex(0,0);
A24=complex(0,0);
A34=complex(0,0);

long int n2,n4;

if(imin==2){ n2=n0; n4=n1;}else{ n2=n1; n4=n0; }

/*n2  */

for(a=0;a<=3;a++){  A32+=m1m(a+1)*BB[3-a][3]*BB[a][0];/*-1 x A3*/
  
  for(n=0;n<=n2;n++){  A12+=m1m(a+1)*BB[3-a][3]*BB[a][1]*b[3][n*4+(3-a)]*b[1][(n2-n)*4+a]; 
                     }//n
  
for(n=0;n<=n2-n4;n++){ A22+=m1m(a+1)*BB[3-a][2]*BB[a][1]*b[2][n*4+(3-a)]*b[1][(n2-n4-n)*4+a];
                      }//n  
  
  
}//a

/* n4 */

for(a=0;a<=3;a++){  A24+=m1m(a+1)*BB[3-a][2]*BB[a][1];/* -1 x A3*/
  
  for(n=0;n<=n4;n++){  A14+=m1m(a+1)*BB[3-a][3]*BB[a][1]*b[3][n*4+(3-a)]*b[1][(n4-n)*4+a]; 
                     }//n
  
for(n=0;n<=n4-n2;n++){ A34+=m1m(a+1)*BB[3-a][3]*BB[a][0]*b[3][n*4+(3-a)]*b[0][(n2-n4-n)*4+a];
                      }//n  
  
  
}//a

/* final values of the constants*/
cl_N alfa2,alfa4;

alfa2=(A14*A22-A12*A24)/(A24*A32-A22*A34);     

alfa4=(A14*A32-A12*A34)/(A22*A34-A24*A32);
                   
/**********************************************************/

cl_N * Q[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

/*Q_{a|i}[u] at the large u cutoff line */

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  Q[a][i][k]=complex(0,0);
  for(n=0;n<=NQ0[i];n++){  Q[a][i][k]+=b[i][(n*4+a)]*TuANIn[i][n*lc+k]; }/*n*/
                       Q[a][i][k]*=(BB[a][i]*TuMiMaNI[(k*4+a)*4+i]);   
}/*k*/
}/*a*/
}/*i*/


for(i=0;i<=3;i++) { delete [] b[i]; }/*i*/



/* Arrays to store  P_a(u_A+I*n)   */

cl_N * Puj[4];

for(a=0;a<=3;a++){ Puj[a]=new cl_N[NImax*lc];       }

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NImax-1;n++){ Puj[a][n*lc+k]=complex(0,0);
                          
        for(m=0;m<=N00;m++){  Puj[a][n*lc+k]+=c[a][m]*PujT[((m*NImax+n)*lc+k)*4+a];
	                     
	  	}/*m*/			  
    }/*n*/
  }/*k*/
}/*a*/		   



cl_N vector2[4];/*auxiliary vector */


/**************************************************/
/*** Computing Q_{a|i}[u+I/2] ***/
/***************************************************/


for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){ 
    for(n=(NI[i]-1);n>=0;n--){
    for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q[a][i][k]+=m1m(b0+1)*Puj[3-b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q[a][i][k]*=Puj[a][n*lc+k];
    Q[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
    }/*k*/
    
for(a=0;a<=3;a++){ delete [] Puj[a];  }/*a*/

  
/* P_a and Ptilde_a on the real cut: P_a(u_k), Ptilde_a(u_k). */

cl_N * P[4];
cl_N * Pt[4];


for(a=0;a<=3;a++){ P[a]=new cl_N[lc];
                   Pt[a]=new cl_N[lc]; 
}/*a*/
		   
for(a=0;a<=3;a++){ 
  for(k=0;k<=lc-1;k++){  P[a][k]=complex(0,0);  
                       Pt[a][k]=complex(0,0);
		       
	for(n=0;n<=N00;n++){  P[a][k]+=c[a][n]*PaT[(n*lc+k)*4+a];
	                     Pt[a][k]+=c[a][n]*PtaT[(n*lc+k)*4+a];			     
    }/*n*/	       
  }/*k*/
}/*a*/		   
		   
/* Q_i[u] on the cut*/    

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qlower[k*4+i]-=m1m(a+1)*P[3-a][k]*Q[a][i][k];  }
}}

/* Qtilde_i[u] on the cut*/     

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qtlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qtlower[k*4+i]-=m1m(a+1)*Pt[3-a][k]*Q[a][i][k];  }
}}
 


/* Final fixing of the ambiguity: computation of new Q_i, Qtilde_i and Q_{a|i} */

for(n=0;n<=lc-1;n++){  Qlower[n*4+1]=Qlower[n*4+1]+alfa2*Qlower[n*4+0]; 
                       Qlower[n*4+3]=Qlower[n*4+3]+alfa4*Qlower[n*4+2];
		       
		       Qtlower[n*4+1]=Qtlower[n*4+1]+alfa2*Qtlower[n*4+0]; 
                       Qtlower[n*4+3]=Qtlower[n*4+3]+alfa4*Qtlower[n*4+2];
		       
		     for(a=0;a<=3;a++){ Q[a][1][n]=Q[a][1][n]+alfa2*Q[a][0][n];
		                        Q[a][3][n]=Q[a][3][n]+alfa4*Q[a][2][n];
		       }//a  
		       
}//n



/* Gluing constant defined using all Q_i s*/

cl_N alfaQ=complex(0,0);

for(k=0;k<=lc-1;k++){ alfaQ+=(Qlower[k*4+0]/conjugate(Qlower[k*4+2])+Qtlower[k*4+0]/conjugate(Qtlower[k*4+2]) 
                             -Qlower[k*4+1]/conjugate(Qlower[k*4+3])-Qtlower[k*4+1]/conjugate(Qtlower[k*4+3]) );  }/*k*/

alfaQ/=(4*lc);

alfaQ=complex(realpart(alfaQ),0);


/* Computing the differences: P_a[u]-(P_a[u] computed back from gluing) and Ptilde_a[u]-(Ptilde_a[u] computed back from gluing)*/

for(a=0;a<=3;a++){
for(n=0;n<=lc-1;n++){ 
  
  deltaP[n*4+a]=Q[a][0][n]*(Qlower[n*4+3]+conjugate(Qlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qlower[n*4+2]-conjugate(Qlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qlower[n*4+1]+conjugate(Qlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qlower[n*4+0]-conjugate(Qlower[n*4+2])*alfaQ);
  
  deltaPt[n*4+a]=Q[a][0][n]*(Qtlower[n*4+3]+conjugate(Qtlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qtlower[n*4+2]-conjugate(Qtlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qtlower[n*4+1]+conjugate(Qtlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qtlower[n*4+0]-conjugate(Qtlower[n*4+2])*alfaQ);
		      
}/*n*/
}/*a*/

		   
for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q[a][i]; }/*i*/
}/*a*/		   

for(a=0;a<=3;a++){ delete [] P[a];
                   delete [] Pt[a];} /*a*/

}/* QsingleLRuji ends */

/************************************************************************************************************/


/* The function computing the Q_i[u] functions and the difference of P_a[u] and P_a[u] computed back from the gluing equations at the sampling points of the cut : the double case */

void QdoubleLRuji(cl_N * Qlower,cl_N * Qtlower,cl_N * deltaP,cl_N * deltaPt,cl_N A[4],cl_N AA[4][4],cl_N BB[4][4],cl_N alfa[4][4],cl_N * c[4],cl_N * sigmasub[4],cl_N * PaT,cl_N * PtaT,cl_N * PujT,cl_N * T1,cl_N * T2,cl_N * T41,cl_N * T3[4],cl_N * T5[4],cl_N * Y1,cl_N * Y2,cl_N * Y31,cl_N * Y32[4],cl_N * S1n[4],cl_N * S1,cl_N * S31,cl_N * S32[4],cl_N * TuANIn[4],cl_N * TuMiMaNI,cl_N g,long int N00,long int NQ0[2],long int NQmax,unsigned long int * NI,unsigned long int NImax,long int lc,long int ni[4],cl_N lambda0[4],long int * niminmax/*4*/,cl_N * gvmin,cl_N * gvmax,cl_N trnullmin[4],cl_N trnullmax[4],cl_N * scTm[4],float_format_t prec){
  

  long int i,a,m,n,k,j,b0;
  
  
  long int imin,imax,n0;
  
  /* the level of zero-mode eequtions*/
  n0=niminmax[0];
  /*and the i-indeces corresponing to it*/
  imin=niminmax[1];
  imax=niminmax[2];




cl_N * ksub[4];

for(a=0;a<=3;a++){ ksub[a]=new cl_N[NQmax+1]; } 

/* Computing k_{a,n} upto n<=n0-1: P_a[u]Sum k_{a,n}/u^(n+Mtilde[a])*/	
                   		   		   
kanfuncpart(ksub,c,sigmasub,NQmax,N00,0,min(n0-1,NQmax));
		   
/* Computing q^{a,b}_n=k_{a,n} k^{b}_n upto n<=n0-1: P^a[u]Sum k^{a}_n/u^(1-n-Mtilde[a])*/

cl_N * q=new cl_N[4*4*(1+NQmax)];


for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=0;n<=min(n0-1,NQmax);n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   

/* Initialize b_{a|i}_n */
cl_N * b[4];

for(i=0;i<=3;i++) { b[i]=new cl_N[4*(NQ0[i]+1)]; }/*i*/

for(i=0;i<=3;i++){
for(m=0;m<=NQ0[i];m++){
for(a=0;a<=3;a++){
 if(m==0){b[i][(m*4+a)]=complex(1,0); } else{b[i][(m*4+a)]=complex(0,0); }  
}/*a*/
}/*m*/
}/*i*/		   
		   
cl_N x4[4];
cl_N v4[4];
cl_N M4[4][4];

cl_N x3[3];
cl_N v3[3];
cl_N M3[3][3];


cl_N * f1=new cl_N[4];
cl_N * f2=new cl_N[4];



 long int l0;
 

 
/* Solving b_{a|i}_n upto n0-1 */
 
for(i=0;i<=3;i++){
  
  for(m=1;m<=n0-1;m++){


  
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);  
  
if(m%2==0){ l0=ldiv(m,2L).quot;  F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax); }
else{ l0=ldiv(m-1,2L).quot; F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  }  

  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
  for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* the matrix  */  

/*Solving M4.x4=v4 */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 


/* Updating b_{a|i,n} */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
      
}/*m*/
  
  
}/*i*/



/*************************************/
/************************************/
/* Procedure for fixing c_{1,n0}   */
/*************************************/
/************************************/

cl_N R[4];

cl_N * Jf=new cl_N[1];

cl_N * f1aux=new cl_N[4];

cl_N dk[4];
cl_N dqhat[4][4];
cl_N dF12[4];


/*initialization*/

for(a=0;a<=3;a++){ dk[a]=complex(0,0); dF12[a]=complex(0,0); R[a]=complex(0,0); 
  for(b0=0;b0<=3;b0++){ dqhat[a][b0]=complex(0,0); }//b0
}//a

/* dk  */

dkanfunc(&dk[0],c,sigmasub,NQmax,N00,n0);


/* dF12  */

dF2mi(n0,imin,&dF12[0],AA,BB,b[imin],q,S1n[imin],S1,S31,S32[imin],NQ0[imin],NQmax);  

if(m%2==0){ l0=ldiv(n0,2L).quot;  F1psiLR(imin,l0,f1aux,BB,alfa,b[imin],T1,T2,T3[imin],T41,T5[imin],NQ0[imin],NQmax); }
else{ l0=ldiv(n0-1,2L).quot; F1ptliLR(imin,l0,f1aux,BB,alfa,b[imin],Y1,Y2,Y31,Y32[imin],NQ0[imin],NQmax);  }  

for(a=0;a<=3;a++){ dF12[a]=f1aux[a]-dF12[a]; }//a

/* dqhat  */

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){   dqhat[a][b0]=complex(0,0);
    
    for(j=1;j<=n0-1;j++){ dqhat[a][b0]+=ksub[a][j]*m1m(b0+1)*ksub[3-b0][n0-j];  }//k
    
    dqhat[a][b0]/=AA[a][b0];
                         
  }//b0
}//a


/* Jf and R_a at n0 */

JfRaconstructor(imin,n0,Jf,R,A,BB,dqhat,dk,dF12,trnullmin,gvmin);


for(a=1;a<=3;a++){   Jf[0]-=c[a][n0]*R[a];  }//a

c[0][n0]=Jf[0]/R[0];

/*****************************************/
/* Procedure for fixing c_{1,n0}  ended */
/***************************************/



/*****************************************************************************************************************/
/****************************************************************************************************************/
/* Computing k_{a,n} and q^{a,b}_n from  n=n0 upto cutoff value  */

kanfuncpart(ksub,c,sigmasub,NQmax,N00,n0,NQmax);
 

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=n0;n<=NQmax;n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   

/*************************************************************************************/

/*Solving b-recursion at n0 */

/*************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n0;m<=n0;m++){
 
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  /*since f1 is known at n0 */

if(m==n0 && i==imin){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;  
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  /*if ends*/}/*m==n0 if ends*/

//for(i=0;i<=3;i++){  
  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
if(m==ni[i]){ /*the zero-mode case */
  
/* 3x3 minor*/
for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
/*3-dimensional source*/  
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  
/*solving minor equation */
linsolvepE(&M3[0][0],&v3[0],x3,3);
/*embedding the result into a 4-dimensional vector */
  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i]; 
/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}/*the zero-mode case ends*/
else{    /* usual case*/

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix  */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} /* else ends*/
      
}/*m*/
   
}/*for-i ends */
/********************************************************************************************/





/* Computation of alfa2, a constant to be fixed as a consequence of the zero-mode ambiguity */

cl_N A1;
cl_N A3;

A1=complex(0,0);
A3=complex(0,0);

for(a=0;a<=3;a++){  A3+=m1m(a)*BB[3-a][3]*BB[a][0];/* -1 x A3*/
  
  for(n=0;n<=n0;n++){  A1+=m1m(a+1)*BB[3-a][3]*BB[a][1]*b[3][n*4+(3-a)]*b[1][(n0-n)*4+a];   }//n
}//a

cl_N alfa2;

alfa2=A1/A3;

/*************************************************************************************/

/*Solving the recursion for b from n0+1 upto the cutoff */

/******************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n0+1;m<=NQ0[i];m++){

F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  /*since at m=n0 f1 has been computed */

if(m==n0 && i==imax){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;   
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  /*if ends*/}/*m==n0 if ends*/


  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
if(m==ni[i]){ /*zero mode case */
  
for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
 
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  

linsolvepE(&M3[0][0],&v3[0],x3,3);
/*embedding the result into a 4-dimensional vector */
  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
 
/*solution for b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}/*zero mode case ends*/
else{    /*usual case*/

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix  */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} /* else ends*/
      
}/*m*/
  
  
}/*for-i ends */
/********************************************************************************************/


delete [] q; 	

delete [] f1;
delete [] f2;

delete [] f1aux;
delete [] Jf;



for(a=0;a<=3;a++){ delete [] ksub[a];  }
                   
/**********************************/		   
/* Process to compute Q_{a|i}[u] */
/********************************/              

cl_N * Q[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

/*Q_{a|i}[u] at the large u cutoff line */
for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  Q[a][i][k]=complex(0,0);
  for(n=0;n<=NQ0[i];n++){  Q[a][i][k]+=b[i][(n*4+a)]*TuANIn[i][n*lc+k]; }/*n*/
                       Q[a][i][k]*=(BB[a][i]*TuMiMaNI[(k*4+a)*4+i]);   
}/*k*/
}/*a*/
}/*i*/


for(i=0;i<=3;i++) { delete [] b[i]; }/*i*/

/* Arrays to store  P_a(u_A+I*n)   */

cl_N * Puj[4];


for(a=0;a<=3;a++){ Puj[a]=new cl_N[NImax*lc];
                   
  
}

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NImax-1;n++){ Puj[a][n*lc+k]=complex(0,0);
                          
        for(m=0;m<=N00;m++){  Puj[a][n*lc+k]+=c[a][m]*PujT[((m*NImax+n)*lc+k)*4+a];
	                     
	  	}/*m*/			  
    }/*n*/
  }/*k*/
}/*a*/		   


cl_N vector2[4];/*auxiliary vector*/


/**************************************************/
/*** Computing Q_{a|i}[u+I/2] ***/
/***************************************************/

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){ 
    for(n=(NI[i]-1);n>=0;n--){
    for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q[a][i][k]+=m1m(b0+1)*Puj[3-b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q[a][i][k]*=Puj[a][n*lc+k];
    Q[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
    }/*for-k ends*/
    
for(a=0;a<=3;a++){ delete [] Puj[a];                   
}/*a*/
/********************************************/
/*********************************************/
/**********************************************/

/* P_a and Ptilde_a on the real cut: P_a(u_k), Ptilde_a(u_k). */

cl_N * P[4];
cl_N * Pt[4];


for(a=0;a<=3;a++){ P[a]=new cl_N[lc];
                   Pt[a]=new cl_N[lc];		    
}/*a*/
		   
for(a=0;a<=3;a++){ 
  for(k=0;k<=lc-1;k++){  P[a][k]=complex(0,0);  
                       Pt[a][k]=complex(0,0);
		       
	for(n=0;n<=N00;n++){  P[a][k]+=c[a][n]*PaT[(n*lc+k)*4+a];
	                     Pt[a][k]+=c[a][n]*PtaT[(n*lc+k)*4+a];
			     
                        }/*n*/	       
  }/*k*/
}/*a*/		   
		   
/* Q_i[u] on the cut*/      

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qlower[k*4+i]-=m1m(a+1)*P[3-a][k]*Q[a][i][k];  }
}}

/* Qtilde_i[u] on the cut*/     

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qtlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qtlower[k*4+i]-=m1m(a+1)*Pt[3-a][k]*Q[a][i][k];  }
}}



/*Computing the glunig constant */

cl_N alfaQ=complex(0,0);

for(k=0;k<=lc-1;k++){ alfaQ+=(Qlower[k*4+0]/conjugate(Qlower[k*4+2])+Qtlower[k*4+0]/conjugate(Qtlower[k*4+2]));  }/*k*/

alfaQ/=(2*lc);

alfaQ=complex(realpart(alfaQ),0);



/* Final fixing of the ambiguity: computation of new Q_i, Qtilde_i and Q_{a|i} */

for(n=0;n<=lc-1;n++){  Qlower[n*4+1]=Qlower[n*4+1]+alfa2*Qlower[n*4+0];                        
		       
		       Qtlower[n*4+1]=Qtlower[n*4+1]+alfa2*Qtlower[n*4+0]; 
                       		       
		     for(a=0;a<=3;a++){ Q[a][1][n]=Q[a][1][n]+alfa2*Q[a][0][n];
		                        
		       }//a  
		       
}//n


/* Computing the differences: P_a[u]-(P_a[u] computed back from gluing) and Ptilde_a[u]-(Ptilde_a[u] computed back from gluing)*/

for(a=0;a<=3;a++){
for(n=0;n<=lc-1;n++){ 
  
  deltaP[n*4+a]=Q[a][0][n]*(Qlower[n*4+3]+conjugate(Qlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qlower[n*4+2]-conjugate(Qlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qlower[n*4+1]+conjugate(Qlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qlower[n*4+0]-conjugate(Qlower[n*4+2])*alfaQ);
  
  deltaPt[n*4+a]=Q[a][0][n]*(Qtlower[n*4+3]+conjugate(Qtlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qtlower[n*4+2]-conjugate(Qtlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qtlower[n*4+1]+conjugate(Qtlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qtlower[n*4+0]-conjugate(Qtlower[n*4+2])*alfaQ);
  		      
}/*n*/
}/*a*/


for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q[a][i]; }/*i*/
}/*a*/		   

for(a=0;a<=3;a++){ delete [] P[a];
                   delete [] Pt[a];} /*a*/

}/* QdoubleLRuji ends */

/**********************************************************************/

/*The same as the QsingleLRuji() function, but it does the computations two times with NI imaginary cutoff and with NI-1 imaginary cutoff */

void QsingleLRujicm(cl_N * Qlower,cl_N * Qtlower,cl_N * deltaP,cl_N * deltaPt,cl_N * deltaPv,cl_N * deltaPtv,cl_N A[4],cl_N AA[4][4],cl_N BB[4][4],cl_N alfa[4][4],cl_N * c[4],cl_N * sigmasub[4],cl_N * PaT,cl_N * PtaT,cl_N * PujT,cl_N * T1,cl_N * T2,cl_N * T41,cl_N * T3[4],cl_N * T5[4],cl_N * Y1,cl_N * Y2,cl_N * Y31,cl_N * Y32[4],cl_N * S1n[4],cl_N * S1,cl_N * S31,cl_N * S32[4],cl_N * TuANIn[4],cl_N * TuMiMaNI,cl_N g,long int N00,long int NQ0[2],long int NQmax,unsigned long int * NI,unsigned long int NImax,long int lc,long int ni[4],cl_N lambda0[4],long int * niminmax/*4*/,cl_N * gvmin,cl_N * gvmax,cl_N trnullmin[4],cl_N trnullmax[4],cl_N * scTm[4],float_format_t prec,cl_N * TuANInshift[4],cl_N * TuMiMaNIshift,cl_N * Qalfadev){
  
 
  long int i,a,m,n,k,j,b0;
  
  
  long int imin,imax,n0,n1;
  
  n0=niminmax[0];
  n1=niminmax[1];
  imin=niminmax[2];
  imax=niminmax[3];


cl_N * ksub[4];

for(a=0;a<=3;a++){ ksub[a]=new cl_N[NQmax+1];  }
                    		   
		   
kanfuncpart(ksub,c,sigmasub,NQmax,N00,0,min(n0-1,NQmax));
		   

cl_N * q=new cl_N[4*4*(1+NQmax)];


for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=0;n<=min(n0-1,NQmax);n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   


cl_N * b[4];

for(i=0;i<=3;i++) { b[i]=new cl_N[4*(NQ0[i]+1)]; }/*i*/

for(i=0;i<=3;i++){
for(m=0;m<=NQ0[i];m++){
for(a=0;a<=3;a++){
 if(m==0){b[i][(m*4+a)]=complex(1,0); } else{b[i][(m*4+a)]=complex(0,0); }  
}/*a*/
}/*m*/
}/*i*/		   
		   
cl_N x4[4];
cl_N v4[4];
cl_N M4[4][4];

cl_N x3[3];
cl_N v3[3];
cl_N M3[3][3];


cl_N * f1=new cl_N[4];
cl_N * f2=new cl_N[4];


 long int l0;
 
for(i=0;i<=3;i++){
  
  for(m=1;m<=n0-1;m++){
  
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);  
  
if(m%2==0){ l0=ldiv(m,2L).quot;  F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax); }
else{ l0=ldiv(m-1,2L).quot; F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  }  

  
  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ 
   
  for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix */  
      

linsolvepE(&M4[0][0],&v4[0],x4,4); 


for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
      
}/*m*/
  
}/*i*/

/* c_{1,n0} fixing */
cl_N R[4];

cl_N * Jf=new cl_N[1];

cl_N * f1aux=new cl_N[4];

cl_N dk[4];
cl_N dqhat[4][4];
cl_N dF12[4];

/*inicialization*/

for(a=0;a<=3;a++){ dk[a]=complex(0,0); dF12[a]=complex(0,0); R[a]=complex(0,0); 
  for(b0=0;b0<=3;b0++){ dqhat[a][b0]=complex(0,0); }//b0
}//a

/* dk  */

dkanfunc(&dk[0],c,sigmasub,NQmax,N00,n0);

/* dF12 */

dF2mi(n0,imin,&dF12[0],AA,BB,b[imin],q,S1n[imin],S1,S31,S32[imin],NQ0[imin],NQmax);  

if(m%2==0){ l0=ldiv(n0,2L).quot;  F1psiLR(imin,l0,f1aux,BB,alfa,b[imin],T1,T2,T3[imin],T41,T5[imin],NQ0[imin],NQmax); }
else{ l0=ldiv(n0-1,2L).quot; F1ptliLR(imin,l0,f1aux,BB,alfa,b[imin],Y1,Y2,Y31,Y32[imin],NQ0[imin],NQmax);  }  

for(a=0;a<=3;a++){ dF12[a]=f1aux[a]-dF12[a]; }//a

/* dqhat  */

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){   dqhat[a][b0]=complex(0,0);
    
    for(j=1;j<=n0-1;j++){ dqhat[a][b0]+=ksub[a][j]*m1m(b0+1)*ksub[3-b0][n0-j];  }//k
    
    dqhat[a][b0]/=AA[a][b0];
                         
  }//b0
}//a

/* Jf and R_a at n0 */

JfRaconstructor(imin,n0,Jf,R,A,BB,dqhat,dk,dF12,trnullmin,gvmin);


/* fixing c_{1,n0}  */

c[0][n0]=Jf[0];

for(a=1;a<=3;a++){ c[0][n0]-=c[a][n0]*R[a]; }//a

c[0][n0]/=R[0];

/*****************************************************************************************************************/
/****************************************************************************************************************/


kanfuncpart(ksub,c,sigmasub,NQmax,N00,n0,min(n1-1,NQmax));
 


for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=n0;n<=min(n1-1,NQmax);n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   


/******************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n0;m<=min(n1-1,NQmax);m++){
 
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  

if(m==n0 && i==imin){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;  
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  }


  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/
    
if(m==ni[i]){ 
  

for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 

for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  

linsolvepE(&M3[0][0],&v3[0],x3,3);

  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
 
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}else{    

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix  */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 


for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} 
      
}/*m*/
  
  
}/*i */
/********************************************************************************************/

/* solving zero-mode problem at n1  */

/*Re-initialization*/

for(a=0;a<=3;a++){ dk[a]=complex(0,0); dF12[a]=complex(0,0); R[a]=complex(0,0); 
  for(b0=0;b0<=3;b0++){ dqhat[a][b0]=complex(0,0); }//b0
}//a


/* dk  */

dkanfunc(&dk[0],c,sigmasub,NQmax,N00,n1);


/* dF12  */

dF2mi(n1,imax,&dF12[0],AA,BB,b[imax],q,S1n[imax],S1,S31,S32[imax],NQ0[imax],NQmax);  

if(m%2==0){ l0=ldiv(n1,2L).quot;  F1psiLR(imax,l0,f1aux,BB,alfa,b[imax],T1,T2,T3[imax],T41,T5[imax],NQ0[imax],NQmax); }
else{ l0=ldiv(n1-1,2L).quot; F1ptliLR(imax,l0,f1aux,BB,alfa,b[imax],Y1,Y2,Y31,Y32[imax],NQ0[imax],NQmax);  }  

for(a=0;a<=3;a++){ dF12[a]=f1aux[a]-dF12[a]; }//a

/* dqhat  */

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){   dqhat[a][b0]=complex(0,0);
    
    for(j=1;j<=n1-1;j++){ dqhat[a][b0]+=ksub[a][j]*m1m(b0+1)*ksub[3-b0][n1-j];  }//k
    
    dqhat[a][b0]/=AA[a][b0];
                         
  }//b0
}//a

/*  Jf es R_a at n1 */

JfRaconstructor(imax,n1,Jf,R,A,BB,dqhat,dk,dF12,trnullmax,gvmax);



/*  c_{1,n1} fixing */

c[0][n1]=Jf[0];

for(a=1;a<=3;a++){ c[0][n1]-=c[a][n1]*R[a]; }//a

c[0][n1]/=R[0];

/*****************************************************************************************************************/
/****************************************************************************************************************/


kanfuncpart(ksub,c,sigmasub,NQmax,N00,n1,NQmax);
 


for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=n1;n<=NQmax;n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   



/*solving b-recursion */

/******************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n1;m<=NQ0[i];m++){
 
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  

if(m==n1 && i==imax){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;  
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
} }


  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ 
    
if(m==ni[i]){ 
  

for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
  
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  

linsolvepE(&M3[0][0],&v3[0],x3,3);

  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
 
/* solution for b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}else{    

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} 
      
}/*m*/
  
  
}/*i*/
/********************************************************************************************/



delete [] q; 		

delete [] f1;
delete [] f2;

delete [] f1aux;
delete [] Jf;

		   
/*Producing Q_{a|i}[u] */		   
		   
for(a=0;a<=3;a++){ delete [] ksub[a];  }
                  
                   
/*** computing alfa2 and alfa4 **/                   
                   
cl_N A12,A14;
cl_N A22,A24;
cl_N A32,A34;

A12=complex(0,0);
A22=complex(0,0);
A32=complex(0,0);
A14=complex(0,0);
A24=complex(0,0);
A34=complex(0,0);

long int n2,n4;

if(imin==2){ n2=n0; n4=n1;}else{ n2=n1; n4=n0; }

/*n2 */

for(a=0;a<=3;a++){  A32+=m1m(a+1)*BB[3-a][3]*BB[a][0];/*-1 x A3*/
  
  for(n=0;n<=n2;n++){  A12+=m1m(a+1)*BB[3-a][3]*BB[a][1]*b[3][n*4+(3-a)]*b[1][(n2-n)*4+a]; 
                     }//n
  
for(n=0;n<=n2-n4;n++){ A22+=m1m(a+1)*BB[3-a][2]*BB[a][1]*b[2][n*4+(3-a)]*b[1][(n2-n4-n)*4+a];
                      }//n  
  
  
}//a

/******** n4 ******/

for(a=0;a<=3;a++){  A24+=m1m(a+1)*BB[3-a][2]*BB[a][1];/* -1 x A3*/
  
  for(n=0;n<=n4;n++){  A14+=m1m(a+1)*BB[3-a][3]*BB[a][1]*b[3][n*4+(3-a)]*b[1][(n4-n)*4+a]; 
                     }//n
  
for(n=0;n<=n4-n2;n++){ A34+=m1m(a+1)*BB[3-a][3]*BB[a][0]*b[3][n*4+(3-a)]*b[0][(n2-n4-n)*4+a];
                      }//n  
  
  
}//a


cl_N alfa2,alfa4;

alfa2=(A14*A22-A12*A24)/(A24*A32-A22*A34);     

alfa4=(A14*A32-A12*A34)/(A22*A34-A24*A32);
                   
/**********************************************************/

cl_N * Q[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

cl_N * Q1[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q1[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/


cl_N * Q2[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q2[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/


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


/* P_a(u_A+I*n)  */

cl_N * Puj[4];


for(a=0;a<=3;a++){ Puj[a]=new cl_N[NImax*lc];
                   
  
}

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NImax-1;n++){ Puj[a][n*lc+k]=complex(0,0);
                         
        for(m=0;m<=N00;m++){  Puj[a][n*lc+k]+=c[a][m]*PujT[((m*NImax+n)*lc+k)*4+a];	                     
	  	}/*m*/			  
    }/*n*/
  }/*k*/
}/*a*/		   


cl_N vector2[4];

/*********************************************/
/***************************************************/

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){ 
    for(n=(NI[i]-1);n>=(NI[i]-1);n--){
    for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q[a][i][k]+=m1m(b0+1)*Puj[3-b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q[a][i][k]*=Puj[a][n*lc+k];
    Q[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
    }/*k*/

int s;

for(k=0;k<lc;k++){ 


/********************************************/
  
  for(i=0;i<=3;i++){ 
    for(n=(NI[i]-2);n>=0;n--){
    for(j=0;j<=3;j++){ vector2[j]=Q2[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q2[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q2[a][i][k]+=m1m(b0+1)*Puj[3-b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q2[a][i][k]*=Puj[a][n*lc+k];
    Q2[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
  
 /***************************** *******/ 
  
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
 

 
  /*****************************************************************************/
    }/*for-k ends*/
 /****************************************************************************/
 
 for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q2[a][i]; }/*i*/
}/*a*/	

    
for(a=0;a<=3;a++){ delete [] Puj[a]; }/*a*/
/********************************************/
/*********************************************/
/**********************************************/

/*  P_a(u_k)  */

cl_N * P[4];
cl_N * Pt[4];

for(a=0;a<=3;a++){ P[a]=new cl_N[lc];
                   Pt[a]=new cl_N[lc];  }/*a*/
		   
for(a=0;a<=3;a++){ 
  for(k=0;k<=lc-1;k++){  P[a][k]=complex(0,0);  
                       Pt[a][k]=complex(0,0);
		       
	for(n=0;n<=N00;n++){  P[a][k]+=c[a][n]*PaT[(n*lc+k)*4+a];
	                     Pt[a][k]+=c[a][n]*PtaT[(n*lc+k)*4+a];			     
                        }/*n*/	       
  }/*k*/
}/*a*/		   
		  

cl_N * Qlowerm1=new cl_N[4*lc];     


for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qlower[k*4+i]=complex(0,0); Qlowerm1[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qlower[k*4+i]-=m1m(a+1)*P[3-a][k]*Q[a][i][k]; Qlowerm1[k*4+i]-=m1m(a+1)*P[3-a][k]*Q1[a][i][k]; }
}}

/* Computing the norm of the differneces of Q_i with NI and Q_i with NI-1 */

Qalfadev[0]=complex(0,0);

for(i=0;i<=3;i++){
  for(k=0;k<=lc-1;k++){ Qalfadev[0]+=(Qlower[k*4+i]-Qlowerm1[k*4+i])*conjugate(Qlower[k*4+i]-Qlowerm1[k*4+i]);   }//k
}//i

Qalfadev[0]=sqrt(Qalfadev[0]);
    

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qtlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qtlower[k*4+i]-=m1m(a+1)*Pt[3-a][k]*Q[a][i][k];  }
}}



/* Giving new Q_i and Q_{a|i} */

for(n=0;n<=lc-1;n++){  Qlower[n*4+1]=Qlower[n*4+1]+alfa2*Qlower[n*4+0]; 
                       Qlower[n*4+3]=Qlower[n*4+3]+alfa4*Qlower[n*4+2];
		       
		       Qtlower[n*4+1]=Qtlower[n*4+1]+alfa2*Qtlower[n*4+0]; 
                       Qtlower[n*4+3]=Qtlower[n*4+3]+alfa4*Qtlower[n*4+2];
		       
		     for(a=0;a<=3;a++){ Q[a][1][n]=Q[a][1][n]+alfa2*Q[a][0][n];
		                        Q[a][3][n]=Q[a][3][n]+alfa4*Q[a][2][n];
		       }//a  
		       
}//n



/* computing the gluing constant*/

cl_N alfaQ=complex(0,0);

for(k=0;k<=lc-1;k++){ alfaQ+=(Qlower[k*4+0]/conjugate(Qlower[k*4+2])+Qtlower[k*4+0]/conjugate(Qtlower[k*4+2]) 
                             -Qlower[k*4+1]/conjugate(Qlower[k*4+3])-Qtlower[k*4+1]/conjugate(Qtlower[k*4+3]) );  }/*k*/

alfaQ/=(4*lc);

alfaQ=complex(realpart(alfaQ),0);

/* Computing the differences: P_a[u]-(P_a[u] computed back from gluing) and Ptilde_a[u]-(Ptilde_a[u] computed back from gluing)*/

for(a=0;a<=3;a++){
for(n=0;n<=lc-1;n++){ 
  
  deltaP[n*4+a]=Q[a][0][n]*(Qlower[n*4+3]+conjugate(Qlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qlower[n*4+2]-conjugate(Qlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qlower[n*4+1]+conjugate(Qlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qlower[n*4+0]-conjugate(Qlower[n*4+2])*alfaQ);
  
  deltaPt[n*4+a]=Q[a][0][n]*(Qtlower[n*4+3]+conjugate(Qtlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qtlower[n*4+2]-conjugate(Qtlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qtlower[n*4+1]+conjugate(Qtlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qtlower[n*4+0]-conjugate(Qtlower[n*4+2])*alfaQ);
  
  /*only with the gluing terms  */
  
  deltaPv[n*4+a]=Q[a][0][n]*(conjugate(Qlower[n*4+1])/alfaQ)+Q[a][1][n]*(conjugate(Qlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(conjugate(Qlower[n*4+3])*alfaQ)+Q[a][3][n]*(conjugate(Qlower[n*4+2])*alfaQ);
  
  deltaPtv[n*4+a]=Q[a][0][n]*(conjugate(Qtlower[n*4+1])/alfaQ)+Q[a][1][n]*(conjugate(Qtlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(conjugate(Qtlower[n*4+3])*alfaQ)+Q[a][3][n]*(conjugate(Qtlower[n*4+2])*alfaQ);
    		      
}/*n*/
}/*a*/


Qalfadev[1]=complex(0,0);

for(k=0;k<=lc-1;k++){ Qalfadev[1]+=(alfaQ-(Qlower[k*4+0]/conjugate(Qlower[k*4+2])))*conjugate(alfaQ-(Qlower[k*4+0]/conjugate(Qlower[k*4+2]))); 
    Qalfadev[1]+=(alfaQ-(Qtlower[k*4+0]/conjugate(Qtlower[k*4+2])))*conjugate(alfaQ-(Qtlower[k*4+0]/conjugate(Qtlower[k*4+2])));
    Qalfadev[1]+=(alfaQ+(Qlower[k*4+1]/conjugate(Qlower[k*4+3])))*conjugate(alfaQ+(Qlower[k*4+1]/conjugate(Qlower[k*4+3])));
    Qalfadev[1]+=(alfaQ+(Qtlower[k*4+1]/conjugate(Qtlower[k*4+3])))*conjugate(alfaQ+(Qtlower[k*4+1]/conjugate(Qtlower[k*4+3])));
}//k

Qalfadev[1]=sqrt(Qalfadev[1]);
Qalfadev[1]/=(4*lc*alfaQ);


delete [] Qlowerm1;

for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q[a][i]; delete [] Q1[a][i];  }/*i*/
}/*a*/
		   

for(a=0;a<=3;a++){ delete [] P[a];
                   delete [] Pt[a];} /*a*/

}/* QsingleLRujicm ends */

/*************************************************************************************************************************************/


/*The same as the QdoubleLRuji() function, but it does the computations two times with NI imaginary cutoff and with NI-1 imaginary cutoff */

void QdoubleLRujicm(cl_N * Qlower,cl_N * Qtlower,cl_N * deltaP,cl_N * deltaPt,cl_N * deltaPv,cl_N * deltaPtv,cl_N A[4],cl_N AA[4][4],cl_N BB[4][4],cl_N alfa[4][4],cl_N * c[4],cl_N * sigmasub[4],cl_N * PaT,cl_N * PtaT,cl_N * PujT,cl_N * T1,cl_N * T2,cl_N * T41,cl_N * T3[4],cl_N * T5[4],cl_N * Y1,cl_N * Y2,cl_N * Y31,cl_N * Y32[4],cl_N * S1n[4],cl_N * S1,cl_N * S31,cl_N * S32[4],cl_N * TuANIn[4],cl_N * TuMiMaNI,cl_N g,long int N00,long int NQ0[2],long int NQmax,unsigned long int * NI,unsigned long int NImax,long int lc,long int ni[4],cl_N lambda0[4],long int * niminmax/*4*/,cl_N * gvmin,cl_N * gvmax,cl_N trnullmin[4],cl_N trnullmax[4],cl_N * scTm[4],float_format_t prec,cl_N * TuANInshift[4],cl_N * TuMiMaNIshift,cl_N * Qalfadev){
  
  long int i,a,m,n,k,j,b0;
  
  long int imin,imax,n0;
  
  n0=niminmax[0];
  imin=niminmax[1];
  imax=niminmax[2];


cl_N * ksub[4];

for(a=0;a<=3;a++){ ksub[a]=new cl_N[NQmax+1];     }
		   
kanfuncpart(ksub,c,sigmasub,NQmax,N00,0,min(n0-1,NQmax));


cl_N * q=new cl_N[4*4*(1+NQmax)];


for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=0;n<=min(n0-1,NQmax);n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   

/**/

/* solving for b*/
cl_N * b[4];

for(i=0;i<=3;i++) { b[i]=new cl_N[4*(NQ0[i]+1)]; }/*i*/

for(i=0;i<=3;i++){
for(m=0;m<=NQ0[i];m++){
for(a=0;a<=3;a++){
 if(m==0){b[i][(m*4+a)]=complex(1,0); } else{b[i][(m*4+a)]=complex(0,0); } 
}/*a*/
}/*m*/
}/*i*/		   
		   
cl_N x4[4];
cl_N v4[4];
cl_N M4[4][4];

cl_N x3[3];
cl_N v3[3];
cl_N M3[3][3];


cl_N * f1=new cl_N[4];
cl_N * f2=new cl_N[4];


 long int l0;
 

for(i=0;i<=3;i++){
  
  for(m=1;m<=n0-1;m++){

  
F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);  
  
if(m%2==0){ l0=ldiv(m,2L).quot;  F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax); }
else{ l0=ldiv(m-1,2L).quot; F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  }  
 
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
  for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix */  
    

linsolvepE(&M4[0][0],&v4[0],x4,4); 


/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
      
}/*m*/
  
  
}/*for-i ends */


/*  c_{1,n0} fixing process */
cl_N R[4];

cl_N * Jf=new cl_N[1];

cl_N * f1aux=new cl_N[4];

cl_N dk[4];
cl_N dqhat[4][4];
cl_N dF12[4];


/*initialization */

for(a=0;a<=3;a++){ dk[a]=complex(0,0); dF12[a]=complex(0,0); R[a]=complex(0,0); 
  for(b0=0;b0<=3;b0++){ dqhat[a][b0]=complex(0,0); }//b0
}//a

/* dk */

dkanfunc(&dk[0],c,sigmasub,NQmax,N00,n0);

/* dF12  */

dF2mi(n0,imin,&dF12[0],AA,BB,b[imin],q,S1n[imin],S1,S31,S32[imin],NQ0[imin],NQmax);  

if(m%2==0){ l0=ldiv(n0,2L).quot;  F1psiLR(imin,l0,f1aux,BB,alfa,b[imin],T1,T2,T3[imin],T41,T5[imin],NQ0[imin],NQmax); }
else{ l0=ldiv(n0-1,2L).quot; F1ptliLR(imin,l0,f1aux,BB,alfa,b[imin],Y1,Y2,Y31,Y32[imin],NQ0[imin],NQmax);  }  /*if ends*/

for(a=0;a<=3;a++){ dF12[a]=f1aux[a]-dF12[a]; }//a

/* dqhat  */

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){   dqhat[a][b0]=complex(0,0);
    
    for(j=1;j<=n0-1;j++){ dqhat[a][b0]+=ksub[a][j]*m1m(b0+1)*ksub[3-b0][n0-j];  }//k
    
    dqhat[a][b0]/=AA[a][b0];
                         
  }//b0
}//a


/* Jf and R_a at n0  */

JfRaconstructor(imin,n0,Jf,R,A,BB,dqhat,dk,dF12,trnullmin,gvmin);


for(a=1;a<=3;a++){   Jf[0]-=c[a][n0]*R[a];  }//a

/* fixing c_{1,n0} */

c[0][n0]=Jf[0]/R[0];

/****************************************************************************************************************/


kanfuncpart(ksub,c,sigmasub,NQmax,N00,n0,NQmax);
 

for(a=0;a<=3;a++){
  for(b0=0;b0<=3;b0++){ 
    for(n=n0;n<=NQmax;n++){  q[(n*4+a)*4+b0]=complex(0,0);
    
    for(m=0;m<=n;m++){  q[(n*4+a)*4+b0]+=ksub[a][m]*m1m(b0+1)*ksub[3-b0][n-m];  }/*m*/ 
      
      q[(n*4+a)*4+b0]/=AA[a][b0];
      
    }/*n*/
  }/*b0*/
}/*a*/   

/****************************************************************************************/

/* solving the remaining recursion for b from n0 to ...*/

/******************************************************************************************/
/*n0 case*/
for(i=0;i<=3;i++){ 
for(m=n0;m<=n0;m++){

F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);
  

if(m==n0 && i==imin){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;   
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  }

  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source */
    
if(m==ni[i]){ /*zero-mode case */
  
/*minor*/
for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
/*minor source*/  
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  
/*solving minor equation */
linsolvepE(&M3[0][0],&v3[0],x3,3);
/*embedding the resul into 4-dimensional vector */
  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
/* solution for b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}else{    

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix  */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* updating b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} /* else ends*/
      
}/*m*/
  
  
}/*for-i ends */
/********************************************************************************************/

/* computing  alfa2 for fixing the ambiguity  */

cl_N A1;
cl_N A3;

A1=complex(0,0);
A3=complex(0,0);

for(a=0;a<=3;a++){  A3+=m1m(a)*BB[3-a][3]*BB[a][0];/* -1 x A3*/
  
  for(n=0;n<=n0;n++){  A1+=m1m(a+1)*BB[3-a][3]*BB[a][1]*b[3][n*4+(3-a)]*b[1][(n0-n)*4+a];   }//n
}//a

cl_N alfa2;

alfa2=A1/A3;


/*solving b-recusrion for the remaining indexes*/

/******************************************************************************************/
for(i=0;i<=3;i++){ 
for(m=n0+1;m<=NQ0[i];m++){

F2miLR(i,m,f2,AA,BB,b[i],q,S1n[i],S1,S31,S32[i],NQ0[i],NQmax);

if(m==n0 && i==imax){ for(a=0;a<=3;a++){ f1[a]=f1aux[a]; } }else{  
  if(m%2==0){ l0=ldiv(m,2L).quot;  
                                   F1psiLR(i,l0,f1,BB,alfa,b[i],T1,T2,T3[i],T41,T5[i],NQ0[i],NQmax);   
  }
else{ l0=ldiv(m-1L,2L).quot; 
                             F1ptliLR(i,l0,f1,BB,alfa,b[i],Y1,Y2,Y31,Y32[i],NQ0[i],NQmax);  
}  }


  
  for(a=0;a<=3;a++){ v4[a]=f1[a]-f2[a];}/*a*/ /* source vector */
    
if(m==ni[i]){ /*zero mode case */
  
/* 3x3 minor*/
for(a=0;a<=2;a++){for(b0=0;b0<=2;b0++){  M3[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ 
/*minor source*/  
for(a=0;a<=2;a++){ v3[a]=v4[a]-scTm[i][((m)*4+a)*4+3]*lambda0[i];    }/*a*/  
/*solving minor problem */
linsolvepE(&M3[0][0],&v3[0],x3,3);
/*embedding the result into a 4-dimensional vector */
  for(a=0;a<=2;a++){ x4[a]=x3[a];    }/*a*/
    x4[3]=lambda0[i];
 
/* solution for b */
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
}else{    

for(a=0;a<=3;a++){for(b0=0;b0<=3;b0++){  M4[a][b0]=scTm[i][((m)*4+a)*4+b0];   }/*b0*/ }/*a*/ /* matrix */
    
linsolvepE(&M4[0][0],&v4[0],x4,4); 

/* updating b*/
for(a=0;a<=3;a++){ b[i][(m*4+a)]=x4[a];}
} /* else ends*/
      
}/*m*/
 
}/*for-i ends */
/********************************************************************************************/


delete [] q; 		

delete [] f1;
delete [] f2;

delete [] f1aux;
delete [] Jf;



		   


/* Computing Q_{a|i}[u] */

for(a=0;a<=3;a++){ delete [] ksub[a];  }
                               
cl_N * Q[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/

cl_N * Q1[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q1[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/


cl_N * Q2[4][4];

for(i=0;i<=3;i++){
for(a=0;a<=3;a++){  Q2[a][i]=new cl_N[lc];
 }/*a*/
}/*i*/


/* Q_{a|i}[u+I*(NI+1/2)]*/
for(i=0;i<=3;i++){
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){  Q[a][i][k]=complex(0,0);
  for(n=0;n<=NQ0[i];n++){  Q[a][i][k]+=b[i][(n*4+a)]*TuANIn[i][n*lc+k]; }/*n*/
                       Q[a][i][k]*=(BB[a][i]*TuMiMaNI[(k*4+a)*4+i]);   
}/*k*/
}/*a*/
}/*i*/

/* Q_{a|i}[u+I*(NI+1/2-1)]*/

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
                         

/*  P_a(u_A+I*n)  and Ptilde_a(u_A+I*n)  */

cl_N * Puj[4];
//cl_N * Pfuj[4];

for(a=0;a<=3;a++){ Puj[a]=new cl_N[NImax*lc];           }

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NImax-1;n++){ Puj[a][n*lc+k]=complex(0,0);
                          
        for(m=0;m<=N00;m++){  Puj[a][n*lc+k]+=c[a][m]*PujT[((m*NImax+n)*lc+k)*4+a];
	                     
	  	}/*m*/			  
    }/*n*/
  }/*k*/
}/*a*/		   


cl_N vector2[4];/*auxiliary vector*/

/*********************************************/
/***************************************************/

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){ 
    for(n=(NI[i]-1);n>=(NI[i]-1);n--){
    for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q[a][i][k]+=m1m(b0+1)*Puj[3-b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q[a][i][k]*=Puj[a][n*lc+k];
    Q[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
    }/*for-k ends*/


int s;


for(k=0;k<lc;k++){ 


/********************************************/
  
  for(i=0;i<=3;i++){ 
    for(n=(NI[i]-2);n>=0;n--){
    for(j=0;j<=3;j++){ vector2[j]=Q2[j][i][k]; } 
  for(a=0;a<=3;a++){
      Q2[a][i][k]=complex(0,0);
      for(b0=0;b0<=3;b0++){ Q2[a][i][k]+=m1m(b0+1)*Puj[3-b0][n*lc+k]*vector2[b0];  }/*b0*/  
	Q2[a][i][k]*=Puj[a][n*lc+k];
    Q2[a][i][k]+=vector2[a];  
      }/*a*/
    }/*n*/
  }/*i*/
  
 /*** Q_{a|i}[u+I*1/2] with NI cutoff **/ 
  
 for(i=0;i<=3;i++){  for(j=0;j<=3;j++){ vector2[j]=Q[j][i][k]; }
  for(a=0;a<=3;a++){ 
    
    Q[a][i][k]=complex(0,0);
  
  for(s=0;s<=3;s++){  Q[a][i][k]+=Q2[a][s][k]*vector2[s];  }//  
  }//i  
}//a


/*** Q_{a|i}[u+I*1/2] with NI-1 cutoff **/ 

for(i=0;i<=3;i++){  for(j=0;j<=3;j++){ vector2[j]=Q1[j][i][k]; }
  for(a=0;a<=3;a++){ 
    
    Q1[a][i][k]=complex(0,0);
  
  for(s=0;s<=3;s++){  Q1[a][i][k]+=Q2[a][s][k]*vector2[s];  }//  
  }//i  
}//a

  /*****************************************************************************/
    }/*for-k ends*/
 /****************************************************************************/
 
 for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q2[a][i]; }/*i*/
}/*a*/	


    
for(a=0;a<=3;a++){ delete [] Puj[a];  }/*a*/

/**********************************************/

/* P_a(u_k) and Ptilde_a(u_k) */

cl_N * P[4];
cl_N * Pt[4];

for(a=0;a<=3;a++){ P[a]=new cl_N[lc];
                   Pt[a]=new cl_N[lc];		  
}/*a*/
		   
for(a=0;a<=3;a++){ 
  for(k=0;k<=lc-1;k++){  P[a][k]=complex(0,0);  
                       Pt[a][k]=complex(0,0);
		      
	for(n=0;n<=N00;n++){  P[a][k]+=c[a][n]*PaT[(n*lc+k)*4+a];
	                     Pt[a][k]+=c[a][n]*PtaT[(n*lc+k)*4+a];
			     
                        }/*n*/	       
  }/*k*/
}/*a*/		   

    

/* Q_i[u] computed with NI-1*/
cl_N * Qlowerm1=new cl_N[4*lc];     

/* Computing Q_i[u] at the cut*/
for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qlower[k*4+i]=complex(0,0); Qlowerm1[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qlower[k*4+i]-=m1m(a+1)*P[3-a][k]*Q[a][i][k]; Qlowerm1[k*4+i]-=m1m(a+1)*P[3-a][k]*Q1[a][i][k]; }
}}

/* Computing the norm of the differences of Q_i[u] with NI and Q_i[u] with NI-1 */

Qalfadev[0]=complex(0,0);

for(i=0;i<=3;i++){
  for(k=0;k<=lc-1;k++){ Qalfadev[0]+=(Qlower[k*4+i]-Qlowerm1[k*4+i])*conjugate(Qlower[k*4+i]-Qlowerm1[k*4+i]);   }//k
}//i

Qalfadev[0]=sqrt(Qalfadev[0]);

/* Computing Qtilde_i[u] at the cut*/

for(k=0;k<lc;k++){
  for(i=0;i<=3;i++){  Qtlower[k*4+i]=complex(0,0);
for(a=0;a<=3;a++){ Qtlower[k*4+i]-=m1m(a+1)*Pt[3-a][k]*Q[a][i][k];  }
}}


/*Computing the gluing constant */

cl_N alfaQ=complex(0,0);

for(k=0;k<=lc-1;k++){ alfaQ+=(Qlower[k*4+0]/conjugate(Qlower[k*4+2])+Qtlower[k*4+0]/conjugate(Qtlower[k*4+2]));  }/*k*/

alfaQ/=(2*lc);

alfaQ=complex(realpart(alfaQ),0);


/* Computing the new Q_i[u] and Q_{a|i}[u+I/2]*/

for(n=0;n<=lc-1;n++){  Qlower[n*4+1]=Qlower[n*4+1]+alfa2*Qlower[n*4+0]; 
                       	
		       Qtlower[n*4+1]=Qtlower[n*4+1]+alfa2*Qtlower[n*4+0]; 
                      	       
		     for(a=0;a<=3;a++){ Q[a][1][n]=Q[a][1][n]+alfa2*Q[a][0][n];
		                        
		       }//a  
		       
}//n


/* Computing the differences: P_a[u]-(P_a[u] computed back from gluing) and Ptilde_a[u]-(Ptilde_a[u] computed back from gluing)*/

for(a=0;a<=3;a++){
for(n=0;n<=lc-1;n++){ 
  
  deltaP[n*4+a]=Q[a][0][n]*(Qlower[n*4+3]+conjugate(Qlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qlower[n*4+2]-conjugate(Qlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qlower[n*4+1]+conjugate(Qlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qlower[n*4+0]-conjugate(Qlower[n*4+2])*alfaQ);
  
  deltaPt[n*4+a]=Q[a][0][n]*(Qtlower[n*4+3]+conjugate(Qtlower[n*4+1])/alfaQ)-Q[a][1][n]*(Qtlower[n*4+2]-conjugate(Qtlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(Qtlower[n*4+1]+conjugate(Qtlower[n*4+3])*alfaQ)-Q[a][3][n]*(Qtlower[n*4+0]-conjugate(Qtlower[n*4+2])*alfaQ);
  
  
    /*only with the gluing terms */
  
  deltaPv[n*4+a]=Q[a][0][n]*(conjugate(Qlower[n*4+1])/alfaQ)+Q[a][1][n]*(conjugate(Qlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(conjugate(Qlower[n*4+3])*alfaQ)+Q[a][3][n]*(conjugate(Qlower[n*4+2])*alfaQ);
  
  deltaPtv[n*4+a]=Q[a][0][n]*(conjugate(Qtlower[n*4+1])/alfaQ)+Q[a][1][n]*(conjugate(Qtlower[n*4+0])/alfaQ)+
  Q[a][2][n]*(conjugate(Qtlower[n*4+3])*alfaQ)+Q[a][3][n]*(conjugate(Qtlower[n*4+2])*alfaQ);
		      
}/*n*/
}/*a*/

/* Computing the norm of the deviations of the Q/Qtilde ratios from the gluing constant */

Qalfadev[1]=complex(0,0);

for(k=0;k<=lc-1;k++){ Qalfadev[1]+=(alfaQ-(Qlower[k*4+0]/conjugate(Qlower[k*4+2])))*conjugate(alfaQ-(Qlower[k*4+0]/conjugate(Qlower[k*4+2]))); 
    Qalfadev[1]+=(alfaQ-(Qtlower[k*4+0]/conjugate(Qtlower[k*4+2])))*conjugate(alfaQ-(Qtlower[k*4+0]/conjugate(Qtlower[k*4+2])));
    Qalfadev[1]+=(alfaQ+(Qlower[k*4+1]/conjugate(Qlower[k*4+3])))*conjugate(alfaQ+(Qlower[k*4+1]/conjugate(Qlower[k*4+3])));
    Qalfadev[1]+=(alfaQ+(Qtlower[k*4+1]/conjugate(Qtlower[k*4+3])))*conjugate(alfaQ+(Qtlower[k*4+1]/conjugate(Qtlower[k*4+3])));
}//k

Qalfadev[1]=sqrt(Qalfadev[1]);
Qalfadev[1]/=(4*lc*alfaQ);


delete [] Qlowerm1;

for(a=0;a<=3;a++){
 for(i=0;i<=3;i++){ delete [] Q[a][i]; delete [] Q1[a][i];  }/*i*/
}/*a*/	

	   

for(a=0;a<=3;a++){ delete [] P[a];
                   delete [] Pt[a];} /*a*/

}/* QdoubleLRujicm ends */

/**********************************************************************/



/*****************************************/

/*A function to determine the set of equations obtained after the Fourier-transformation: the double case */ 

/*****************************************/

void doubleQtoEtypeIInewton(cl_N * E,cl_N * deltaP,cl_N * deltaPt/*,cl_N * deltaPup,cl_N * deltaPtup*/,cl_N * CT,cl_N * CU,int Nas[4][2],cl_N * suA/*lc*/,long int lc,long int N00,cl_N g,int ncmax[4],long int Mtint[4],int nmin/*double zerusmodus helye */){
  
  long int k,a,n;
  
  /*symmetric and antisymmetric combinations */
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  cl_N * cS[4];
  
  cl_N * cA[4];
  
  for(a=0;a<=3;a++){ cS[a]=new cl_N[ncmax[a]+1]; cA[a]=new cl_N[ncmax[a]+1]; }//a
  
  
  
  for(a=0;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaP[k*4+a]+deltaPt[k*4+a])/2; 
                          fA[k*4+a]=(deltaP[k*4+a]-deltaPt[k*4+a])/2/suA[k];			  			
    }/*k*/
  }/*a*/
  


for(a=0;a<=3;a++){ 
  for(n=0;n<=ncmax[a];n++){  cS[a][n]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cS[a][n]+=fS[k*4+a]*CT[(lc-1-k)*lc+n]; 
    
  }/*k*/
  
  cS[a][n]/=lc;    
  }/*n*/
}/*a*/


for(a=0;a<=3;a++){ 
  
  cA[a][0]=complex(0,0); 
  
  
  /*n from 1 */ 
  for(n=1;n<=ncmax[a];n++){  cA[a][n]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cA[a][n]+=fA[k*4+a]*CU[(lc-1-k)*lc+n-1]; 
    
  }/*k*/
  
  cA[a][n]/=lc;
  
  cA[a][n]*=2*complex(0,1)*g;
    
  }/*n*/
}/*a*/




/* the c_{a,n} coefficients */
a=0;
k=0;
for(n=0;n<=N00;n++){ if( n==nmin){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[n-k]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[n-k]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

a=1;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[1]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[N00-1+n-k]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[N00-1+n-k]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=2;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[2])  ||  n==(Mtint[1]-Mtint[2]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[2*N00-2+n-k]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[2*N00-2+n-k]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

a=3;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[3]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[3*N00-4+n-k]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[3*N00-4+n-k]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/



delete [] fS;
delete [] fA;


for(a=0;a<=3;a++){ delete [] cS[a]; 
                   delete [] cA[a]; }//a

  
}

/****************************/
/*A function to determine the set of equations obtained after the Fourier-transformation: the single case */
/***************************/

void singleQtoEtypeIInewton(cl_N * E,cl_N * deltaP,cl_N * deltaPt,cl_N * CT,cl_N * CU,int Nas[4][2],cl_N * suA/*lc*/,long int lc,long int N00,cl_N g,int ncmax[4],long int Mtint[4],int nmin,int nmax){
  
  long int k,a,n;
  
  
  /* symmetric and antisymmetric combinations */
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  
  cl_N * cS[4];
  
  cl_N * cA[4];
  
  for(a=0;a<=3;a++){ cS[a]=new cl_N[ncmax[a]+1]; cA[a]=new cl_N[ncmax[a]+1]; }//a
  
  
  for(a=0;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaP[k*4+a]+deltaPt[k*4+a])/2; 
                          fA[k*4+a]=(deltaP[k*4+a]-deltaPt[k*4+a])/2/suA[k];			  
    }/*k*/
  }/*a*/
  

for(a=0;a<=3;a++){ 
  for(n=0;n<=ncmax[a];n++){  cS[a][n]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cS[a][n]+=fS[k*4+a]*CT[(lc-1-k)*lc+n];  
    
  }/*k*/
  
  cS[a][n]/=lc;
    
  }/*n*/
}/*a*/


for(a=0;a<=3;a++){ 
  
  cA[a][0]=complex(0,0); 
  
  
  /*n 1-tol megy */ 
  for(n=1;n<=ncmax[a];n++){  cA[a][n]=complex(0,0); ;
    
  for(k=0;k<=lc-1;k++){   cA[a][n]+=fA[k*4+a]*CU[(lc-1-k)*lc+n-1]; 
    
  }/*k*/
  
  cA[a][n]/=lc;
  
  cA[a][n]*=2*complex(0,1)*g;
    
  }/*n*/
}/*a*/

/* the c_{a,n} coefficients */

a=0;
k=0;
for(n=0;n<=N00;n++){ if( n==nmin || n==nmax ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[n-k]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[n-k]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=1;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[1]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[N00-2+n-k]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[N00-2+n-k]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

a=2;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[2])  ||  n==(Mtint[1]-Mtint[2]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[2*N00-3+n-k]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[2*N00-3+n-k]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=3;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[3]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ E[3*N00-5+n-k]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ E[3*N00-5+n-k]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

  
delete [] fS;
delete [] fA;


for(a=0;a<=3;a++){ delete [] cS[a]; delete [] cA[a]; }//a
  
}

/****************************/


/*A function to determine the derivative with respect to the j-th variable of the equations 
 *obtained after the Fourier-transformation: double case */  

void doubleQtoDEtypeIInewton(cl_N * DE,cl_N * deltaP,cl_N * deltaPt,cl_N * CT,cl_N * CU,int Nas[4][2],cl_N * suA/*lc*/,long int lc,long int N00,cl_N g,int ncmax[4],long int Mtint[4],int nmin,long int dimV,long int j){
  
  long int k,a,n;
  
  
  /*symmetric and antisymmetric combinations*/
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  cl_N * cS[4];
  
  cl_N * cA[4];
  
  for(a=0;a<=3;a++){ cS[a]=new cl_N[ncmax[a]+1]; cA[a]=new cl_N[ncmax[a]+1]; }//a

  
  for(a=0;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaP[k*4+a]+deltaPt[k*4+a])/2; 
                          fA[k*4+a]=(deltaP[k*4+a]-deltaPt[k*4+a])/2/suA[k];			  
    }/*k*/
  }/*a*/
 
for(a=0;a<=3;a++){ 
  for(n=0;n<=ncmax[a];n++){  cS[a][n]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cS[a][n]+=fS[k*4+a]*CT[(lc-1-k)*lc+n]; 
    
  }/*k*/
  
  cS[a][n]/=lc;
    
  }/*n*/
}/*a*/


for(a=0;a<=3;a++){ 
  
  cA[a][0]=complex(0,0); 
  
  
  /*n runs from 1*/ 
  for(n=1;n<=ncmax[a];n++){  cA[a][n]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cA[a][n]+=fA[k*4+a]*CU[(lc-1-k)*lc+n-1];  
    
  }/*k*/
  
  cA[a][n]/=lc;
  
  cA[a][n]*=2*complex(0,1)*g;
    
  }/*n*/
}/*a*/

/* "c_{a,n}" */
a=0;
k=0;
for(n=0;n<=N00;n++){ if( n==nmin){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(n-k)*dimV+j]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(n-k)*dimV+j]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=1;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[1]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(N00-1+n-k)*dimV+j]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(N00-1+n-k)*dimV+j]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

a=2;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[2])  ||  n==(Mtint[1]-Mtint[2]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(2*N00-2+n-k)*dimV+j]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(2*N00-2+n-k)*dimV+j]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=3;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[3]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(3*N00-4+n-k)*dimV+j]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(3*N00-4+n-k)*dimV+j]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


delete [] fS;
delete [] fA;

for(a=0;a<=3;a++){ delete [] cS[a]; delete [] cA[a]; }//a
  
}

/****************************/
/*************************/
/*A function to determine the derivative with respect to the j-th variable of the equations 
 *obtained after the Fourier-transformation: single case */ 

void singleQtoDEtypeIInewton(cl_N * DE,cl_N * deltaP,cl_N * deltaPt,cl_N * CT,cl_N * CU,int Nas[4][2],cl_N * suA/*lc*/,long int lc,long int N00,cl_N g,int ncmax[4],long int Mtint[4],int nmin,int nmax/*single zerusmodusok helye */,long int dimV,long int j){
  
  long int k,a,n;
 
  
  /*symmetric and antisymmetric combinations */
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  cl_N * cS[4];
  
  cl_N * cA[4];
  
  for(a=0;a<=3;a++){ cS[a]=new cl_N[ncmax[a]+1]; cA[a]=new cl_N[ncmax[a]+1]; }//a
  
  for(a=0;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaP[k*4+a]+deltaPt[k*4+a])/2; 
                          fA[k*4+a]=(deltaP[k*4+a]-deltaPt[k*4+a])/2/suA[k];		  
    }/*k*/
  }/*a*/
  

for(a=0;a<=3;a++){ 
  for(n=0;n<=ncmax[a];n++){  cS[a][n]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cS[a][n]+=fS[k*4+a]*CT[(lc-1-k)*lc+n]; 
    
  }/*k*/
  
  cS[a][n]/=lc;
  
  }/*n*/
}/*a*/


for(a=0;a<=3;a++){ 
  
  cA[a][0]=complex(0,0); 
  
  /*n runs from 1 */ 
  for(n=1;n<=ncmax[a];n++){  cA[a][n]=complex(0,0);
    
  for(k=0;k<=lc-1;k++){   cA[a][n]+=fA[k*4+a]*CU[(lc-1-k)*lc+n-1]; 
    
  }/*k*/
  
  cA[a][n]/=lc;
  
  cA[a][n]*=2*complex(0,1)*g;
    
  }/*n*/
}/*a*/

/* "c_{a,n}" */
a=0;
k=0;
for(n=0;n<=N00;n++){ if( n==nmin || n==nmax ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(n-k)*dimV+j]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(n-k)*dimV+j]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=1;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[1]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(N00-2+n-k)*dimV+j]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(N00-2+n-k)*dimV+j]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

a=2;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[2])  ||  n==(Mtint[1]-Mtint[2]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(2*N00-3+n-k)*dimV+j]=imagpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(2*N00-3+n-k)*dimV+j]=imagpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/


a=3;
k=0;
for(n=1;n<=N00;n++){ if( n==(Mtint[0]-Mtint[3]) ){k++; }else{ 
    
    if(n>=Nas[a][0]){ DE[(3*N00-5+n-k)*dimV+j]=realpart(cS[a][(-Nas[a][0]+n)]+cA[a][(-Nas[a][0]+n)]); }else{ DE[(3*N00-5+n-k)*dimV+j]=realpart(cS[a][(Nas[a][0]-n)]-cA[a][(Nas[a][0]-n)]); } 
}
}/*n*/

delete [] fS;
delete [] fA;

for(a=0;a<=3;a++){ delete [] cS[a]; delete [] cA[a]; }//a

}

/****************************/


/* Function computing c_{a,n0} from P and Ptilde: n0=given by an argument */  

void QtoEtypeIIcmgen(long int n0/*-1,..,N0!*/,cl_N * Ev,cl_N * deltaPv,cl_N * deltaPtv,cl_N * CT,cl_N * CU,int Nas[4][2],cl_N * suA/*lc*/,long int lc,cl_N g){
  
  long int k,a,n;
  
  
  /*index vector*/
  long int na[4];
  
  for(a=0;a<=3;a++){ na[a]=labs(n0-Nas[a][0]);  }//a
  
  /*symmetric and antisymmetric combintions*/
  
  cl_N * fS=new cl_N[4*lc];
  cl_N * fA=new cl_N[4*lc];
  
  
  cl_N * cS=new cl_N[4];
  
  cl_N * cA=new cl_N[4];
  
  for(a=0;a<=3;a++){
    for(k=0;k<=lc-1;k++){ fS[k*4+a]=(deltaPv[k*4+a]+deltaPtv[k*4+a])/2; 
                          fA[k*4+a]=(deltaPv[k*4+a]-deltaPtv[k*4+a])/2/suA[k];			  
    }/*k*/
  }/*a*/
  

for(a=0;a<=3;a++){ 
   
  cS[a]=complex(0,0); 
    
  for(k=0;k<=lc-1;k++){   cS[a]+=fS[k*4+a]*CT[(lc-1-k)*lc+na[a]];
    
  }/*k*/
  
  cS[a]/=lc;
 
}/*a*/


for(a=0;a<=3;a++){ 
    
    cA[a]=complex(0,0); 
    
    if(na[a]==0){}else{
    
  for(k=0;k<=lc-1;k++){   cA[a]+=fA[k*4+a]*CU[(lc-1-k)*lc+na[a]-1]; 
    
  }/*k*/
  
  cA[a]/=lc;
  
  
  cA[a]*=2*complex(0,1)*g;
    }//else ends
}/*a*/

  
for(a=0;a<=3;a++){ 
    
    if(n0>=Nas[a][0]){ Ev[a]=cS[a]+cA[a]; }else{ Ev[a]=cS[a]-cA[a]; } 
   
}/*a*/

delete [] fS;
delete [] fA;
delete [] cS;
delete [] cA;

  
}



/*************************************************************************/
/*************************************************************************/
/*************************************************************************/ 
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/ 
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/ 
 
 
int main(int argc, char * argv[]) {
  
long int i,j,k,n,a,m,b0;  
  
 long int workingprecision=atol(argv[1]); /*working precision*/

float_format_t precision=float_format(workingprecision);

long int N0=atol(argv[2]);  /* cutoff parameter c_ {a, n} n = 0, ..., N0 */

long int NCb=atol(argv[3]);
/*cutoff parameter, the last kept term in the series of Q_ {a | i} is u^(-Mhat[i])/u^(NCb) */

long int NIcutoff=atol(argv[4]); /* pull down cutoff */

long int lc=atol(argv[5]);  /*number of Chebyshev sampling points on the cut */

long int DH=atol(argv[6]);  /* integer derivative parameter \[Epsilon] = 10^(-DH) */

long int precssf=atol(argv[7]); /* integer quit condition parameter */

long int precDelta=atol(argv[8]); /* integer quit condition parameter */

long int maxiter=atol(argv[9]);/* iteration quit condition parameter, if the number of iterations is 
larger than maxiter, than the programm finishes running, and returns
the value, it obtained upto maxiter iterations */

/* The quit condition is as follows: either "number of iterations" > maxiter, 
                                      
 or the simultaneous fullfillment of two conditions is required:
                                         
|Delta^{(n)}-Delta^{(n-1)}|<10^(-precDelta) and 
"the sum of |c_{a,n}-(c_{a,n} from gluing)|^2"<10^(-precssf), 

where Delta^{(n)} is the value of the anomalous part dimension after the n-th iteration.*/


long int nb[2];
long int na[2];
long int nf[4];

nb[0]=atol(argv[10]);//n_{b_1}
nb[1]=atol(argv[11]);//n_{b_2}

nf[0]=atol(argv[12]);//n_{f_1}
nf[1]=atol(argv[13]);//n_{f_2}
nf[2]=atol(argv[14]);//n_{f_3}
nf[3]=atol(argv[15]);//n_{f_4}

na[0]=atol(argv[16]);//n_{a_1}
na[1]=atol(argv[17]);//n_{a_2}

long int multiplicity=atol(argv[18]); /*multiplicity*/


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

long int N0be=(argc-25)/4;

/*initial values of c_{a,n}*/

cl_R * c0[4];

for(a=0;a<=3;a++){ c0[a]=new cl_R[N0be+1]; }/*a*/

for(j=0;j<=N0be;j++){    
  
                    stext<<argv[21+j]; c0[0][j]=cl_float((cl_R) stext.str().c_str(),precision); 
                         stext.str(std::string());
		       
		    stext<<argv[22+N0be+j]; c0[1][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                         stext.str(std::string());   
		       
		       stext<<argv[23+2*N0be+j]; c0[2][j]=cl_float((cl_R) stext.str().c_str(),precision);  
                        stext.str(std::string());
		       
		       stext<<argv[24+3*N0be+j]; c0[3][j]=cl_float((cl_R) stext.str().c_str(),precision);
                     stext.str(std::string());	       
}/*j*/

/* some useful constants */

cl_N Pi=pi(precision);
cl_R two=cl_float(2,precision);
cl_R egy=cl_float(1,precision);
cl_R half=egy/two;
cl_N II=complex(0,1);



/* "lattice constant" for differentiation */ 
 const cl_N h=expt(cl_float(10,precision),(-1)*DH);
 


/* Computing the levels of the zero-mode equations */

long int ni[4]; 

ni[0]=nb[0]-nb[1]-1;
ni[1]=nb[1]-nb[0]+1;
ni[2]=na[1]-na[0]-1;
ni[3]=na[0]-na[1]+1; 


/* (nmin,imin) and (namx,imax) are the two (n,i) pairs, where the zero-mode problem arises */

long int nmin;
long int nmax;

long int imin;
long int imax;
/* auxiliary index */
long int index;

/* only the 2 positive elements should be found */

nmax=0;
imax=0;
index=0;

for(i=0;i<=3;i++) { if(ni[i]>0 && index==0){imin=i; nmin=ni[i]; index++; }else{}/*1. if ends */
                    if(ni[i]>0 && index==1){imax=i; nmax=ni[i];  }else{}/*2. if ends */}/*i*/

if(nmin>nmax){ index=nmax; nmax=nmin; nmin=index;
               index=imax; imax=imin; imin=index; /* ordering */
                
}else{}


/* an indicator for selecting between the double and single cases */
long int Lniminmax;

if(nmax==nmin){ Lniminmax=3; /*double case: nmin=nmax*/}
else{ Lniminmax=4; /*single case: nmin < nmax */ }


/* Storing the (n,i) pairs in appropriate vectors*/

long int * niminmax=new long int[Lniminmax]; 

if(Lniminmax==3){//double case
niminmax[0]=nmin;
niminmax[1]=imin;
niminmax[2]=imax;}else{ /*single case*/
  niminmax[0]=nmin;
niminmax[1]=nmax;
niminmax[2]=imin;
niminmax[3]=imax;
}/*if-else ends*/



/* An arbitrary 4-dimensional vector for treating the zero-mode ambiguity for Q_i[i] */

cl_N lambda0[4];
/* we fix it to be zero for the entire code, but in principle any value could be good*/
for(i=0;i<=3;i++){ lambda0[i]=complex(0,0);}


/* Filling Mtilde vector with its integer values */

long int Mtint[4];

Mtint[0]=nf[0]+2;
Mtint[1]=nf[1]+1;
Mtint[2]=nf[2];
Mtint[3]=nf[3]-1;

/* for the LR-symmetric formulation the double of the scaling parameter */

long int kettoLAMBDA=(1-Mtint[0]-Mtint[3]);


/* Constructing Mtilde for the LR-symmetric formulation */

cl_N Mt[4];

for(a=0;a<=3;a++){ Mt[a]=cl_float((cl_I) Mtint[a],precision)+cl_float((cl_I) kettoLAMBDA,precision)/cl_float(2,precision);  }

/* The double of the LR-symmetric Mtilde is always integer */

long int twiceMt[4];

for(a=0;a<=3;a++){ twiceMt[a]=2*Mtint[a]+kettoLAMBDA; }/*a*/

/*The value of L */

cl_N L=cl_float((cl_I) nf[0]+nf[1]+nf[2]+nf[3]+na[0]+na[1]-nb[0]-nb[1],precision)/two;
long int Lint=(nf[0]+nf[1]+nf[2]+nf[3]+na[0]+na[1]-nb[0]-nb[1])/2;



/* initializing Delta*/
cl_N Delta=cl_float(Delta0,precision);



/* Mhat at g=0 */
cl_N Mhat0[4];

Mhat0[0]=L+cl_float((cl_I) nb[0]+1,precision)+cl_float((cl_I) kettoLAMBDA,precision)/two;
Mhat0[1]=L+cl_float((cl_I) nb[1]+2,precision)+cl_float((cl_I) kettoLAMBDA,precision)/two;
Mhat0[2]=cl_float((cl_I) -1-na[0],precision)+cl_float((cl_I) kettoLAMBDA,precision)/two;
Mhat0[3]=cl_float((cl_I) -na[1],precision)+cl_float((cl_I) kettoLAMBDA,precision)/two;

/* Mhat at the begining. I.e initial value. */
cl_N Mhat[4];

Mhat[0]=L+cl_float((cl_I) nb[0]+1,precision)+Delta/two+cl_float((cl_I) kettoLAMBDA,precision)/two;
Mhat[1]=L+cl_float((cl_I) nb[1]+2,precision)+Delta/two+cl_float((cl_I) kettoLAMBDA,precision)/two;
Mhat[2]=cl_float((cl_I) -1-na[0],precision)-Delta/two+cl_float((cl_I) kettoLAMBDA,precision)/two;
Mhat[3]=cl_float((cl_I) -na[1],precision)-Delta/two+cl_float((cl_I) kettoLAMBDA,precision)/two;



/* Filling up c_{a,n}'s with their initial values */

long int N0c=max(nmax,N0);

cl_N * c[4];
cl_N * cf[4];


for(a=0;a<=3;a++){ c[a]=new cl_N[N0c+1];
                   cf[a]=new cl_N[N0c+1];   
}


for(a=0;a<=3;a++){
  for(n=0;n<=N0;n++){  
    
    if(n<=N0be){/*if-n begins*/
    
    if(a==0 || a==2){  c[a][n]=II*cl_float(c0[a][n],precision);
                        }else{
                                                   c[a][n]=cl_float(c0[a][n],precision);
                                                         }/*if-a ends*/
}else{ c[a][n]=complex(0,0);
       }/*if -n ends*/
  }/*n*/
for(n=N0+1;n<=N0c;n++){
    c[a][n]=complex(0,0);
    
  }/*n*/
}/*a*/

/* c0 can be deleted */
for(a=0;a<=3;a++){ delete [] c0[a];  }/*a*/
  
for(a=0;a<=3;a++){
  for(n=0;n<=N0;n++){ cf[a][n]=c[a][n]; }//i
}//a  

/*update N0 if it is necessary*/
N0=N0c;



cl_N A[4]; /* A_a */
cl_N Af[4]; /* A^a */
cl_N AAproduct[4];

AAfunc(Mt,Mhat,AAproduct);

/* Giving the initial values for A_a and c_{a,0} */
/*for a=1,2 */
for(a=0;a<=1;a++){
A[a]=VolinAfunc(a,Mt,Mhat); 
Af[a]=AAproduct[a]/A[a];
c[a][0]=A[a]/expt(g,Mt[a]);
}/*a*/
 
/* for a=3,4 */
 
A[2] = AAproduct[1]/A[1];
c[2][0] = A[2]/expt(g,Mt[2]); 
Af[2] = AAproduct[2]/A[2];

 
A[3] =cl_float(-1,precision)*AAproduct[0]/A[0];
c[3][0] = A[3]/expt(g,Mt[3]); 
Af[3] = AAproduct[3]/A[3];



/*  4x4 matrix: A_a*A^b */

cl_N AA[4][4];
for(a=0;a<=3;a++){ 
  for(b0=0;b0<=3;b0++){  AA[a][b0]=A[a]*Af[b0]; 
  }/*b0*/
  }/*a*/


cl_N B[4]; /*The B_i vector and its initialization*/

/*B_1 and B_2*/


for(i=0;i<=1;i++){ B[i]=Bnikafunc(i,Mt,Mhat);  }

/*B_3 and B_4*/

B[2]=BBfunc(1,Mt,Mhat)/B[1];

B[3]=BBfunc(3,Mt,Mhat)/B[0];


/*  The B_a|i matrix */

cl_N BB[4][4];
/* and two other useful auxiliary matrices */
cl_N alfa[4][4];
cl_N beta[4][4];

for(a=0;a<=3;a++){ 
  for(i=0;i<=3;i++){  alfa[a][i]=Mhat[i]-Mt[a];
                      beta[a][i]=Mt[a]-Mt[i];
  BB[a][i]=II*A[a]*B[i]/(Mt[a]-Mhat[i]);
 }/*i*/
}/*a*/

/*two useful auxiliary vectors */

cl_N * gvmin;
gvmin=new cl_N[4];

for(b0=0;b0<=3;b0++){
  
gvmin[b0]=expt(g,Mt[b0]+nmin); }


cl_N * gvmax;
gvmax=new cl_N[4];

for(b0=0;b0<=3;b0++){
  
gvmax[b0]=expt(g,Mt[b0]+nmax); }




/* Important auxiliary array*/
int Nas[4][2]; /*Na*/

for(a=0;a<=3;a++){ 
  
  
  if(kettoLAMBDA%2==0) { Nas[a][0]=(-1)*Mtint[a]-kettoLAMBDA/2;  Nas[a][1]=Mtint[a]+kettoLAMBDA/2-1; }else{ Nas[a][0]=(-1)*Mtint[a]+(Mtint[0]+Mtint[3])/2;  Nas[a][1]=Mtint[a]-1-(Mtint[0]+Mtint[3])/2; } }/*a*/




/* correcting lc, if it is too small for Fourier-transformation (to avoid segmentation fault) */

int lcv[8];
 
lcv[0]=1+abs(1+Nas[0][0]);
lcv[1]=1+abs(1+Nas[1][0]);
lcv[2]=1+abs(1+Nas[2][0]);
lcv[3]=1+abs(1+Nas[3][0]);

lcv[4]=1+abs(N0-Nas[0][0]);
lcv[5]=1+abs(N0-Nas[1][0]);
lcv[6]=1+abs(N0-Nas[2][0]);
lcv[7]=1+abs(N0-Nas[3][0]);

/*lc must be greater, than the maximum */

for(j=0;j<=7;j++){ if(lc<lcv[j]){ lc=lcv[j]; }   }//j

/***************************************************/
/* chat - nmax index set for FFT*/
int ncmax[4];

for(a=0;a<=3;a++){ ncmax[a]=abs(N0-Nas[a][0]); if(ncmax[a]<abs(1+Nas[a][0])){ ncmax[a]=abs(1+Nas[a][0]); }  }//a


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

for(i=0;i<=3;i++){ NQT[i]=NCb; }/*i*/

long int NQ;


NQ=NQT[0];
  
for(j=1;j<=3;j++){ if(NQT[j]>NQ){ NQ=NQT[j]; } }/*j*/



/* (-I*2)^(j) vector j=0,...,NQ-1 */

cl_N * m2ij=new cl_N[NQ];

for(j=0;j<=NQ-1;j++){ m2ij[j]=expt(complex(0,cl_float(-2,precision)),j);
}


/* useful constants*/

cl_R fel=cl_float(1,precision)/cl_float(2,precision);
cl_R mnegyed=cl_float(-1,precision)/cl_float(4,precision);

/* (I/2)^(n) vector n=0,..,NQ */

cl_N * ip2n=new cl_N[NQ+1];

for(j=0;j<=NQ;j++){ ip2n[j]=expt(complex(0,fel),j);
}

/* definition of lmax  */

long int lmax=ldiv(NQ,2L).quot;

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


/* array size constant */

long int Lmax=ldiv(NQ-1L,2L).quot;//floor((NQ-1)/2.0);

/*auxiliary array fo computing b_{a|i,n}*/

cl_N * y1=new cl_N[(Lmax+1)*(Lmax+1)];


Y1func(y1,m1p4k,Lmax,precision);


/*auxiliary array fo computing b_{a|i,n}*/

cl_N * y2=new cl_N[(Lmax+1)*(Lmax+1)];

Y2func(y2,m1p4k,Lmax,precision);


/*auxiliary array fo computing b_{a|i,n}*/

cl_N * y31=new cl_N[(2*Lmax+1)*(Lmax+1)];

Y31func(y31,m4k,Lmax,precision);


/*auxiliary array fo computing b_{a|i,n}*/

cl_N * s1=new cl_N[(NQ-1)*(NQ+1)];


S1func(s1,ip2n,NQ,precision);


/*auxiliary array fo computing b_{a|i,n}*/

cl_N * s31=new cl_N[(NQ-1)*(NQ-1)];

S31func(s31,m2ij,NQ,precision);



/* array size constant*/
long int LN0=ldiv(N0c,2L).quot;


/* Auxiliary array for computing the large u-expansion of P_a[u] */
cl_N * SigmaSub[4];
for(a=0;a<=3;a++) { SigmaSub[a]=new cl_N[(1+NQ)*(1+LN0)]; }


for(a=0;a<=3;a++){
sigmasubfunc(SigmaSub[a],twiceMt[a],N0c,NQ,g,precision);
}/*a*/



/* Integer cutoffs for imaginary "pull down process" to determine Q_i[u] at the sampling points of the cut*/

unsigned long int * NIT=new unsigned long int[4];

for(i=0;i<=3;i++){ NIT[i]=NIcutoff; }//a


/*(Delta-dependent) nullvectors of the zero-mode problems*/
cl_N trnullmin[4];/*belonging to imin*/
cl_N trnullmax[4];/*belonging to imax*/

nullTRvectorfunc(imin,A,alfa,beta,trnullmin);
nullTRvectorfunc(imax,A,alfa,beta,trnullmax);


/*cutoff parameter used for all i for the large u-expansions of Q_{a|i}[u] */

NQ=NQT[0];
  
for(j=1;j<=3;j++){ if(NQT[j]>NQ){ NQ=NQT[j]; } }
  
/* Integer cutoff used for all i values as the imaginary pull down cutoff parameter for the determination of Q_i[u] at the sampling points of the cut*/

unsigned long int NI=NIT[0];
  
for(j=1;j<=3;j++){ if(NIT[j]>NI){ NI=NIT[j]; } }  
  

/* integers for certain array sizes*/

long int lmaxT[4];
long int LmaxT[4];

for(i=0;i<=3;i++) { lmaxT[i]=ldiv(NQT[i],2L).quot;
  if(NQT[i]==0){ LmaxT[i]=-1; }else{ LmaxT[i]=ldiv(NQT[i]-1,2L).quot; } }/*i*/
 



/**************************/
/*************************/
/***************************/

/* Auxiliary arrays for P_a[u] and Ptilde_a[u] at the sampling points of the cut */
cl_N * TXm=new cl_N[(N0+1)*lc];
cl_N * TXmi=new cl_N[(N0+1)*lc];


cl_N * Xrealaux=new cl_N[lc]; 
cl_N * Xrealauxhalf=new cl_N[lc]; 


for(k=0;k<=lc-1;k++){  Xrealaux[k]=X(uA[k]/g);  Xrealauxhalf[k]=sqrt(Xrealaux[k]); }



for(n=0;n<=N0;n++){
  for(k=0;k<=lc-1;k++){ TXm[n*lc+k]=expt(Xrealaux[k],n); 
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

/* filling them up with values */

for(a=0;a<=3;a++){
 for(k=0;k<=lc-1;k++){
   for(n=0;n<=N0;n++){  PaT[(n*lc+k)*4+a]=Txmsub[k*4+a]*TXm[n*lc+k];
                        PtaT[(n*lc+k)*4+a]=Txmsubi[k*4+a]*TXmi[n*lc+k]; 
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
    for(n=0;n<=NI-1;n++){  TxAnm[(m*NI+n)*lc+k]=egy/expt(XAn[n*lc+k],m);  }/*n*/
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

for(a=0;a<=3;a++){
  for(k=0;k<=lc-1;k++){
    for(n=0;n<=NI-1;n++){
      for(m=0;m<=N0;m++){  PujT[((m*NI+n)*lc+k)*4+a]=TxMt[(n*lc+k)*4+a]*TxAnm[(m*NI+n)*lc+k];
	                    
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
  for(n=0;n<=NQT[i];n++){  TuANIn[i][(n*lc+k)]=egy/expt(uA[k]+II*(NIT[i]+half),n);
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

/* For error estimation: for the Qconstructor...cm() function , the previous arrays corresponding to NI-1 cutoff*/

cl_N * TuANInshift[4];
for(i=0;i<=3;i++){ TuANInshift[i]=new cl_N[lc*(NQT[i]+1)]; }/*i*/

for(i=0;i<=3;i++){
for(k=0;k<=lc-1;k++){ 
  for(n=0;n<=NQT[i];n++){  TuANInshift[i][n*lc+k]=egy/expt(uA[k]+II*(((int) NIT[i])-1+half),n);
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


/**********************************************************************************************/
/* Memory allocatioin and initialization of some important Delta-dependent auxiliary arrays */
/*******************************************************************************************/
cl_N * alfaais[4];  
cl_N * T5[4]; 
cl_N * T3[4]; 
cl_N * Y32[4]; 
cl_N * S1n[4];  
cl_N * S32[4]; 

for(i=0;i<=3;i++){
 alfaais[i]=new cl_N[4*(NQT[i]+2)];
 T5[i]=new cl_N[4*(lmaxT[i]+1)*(2*lmaxT[i]-2)];
 T3[i]=new cl_N[4*(lmaxT[i]+1)];
 Y32[i]=new cl_N[4*(2*LmaxT[i]+1)*(LmaxT[i]+1)];
 S1n[i]=new cl_N[4*(NQT[i]+1)];
 S32[i]=new cl_N[4*(NQT[i]+1)*(NQT[i]-1)];  }/*i*/
  
  
  /********************************/
  /********************************/
  /********************************/


/*Arrays for Q_i[u] and Qtilde_i[u] at the sampling points: Qlower[n*4+i-1]=Q_i[u_n] etc.*/
cl_N * Qlower=new cl_N[4*lc]; 
cl_N * Qtlower=new cl_N[4*lc];  

/* P_a[u_n]-(P_a[u_n] from gluing)=deltaP[n*4+a-1] , and its tilded version as well*/
cl_N * deltaP=new cl_N[4*lc];  
cl_N * deltaPt=new cl_N[4*lc];  

/*Arrays for values of P_a[u] computed back from gluing*/
cl_N * deltaPv;
cl_N * deltaPtv;

/* array to store c_{a,-1} computed back from gluing condition */
cl_N * Ev=new cl_N[4];



cl_N * scT[4];

for(i=0;i<=3;i++){

scT[i]=new cl_N[4*4*(NQT[i]+1)];

/*computation of the 4x4 matrices for the linear problems for b_{a|i,n}*/

totalscTmakerLRi(i,AA,BB,alfa,NQT[i],scT[i]);

/*computation of Delta dependent auxiliary arrays*/
AlfatoTLRi(i,alfa,alfaais[i],T3[i],T5[i],Y32[i],S1n[i],S32[i],ip2n,m1p4k,lmaxT[i],LmaxT[i],NQT[i],precision);

}/*i*/


long int l;


/* Arrays for storing some numbers for running diagnostics */
cl_N * Qalfadev=new cl_N[2] ;

cl_R * Qalfadevout=new cl_R[(maxiter+1)*2];


for(i=0;i<=(maxiter+1)*2-1;i++){ Qalfadevout[i]=0; }//i

   
   
/* comuting left nullvectors*/
nullTRvectorfunc(imin,A,alfa,beta,trnullmin);
nullTRvectorfunc(imax,A,alfa,beta,trnullmax);


/*variable to store the values of P_a[u] computed back from gluing condition*/
deltaPv=new cl_N[4*lc];
deltaPtv=new cl_N[4*lc];

/* Computing the Q_i[u], P_a[u]-(P_a[u] from gluing) and the values of P_a[u] computed back from gluing condition and their tilded versions as well*/

if(Lniminmax==3){ //double case
QdoubleLRujicm(Qlower,Qtlower,deltaP,deltaPt,deltaPv,deltaPtv,A,AA,BB,alfa,c,SigmaSub,PaT,PtaT,PujT,t1,t2,t41,T3,T5,y1,y2,y31,Y32,S1n,s1,s31,S32,TuANIn,TuMiMaNI,g,N0,NQT,NQ,NIT,NI,lc,ni,lambda0,niminmax,gvmin,gvmax,trnullmin,trnullmax,scT,/*aprime,*/precision,TuANInshift,TuMiMaNIshift,Qalfadev);   
  }
else{ //single case
  QsingleLRujicm(Qlower,Qtlower,deltaP,deltaPt,deltaPv,deltaPtv,A,AA,BB,alfa,c,SigmaSub,PaT,PtaT,PujT,t1,t2,t41,T3,T5,y1,y2,y31,Y32,S1n,s1,s31,S32,TuANIn,TuMiMaNI,g,N0,NQT,NQ,NIT,NI,lc,ni,lambda0,niminmax,gvmin,gvmax,trnullmin,trnullmax,scT,precision,TuANInshift,TuMiMaNIshift,Qalfadev); }
  
Qalfadevout[0*2+0]=realpart(Qalfadev[0]);
Qalfadevout[0*2+1]=realpart(Qalfadev[1]);

delete [] deltaPv;
delete [] deltaPtv;

long int n0=7;
long int i0=1;


/* the dimension of independent variables */
long int dimV;

if(Lniminmax==3){ dimV=4*N0-4; }else{ dimV=4*N0-5; };


/*Auxiliary array for storing the correction coming from Newton's method*/
cl_N * dc=new cl_N[dimV];

int * errorcode=new int[1];
/*Setting inital value for the errocode */
errorcode[0]=0;

/*Initialization of the independent variables of the problem */
cl_N * V=new cl_N[dimV];


if(Lniminmax==3){ doubleCtoVLR(Delta,nmin,nmax,c,Mtint,V,N0);
  }else{ singleCtoVLR(Delta,nmin,nmax,c,Mtint,V,N0); }
  


/* A new copy of the true variables*/

cl_N * V0=new cl_N[dimV];

for(j=0;j<=dimV-1;j++){ V0[j]=V[j]; }/*j*/
  

/* Putting the inital values of the equations into a vector denoted by E */

long int dimE=dimV;

cl_N * E=new cl_N[dimE];


if(kettoLAMBDA%2==1 || kettoLAMBDA%2==-1){
for(a=0;a<=3;a++){ for(k=0;k<=lc-1;k++){ deltaP[k*4+a]/=Xrealauxhalf[k]; deltaPt[k*4+a]*=Xrealauxhalf[k]; }/*k*/  }/*a*/ }



if(Lniminmax==3){ //double
  doubleQtoEtypeIInewton(E,deltaP,deltaPt,CT,CU,Nas,suA,lc,N0,g,ncmax,Mtint,nmin);  }
else{ //single
  singleQtoEtypeIInewton(E,deltaP,deltaPt,CT,CU,Nas,suA,lc,N0,g,ncmax,Mtint,nmin,nmax); }
  
    
/*Making a copy of the equations*/  

cl_N * E0=new cl_N[dimE];

for(j=0;j<=dimE-1;j++){ E0[j]=E[j]; }/*j*/


  
/* An auxiliary vector for bookkeeping the realness or imaginariness of the variables*/

cl_N * PhiV=new cl_N[dimV]; 

 PhiV[0]=complex(1,0); 
 
if(Lniminmax==3){ 

for(k=1;k<=N0-1;k++){ PhiV[k]=complex(0,1); }

for(k=N0;k<=2*N0-2;k++){ PhiV[k]=complex(1,0); }

for(k=2*N0-1;k<=3*N0-4;k++){ PhiV[k]=complex(0,1); }

for(k=3*N0-3;k<=4*N0-5;k++){ PhiV[k]=complex(1,0); } }else{
  
for(k=1;k<=N0-2;k++){ PhiV[k]=complex(0,1); }

for(k=N0-1;k<=2*N0-3;k++){ PhiV[k]=complex(1,0); }

for(k=2*N0-2;k<=3*N0-5;k++){ PhiV[k]=complex(0,1); }

for(k=3*N0-4;k<=4*N0-6;k++){ PhiV[k]=complex(1,0); }  
  
  
}//else ends





/*Norm square of "equations" */
cl_R Enormainf=realpart(norma2(E0,dimE));

/*Norm square of "equations" */
cl_R Enormainfold;

/* Variable to store the norm square of c_{a,-1} computed back from gluing*/
cl_R Evnorma=888;

/* Actual value of Delta*/
cl_R Deltainf=realpart(V0[0]);



             

/* iteration number  */
long int * it=new long int[1];


/*Derivative matrix*/  
cl_N * DEplus=new cl_N[dimE*dimV];


/* Variable to store |Delta^{(n)}-Delta^{(n-1)}| */
cl_R DDelta=cl_float(9,precision);

/* Upper bound for the norm squared of the equations in the quit condition*/
cl_R ssfmax=expt(cl_float(10,precision),(-1)*precssf);
 

it[0]=0L;
while( !( Enormainf<ssfmax && DDelta<expt(cl_float(10,precision),(-1)*precDelta) ) && it[0]<maxiter){// <- quit conditions

  
 /* Making a copy of the norm squared of the equations obtained after the previous ietration */
Enormainfold=Enormainf; 



/*********************************************************/
/* Here starts the compuation of the derivative matrix */
/********************************************************/

for(n=0;n<=dimV-1;n++){//cycle for shifthing the n-th varible with the "lattice constant"
  
V[n]+=PhiV[n]*h;


/*Delta and  the c's from the shifted variables*/
if(Lniminmax==3){ 
doubleVtoCLR(&Delta,nmin,nmax,c,Mtint,V,N0); }
else{ singleVtoCLR(&Delta,nmin,nmax,c,Mtint,V,N0); }


/* computing the Delta dependent auxiliary arrays and the necessary nullvectors */
if(n==0 || n==1) {   Ca0funcLR(Delta,c,Mt,Mhat0,Mhat,A,AA,B,BB,alfa,g,precision);
  
for(i=0;i<=3;i++){  
  
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){ 
    TuMiMaNI[(k*4+a)*4+i]=expt(uA[k]+II*(NIT[i]+half),Mhat[i]-Mt[a]);   
}/*a*/
  }/*k*/   

totalscTmakerLRi(i,AA,BB,alfa,NQT[i],scT[i]);

AlfatoTLRi(i,alfa,alfaais[i],T3[i],T5[i],Y32[i],S1n[i],S32[i],ip2n,m1p4k,lmaxT[i],LmaxT[i],NQT[i],precision);

}/*i*/

nullTRvectorfunc(imin,A,alfa,beta,trnullmin);
nullTRvectorfunc(imax,A,alfa,beta,trnullmax);

}/*if -n ends*/


/*computing Q_i[u] and Qtilde_i[u] from the auxiliary arrays */
if(Lniminmax==3){ //double case
QdoubleLRuji(Qlower,Qtlower,deltaP,deltaPt,A,AA,BB,alfa,c,SigmaSub,PaT,PtaT,PujT,t1,t2,t41,T3,T5,y1,y2,y31,Y32,S1n,s1,s31,S32,TuANIn,TuMiMaNI,g,N0,NQT,NQ,NIT,NI,lc,ni,lambda0,niminmax,gvmin,gvmax,trnullmin,trnullmax,scT,precision);   
  }
else{ //single case
  QsingleLRuji(Qlower,Qtlower,deltaP,deltaPt,A,AA,BB,alfa,c,SigmaSub,PaT,PtaT,PujT,t1,t2,t41,T3,T5,y1,y2,y31,Y32,S1n,s1,s31,S32,TuANIn,TuMiMaNI,g,N0,NQT,NQ,NIT,NI,lc,ni,lambda0,niminmax,gvmin,gvmax,trnullmin,trnullmax,scT,precision); }
		 
/*computing the "shifted" equations and puting them to a matrix */

if(kettoLAMBDA%2==1 || kettoLAMBDA%2==-1){ 
for(a=0;a<=3;a++){ for(k=0;k<=lc-1;k++){ deltaP[k*4+a]/=Xrealauxhalf[k]; deltaPt[k*4+a]*=Xrealauxhalf[k]; }/*k*/  }/*a*/ }


if(Lniminmax==3){ //double case
   doubleQtoDEtypeIInewton(DEplus,deltaP,deltaPt,CT,CU,Nas,suA,lc,N0,g,ncmax,Mtint,nmin,dimV,n);   }
else{ //single case
  singleQtoDEtypeIInewton(DEplus,deltaP,deltaPt,CT,CU,Nas,suA,lc,N0,g,ncmax,Mtint,nmin,nmax,dimV,n); }

/*shifting back the variables into their initial values*/
V[n]-=PhiV[n]*h;
  
}/* cycle ends*/


/* computation of the derivative matrix out of the "shifted equations" */

for(j=0;j<=dimE-1;j++){
  for(k=0;k<=dimV-1;k++){  DEplus[j*dimV+k]-=E0[j];  DEplus[j*dimV+k]/=(h*PhiV[k]);
     }/**/
}/*j*/


/*********************************************************/
/* Here ends the compuation of the derivative matrix */
/********************************************************/


/* Solving the linear problem to get the change of variables in the Newton-method*/
linsolvepE(DEplus,E0,dc,dimV);


/* Updating the variables with their new values */

for(j=0;j<=dimV-1;j++){      V[j]=V0[j]-dc[j]; }/*j*/


/*The computation of the new values of Delta and the c's from the updated variables*/
  
if(Lniminmax==3){ 
doubleVtoCLR(&Delta,nmin,nmax,c,Mtint,V,N0); }
else{ singleVtoCLR(&Delta,nmin,nmax,c,Mtint,V,N0); }

/* Forcing them to be purely real or imaginary */


  for(n=0;n<=N0;n++){  c[0][n]=complex(0,imagpart(c[0][n]));
                       c[2][n]=complex(0,imagpart(c[2][n]));
		       c[1][n]=complex(realpart(c[1][n]),0);
		       c[3][n]=complex(realpart(c[3][n]),0);
  }//n

Delta=complex(realpart(Delta),0);  





/*Computing Delta-dependent auxiliary quantities*/
Ca0funcLR(Delta,c,Mt,Mhat0,Mhat,A,AA,B,BB,alfa,g,precision);

for(i=0;i<=3;i++){  
  
for(a=0;a<=3;a++){
for(k=0;k<=lc-1;k++){ 
    TuMiMaNI[(k*4+a)*4+i]=expt(uA[k]+II*(NIT[i]+half),Mhat[i]-Mt[a]); 
    TuMiMaNIshift[(k*4+a)*4+i]=expt(uA[k]+II*(((int) NIT[i])-1+half),Mhat[i]-Mt[a]);
}/*a*/
  }/*k*/   


totalscTmakerLRi(i,AA,BB,alfa,NQT[i],scT[i]);

AlfatoTLRi(i,alfa,alfaais[i],T3[i],T5[i],Y32[i],S1n[i],S32[i],ip2n,m1p4k,lmaxT[i],LmaxT[i],NQT[i],precision);

}/*i*/


/* updating the nullvectors*/

nullTRvectorfunc(imin,A,alfa,beta,trnullmin);
nullTRvectorfunc(imax,A,alfa,beta,trnullmax);

/*variable to store the values of P_a[u] computed back from gluing condition*/
deltaPv=new cl_N[4*lc];
deltaPtv=new cl_N[4*lc];


/* Computing the Q_i[u], Qtilde_i[u], and the values of P_a[u] computed back from gluing condition*/
if(Lniminmax==3){ //double case
QdoubleLRujicm(Qlower,Qtlower,deltaP,deltaPt,deltaPv,deltaPtv,A,AA,BB,alfa,c,SigmaSub,PaT,PtaT,PujT,t1,t2,t41,T3,T5,y1,y2,y31,Y32,S1n,s1,s31,S32,TuANIn,TuMiMaNI,g,N0,NQT,NQ,NIT,NI,lc,ni,lambda0,niminmax,gvmin,gvmax,trnullmin,trnullmax,scT,/*aprime,*/precision,TuANInshift,TuMiMaNIshift,Qalfadev);   
  }
else{ //single case
  QsingleLRujicm(Qlower,Qtlower,deltaP,deltaPt,deltaPv,deltaPtv,A,AA,BB,alfa,c,SigmaSub,PaT,PtaT,PujT,t1,t2,t41,T3,T5,y1,y2,y31,Y32,S1n,s1,s31,S32,TuANIn,TuMiMaNI,g,N0,NQT,NQ,NIT,NI,lc,ni,lambda0,niminmax,gvmin,gvmax,trnullmin,trnullmax,scT,precision,TuANInshift,TuMiMaNIshift,Qalfadev); }
  
Qalfadevout[(it[0]+1)*2+0]=realpart(Qalfadev[0]);
Qalfadevout[(it[0]+1)*2+1]=realpart(Qalfadev[1]);  
  


/* Computing the equations out of the new data */

if(kettoLAMBDA%2==1 || kettoLAMBDA%2==-1){
for(a=0;a<=3;a++){ for(k=0;k<=lc-1;k++){ deltaP[k*4+a]/=Xrealauxhalf[k]; deltaPt[k*4+a]*=Xrealauxhalf[k]; 
                                         deltaPv[k*4+a]/=Xrealauxhalf[k]; deltaPtv[k*4+a]*=Xrealauxhalf[k];  
       
}/*k*/  }/*a*/ }



if(Lniminmax==3){ //double case
  doubleQtoEtypeIInewton(E,deltaP,deltaPt,CT,CU,Nas,suA,lc,N0,g,ncmax,Mtint,nmin);  }
else{ //single case
  singleQtoEtypeIInewton(E,deltaP,deltaPt,CT,CU,Nas,suA,lc,N0,g,ncmax,Mtint,nmin,nmax); }

  
  
/*Computing c_{a,-1} from gluing condition*/
QtoEtypeIIcmgen(-1L,Ev,deltaPv,deltaPtv,CT,CU,Nas,suA,lc,g);



delete [] deltaPv;
delete [] deltaPtv;


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
    
  /* Computing the norm squared of c_{a,2} obtained after the iteration */  
    Evnorma=realpart(norma2(Ev,4));

it[0]+=1L;}/*The "big" iteration ends*/


/*If the required convergence and accuracy couldn't be reach within maxiter number of iterations, the value of errorcode variable changes from 0 to 1  */

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
    
  
/* running parameters */   

strm<<it[0]<<","<<realpart(Enormainf)<<","<<realpart(DDelta)<<","<<N0<<","<<NCb<<","<<NIcutoff<<","<<lc<<","<<DH<<","<<workingprecision<<","
<<precDelta<<","<<precssf<<","<<maxiter<<","<<Evnorma<<"},{";


/*Q shift deviations */

strm<<Qalfadevout[0*2+0];

/*Deviations of Q/Qtilde ratios from the gluing constant */

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

delete [] deltaP;
delete [] deltaPt;

delete [] errorcode; 
 
delete [] Ev;

delete [] Qalfadev;
delete [] Qalfadevout; 


for(a=0;a<=3;a++){ delete [] c[a];                  
		   delete [] SigmaSub[a];		   
}/*a*/



delete [] m2ij;
delete [] ip2n;
delete [] m4k;
delete [] m1p4k;





delete [] t41;
delete [] t2;
delete [] t1;
delete [] y1;
delete [] y2;
delete [] y31;
delete [] s1;
delete [] s31;

for(i=0;i<=3;i++){
delete [] alfaais[i];
delete [] T5[i];
delete [] T3[i];
delete [] Y32[i];
delete [] S1n[i];
delete [] S32[i];  }




delete [] TXm;
delete [] TXmi;

delete [] PaT;
delete [] PtaT;


delete [] PujT;



delete [] TuMiMaNI;
delete [] TuMiMaNIshift;

delete [] Qlower;

delete [] Qtlower;


delete [] V;
delete [] V0;

delete [] E;
delete [] E0;

delete [] PhiV;

delete [] DEplus;


delete [] gvmin;
delete [] gvmax;

for(i=0;i<=3;i++){ delete [] TuANIn[i]; delete [] TuANInshift[i]; }

delete [] niminmax;

delete [] NIT;

/*****/

delete [] CU;
delete [] CT;
delete [] suA;

delete [] Xrealauxhalf;



return 0;
}
/******************************************************************************************/
/******************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
