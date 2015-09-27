#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "circuit.h"
#include "hash.h"
#include "mna.h"
#include "mna_sparse.h"
#include "options.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "transient.h"
#include "csparse.h"
#define pi 3.14159265359

void transient(){

	double TRAN_time_step;
	gsl_vector *temp_e;
	int i,j;
	double *e_sparse;
	double *e_prev_sparse;
	gsl_vector *e;
	gsl_vector *e_prev;

	if(SPARSE){
	  e_sparse=calloc(sizeB,sizeof(double));
	  e_prev_sparse=calloc(sizeB,sizeof(double));
	  initPlotFiles("Sparse");
	}else{
	  e = gsl_vector_calloc(sizeB);
	  e_prev = gsl_vector_calloc(sizeB);
	  initPlotFiles("Dense");
	}

	if(SPARSE){
	  for(i = 0; i < sizeB; ++i){
	    e_prev_sparse[i] = e_sparse[i] = B_sparse[i];
	  }
	}else{
	  gsl_blas_dcopy(B,e);
	  gsl_blas_dcopy(B,e_prev);
	}

	if(SPARSE){
	  
	  //TRAPEZODIAL
	  if(METHOD==0){
	    C_sparse = cs_add(C_sparse, E_sparse,1, 2/time_step);	//E_sparse=G-2/h*E=G-2/H*E
	    E_sparse =cs_add(C_sparse,E_sparse,1,-4/time_step);
	  }else{
	    E_sparse =cs_add(E_sparse,E_sparse,0,1/time_step); 	//E_sparse=1/h*E=1/h*C
	    C_sparse = cs_add (C_sparse, E_sparse, 1, 1);
	  }
	}else{

	  //TRAPEZODIAL
	  if(METHOD==0){
	    for(i = 0; i < sizeA; i++)
	      for(j = 0; j < sizeA; j++){
		gsl_matrix_set(A,i,j,gsl_matrix_get(A,i,j) + (2/time_step)*gsl_matrix_get(C,i,j));					
		gsl_matrix_set(C,i,j,gsl_matrix_get(A,i,j) - (4/time_step)*gsl_matrix_get(C,i,j));	//C=G-(2/h)*C
	      }
	  }else{
	    for(i = 0; i < sizeA; i++)
	      for(j = 0; j < sizeA; j++){
		gsl_matrix_set(C,i,j,(1/time_step)*gsl_matrix_get(C,i,j));		//C=1/h*C
		gsl_matrix_set(A,i,j,gsl_matrix_get(A,i,j) + gsl_matrix_get(C,i,j));	//A=G+1/h*C=A+1/h*C	
	      }
	  }
	}

	for(TRAN_time_step=0; TRAN_time_step<=end_time; TRAN_time_step+=time_step){
	  //printf("time:%g\n",TRAN_time_step);
	  if(SPARSE){
	    for(i = 0; i < sizeB; ++i){
	      e_prev_sparse[i] = e_sparse[i];
	    }
	  }else{
	    gsl_blas_dcopy(e,e_prev);
	  }

	  //COMPUTE NEW e(tk)
	  temp_e=compute_E(TRAN_time_step);

	  if(SPARSE){
	    free(e_sparse);
	    e_sparse=gslvector2double(temp_e,sizeB);
	  }else{
	    gsl_blas_dcopy(temp_e,e);
	  }

	  gsl_vector_free(temp_e);
	  
	  if(METHOD == 1){
	    //printf("Backward Euler\n");		
	    if(SPARSE)
	      TRAN_backward_euler_sparse(x_sparse,E_sparse,e_sparse);
	    else
	      TRAN_backward_euler(x, C,e);
	  }else{
	    //printf("Trapezoidial\n");
	    if(SPARSE){
	      TRAN_trapezoidial_sparse(x_sparse,E_sparse,e_sparse,e_prev_sparse);
	    }else{
	      TRAN_trapezoidial(x, C, e, e_prev);
	    }
	  }
	  
	  /*EPILUSH SUSTHMATOS*/		///################################
	  if(SPARSE==0){
	    if(SPD==0)
	      solve(TRAN_time_step);
	    else
	      solve_spd(TRAN_time_step);
	  }else{
	    if(SPD==0)
	      solveSparse(TRAN_time_step);
	    else	
	      solve_spdSparse(TRAN_time_step); 
	  }					///################################

	}
	if(SPARSE){
	  free(e_sparse);
	  free(e_prev_sparse);	
	}
	else{
	  gsl_vector_free(e);
	  gsl_vector_free(e_prev);
	}

}

void TRAN_backward_euler_sparse(double *x_sparse, cs *E_sparse,double *e){

	//printf("SPARSE\n");
	int i;
	double temp[sizeA];
	
	memset(temp, 0, sizeA*sizeof(double));
	
	cs_gaxpy (E_sparse, x_sparse, temp);
	
	for(i=0;i<sizeB;i++){
		B_sparse[i]=e[i]+temp[i];
	}

}

void TRAN_backward_euler(gsl_vector *x, gsl_matrix *C,gsl_vector *e){
	
	//printf("DENSE\n");

	gsl_vector *temp;
	
	temp=gsl_vector_calloc(sizeB);

	gsl_blas_dgemv(CblasNoTrans,1.0,C,x,0.0,temp);	//C=C*X=(1/h)*C*x(tk-1)
	gsl_vector_add(temp,e);				//C=e(tk) + (1/h)*C*x(tk-1)
	gsl_blas_dcopy(temp,B);
	

}

void TRAN_trapezoidial_sparse(double *x_sparse, cs *E_sparse,double *e,double *e_prev){

	//printf("SPARSE\n");

	int i;	
	double temp[sizeA];
	memset(temp, 0, sizeA*sizeof(double));
	cs_gaxpy (E_sparse, x_sparse, temp);
	
	for(i=0;i<sizeB;i++){
		B_sparse[i]=e[i]+e_prev[i]-temp[i];
	}


}

void TRAN_trapezoidial(gsl_vector *x, gsl_matrix *C,gsl_vector *e,gsl_vector *e_prev){

	//printf("DENSE\n");

	gsl_vector *temp_e;
	temp_e=gsl_vector_calloc(sizeB);
	gsl_vector *temp;
	temp=gsl_vector_calloc(sizeB);
		
	gsl_blas_dgemv(CblasNoTrans,1.0,C,x,0.0,temp);	//C=C*x=(G-(2/h)*C)*x(tk-1)
	gsl_vector_add(temp_e,e_prev);
	gsl_vector_add(temp_e,e);
	gsl_vector_sub(temp_e,temp);
	gsl_blas_dcopy(temp_e,B);	
	

}

gsl_vector *compute_E(double t){
  VoltageSourceT *runnerV;
  CurrentSourceT *runnerI;
  struct PWL *currPwl;
  struct PWL *prevPwl;
  gsl_vector *E;
  int i,j,k=0,b=1;
  double value,t_temp;


  E = gsl_vector_calloc(sizeA);

  runnerI = rootI;
  while(runnerI!=NULL){
    i = atoi(ht_get(hashtable,runnerI->pos_term));
    j = atoi(ht_get(hashtable,runnerI->neg_term));
    
    if(strcmp(runnerI->transient_spec,"EXP")==0){
      if(t<=runnerI->exp->td1){
	value = runnerI->exp->i1;
      }else if(t<=runnerI->exp->td2){
	value = runnerI->exp->i1 + (runnerI->exp->i2 - runnerI->exp->i1)*(1 - exp(-(t-runnerI->exp->td1)/runnerI->exp->tc1));
      }else if(t<=end_time){
	value = runnerI->exp->i1 + (runnerI->exp->i2 - runnerI->exp->i1)*(exp(-(t-runnerI->exp->td2)/runnerI->exp->tc2) - exp(-(t-runnerI->exp->td1)/runnerI->exp->tc1));
      }
    }else if(strcmp(runnerI->transient_spec,"SIN")==0){
      if(t<=runnerI->sin->td){
	value = runnerI->sin->i1 + runnerI->sin->ia*sin(2*pi*runnerI->sin->ph/360);
      }else if(t<=end_time){
	value = runnerI->sin->i1 + runnerI->sin->ia*sin(2*pi*runnerI->sin->fr*(t-runnerI->sin->td)+2*pi*runnerI->sin->ph/360);
      }
    }else if(strcmp(runnerI->transient_spec,"PULSE")==0){
      k=(t-runnerI->pulse->td)/runnerI->pulse->per;
      if(k<0)k=0;
      
      t_temp=t-k*runnerI->pulse->per;
      
      if(t_temp<=runnerI->pulse->td){
	value = runnerI->pulse->i1;
      }else if(t_temp<=runnerI->pulse->td+runnerI->pulse->tr){
	//linear i1->i2
	value = linear_value(runnerI->pulse->td, runnerI->pulse->i1, runnerI->pulse->td+runnerI->pulse->tr, runnerI->pulse->i2, t_temp);
      }else if(t_temp<=runnerI->pulse->td+runnerI->pulse->tr+runnerI->pulse->pw){
	value = runnerI->pulse->i2;
      }else if(t_temp<=runnerI->pulse->td+runnerI->pulse->tr+runnerI->pulse->pw+runnerI->pulse->tf){
	//linear i2->i1
	value = linear_value(runnerI->pulse->td+runnerI->pulse->tr+runnerI->pulse->pw, runnerI->pulse->i2, runnerI->pulse->td+runnerI->pulse->tr+runnerI->pulse->pw+runnerI->pulse->tf, runnerI->pulse->i2, t_temp);
      }else if(t_temp<=runnerI->pulse->td+runnerI->pulse->per){
	value = runnerI->pulse->i1;
      }
    }else if(strcmp(runnerI->transient_spec,"PWL")==0){
      currPwl=runnerI->pwl;
      if(t<=currPwl->t){
	value = currPwl->i;
      }else{
	prevPwl=currPwl;
	currPwl=currPwl->next;
	while(currPwl!=NULL){
	  if(t<=currPwl->t){
	    //linear
	    value = linear_value(prevPwl->t, prevPwl->i, currPwl->t, currPwl->i, t);
	    break;
	  }
	  prevPwl=currPwl;
	  currPwl=currPwl->next;
	}
	if(currPwl==NULL){
	  value = prevPwl->i;
	}
      }
      
    }else{
      value = runnerI->value;
    }
    if(i!=0){
      gsl_vector_set (E, i-1, gsl_vector_get(E,i-1) - value);		
    }
    if(j!=0){

      gsl_vector_set (E, j-1, gsl_vector_get(E,j-1) + value);
    }
    runnerI=runnerI->next;
  }

  runnerV = rootV;
  while(runnerV!=NULL){
    if(strcmp(runnerV->transient_spec,"EXP")==0){
      if(t<=runnerV->exp->td1){
	value = runnerV->exp->i1;
      }else if(t<=runnerV->exp->td2){
	value = runnerV->exp->i1 + (runnerV->exp->i2 - runnerV->exp->i1)*(1 - exp(-(t-runnerV->exp->td1)/runnerV->exp->tc1));
      }else if(t<=end_time){
	value = runnerV->exp->i1 + (runnerV->exp->i2 - runnerV->exp->i1)*(exp(-(t-runnerV->exp->td2)/runnerV->exp->tc2) - exp(-(t-runnerV->exp->td1)/runnerV->exp->tc1));
      }
    }else if(strcmp(runnerV->transient_spec,"SIN")==0){
      if(t<=runnerV->sin->td){
	value = runnerV->sin->i1 + runnerV->sin->ia*sin(2*pi*runnerV->sin->ph/360);
      }else if(t<=end_time){
	value = runnerV->sin->i1 + runnerV->sin->ia*sin(2*pi*runnerV->sin->fr*(t-runnerV->sin->td)+2*pi*runnerV->sin->ph/360);
      }
    }else if(strcmp(runnerV->transient_spec,"PULSE")==0){
      k=(t-runnerV->pulse->td)/runnerV->pulse->per;
      if(k<0)k=0;
      
      t_temp=t-k*runnerV->pulse->per;
      
      if(t_temp<=runnerV->pulse->td){
	value = runnerV->pulse->i1;
      }else if(t_temp<=runnerV->pulse->td+runnerV->pulse->tr){
	//linear i1->i2
	value = linear_value(runnerV->pulse->td, runnerV->pulse->i1, runnerV->pulse->td+runnerV->pulse->tr, runnerV->pulse->i2, t_temp);
      }else if(t_temp<=runnerV->pulse->td+runnerV->pulse->tr+runnerV->pulse->pw){
	value = runnerV->pulse->i2;
      }else if(t_temp<=runnerV->pulse->td+runnerV->pulse->tr+runnerV->pulse->pw+runnerV->pulse->tf){
	//linear i2->i1
	value = linear_value(runnerV->pulse->td+runnerV->pulse->tr+runnerV->pulse->pw, runnerV->pulse->i2, runnerV->pulse->td+runnerV->pulse->tr+runnerV->pulse->pw+runnerV->pulse->tf, runnerV->pulse->i2, t_temp);
      }else if(t_temp<=runnerV->pulse->td+runnerV->pulse->per){
	value = runnerV->pulse->i1;
      }
    }else if(strcmp(runnerV->transient_spec,"PWL")==0){
      currPwl=runnerV->pwl;
      if(t<=currPwl->t){
	value = currPwl->i;
      }else{
	prevPwl=currPwl;
	currPwl=currPwl->next;
	while(currPwl!=NULL){
	  if(t<=currPwl->t){
	    //linear
	    value = linear_value(prevPwl->t, prevPwl->i, currPwl->t, currPwl->i, t);
	    break;
	  }
	  prevPwl=currPwl;
	  currPwl=currPwl->next;
	}
	if(currPwl==NULL){
	  value = prevPwl->i;
	}
      }
      
    }else{
      value = runnerV->value;
    }
    gsl_vector_set (E, hashNode_num-2+b, value);

    b++;
    runnerV=runnerV->next;
  }

  return E;
}

double linear_value(double x1, double y1, double x2, double y2, double x){
  double y;
  double m;
  
  m = (y2-y1)/(x2-x1);
  
  y = m*(x-x1) + y1;
  
  return y ;
}


