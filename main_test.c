/*
 * FINAL IBM 
 * Μυρτώ Ζάββου 989
 * Πέτρος Κούλαλης 946
 * Χρήστος Κυρίτσης 935
 */



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "circuit.h"
#include "hash.h"
#include "mna.h"
#include "options.h"
#include "free.h"
#include "csparse.h"
#include "mna_sparse.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include "transient.h"

int main(int argc, char *argv[]){
  
  /*open file*/
  FILE *input_file_ptr;
  input_file_ptr=fopen( argv[1], "r" );
  printf( "**** opening file %s\n",argv[1] );
    
  if( !input_file_ptr )
  {
    printf( "Error: can't find the specified file\n " );
  }
  circuit_init();
  readCircuit(input_file_ptr);
  
  printf("SPARSE=%d\n",SPARSE);
  printf("SPD=%d\n",SPD);
  printf("ITER=%d\n",ITER);
  printf("METHOD=%d\n",METHOD);
  printf("TRAN=%d\n",TRAN);
  printf("plot=%d\n",plot);
  
  //An den yparxei geiwsi to programma termatizei
  if(ground==0)
  {
    printf("\nERROR: No ground. Exit program.\n");
    return(0);
  }

  //printList();
  //printHash();
  
  if(SPARSE==0)
  {	
    computeMNA();	
  }
  else
  {	
    computeMNASparse();	   
  }
  
  double time = -1;
  
  if(SPARSE==0)
  {
    if(SPD==0)
    {    
      solve(time);
    }
    else
    {
      solve_spd(time);
    }
  }
  else
  {
    if(SPD==0)
    {		
      solveSparse(time);	  
    }
    else
    {		
      solve_spdSparse(time);		
    }
  }
  
  if(TRAN==1) transient();
  
  freeMemory();
  return(0);

}