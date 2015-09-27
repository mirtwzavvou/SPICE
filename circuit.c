#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "circuit.h"
#include "hash.h"
#include "mna.h"
#include "options.h"
#include "mna_sparse.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>

void circuit_init(){

	rootV=NULL;
	rootI=NULL;
	rootR=NULL;
	rootC=NULL;
	rootL=NULL;
	rootD=NULL;
	rootM=NULL;
	rootQ=NULL;
	ground=0;
	SPD=0;
	ITER=0;
	itol_value=1e-3;
	SPARSE=0;
	TRAN=0;
	METHOD=0;
	plot=0;
	m2=0;
	sizeA_sparse=0;
	sizeD_sparse=0;
	hashNode_num=1;
	hashtable = ht_create( 10000 );	//Create Hash table,oso ligotero tosa ligotera collision 
	
	dc_sweep=0;
	sweep_value_I=0;

} 
 
void readCircuit(FILE *input_ptr) {
  
  char *input_line;					//input line we read
  input_line=( char * )malloc( sizeof( char ) * 1024 ); //allocate memory for variable input_line
  int i=0;
  
  while( fgets(input_line, sizeof(input_line)*1024, input_ptr) != NULL) //read file line by line
  {


    //elegxos gia keni grammi
    if(input_line[i]=='\n')
    {
      continue;
    }
    
    //elegxos gia keno i tab stin arxi tis grammis
    if(input_line[i]==' ' || input_line[i]=='\t')
    {
      //sunexizei mexri na vrei != apo keno i tab
      while(input_line[i]== ' ' || input_line[i]=='\t')
      {
	i++;
      } 
    }
    
    //diavazei ta sxolia
    if(input_line[i] == '*')
    {
      continue;
    }
    
    if( (input_line[i]== 'v') || (input_line[i] == 'V'))
    {
      insertV(input_line);
      m2++;
    }
    
    if( (input_line[i]== 'r') || (input_line[i] == 'R'))
    {
      insertR(input_line);
    }
    
    if( (input_line[i]== 'i') || (input_line[i] == 'I'))
    {
      insertI(input_line);
    }
    
    if( (input_line[i]== 'c') || (input_line[i] == 'C'))
    {
      insertC(input_line);
    }
    
    if( (input_line[i]== 'l') || (input_line[i] == 'L'))
    {
      insertL(input_line);
      m2++;
    }
    
    if( (input_line[i]== 'd') || (input_line[i] == 'D'))
    {
      insertD(input_line);
    }
    
    if( (input_line[i]== 'm') || (input_line[i] == 'M'))
    {
      insertM(input_line);
    }
    
    if( (input_line[i]== 'q') || (input_line[i] == 'Q'))
    {
      insertQ(input_line);
    }
    
    if(input_line[i] == '.') 
    {
      options(input_line);
    }
    
    i=0;
    
  }
  
}
  
void insertV(char *line){
  
  char *readElement;
  const char delim[] = "\n ,\t()=";
  struct PWL* new2;
  struct PWL* temp;
  
  VoltageSourceT *newElement;

  newElement = (VoltageSourceT*) malloc(sizeof(VoltageSourceT));
  
  readElement=strtok (line, delim);

  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou
  
  readElement = strtok (NULL, delim);
  
  newElement->pos_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->pos_term,readElement);	//opou pos_term einai to positive
  
  AddToHashtable(readElement);
  
  readElement = strtok (NULL, delim);
  
  newElement->neg_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->neg_term,readElement);	//opou neg_term einai to negative
  
  AddToHashtable(readElement);
  
  readElement = strtok (NULL, delim);
  
  //opou value einai to value
  newElement->value=atof(readElement);
  
  readElement = strtok (NULL, delim);
   
  if (readElement!=NULL) {
    if(strcmp(readElement,"EXP")==0 || strcmp(readElement, "exp") ==0 )  
    {
	strcpy(newElement->transient_spec,"EXP");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	
	newElement->exp=(struct EXP*)malloc(sizeof(struct EXP));
	newElement->exp->i1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->i1: %lf\n",newElement->exp->i1);
	newElement->exp->i2 = atof(strtok (NULL,delim));
//	  printf("\n\t\t\t\texp->i2: %lf\n",newElement->exp->i2);
	newElement->exp->td1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->td1: %lf\n",newElement->exp->td1);
	newElement->exp->tc1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->tc1: %lf\n",newElement->exp->tc1);
	newElement->exp->td2 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->td2: %lf\n",newElement->exp->td2);
	newElement->exp->tc2 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->tc2: %lf\n",newElement->exp->tc2);
    }
    else if(strcmp(readElement,"SIN")==0||strcmp(readElement,"sin")==0)
    {
	strcpy(newElement->transient_spec,"SIN");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	newElement->sin=(struct SIN*)malloc(sizeof(struct SIN));
	newElement->sin->i1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->i1: %lf\n",newElement->sin->i1);
	newElement->sin->ia = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->ia: %lf\n",newElement->sin->ia);
	newElement->sin->fr = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->fr: %lf\n",newElement->sin->fr);
	newElement->sin->td = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->td: %lf\n",newElement->sin->td);
	newElement->sin->df = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->df: %lf\n",newElement->sin->df);
	newElement->sin->ph = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->ph: %lf\n",newElement->sin->ph);
    }
    else if(strcmp(readElement,"PULSE")==0||strcmp(readElement,"pulse")==0)
    {
	strcpy(newElement->transient_spec,"PULSE");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	newElement->pulse=(struct PULSE*)malloc(sizeof(struct PULSE));
	newElement->pulse->i1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->i1: %lf\n",newElement->pulse->i1);
	newElement->pulse->i2 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->i2: %lf\n",newElement->pulse->i2);
	newElement->pulse->td = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->td: %lf\n",newElement->pulse->td);
	newElement->pulse->tr = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->tr: %lf\n",newElement->pulse->tr);
	newElement->pulse->tf = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->tf: %lf\n",newElement->pulse->tf);
	newElement->pulse->pw = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->pw: %lf\n",newElement->pulse->pw);
	newElement->pulse->per = atof(strtok (NULL, " ,\t()"));
//	  printf("\n\t\t\t\tpulse->per: %lf\n",newElement->pulse->per);
    }
    else if(strcmp(readElement,"PWL")==0||strcmp(readElement,"pwl")==0)
    {
	strcpy(newElement->transient_spec,"PWL");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	newElement->pwl=(struct PWL*)malloc(sizeof(struct PWL));
	newElement->pwl->t = atof(strtok (NULL, " ,\t()\n"));
//	  printf("\n\t\t\t\tpwl->t: %lf\n",newElement->pwl->t);
	newElement->pwl->i = atof(strtok (NULL, " ,\t()\n"));
//	  printf("\n\t\t\t\tpwl->i: %lf\n",newElement->pwl->i);
	newElement->pwl->next=NULL;
	temp=newElement->pwl;
	readElement = strtok (NULL, " ,\t()\n");
	
	while(readElement!=NULL){
	  new2=(struct PWL*)malloc(sizeof(struct PWL));
	  temp->next=new2;
	  new2->t=atof(readElement);
//	  printf("\n\t\t\t\tpwl->t: %lf\n",new2->t);
	  new2->i=atof(strtok (NULL, " ,\t()\r\n"));
//	  printf("\n\t\t\t\tpwl->i: %lf\n",new2->i);
	  readElement = strtok (NULL, " ,\t()\r\n");
//	  printf("\n\t\t\t\treadElement: %s\n",readElement);
	  temp=new2;
	}
    }
  }
  
  newElement->next = rootV;
  rootV = newElement;
  
  
}


void insertR(char *line){
  
  char *readElement;
  const char delim[] = "\n ,\t()=";
  
  ResistorT *newElement;

  newElement = (ResistorT*) malloc(sizeof(ResistorT));
  
  readElement=strtok (line, delim);
      
  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou
  
  readElement = strtok (NULL, delim);
  
  newElement->pos_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->pos_term,readElement);	//opou pos_term einai to positive

  AddToHashtable(readElement);
  
  readElement = strtok (NULL, delim);
  
  newElement->neg_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->neg_term,readElement);	//opou neg_term einai to negative
  
  AddToHashtable(readElement);
  
  readElement = strtok (NULL, delim);
  
  //opou value einai to value
  newElement->value=atof(readElement);
  
  //tsekaroume an uparxei toulaxiston enas kombos 0
  if( (strlen(newElement->pos_term) == 1 && newElement->pos_term[0] == '0') || (strlen(newElement->neg_term) == 1 && newElement->neg_term[0] == '0'))
  {
    ground=1;
  } 

  newElement->next = rootR;
  rootR = newElement;
  
}


void insertI(char *line){

  char *readElement;
  const char delim[] = "\n ,\t()=";\
  struct PWL* new2;
  struct PWL* temp;
  
  CurrentSourceT *newElement;

  newElement = (CurrentSourceT*) malloc(sizeof(CurrentSourceT));

  readElement=strtok (line, delim);
  
  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou
 
  readElement = strtok (NULL, delim);
  
  newElement->pos_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->pos_term,readElement);	//opou pos_term einai to positive

  AddToHashtable(readElement);
 
  readElement = strtok (NULL, delim);  
  
  newElement->neg_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->neg_term,readElement);	//opou neg_term einai to negative
  
  AddToHashtable(readElement);
  
  readElement = strtok (NULL, delim);  
  
  //opou value einai to value
  newElement->value=atof(readElement);
  
  //tsekaroume an uparxei toulaxiston enas kombos 0
  if( (strlen(newElement->pos_term) == 1 && newElement->pos_term[0] == '0') || (strlen(newElement->neg_term) == 1 && newElement->neg_term[0] == '0'))
  {
    ground=1;
  } 
  
 readElement = strtok (NULL, delim);
   
  if (readElement!=NULL) {
    if(strcmp(readElement,"EXP")==0 || strcmp(readElement, "exp") ==0 )  
    {
	strcpy(newElement->transient_spec,"EXP");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	
	newElement->exp=(struct EXP*)malloc(sizeof(struct EXP));
	newElement->exp->i1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->i1: %lf\n",newElement->exp->i1);
	newElement->exp->i2 = atof(strtok (NULL,delim));
//	  printf("\n\t\t\t\texp->i2: %lf\n",newElement->exp->i2);
	newElement->exp->td1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->td1: %lf\n",newElement->exp->td1);
	newElement->exp->tc1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->tc1: %lf\n",newElement->exp->tc1);
	newElement->exp->td2 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->td2: %lf\n",newElement->exp->td2);
	newElement->exp->tc2 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\texp->tc2: %lf\n",newElement->exp->tc2);
    }
    else if(strcmp(readElement,"SIN")==0||strcmp(readElement,"sin")==0)
    {
	strcpy(newElement->transient_spec,"SIN");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	newElement->sin=(struct SIN*)malloc(sizeof(struct SIN));
	newElement->sin->i1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->i1: %lf\n",newElement->sin->i1);
	newElement->sin->ia = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->ia: %lf\n",newElement->sin->ia);
	newElement->sin->fr = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->fr: %lf\n",newElement->sin->fr);
	newElement->sin->td = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->td: %lf\n",newElement->sin->td);
	newElement->sin->df = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->df: %lf\n",newElement->sin->df);
	newElement->sin->ph = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tsin->ph: %lf\n",newElement->sin->ph);
    }
    else if(strcmp(readElement,"PULSE")==0||strcmp(readElement,"pulse")==0)
    {
	strcpy(newElement->transient_spec,"PULSE");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	newElement->pulse=(struct PULSE*)malloc(sizeof(struct PULSE));
	newElement->pulse->i1 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->i1: %lf\n",newElement->pulse->i1);
	newElement->pulse->i2 = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->i2: %lf\n",newElement->pulse->i2);
	newElement->pulse->td = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->td: %lf\n",newElement->pulse->td);
	newElement->pulse->tr = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->tr: %lf\n",newElement->pulse->tr);
	newElement->pulse->tf = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->tf: %lf\n",newElement->pulse->tf);
	newElement->pulse->pw = atof(strtok (NULL, delim));
//	  printf("\n\t\t\t\tpulse->pw: %lf\n",newElement->pulse->pw);
	newElement->pulse->per = atof(strtok (NULL, " ,\t()"));
//	  printf("\n\t\t\t\tpulse->per: %lf\n",newElement->pulse->per);
    }
    else if(strcmp(readElement,"PWL")==0||strcmp(readElement,"pwl")==0)
    {
	strcpy(newElement->transient_spec,"PWL");
//	  printf("\n\t\t\t\tnewElement->transient_spec: %s\n",newElement->transient_spec);
	newElement->pwl=(struct PWL*)malloc(sizeof(struct PWL));
	newElement->pwl->t = atof(strtok (NULL, " ,\t()\n"));
//	  printf("\n\t\t\t\tpwl->t: %lf\n",newElement->pwl->t);
	newElement->pwl->i = atof(strtok (NULL, " ,\t()\n"));
//	  printf("\n\t\t\t\tpwl->i: %lf\n",newElement->pwl->i);
	newElement->pwl->next=NULL;
	temp=newElement->pwl;
	readElement = strtok (NULL, " ,\t()\n");
	
	while(readElement!=NULL){
	  new2=(struct PWL*)malloc(sizeof(struct PWL));
	  temp->next=new2;
	  new2->t=atof(readElement);
//	  printf("\n\t\t\t\tpwl->t: %lf\n",new2->t);
	  new2->i=atof(strtok (NULL, " ,\t()\r\n"));
//	  printf("\n\t\t\t\tpwl->i: %lf\n",new2->i);
	  readElement = strtok (NULL, " ,\t()\r\n");
//	  printf("\n\t\t\t\treadElement: %s\n",readElement);
	  temp=new2;
	}
    }
  }  

  newElement->next = rootI;
  rootI = newElement;

}

void insertC(char *line){

  char *readElement;
  const char delim[] = "\n ,\t()=";
  
  CapacitorT *newElement;

  newElement = (CapacitorT*) malloc(sizeof(CapacitorT));
  
  readElement=strtok (line, delim);
  
  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou
  
  readElement = strtok (NULL, delim);
  
  newElement->pos_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->pos_term,readElement);	//opou pos_term einai to positive

  AddToHashtable(readElement);

  readElement = strtok (NULL, delim);  
  
  newElement->neg_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->neg_term,readElement);	//opou neg_term einai to negative
  
  AddToHashtable(readElement);
  
  readElement = strtok (NULL, delim);  
  
  //opou value einai to value
  newElement->value=atof(readElement);
  
  //tsekaroume an uparxei toulaxiston enas kombos 0
  if( (strlen(newElement->pos_term) == 1 && newElement->pos_term[0] == '0') || (strlen(newElement->neg_term) == 1 && newElement->neg_term[0] == '0'))
  {
    ground=1;
  } 

  newElement->next = rootC;
  rootC = newElement;  
  
}

void insertL(char *line){
 
  char *readElement;
  const char delim[] = "\n ,\t()=";
  
  InductorT *newElement;

  newElement = (InductorT*) malloc(sizeof(InductorT));
  
  readElement=strtok (line, delim);

  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou

  readElement = strtok (NULL, delim);   
  
  newElement->pos_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->pos_term,readElement);	//opou pos_term einai to positive

  AddToHashtable(readElement);

  readElement = strtok (NULL, delim); 
  
  newElement->neg_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->neg_term,readElement);	//opou neg_term einai to negative
  
  AddToHashtable(readElement);

  readElement = strtok (NULL, delim); 
  
  //opou value einai to value
  newElement->value=atof(readElement);
  
  //tsekaroume an uparxei toulaxiston enas kombos 0
  if( (strlen(newElement->pos_term) == 1 && newElement->pos_term[0] == '0') || (strlen(newElement->neg_term) == 1 && newElement->neg_term[0] == '0'))
  {
    ground=1;
  }  
 
  newElement->next = rootL;
  rootL = newElement; 

}

void insertD(char *line){
  
  char *readElement;
  const char delim[] = "\n ,\t()=";  
  
  DiodeT *newElement;

  newElement = (DiodeT*) malloc(sizeof(DiodeT));
  
  readElement=strtok (line, delim);

  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou

  readElement = strtok (NULL, delim); 
  
  newElement->pos_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->pos_term,readElement);	//opou pos_term einai to positive

  readElement = strtok (NULL, delim); 
  
  newElement->neg_term=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->neg_term,readElement);	//opou neg_term einai to negative

  readElement = strtok (NULL, delim); 
  
  newElement->model_name=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->model_name,readElement);	//opou model_name einai to model_name
  
  readElement = strtok (NULL, delim);   
  
  //opou area einai to area
  newElement->area=atof(readElement);
  
    //tsekaroume an uparxei toulaxiston enas kombos 0
  if( (strlen(newElement->pos_term) == 1 && newElement->pos_term[0] == '0') || (strlen(newElement->neg_term) == 1 && newElement->neg_term[0] == '0'))
  {
    ground=1;
  } 
  
  newElement->next = rootD;
  rootD = newElement;  
  
}

void insertM(char *line){
  
  char *readElement;
  const char delim[] = "\n ,\t()="; 
  
  MOST *newElement;

  newElement = (MOST*) malloc(sizeof(MOST));
  
  readElement=strtok (line, delim);

  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou
  
  readElement = strtok (NULL, delim);   
  
  newElement->drain=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->drain,readElement);	//opou drain einai to drain
  
  readElement = strtok (NULL, delim);
  
  newElement->gate=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->gate,readElement);	//opou gate einai to gate
  
  readElement = strtok (NULL, delim);
  
  newElement->source=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->source,readElement);	//opou source einai to source
  
  readElement = strtok (NULL, delim);
  
  newElement->body=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->body,readElement);	//opou body einai to body
  
  readElement = strtok (NULL, delim);
  
  newElement->model_name=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->model_name,readElement);	//opou model_name einai to model_name
  
  readElement = strtok (NULL, delim);
  
  newElement->Lvalue=atof(readElement);		//opou Lvalue to length
  
  readElement = strtok (NULL, delim);
  
  newElement->Wvalue=atof(readElement);		//opou Wvalue to width
  
  newElement->next = rootM;
  rootM = newElement;  
  
}

void insertQ(char *line){
  
  char *readElement;
  const char delim[] = "\n ,\t()="; 
  
  BJTT *newElement;

  newElement = (BJTT*) malloc(sizeof(BJTT));
  
  readElement=strtok (line, delim);

  newElement->name=( char * )malloc( ( strlen(readElement) + 1 ) * sizeof( char ) );
  strcpy(newElement->name,readElement);		//name tou stoixeiou
  
  readElement = strtok (NULL, delim);
  
  newElement->collector=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->collector,readElement);	//opou collector einai to collector
  
  readElement = strtok (NULL, delim);
  
  newElement->base=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->base,readElement);	//opou base einai to base
  
  readElement = strtok (NULL, delim);
  
  newElement->emitter=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->emitter,readElement);	//opou emitter einai to emitter
  
  readElement = strtok (NULL, delim);
  
  newElement->model_name=(char *)malloc( ( strlen(readElement) + 1 ) * sizeof(char) );
  strcpy(newElement->model_name,readElement);	//opou model_name einai to model_name
  
  readElement = strtok (NULL, delim);
  
  //opou area einai to area
  newElement->area=atof(readElement);
 
  newElement->next = rootQ;
  rootQ = newElement;  
  
}

void AddToHashtable(char *string){

	char str[15];
	
	if(strcmp(string,"0")==0){
		ht_set(hashtable, string, "0");
		return;
	}

	sprintf(str, "%d", hashNode_num);
	ht_set(hashtable, string, str);
}

void printList(){

	VoltageSourceT *nodeV;
	ResistorT *nodeR;
	CurrentSourceT *nodeI;
	CapacitorT *nodeC;
	InductorT *nodeL;
	DiodeT *nodeD;
	MOST *nodeM;
	BJTT *nodeQ;

	nodeV=rootV;
	printf("\n---Voltage Sources---\n");
	while(nodeV!=NULL){
		printf("name = %s  ",nodeV->name);
		printf("pos_term = %s  ",nodeV->pos_term);
		printf("neg_term = %s  ",nodeV->neg_term);
		printf("value = %e\n",nodeV->value);
	
		nodeV=nodeV->next;	
	}

	nodeR=rootR;
	printf("\n---Resistors---\n");
	while(nodeR!=NULL){
		printf("name = %s ",nodeR->name);
		printf("pos_term = %s  ",nodeR->pos_term);
		printf("neg_term = %s  ",nodeR->neg_term);
		printf("value = %e\n",nodeR->value);
	
		nodeR=nodeR->next;	
	}
	
	nodeI=rootI;
	printf("\n---Current Sources---\n");
	while(nodeI!=NULL){
		printf("name = %s  ",nodeI->name);
		printf("pos_term = %s  ",nodeI->pos_term);
		printf("neg_term = %s  ",nodeI->neg_term);
		printf("value = %e\n",nodeI->value);
	
		nodeI=nodeI->next;	
	}

	nodeC=rootC;
	printf("\n---Capacitors---\n");
	while(nodeC!=NULL){
		printf("name = %s ",nodeC->name);
		printf("pos_term = %s  ",nodeC->pos_term);
		printf("neg_term = %s  ",nodeC->neg_term);
		printf("value = %e\n",nodeC->value);
	
		nodeC=nodeC->next;	
	}

	nodeL=rootL;
	printf("\n---Inductors---\n");
	while(nodeL!=NULL){
		printf("name = %s  ",nodeL->name);
		printf("pos_term = %s  ",nodeL->pos_term);
		printf("neg_term = %s  ",nodeL->neg_term);
		printf("value = %e\n",nodeL->value);
	
		nodeL=nodeL->next;	
	}

	nodeD=rootD;
	printf("\n---Diodes---\n");
	while(nodeD!=NULL){
		printf("name = %s  ",nodeD->name);
		printf("pos_term = %s  ",nodeD->pos_term);
		printf("neg_term = %s  ",nodeD->neg_term);
		printf("model_name = %s ",nodeD->model_name);
		printf("area = %.6lf\n",nodeD->area);
	
		nodeD=nodeD->next;	
	}

	nodeM=rootM;
	printf("\n---MOSFET Transistors---\n");
	while(nodeM!=NULL){
		printf("name = %s  ",nodeM->name);
		printf("drain = %s  ",nodeM->drain);
		printf("gate = %s  ",nodeM->gate);
		printf("source = %s  ",nodeM->source);
		printf("body = %s  ",nodeM->body);
		printf("model_name = %s ",nodeM->model_name);
		printf("lenght = %e  ",nodeM->Lvalue);
		printf("width = %e\n",nodeM->Wvalue);
		nodeM=nodeM->next;	
	}

	nodeQ=rootQ;
	printf("\n---BJT Transistors---\n");
	while(nodeQ!=NULL){
		printf("name = %s  ",nodeQ->name);
		printf("collector = %s  ",nodeQ->collector);
		printf("base = %s  ",nodeQ->base);
		printf("emitter = %s  ",nodeQ->emitter);
		printf("model_name = %s ",nodeQ->model_name);		
		printf("area = %.6lf\n",nodeQ->area);
	
		nodeQ=nodeQ->next;	
	}
	
}

void printHash(){
	int j;
	entry_t *next = NULL;
	for(j=0;j<hashtable->size;j++){
		next = hashtable->table[j];
		while(next != NULL) {
			printf("node (%s) = hashtable node (%s)\n",next->key,next->value);
			next = next->next;
		}
	}
}

//Arxikopoiisi twn arxeiwn pou tha apothikeutoun ta apotelesmata twn ploted komvwn
void initPlotFiles(char *str){
  
  int i;
  char filename[100];
  FILE *fp;
  
  for(i=0;i<plot_size;i++){	//Adeiasma twn arxeiwn
	    strcpy(filename,"./PlotFiles/");
	    strcat(filename,str);
	    strcat(filename,"-Node ");
	    strcat(filename,plot_names[i]);
	    fp = fopen(filename, "w");
	    fflush(fp);
	    fclose(fp);
	  }
}

//Apothikeusi sta arxeia apotelesmatwn
void plotFiles(char *str, double *plot_table, double current_value, char *msg){
  
  int i;
  char filename[100];
  FILE *fp;
  
  for(i=0;i<plot_size;i++){	//Gemisma twn arxeiwn
 	strcpy(filename,"./PlotFiles/");
	strcat(filename,str);
	strcat(filename,"-Node ");
 	strcat(filename,plot_names[i]);
 	fp = fopen(filename, "a");
 	if (fp == NULL) {
	  printf("Can't open output file %s!\n",filename);
	  return;
 	}
 	//fprintf(fp,"%s\n",msg);
 	if(current_value!=-1){
	  fprintf(fp,"%lf %lf\n",current_value, plot_table[plot_nodes[i]-1]);
	}else{
	  fprintf(fp,"%lf\n",plot_table[plot_nodes[i]-1]);
	}
 	fflush(fp);
 	fclose(fp);
 }
  
  
}


