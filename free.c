#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "circuit.h"
#include "hash.h"
#include "mna.h"
#include "options.h"
#include "mna_sparse.h"
#include "csparse.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

//FREE MEMORY
void freeMemory(){


	//FREE lists
	VoltageSourceT *nodeV, *freeNodeV;
	CurrentSourceT *nodeI, *freeNodeI;
	ResistorT *nodeR, *freeNodeR;
	CapacitorT *nodeC, *freeNodeC;
	InductorT *nodeL, *freeNodeL;
	DiodeT *nodeD, *freeNodeD;
	MOST *nodeM, *freeNodeM;
	BJTT *nodeQ, *freeNodeQ;
	struct PWL *Pwl,*freePwl;

	nodeV=rootV;
	while(nodeV!=NULL){
		freeNodeV=nodeV;
		nodeV=nodeV->next;
		if(freeNodeV->exp!=NULL)free(freeNodeV->exp);
		if(freeNodeV->sin!=NULL)free(freeNodeV->sin);
		if(freeNodeV->pulse!=NULL)free(freeNodeV->pulse);
		Pwl=freeNodeV->pwl;
		while(Pwl!=NULL){
		  freePwl=Pwl;
		  Pwl=Pwl->next;
		  free(freePwl);
		}
		free(freeNodeV->name);
		free(freeNodeV->pos_term);
		free(freeNodeV->neg_term);
		free(freeNodeV);
	}

	nodeI=rootI;
	while(nodeI!=NULL){
		freeNodeI=nodeI;
		nodeI=nodeI->next;
		if(freeNodeI->exp!=NULL)free(freeNodeI->exp);
		if(freeNodeI->sin!=NULL)free(freeNodeI->sin);
		if(freeNodeI->pulse!=NULL)free(freeNodeI->pulse);
		Pwl=freeNodeI->pwl;
		while(Pwl!=NULL){
		  freePwl=Pwl;
		  Pwl=Pwl->next;
		  free(freePwl);
		}
		free(freeNodeI->name);
		free(freeNodeI->pos_term);
		free(freeNodeI->neg_term);
		free(freeNodeI);
	}

	nodeR=rootR;
	while(nodeR!=NULL){
		freeNodeR=nodeR;
		nodeR=nodeR->next;
		free(freeNodeR->name);
		free(freeNodeR->pos_term);
		free(freeNodeR->neg_term);
		free(freeNodeR);
	}

	nodeC=rootC;
	while(nodeC!=NULL){
		freeNodeC=nodeC;
		nodeC=nodeC->next;
		free(freeNodeC->name);
		free(freeNodeC->pos_term);
		free(freeNodeC->neg_term);
		free(freeNodeC);
	}

	nodeL=rootL;
	while(nodeL!=NULL){
		freeNodeL=nodeL;
		nodeL=nodeL->next;
		free(freeNodeL->name);
		free(freeNodeL->pos_term);
		free(freeNodeL->neg_term);
		free(freeNodeL);
	}

	nodeD=rootD;
	while(nodeD!=NULL){
		freeNodeD=nodeD;
		nodeD=nodeD->next;
		free(freeNodeD->name);
		free(freeNodeD->pos_term);
		free(freeNodeD->neg_term);
		if(freeNodeD->model_name!=NULL)free(freeNodeD->model_name);
		free(freeNodeD);
	}

	nodeM=rootM;
	while(nodeM!=NULL){
		freeNodeM=nodeM;
		nodeM=nodeM->next;
		free(freeNodeM->name);
		free(freeNodeM->drain);
		free(freeNodeM->gate);
		free(freeNodeM->source);
		free(freeNodeM->body);
		if(freeNodeM->model_name!=NULL)free(freeNodeM->model_name);
		free(freeNodeM);
	}

	nodeQ=rootQ;
	while(nodeQ!=NULL){
		freeNodeQ=nodeQ;
		nodeQ=nodeQ->next;
		free(freeNodeQ->name);
		free(freeNodeQ->collector);
		free(freeNodeQ->base);
		free(freeNodeQ->emitter);
		if(freeNodeQ->model_name!=NULL)free(freeNodeQ->model_name);
		free(freeNodeQ);
	};
	
	
	//free MNA
	gsl_permutation_free(p);
	gsl_vector_free(B);
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_matrix_free(C);
	

	//free csparse matrices 
	cs_spfree(C_sparse);
//	cs_sfree(S);
//	cs_nfree(N);
	cs_spfree(E_sparse);
	
	//free double arrays
	free(B_sparse);
	free(x_sparse);
	
	//free hash tables
	
	int i;
	entry_t *temp;
	entry_t *current = NULL;
	for(i=0;i<hashtable->size;i++){
		current = hashtable->table[i];
		while(current != NULL) {
			temp=current;
			current = current->next;
			free(temp);
		}
	}
	free(hashtable->table);
	free(hashtable);

	
	
	//free all variables
	free(plot_nodes);
	for(i=0;i<plot_size;i++){
	  free(plot_names[i]);
	}

}
