#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "circuit.h"
#include "hash.h"
#include "mna.h"
#include "options.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>

//Sinartisi pou diavazei ta options, diladi diavazei o,ti yparxei meta apo '.'
//sto netlist kai kathorizei ti tha ginei.
//
//Diavazei: OPTIONS, DC, TRAN, PLOT
//
//Diavazei prwta olokliri ti grammi apo to netlist kai meta vlepei se poia apo
//tis panw periptwseis eimaste opote prattei ta analoga...

void options(char *input_line){

  char input_line2[1024];
  char *readOption;
  const char delim[] = "\n ,\t()";
  VoltageSourceT *nodeV;
  CurrentSourceT *nodeI;
  int n,fid;
  char *gt_id;

  input_line++;						//olokliri i grammi meta to '.'
  strcpy(input_line2,input_line);
  readOption = strtok (input_line, delim);	//diavazei thn prwth leksi tis grammis
 
  //Diavase ta OPTIONS
  if(!(strcmp(readOption,"OPTIONS"))){	
	  while((readOption = strtok (NULL, delim))!=NULL){
	    if ((strcmp(readOption,"SPD"))==0)
	    {
	      SPD=1;
	    }
	    if (!(strcmp(readOption,"ITER")))
	    {
	      ITER=1;
	    }
	    if (!(strcmp(readOption,"METHOD=TR")))
	    {
	      //printf("METHOD=TR");
		METHOD=0;
	    }
	    if (!(strcmp(readOption,"METHOD=BE")))
	    {
		//printf("METHOD=BE");
		METHOD=1;
	    }
	    if (!(strcmp(readOption,"ITOL")))
	    {
	      if (strtok (NULL, delim) != NULL)
	      {
		itol_value=atof(strtok (NULL, delim));
	      }
	    }
	    if (!(strcmp(readOption,"SPARSE")))
	    {
	      SPARSE=1;
	    }
	  }
	  return;
  } //diavase to DC
  else if(!(strcmp(readOption,"DC")))
  {	
    readOption = strtok (NULL, delim);
    
    if(readOption==NULL)
    {
      return;     //an meta to DC exw NULL
    } 
    
    // an vrw V meta to DC,DC Sweep gia voltage source
    if(readOption[0]=='V')
    {
      nodeV=rootV;
      n=hashNode_num;
  
      while(nodeV!=NULL){
	if(!(strcmp(nodeV->name,readOption))){	  
	  sweep_source=n;
	  break;
	}
	n++;
	nodeV=nodeV->next;	
      }
      start_value = atof(strtok (NULL, delim));
      end_value = atof(strtok (NULL, delim));
      sweep_step = atof(strtok (NULL, delim));   
      dc_sweep=1;
      //printf("\n\n");
      //printf("DC Sweep \n\n");
      //printf("Voltage Source:\t\t%s\nSource Array B: %d\n",readOption,sweep_source);
      //printf("Start Value:\t\t%.3lf\n",start_value);
      //printf("End Value:\t\t%.3lf\n",end_value);
      //printf("Sweep Step:\t\t%.3lf\n",sweep_step);
      //printf("\n\n");
      return;
    } //an vrw I meta to DC, tote DC Sweep gia current source
    else if(readOption[0]=='I')
    {

      nodeI=rootI;

      while(nodeI!=NULL){
	if(!(strcmp(nodeI->name,(readOption)))){
	  sweep_source=-1;						//DC Sweep gia pigi reumatos
	  sweep_posNode=atoi(ht_get(hashtable, nodeI->pos_term));	//pos_term tis pigis reumatos
	  sweep_negNode=atoi(ht_get(hashtable, nodeI->neg_term));	//neg_term tis pigis reumatos
	  sweep_value_I=nodeI->value;
	  break;
	}
	nodeI=nodeI->next;	
      }
      start_value = atof(strtok (NULL, delim));
      end_value = atof(strtok (NULL, delim));
      sweep_step = atof(strtok (NULL, delim));
      dc_sweep=1;
      //printf("\n\n");
      //printf("DC SWEEP\n\n");
      //printf("Current Source:\t%s\nPositive node: %d\nNegative node: %d\n", readOption, sweep_posNode, sweep_negNode);
      //printf("Start Value:\t%.1e\n",start_value);
      //printf("End Value:\t%.1e\n",end_value);
      //printf("Sweep Step:\t%.1e\n",sweep_step);
      //printf("\n\n");
      return;
    }
    
    
  }
  else if(!(strcmp(readOption,"TRAN")))
  {
    TRAN=1;
    time_step = atof(strtok (NULL, delim));
    end_time = atof(strtok (NULL, delim));
    return;
  }
  else if(!(strcmp(readOption,"PLOT")))
  {

	  readOption = strtok (NULL, delim);
	  if((readOption==NULL) || (strcmp(readOption,"V"))!=0)
	  {
	    return;	    
	  }
	  readOption = strtok (NULL, delim);
	  plot_size=1;
	  while(readOption!=NULL){
		readOption = strtok (NULL, delim);
		if((readOption==NULL) || (strcmp(readOption,"V"))!=0)
		{
		  break;
		}
		readOption = strtok (NULL, delim);
		plot_size++;
	 }

         //printf("PLOT_SIZE=%d\n",plot_size);
	 plot_nodes= (int *)calloc(plot_size+1,sizeof(int)) ;
	 plot_names= (char **)calloc(plot_size,sizeof(char *));
	 
	 readOption = strtok (input_line2, delim);
	 readOption = strtok (NULL, delim);
	 if((readOption==NULL) || (strcmp(readOption,"V"))!=0){return;}
	 readOption = strtok (NULL, delim);
	 n=0;
	 
	 while(readOption!=NULL){
	      gt_id=ht_get(hashtable, readOption);
	      if(gt_id==NULL){printf("PLOT NODE DOESN'T EXIST\nPROGRAM STOPPED!!!\n\n");exit(1);}
	      fid=atoi (gt_id);
	      plot_nodes[n]=fid;
	      plot_names[n]=strdup(readOption);
	      readOption = strtok (NULL, delim);
	      if((readOption==NULL) || (strcmp(readOption,"V"))!=0){break;}
	      readOption = strtok (NULL, delim);
	      n++;
	 }
	 plot=1;
	 return;
  }
}

