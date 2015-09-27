#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef STRUCTS_H
#define STRUCTS_H

#define EPS 1e-14

//transient_spec EXP
struct EXP{
  double i1,i2,td1,tc1,td2,tc2;
};

//transient_spec SIN
struct SIN{
  double i1,ia,fr,td,df,ph;
};

//transient_spec PULSE
struct PULSE{
  double i1,i2,td,tr,tf,pw,per;
};

//transient_spec PWL
struct PWL{
  double t,i;
  struct PWL *next;
};

//Voltage Source
typedef struct VoltageSource{
  
  char *name;
  char *pos_term;
  char *neg_term;
  double value;
  char transient_spec[6];
  struct EXP *exp;
  struct SIN *sin;
  struct PULSE *pulse;
  struct PWL *pwl;
  struct VoltageSource *next;
  
}VoltageSourceT;

//Current Source
typedef struct CurrentSource{
  
  char *name;
  char *pos_term;
  char *neg_term;
  double value;
  char transient_spec[6];
  struct EXP *exp;
  struct SIN *sin;
  struct PULSE *pulse;
  struct PWL *pwl;
  struct CurrentSource *next;
  
}CurrentSourceT;

//Resistor
typedef struct Resistor{
  
  char *name;
  char *pos_term;
  char *neg_term;
  double value;
  struct Resistor *next;
  
}ResistorT;

//Capacitor
typedef struct Capacitor{
  
  char *name;
  char *pos_term;
  char *neg_term;
  double value;
  struct Capacitor *next;
  
}CapacitorT;

//Inductor
typedef struct Inductor{
  
  char *name;
  char *pos_term;
  char *neg_term;
  double value;
  struct Inductor *next;
  
}InductorT;

//Diode
typedef struct Diode{
  
  char *name;
  char *pos_term;
  char *neg_term;
  char *model_name;
  double area;
  struct Diode *next;
  
}DiodeT;

//MOSFET
typedef struct MOS{
  
  char *name;
  char *drain;
  char *gate;
  char *source;
  char *body;
  char *model_name;
  double Lvalue;
  double Wvalue; 
  
  struct MOS *next;
  
}MOST;

//BJT
typedef struct BJT{
  
  char *name;
  char *collector;
  char *base;
  char *emitter;
  char *model_name;
  double area;
  struct BJT *next;
  
}BJTT;


//Root twn listwn opou apothikeuetai kathe stoixeio tou kiklwmatos
VoltageSourceT *rootV;
CurrentSourceT *rootI;
ResistorT *rootR;
CapacitorT *rootC;
InductorT *rootL;
DiodeT *rootD;
MOST *rootM;
BJTT *rootQ;

int ground;	//arxikopoieitai me 0 kai ginetai isi me 1 otan diavasoume komvo '0' (geiwsi).
int SPD;	//einai isi me 1 an exoume diavasei SPD sta OPTIONS. Alliws =0
int ITER;	//einai isi me 1 an exoume diavasei ITER sta OPTIONS. Alliws =0
int SPARSE;	//einai isi me 1 an exoume diavasei SPARSE sta OPTIONS. Alliws =0
int TRAN;		//einai isi me 1 an exoume diavasei .TRAN .Alliws =0
double time_step;	//Vhma metavatikhs analushs
double end_time;	//final time metavatikhs analushs
int METHOD;		//0:TR(DEFAULT) , 1:BE.
double itol_value;	//timi tis itol value an tin diavasoume
int plot;	//einai isi me 1 otan tha kanoume PLOT kapoia dinamika komvwn
int m2; 	//plh8os autepagwgwn kai phgwn tashs
int hashNode_num;	//plithos komvwn -> voithaei stin antistoixisi twn komvwn me akeraies times
int dc_sweep;	// einai isi me 1 an tha kanoume sweep. Alliws =0
int sweep_source;	//An kanoume sweep me pigi tasis, mas deixnei ti thesis tis pigis ston pinaka B. An kanoume sweep me pigi reumatos einai isi me -1
int sweep_posNode;	//An kanoume sweep me pigi reumatos mas deixnei ton positive komvo pou sindeetai auti i pigi(tin akeraia timi tou komvou apo to hashtable)
int sweep_negNode;	//An kanoume sweep me pigi reumatos mas deixnei ton negative komvo pou sindeetai auti i pigi(tin akeraia timi tou komvou apo to hashtable)
double start_value;	//krata to start value gia to sweep(an ginetai)
double end_value;	//krata to end value gia to sweep(an ginetai)
double sweep_step;	//krata to sweep step gia to sweep(an ginetai)
int *plot_nodes;	//einai o pinakas pou apothikeuontai ta ids ton komvwn pou ginontai PLOT
char **plot_names;	//einai o pinakas pou apothikeuontai ta onomata ton komvwn pou ginontai PLOT
int plot_size;		//megethos tou pinaka pou apothikeuontai ta ids ton komvwn pou ginontai PLOT
double sweep_value_I;	//krataei to value tis pigis reumatos pou tha ginei sweep to opoio diavazetai apo to netlist

//Declaration twn sinartisewn tou arxeiou circuit.c
void circuit_init();
void readCircuit(FILE *input_ptr);
void insertV(char *line);
void insertI(char *line);
void insertR(char *line);
void insertC(char *line);
void insertL(char *line);
void insertD(char *line);
void insertM(char *line);
void insertQ(char *line);
void AddToHashtable(char *string);
void printList();
void printHash();
void initPlotFiles(char *str);
void plotFiles(char *str, double *plot_table, double current_value, char *msg);

#endif
