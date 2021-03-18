#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include <stdbool.h>
#include "my_functions.h"
#include "RKF78.h" 
#include <float.h>
#define	MAXDOUBLE	DBL_MAX

//////////////////////////////////////////////////////////////////////////// Genetic

// Create individual
Individual create_Individual(){   
    Individual ind; 
      
    ind.IC[0]    = floor(uniform()*pow(2,30));
    ind.IC[1]    = floor(uniform()*pow(2,30));  
    ind.IC[2]    = floor(uniform()*pow(2,30)); 
    
    ind.Pars[0]  = floor(uniform()*pow(2,40));     
    ind.Pars[1]  = floor(uniform()*pow(2,20)); 
    ind.Pars[2]  = floor(uniform()*pow(2,10)); 
    ind.Pars[3]  = floor(uniform()*pow(2,10)); 
    ind.Pars[4]  = floor(uniform()*pow(2,20)); 
    ind.Pars[5]  = floor(uniform()*pow(2,20)); 
    ind.Pars[6]  = floor(uniform()*pow(2,20));  
    ind.Pars[7]  = floor(uniform()*pow(2,10)); 
    ind.Pars[8]  = floor(uniform()*pow(2,20)); 
    ind.Pars[9]  = floor(uniform()*pow(2,40));   
    ind.Pars[10] = floor(uniform()*pow(2,40));   
     
    ind.fitness  = MAXDOUBLE; 
    
    return ind;
}

// Create initial population
Individual * create_Population(int population_size){
    Individual * population = malloc(sizeof(Individual) * population_size);
    for(int i = 0; i < population_size; i++){
        *(population + i ) = create_Individual();
    }
    return population;
}

// Function for competition
int compare_fitness (const void * a, const void * b) {
    Individual *ind1 = (Individual *)a;
    Individual *ind2 = (Individual *)b;
    return ( ind1->fitness - ind2->fitness );
}

// Selection
int tournament_selection(Individual * population,int population_size){          
    int first_fighter  = floor(uniform()*population_size); 
    int second_fighter = floor(uniform()*population_size);
    return population[first_fighter].fitness < population[second_fighter].fitness ? first_fighter : second_fighter;
}

//Mutation
void Mutation(Individual* ind, double prob) {
  if(uniform() < prob){
    int i = random_in(0,PARAMETERS_GENES_NUMBER-1);   
    int x;          
    if (i==10 || i == 9 || i == 0)                            {x = random_in(0,pow(2,20))*random_in(0,pow(2,20));} 
    else if (i == 1 || i == 4 || i == 5 || i == 6 || i == 8)  {x = random_in(0,pow(2,20));}
    else if (i == 2 || i == 3 || i == 7)                      {x = random_in(0,pow(2,10));}
    else {x = 0; Exit_Error("When choosing gen to mutate", 7);}
    ind->Pars[i] ^= (1U << x);
  }
}

// Mutation to bitflip an individual's parameters
void BitFlipMutation(Individual* ind, double prob){
  int par = random_in(0,PARAMETERS_GENES_NUMBER-1);
  unsigned char len = 8*sizeof(ind->Pars[par]);
  register unsigned char i;
  for(i=0; i < len; i++) if(uniform() < prob) ind->Pars[par] = (ind->Pars[par])^(1U << i);
}



// Crossover
void OnePointCrossover(Individual* parent1, Individual* parent2, Individual* child1, Individual* child2, double prob){
    Copy_Individual(parent1, child1);
    Copy_Individual(parent2, child2);
    
    if(uniform() < prob){
      unsigned char len = 8*sizeof(unsigned int);
      unsigned char d = uniform()*(len-1) + 1, di = len - d; 
      int par = random_in(0,PARAMETERS_GENES_NUMBER-1);
      child1->Pars[par] = (( parent1->Pars[par] >> d) << d) | (( parent2->Pars[par] << di) >> di);
      child2->Pars[par] = (( parent2->Pars[par] >> d) << d) | (( parent1->Pars[par] << di) >> di);
    }
}


////////////////////////////////////////////////////////////////////////// Debuging

// Print individual's parameters and fitness
void print_individual(Individual ind){
    printf("\nPrint of the individual values\nInitial conditions for ODE values: \n ");
    printf("E(0) | I1(0) | A(0)\n");
    for (int i = 0; i < IC_GENES_NUMBER ; i++ ){
        printf("%lu | ", ind.IC[i]);
    }
    printf("\nGenotype values: \n ");
    printf("beta  |  phi  |  epsI  |  epsY  |  sigma  |  gamma1  |  gamma2  |  kappa  |  p  |  alpha  |  delta\n");
    for (int i = 0; i < PARAMETERS_GENES_NUMBER ; i++ ){
        printf("%lu | ", ind.Pars[i]);
    }
    printf("\nFitness value: %f\n", ind.fitness);
}



//////////////////////////////////////////////////////////////////////////// Random 

double uniform(){  
    return (double)rand() / (double)RAND_MAX ;
}
int random_in(int lower, int upper) {  
    int rnd = (rand() % (upper - lower + 1)) + lower; 
    return rnd; 
} 

//////////////////////////////////////////////////////////////////////////// RKF78 

#define CoreModelDIM 8
void CoreModel(double t, double *x, unsigned dim, double *der, void *Params){
	ODE_Parameters *par = (ODE_Parameters *) Params;                                        // To simplify the usage of Params (void pointer)
	double sigmae = par->sigma*x[1],
	gamma1i1 = par->gamma1*x[2],
	kappaA = par->kappa*x[3],
	alphai2 = par->alpha*x[5];
	der[0] = par->phi*x[2] + x[3] + (1-par->epsI)*(x[4]+x[5]) + (1-par->epsY)*x[6];
	der[0] = - par->beta * (x[0] * der[0])/par->PopSize;                                    // S
	der[1] = -der[0] - sigmae;                                                              // E
	der[2] = sigmae - gamma1i1;                                                             // I_1
	der[3] = (1-par->p)*gamma1i1 - kappaA - par->gamma2*x[3] ;                              // A
	der[4] = kappaA - par->gamma2*x[4];                                                     // A_d
	der[5] = par->p*gamma1i1 - par->gamma2*x[5] - alphai2;                                  // I_2
	der[6] = alphai2 - (par->gamma2+par->delta)*x[6];                                       // Y
	der[7] = par->gamma2*(x[3] + x[4] + x[5] + x[6]);                                       // R
}

#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

int GeneratePredictionFromIndividual(double *xt, void *ODE_pars, DataForFitting *Pred) {
	register unsigned ndays; 
	double t = 0.0, err, h = 1.e-3;
	for(ndays=1; ndays <= Pred->N_days; ndays++) {  
		int status;
		while(t+h < ndays) {
			status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
			if(status) return status;
		}	
		h = ndays - t;
		status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
		if(status) return status;
		Pred->Data_Time_Series[ndays][0] = xt[4]; 
		Pred->Data_Time_Series[ndays][1] = xt[5];
		Pred->Data_Time_Series[ndays][2] = xt[6]; 
		Pred->Data_Time_Series[ndays][3] = xt[7]; 
		Pred->Data_Time_Series[ndays][4] = Pred->PopSize - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]);
	}
	return 0;
}

#define crom2IC(c) (((double) c)/1000)
#define crom2HSPar(c) (((double) c)/1099511627776UL)
#define crom2Par(c) (((double) c)/1048576U)
#define crom2LSPar(c) (((double) c)/1024U)

int CoreModelVersusDataQuadraticError(Individual *ind, void *TheData) {
	
  DataForFitting *TDfF = (DataForFitting *) TheData;
  
	DataForFitting ThePrediction = { TDfF->PopSize, TDfF->N_days, 
  { { TDfF->Data_Time_Series[0][0], TDfF->Data_Time_Series[0][1], TDfF->Data_Time_Series[0][2], TDfF->Data_Time_Series[0][3], TDfF->Data_Time_Series[0][4] } } };
	
  double xt[CoreModelDIM] = { TDfF->PopSize, crom2IC(ind->IC[0]), crom2IC(ind->IC[1]), crom2IC(ind->IC[2]), TDfF->Data_Time_Series[0][0],
		TDfF->Data_Time_Series[0][1], TDfF->Data_Time_Series[0][2], TDfF->Data_Time_Series[0][3] };
	
  xt[0] -= (xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]);
	
  ODE_Parameters ODE_pars = { 
    crom2HSPar(ind->Pars[0]), crom2Par(ind->Pars[1]), crom2LSPar(ind->Pars[2]), crom2LSPar(ind->Pars[3]),
		crom2Par(ind->Pars[4]), crom2Par(ind->Pars[5]), crom2Par(ind->Pars[6]), crom2LSPar(ind->Pars[7]),
		crom2Par(ind->Pars[8]), crom2HSPar(ind->Pars[9]), crom2HSPar(ind->Pars[10]), TDfF->PopSize };
  
  if(GeneratePredictionFromIndividual(xt, &ODE_pars, &ThePrediction)) { ind->fitness = MAXDOUBLE;  return -1; }
	
	double fitn=0; double sum=0; 
	for (int i=0; i<Ndays_Time_Series; i++){
		for (int j=0; j<Nvar_Time_Series; j++){
        sum = sum + pow(ThePrediction.Data_Time_Series[i][j] - TDfF->Data_Time_Series[i][j],2);
    }

		fitn = fitn + (i+1)*sum; 
		sum = 0;
	}
	ind->fitness = fitn;
  return 0;
}

/////////////////////////////////////////////////////////////////// CONJUGATE GRADIENT DESCENT

void evaluate(Individual *ind, void *TheData, ODE_Parameters ODE_pars) {
	
  DataForFitting *TDfF = (DataForFitting *) TheData;
  
	DataForFitting ThePrediction = { TDfF->PopSize, TDfF->N_days, 
  { { TDfF->Data_Time_Series[0][0], TDfF->Data_Time_Series[0][1], TDfF->Data_Time_Series[0][2], TDfF->Data_Time_Series[0][3], TDfF->Data_Time_Series[0][4] } } };
	
  double xt[CoreModelDIM] = { TDfF->PopSize, crom2IC(ind->IC[0]), crom2IC(ind->IC[1]), crom2IC(ind->IC[2]), TDfF->Data_Time_Series[0][0],
		TDfF->Data_Time_Series[0][1], TDfF->Data_Time_Series[0][2], TDfF->Data_Time_Series[0][3] };
	
  xt[0] -= (xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]);
  
  if(GeneratePredictionFromIndividual(xt, &ODE_pars, &ThePrediction)) { ind->fitness = MAXDOUBLE;}
	
	double fitn=0; double sum=0; 
	for (int i=0; i<Ndays_Time_Series; i++){
		for (int j=0; j<Nvar_Time_Series; j++){
        sum = sum + pow(ThePrediction.Data_Time_Series[i][j] - TDfF->Data_Time_Series[i][j],2);
    }
		fitn = fitn + (i+1)*sum; 
		sum = 0;
	}
	ind->fitness = fitn;
}

void gradient(Individual *ind, void *TheData, double grad[]){
  Individual *aux_ind; 
  if ((aux_ind = (Individual*) malloc(sizeof(Individual))) == NULL) Exit_Error("When allocating memory for individual", 9);
  DataForFitting *TDfF = (DataForFitting *) TheData;
  double h = 0.000001;
  
  for (int par = 0;par<PARAMETERS_GENES_NUMBER;par++){
    Copy_Individual(ind,aux_ind);
    ODE_Parameters ODE_pars = { crom2HSPar(aux_ind->Pars[0]), crom2Par(aux_ind->Pars[1]),   crom2LSPar(aux_ind->Pars[2]),  crom2LSPar(aux_ind->Pars[3]),
 	                              crom2Par(aux_ind->Pars[4]),   crom2Par(aux_ind->Pars[5]),   crom2Par(aux_ind->Pars[6]),    crom2LSPar(aux_ind->Pars[7]),
		                            crom2Par(aux_ind->Pars[8]),   crom2HSPar(aux_ind->Pars[9]), crom2HSPar(aux_ind->Pars[10]), TDfF->PopSize };
    
    if (par == 0) ODE_pars.beta += h;        else if (par == 1) ODE_pars.phi += h;    else if (par == 2) ODE_pars.epsI += h;
    else if (par == 3) ODE_pars.epsY += h;   else if (par == 4) ODE_pars.sigma += h;  else if (par == 5) ODE_pars.gamma1 += h;
    else if (par == 6) ODE_pars.gamma2 += h; else if (par == 7) ODE_pars.kappa += h;  else if (par == 8) ODE_pars.p += h;
    else if (par == 9) ODE_pars.alpha += h;  else if (par == 10) ODE_pars.delta += h;
    
    evaluate(aux_ind, TheData, ODE_pars);
    grad[par] = (ind->fitness - aux_ind->fitness)/(h);
  }
  free(aux_ind);
}

double alpha_search(Individual *ind, void *TheData, double d[]){
  Individual *aux_ind; 
  if ((aux_ind = (Individual*) malloc(sizeof(Individual))) == NULL) Exit_Error("When allocating memory for individual", 10);
  DataForFitting *TDfF = (DataForFitting *) TheData;
  
  double f_value = MAXDOUBLE;
  double alpha = 0, min_alpha = MAXDOUBLE;
  double step = 0.000001;
  for (int k = 0;k<200;k++){
    Copy_Individual(ind,aux_ind);
    ODE_Parameters ODE_pars = { crom2HSPar(aux_ind->Pars[0]), crom2Par(aux_ind->Pars[1]),   crom2LSPar(aux_ind->Pars[2]),  crom2LSPar(aux_ind->Pars[3]),
 	                              crom2Par(aux_ind->Pars[4]),   crom2Par(aux_ind->Pars[5]),   crom2Par(aux_ind->Pars[6]),    crom2LSPar(aux_ind->Pars[7]),
		                            crom2Par(aux_ind->Pars[8]),   crom2HSPar(aux_ind->Pars[9]), crom2HSPar(aux_ind->Pars[10]), TDfF->PopSize };
    alpha = k*step; 
    for (int par = 0;par<PARAMETERS_GENES_NUMBER;par++){ 
      if (par == 0)       ODE_pars.beta   += alpha*crom2HSPar(d[par]);      else if (par == 1) ODE_pars.phi    += alpha*crom2Par(d[par]);   
      else if (par == 2)  ODE_pars.epsI   += alpha*crom2LSPar(d[par]);      else if (par == 3) ODE_pars.epsY   += alpha*crom2LSPar(d[par]);   
      else if (par == 4)  ODE_pars.sigma  += alpha*crom2Par(d[par]);        else if (par == 5) ODE_pars.gamma1 += alpha*crom2Par(d[par]);
      else if (par == 6)  ODE_pars.gamma2 += alpha*crom2Par(d[par]);        else if (par == 7) ODE_pars.kappa  += alpha*crom2LSPar(d[par]); 
      else if (par == 8)  ODE_pars.p      += alpha*crom2Par(d[par]);        else if (par == 9) ODE_pars.alpha  += alpha*crom2HSPar(d[par]);  
      else if (par == 10) ODE_pars.delta  += alpha*crom2HSPar(d[par]);
    
    evaluate(aux_ind, TheData, ODE_pars);
    if (aux_ind->fitness < f_value) {f_value = aux_ind->fitness; min_alpha = alpha;}
    }
  }
  free(aux_ind);
  return min_alpha;  
}

double get_beta(double grad[], double prev_grad[]){
  double beta = 0;
  for (int gi = 0; gi < PARAMETERS_GENES_NUMBER; gi++){
    beta += grad[gi]*(grad[gi]-prev_grad[gi])/(grad[gi]*prev_grad[gi]); 
  }
  if (beta < 0) beta = 0;
   
  return beta;
}

#define Par2HSDiscPar(c) (c*1099511627776UL)
#define Par2DiscPar(c) (c*1048576U)
#define Par2LSDiscPar(c) (c*1024U)

void CGDescent(Individual *ind, void *TheData){
  Individual *aux_ind; 
  if ((aux_ind = (Individual*) malloc(sizeof(Individual))) == NULL) Exit_Error("When allocating memory in CGDescent", 11);
  Copy_Individual(ind,aux_ind);
  
  double beta, alpha = 1;
  double grad[PARAMETERS_GENES_NUMBER], prev_grad[PARAMETERS_GENES_NUMBER], d[PARAMETERS_GENES_NUMBER], s[PARAMETERS_GENES_NUMBER];
   
  gradient(ind,TheData,grad);
  
  
  for (int par = 0;par<PARAMETERS_GENES_NUMBER;par++){ d[par] = -grad[par];}
  alpha = alpha_search(ind, TheData, d); 
  for (int par = 0;par<PARAMETERS_GENES_NUMBER;par++){ ind->Pars[par] += alpha*d[par]; s[par] = d[par];}
  
  int iter_max = 10; 
  for (int iter = 0; iter < iter_max; iter++){
    for (int par = 0;par<PARAMETERS_GENES_NUMBER;par++){ prev_grad[par] = grad[par];}
    gradient(ind, TheData,grad);
    beta = get_beta(grad,prev_grad); 
    for (int par = 0;par<PARAMETERS_GENES_NUMBER;par++){ s[par] = d[par] + beta*s[par];}
    alpha = alpha_search(ind, TheData, d);
    ind->Pars[0] += alpha*Par2HSDiscPar(s[0]);  ind->Pars[1] += alpha*Par2DiscPar(s[1]);       ind->Pars[2] += alpha*Par2LSDiscPar(s[2]);
    ind->Pars[3] += alpha*Par2HSDiscPar(s[3]);     ind->Pars[4] += alpha*Par2DiscPar(s[4]);       ind->Pars[5] += alpha*Par2DiscPar(s[5]);
    ind->Pars[6] += alpha*Par2DiscPar(s[6]);    ind->Pars[7] += alpha*Par2LSDiscPar(s[7]);     ind->Pars[8] += alpha*Par2DiscPar(s[8]);
    ind->Pars[9] += alpha*Par2HSDiscPar(s[9]);  ind->Pars[10] += alpha*Par2HSDiscPar(s[10]);
  }
  free(aux_ind);
}


////////////////////////////////////////////////////////////////// SIMULATED ANNEALING

void SA(Individual *ind, void *TheData){
  Individual *aux_ind; 
  if ((aux_ind = (Individual*) malloc(sizeof(Individual))) == NULL) Exit_Error("When allocating memory in SA", 12);
  DataForFitting *TDfF = (DataForFitting *) TheData;
  
  double accept_perturbation; 
  double goal = ind->fitness*7/10; 
  int iters = 0, iters_max = 10;

  while (ind->fitness > goal && iters < iters_max){
      for (int par = 0; par<PARAMETERS_GENES_NUMBER;par++){ 
          Copy_Individual(ind,aux_ind);
          
          ODE_Parameters ODE_aux_pars = { crom2HSPar(aux_ind->Pars[0]), crom2Par(aux_ind->Pars[1]),   crom2LSPar(aux_ind->Pars[2]),  crom2LSPar(aux_ind->Pars[3]),
 	                                        crom2Par(aux_ind->Pars[4]),   crom2Par(aux_ind->Pars[5]),   crom2Par(aux_ind->Pars[6]),    crom2LSPar(aux_ind->Pars[7]),
		                                      crom2Par(aux_ind->Pars[8]),   crom2HSPar(aux_ind->Pars[9]), crom2HSPar(aux_ind->Pars[10]), TDfF->PopSize };
          
          if      (par == 0)  ODE_aux_pars.beta   += random_in(-1,1)*uniform();     else if (par == 1) ODE_aux_pars.phi    += random_in(-1,1)*uniform();   
          else if (par == 2)  ODE_aux_pars.epsI   += random_in(-1,1)*uniform();     else if (par == 3) ODE_aux_pars.epsY   += random_in(-1,1)*uniform();   
          else if (par == 4)  ODE_aux_pars.sigma  += random_in(-1,1)*uniform();     else if (par == 5) ODE_aux_pars.gamma1 += random_in(-1,1)*uniform();
          else if (par == 6)  ODE_aux_pars.gamma2 += random_in(-1,1)*uniform();     else if (par == 7) ODE_aux_pars.kappa  += random_in(-1,1)*uniform();
          else if (par == 8)  ODE_aux_pars.p      += random_in(-1,1)*uniform();     else if (par == 9) ODE_aux_pars.alpha  += random_in(-1,1)*uniform();  
          else if (par == 10) ODE_aux_pars.delta  += random_in(-1,1)*uniform();
           
          evaluate(aux_ind, TheData, ODE_aux_pars);
      
          accept_perturbation = pow(euler, -(aux_ind->fitness - ind->fitness)/ind->fitness);     
          
          if ((aux_ind->fitness < infeasible_value) && ((aux_ind->fitness < ind->fitness) || (uniform() < accept_perturbation) )){     
            Copy_Individual(aux_ind,ind);
          }   
      }      
      iters++;  
  } 
  free(aux_ind);
}

void Copy_Individual(Individual *ind,Individual *aux_ind){
    for (int ic = 0; ic < IC_GENES_NUMBER; ic++){aux_ind->IC[ic] = ind->IC[ic];}
    for (int gen = 0; gen < PARAMETERS_GENES_NUMBER; gen++){aux_ind->Pars[gen] = ind->Pars[gen];}  
    aux_ind->fitness = ind->fitness;
}



//////////////////////////////////////////////////////////////////////////// Errors
void Exit_Error(const char *miss, int errorcode) {
    fprintf (stderr, "\n ERROR: %s.\n Stopping...\n\n", miss); exit(errorcode);
}


