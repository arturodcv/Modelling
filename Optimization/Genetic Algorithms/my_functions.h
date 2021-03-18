typedef struct{   //Store the 11 ODE parameters
  double beta, phi, epsI, epsY;
  double sigma;
  double gamma1, gamma2;
  double kappa;
  double p;
  double alpha;
  double delta;
  double PopSize;
} ODE_Parameters;

#define IC_GENES_NUMBER 3
#define PARAMETERS_GENES_NUMBER 11
typedef struct {  //Store the individuals data type
  unsigned long IC[IC_GENES_NUMBER];  // Initial conditions. In this order: E(0), I1(0), A(0).
  unsigned long Pars[PARAMETERS_GENES_NUMBER]; //Genotype. In this order: beta,phi, epsI, epsY, sigma, gamma1, gamma2, kappa, p, alpha, delta.
double fitness;
} Individual; 

// Data
#define Ndays_Time_Series 101
#define Nvar_Time_Series 5
typedef struct{
	double PopSize;
	unsigned long N_days;
	double Data_Time_Series[Ndays_Time_Series][Nvar_Time_Series];
} DataForFitting;

int Data_from_txt(FILE *ftxt,DataForFitting *TheData);
void Get_data(DataForFitting *TheData); 
void print_data(DataForFitting *TheData);

//Genetic stuff
Individual * create_Population(int population_size);
double uniform(void); 
int random_in(int lower, int upper);
void print_individual(Individual ind);
int compare_fitness (const void * a, const void * b);
int tournament_selection(Individual * population, int population_size);
void Mutation(Individual* ind, double prob);
void BitFlipMutation(Individual* ind, double prob);
void OnePointCrossover(Individual* parent1, Individual* parent2, Individual* child1, Individual* child2, double prob);


// RKF78
#define CoreModelDIM 8
void CoreModel(double t, double *x, unsigned dim, double *der, void *Params);
int GeneratePredictionFromIndividual(double *xt, void *ODE_pars, DataForFitting *Pred);
int CoreModelVersusDataQuadraticError(Individual *ind, void *TheData);
double fitness_function(Individual *ind, DataForFitting *TheData);

//CGDescent
void CGDescent(Individual *ind, void *TheData);
void gradient(Individual *ind, void *TheData,double grad[]);
double alpha_search(Individual *ind, void *TheData,double d[]);
double get_beta(double grad[],double prev_grad[]);
void evaluate(Individual *ind, void *TheData, ODE_Parameters ODE_pars);

// Simulated Annealing
void evalf(Individual *ind,  void *TheData);
void SA(Individual *ind, void *TheData);
void Copy_Individual(Individual *ind,Individual *aux_ind);

// Errors
void Exit_Error(const char *miss, int errorcode);

#define infeasible_value 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.000000

#define euler 2.71828182845904523536
