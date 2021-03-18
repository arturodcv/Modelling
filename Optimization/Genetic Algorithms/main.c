#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include <stdbool.h>
#include <float.h>
#include "my_functions.h"

////////////////////////////////////////////////////////////////////////////////////
////////////
////////////  to compile use
////////////  gcc -g -Wall -O3 -o main main.c my_functions.c RKF78.c read_data.c -lm
////////////
////////////////////////////////////////////////////////////////////////////////////

 
int main(){
    DataForFitting *TheData;  
    if ((TheData = malloc (sizeof(DataForFitting))) == NULL){
        Exit_Error("when allocating memory for the data for fitting matrix", 1);}
    Get_data(TheData);
    
    //time_t t; int seed = (unsigned) time(&t); srand(seed); printf("\nSeed = %d\n", seed);
    time_t t; srand((unsigned) time(&t)); 
    
    int population_size, max_num_generation;

    printf("\nSize of the population: "); 
    if (scanf("%d", &population_size) != 1) Exit_Error("when reading the size of the population", 2);
    while ( population_size <= 0 ) {
      printf("\nInvalid choice of population size. The number must be a positive integer. Please enter the desired size of the population again: ");
        if (scanf("%d", &population_size) != 1) Exit_Error("when reading the size of the population", 3);
    }
    
    printf("\nNumber of generations: "); 
    if (scanf("%d", &max_num_generation) != 1) Exit_Error("when reading the number of the generations", 4);
    while ( max_num_generation <= 0 ) {
      printf("\nInvalid choice of number of generations. The number must be a positive integer. Please enter number of generations again: ");
        if (scanf("%d", &max_num_generation) != 1) Exit_Error("when reading the number of the generations", 5);
    }
    
    if (population_size % 2 == 1) population_size += 1;
    time_t GA_time;
    
    
    FILE * output;
    output=fopen("output_genetic.txt","w");
    fprintf(output, "Genetic Algorithm started for a population with %d individuals and %d generations\n",population_size,max_num_generation);
     
    GA_time = clock(); 
    
    printf("Creating a population with %i individuals\n", population_size);
    Individual * population = create_Population(population_size);
    
    int generation = 0;
    double fitness_goal = 100, fittest = DBL_MAX;
    
    while (generation <= max_num_generation && population->fitness > fitness_goal){
      Individual * next_generation = NULL;
      if ((next_generation = (Individual*) malloc(population_size*sizeof(Individual))) == NULL)  
          Exit_Error("When allocating memory for the next generation", 6);
      
      printf("\n\n\n------------- Generation Number %i ------------- \n",generation);  
      fprintf(output,"\n\n------------- Generation Number %i ------------- \n",generation);
      
      double infeasibles = 0;      

      int best = random_in(0,population_size-1),second_best = random_in(0,population_size-1);  //Elitism: the two best individuals go to the next generation
      double fittest_val = infeasible_value;

      for (int i=0;i < population_size; i++){
          if (((population+i)->fitness) < fittest_val) {fittest_val = (population+i)->fitness; second_best = best; best = i;}
      } 
      
      OnePointCrossover(population+best,(population+second_best),(next_generation),(next_generation+1),0);
      CoreModelVersusDataQuadraticError((next_generation),TheData); 
      CoreModelVersusDataQuadraticError((next_generation+1),TheData); 
      
      for (int i = 1;i < population_size/2; i++){ 
          int is_feasible = 1;
          double prob_mutation;
          while (is_feasible){     // To control the proportion of infeasible individuals 
            
            Individual *parent1 = population + tournament_selection(population, population_size);   // Selection for crossover
            Individual *parent2 = population + tournament_selection(population, population_size);
            OnePointCrossover(parent1, parent2, (next_generation + 2*i), (next_generation + 2*i + 1),0.85); 
            
            for (int child = 0; child <= 1; child++) {  
              
              prob_mutation = (next_generation + 2*i + child)->fitness/(pow(10,12) + (next_generation + 2*i + child)->fitness); // Dynamic decreasing mutation
              Mutation((next_generation + 2*i + child), prob_mutation); 
               
              CoreModelVersusDataQuadraticError((next_generation + 2*i + child),TheData); 
            
              if (((next_generation + 2*i + child)->fitness < infeasible_value)) {
                //CGDescent((next_generation + 2*i + child),TheData);    // Conjugate Gradient Descent
                //SA((next_generation + 2*i + child),TheData);             // Simulated Annealing 
              }
            }  
            
            if (((next_generation + 2*i)->fitness < infeasible_value) || ((next_generation + 2*i + 1)->fitness < infeasible_value)) is_feasible = 0; 
            
          }
      }
      
      Individual * aux = population;
      population = next_generation;
      free(aux); 
      
      fittest = DBL_MAX;
      double mean = 0;
      for (int j = 0;j<population_size;j++){
          if ((population+j)->fitness < infeasible_value){mean  = mean + (population+j)->fitness;}
          else infeasibles++;
          if ((population+j)->fitness < fittest) {fittest = (population+j)->fitness;}
          }
          
      printf("\nFittest individual fitness value                      = %f\n",fittest);
      fprintf(output,"\nFittest individual fitness value                      = %f\n",fittest);
      //fprintf(output,"%f\n", fittest);
      printf("Mean fitness of the feasibles in population           = %f\n",mean/(population_size-infeasibles));
      fprintf(output,"Mean fitness of the feasibles in population           = %f\n",mean/(population_size-infeasibles));
      //fprintf(output,"%f\n", mean/(population_size-infeasibles));
      printf("Proportion of infeasibles individuals in generation %i = %f\n",generation,infeasibles/population_size);
      fprintf(output,"Proportion of infeasibles individuals                 = %f\n",infeasibles/population_size);
      //fprintf(output,"%f\n",infeasibles/population_size );
       
      generation++;
    }
    GA_time = clock()- GA_time;
    if (generation == max_num_generation+1){ 
        printf("\n\nGenetic algorithm has ended because it has reached the maximum number of generations. Total time spents is %f seconds \nThe fittest individual achieved a score of %f \n", ((double)GA_time)/CLOCKS_PER_SEC,  fittest);
        fprintf(output,"\n\nGenetic algorithm has ended because it has reached the maximum number of generations. Total time spent is %f seconds \nThe fittest individual achieved a score of %f \n",    ((double)GA_time)/CLOCKS_PER_SEC, fittest);
    }
    
    else {printf("\nGenetic algorithm has converged to a solution in %f seconds! The individual thought to be the one that started the COVID-19 pandemia is\n", ((double)GA_time)/CLOCKS_PER_SEC);
          print_individual(population[0]);
          fprintf(output,"\nGenetic algorithm has converged to a solution in %f seconds! Check the compiling screen to get the parameters of the individual thought to be the one that started the COVID-19 pandemia.\n", ((double)GA_time)/CLOCKS_PER_SEC);
    }
    
    free(population);
    return 0;
} 












