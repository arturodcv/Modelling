#define PI 4*atan(1)
#define R 6371
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include "functions.h"

Node * nodes;
Node_status* node_status;




int main (){	
		time_t Astar_time; // A* algorithm time
		time_t reading_time; // reading function time
		time_t total_time; // total time of code
    time_t output_time; // time in writing the output in a txt file
		total_time = clock();

		unsigned long num_nodes = 23895681UL; 
    if ((node_status = malloc ( sizeof(Node_status)*num_nodes)) == NULL){
        Exit_Error("when allocating memory for the nodes vector", 1);}
	    
    
    printf("Reading the binary file... ...");    
		//Read file
    reading_time = clock();
		reading_from_file();
		reading_time = clock()-reading_time;
    printf(" ... Binary file reading finalized \n \n");
   
    //Initializate values
    Init_values(node_status, nodes);
		free(nodes);
   
    int start_node, goal_node, start_node_id, goal_node_id;
    int decission;
    
    printf("Do you want the default starting and goal nodes? Or do you want to customize them? \n Default: From Santa Maria del Mar, Barcelona to Giralda, Sevilla            \n");
    printf("\n Default: 0 \n Custom: 1 \n \n Please select:");
    if (scanf("%d", &decission) != 1) Exit_Error("when reading the id of starting node", 2);
    if (decission == 0){ 
      start_node_id = 240949599;
      goal_node_id  = 195977239;
      start_node = binary_search(start_node_id);
      goal_node = binary_search(goal_node_id);
      }
    
    else if (decission == 1){
      printf("\nChoose IDs for start and goal nodes:\n");  
      printf("ID of start node: "); 
      if (scanf("%d", &start_node_id) != 1) Exit_Error("when reading the id of starting node", 3);
      while ( (start_node = binary_search(start_node_id)) == -1 ) {
          printf("Invalid choice of starting node. Please enter starting node ID again: ");
          if (scanf("%d", &start_node_id) != 1) Exit_Error("when reading the ID of starting node", 4);
      }
    
      printf("ID of goal node: "); 
      if (scanf("%d", &goal_node_id) != 1) Exit_Error("when reading the id of starting node", 5);
      while ( (goal_node = binary_search(goal_node_id)) == -1 ) {
          printf("Invalid choice of starting node. Please enter starting node ID again: ");
          if (scanf("%d", &goal_node_id) != 1) Exit_Error("when reading the ID of starting node", 6);
      }
    }
    else {Exit_Error("when reading the custom or default start and goal nodes", 7);
          return 0;}
    
		printf("Index of start_node node:  %d \n", start_node);
		printf("Index of goal_node node: %d \n", goal_node );
   
    printf("Choose the evaluation function:\n");
    printf("\t1 - Default:              f = g + h\n");
    printf("\t2 - Weighted:             f = (1-w)*g + w*h\n");
    printf("\t3 - Dynamic weighting:    f = g + h + e*(1-d/N)*h\n");
    int evaluation = 0;
    if (scanf("%d", &evaluation) != 1) Exit_Error("when reading the evaluation function", 8);
    while ( evaluation < 1 || evaluation > 3) {
        printf("Invalid choice of evaluation function. Please enter your choice of evaluation again: ");
        if (scanf("%d", &evaluation) != 1) Exit_Error("when reading the evaluation function\n", 9);
    }
    double param = 0;
    if (evaluation == 1) printf("You have chosen default evaluation function.\n");
    else if (evaluation == 2) {
        printf("You have chosen weighted evaluation. Please input parameter w between 0 and 1: ");
        if (scanf("%lf", &param) != 1) Exit_Error("when reading the parameter for weighted evaluation function \n", 10);
        while ( param < 0 || param > 1) {
            printf("Invalid choice of w. Please enter w again: ");
            if (scanf("%lf", &param) != 1) Exit_Error("when reading the parameter for weighted evaluation", 11);
        }
    }
    else if (evaluation == 3) {
        printf("You have chosen dynamic weighting evaluation. Please select parameter epsilon (e): \n");
        if (scanf("%lf", &param) != 1) Exit_Error("when reading the parameter for weighted evaluation", 12);
    }  
   
		printf("Starting the A* algorithm\n");
		Astar_time=clock();	 
   
    // A* algorithm
		if ( Astar_algorithm(start_node,goal_node, evaluation, param) == 1){
      Astar_time=clock()-Astar_time;
	  	total_time = clock() - total_time;
	  	printf("A* has found a path \n");

      
      
	  	
      // Print the path from the start node to the goal node
      output_time = clock();           
      print_path(goal_node, reading_time, Astar_time, total_time);
      output_time = clock() - output_time;
		  
      printf("Time elapsed by A* Algorithm %fs\n", ((double)Astar_time)/CLOCKS_PER_SEC );
      printf("Time elapsed by writing the output %fs \n",((double)output_time)/CLOCKS_PER_SEC);

		}
   
    else {printf("Astar algorithm did not found a path from start node to goal node \n"); 
          Exit_Error("when doing the Astar_algorithm", 13); 
          total_time = clock() - total_time;}
		return 1;
	}
