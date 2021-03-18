#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include <stdbool.h>
#include "functions.h" 


#define PI 4*atan(1)
#define R 6371 // mean earth ratio in kms

unsigned long num_nodes = 23895681UL;

Node * nodes;
Node_status* node_status;

//Read file
void reading_from_file(){
	int i;
	static unsigned long * allsuccesors;
	unsigned long ntotnum_succesors=0UL;
	FILE * file_in;
  
  // Open 
  if ((file_in = fopen("binary.bin","rb")) == NULL)
		Exit_Error("The data file does not exist or cannot be opened\n",14);
	
 // Header
  if ((fread(&num_nodes,sizeof(unsigned long),1,file_in)+fread(&ntotnum_succesors,sizeof(unsigned long),1,file_in))!=2) 
    Exit_Error("when reading the header oh the binary file\n ", 15);
	
  // Memory
  if ((nodes=malloc(sizeof(Node)*num_nodes))==NULL)  
    Exit_Error("when allocating memory for the nodes vector\n ", 16);
	
  if((allsuccesors=malloc(sizeof(unsigned long)*ntotnum_succesors))==NULL) 
    Exit_Error("when allocating memory for the edges vector\n ", 17);
 
  // Data  
  if (fread(nodes,sizeof(Node),num_nodes,file_in)!=num_nodes) 
    Exit_Error("when reading nodes from the binary data file\n ", 18);
	
  if (fread(allsuccesors,sizeof(unsigned long),ntotnum_succesors,file_in)!=ntotnum_succesors) 
    Exit_Error("when reading succesors from the binary data file\n ", 19);
	
  fclose(file_in);
	for (i=0;i<num_nodes; i++) 
    if (nodes[i].num_succesors){
		  nodes[i].succesors=allsuccesors;
		  allsuccesors+=nodes[i].num_succesors;}  
  }

//Initializate list
void Init_list ( List * list){  
 list->head = list->tail = NULL; 
}

//Initializate values
void Init_values (Node_status * node_status, Node * nodes){ 
  int i;
  for ( i = 0; i < num_nodes; i ++){
			node_status[i].latitude = nodes[i].latitude;
			node_status[i].longitude =	nodes[i].longitude;
			node_status[i].num_succesors = nodes[i].num_succesors;
			node_status[i].succesors = nodes[i].succesors;
			node_status[i].g = INFINITY;
      node_status[i].queue_status = NONE;
			node_status[i].h = INFINITY;
			node_status[i].id = nodes[i].id;
			node_status[i].parent = NULL;
		}
}


// Equirectangular distance function
double h_equir (Node_status * node_1, Node_status * node_2) { 
	double phi1 = node_1 -> latitude/180*PI;
	double phi2 = node_2 -> latitude/180*PI;
	double lambda1 = node_1 -> longitude/180*PI;
	double lambda2 = node_2 -> longitude/180*PI;
 
	double dif_long = fabs(lambda1-lambda2);
	double dif_lat = fabs(phi1-phi2);
 
	return R*sqrt(pow(dif_long*cos((phi1+phi2)/2),2) + pow(dif_lat,2));
}


// Haversine distance function
double h_haversine(Node_status * node_1, Node_status * node_2 ) { 	 
	double phi1= node_1->latitude/180*PI;
	double phi2= node_2->latitude/180*PI;
	double lambda1=node_1->longitude/180*PI;
	double lambda2=node_2->longitude/180*PI;
  
  double dif_long = fabs(lambda1-lambda2);
	double dif_lat = fabs(phi1-phi2);
	
	double a = sin(dif_lat/2)*sin(dif_lat/2)+cos(phi1)*cos(phi2)*sin(dif_long/2)*sin(dif_long/2);
	return R*2*atan2(sqrt(a),sqrt(1-a));
}

// Spherical distance function
double h_spherical (Node_status * node_1, Node_status * node_2 ){ 
	double phi1= node_1->latitude/180*PI;
	double phi2= node_2->latitude/180*PI;
	double lambda1=node_1->longitude/180*PI;
	double lambda2=node_2->longitude/180*PI;
 
	double dif_long= fabs(lambda1-lambda2);
	return R*acos( sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(dif_long));
}


double evaluation_function (Node_status * succesor_node, int evaluation, double param, unsigned int start_node, unsigned int goal_node) {
    if (evaluation == 1){                                                                    // default  
      return succesor_node -> g + succesor_node -> h;}                    
    
    else if (evaluation == 2){
      return (1-param)*(succesor_node -> g) + param*(succesor_node -> h);}                   // weighted 
    
    else if (evaluation == 3) {                                                              // dynamic weighting
        double extra = 1 - h_equir( succesor_node, &(node_status[start_node]) ) / ( h_equir(succesor_node,&(node_status[goal_node])) + h_equir(succesor_node , &(node_status[start_node])) );
        return succesor_node -> g + succesor_node -> h + param*extra*(succesor_node -> h);}
    Exit_Error("when evaluating function", 22);
    return 0.0f;
}

// Binary search function
long binary_search(unsigned long key)
{  
	unsigned long start=0UL, afterend=num_nodes, middle;
	unsigned long try;
	while(afterend > start)
	{ 
		middle = start + ((afterend-start-1)>>1); try = node_status[middle].id; 
		if (key == try){
        return middle;}
		else if ( key > try ) start=middle+1;
		else afterend=middle;
	}
	return -1;
}

// Remove node from the open list
void pop_node(List * list, Node_status * target){ 
  if ( target == list->head){ 
		list->head = target->next; 
		(target->next)->prev = NULL;}
  else if ( target == list->tail){ 
		list->tail = target->prev; 
		(target->prev)->next = NULL;}
  else {
		(target->next) -> prev = target -> prev; 
		(target->prev) -> next = target -> next ; }
}

// Insert a node at the end of the open list
void insert_node_at_end(List * list, Node_status* target , Node_status* to_insert){
		list-> tail = to_insert; 
		target-> next = to_insert;
		to_insert-> prev = target ;
		to_insert->next = NULL;
}

// Insert a node in the open list
void insert_node_before(List * list, Node_status * target , Node_status * to_insert){
	if ( list->head == target){ 
		target->prev = to_insert; 
		list->head = to_insert; 
		to_insert->prev = NULL;
		to_insert->next = target; }
	else {
		to_insert-> prev = (target-> prev);
		to_insert-> next = target ;
		(target->prev)-> next = to_insert ;
		target-> prev = to_insert;} 
} 

// Add a node to the open_list
void add_to_open_list (List * list, Node_status * succesor){
  Node_status * aux; 
  aux = list->head; 
  
  if ( list->head == NULL ){ 
    list->head = list->tail = succesor; 
    succesor->prev = succesor->next = NULL;} 
  while ( (aux->next != NULL) && ((aux->next)->f <= succesor->f ) ) {  
    aux = aux->next;}  
    
  if (aux->next == NULL) { 
    aux->next = succesor;
    insert_node_at_end(list, list->tail, succesor);}   
  else if ((aux->next)->f > succesor->f) {                                            
    insert_node_before (list, aux, succesor);}
}   

// A* Algorithm 
int Astar_algorithm (unsigned int start_node, unsigned int goal_node, int evaluation, double param)
{
	Node_status *current_node , *succesor_node;

	int i; 
	double succesor_current_cost; 
	List open_list;
	Init_list(&open_list);

	open_list.head = open_list.tail = &(node_status[start_node]); 
	node_status[start_node].g = 0;
	node_status[start_node].h = h_equir(&(node_status[start_node]),&(node_status[goal_node]));
  
  while (open_list.head != NULL){ 

		current_node = open_list.head ; 
		if ( current_node->id == node_status[goal_node].id) {
        return 1;} 
        for ( i = 0; i < (current_node->num_succesors) ; i++){
		      succesor_node = & node_status[*(current_node->succesors+i)]; 
			    succesor_current_cost = current_node -> g + h_equir(succesor_node,current_node); 
          if ( succesor_node -> queue_status == OPEN){
			      if ( succesor_node -> g <= succesor_current_cost){ continue;}			
			        pop_node(&open_list, succesor_node);
		      }
          else if ( succesor_node -> queue_status == CLOSED){		
			      if ( succesor_node -> g <= succesor_current_cost) {continue;}			
              succesor_node -> queue_status = OPEN; 
		      }
          else{ 
            succesor_node -> queue_status = OPEN;}              				 
       succesor_node -> g = succesor_current_cost;
		   succesor_node -> h = h_equir (succesor_node, &(node_status[goal_node])); 
       succesor_node -> f = evaluation_function(succesor_node, evaluation, param, start_node, goal_node); 
       add_to_open_list ( &open_list, succesor_node);				   
		   succesor_node -> parent = current_node;
		   }      
   current_node -> queue_status = CLOSED;
   pop_node (&open_list,current_node); 
	}
	return -1;
}

void print_path(int goal_node, double reading_time, double Astar_time, double total_time){
      FILE * output;
   		output=fopen("output_Astar.txt","w");
 
      fprintf(output, "Reading binary file done. Elapsed time: %fs CPU seconds \n",((double)reading_time)/CLOCKS_PER_SEC);
	  	fprintf(output,"\nA* algorithm has finished. Elapsed time: %fs  CPU seconds \n",((double)Astar_time)/CLOCKS_PER_SEC);
      
      fprintf(output,"\nThis is the path from start node to goal node A* found:\n"); 
	  		
	  	fprintf(output, " Node Id | Distance | Latitude | Longitude \n");
      //printf(" Iteration | Node Id | Distance | Latitude | Longitude \n");
	    
      Node_status* path_node;
      Node_status* aux_node; 
      if ((aux_node = malloc ( sizeof(Node_status)*num_nodes)) == NULL){
        Exit_Error("when allocating memory for the nodes vector", 20);}

      int i;  // We have saved the path from goal to start. We use an aux structure to print the path correctly.  
	  	for ( i=0,path_node = &(node_status[goal_node]); path_node; i++, path_node = path_node -> parent){
        aux_node[i].id = path_node -> id; 
			  aux_node[i].g  = path_node->g;  
        aux_node[i].latitude = path_node->latitude;
        aux_node[i].longitude = path_node->longitude;                     
		  }    
      
      int length_path = i;      
      for ( i=1; i <= length_path; i++){
        // Uncomment this to see the output in the terminal
		    //printf(" %d | %10.lu | %3.4f | %f | %f   \n", i-1,aux_node[length_path-i].id, aux_node[length_path-i].g , aux_node[length_path-i].latitude, aux_node                                                                 [length_path-i].longitude);
        fprintf(output,"%10.lu | %3.4f | %f | %f\n", aux_node[length_path-i].id, aux_node[length_path-i].g,aux_node[length_path-i].latitude, aux_node                                                                 [length_path-i].longitude);                                                       
		  }  
} 
 
// Errors
void Exit_Error(const char *miss, int errorcode) {
    fprintf (stderr, "\n ERROR: %s.\n Stopping...\n\n", miss); exit(errorcode);
}
