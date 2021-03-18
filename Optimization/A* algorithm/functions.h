// Structure definitions
typedef struct {
	unsigned long id;
  unsigned short name_lenght;
	char *name;
	double latitude, longitude;
	unsigned short num_succesors;
	unsigned long *succesors;
}Node;

typedef char Queue;
  enum whichQueue { NONE, OPEN, CLOSED };

typedef struct Node{
	unsigned long id;
	double latitude, longitude;
	unsigned short num_succesors;
	unsigned long *succesors;
	struct Node *prev, *next,*parent; 
  Queue queue_status;
	double f, h, g;   
}Node_status;
 
typedef struct{
  Node_status *head, *tail; 
} List; 


// Init functions
void Init_list ( List * list);
void Init_values (Node_status * node_status, Node * nodes);

// Heuristic functions 
double h_equir(Node_status *node_1, Node_status *node_2 );
double h_haversine(Node_status * node_1, Node_status * node_2 );
double h_spherical (Node_status * node_1, Node_status * node_2 ); 

double evaluation_function (Node_status * succesor_node, int mode, double param, unsigned int start_node, unsigned int goal_node);

 // Binary Search
long binary_search(unsigned long key);

//List functions
void pop_node(List * list, Node_status * target);
void insert_node_at_end(List * list, Node_status* target , Node_status* to_insert);
void insert_node_before(List * list, Node_status * target , Node_status * to_insert);
void add_to_open_list (List * list, Node_status * succesor_node);

//Other functions
void reading_from_file();
int Astar_algorithm (unsigned int start_node, unsigned int goal_node, int evaluation, double param);
void print_path(int goal_node, double reading_time, double Astar_time, double total_time);
void Exit_Error(const char *miss, int errorcode);