#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

// structure for the nodes vector
typedef struct {
	unsigned long id;  //node identification
  unsigned short name_lenght;
	char *name;
	double latitude, longitude;  // node position
	unsigned short num_succesors; //number of node successors
	unsigned long *succesors;
}Node;

// structure to store the successors that will be connected
typedef struct{
	unsigned long id;
	int position;
}succesor;


// global declarations
unsigned long nnodes = 23895681UL;

Node * nodes;

succesor * stored = NULL;

short num = 0; //number of nodes in the successor vector

unsigned int k = 0; // index for the nodes vector


// Returns 0 if is type node, 1 if is type way and 2 if relation. -1 if none of the above
int line_type(char *piece)
{
	if(strcmp(piece,"node")==0) {return 0;}
	else if (strcmp(piece,"way")==0) {return 1;}
	else if (strcmp(piece,"relation")==0) {return 2;}
	else {return -1;}
}

// binary search algorithm taken from alseda's notes
long binary_search (unsigned long key, Node *list, unsigned long ListLength)
{ 
	unsigned long start=0UL, AfterEnd=ListLength, middle;
	unsigned long try;

	while(AfterEnd > start)
	{ 
		middle = start + ((AfterEnd-start-1)>>1); try = list[middle].id; 
		if (key == try) return middle;
		else if ( key > try ) start=middle+1;
		else AfterEnd=middle;
	}

	return -1;
}

// stores the node on the way line in the stored vector
void store_node(char *piece) 
{
	long position = binary_search(atol(piece),nodes,nnodes);
 
	if (position != -1 ) {
		num++;
		stored = realloc(stored,sizeof(succesor)*(num));
		stored[num-1].id = atol(piece);
		stored[num-1].position = position;
	}
}

// takes the nodes in the stored vector and links them
void link_succ(int oneway)
{
	int i;
	for( i = 0; i < num-1; i++)
    {	    	
    	nodes[stored[i].position].num_succesors+=1;
    	nodes[stored[i].position].succesors = realloc(nodes[stored[i].position].succesors,sizeof(unsigned long)*nodes[stored[i].position].num_succesors);
    	nodes[stored[i].position].succesors[nodes[stored[i].position].num_succesors-1] = stored[i+1].position;

    	if (oneway==0)
    	{
    	nodes[stored[i+1].position].num_succesors++;
    	nodes[stored[i+1].position].succesors = realloc(nodes[stored[i+1].position].succesors,sizeof(unsigned long)*nodes[stored[i+1].position].num_succesors);
    	nodes[stored[i+1].position].succesors[nodes[stored[i+1].position].num_succesors-1] = stored[i].position;
    	}
    }
}

// It registers the nodes and their parameters in a node array  
void get_node (int n, char *piece)
{
	if(n == 1) {nodes[k].id =strtoul(piece, NULL, 0);}
	else if(n == 2) {nodes[k].name = piece;}
    else if (n == 9) {nodes[k].latitude = atof(piece);}
    else if (n == 10) {nodes[k].longitude = atof(piece);}
}


void reading (FILE *fcsv)
{
	int type = 0;
  	int n = 0;
	unsigned short oneway = 0;
	unsigned long third = 0;
	char *line = NULL;
	size_t size_ = 0;
	ssize_t bytes_read;
	char *piece = NULL;
  	unsigned long id_prev = 0;

	if (fcsv == NULL) exit(EXIT_FAILURE);


  int progress = 0;
	while((bytes_read = getline(&line, &size_, fcsv))!=-1)
    { 
    	if (third > 2){
		while((piece=strsep(&line,"|"))!=NULL)
    	{
    		if (n ==0)  
    		{
    			type = line_type(piece);  
    			
    			if (type!=-1) {
					if(num>1) {
						link_succ(oneway);}}
						free(stored); 
						stored = NULL; 
						num = 0; 
						oneway = 0;
			}
			if (type==0 ) {get_node(n, piece);}

			else if(type==1 ) 
			{ 
				if (n==7 && strcmp(piece,"oneway")==0) oneway=1;  
				if (n>=9) {
          			if (atol(piece) != id_prev){store_node(piece); id_prev = atol(piece);} // checks for repeating nodes
          			else {type = -1;}
          		} 
			}
			n++; 
		}
    	id_prev = 0;
		k++;
		n = 0;}
		third++;
    if (third % 2535379 == 0){
      printf ("%d %% of the csv read \n", progress);
      progress = progress + 10;}
	}
  printf("100%% of the csv read \n");
	printf("Number of lines read: %lu\n", third);
}

void writing(){
	FILE *fin;
	int i;
	
	unsigned long ntotnsucc=0UL;
	for(i=0;i<nnodes; i++) ntotnsucc+=nodes[i].num_succesors;


	if((fin=fopen("binary.bin","wb"))==NULL) printf("The binary Output file can't be opened\n");

	/* Global data---header */
	 if((fwrite(&nnodes,sizeof(unsigned long),1,fin)+
	 	fwrite(&ntotnsucc,sizeof(unsigned long),1,fin))!=2) printf("When initializing the output binary data file\n");

	 /* Writting all nodes */

	 if( fwrite(nodes,sizeof(Node),nnodes,fin)!=nnodes) printf ("When writing nodes to the output binary data file\n");

	//Writting successors in blocks
	for (i=0;i<nnodes; i++) 
		if(nodes[i].num_succesors)
		{
			if(fwrite(nodes[i].succesors,sizeof(unsigned long),nodes[i].num_succesors,fin)!=nodes[i].num_succesors)
			printf("When writing edges to the output binary data file\n");	
		}


	fclose(fin);
}



int main (){

	int i;
	FILE *fcsv;
	
	if((nodes = (Node*)malloc(nnodes*sizeof(Node)))==NULL) perror("When allocating memory for the nodes vector");
	

	// Initialization

	for(i = 0 ; i < nnodes ; i++)
	{
		nodes[i].id = 0;
		nodes[i].num_succesors = 0;
		nodes[i].succesors = malloc(sizeof(unsigned long)*1);
	}
    
	printf("Initialization finished. Please wait while the reading of the csv file and the writing of the binary file takes place \n");
	
	time_t time;
	time = clock();
	fcsv = fopen("spain.csv","r"); 
	reading(fcsv); 
	printf("Reading of the csv file finished.\n");	
	fclose(fcsv);  	
	writing();
	printf("Binary file written.\n");	
	time = clock()-time;
	printf("The elapsed time by reading and writing functions is: %fs \n",((double)time)/CLOCKS_PER_SEC);	

	return 0; 
}

