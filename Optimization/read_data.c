#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include "my_functions.h"

int Data_from_txt(FILE *ftxt,DataForFitting *TheData){
	char *line = NULL;
	size_t size_ = 0;
	ssize_t bytes_read;
	char *piece = NULL;
	int nline = 0; 
	int npiece = 0;
	
	while((bytes_read = getline(&line, &size_, ftxt))!=-1)
    { 
    	while((piece=strsep(&line," "))!=NULL)
    	{
    		if (npiece>0){ TheData->Data_Time_Series[nline][npiece-1] = atof(piece); }
    		npiece++;
    	}
    	nline++;
    	npiece = 0;
	}
  TheData->N_days = nline-1; 
  TheData->PopSize = 100; 
  return 1;
}

void Get_data(DataForFitting *TheData){
   FILE *ftxt;
   if ((ftxt = fopen( "Genetic_data.txt" , "r")) == NULL)
		 Exit_Error("The genetic_data.txt does not exist or cannot be opened\n",2);
	 printf("File opened, starting to read the txt file.\n");
	 if (Data_from_txt(ftxt,TheData) != 1) Exit_Error("Data was not read correctly.\n",2);
   printf("Data from txt file was read correctly.\n ");
}

void print_data(DataForFitting *TheData){
   int i, j; printf("Printing data:\n");
 	 for (i=0; i<Ndays_Time_Series; i++){
	 	  for (j=0; j<Nvar_Time_Series; j++){
	 		  printf("%f, ", TheData->Data_Time_Series[i][j]);}
   printf("\n");
   }
   printf("N_days = %lu\n", TheData->N_days);
   printf("PopSize = %f\n", TheData->PopSize);
}





