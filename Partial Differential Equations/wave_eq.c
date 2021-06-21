#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include<string.h>
#include <stdbool.h>

/*   
/		This C program is created by Arturo del Cerro (UAB) in order to solve the wave equation problem 
/		for different values the user can select such as the method used to solve the system of
/		linear equations.
/
/		Contact: arturo.delcerro@e-campus.uab.cat
/       More projects: https://github.com/arturodcv
/
/
/		To execute this program execute the following commands
/					
/								gcc -Wall -O3 -o wave_eq wave_eq.c -lm       <---- To compile
/								./wave_eq                                    <---- To execute
/
/
*/


// If desired the user can modify this function values in order to see different behaviours
float  xini(float x){ //Function of x
    return sin(x);
}

float derini(float x){ // This function must be the derivative of the xini function
    return cos(x);
}

float tini(float x){ //Function of x
    return x;
}
float tfin(float x){ //Function of x
    return 0;
}

float a(float x,float t){ //Function of t and x. It could return values as x*t or x^2 + 2*t + 2
    return 1;
}

float g(float x,float t){ //Function of t and x
    return x*t;
}

void tridiagonal_matrix_algorithm( float *A, float *b,float *x, int size_){
	/* 
	This function produces a simplified Gaussian elimination algorithm in order to solve a tridiagonal system of equations 
	of the form A*x = b. 
	A tridiagonal system for n unknowns can be written as:
			up_{i}x_{i-1} + diag_{i}x_{i} + low_{i}x_{i+1} = b_{i}
	and the solution can be obtained in O(n) instead of O(n^3) that produces the Gauss -Jordan algorithm.
	As we know that our matrix is diagonal dominant, then the stability of this algorithm is garantized.
	
	More information: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	*/
	int i;
	float w;
	float *diag_, *low_, *up_, *b_;
	diag_     = (float* ) malloc( (size_) * sizeof(float) );
	low_      = (float* ) malloc( (size_ - 1) * sizeof(float) );
	up_       = (float* ) malloc( (size_ - 1) * sizeof(float) );
	b_        = (float* ) malloc( (size_) * sizeof(float) );


	for (i=0; i < size_ ;i++){						// We make copies of the main diagonal and the b term 
		diag_[i] = A[i*(size_) + i]; 
		b_[i] = b[i]; 
	}

	for (i=0;i < size_ - 1;i++){					// We make copies of the other two main diagonals
		low_[i] = A[(i+1)*(size_) + i];  
		up_[i]  = A[i*(size_) + i+1];  
	}

	for (i=1;i<size_;++i){							// Here the algorithm starts
        w = low_[i-1]/diag_[i-1]; 					// In a first step the lower diagonal is eliminated by 0s.
        diag_[i] = diag_[i] - w*up_[i-1]; 
        b_[i] = b_[i] - w*b_[i-1]; 
    }

    for (i = 0; i < size_; ++i){
    	x[i] = diag_[i];
    }

    x[size_-1] = b_[size_-1]/diag_[size_-1];

    for (i=size_-2;i>-1;--i){						// Then the backward substitution on x produces the solution
        x[i] = (b_[i] - up_[i]*x[i+1])/diag_[i];	
    }
}

void gauss_jordan(float *A, float *B, float *x, int size){
	/* This function produces the Gauss - Jordan algorithm for solving systems of linear equations 
	In our case, we are looking for a solution, x if exists, for A*x = B. 
	*/
	float *C;
	C = (float*) malloc( (size)*(size + 1)*sizeof(float) );
	float ratio;

	for (int i = 0; i < size; ++i){					// In a first term the matrix A and B are joined into the 
		for (int j = 0; j < size; ++j){				// extended matrix C of size N * (N+1)
			C[i*(size + 1) + j] = A[i*size + j]; 
		}
		C[i*(size + 1) + size ] = B[i];
	}

	for (int i = 0; i < size; ++i){
		for (int j = i+1; j < size; ++j){
			ratio = C[j*(size+1) + i] / C[i*(size+1) + i];
			for (int k = 0; k < size + 1; ++k)
			{
				C[j*(size +1) + k] = C[j*(size +1) + k] - ratio * C[i*(size + 1) + k]; //0s are made below the diagonal
			}
		}
	}  
	x[size-1] = C[(size-1)*(size+1) + size]/C[(size-1)*(size +1) + size-1] ;
	for (int i = size - 2; i > -1; --i){
		x[i] = C[i*(size +1) + size];
		for (int j = i + 1; j < size; ++j){
			x[i] = x[i] - C[i*(size + 1) + j] * x[j]; //Row operations are computed in order to give the solution
		}
		x[i] = x[i] / C[i*(size + 1) + i];
	}
}

void Exit_Error(const char *miss, int errorcode) {
    fprintf (stderr, "\n ERROR: %s.\n Stopping...\n\n", miss); exit(errorcode);
}


int main(){

int values_decission, algorithm_decission; 
float x0 ,x1 ,t0, t1, dx, dt;

printf("Do you want a set of default wave equation values? Or do you want to customize them?\n"); 
printf("\n'Easy' default values (Results in a 60x25 grid): \n\n    x0 = -3.0, x1 = 3.0 \n    t0 = 0.0, t1 = 13.0/10.0\n    dx = 1.0/10.0, dt = 1.0/20.0   \n"); 
printf("\n'More difficult' default values (Results in a 260x600 grid) :\n\n    x0 = -3.0, x1 = 3.0 \n    t0 = 0.0, t1 = 2.0\n    dx = 23.0/1000.0, dt = 0.00333   \n");  
printf("\nEasy values: 0 \nMore difficult values: 1 \nCustom: 2 \n \nPlease select: \n");

if (scanf("%d", &values_decission) != 1) Exit_Error("when reading the default or custom values for the wave equation", 0);
if (values_decission == 0){ 

	/* In case you need to hardcode the default values and the scanf isn't working, please change this below values */
	x0 = -3.0; x1 = 3.0; 
	t0 = 0.0; t1 = 1.3; 
	dx = 0.1; dt = 0.05;

}
if (values_decission == 1){ 

	x0 = -3.0; x1 = 3.0; 
	t0 = 0.0; t1 = 2.0; 
	dx = 0.023; dt = 0.00333333;

}
else if (values_decission == 2){
	printf("\n\nPlease enter the values for the next variables:");
	printf("\n(Remember that x0,x1,t0,t1,dx and dt are float types (You must write them as 3.0 or 1.0/4.0)"); 

	printf("\nx0: "); if (scanf("%f", &x0) != 1) Exit_Error("when reading the x0 value", 1);
	printf("\nx1: "); if (scanf("%f", &x1) != 1) Exit_Error("when reading the x1 value", 2);
	printf("\nt0: "); if (scanf("%f", &t0) != 1) Exit_Error("when reading the t0 value", 3);
	printf("\nt1: "); if (scanf("%f", &t1) != 1) Exit_Error("when reading the t1 value", 4);
	printf("\ndx: "); if (scanf("%f", &dx) != 1) Exit_Error("when reading the dx value", 5);
	printf("\ndt: "); if (scanf("%f", &dt) != 1) Exit_Error("when reading the dt value", 4);
	printf("\n");
}

printf("\nWhich algorithm do you want to use for solving the system of linear equations?");
printf("\nGaussian elimination method (Gauss - Jordan algorithm): 0 \nTridiagonal matrix algorithm method:                    1");
printf("\n\nPlease select:"); 
if (scanf("%d", &algorithm_decission) != 1) Exit_Error("when reading the selected algorithm", 1);
printf("\nProgram is now running! \n\n");



int M = (x1-x0)/dx, N1 = (t1-t0)/dt; 
float lamda = dt/dx*1;


// Allocate memory for the matrices we need
float *R, *V, *S, *sol; 
R    = (float* ) malloc( 1024 * (N1 +1) * (M+1) * sizeof(float) ); //This 1024 is for extra memory in case the matrix is extremly large
V    = (float* ) malloc( 1024 * (M-1) * (M-1) * sizeof(float) );
S    = (float* ) malloc( 1024 * (M-1) * sizeof(float) );
sol  = (float* ) malloc( 1024 * (M-1) * sizeof(float) );

// Prepare the output file
FILE * output;
output=fopen("output.txt","w");

//Time counters
double total_solver_time,total_output_time;
time_t total_time, solver_time, output_time; 
total_solver_time = 0.0;
total_output_time = 0.0;
total_time = clock();


//Here the wave equation program starts
int i,j,n;

for (i = 0;i < N1 + 1; ++i){ 
	for (j = 0;j < M+1; ++j){
		R[i*(N1+1)+j] = 0.0f;
	}
} 

for (i = 0;i < M-1;++i ){
	for (j=0;j<M-1;++j){
		V[i*(M-1)+j] = 0.0;
	}
	S[i] = 0.0;
}

for (i=0;i<M+1;++i){
	R[i] = xini(x0 + (i-1)*dx); 
    R[M+1+i] = R[i] + derini(x0+(i-1)*dx)*dt; 
}

for (i=0;i<N1+1;++i){
    R[i*(M+1)] = tini(t0+(i-1)*dt); 
    R[i*(M+1) + M] = tfin(t0+(i-1)*dt);
} 

for (n =1;n<N1+1;++n){
	for (i =0;i<M-1;++i){
		for (j =0;j<M-1;++j){
			V[i*(M-1) + j] = 0.0;
		}
		V[i*(M-1) + i] = 1.0;
		V[i*(M-1) + i-1] = -(pow(lamda,2)*a(x0+(-1+i)*dx,t0+(-1+n)*dt))/(1+2*pow(lamda,2)*a(x0+(-1+i)*dx,t0+(-1+n)*dt));
		V[i*(M-1) + i+1] = -(pow(lamda,2)*a(x0+(-1+i)*dx,t0+(-1+n)*dt))/(1+2*pow(lamda,2)*a(x0+(-1+i)*dx,t0+(-1+n)*dt));
	}

	for (i =1;i<M-1;++i){
		S[i-1] = -(-pow(dt,2)*g(x0+(-1+i)*dx,t0+(-1+n)*dt) + R[(-2+n)*(M+1) + i] -2*R[(-1+n)*(M+1) + i])/(1+2*pow(lamda,2)*a(x0+(-1+i)*dx,t0+(-1+n)*dt));
	}

	S[0] = S[0] +  (pow(lamda,2)*a(x0+(-1+1)*dx,t0+(-1+n)*dt) * R[n*(M+1)])/(1+2*pow(lamda,2)*a(x0+(-1+1)*dx,t0+(-1+n)*dt));
	S[M-2] = S[M-2] + (pow(lamda,2)*a(x0+(-1+1)*dx,t0+(-1+n)*dt) * R[n*(M+1) + M]) / (1+2*pow(lamda,2)*a(x0+(-1+1)*dx,t0+(-1+n)*dt));

	solver_time = clock();
	if (algorithm_decission == 0){
		gauss_jordan(V,S,sol,M-1);
	}
	else if (algorithm_decission == 1){
		tridiagonal_matrix_algorithm(V, S, sol, M-1);
	}

	solver_time = clock() - solver_time;
	total_solver_time = total_solver_time + ((double)solver_time) / CLOCKS_PER_SEC ;

	output_time = clock();
	for (i=0;i<M-2;i++){
		R[n*(M+1) + i+1] = sol[i]; 
		fprintf(output,"%f;", sol[i]);
	}
	fprintf(output, "\n" );
	output_time = clock() - output_time;
	total_output_time = total_output_time + ((double)output_time) / CLOCKS_PER_SEC ;

}

total_time = clock() - total_time;

printf("Time spent executing the program: %fs \n",  ((double)total_time) / CLOCKS_PER_SEC );
printf("Time spent solving the linear system of equations: %fs \n", total_solver_time);
printf("Time spent writing the output file: %fs \n",total_output_time );

printf("\n\nProgram finished correctly! Check the output.txt file for results or execute the 'python3 results.py' command to see an animation of this wave equation. \n\n");

//Close our output file
fclose(output);

return 0;
}