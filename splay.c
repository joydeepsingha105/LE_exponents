/*
 *  splay.c
 *  
 *
 *  Created by Joydeep Singha on 31/03/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 This code calculates teh lyapunov exponent for a system of two coupled sine circle map lattices.
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

double matrix(double *x_1, double *x_2, gsl_matrix * data, double K, double omega, double eps_1, int n);	/*function prototype that calculate the jacobian*/
double iterate (double K, double omega, double eps_1, double *x_1, double *x_2, int n);						/*function prototype iterates the system*/
double modulo_1(double *phase);																				/*function prototype to perform the modulus 1 operation*/
double ellipse( gsl_matrix * data, gsl_matrix * U, int n);													/*function prototype to multiply the jacobian with a unit orthogonal set of vectors*/
double normalize(gsl_matrix * U, int n);																	/*function prototype to normalize a set of vector*/
double gram_schmidt(gsl_matrix * U, int n);																	/*function prototype to perform the gram schmidt orthogonalization on a set of vectors*/
double norm(double *vector, int n);																			/*function prototype to finc out a unit vector given a vector*/
//double maximum(double *array, int n);                                                        /*function prototype to find the largest number from a given array of numbers*/
float ran2(long *idum);
int main (void)
{
    FILE *lyap;                                                                                         /*file pointer to write the largest lyapunov exponent for a value of the parameter K*/
	FILE *group_1;																							/*file pointer to write the intial condition*/
//	double max;
	int n, i, j, k, l;																						/* n = number of maps in each group. i, j, k, l= intergers to run loops*/
	double *x_1, *ini_1, *ini_2;																							/*holds the phases of group 1 maps*/
	double *x_2;																							/*holds the phases of group 2 maps*/
	double K, omega, eps_1, eps_2;																					/*parameterers in our system*/
	int max_iter;																							/*maximum number of iterations*/
    int trans;                                                                                              /*number of iterations at the begining upto which we allow the transients to go away*/
	gsl_matrix * data;																						/*gsl matrix which holds the jacobian matrix for a state*/	
	gsl_matrix * U;																							/*gsl matrix which holds the unit orthogonal vectors*/
	double U_vec;																							/*variables needed for intermediate calculations*/
	double value;
    double eps_2_max, eps_2_min;                                                                                    /*maximum and minimum values of the parameter K*/
    double d_eps_2;                                                                                              /*minimum interval to increment the values of K*/
	double *lyapunov;																						/*an array hold the values of lyapunov exponent*/
	group_1 = fopen("group_1.txt","w");																		/*opening files*/
	lyap = fopen("lyapunov_eps_2.txt","w");
	max_iter = 3010000;
	trans = 3000000;
	n = 150;
	data = gsl_matrix_alloc ((2 * n), (2 * n));																/*allocating the gsl matrices*/
	U = gsl_matrix_alloc ((2 * n), (2 * n));
	x_1 = (double *) malloc(n * sizeof(double));															/*dynamic allocation of the arrays*/
	x_2 = (double *) malloc(n * sizeof(double));
	ini_1 = (double *) malloc(n * sizeof(double));															/*dynamic allocation of the arrays*/
	ini_2 = (double *) malloc(n * sizeof(double));
	lyapunov = (double *) malloc((2 * n) * sizeof(double));
    
    eps_2_min = 0.2;
    eps_2_max = 0.4;
	
	long idum;
	double ran;
	idum = -999000000;

	
	if(x_1 == NULL)
	{
		printf("cannot allocate memory for x_1\n");
		exit(1);
	}
	if(x_2 == NULL)
	{
		printf("cannot allocate memroy for x_2\n");
		exit(1);
	}
	if(ini_1 == NULL)
	{
		printf("cannot allocate memory for x_1\n");
		exit(1);
	}
	if(ini_2 == NULL)
	{
		printf("cannot allocate memroy for x_2\n");
		exit(1);
	}
	if(lyapunov == NULL)
	{
		printf("cannot allocate memory for lyapunov\n");
	}
	for(i = 0;i < (4 * n);i++)																					/*setting the initial condition*/
	{
		ran = ran2(&idum);
	}
	for(i = 0;i < n;i++)																					/*setting the initial condition*/
	{
		ran = ran2(&idum);
		ini_1[i] = ran; //(double)i/(2 * n);
			
		ran = ran2(&idum);
		ini_2[i] = ran;
	}

    K = pow(10.0,-6);
	omega = 0.27;
    d_eps_2 = (eps_2_max - eps_2_min)/30;
    for(l = 0;l < 30;l++)                                                                                      /*running loop over the range of values of K*/
    {
//		max = 0.0;
        eps_2 = eps_2_min + (l * d_eps_2);
		eps_1 = 1.0 - eps_2;
		for(i = 0;i < (2 * n);i++)																				/*running loops to create the unit orthogonal set of vectors*/
		{
			lyapunov[j] = 0.0;																					/*intiallizing the array which holds the values of lyapunov exponents*/
			for(j = 0;j < (2 * n);j++)
			{
				if(i == j)
				{
					gsl_matrix_set(U, i, j, 1.0);																/*setting values to the elements of the matrix U which holds the unit orhogonal vectors columnwise*/
				}
				else if( i != j)
				{
					gsl_matrix_set(U, i, j, 0.0);
				}
			}
		}

        for(i = 0;i < n;i++)																					/*setting the initial condition*/
        {
            x_1[i] = ini_1[i]; //(double)i/(2 * n);
			
            x_2[i] = ini_2[i];
        }
		
        for(j = 0;j < max_iter;j++)																				/*iterating the system with a given initial condition and parameter values*/
        {
            iterate(K, omega, eps_1, x_1, x_2, n);																/*calling the function which iterates the system to one step*/
        //    printf("K = %e\t%d\n",K, j);
            if(j >= trans)
            {
                matrix(x_1, x_2, data, K, omega, eps_1, n);															/*calling the function to calculate the jacobian*/
                ellipse(data, U, n);																				/*multiplying the jacobian to U which contains normalized vectors in its column*/
                gram_schmidt(U, n);																					/*finding a new set of orthogonal vectors using gram schmidt orthogonalization process*/
                /*things need to done here*/
                for(i = 0;i < (2 * n);i++)																			/*calculating the lyapunov exponent by reading each column from the matrix U and */
                {
                    value = 0.0;
                    for(k = 0;k < (2 * n);k++)																		/*reading out a column the matrix in each run of the loop and calculating the norm*/
                    {
                        U_vec = gsl_matrix_get(U, k, i);
                        value = value + (U_vec * U_vec);
                    }
                    value = pow(value, 0.5);
                    lyapunov[i] =  lyapunov[i] + log(value);														/*calculating the ith lyapunov exponent*/
                }
                normalize(U, n);																					/*calling a normalizing the column vectors in the matrix U*/
            }
        }
        for(j = 0;j < (2 * n);j++)
        {
            lyapunov[j] = lyapunov[j]/(max_iter-trans);
		//	printf("%e\n",lyapunov[j]);
        }
//		max = lyapunov[0];
//		for(j = 0;j < (2 * n);j++)
//		{
//			if((max < lyapunov[j])||(max == lyapunov[j]))
//			{
//				max = lyapunov[j];
//			}
//		}
//        value = maximum(lyapunov, n);                                                                          /*calling a function to find the largest lyapunov exponent*/
        printf("l = %d \t eps_1 = %e \n", l, eps_1); 													 /*printing the largest lyapunov exponent terminal and file*/
        for(j = 0;j < (2 * n);j++)
		{
			fprintf(lyap,"%e\t%e\n",eps_1, lyapunov[j]);
		}
		

    }
		
		
	for(i = 0;i < n;i++)																					/*writing finals states to an output file after all iterations*/
	{		
		fprintf(group_1,"%d\t%d\t%e\t%e\n",j,i,x_1[i],x_2[i]);
	}
	fprintf(group_1,"\n");
	

	free(x_1);																								/*freeing all the variables*/
	free(x_2);
	free(ini_1);
	free(ini_2);
	free(lyapunov);
	fclose(group_1);
	fclose(lyap);
	gsl_matrix_free(data);
	gsl_matrix_free(U);
	return(0);
}
//double maximum(double *array, int n)                                                               /*function definition to find the largest number from a given array of numbers*/
//{
//    double max;
//    int i;
//    max = array[0];
//    for(i = 0;i < (2 * n);i++)
//    {
//        if(max < array[i]);
//        {
//            max = array[i];
//       }
//    }
//    return(max);
//}
