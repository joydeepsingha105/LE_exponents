/*
 *  ellipse.c
 *  
 *
 *  Created by Joydeep Singha on 31/03/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_matrix.h>
double matrix(double *x_1, double *x_2, gsl_matrix * data, double K, double omega, double eps_1, int n);
double matrix(double *x_1, double *x_2, gsl_matrix * data, double K, double omega, double eps_1, int n)
{
	int i, j;
	const double pi = 4.0 * atan(1.0);
	double eps_2;
	
	gsl_matrix * A;
	gsl_matrix * B;
	gsl_matrix * C;
	gsl_matrix * D;
	
	double *f;
	double *g;

	A = gsl_matrix_alloc (n, n);
	B = gsl_matrix_alloc (n, n);
	C = gsl_matrix_alloc (n, n);
	D = gsl_matrix_alloc (n, n);
	eps_2 = 1.0 - eps_1;	
	f = (double *) malloc(n * sizeof(double));
	if(f == NULL)
	{
		printf("cannot allocate memory for f \n");
		exit(1);
	}
	
	g = (double *) malloc(n * sizeof(double));
	if(g == NULL)
	{
		printf("cannot allocate memory for g \n");
		exit(1);
	}
	
	for(i = 0;i < n;i++)
	{
		f[i] = 1.0 - (K * cos(2.0 * pi * x_1[i]));
		g[i] = 1.0 - (K * cos(2.0 * pi * x_2[i]));
	}

	for(i = 0;i < n;i++)
	{
		for(j = 0;j < n;j++)
		{
			if(i == j)
			{
				//diagonal elements of A
				gsl_matrix_set (A, i, j, ((2 - ((n - 1) * eps_1/n) - eps_2) * f[i]));
				
				//diagonal elements of D
				gsl_matrix_set (D, i, j, ((2 - ((n - 1) * eps_1/n) - eps_2) * g[i]));
			}
			
			if(i != j)
			{
				//off diagonal elements of A
				gsl_matrix_set (A, i, j, (eps_1 * f[j]/n));
				
				//off diagonal elements of D
				gsl_matrix_set (D, i, j, (eps_1 * g[j]/n));
			}
			
			//elements of B
			gsl_matrix_set (B, i, j, (eps_2 * g[j]/n));
			
			//elements of C
			gsl_matrix_set (C, i, j, (eps_2 * f[j]/n));
		}
	}

	for(i = 0;i < (2 * n);i++)
	{
		for(j = 0;j < (2 * n);j++)
		{
			//block A
			if((i < n) && (j < n))
			{
				gsl_matrix_set (data, i, j, gsl_matrix_get (A, i, j));
			}
			
			//block B
			if((i < n) && (j > (n - 1)) && (j < (2 * n)))
			{
				gsl_matrix_set (data, i, j, gsl_matrix_get (B, i, j - n));
			}
			
			//block C
			if((i > (n - 1)) && (i < (2 * n)) && (j < n))
			{
				gsl_matrix_set (data, i, j, gsl_matrix_get (C, i - n, j));
			}
			
			//block D
			if((i > (n - 1)) && (j > (n - 1)) && (i < (2 * n)) && (i < (2 * n)))
			{
				gsl_matrix_set (data, i, j, gsl_matrix_get (D, i - n, j - n));
			}
		}
	}
	
	free(f);
	free(g);
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(C);
	gsl_matrix_free(D);

}

