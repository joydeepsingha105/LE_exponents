/*
 *  gram_schmidt.c
 *  
 *
 *  Created by Joydeep Singha on 01/04/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 given a set of vectors along the columns of a matrix U of order 2n this program performs a gram schmidt orthogonalization process on the set of column vectors and thus replaces the columns with a set of orthogonal vectors
 */

 #include<stdio.h>
 #include<stdlib.h>
 #include<math.h>
 #include<gsl/gsl_math.h>
 #include<gsl/gsl_matrix.h>
 #include<gsl/gsl_vector.h>
 double gram_schmidt(gsl_matrix * U, int n);												/*function prototypes*/
 double norm(double *vector, int j, int n);
 double component(gsl_matrix * U, gsl_matrix * temp, double *comp, double *o_vec, int n);
 double gram_schmidt(gsl_matrix * U, int n)													/*function definition*/
 {
	int i, j, k;																			/*integers to run loops*/
	gsl_vector * U_vec;																		/*holds the non ortogonal columns*/
	gsl_vector * temp_vec;																	
	gsl_vector * new_vec;																	/*holds the orthogonal columns*/
	double absolute;																		/*variables to calculate intermediate values*/
	double value;
	U_vec = gsl_vector_alloc(2 * n);														/*allocating the gsl vectors*/
	temp_vec = gsl_vector_alloc(2 * n);
	new_vec = gsl_vector_alloc(2 * n);
	for(i = 0;i < (2 * n);i++)																/*calculating orthogonal sets of vectors*/
	{
		if(i > 0)																			/*we keep i = 0 column unchanged*/
		{
			for(j = 0;j < (2 * n);j++)														/*loading a vector or i > 0 & we will make this vector perpendicular to the previous vectors*/
			{
				gsl_vector_set(U_vec, j, gsl_matrix_get(U, j, i));
			}
			gsl_vector_memcpy(new_vec, U_vec);												/*copying the U_vec to new_vec to keep this vector safe for further calculations*/
			for(k = 0;k < i;k++)															/*looping through the vectors whose column numbers are less than that of the chosen vectors*/
			{	
				value = 0.0;																/*initializing the variables*/
				absolute = 0.0;
				for(j = 0;j < (2 * n);j++)													/*looding the vectors elements from the kth column*/
				{
					gsl_vector_set(temp_vec, j, gsl_matrix_get(U, j, k));											/*loading the kth vector*/		
					value = value + (gsl_vector_get(U_vec, j) * gsl_vector_get(temp_vec, j));						/*calculating the dot product between the kth vector and the ith vector*/
					absolute = absolute + (gsl_vector_get(temp_vec, j) * gsl_vector_get(temp_vec, j));			/*calculating the square of the norm of the kth vector*/
				}
				gsl_vector_scale(temp_vec, ((value)/(absolute)));													/*component of the ith vector along the kth vector*/
				gsl_vector_sub(new_vec, temp_vec);																	/*subtracting the component of the ith vector along the kth vector from the ith vector*/
			}
			for(j = 0;j < (2 * n);j++)																				/*replacing the vector in the ith column with the orthogonal vector*/
			{
				gsl_matrix_set(U, j, i, gsl_vector_get(new_vec, j));
			}
		}
	}
	gsl_vector_free(U_vec);																							/*freeing the variables*/
	gsl_vector_free(new_vec);
	gsl_vector_free(temp_vec);
//	free(absolute);
//	free(value);
 }
		
