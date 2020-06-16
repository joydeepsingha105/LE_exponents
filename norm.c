/*
 *  norm.c
 *  
 *
 *  Created by Joydeep Singha on 02/04/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 this code normalizes a vector array
 */
 #include<stdio.h>
 #include<stdlib.h>
 #include<math.h>
 double norm(double *vector, int n);				/*function prototype*/
double norm(double *vector, int n)					/*function definition*/
 {	
	int i;											/*variable to run loop*/
	double value;									/*variable to calculate intermediate values*/
	value = 0.0;
	for(i = 0;i < (2 * n);i++)						/*squaring each element and summing*/
	{
		value = value + (vector[i] * vector[i]);
	}
	value = pow(value, 0.5);						/*norm of the vector*/
	for(i = 0;i < (2 * n);i++)						/*deving each element of the vector by the norm*/
	{
		vector[i] = vector[i]/(value);
	}
//	free(value);									/*freeing the variables*/
 }

