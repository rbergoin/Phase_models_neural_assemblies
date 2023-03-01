#include "graph.h"


/*
* Change the input vector by null values
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
*/
void addNullInputs(double *inputs, int n)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		inputs[i] = 0.0;
	}
}


/*
* Change the input vector by binary values on a specific cluster
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
*/
void addBinaryLocalized(double *inputs, int n, int index, int nbInput)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)))
		{
			inputs[i] = 3.0;
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}


/*
* Change the input vector by binary values on a specific cluster with random activation
*
* @param	inputs			the input vector
* @param	n				the number of input in the vector
* @param	index			the index of the center position the cluster
* @param	nbInput			the number of input to activate
*/
void addBinaryLocalizedRandom(double *inputs, int n, int index, int nbInput)
{
	int i;
	
	//For all inputs
	for(i=0; i < n; i++)
	{
		//If we are on the cluster
		if((i>=(index-nbInput/2)) && (i<(index+nbInput/2)))
		{
			if((rand()%2)!=0)
			{
				inputs[i] = 3.0;
			}
			else
			{
				inputs[i] = 0.0;
			}
		}
		else
		{
			inputs[i] = 0.0;
		}
	}
}

