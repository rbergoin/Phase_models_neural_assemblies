#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>


/*
* Gaussian function
*
* @param	x		variable
* @param	mu		mean of the gaussian
* @param	sigma	standard deviation of the gaussian
* @param	ampl	amplitude of the gaussian
*/
double gaussian(double x, float mu, float sigma, float ampl)
{
	return ampl * 1.0/(sigma*sqrt(2.0*M_PI)) * exp(-0.5*pow((x-mu)/sigma, 2.0));
}


/*
* Calculate the euclidean distance between two points
*
* @param	x1		abscissa of point 1
* @param	y1		ordinate of point 1
* @param	x2		abscissa of point 2
* @param	y2		ordinate of point 2
*/
double calculate_distance(float x1, float y1, float x2, float y2) 
{ 
	return sqrt(pow((x2-x1), 2.0) + pow((y2-y1), 2.0));
}


/*
* Calculate the euclidean distance between two 3D points
*
* @param	x1		abscissa of point 1
* @param	y1		ordinate of point 1
* @param	z1		ordinate of point 1
* @param	x2		abscissa of point 2
* @param	y2		ordinate of point 2
* @param	z2		ordinate of point 2
*/
double calculate_distance_3D(float x1, float y1, float z1, float x2, float y2, float z2) 
{ 
	return sqrt(pow((x2-x1), 2.0) + pow((y2-y1), 2.0) + pow((z2-z1), 2.0));
}


/*
* Generate an random number uniformly
*
* @param	min		range min
* @param	max		range max
*/
double get_random(double min, double max) 
{ 
	return ((double) rand() / ((double)RAND_MAX)) * (max-min) + min;
}


/*
* Generate an random number distributed with a normal/gaussian law
*
* @param	mu		mean
* @param	sigma	standard deviation
* @param	min		range min
* @param	max		range max
*/
double get_random_normal(double mu, double sigma, double min, double max)
{
	double u, r, theta;           //Variables for Box-Muller method
	double x;                     //Normal(0, 1) rv
	double norm_rv;               //The adjusted normal rv

    //Generate u
	u = 0.0;
	while (u == 0.0)
	{
		u = get_random(0.0, 1.0);
	}

	//Compute r
	r = sqrt(-2.0 * log(u));

	//Generate theta
	theta = 0.0;
	while (theta == 0.0)
	{
		theta = 2.0 * M_PI * get_random(0.0, 1.0);
	}

    //Generate x value
    x = r * cos(theta);

    //Adjust x value for specified mean and variance
    norm_rv = (x * sigma) + mu;
	
	if ((norm_rv<=max) && (norm_rv>=min))
	{
	    //Return the normally distributed RV value
	    return norm_rv;
	}
	else
	{
		return get_random_normal(mu, sigma, min, max);
	}
}


/*
* Calculate the degree of a graph 
*
* @param	a		adjacency matrix
* @param	n	 	number total of node
*/
int graphDegree(int **a, int n)
{
	int i, j, degree = 0;
	
	for(i=0; i<n; i++)		//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(a[i][j]==1)		//if there is a link
			{
				degree++;
			}
		}
	}
	
	return degree;
}


/*
* Calculate the degree of a node 
*
* @param	i		index of the node
* @param	a		adjacency matrix
* @param	n	 	number total of node
*/
int nodeDegree(int i, int **a, int n)
{
	int j, degree = 0;

	for(j=0; j<n; j++)		//for each nodes
	{
		if(a[i][j]==1)		//if there is a link
		{
			degree++;
		}
	}
	
	return degree;
}


/*
* Initialize the adjacency matrix with random regular connections
*
* @param	a				adjacency matrix
* @param	K				degree of connection of the a(incoming degree)
* @param	n	 			number of node
* @param	ratio			ratio of excitatory neurons
*/
void createRandomRegular(int **a, int K, int n, int ratio)
{
	int i, j, k;
	
	//excitatory to excitatory connections
	for(i=0; i<n/100.0*ratio; i++)		//for each nodes
	{
		k=1;
		do
		{
			j = rand()%(int)(n/100.0*ratio);	//get random node
			if((a[i][j]==0) && (i!=j))	//if it is not already connected
			{
				a[i][j] = 1;
				k++;
			}
		} while(k<= (K*ratio/100.0) );
	}
	
	//inhibitory to excitatory connections
	for(i=0; i<n/100.0*ratio; i++)		//for each nodes
	{
		/*for(j=n/100.0*ratio; j<n; j++)		//for each nodes
		{
			if(i!=j)	
			{
				a[i][j] = 1;
			}
		}*/
		
		k=1;
		do
		{
			j = n/100.0*ratio + rand()%(int)(n/100.0*(100-ratio));	//get random node
			if((a[i][j]==0) && (i!=j))	//if it is not already connected
			{
				a[i][j] = 1;
				k++;
			}
		} while(k<=(K*(100-ratio)/100.0));
	}
	
	//excitatory to inhibitory connections
	for(i=n/100.0*ratio; i<n; i++)		//for each nodes
	{
		/*for(j=0; j<n/100.0*ratio; j++)		//for each nodes
		{
			if(i!=j)	
			{
				a[i][j] = 1;
			}
		}*/
		
		k=1;
		do
		{
			j = rand()%(int)(n/100.0*ratio);	//get random node
			if((a[i][j]==0) && (i!=j))	//if it is not already connected
			{
				a[i][j] = 1;
				k++;
			}
		} while(k<=(K*ratio/100.0));
	}
	
	//inhibitory to inhibitory connections
	for(i=n/100.0*ratio; i<n; i++)		//for each nodes
	{
		/*for(j=n/100.0*ratio; j<n; j++)		//for each nodes
		{
			if(i!=j)	
			{
				a[i][j] = 1;
			}
		}*/
		
		k=1;
		do
		{
			j = n/100.0*ratio + rand()%(int)(n/100.0*(100-ratio));	//get random node
			if((a[i][j]==0) && (i!=j))	//if it is not already connected
			{
				a[i][j] = 1;
				k++;
			}
		} while(k<=(K*(100-ratio)/100.0));
	}
}


/*
* Initialize the adjacency matrix with a small world connection
*
* @param	a				adjacency matrix
* @param	K				degre of connection of the a(incoming degree)
* @param	p				probability to change a link(rewiring)
* @param	n	 			number of node
*/
void createSmallWorld(int **a, int K, float p, int n)
{
	int i, j, newJ;
	
	//Connect to nearest K nodes
	for(i=0; i<n; i++)		//for each nodes
	{
		for(j=(i/K)*K; j<(i/K)*K+K; j++)
		{
			if((i!=j) && (i<(n/K)*K))
			{
				a[i][j] = 1;
			}
		}
	}
	
	//Create new links
	for(i=0; i<n; i++)		//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(a[i][j] == 1)		//if theres is a link
			{
				if (get_random(0, 1) < p)
				{
					newJ = rand()%((n/K)*K);
					if((a[i][newJ] == 0) && (i!=newJ))		//if it is not already connected
					{
						a[i][j] = 0;
						a[i][newJ] = 1;
					}
				}
			}
		}
	}
}


/*
* Initialize the adjacency matrix with a random small world connection
*
* @param	a				adjacency matrix
* @param	K				degre of connection of the a(incoming degree)
* @param	p				probability to change a link(rewiring)
* @param	n	 			number of node
*/
void createSmallWorldRandom(int **a, int K, float p, int n)
{
	int i, j, newJ;
	
	//Connect to nearest nodes randomly
	for(i=0; i<n; i++)		//for each nodes
	{
		for(j=(i/K)*K; j<(i/K)*K+K; j++)
		{
			if((i!=j) && (i<(n/K)*K))
			{
				a[i][j] = rand()%2;
			}
		}
	}
	
	//Create new links
	for(i=0; i<n; i++)		//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(a[i][j] == 1)		//if theres is a link
			{
				if (get_random(0, 1) < p)
				{
					newJ = rand()%((n/K)*K);
					if((a[i][newJ] == 0) && (i!=newJ))		//if it is not already connected
					{
						a[i][j] = 0;
						a[i][newJ] = 1;
					}
				}
			}
		}
	}
}


/*
* Initialize the adjacency matrix with a scale free connection (Barabási–Albert model)
*
* @param	a				adjacency matrix
* @param	K				degre of connection of the a(incoming degree)
* @param	m0				the initial number of nodes
* @param	n	 			number of node
*/
void createScaleFree(int **a, int K, int m0, int n)
{
	int i, j;
	double p;
	
	//Connect to nearest K nodes
	for(i=0; i<m0; i++)		//for each initial nodes
	{
		for(j=0; j<K; j++)		//for nearest K nodes
		{
			a[i][j] = 1;
		}
	}
	
	for(j=m0; j<n; j++)		//for each new node/hub
	{
		for (i=0; i<j; i++) 
		{
			p = (double)nodeDegree(i, a, n)/graphDegree(a, n);
			if (get_random(0, 1) < p)
			{
				a[j][i] = 1;
				a[i][j] = 1;
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a scale free random connection (Barabási–Albert model)
*
* @param	a				adjacency matrix
* @param	m0				the initial number of nodes
* @param	n	 			number of node
*/
void createScaleFreeRandom(int **a, int m0, int n)
{
	int i, j;
	double p;
	
	//Connect to random nearest nodes
	for(i=0; i<m0; i++)		//for each initial nodes
	{
		for(j=0; j<m0; j++)
		{
			a[i][j] = rand()%2;
		}
	}
	
	for(j=m0; j<n; j++)		//for each new node/hub
	{
		for (i=0; i<j; i++) 
		{
			p = (double) nodeDegree(i, a, n)/graphDegree(a, n);
			if (get_random(0, 1) < p)
			{
				a[j][i] = 1;
				a[i][j] = 1;
			}
		}
	}	
}

/*
* Initialize the adjacency matrix with a modular connection
*
* @param	a				adjacency matrix
* @param	K				degre of connection of the a(incoming degree)
* @param	n	 			number of node
*/
void createModular(int **a, int K, int n)
{
	int i, j;

	
	//Connect to nearest nodes randomly
	for(i=0; i<n; i++)		//for each nodes
	{
		for(j=(i/K)*K; j<(i/K)*K+K; j++)
		{
			if((i!=j) && (i<(n/K)*K))
			{
				a[i][j] = rand()%2;
			}
		}
	}
	
	for(i=1; i<n; i+=(K/2)+1)		//connect hubs from each other
	{
		for (j=1; j<n; j+=(K/2)+1) 			
		{
			a[i][j] = 1;
		}
	}	
}


/*
* Initialize the adjacency matrix with a multi layer connection
*
* @param	a				adjacency matrix
* @param	nbLayer			number of layer
* @param	n	 			number of node
*/
void createMultiLayer(int **a, int nbLayer, int n)
{
	int i, j;

	for(j=0; j<n; j++)		//for each nodes
	{
		for(i=(j/(n/nbLayer)+1)*(n/nbLayer); i<(j/(n/nbLayer)+1)*(n/nbLayer)+(n/nbLayer); i++)
		{
			if((i<n))
			{
				a[i][j] = 1;
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a multi layer connection with feedback
*
* @param	a				adjacency matrix
* @param	nbLayer			number of layer
* @param	n	 			number of node
*/
void createMultiLayerFeedback(int **a, int nbLayer, int n)
{
	int i, j;

	for(j=0; j<n; j++)		//for each nodes
	{
		for(i=(j/(n/nbLayer)+1)*(n/nbLayer); i<(j/(n/nbLayer)+1)*(n/nbLayer)+(n/nbLayer); i++)
		{
			if((i<n))
			{
				a[i][j] = 1;
				a[j][i] = 1;
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with all input linked to a fully connected topologie
*
* @param	a				adjacency matrix
* @param	nbInput			number of input node
* @param	n	 			number of node
*/
void createAllToFull(int **a, int nbInput, int n)
{
	int i, j;
	
	for(i=0; i<n; i++)			//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(i!=j)
			{
				if((i<nbInput) && (j>=nbInput))
				{
					a[i][j] = 1;
				}
				else if((j<nbInput) && (i>=nbInput))
				{
					a[i][j] = 1;
				}
				else if((j>=nbInput) && (i>=nbInput))
				{
					a[i][j] = 1;
				}
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a map topology number of neurons n = h*w
*
* @param	a				adjacency matrix
* @param	h	 			height of the map
* @param	w	 			width of the map
*/
void createMapTopology(int **a, int h, int w)
{
	int i, j;
	
	for(i=0; i<h; i++)			//for each nodes
	{
		for(j=0; j<w; j++)		//for each nodes
		{
			if((j+1)<w)
			{
				a[w*i + j][w*i + (j+1)] = 1;
			}
			if((j-1)>=0)
			{
				a[w*i + j][w*i + (j-1)] = 1;
			}
			if((i+1)<h)
			{
				a[w*i + j][w*(i+1) + j] = 1;
			}
			if((i-1)>=0)
			{
				a[w*i + j][w*(i-1) + j] = 1;
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a reservoir computing connection
*
* @param	a				adjacency matrix
* @param	nbInput			number of input node
* @param	nbOutput		number of output node
* @param	n	 			number of node total
*/
void createReservoir(int **a, int nbInput, int nbOutput, int n)
{
	int i, j;

	for(j=0; j<n; j++)		//for each nodes
	{
		for(i=0; i<n; i++)		//for each nodes
		{
			if(i<nbInput)				//If node i is an input
			{
				if((j>=nbInput) && (j < (n-nbOutput)))		//If node j is in the reservoir
				{
					a[i][j] = 1;
				}
			}
			else if(i>=(n-nbOutput))	//If node i is an output
			{
				if((j>=nbInput) && (j < (n-nbOutput)))		//If node j is in the reservoir
				{
					a[i][j] = 1;
				}
			}
			else						//If node i is in the reservoir
			{
				a[i][j] = 1;
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a reservoir computing connection
*
* @param	a				adjacency matrix
* @param	nbInput			number of input node
* @param	nbOutput		number of output node
* @param	n	 			number of node total
*/
void createReservoir2(int **a, int nbInput, int nbOutput, int n)
{
	int i, j;
	int sizeGroup = (n-(nbInput + nbOutput))/(nbInput + nbOutput);

	for(j=0; j<n; j++)		//for each nodes
	{
		for(i=0; i<n; i++)		//for each nodes
		{
			if(i<nbInput)				//If node i is an input
			{
				if((j>= (nbInput+i*sizeGroup)) && (j < (nbInput+(i+1)*sizeGroup)))		//If node j is in his group
				{
					a[i][j] = 1;
					a[j][i] = 1;
				}
			}
			else if(i>=(n-nbOutput))	//If node i is an output/readout
			{
				if((j>= (nbInput+nbInput*sizeGroup +(i-n+2)*sizeGroup)) && (j < (nbInput+nbInput*sizeGroup +(i+1-n+2)*sizeGroup)))		//If node j is in his group
				{
					a[i][j] = 1;
					a[j][i] = 1;
				}
			}
			else if((i>=nbInput) && (i<nbInput+nbInput*sizeGroup))						//If node i is in the hidden layer
			{
				if((j>=nbInput+nbInput*sizeGroup) && (j<(n-nbOutput)))						//If node j is in the analog ouput layer
				{
					a[i][j] = 1;
					a[j][i] = 1;
				}
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a random topology depending on physical distances
*
* @param	a				adjacency matrix
* @param	n	 			number of node
* @param	w				width of the 2D array
* @param	inhibitory		if neurons are inhibitory or excitatory
*/
void createPhysicalTopologyDependant(int **a, int n, int w, int *type_neuron)
{
	int i, j;
	
	double distance;
	double proba;
	
	for(i=0; i<n; i++)			//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(i!=j)
			{
				distance = calculate_distance(i%w, i/w, j%w, j/w);	//distance max 12.72 for 10x10
				proba = get_random(0.0, 1.0);
				if(type_neuron[j]!=0)
				{
					if(proba<gaussian(distance, 1.0, 1.0, 2.5))
					{
						a[i][j] = 1;
					}
				}
				else
				{
					if(proba<gaussian(distance, 2.0, 2.0, 5.0))
					{
						a[i][j] = 1;
					}
				}
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a random topology depending on physical distances
*
* @param	a				adjacency matrix
* @param	n	 			number of node
* @param	w				width of the 2D array
* @param	inhibitory		if neurons are inhibitory or excitatory
*/
void createPhysicalTopologyDependant_2(int **a, int n, int w, int *type_neuron)
{
	int i, j;
	
	double distance;
	double proba;
	
	for(i=0; i<n; i++)			//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(i!=j)
			{
				proba = get_random(0.0, 1.0);
				if((type_neuron[j]!=0) || (type_neuron[i]!=0))
				{
					distance = 0;
					a[i][j] = 1;
				}
				else
				{
					distance = calculate_distance(i%w, i/w, j%w, j/w);
					if(proba<gaussian(distance, 1.0, 4.0, 10.0))
					{
						a[i][j] = 1;
					}
				}
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a random topology depending on physical distances
*
* @param	a				adjacency matrix
* @param	n	 			number of node
* @param	w				width of the 2D array
* @param	type_neuron		type of the neuron
*/
void createPhysicalTopologyDependant_3(int **a, int n, int w, int *type_neuron)
{
	int i, j;
	
	double distance;
	double proba;
	
	for(i=0; i<n; i++)			//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(i!=j)
			{
				distance = calculate_distance(i%w, i/w, j%w, j/w);	//distance max 12.72 for 10x10
				proba = get_random(0.0, 1.0);
				if(type_neuron[j]!=0)
				{
					if(proba<gaussian(distance, 1.0, 2.0, 5.0)) // 3 voisinage
					{
						a[i][j] = 1;
					}
				}
				else
				{
					if(proba<=0.5)
					{
						a[i][j] = 1;
					}
				}
			}
		}
	}	
}


/*
* Initialize the adjacency matrix with a random topology depending on physical distances in 3D
*
* @param	a				adjacency matrix
* @param	n	 			number of node
* @param	w				width of the 3D array
* @param	h				height of the 3D array
* @param	inhibitory		if neurons are inhibitory or excitatory
*/
void createPhysicalTopologyDependant3D(int **a, int n, int w, int h, bool *inhibitory)
{
	int i, j;
	
	double distance;
	double proba;
	
	for(i=0; i<n; i++)			//for each nodes
	{
		for(j=0; j<n; j++)		//for each nodes
		{
			if(i!=j)
			{
				distance = calculate_distance_3D(i%w, (i/w)%h, i/(w*h), j%w, (j/w)%h, j/(w*h));	//distance max 6.92 for 5x5x5
				proba = get_random(0.0, 1.0);
				if(inhibitory[j])
				{
					if(proba<gaussian(distance, 1.0, 1.0, 2.5))
					{
						a[i][j] = 1;
					}
				}
				else
				{
					if(proba<gaussian(distance, 2.0, 2.0, 5.0))
					{
						a[i][j] = 1;
					}
				}
			}
		}
	}	
}
