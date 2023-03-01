#include "input.h"


/*
* Struct to represent a network of N theta neurons 
*
*/
struct neurons 
{
    int n;								//number of neurons
	float ratio;						//ratio of excitatory neurons
    float epsilon1;						//the dynamic of the neurons for slow adaptation
	float epsilon2;						//the dynamic of the neurons for fast adaptation
	float dt;							//integration time step
	
	int **a;							//adjacency matrix
	long double **k;					//coupling weights matrix
	double g; 							//coupling strength
	bool *inhibitory;					//if neurons are inhibitory or excitatory
	
	long double *theta;					//phases of the neurons
	float *eta;							//bifurcation parameter of the neurons
	double *inputs;						//inputs of the neurons
};


/*
* Copy values of a table in a new one
*
* @param	table	table to copy
* @param	n		size of the table
*/
long double *copyTable(long double *table, int n) 
{ 
	int i;
	long double *copy = (long double *) malloc(sizeof(long double) * n);
	
	for(i=0; i<n; i++)
	{
		copy[i] = table[i];
	}
	
	return copy;
}


/*
* Copy values of a matrix in a new one
*
* @param	matrix	matrix to copy
* @param	n		size of the matrix
*/
long double **copyMatrix(long double **matrix, int n) 
{ 
	int i, j;
	long double **copy = (long double **) malloc(sizeof(long double *) * n);

	for (i=0; i<n; i++) 
	{
		copy[i] = (long double *) malloc(sizeof(long double) * n);
		for (j=0; j<n; j++) 
		{
			copy[i][j] = matrix[i][j];
		}
	}
	
	return copy;
}


/*
* Copy values of a matrix in a new one
*
* @param	matrix	matrix to copy
* @param	n		size of the matrix
*/
void freeMatrix(long double **matrix, int n)  
{ 
	int i;
	
	//Free matrix
	for (i=0; i<n; i++) 
	{
		free(matrix[i]);
	}
	free(matrix);
}


/*
* Function to save a spike of a neuron
*
* @param	i			index of the neuron
* @param	t_spike		time of the spike in second
*/
void saveSpike(int i, long double t_spike)
{
	/**** Save spikes  through the time ****/
	FILE *fptr = fopen("spikes.txt","a");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	fprintf(fptr,"%d %3.5Lf\n", i, t_spike);
	
	fclose(fptr);
}


/*
* Initialize a network of neurons
*
* @param	n					number of neurons
* @param	epsilon1			the dynamic of the neurons for slow adaptation
* @param	epsilon2			the dynamic of the neurons for fast adaptation
* @param	dt					integration time step
* @param	g 		 			the coupling strength
* @param	adjacencyPolicy		type of policy for the adjacency matrix creation
* @param	weightPolicy		type of policy for the coupling weights matrix creation
* @param	etaPolicy			type of policy for the bifurcation parameter
* @param	phasePolicy			type of policy for the phases array creation
* @param	ratio				ratio of excitatory neurons
*/
struct neurons initneurons(int n, float epsilon1, float epsilon2, float dt, double g, char* adjacencyPolicy, char weightPolicy, char etaPolicy, char phasePolicy, float ratio)
{    	
	int i, j, a;
	struct neurons neurons;
	neurons.n = n;
	neurons.epsilon1 = epsilon1;
	neurons.epsilon2 = epsilon2;
	neurons.dt = dt;
	neurons.g = g;
	neurons.ratio = ratio;
	
	//Allocate and initialize adjacency matrix
	neurons.a = (int **) malloc(sizeof(int *) * n);
	
	for (i=0; i<n; i++) 
	{
		neurons.a[i] = (int *) malloc(sizeof(int) * n);
		for (j=0; j<n; j++) 
		{
			neurons.a[i][j] = 0;							//By default no connected
			
			if(strcmp(adjacencyPolicy, "f")==0)					//Fully connected, no recurrent links
			{
				if(i!=j)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "fr")==0)			//Fully connected with recurrent links (all)
			{
				neurons.a[i][j] = 1;
			}
			else if(strcmp(adjacencyPolicy, "fe")==0)			//Fully connected except one neuron with no connection
			{
				if(i!=0)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "re")==0)			//Randomly uniform connected except one neuron with no connection
			{
				if(i!=0)
				{
					neurons.a[i][j] = rand()%2;
				}
			}
			else if(strcmp(adjacencyPolicy, "ru")==0)			//Randomly uniform connected
			{
				neurons.a[i][j] = rand()%2;
			}
			else if(strcmp(adjacencyPolicy, "r23")==0)			//Randomly 2/3 connected
			{
				if((rand()%3)%2 == 0)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r13")==0)			//Randomly 1/3 connected
			{
				if((rand()%3)%2 == 1)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r14")==0)			//Randomly 1/4 connected
			{
				if((rand()%4) == 0)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r34")==0)			//Randomly 3/4 connected
			{
				if((rand()%4) != 0)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r120")==0)			//Randomly 1/20 connected
			{
				if((rand()%20) == 0)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "o")==0)			//One to one connected (in index order)
			{
				if(i == j+1)
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "t")==0)			//Two by two connected (in index order)
			{
				if((i == j+1) || (i == j-1))
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "tr")==0)			//Two by two and recurrently connected (in index order)
			{
				if((i == j+1) || (i == j-1) || (i == j))
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "c")==0)			//Two by two in circle connected (in index order)
			{
				if((i == j+1) || (i == j-1))
				{
					neurons.a[i][j] = 1;
				}
				if(i == 0)
				{
					neurons.a[i][n-1] = 1;
				}
				if(i == (n-1))
				{
					neurons.a[i][0] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "oa")==0)			//One(0) to all(1-N) connected (connected in index order)
			{
				if(i!=0)
				{
					neurons.a[i][0] = 1;
				}
				if((i == j+1) || (i == j-1))
				{
					neurons.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "oaf")==0)			//One(0) to all(1-N) connected with feedback (connected in index order)
			{
				if(i!=0)
				{
					neurons.a[i][0] = 1;
				}
				if((i == j+1) || (i == j-1))
				{
					neurons.a[i][j] = 1;
					if(j!=0)
					{
						neurons.a[0][j] = 1;
					}
				}
			}
			else if(strcmp(adjacencyPolicy, "ao")==0)			//All(1-N) to one connected (0)
			{
				if((i==0) && (j!=0))
				{
					neurons.a[i][j] = 1;
				}
				else if((j==0) && (i!=0))
				{
					neurons.a[i][j] = 1;
				}
			}
		}
	}
	
	if(strcmp(adjacencyPolicy, "rre")==0)						
	{
		createRandomRegular(neurons.a, n/4.0, n, ratio);			//Random regular connected, K=n/4 
	}
	else if(strcmp(adjacencyPolicy, "sw")==0)			
	{
		createSmallWorld(neurons.a, n/4.0, 0.1, n);					//Small world connected, K=n/4, p=0.1
	}
	else if(strcmp(adjacencyPolicy, "sf")==0)
	{
		createScaleFree(neurons.a, (0.05*n), 0.05*n, n);			//Scale free connected, K=m0, m0=0.05n
	}
	else if(strcmp(adjacencyPolicy, "swr")==0)			
	{
		createSmallWorldRandom(neurons.a, n/4.0, 0.1, n);			//Small world random connected, K=n/4, p=0.1
	}
	else if(strcmp(adjacencyPolicy, "sfr")==0)
	{
		createScaleFreeRandom(neurons.a, 0.05*n, n);	//Scale free random connected, m0=0.05n
	}
	else if(strcmp(adjacencyPolicy, "mod")==0)
	{
		createModular(neurons.a, n/4.0, n);							//Modular connected, K=n/4
	}
	else if(strcmp(adjacencyPolicy, "ml")==0)
	{
		createMultiLayer(neurons.a, 4, n);							//Multi layer connected, nb layer = 4
	}
	else if(strcmp(adjacencyPolicy, "mlf")==0)
	{
		createMultiLayerFeedback(neurons.a, 4, n);					//Multi layer connected, nb layer = 4
	}
	else if(strcmp(adjacencyPolicy, "af")==0)
	{
		createAllToFull(neurons.a, 10, n);							//All to full connected, nb input = 10
	}
	else if(strcmp(adjacencyPolicy, "res")==0)
	{
		createReservoir(neurons.a, 6, 2, n);						//Reservoirconnected, nb input = 6, nb output = 2
	}
	else if(strcmp(adjacencyPolicy, "map")==0)			
	{
		createMapTopology(neurons.a, 5, n/5);						//Map topology, n = w*h, ratio /5
	}
	
	
	//Allocate and initialize inhibitory array
	neurons.inhibitory = (bool *) malloc(sizeof(bool) * n);
	for(i=0; i<n; i++)
	{
		//The first ratio% are excitatory neurons
		if(i<(n/100.0*ratio))
		{
			neurons.inhibitory[i] = false;
		}
		else	//The % remaining are inhibitory neurons
		{
			neurons.inhibitory[i] = true; 
		}
	}	
	
	
	//Allocate and initialize coupling weights array
	neurons.k = (long double **) malloc(sizeof(long double *) * n);
	a = n - (n/100.0*ratio); //number of inhibitory neurons
	for (i=0; i<n; i++) 
	{
		neurons.k[i] = (long double *) malloc(sizeof(long double) * n);
		for (j=0; j<n; j++) 
		{
			if(neurons.a[i][j]==1)  //if the link between node i and j exists
			{
				switch(weightPolicy)
				{
					case 'g' :	//Random value distributed with a normal/gaussian law
						if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = get_random_normal(0.0, 0.001, -0.99999, -0.00001);
						}
						else
						{
							neurons.k[i][j] = get_random_normal(0.0, 0.001, 0.00001, 0.99999);
						}
						break;
					case 'r' :	//Random value uniform between -1 and 1 (depending of type of neuron)
						if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = get_random(-0.99999, -0.00001);
						}
						else
						{
							neurons.k[i][j] = get_random(0.00001, 0.99999);
						}
						//neurons.k[i][j] = get_random(-0.99999, 0.99999);
						break;
					case 'd' :	//Double cluster intialisation
						if(neurons.inhibitory[j] && neurons.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/2 && j<(n/100.0*ratio)+a/2) || (i>=(n/100.0*ratio)+a/2 && j>=(n/100.0*ratio)+a/2))
							{
								neurons.k[i][j] = -0.0001;
							}
							else
							{
								neurons.k[i][j] = -0.9999;
							}
						}
						else if(neurons.inhibitory[j])
						{
							if((i<(n/100.0*ratio)/2 && j<(n/100.0*ratio)+a/2) || (i>=(n/100.0*ratio)/2 && j>=(n/100.0*ratio)+a/2))
							{
								neurons.k[i][j] = -0.0001;
							}
							else
							{
								neurons.k[i][j] = -0.9999;
							}
						}
						else if(neurons.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/2 && j<(n/100.0*ratio)/2) || (i>=(n/100.0*ratio)+a/2 && j>=(n/100.0*ratio)/2))
							{
								neurons.k[i][j] = 0.9999;
							}
							else
							{
								neurons.k[i][j] = 0.0001;
							}
						}
						else
						{
							if((i<(n/100.0*ratio)/2 && j<(n/100.0*ratio)/2) || (i>=(n/100.0*ratio)/2 && j>=(n/100.0*ratio)/2))
							{
								neurons.k[i][j] = 0.9999;
							}
							else
							{
								neurons.k[i][j] = 0.0001;
							}
						}
						break;
					case 'q' :	//Four clusters initialisation
						if(neurons.inhibitory[j] && neurons.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/4 && j<(n/100.0*ratio)+a/4) || (i>=(n/100.0*ratio)+a/4 && j>=(n/100.0*ratio)+a/4 && i<(n/100.0*ratio)+2*a/4 && j<(n/100.0*ratio)+2*a/4) || (i>=(n/100.0*ratio)+2*a/4 && j>=(n/100.0*ratio)+2*a/4 && i<(n/100.0*ratio)+3*a/4 && j<(n/100.0*ratio)+3*a/4) || (i>=(n/100.0*ratio)+3*a/4 && j>=(n/100.0*ratio)+3*a/4))
							{
								neurons.k[i][j] = -0.0001;
							}
							else
							{
								neurons.k[i][j] = -0.9999;
							}
						}
						else if(neurons.inhibitory[j])
						{
							if((i<(n/100.0*ratio)/4 && j<(n/100.0*ratio)+a/4) || (i>=(n/100.0*ratio)/4 && j>=(n/100.0*ratio)+a/4 && i<(n/100.0*ratio)*2/4 && j<(n/100.0*ratio)+2*a/4) || (i>=(n/100.0*ratio)*2/4 && j>=(n/100.0*ratio)+2*a/4 && i<(n/100.0*ratio)*3/4 && j<(n/100.0*ratio)+3*a/4) || (i>=(n/100.0*ratio)*3/4 && j>=(n/100.0*ratio)+3*a/4))
							{
								neurons.k[i][j] = -0.0001;
							}
							else
							{
								neurons.k[i][j] = -0.9999;
							}
						}
						else if(neurons.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/4 && j<(n/100.0*ratio)/4) || (i>=(n/100.0*ratio)+a/4 && j>=(n/100.0*ratio)/4 && i<(n/100.0*ratio)+2*a/4 && j<(n/100.0*ratio)*2/4) || (i>=(n/100.0*ratio)+2*a/4 && j>=(n/100.0*ratio)*2/4 && i<(n/100.0*ratio)+3*a/4 && j<(n/100.0*ratio)*3/4) || (i>=(n/100.0*ratio)+3*a/4 && j>=(n/100.0*ratio)*3/4))
							{
								neurons.k[i][j] = 0.9999;
							}
							else
							{
								neurons.k[i][j] = 0.0001;
							}
						}
						else
						{
							if((i<(n/100.0*ratio)/4 && j<(n/100.0*ratio)/4) || (i>=(n/100.0*ratio)/4 && j>=(n/100.0*ratio)/4 && i<(n/100.0*ratio)*2/4 && j<(n/100.0*ratio)*2/4) || (i>=(n/100.0*ratio)*2/4 && j>=(n/100.0*ratio)*2/4 && i<(n/100.0*ratio)*3/4 && j<(n/100.0*ratio)*3/4) || (i>=(n/100.0*ratio)*3/4 && j>=(n/100.0*ratio)*3/4))
							{
								neurons.k[i][j] = 0.9999;
							}
							else
							{
								neurons.k[i][j] = 0.0001;
							}
						}
						break;
					case 'p' :	//Constant value of 0 or +1 (depending of type of neuron)
						if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = -0.0001;
						}
						else
						{
							neurons.k[i][j] = 0.9999;
						}
						break;
					case 'm' :	//Constant value of -1 or 0 (depending of type of neuron)
						if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = -0.9999;
						}
						else
						{
							neurons.k[i][j] = 0.0001;
						}
						break;
					case 'n' :	//Constant value null
						neurons.k[i][j] = 0.0;
						break;
					case 'i' :	//all indep
						if(neurons.inhibitory[j] && (j == i+50))
						{
							neurons.k[i][j] = -0.00001;
						}
						else if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = -0.99999;
						}
						else if(neurons.inhibitory[i] && (j==i-50))
						{
							neurons.k[i][j] = 0.99999;
						}
						else
						{
							neurons.k[i][j] = 0.00001;
						}
						break;
					case 't' :	//Constant value of -1 or +1 (depending of type of neuron)
						if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = -0.9999;
						}
						else
						{
							neurons.k[i][j] = 0.9999;
						}
						break;
					default : //Random value uniform between -1 and 1 (depending of type of neuron)
						if(neurons.inhibitory[j])
						{
							neurons.k[i][j] = get_random(-0.9999, -0.0001);
						}
						else
						{
							neurons.k[i][j] = get_random(0.0001, 0.9999);
						}
						break;
				}
			}
			else
			{
				neurons.k[i][j] = 0.0000001;
			}
		}
	}	
	
	
	//Allocate and initialize inputs array
	neurons.inputs =  (double *) calloc(n, sizeof(double));
	
	
	//Allocate and initialize bifurcation parameter array
	neurons.eta = (float *) malloc(sizeof(float) * n);
	for(i=0; i<n; i++)
	{
		switch(etaPolicy)
		{
			case 'm' :	//Identical value of -0.2
				neurons.eta[i] = -0.2; 
				break;
			case 'n' :	//Identical value of 0.0
				neurons.eta[i] = 0.0; 
				break;
			case 'p' :	//Identical value
				neurons.eta[i] = 1.5; 
				break;
			case 'r' :	//Non-identical uniformaly random value between -1.0 and 1.0
				neurons.eta[i] = get_random(-1.0, 1.0); 
				break;
			case 'g' :	//Non-identical normallly random value between 0.0 and 1.0, mean=1.5, std=0.01
				neurons.eta[i] = get_random_normal(1.5, 0.01, 1.4, 1.6);
				break;
			default :	//Identical value of 0.0
				neurons.eta[i] = 0.0; 
				break;
		}
	}
	
	
	//Allocate and initialize phases array
	neurons.theta = (long double *) malloc(sizeof(long double) * n);
	for(i=0; i<n; i++)
	{
		switch(phasePolicy)
		{
			case 'r' :	//Random value between -M_PI and M_PI
				neurons.theta[i] = get_random(-M_PI, M_PI); 
				break;
			case 'c' :	//Constant value of 0
				neurons.theta[i] = 0.0;
				break;
			case '2' :	//Two phases (0 and PI)
				neurons.theta[i] = (rand()%2)*M_PI;
				break;
			case '3' :	//Three phases (0, PI/2, PI)
				neurons.theta[i] = (rand()%3)*M_PI/2.0;
				break;
			case 'o' :	//n phases ordered
				neurons.theta[i] = ((double)i/n)*M_PI;
				break;
			default :	//Random value between -M_PI and PI
				neurons.theta[i] = get_random(-M_PI, M_PI); 
				break;
		}
	}
	
	return neurons;
}



/*
* Initialize a network of neurons with a save
*
* @param	n			number of neurons
* @param	epsilon1	the dynamic of the neurons for slow adaptation
* @param	epsilon2	the dynamic of the neurons for fast adaptation
* @param	dt			integration time step
* @param	g 		 	the coupling strength
* @param	ratio		ratio of excitatory neurons
*/
struct neurons initneuronsSaved(int n, float epsilon1, float epsilon2, float dt,  double g, float ratio)
{    	
	FILE *fp;
	int i, j;
	struct neurons neurons;
	neurons.n = n;
	neurons.epsilon1 = epsilon1;
	neurons.epsilon2 = epsilon2;
	neurons.dt = dt;
	neurons.g = g;
	neurons.ratio = ratio;
	
	
	fp = fopen("save/adjacency.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	//Allocate and initialize adjacency matrix
	neurons.a = (int **) malloc(sizeof(int *) * n);
	for (i=0; i<n; i++) 
	{
		neurons.a[i] = (int *) malloc(sizeof(int) * n);
		for (j=0; j<n; j++) 
		{
			fscanf(fp, "%d ", &neurons.a[i][j]);
		}
	}
	fclose(fp);
	
	
	fp = fopen("save/weights_matrix_3.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	//Allocate and initialize coupling weights array
	neurons.k = (long double **) malloc(sizeof(long double *) * n);
	for (i=0; i<n; i++) 
	{
		neurons.k[i] = (long double *) malloc(sizeof(long double) * n);
		for (j=0; j<n; j++) 
		{
			fscanf(fp, "%Lf ", &neurons.k[i][j]);
		}
	}
	fclose(fp);	
	
	
	//Allocate and initialize inputs array
	neurons.inputs =  (double *) calloc(n, sizeof(double));
	
	
	fp = fopen("save/eta.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	//Allocate and initialize  bifurcation parameter array
	neurons.eta = (float *) malloc(sizeof(float) * n);
	for(i=0; i<n; i++)
	{
		fscanf(fp, "%f ", &neurons.eta[i]);
	}
	fclose(fp);
	
	
	fp = fopen("save/phases.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	fseek(fp, -((18*n)+1), SEEK_END);
	//Allocate and initialize phases array
	neurons.theta = (long double *) malloc(sizeof(long double) * n);
	for(i=0; i<n; i++)
	{
		fscanf(fp, "%Lf ", &neurons.theta[i]);
	}
	fclose(fp);
	
	
	//Allocate and initialize inhibitory array
	neurons.inhibitory = (bool *) malloc(sizeof(bool) * n);
	for(i=0; i<n; i++)
	{
		//The first ratio% are excitatory neurons
		if(i<(n/100*ratio))
		{
			neurons.inhibitory[i] = false;
		}
		else	//The 20% remaining are inhibitory neurons
		{
			neurons.inhibitory[i] = true; 
		}
	}
	
	return neurons;
}


/*
* Free the memory of the neurons
*
* @param	neurons			pointer on the current network
*/
void freeNeurons(struct neurons *neurons)
{    	
	int i;
	
	//Free coupling weights array
	for (i=0; i<neurons->n; i++) 
	{
		free(neurons->k[i]);
	}
	free(neurons->k);
	
	//Free adjacency matrix
	for (i=0; i<neurons->n; i++) 
	{
		free(neurons->a[i]);
	}
	free(neurons->a);
	
	//Free input array
	free(neurons->inputs);
	
	//Free bifurcation parameter array
	free(neurons->eta);
	
	//Free phases array
	free(neurons->theta);
	
	//Free inhibitory array
	free(neurons->inhibitory);
}



/*
* Euler method 2 for the phase
*
* @param	theta_i		the phase for the neuron i
* @param	theta_j		the phases vector for the neighbor of i
* @param	dt			integration time step
* @param	g 		 	the coupling strength
* @param	a			the adjacency vector for i and its neighbors
* @param	k			the couping vector for i and its neighbors
* @param	eta			the bifurcation parameter of the neuron i
* @param	n			number of neighbors
* @param	input		input/noise added in the neuron
*/ 
long double euler_theta(long double theta_i, long double *theta_j, float dt, double g, int *a, long double *k, float eta, int n, int degree, double input)
{   
	int j;
	
	long double sum = 0.0;
	//for each neighbors of neuron i
	for(j=0; j < n; j++)
    {	
		sum += a[j]*k[j]*sin(theta_j[j]-theta_i); //-0.2*M_PI //delay
	}
  	
	//Stratonovich method
	double noise = get_random_normal(0.0, 0.1, -1.0, 1.0);
	return theta_i + dt*((1.0-cos(theta_i)) + (1.0+cos(theta_i))* (eta + g/degree * sum +input)) -0.5*pow(0.1, 2.0)*(1.0+cos(theta_i))*sin(theta_i)*dt + (1.0+cos(theta_i))*sqrt(dt)*noise;
    
    //No noise
	//return theta_i + dt*((1.0-cos(theta_i)) + (1.0+cos(theta_i))* (eta + g/degree * sum +input));
}  


/*
* Update neurons' phases of the network
*
* @param	neurons			pointer on the current network
* @param	t				iteration time
*/   
void update_phases(struct neurons *neurons, int t)
{
    int i;
	
	long double *phases = copyTable(neurons->theta, neurons->n);	//copy phases before updating
	
	//for each neurons of the network
	for(i=0; i < neurons->n; i++)
    {   
		neurons->theta[i] = euler_theta(neurons->theta[i], phases, neurons->dt, neurons->g, neurons->a[i], neurons->k[i], neurons->eta[i], neurons->n, nodeDegree(i, neurons->a, neurons->n), neurons->inputs[i]);
		
		if(neurons->theta[i] >= M_PI)	//phases must be between -M_PI and PI (periodic)
		{
			neurons->theta[i] = neurons->theta[i] - (2.0*M_PI);
			
			saveSpike(i, t*neurons->dt);
		}
		else if(neurons->theta[i] < -M_PI)
		{
			neurons->theta[i] = -M_PI;
		}
    }
	
	free(phases); 	//free memory 
}


/*
* Heaviside function
*
* @param	value	value of the function
*
* @return 1 if up to 0.1, 0.0 otherwise
*/ 
double Heaviside(double value)
{
	if(value>0.1)
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}

/*
* Runge-kutta method 4th order for the weight (k)
*
* @param	k_i_j		the couping value for link from node j and to i
* @param	theta_i		the phase for the neuron i
* @param	theta_j		the phase for the neuron j
* @param	dt			the stepsize of the RK method
* @param	epsilon		the dynamic of the neurons
* @param	plasticityPolicy	type of policy for the plasticity function
*/ 
double RK_k_Berner(long double k_i_j, long double theta_i, long double theta_j, float dt, float epsilon1, float epsilon2, double input)
{   
	double k1, k2, k3, k4;
	
	k1 = (epsilon1 + epsilon2*Heaviside(input)) * -(k_i_j + sin((theta_i - theta_j) -0.5*M_PI));
	
	k2 = (epsilon1 + epsilon2*Heaviside(input)) * -((k_i_j+k1*dt/2.0) + sin((theta_i - theta_j) -0.5*M_PI));
	
	k3 = (epsilon1 + epsilon2*Heaviside(input)) * -((k_i_j+k2*dt/2.0) + sin((theta_i - theta_j) -0.5*M_PI));
	
	k4 = (epsilon1 + epsilon2*Heaviside(input)) * -((k_i_j+k3*dt) + sin((theta_i - theta_j) -0.5*M_PI));
  
    return (k_i_j + dt*((k1 + 2.0*k2 + 2.0*k3 + k4)/6.0));
}


/*
* Runge-kutta method 4th order for the weight (k) symetric hebb
*
* @param	k_i_j			the couping value for link from node j and to i
* @param	theta_i			the phase for the neuron i
* @param	theta_j			the phase for the neuron j
* @param	dt				the stepsize of the RK method
* @param	epsilon1		the dynamic of the neurons for slow adaptation
* @param	epsilon2		the dynamic of the neurons for fast adaptation
* @param	input			input/noise added in the neuron j
* @param	inhibitory_i	if the neuron i inhibitory or not
* @param	inhibitory_j	if the neuron j inhibitory or not
*/ 
long double RK_k_symmetric(long double k_i_j, long double theta_i, long double theta_j, float dt, float epsilon1, float epsilon2, double input, bool inhibitory_i, bool inhibitory_j)
{   
	long double k1, k2, k3, k4;
	long double theta_difference = fabsl(theta_j - theta_i);
	
	if(inhibitory_j)
	{			
		if(theta_difference<M_PI)
		{
			k1 = epsilon1 * (0.0 - k_i_j) * (1.0 + k_i_j) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k2 = epsilon1 * (0.0 - (k_i_j+k1*dt/2.0)) * (1.0 + (k_i_j+k1*dt/2.0)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k3 = epsilon1 * (0.0 - (k_i_j+k2*dt/2.0)) * (1.0 + (k_i_j+k2*dt/2.0)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k4 = epsilon1 * (0.0 - (k_i_j+k3*dt)) * (1.0 + (k_i_j+k3*dt)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
		}
		else
		{	
			k1 = epsilon1 * (0.0 - k_i_j) * (1.0 + k_i_j) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k2 = epsilon1 * (0.0 - (k_i_j+k1*dt/2.0)) * (1.0 + (k_i_j+k1*dt/2.0)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k3 = epsilon1 * (0.0 - (k_i_j+k2*dt/2.0)) * (1.0 + (k_i_j+k2*dt/2.0)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k4 = epsilon1 * (0.0 - (k_i_j+k3*dt)) * (1.0 + (k_i_j+k3*dt)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));		
		}
	}
	else if(inhibitory_i)
	{	
		if(theta_difference<M_PI)
		{
			k1 = epsilon1 * (k_i_j - 0.0) * (1.0 - k_i_j) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k2 = epsilon1 * ((k_i_j+k1*dt/2.0) - 0.0) * (1.0 - (k_i_j+k1*dt/2.0)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k3 = epsilon1 * ((k_i_j+k2*dt/2.0) - 0.0) * (1.0 - (k_i_j+k2*dt/2.0)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k4 = epsilon1 * ((k_i_j+k3*dt) - 0.0) * (1.0 - (k_i_j+k3*dt)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));	
		}
		else
		{	
			k1 = epsilon1 * (k_i_j - 0.0) * (1.0 - k_i_j) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k2 = epsilon1 * ((k_i_j+k1*dt/2.0) - 0.0) * (1.0 - (k_i_j+k1*dt/2.0)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k3 = epsilon1 * ((k_i_j+k2*dt/2.0) - 0.0) * (1.0 - (k_i_j+k2*dt/2.0)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k4 = epsilon1 * ((k_i_j+k3*dt) - 0.0) * (1.0 - (k_i_j+k3*dt)) *(exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));	
		}
	}
	else
	{
		if(theta_difference<M_PI)
		{
			k1 = (epsilon1 + epsilon2*Heaviside(input)) * (k_i_j - 0.0) * (1.0 - k_i_j) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k2 = (epsilon1 + epsilon2*Heaviside(input)) * ((k_i_j+k1*dt/2.0) - 0.0) * (1.0 - (k_i_j+k1*dt/2.0)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k3 = (epsilon1 + epsilon2*Heaviside(input)) * ((k_i_j+k2*dt/2.0) - 0.0) * (1.0 - (k_i_j+k2*dt/2.0)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));
	
			k4 = (epsilon1 + epsilon2*Heaviside(input)) * ((k_i_j+k3*dt) - 0.0) * (1.0 - (k_i_j+k3*dt)) * (exp2l(-theta_difference/0.1)-exp2l((theta_difference-M_PI)/0.5));	
		}
		else
		{	
			k1 = (epsilon1 + epsilon2*Heaviside(input)) * (k_i_j - 0.0) * (1.0 - k_i_j) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k2 = (epsilon1 + epsilon2*Heaviside(input)) * ((k_i_j+k1*dt/2.0) - 0.0) * (1.0 - (k_i_j+k1*dt/2.0)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k3 = (epsilon1 + epsilon2*Heaviside(input)) * ((k_i_j+k2*dt/2.0) - 0.0) * (1.0 - (k_i_j+k2*dt/2.0)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));
	
			k4 = (epsilon1 + epsilon2*Heaviside(input)) * ((k_i_j+k3*dt) - 0.0) * (1.0 - (k_i_j+k3*dt)) * (exp2l((theta_difference-2.0*M_PI)/0.1)-exp2l(-(theta_difference-M_PI)/0.5));	
		}
	}
  	
    return (k_i_j + dt*((k1 + 2.0*k2 + 2.0*k3 + k4)/6.0));
}


/*
* Update coupling weights of the network
*
* @param	neurons			pointer on the current network
*/  
void update_weights(struct neurons *neurons) 
{
    int i, j;
	//for each neurons of the network
	for(i=0; i < neurons->n; i++)
    {
		for(j=0; j < neurons->n; j++)
	    {
			if(neurons->a[i][j]==1)  //if the link between from node j to i exists
			{
				//neurons->k[i][j] = RK_k_Berner(neurons->k[i][j], neurons->theta[i], neurons->theta[j], neurons->dt, neurons->epsilon1, neurons->epsilon2, neurons->inputs[j]);
				neurons->k[i][j] =  RK_k_symmetric(neurons->k[i][j], neurons->theta[i], neurons->theta[j], neurons->dt, neurons->epsilon1, neurons->epsilon2, neurons->inputs[j], neurons->inhibitory[i], neurons->inhibitory[j]);
			}
		}
	}
}


/*
* Reset exitatory weights with random values
*
* @param	neurons			pointer on the current network
*/ 
void reset_excitatory_weights_random(struct neurons *neurons)
{
    int i, j;
	//for each neurons of the network
	for(i=0; i < neurons->n; i++)
    {
		for(j=0; j < neurons->n; j++)
	    {
			if(!neurons->inhibitory[j] && !neurons->inhibitory[i])
			{
				neurons->k[i][j] = get_random(0.00001, 0.99999);
			}
		}
	}
}