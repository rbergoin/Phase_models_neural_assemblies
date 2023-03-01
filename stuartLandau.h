#include "input.h"


/*
* Struct to represent a network of N Stuart-Landau oscillators 
*
*/
struct oscillators 
{
    int n;								//number of oscillators
	float ratio;						//ratio of excitatory oscillators
    float epsilon1;						//the dynamic of the oscillators for slow adaptation
	float epsilon2;						//the dynamic of the oscillators for fast adaptation
	float dt;							//integration time step
	
	int **a;							//adjacency matrix
	long double **k;					//coupling weights matrix
	double g; 							//coupling strength
	bool *inhibitory;					//if oscillators are inhibitory or excitatory
	
	long double complex *z;				//state of the oscillators
	double *omega;						//natural frequencies of the oscillators
	float *alpha;						//bifurcation parameter of the oscillators
	double *inputs;						//inputs of the oscillators
};


/*
* Copy values of a table in a new one
*
* @param	table	table to copy
* @param	n		size of the table
*/
long double complex *copyTable(long double complex *table, int n) 
{ 
	int i;
	long double complex *copy = (long double complex *) malloc(sizeof(long double complex) * n);
	
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
* Function to save a spike of an oscillator
*
* @param	i			index of the oscillator
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
* Initialize a network of oscillators
*
* @param	n					number of oscillators
* @param	epsilon1			the dynamic of the oscillators for slow adaptation
* @param	epsilon2			the dynamic of the oscillators for fast adaptation
* @param	dt					integration time step
* @param	g 		 			the coupling strength
* @param	alphaPolicy			type of policy for the bifurcation parameter
* @param	adjacencyPolicy		type of policy for the adjacency matrix creation
* @param	weightPolicy		type of policy for the coupling weights matrix creation
* @param	frequencyPolicy		type of policy for the natural frequencies array creation
* @param	statePolicy			type of policy for the states array creation
* @param	ratio				ratio of excitatory oscillators
*/
struct oscillators initOscillators(int n, float epsilon1, float epsilon2, float dt, double g, char alphaPolicy, char* adjacencyPolicy, char weightPolicy, char frequencyPolicy, char statePolicy, float ratio)
{    	
	int i, j;
	long double a;
	struct oscillators oscillators;
	oscillators.n = n;
	oscillators.epsilon1 = epsilon1;
	oscillators.epsilon2 = epsilon2;
	oscillators.dt = dt;
	oscillators.g = g;
	oscillators.ratio = ratio;
	
	//Allocate and initialize  bifurcation parameter array
	oscillators.alpha = (float *) malloc(sizeof(float) * n);
	for(i=0; i<n; i++)
	{
		switch(alphaPolicy)
		{
			case 'm' :	//Identical value of -1.0
				oscillators.alpha[i] = -1.0; 
				break;
			case 'n' :	//Identical value of 0.0
				oscillators.alpha[i] = 0.0; 
				break;
			case 'p' :	//Identical value of 1.0
				oscillators.alpha[i] = 1.0; 
				break;
			case 'r' :	//Non-identical uniformaly random value between 0.0 and 1.0
				oscillators.alpha[i] = get_random(-1.0, 1.0); 
				break;
			case 'g' :	//Non-identical normallly random value between 0.0 and 1.0, mean=0, std=0.1
				oscillators.alpha[i] = get_random_normal(0.0, 0.2, -1.0, 1.0);
				break;
			default :	//Identical value of 0.0
				oscillators.alpha[i] = 0.0; 
				break;
		}
	}
	
	
	//Allocate and initialize adjacency matrix
	oscillators.a = (int **) malloc(sizeof(int *) * n);
	for (i=0; i<n; i++) 
	{
		oscillators.a[i] = (int *) malloc(sizeof(int) * n);
		for (j=0; j<n; j++) 
		{
			oscillators.a[i][j] = 0;							//By default no connected
			
			if(strcmp(adjacencyPolicy, "f")==0)					//Fully connected, no recurrent links
			{
				if(i!=j)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "fr")==0)			//Fully connected with recurrent links (all)
			{
				oscillators.a[i][j] = 1;
			}
			else if(strcmp(adjacencyPolicy, "fe")==0)			//Fully connected except one oscillator with no connection
			{
				if(i!=0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "re")==0)			//Randomly uniform connected except one oscillator with no connection
			{
				if(i!=0)
				{
					oscillators.a[i][j] = rand()%2;
				}
			}
			else if(strcmp(adjacencyPolicy, "ru")==0)			//Randomly uniform connected
			{
				oscillators.a[i][j] = rand()%2;
			}
			else if(strcmp(adjacencyPolicy, "r23")==0)			//Randomly 2/3 connected
			{
				if((rand()%3)%2 == 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r13")==0)			//Randomly 1/3 connected
			{
				if((rand()%3)%2 == 1)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r14")==0)			//Randomly 1/4 connected
			{
				if((rand()%4) == 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r34")==0)			//Randomly 3/4 connected
			{
				if((rand()%4) != 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "r120")==0)			//Randomly 1/20 connected
			{
				if((rand()%20) == 0)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "o")==0)			//One to one connected (in index order)
			{
				if(i == j+1)
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "t")==0)			//Two by two connected (in index order)
			{
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "tr")==0)			//Two by two and recurrently connected (in index order)
			{
				if((i == j+1) || (i == j-1) || (i == j))
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "c")==0)			//Two by two in circle connected (in index order)
			{
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
				}
				if(i == 0)
				{
					oscillators.a[i][n-1] = 1;
				}
				if(i == (n-1))
				{
					oscillators.a[i][0] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "oa")==0)			//One(0) to all(1-N) connected (connected in index order)
			{
				if(i!=0)
				{
					oscillators.a[i][0] = 1;
				}
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
				}
			}
			else if(strcmp(adjacencyPolicy, "oaf")==0)			//One(0) to all(1-N) connected with feedback (connected in index order)
			{
				if(i!=0)
				{
					oscillators.a[i][0] = 1;
				}
				if((i == j+1) || (i == j-1))
				{
					oscillators.a[i][j] = 1;
					if(j!=0)
					{
						oscillators.a[0][j] = 1;
					}
				}
			}
			else if(strcmp(adjacencyPolicy, "ao")==0)			//All(1-N) to one connected (0)
			{
				if((i==0) && (j!=0))
				{
					oscillators.a[i][j] = 1;
				}
				else if((j==0) && (i!=0))
				{
					oscillators.a[i][j] = 1;
				}
			}
		}
	}
	if(strcmp(adjacencyPolicy, "rre")==0)						
	{
		createRandomRegular(oscillators.a, n/4.0, n, ratio);			//Random regular connected, K=n/4 
	}
	else if(strcmp(adjacencyPolicy, "sw")==0)			
	{
		createSmallWorld(oscillators.a, n/4.0, 0.1, n);					//Small world connected, K=n/4, p=0.1
	}
	else if(strcmp(adjacencyPolicy, "sf")==0)
	{
		createScaleFree(oscillators.a, (0.05*n), 0.05*n, n);			//Scale free connected, K=m0, m0=0.05n
	}
	else if(strcmp(adjacencyPolicy, "swr")==0)			
	{
		createSmallWorldRandom(oscillators.a, n/4.0, 0.1, n);			//Small world random connected, K=n/4, p=0.1
	}
	else if(strcmp(adjacencyPolicy, "sfr")==0)
	{
		createScaleFreeRandom(oscillators.a, 0.05*n, n);	//Scale free random connected, m0=0.05n
	}
	else if(strcmp(adjacencyPolicy, "mod")==0)
	{
		createModular(oscillators.a, n/4.0, n);							//Modular connected, K=n/4
	}
	else if(strcmp(adjacencyPolicy, "ml")==0)
	{
		createMultiLayer(oscillators.a, 4, n);							//Multi layer connected, nb layer = 4
	}
	else if(strcmp(adjacencyPolicy, "mlf")==0)
	{
		createMultiLayerFeedback(oscillators.a, 4, n);					//Multi layer connected, nb layer = 4
	}
	else if(strcmp(adjacencyPolicy, "af")==0)
	{
		createAllToFull(oscillators.a, 10, n);							//All to full connected, nb input = 10
	}
	else if(strcmp(adjacencyPolicy, "res")==0)
	{
		createReservoir(oscillators.a, 6, 2, n);						//Reservoirconnected, nb input = 6, nb output = 2
	}
	else if(strcmp(adjacencyPolicy, "res2")==0)
	{
		createReservoir2(oscillators.a, 3, 2, n);						//Reservoirconnected2, nb input = 3, nb output = 2 (for 20)
	}
	else if(strcmp(adjacencyPolicy, "map")==0)			
	{
		createMapTopology(oscillators.a, 5, n/5);						//Map topology, n = w*h, ratio /5
	}
	
	
	//Allocate and initialize inhibitory array
	oscillators.inhibitory = (bool *) malloc(sizeof(bool) * n);
	for(i=0; i<n; i++)
	{
		//The first ratio% are excitatory oscillators
		if(i<(n/100.0*ratio))
		{
			oscillators.inhibitory[i] = false;
		}
		else	//The % remaining are inhibitory oscillators
		{
			oscillators.inhibitory[i] = true; 
		}
	}	
	
	
	//Allocate and initialize coupling weights array
	oscillators.k = (long double **) malloc(sizeof(long double *) * n);
	a = n - (n/100.0*ratio); //number of inhibitory oscillators
	for (i=0; i<n; i++) 
	{
		oscillators.k[i] = (long double *) malloc(sizeof(long double) * n);
		for (j=0; j<n; j++) 
		{
			if(oscillators.a[i][j]==1)  //if the link between node i and j exists
			{
				switch(weightPolicy)
				{
					case 'g' :	//Random value distributed with a normal/gaussian law
						if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = get_random_normal(0.0, 0.001, -0.99999, -0.00001);
						}
						else
						{
							oscillators.k[i][j] = get_random_normal(0.0, 0.001, 0.00001, 0.99999);
						}
						break;
					case 'r' :	//Random value uniform between -1 and 1 (depending of type of oscillator)
						if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = get_random(-0.99999, -0.00001);
						}
						else
						{
							oscillators.k[i][j] = get_random(0.00001, 0.99999);
						}
						break;
					case 'd' :	//Double cluster intialisation
						if(oscillators.inhibitory[j] && oscillators.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/2 && j<(n/100.0*ratio)+a/2) || (i>=(n/100.0*ratio)+a/2 && j>=(n/100.0*ratio)+a/2))
							{
								oscillators.k[i][j] = -0.0001;
							}
							else
							{
								oscillators.k[i][j] = -0.9999;
							}
						}
						else if(oscillators.inhibitory[j])
						{
							if((i<(n/100.0*ratio)/2 && j<(n/100.0*ratio)+a/2) || (i>=(n/100.0*ratio)/2 && j>=(n/100.0*ratio)+a/2))
							{
								oscillators.k[i][j] = -0.0001;
							}
							else
							{
								oscillators.k[i][j] = -0.9999;
							}
						}
						else if(oscillators.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/2 && j<(n/100.0*ratio)/2) || (i>=(n/100.0*ratio)+a/2 && j>=(n/100.0*ratio)/2))
							{
								oscillators.k[i][j] = 0.9999;
							}
							else
							{
								oscillators.k[i][j] = 0.0001;
							}
						}
						else
						{
							if((i<(n/100.0*ratio)/2 && j<(n/100.0*ratio)/2) || (i>=(n/100.0*ratio)/2 && j>=(n/100.0*ratio)/2))
							{
								oscillators.k[i][j] = 0.99;
							}
							else
							{
								oscillators.k[i][j] = 0.01;
							}
						}
						break;
					case 'q' :	//Four cluster intialisation
						if(oscillators.inhibitory[j] && oscillators.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/4 && j<(n/100.0*ratio)+a/4) || (i>=(n/100.0*ratio)+a/4 && j>=(n/100.0*ratio)+a/4 && i<(n/100.0*ratio)+2*a/4 && j<(n/100.0*ratio)+2*a/4) || (i>=(n/100.0*ratio)+2*a/4 && j>=(n/100.0*ratio)+2*a/4 && i<(n/100.0*ratio)+3*a/4 && j<(n/100.0*ratio)+3*a/4) || (i>=(n/100.0*ratio)+3*a/4 && j>=(n/100.0*ratio)+3*a/4))
							{
								oscillators.k[i][j] = -0.0001;
							}
							else
							{
								oscillators.k[i][j] = -0.9999;
							}
						}
						else if(oscillators.inhibitory[j])
						{
							if((i<(n/100.0*ratio)/4 && j<(n/100.0*ratio)+a/4) || (i>=(n/100.0*ratio)/4 && j>=(n/100.0*ratio)+a/4 && i<(n/100.0*ratio)*2/4 && j<(n/100.0*ratio)+2*a/4) || (i>=(n/100.0*ratio)*2/4 && j>=(n/100.0*ratio)+2*a/4 && i<(n/100.0*ratio)*3/4 && j<(n/100.0*ratio)+3*a/4) || (i>=(n/100.0*ratio)*3/4 && j>=(n/100.0*ratio)+3*a/4))
							{
								oscillators.k[i][j] = -0.0001;
							}
							else
							{
								oscillators.k[i][j] = -0.9999;
							}
						}
						else if(oscillators.inhibitory[i])
						{
							if((i<(n/100.0*ratio)+a/4 && j<(n/100.0*ratio)/4) || (i>=(n/100.0*ratio)+a/4 && j>=(n/100.0*ratio)/4 && i<(n/100.0*ratio)+2*a/4 && j<(n/100.0*ratio)*2/4) || (i>=(n/100.0*ratio)+2*a/4 && j>=(n/100.0*ratio)*2/4 && i<(n/100.0*ratio)+3*a/4 && j<(n/100.0*ratio)*3/4) || (i>=(n/100.0*ratio)+3*a/4 && j>=(n/100.0*ratio)*3/4))
							{
								oscillators.k[i][j] = 0.9999;
							}
							else
							{
								oscillators.k[i][j] = 0.0001;
							}
						}
						else
						{
							if((i<(n/100.0*ratio)/4 && j<(n/100.0*ratio)/4) || (i>=(n/100.0*ratio)/4 && j>=(n/100.0*ratio)/4 && i<(n/100.0*ratio)*2/4 && j<(n/100.0*ratio)*2/4) || (i>=(n/100.0*ratio)*2/4 && j>=(n/100.0*ratio)*2/4 && i<(n/100.0*ratio)*3/4 && j<(n/100.0*ratio)*3/4) || (i>=(n/100.0*ratio)*3/4 && j>=(n/100.0*ratio)*3/4))
							{
								oscillators.k[i][j] = 0.9999;
							}
							else
							{
								oscillators.k[i][j] = 0.0001;
							}
						}
						break;
					case 'p' :	//Constant value of 0 or +1 (depending of type of oscillator)
						if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = -0.0001;
						}
						else
						{
							oscillators.k[i][j] = 0.9999;
						}
						break;
					case 'm' :	//Constant value of -1 or 0 (depending of type of oscillator)
						if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = -0.9999;
						}
						else
						{
							oscillators.k[i][j] = 0.0001;
						}
						break;
					case 'n' :	//Constant value null
						oscillators.k[i][j] = 0.0001;
						break;
					case 'i' :	//all indep
						if(oscillators.inhibitory[j] && (j == i+50))
						{
							oscillators.k[i][j] = -0.0001;
						}
						else if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = -0.9999;
						}
						else if(oscillators.inhibitory[i] && (j==i-50))
						{
							oscillators.k[i][j] = 0.9999;
						}
						else
						{
							oscillators.k[i][j] = 0.0001;
						}
						break;
					case 't' :	//Constant value of -1 or +1 (depending of type of oscillator)
						if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = -0.9999;
						}
						else
						{
							oscillators.k[i][j] = 0.9999;
						}
						break;
					default : //Random value uniform between -1 and 1 (depending of type of oscillator)
						if(oscillators.inhibitory[j])
						{
							oscillators.k[i][j] = get_random(-0.9999, 0.0001);
						}
						else
						{
							oscillators.k[i][j] = get_random(0.0001, 0.9999);
						}
						break;
				}
			}
			else
			{
				oscillators.k[i][j] = 0.0000001;
			}
		}
	}	
	
	
	//Allocate and initialize inputs array
	oscillators.inputs =  (double *) calloc(n, sizeof(double));
	
	
	//Allocate and initialize natural frequencies array
	oscillators.omega = (double *) malloc(sizeof(double) * n);
	for(i=0; i<n; i++)
	{
		switch(frequencyPolicy)
		{
			case 'c' :	//Constant value of 1.0
				oscillators.omega[i] = 1.0;
				break;
			case 'u' :	//Uniformaly random value between 3.0 and 3.2
				oscillators.omega[i] =  1.0 * (2.0 + 0.2*get_random(0.0, 1.0));
				break;	
			case 'v' :	//various values between inhibitory and excitatory
				if(oscillators.inhibitory[i])
				{
					oscillators.omega[i] = 1.0;
				}
				else
				{
					oscillators.omega[i] = get_random_normal(1.0, 1.0, 0.0, 2.0);
				}
				break;
			case 'e' :	//constant value of 2PI except one exception with frequency of 6PI
				if(i!=0)
				{
					oscillators.omega[i] = 1.0;
				}
				else
				{
					oscillators.omega[i] = 2.0;
				}
				break;
			case 'r' :	//Random value between 1 and 2PI
				oscillators.omega[i] = get_random(0.0, 2.0);
				break;
			case 'g' :	//Random value distributed with a normal/gaussian law
				oscillators.omega[i] = get_random_normal(sqrt(1.5)*2.0, 0.02, 1.5, 2.5);
				break;
			case '2' :	//Two frequencies (1.0 and 2.0)
				oscillators.omega[i] = (rand()%2 +1);
				break;
			case '3' :	//Three frequencies (1.0, 2.0, 3.0)
				oscillators.omega[i] = (rand()%3 +1);
				break;
			case '4' :	//Four frequencies (1.0, 2.0, 3.0, 4.0)
				oscillators.omega[i] = (rand()%4 +1);
				break;
			default : //Constant value of 1.0
				oscillators.omega[i] = 1.0;
				break;
		}
	}
	
	
	//Allocate and initialize states array
	oscillators.z = (long double complex *) malloc(sizeof(long double complex) * n);
	for(i=0; i<n; i++)
	{
		switch(statePolicy)
		{
			case 'r' :	//Uniformly random value between -1 and 1
				oscillators.z[i] = CMPLXL(get_random(-1.0, 1.0), get_random(-1.0, 1.0)); 
				break;
			case '0' :	//Synchornised initial conditions 0
				oscillators.z[i] = CMPLXL(0.0, 0.0); 
				break;
			case '1' :	//Synchornised initial conditions 1
				oscillators.z[i] = CMPLXL(1.0, 1.0); 
				break;
			case 'n' :	//Near synchrony initial conditions
				oscillators.z[i] = CMPLXL(0.001, get_random(-1.0, 1.0)); 
				break;
			case 's' :	//Uniformly random, from -sqrt(a) to +sqrt(a)
				a = get_random(0.0, 1.0);
				oscillators.z[i] = CMPLXL(get_random(-sqrt(a), sqrt(a)), get_random(-sqrt(a), sqrt(a))); 
				break;
			default :	//Random uniform value between -1 and 1
				oscillators.z[i] = CMPLXL(get_random(-1.0, 1.0), get_random(-1.0, 1.0)); 
				break;
		}
	}
	
	return oscillators;
}



/*
* Initialize a network of oscillators with a save
*
* @param	n					number of oscillators
* @param	epsilon1			the dynamic of the oscillators for slow adaptation
* @param	epsilon2			the dynamic of the oscillators for fast adaptation
* @param	dt					integration time step
* @param	g 		 			the coupling strength
* @param	ratio				ratio of excitatory oscillators
*/
struct oscillators initOscillatorsSaved(int n, float epsilon1, float epsilon2, float dt, double g, float ratio)
{    	
	FILE *fp;
	FILE *fp2;
	int i, j;
	long double phase, amplitude;
	struct oscillators oscillators;
	oscillators.n = n;
	oscillators.epsilon1 = epsilon1;
	oscillators.epsilon2 = epsilon2;
	oscillators.dt = dt;
	oscillators.g = g;
	oscillators.ratio = ratio;
	
	
	fp = fopen("save/adjacency.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	//Allocate and initialize adjacency matrix
	oscillators.a = (int **) malloc(sizeof(int *) * n);
	for (i=0; i<n; i++) 
	{
		oscillators.a[i] = (int *) malloc(sizeof(int) * n);
		for (j=0; j<n; j++) 
		{
			fscanf(fp, "%d ", &oscillators.a[i][j]);
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
	oscillators.k = (long double **) malloc(sizeof(long double *) * n);
	for (i=0; i<n; i++) 
	{
		oscillators.k[i] = (long double *) malloc(sizeof(long double) * n);
		for (j=0; j<n; j++) 
		{
			fscanf(fp, "%Lf ", &oscillators.k[i][j]);
		}
	}
	fclose(fp);	
	
	
	//Allocate and initialize inputs array
	oscillators.inputs = (double *) malloc(sizeof(double) * n);
	for(i=0; i<n; i++)
	{
		oscillators.inputs[i] = 0.0; 
	}
	
	
	fp = fopen("save/natural_frequencies.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	//Allocate and initialize natural frequencies array
	oscillators.omega = (double *) malloc(sizeof(double) * n);
	for(i=0; i<n; i++)
	{
		fscanf(fp, "%lf ", &oscillators.omega[i]);
	}
	fclose(fp);
	
	
	fp = fopen("save/alpha.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	//Allocate and initialize  bifurcation parameter array
	oscillators.alpha = (float *) malloc(sizeof(float) * n);
	for(i=0; i<n; i++)
	{
		fscanf(fp, "%f ", &oscillators.alpha[i]);
	}
	fclose(fp);
	
	
	fp = fopen("save/phases.txt", "r");
	if(fp == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	fseek(fp, -((18*n)+1), SEEK_END);
	fp2 = fopen("save/amplitudes.txt", "r");
	if(fp2 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	fseek(fp2, -((18*n)+1), SEEK_END);
	//Allocate and initialize states array
	oscillators.z = (long double complex *) malloc(sizeof(long double complex) * n);
	for(i=0; i<n; i++)
	{
		fscanf(fp, "%Lf ", &phase);
		fscanf(fp2, "%Lf ", &amplitude);
		oscillators.z[i] = CMPLXL(amplitude*cos(phase), amplitude*sin(phase)); 
	}
	fclose(fp);
	fclose(fp2);
	
	
	//Allocate and initialize inhibitory array
	oscillators.inhibitory = (bool *) malloc(sizeof(bool) * n);
	for(i=0; i<n; i++)
	{
		//The first ratio% are excitatory oscillators
		if(i<(n/100*ratio))
		{
			oscillators.inhibitory[i] = false;
		}
		else	//The 20% remaining are inhibitory oscillators
		{
			oscillators.inhibitory[i] = true; 
		}
	}
	
	return oscillators;
}


/*
* Free the memory of the oscillators
*
* @param	oscillators			pointer on the current network
*/
void freeOscillators(struct oscillators *oscillators)
{    	
	int i;
	
	//Free coupling weights array
	for (i=0; i<oscillators->n; i++) 
	{
		free(oscillators->k[i]);
	}
	free(oscillators->k);
	
	//Free adjacency matrix
	for (i=0; i<oscillators->n; i++) 
	{
		free(oscillators->a[i]);
	}
	free(oscillators->a);
	
	//Free input array
	free(oscillators->inputs);
	
	//Free natural frequencies array
	free(oscillators->omega);
	
	//Free phases array
	free(oscillators->z);
	
	//Free bifurcation parameter array
	free(oscillators->alpha);
	
	//Free inhibitory array
	free(oscillators->inhibitory);
}
    

/*
* Update oscillators' states of the network
*
* @param	oscillators			pointer on the current network
* @param	t					iteration time
*/   
void update_states(struct oscillators *oscillators, int t)
{
    int i, j;

	long double complex sum;
	long double complex *states = copyTable(oscillators->z, oscillators->n);	//copy phases before updating
	
	//for each oscillators of the network
	for(i=0; i < oscillators->n; i++)
    {   
		sum = CMPLXL(0.0, 0.0);
		
		//for each neighbors of oscillator i
		for(j=0; j < oscillators->n; j++)
	    {
			//sum += oscillators->a[i][j] * oscillators->k[i][j] * (states[j] - oscillators->z[i]);		//Conservative coupling
			sum += oscillators->a[i][j] * oscillators->k[i][j] * states[j];								//Non-conservative coupling
		}
		
		oscillators->z[i] += oscillators->dt * (oscillators->z[i]*(oscillators->alpha[i] + (oscillators->omega[i] + oscillators->inputs[i])*I - cabsl(cpowl(oscillators->z[i], 2))) + oscillators->g/nodeDegree(i, oscillators->a, oscillators->n) * sum) + sqrt(oscillators->dt)* (CMPLXL(get_random_normal(0.0, 0.2, -1.0, 1.0), get_random_normal(0.0, 0.2, -1.0, 1.0))); 
		
		if(cargl(states[i]) - cargl(oscillators->z[i]) >= M_PI)		//if there is a spike
		{				
			saveSpike(i, t*oscillators->dt);
		}
    }
	
	free(states); 	//free memory 
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
* @param	theta_i		the phase for the oscillator i
* @param	theta_j		the phase for the oscillator j
* @param	dt			the stepsize of the RK method
* @param	epsilon		the dynamic of the oscillators
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
* Runge-kutta method 4th order for the weight (k) with symmetric STDP hebbian rule
*
* @param	k_i_j			the couping value for link from node j and to i
* @param	theta_i			the phase for the oscillator i
* @param	theta_j			the phase for the oscillator j
* @param	dt				the stepsize of the RK method
* @param	epsilon1		the dynamic of the oscillators for slow adaptation
* @param	epsilon2		the dynamic of the oscillators for fast adaptation
* @param	input			input/noise added in the oscillator j
* @param	inhibitory_i	if the oscillator i inhibitory or not
* @param	inhibitory_j	if the oscillator j inhibitory or not
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
* @param	oscillators			pointer on the current network
*/  
void update_weights(struct oscillators *oscillators)
{
    int i, j;
	//for each oscillators of the network
	for(i=0; i < oscillators->n; i++)
    {
		for(j=0; j < oscillators->n; j++)
	    {
			if(oscillators->a[i][j]==1)  //if the link between from node j to i exists
			{
				//oscillators->k[i][j] = RK_k_Berner(oscillators->k[i][j], cargl(oscillators->z[i]),cargl(oscillators->z[j]), oscillators->dt, oscillators->epsilon1, oscillators->epsilon2, oscillators->inputs[j]);
				oscillators->k[i][j] = RK_k_symmetric(oscillators->k[i][j], cargl(oscillators->z[i]), cargl(oscillators->z[j]), oscillators->dt, oscillators->epsilon1, oscillators->epsilon2, oscillators->inputs[j], oscillators->inhibitory[i], oscillators->inhibitory[j]);
			}
		}
	}
}


/*
* Reset exitatory weights with random values
*
* @param	oscillators			pointer on the current network
*/ 
void reset_excitatory_weights_random(struct oscillators *oscillators)
{
    int i, j;
	//for each oscillators of the network
	for(i=0; i < oscillators->n; i++)
    {
		for(j=0; j < oscillators->n; j++)
	    {
			if(!oscillators->inhibitory[j] && !oscillators->inhibitory[i])
			{
				oscillators->k[i][j] = get_random(0.00001, 0.99999);
			}
		}
	}
}