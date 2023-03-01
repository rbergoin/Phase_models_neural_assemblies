// gcc -W -Wall -o simulationKuramoto simulationKuramoto.c
// ./simulationKuramoto

#include "kuramoto.h"


/*
* Calculate the rate of change of weights
*
* @param	k		coupling weights matrix at time 1
* @param	k2		coupling weights matrix at time 2
* @param	n		number of oscillators
* @param	d		number of links
*/
long double ChangeRate(long double **k, long double ** k2, int n, int d)
{
    int i, j;
    long double K = 0.0;

	//for each oscillators of the network
	for (i=0; i<n; i++) 
	{
		for (j=0; j<n; j++) 
		{
			K += fabsl(k2[i][j]-k[i][j]);
		}
	}
    
    return K/d;
}


/*
* Save a matrix in a file
*
* @param	fptr	file
* @param	m		matrix
* @param	l		number of line
* @param	c		number of column
*/
void saveMatrix(FILE *fptr, long double **m, int l, int c)
{
	int i, j;
	
	//for each element of the matrix
	for(i=0; i < l; i++)
    {
		for(j=0; j < c; j++)
	    {
			fprintf(fptr,"%3.4Lf ", m[i][j]);
		}
		fprintf(fptr, "\n");
	}
	fprintf(fptr, "\n");
}


int main(int argc, char *argv[])
{	
    srand(time(NULL)); // randomize seed
	
	int t, i, j, degree;
	FILE *fptr;
	FILE *fptr2;
	FILE *fptr3;
	FILE *fptr4;
	long double **k1;
	int nbIterations;
	int n;							//number of oscillators
    float epsilon1;					//the dynamic of the oscillators for slow adaptation
	float epsilon2;					//the dynamic of the oscillators for fast adaptation
	float dt;						//integration time step
	double g; 						//coupling strength
	float ratio;					//ratio of excitatory oscillators
	char* adjacencyPolicy;			//Type of policy for the adjacency matrix creation
	char weightPolicy;				//Type of policy for the coupling weights matrix creation
	char frequencyPolicy;			//Type of policy for the natural frequencies array creation
	char phasePolicy;				//Type of policy for the phases array creation
	int save;						//If we load a save
	float saveData = 0.0;			//Time to save data
	int nbInputOscillators;			//The number of input oscillators
	
	if(argc < 13  || argc > 13) 
	{
		nbIterations = 3000;
		n = 100;
		epsilon1 = 0.00001;
		epsilon2 = 0.1;
		dt = 0.01;
		g = 2.0;
		adjacencyPolicy = "f";
		weightPolicy = 'r';
		frequencyPolicy = 'g';
		phasePolicy = 'r';
		ratio = 80.0;
		save = 0;
	}
	else
	{
		nbIterations = atoi(argv[1]);
		n = atoi(argv[2]);
		epsilon1 = atof(argv[3]);
		epsilon2 = atof(argv[4]);
		dt = atof(argv[5]);	
		g = atof(argv[6]);
		adjacencyPolicy = argv[7];
		weightPolicy = argv[8][0];
		frequencyPolicy = argv[9][0];
		phasePolicy = argv[10][0];
		ratio = atof(argv[11]);
		save = atoi(argv[12]);
	}
	
	
	struct oscillators oscillators;
	if(save)
	{
		oscillators = initOscillatorsSaved(n, epsilon1, epsilon2, dt, g, ratio);	//Create a network of n oscillators from a save
	}
	else
	{
		oscillators = initOscillators(n, epsilon1, epsilon2, dt, g, adjacencyPolicy, weightPolicy, frequencyPolicy, phasePolicy, ratio); 		//Create a network of n oscillators
	}
	
	
	degree =  graphDegree(oscillators.a, n);
	nbInputOscillators = oscillators.n/100.0*ratio;
	
	
	/**** Save phases vector through the time ****/
	fptr = fopen("phases.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each oscillators of the network
	for(i=0; i < oscillators.n; i++)
    {
		fprintf(fptr,"%3.5Lf ",oscillators.theta[i]);
	}
	fprintf(fptr, "\n");
	
	
	/**** Save change rate of weights through the time ****/
	fptr2 = fopen("changeRates.txt","w");
	
	if(fptr2 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	
	/**** Save spikes of oscillators ****/
	fptr3 = fopen("spikes.txt","w");
	
	if(fptr3 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	fclose(fptr3);
	
	
	/**** Save weights matrix ****/
	fptr4 = fopen("weights_matrices.txt" ,"w");
	
	if(fptr4 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	saveMatrix(fptr4, oscillators.k, oscillators.n, oscillators.n);
	k1 = copyMatrix(oscillators.k, oscillators.n);	//copy weights before updating
	
	//Parameters duration (numerical values in time unit)
	nbIterations = nbIterations/dt;
	int periodInput = 20/dt;
	int durationSpontaneousActivity = 200/dt;
	int durationLearning = 1500/dt;	
	int timeConsolidation = 2000/dt;
	
	
	//Simulate for nbIterations iterations
	for(t=0; t < nbIterations; t++)
	{	
		if(((t%periodInput)==0) && (t<durationLearning) && (t>=durationSpontaneousActivity))
		{	
			//Experiment 2 clusters
			j = rand()%2;
			j = (j*nbInputOscillators/2.0)+(nbInputOscillators/2.0)/2;
			addBinaryLocalized(oscillators.inputs, n, j, nbInputOscillators/2.0);
			
			
			//Experiment 2 clusters random
			/*j = rand()%2;
			j = (j*nbInputOscillators/2.0)+(nbInputOscillators/2.0)/2;
			addBinaryLocalizedRandom(oscillators.inputs, n, j, nbInputOscillators/2.0);*/
			
			
			//Experiment 3 clusters
			/*j = rand()%3;
			j = (j*nbInputOscillators/3.0)+(nbInputOscillators/3.0)/2;
			addBinaryLocalized(oscillators.inputs, n, j, nbInputOscillators/3.0);*/
			
			
			//Experiment 3 clusters (2 stimulated)
			/*j = rand()%2;
			j = (j*nbInputOscillators/3.0)+(nbInputOscillators/3.0)/2;
			addBinaryLocalized(oscillators.inputs, n, j, nbInputOscillators/3.0);*/
			
			
			//Experiment 2 clusters with 8 overlaping oscillators
			/*if((t%(2*periodInput))==0)
			{
				addNullInputs(oscillators.inputs, nbInputOscillators);
			}
			else
			{
				if(j>50)
				{
					j = 0;
				}
				else
				{
					j = 1;
				}
				j = (j*(nbInputOscillators/2.0-4))+((nbInputOscillators+8)/2.0)/2;
				addBinaryLocalized(oscillators.inputs, n, j, (nbInputOscillators+8)/2.0);
			}*/
		}
		else if(t==durationLearning)
		{
			addNullInputs(oscillators.inputs, nbInputOscillators);
		}
		else if(t==timeConsolidation)
		{
			oscillators.epsilon1 = 1.0;
		}
		
		//Experiment recall
		/*else if(t==nbIterations)
		{
			oscillators.epsilon1 = 0.00001;
			reset_excitatory_weights_random(&oscillators);
		}
		else if(t==nbIterations+300/dt)
		{
			j = 0;
			j = (j*nbInputOscillators/2.0)+(nbInputOscillators/2.0)/2;
			addBinaryLocalized(oscillators.inputs, n, j, nbInputOscillators/2.0);
		}
		else if(t==nbIterations+305/dt)
		{
			addNullInputs(oscillators.inputs, nbInputOscillators);
		}
		else if(t==nbIterations+600/dt)
		{
			j = 1;
			j = (j*nbInputOscillators/2.0)+(nbInputOscillators/2.0)/2;
			addBinaryLocalized(oscillators.inputs, n, j, nbInputOscillators/2.0);
		}
		else if(t==nbIterations+605/dt)
		{
			addNullInputs(oscillators.inputs, nbInputOscillators);
		}*/
		
		
		update_phases(&oscillators, t);
		update_weights(&oscillators);
		
		
		/***** Save data *****/
		if(t*dt >= saveData)
		{
			saveData += 1.0;	//Save every 0.1 time unit
		
			for(i=0; i < oscillators.n; i++)
	    	{
				fprintf(fptr,"%3.5Lf ", oscillators.theta[i]);
			}
			fprintf(fptr, "\n");
		
			fprintf(fptr2,"%10.15Lf ", ChangeRate(k1, oscillators.k, oscillators.n, degree)); //calculate and register change rate of weights
			
			freeMatrix(k1, oscillators.n);
			k1 = copyMatrix(oscillators.k, oscillators.n);	//copy weights before updating
		}
		
		if(((t+1)==(durationLearning/3)) || ((t+1)==durationLearning) || ((t+1)==nbIterations))
		{
			saveMatrix(fptr4, oscillators.k, oscillators.n, oscillators.n);
		}
	}
	
	fclose(fptr);
	fclose(fptr2);
	fclose(fptr4);
	freeMatrix(k1, oscillators.n);
		
	
	/**** Save adjacency matrix ****/
	fptr = fopen("adjacency.txt" ,"w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each element of the matrix
	for(i=0; i < oscillators.n; i++)
    {
		for(j=0; j < oscillators.n; j++)
	    {
			fprintf(fptr,"%d ", oscillators.a[i][j]);
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);
	
	
	/**** Save natural frequencies ****/
	fptr = fopen("natural_frequencies.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each oscillators of the network
	for(i=0; i < oscillators.n; i++)
    {
		fprintf(fptr,"%f ", oscillators.omega[i]);
	}
	
	fclose(fptr);
	
	
	/**** Save type of oscillators ****/
	fptr = fopen("inhibitory.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each oscillators of the network
	for(i=0; i < oscillators.n; i++)
    {
		if(oscillators.inhibitory[i])
		{
			fprintf(fptr,"-1 ");
		}
		else
		{
			fprintf(fptr,"1 ");
		}
	}
	
	fclose(fptr);
	
	
	/**** Free memory ****/
	freeOscillators(&oscillators);
    
	
    return 0;
}