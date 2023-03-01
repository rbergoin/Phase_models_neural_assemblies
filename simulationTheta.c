// gcc -W -Wall -o simulationTheta simulationTheta.c
// ./simulationTheta

#include "theta.h"


/*
* Calculate the rate of change of weights
*
* @param	k		coupling weights matrix at time 1
* @param	k2		coupling weights matrix at time 2
* @param	n		number of neurons
* @param	d		number of links
*/
long double ChangeRate(long double **k, long double ** k2, int n, int d)
{
    int i, j;
    long double K = 0.0;

	//for each neurons of the network
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
	int n;								//number of neurons
    float epsilon1;						//the dynamic of the neurons for slow adaptation
	float epsilon2;						//the dynamic of the neurons for fast adaptation
	float dt;							//integration time step
	double g; 							//coupling strength
	float ratio;						//ratio of excitatory neurons
	char* adjacencyPolicy;				//Type of policy for the adjacency matrix creation
	char weightPolicy;					//Type of policy for the coupling weights matrix creation
	char bifurcationParameterPolicy;	//Type of policy for the natural frequencies array creation
	char phasePolicy;					//Type of policy for the phases array creation
	int save;							//If we load a save
	float saveData = 0.0;				//Time to save data
	int nbInputNeurons;					//The number of input neurons
	
	if(argc < 13  || argc > 13) 
	{
		nbIterations = 3000;
		n = 100;
		epsilon1 = 0.00001;
		epsilon2 = 0.1;
		dt = 0.01;
		g = 1.0;
		adjacencyPolicy = "f";
		weightPolicy = 'r';
		bifurcationParameterPolicy = 'g';
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
		bifurcationParameterPolicy = argv[9][0];
		phasePolicy = argv[10][0];
		ratio = atof(argv[11]);
		save = atoi(argv[12]);
	}
	
	
	struct neurons neurons;
	if(save)
	{
		neurons =  initneuronsSaved(n, epsilon1, epsilon2, dt, g, ratio);	//Create a network of n neurons from a save
	}
	else
	{
		neurons = initneurons(n, epsilon1, epsilon2, dt, g, adjacencyPolicy, weightPolicy, bifurcationParameterPolicy, phasePolicy, ratio); 		//Create a network of n neurons
	}
	
	
	degree =  graphDegree(neurons.a, n);
	nbInputNeurons = neurons.n/100.0*ratio;
	
	
	/**** Save phases vector through the time ****/
	fptr = fopen("phases.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each neurons of the network
	for(i=0; i < neurons.n; i++)
    {
		fprintf(fptr,"%3.5Lf ",neurons.theta[i]);
	}
	fprintf(fptr, "\n");
	
	
	/**** Save change rate of weights through the time ****/
	fptr2 = fopen("changeRates.txt","w");
	
	if(fptr2 == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	
	/**** Save spikes of neurons ****/
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
	
	saveMatrix(fptr4, neurons.k, neurons.n, neurons.n);
	k1 = copyMatrix(neurons.k, neurons.n);	//copy weights before updating
	
	//Parameters duration (numerical values in time unit)
	nbIterations = nbIterations/dt;
	int periodInput = 20/dt;
	int durationSpontaneousActivity = 200/dt;
	int durationLearning = 1000/dt;	
	int timeConsolidation = 1500/dt;
	
	
	//Simulate for nbIterations iterations
	for(t=0; t < nbIterations; t++)
	{	
		if(((t%periodInput)==0) && (t<durationLearning) && (t>=durationSpontaneousActivity))
		{	
			//Experiment 2 clusters
			j = rand()%2;
			j = (j*nbInputNeurons/2.0)+(nbInputNeurons/2.0)/2;
			addBinaryLocalized(neurons.inputs, n, j, nbInputNeurons/2.0);
			
			
			//Experiment 2 clusters random
			/*j = rand()%2;
			j = (j*nbInputNeurons/2.0)+(nbInputNeurons/2.0)/2;
			addBinaryLocalizedRandom(neurons.inputs, n, j, nbInputNeurons/2.0);*/
			
			
			//Experiment 3 clusters
			/*j = rand()%3;
			j = (j*nbInputNeurons/3.0)+(nbInputNeurons/3.0)/2;
			addBinaryLocalized(neurons.inputs, n, j, nbInputNeurons/3.0);*/
			
			
			//Experiment 3 clusters (2 stimulated)
			/*j = rand()%2;
			j = (j*nbInputNeurons/3.0)+(nbInputNeurons/3.0)/2;
			addBinaryLocalized(neurons.inputs, n, j, nbInputNeurons/3.0);*/
			
			
			//Experiment 2 clusters with 8 overlaping neurons
			/*if((t%(2*periodInput))==0)
			{
				addNullInputs(neurons.inputs, nbInputNeurons);
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
				j = (j*(nbInputNeurons/2.0-4))+((nbInputNeurons+8)/2.0)/2;
				addBinaryLocalized(neurons.inputs, n, j, (nbInputNeurons+8)/2.0);
			}*/
		}
		else if(t==durationLearning)
		{
			addNullInputs(neurons.inputs, nbInputNeurons);
		}
		else if(t==timeConsolidation)
		{
			neurons.epsilon1 = 1.0;
		}
		
		//Experiment recall
		/*else if(t==nbIterations)
		{
			neurons.epsilon1 = 0.00001;
			reset_excitatory_weights_random(&neurons);
		}
		else if(t==nbIterations+300/dt)
		{
			j = 0;
			j = (j*nbInputNeurons/2.0)+(nbInputNeurons/2.0)/2;
			addBinaryLocalized(neurons.inputs, n, j, nbInputNeurons/2.0);
		}
		else if(t==nbIterations+305/dt)
		{
			addNullInputs(neurons.inputs, nbInputNeurons);
		}
		else if(t==nbIterations+600/dt)
		{
			j = 1;
			j = (j*nbInputNeurons/2.0)+(nbInputNeurons/2.0)/2;
			addBinaryLocalized(neurons.inputs, n, j, nbInputNeurons/2.0);
		}
		else if(t==nbIterations+605/dt)
		{
			addNullInputs(neurons.inputs, nbInputNeurons);
		}*/
		
		
		update_phases(&neurons, t);
		update_weights(&neurons);
		
		
		/***** Save data *****/
		if(t*dt >= saveData)
		{
			saveData += 1.0;	//Save every 1.0 time unit
		
			for(i=0; i < neurons.n; i++)
	    	{
				fprintf(fptr,"%3.5Lf ", neurons.theta[i]);
			}
			fprintf(fptr, "\n");
		
			fprintf(fptr2,"%10.15Lf ", ChangeRate(k1, neurons.k, neurons.n, degree)); //calculate and register change rate of weights
			
			freeMatrix(k1, neurons.n);
			k1 = copyMatrix(neurons.k, neurons.n);	//copy weights before updating
		}
		
		if(((t+1)==(durationLearning/2)) || ((t+1)==durationLearning) || ((t+1)==nbIterations))
		{
			saveMatrix(fptr4, neurons.k, neurons.n, neurons.n);
		}
	}
	
	fclose(fptr);
	fclose(fptr2);
	fclose(fptr4);
	freeMatrix(k1, neurons.n);
		
	
	/**** Save adjacency matrix ****/
	fptr = fopen("adjacency.txt" ,"w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each element of the matrix
	for(i=0; i < neurons.n; i++)
    {
		for(j=0; j < neurons.n; j++)
	    {
			fprintf(fptr,"%d ", neurons.a[i][j]);
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);
	
	
	/**** Save bifurcation parameters ****/
	fptr = fopen("eta.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each neurons of the network
	for(i=0; i < neurons.n; i++)
    {
		fprintf(fptr,"%f ", neurons.eta[i]);
	}
	
	fclose(fptr);
	
	
	/**** Save type of neurons ****/
	fptr = fopen("inhibitory.txt","w");
	
	if(fptr == NULL)
	{
		printf("Error!");   
		exit(1);             
	}
	
	//for each neurons of the network
	for(i=0; i < neurons.n; i++)
    {
		if(neurons.inhibitory[i])
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
	freeNeurons(&neurons);
    
	
    return 0;
}