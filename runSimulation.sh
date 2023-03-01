#!/bin/sh

#Compile and execute C program and run visualization with python program
#Developped by Raphael BERGOIN
#Run :  ./runSimulation.sh

model=2 #0 = Kuramamoto model, 1 = Stuart Landau Model, 2 = Theta neuron model


echo "Compilation..."

if [ $model == 0 ]
then
	gcc -W -Wall -o simulationKuramoto simulationKuramoto.c -lm				#Kuramamoto model
elif [ $model == 1 ]
then
	gcc -W -Wall -o simulationStuartLandau simulationStuartLandau.c -lm		#Stuart Landau Model
else
	gcc -W -Wall -o simulationTheta simulationTheta.c -lm 					#Theta neuron model
fi


it=0
number_run=1	#number of times the simulation is run

while [ $it -lt $number_run ]
do

echo "Simulation $it"

echo "Running..."


if [ $model == 0 ]
then
	#Parameters : duration simulation | number of neurons | epsilon 1 | epsilon 2 | time step | coupling strength | adjacency policy | weight policy | frequency policy | phase policy | ratio of excitatory neurons | saved network
	echo "Kuramamoto model"
	./simulationKuramoto 2000 100 0.00001 0.1 0.01 2.0 f r g r 80 0
elif [ $model == 1 ]
then
	#Parameters : duration simulation | number of neurons | epsilon 1 | epsilon 2 | time step | coupling strength | bifurcation parameter policy | adjacency policy | weight policy | frequency policy | state policy | ratio of excitatory neurons | saved network
	echo "Stuart Landau model"
	./simulationStuartLandau 2000 100 0.00001 0.1 0.01 2.0 g f r g r 80 0
else
	#Parameters : duration simulation | number of neurons | epsilon 1 | epsilon 2 | time step | coupling strength | adjacency policy | weight policy | bifurcation parameter policy | phase policy | ratio of excitatory neurons | saved network
	echo "Theta neuron model"
	./simulationTheta 2000 100 0.00001 0.1 0.01 1.0 f r g r 80 0
fi




echo "Plotting..."

save=1 #0 display figures, 1 save figures

python3 plotSimulation.py $save




it=`expr $it + 1`

done