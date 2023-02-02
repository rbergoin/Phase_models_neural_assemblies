# Phase_models_neural_assemblies
Source code of the paper: "Inhibitory neurons control the consolidation of neural assemblies via adaptation to selective stimuli"

To run a simulation execute the script "runSimulation.sh". 
To change the model, change the value of the variable "model" in this script.
To save figures, change the value of the variable "save" in this script. The figures generated are put in the "results" folder.

The parameters of each simulation (i.e. duration, number of neurons...) can be changed by adaptting the corresponding arguments in the execution of the C program.
C code don't require any particular external library.

By default the experiment excuted is the "2 clusters experiment". 
To execute another experiment, comment the lines of the current experiment and uncomment the lines between /* */ of the chosen experiment in the simulation C program.
The ".h" files containing the networks class and the utils functions should not be modified.

During each simulations raw data (i.e. spikes time, weights matrices, phases, parameters...) are generated in the corresponding ".txt" files, notably used to display the plots.

The program "plotSimulation.py" requires external libraries that can be easily installed via the "pip" command.
This python3 program reads the raw data and plot them.
No particular changes are required in this program. However you can adapt the way the neurons are sorted in the "Order" segment.
Additionally, in some plots (such as spikes, phases, or mean firing rate evolution), the x)axis is rescale for better visualization and can therefore be changed.


For any questions and additional requests, contact: raphael.bergoin@upf.edu
