#!/usr/bin/env python

"""
Get and plot data saved

Developped by Raphael BERGOIN

Run : python3 plotSimulation.py
"""

from math import *
import string
import numpy as np
from io import StringIO
from matplotlib.colors import ListedColormap
import codecs
import random
import warnings
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import time
import copy
import cmath
import operator
from scipy import stats
from sklearn.metrics import mutual_info_score
from datetime import datetime
import itertools


warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    

def orderParameter(phi, m, N) :
    """
        Calculate the order parameter
        
        Parameters : 
        phi -- phase vector
        m -- order
        N -- number of neurons
    """
    
    R = 0.0
    for i in range(N) :
        R += cmath.exp(1.0j*m*phi[i])
        
    return abs((1.0/N)*R)
    
    
def localFields(phi, k, N) :
    """
        Calculate the local field of each neuron at time t
        Measure the degree of synchronization of neurons among others
        
        Parameters : 
        phi -- phase vector
        k -- weights matrix
        N -- number of neurons
    """
    
    z = np.zeros(N, dtype=complex)
    for i in range(N) :
        for j in range(N) :
            z[i] += k[i][j]*cmath.exp(1.0j*phi[j])
        
    return abs((1.0/N)*z)
      

def localOrderParameter(phi, N) :
    """
        Calculate the local order parameter between each pairs of neurons
        
        Parameters : 
        phi -- phase vector
        N -- number of neurons
    """
    
    localOrderParameters = np.zeros((N,N))
    
    for i in range(N) :
        for j in range(N) :
            localOrderParameters[i][j] = cos(phi[i]-phi[j])
        
    return localOrderParameters
    
    
def autocorrelation(phi, phi2, N) :
    """
        Autocorrelation function
        
        Parameters : 
        phi -- phase vector at time 1
        phi2 -- phase vector at time 2
        N -- number of neurons
    """
    
    C = 0.0
    for i in range(N) :
        C += cmath.exp(1.0j*phi2[i])*cmath.exp(-1.0j*phi[i])
    
    return np.dot(abs((1.0/N)*C), abs((1.0/N)*C))
    


def getMeanFrequenciesAllTimes(N, T, P, spikesMatrix, timeConstant) :
    """
        Create a matrix of mean frequencies of each neurons at each time step
        
        Parameters : 
        N -- number of neurons
        T -- number of iterations of the simulation
        P -- period for the mean (in time unit)
        spikesMatrix -- matrix for each spikes of each neurons
        timeConstant -- constant factor time per iteration (precision)
    """
    meanFrequencies = np.zeros((N, T))
    
    for i in range(N) :     #For each neurons
        for t in range(int(P/timeConstant), T) :         #For each time
            meanFrequencies[i][t] = sum(map(lambda x : ((x >= t*timeConstant-P) and (x <= t*timeConstant)), spikesMatrix[i]))/P
            #meanFrequencies[i][t] = np.count_nonzero((spikesMatrix[i] >= t*timeConstant-P) & (spikesMatrix[i] <= t*timeConstant))/P
    return meanFrequencies


def getMeanFrequenciesAllTimes_2(phases, N, T, P, spikesMatrix, timeConstant) :
    """
        Create a matrix of mean frequencies of each neurons at each time step (using phases and spikes)
        
        Parameters : 
        phases -- phases of neurons through the time
        N -- number of neurons
        T -- number of iterations of the simulation
        P -- period
        spikesMatrix -- matrix for each spikes of each neurons
    """
    meanFrequencies = np.zeros((N, T))
    
    for i in range(N) :     #For each neurons
        for t in range(int(P/timeConstant), T) :         #For each time
            meanFrequencies[i][t] = (phases[i][int(t)] - phases[i][int(t-P/timeConstant)] + sum(map(lambda x : ((x >= t*timeConstant-P) and (x <= t*timeConstant)), spikesMatrix[i]))*(2.0*pi))/P
    return meanFrequencies



""" Arguments """
    
if len(sys.argv) > 1 :
    if str(sys.argv[1]) == "1" :
        save = True
    else :
        save = False
else :
    save = False
    



"""Get data saved"""   

f = open("weights_matrices.txt", "r")
matrices = []
matrix = []
for x in f:
    if len(x.strip()) == 0 :
        matrices.append(matrix)
        matrix = []
    else :
        lst = x.split()
        matrix.append([float(i) for i in lst])

f.close()
matrices = np.array(matrices)


f = open("adjacency.txt", "r")
adjacency = []
for x in f:
    lst = x.split()
    adjacency.append([int(i) for i in lst])

f.close()
adjacency = np.array(adjacency)


f = open("phases.txt", "r")
phases = []
for x in f:
    lst = x.split()
    phases.append([float(i) for i in lst])

f.close()
phases = np.array(phases)


f = open("changeRates.txt", "r")
changeRates = []
for x in f:
    lst = x.split()
    changeRates = [float(i) for i in lst]

f.close()
changeRates = np.array(changeRates)


f = open("inhibitory.txt", "r")
inhibitory = []
for x in f:
    lst = x.split()
    inhibitory = [int(i) for i in lst]

f.close()
inhibitory = np.array(inhibitory)


f = open("spikes.txt", "r")
spikes = [[] for i in range(len(phases[0]))]
for x in f:
    lst = x.split()
    spikes[int(lst[0])].append(float(lst[1]))

f.close()
spikes = np.array(spikes)



orderParameter1 = []    #list of order parameters (m=1)
orderParameter2 = []    #list of order parameters (m=2)
orderParameter3 = []    #list of order parameters (m=3)
orderParameter4 = []    #list of order parameters (m=4)
localFields1 = []        #list of local fields through the time
autocorrelations = [0.0 for x in range(len(phases)//5)]   #list of autocorrelations


adimensional = True     #If adimensional time
dt = 1.0               #Time step
durationSimulation = 2000
durationLearning = 1000
if adimensional :
    tau = 1.0
else : 
    tau = 0.01
dt = dt*tau
tfinal = (len(phases)-1.0)*dt
tpoints = np.arange(0.0, (tfinal+0.00001), dt)

ta = 0
#Calculate order parameter and autocorrelation from phase evolution
for t in range(0, len(phases)) :
    orderParameter1.append(orderParameter(phases[t], 1, len(phases[t])))    #calculate and register the order parameters (m=1)
    orderParameter2.append(orderParameter(phases[t], 2, len(phases[t])))    #calculate and register the order parameters (m=2)
    orderParameter3.append(orderParameter(phases[t], 3, len(phases[t])))    #calculate and register the order parameters (m=3)
    orderParameter4.append(orderParameter(phases[t], 4, len(phases[t])))    #calculate and register the order parameters (m=4)
    
    ta += 1
    if ta >= len(phases)/5 :
        autocorrelations.append(autocorrelation(phases[t-len(phases)//5], phases[t], len(phases[t])))  #calculate and register autocorrelations

#Moving average
orderParameter1 = np.convolve(orderParameter1, np.ones(200), 'valid') / 200
orderParameter2 = np.convolve(orderParameter2, np.ones(200), 'valid') / 200
orderParameter3 = np.convolve(orderParameter3, np.ones(200), 'valid') / 200
orderParameter4 = np.convolve(orderParameter4, np.ones(200), 'valid') / 200


#Save order parameters 1 and 2 
f = open("results/orderParameters_1.txt", "a")
for item in orderParameter1:
    f.write("%f " % item)
f.write("\n")
f.close()

f = open("results/orderParameters_2.txt", "a")
for item in orderParameter2:
    f.write("%f " % item)
f.write("\n")
f.close()


#Calculate local fields evolution
for t in range(0, len(matrices)) :
    localFields1.append(localFields(phases[int(t*((len(phases)-1)/(len(matrices)-1)))], matrices[t], len(phases[t]))) #calculate and register the local fields of neurons
localFields1 = np.array(localFields1)


#Order according to phases
nbExcit = (inhibitory == 1).sum()
excit = np.arange(nbExcit)
#cluster1 = np.arange(36)
#cluster2 = np.arange(44,nbExcit)
#excit1 = np.arange(56)
order = np.argsort(phases[-1,:])
inhib = [x for x in order if x not in excit]
#hubs = [x for x in order if x not in list(set(cluster1) | set(cluster2) | set(inhib))]
#excit2 = [x for x in order if x not in list(set(excit1) | set(inhib))]
#inhib = inhib[::-1]
if len(inhib) != 0:
    order = np.concatenate((excit, inhib))
    #order = np.concatenate((cluster1, hubs[::-1], cluster2, inhib[::-1]))
    #order = np.concatenate((excit1, excit2[::-1], inhib[::-1]))
else:
    order = excit
    
    
"""
#Order Berner
part1 = np.argsort(phases[-1,0:50])
order = np.argsort(phases[-1,:])
part2 = [x for x in order if x not in part1]
order = np.concatenate((part1, part2))
"""
    
    
#Calculate mean frequencies of the neurons
allMeanFrequencies = getMeanFrequenciesAllTimes(len(phases[0]), len(phases), 1.0, spikes[order], dt)   #mean for periods of 1.0 time unit
#allMeanFrequencies = getMeanFrequenciesAllTimes_2(np.transpose(phases)[order], len(phases[0]), len(phases), 1.0, spikes[order], dt)


"""Plot data""" 

newcmap = np.genfromtxt('./Matplotlib_colourmap_div.csv', delimiter=',')
cmap_div = ListedColormap(newcmap)
 

#Time development of the order parameters
plt.plot(tpoints[0:len(orderParameter1)], orderParameter1, label='order (m=1)') 
plt.plot(tpoints[0:len(orderParameter2)], orderParameter2, label='order (m=2)') 
#plt.plot(tpoints[0:len(orderParameter3)], orderParameter3, label='order (m=3)') 
#plt.plot(tpoints[0:len(orderParameter4)], orderParameter4, label='order (m=4)') 
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#add orders if necessary
plt.gca().legend(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
    plt.savefig('results/order_parameters.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Time development of the global order parameters', fontsize=25)
    plt.show() 
    
    
#Time development of the local fields
for i in range(len(localFields1[0])) :
    plt.plot(np.arange(0, tfinal+0.00001, dt*int((len(phases)-1)/(len(matrices)-1))), localFields1[:,i], label='Neuron %d' % i) 
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
#plt.gca().legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/local_fields.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Time development of the local fields', fontsize=25)
    plt.show()
    

#Local order parameter at the end
localOrderParameters = localOrderParameter(phases[-1,order], len(phases[-1,:]))
plt.matshow(localOrderParameters, cmap=plt.cm.seismic, vmin=0, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
if save : 
    plt.savefig('results/local_order_parameters.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Local order parameters at the end', fontsize=25)
    plt.show() 


#Adjacency matrix
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(adjacency, cmap=plt.cm.gray_r, vmin=0, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/adjacency_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Adjacency matrix', fontsize=25)
    plt.show() 

#Adjacency matrix
#kij : j presynaptic to i postsynaptic neurons (sorted)
plt.matshow(adjacency[order, :][:, order], cmap=plt.cm.gray_r, vmin=0, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/adjacency_matrix_sorted.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Adjacency matrix (sorted)', fontsize=25)
    plt.show() 


#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[0][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T0.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T0', fontsize=25)
    plt.show() 
  
#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[1][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T1.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T1', fontsize=25)
    plt.show()
    
#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[2][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T2.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T2', fontsize=25)
    plt.show()         
    
#Weights matrix after updating (sorted) 
#kij : j presynaptic to i postsynaptic neurons
plt.matshow(matrices[-1][order, :][:, order], cmap=cmap_div, vmin=-1, vmax=1)
plt.gca().set_xlabel('Presynaptic neurons j', fontsize=20)
plt.gca().set_ylabel('Postsynaptic neurons i', fontsize=20)
plt.gca().xaxis.tick_bottom()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/weights_matrix_sorted_T3.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Weight matrix (sorted) T3', fontsize=25)
    plt.show() 


#Distribution of the weights
#plt.hist(matrices[0].flatten(), bins=100, histtype=u'step', density=True, label='before updating')
#plt.hist(matrices[len(phases)//2].flatten(), bins=100, histtype=u'step', density=True, color='green', label='in the middle')
#plt.hist(matrices[-1].flatten(), bins=100, histtype=u'step', density=True, color='red', label='at the end')
plt.hist(matrices[0].flatten(), bins=100, histtype=u'step', density=True, label='before updating', log=True)
plt.hist(matrices[len(matrices)//2].flatten(), bins=100, histtype=u'step', density=True, color='green', label='in the middle', log=True)
plt.hist(matrices[-1].flatten(), bins=100, histtype=u'step', density=True, color='red', label='at the end', log=True)
plt.gca().set_xlabel('Weights', fontsize=20)
plt.gca().set_ylabel('Population', fontsize=20)
plt.gca().legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.gcf().set_size_inches(16, 9)
    plt.savefig('results/distribution_weights.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    #plt.gca().set_title('Distribution of the weights', fontsize=25)
    plt.gca().set_title('Distribution of the weights (log scale)', fontsize=25)
    plt.show() 
    

#Evolution of the change rate of weights
plt.plot(tpoints[1:len(tpoints)], changeRates) 
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/rateChange_weights.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Average evolution of the change rate of weights', fontsize=25)
    plt.show() 
    
   
#Strength of neurons
threshold = 0.0     #If necessary
a = np.absolute(matrices[-1])
#a = np.absolute(matrices[-1][order, :][:, order])
plt.plot(np.sum(np.where(a>threshold, a, 0), axis=1), '_', label='Incoming weights')
plt.plot(np.sum(np.where(a>threshold, a, 0), axis=0), '_', label='Outcoming weights', color='red')
plt.plot(np.sum(np.where(a>threshold, a, 0), axis=0)+np.sum(np.where(a>threshold, a, 0), axis=1), '_', label='Both weights', color='green')
plt.gca().set_ylim([0, len(matrices[-1])*2]) 
plt.gca().legend(fontsize=20)
plt.gca().set_ylabel('Strength', fontsize=20)
plt.gca().set_xlabel('Neurons', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/strength_neurons.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Strength of neurons', fontsize=25)
    plt.show() 
    

#Distribution strength
plt.hist(np.sum(np.where(a>threshold, a, 0), axis=1), bins=20, density=True)
#plt.hist(np.sum(np.where(a>threshold, a, 0), axis=0), bins=20, density=True)
#plt.hist(np.sum(np.where(a>threshold, a, 0), axis=0)+np.sum(np.where(a>threshold, a, 0), axis=1), bins=20, density=True)
plt.gca().set_xlabel('Strength', fontsize=20)
plt.xticks(range(0, len(matrices[-1])+1, 10))
plt.gca().set_ylabel('Population', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/distribution_strength.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Distribution of the strength', fontsize=25)
    plt.show()
   
      
#Distribution of the phase after updating
plt.hist(phases[-1,:], bins=20, density=True)
plt.gca().set_xlabel('$\\theta$', fontsize=20)
plt.xticks([-pi, 0, pi], ['$-\pi$', '0', '$\pi$'])
plt.gca().set_ylabel('Population', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/distribution_phases.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Distribution of the phase', fontsize=25)
    plt.show() 
    

#The phase in function of the number of link
plt.plot(adjacency.sum(axis=1), phases[-1,:], 'o', color='red') 
plt.gca().set_ylim([-pi, pi])
plt.gca().set_ylabel('Phase', fontsize=20)
plt.gca().set_xlabel('Degree', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/nodeDegree_Phases.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Relationship between the node degree and the phases', fontsize=25)
    plt.show() 
    
    
#Evolution of the autocorrelations
plt.plot(tpoints, autocorrelations) 
plt.gca().set_ylim([0, 1])
plt.gca().set_xlim([dt*len(phases)//5, dt*len(phases)])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
if save : 
    plt.savefig('results/autocorrelations.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Autocorrelations of phases', fontsize=25)
    plt.show() 


#Phases evolution
plt.matshow(np.transpose(phases), cmap=plt.cm.viridis, vmin=-pi, vmax=pi, extent=[0, int(len(phases)*dt), 0, len(phases[0])], origin='lower', aspect='auto')
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.gca().set_ylabel('Neurons', fontsize=20)
plt.gca().xaxis.tick_bottom()
cbar = plt.colorbar(ticks=[-pi, 0, pi])
cbar.ax.set_ylabel('$\\theta$', rotation=180)
cbar.ax.set_yticklabels(['$-\pi$', '0', '$\pi$'])  # vertically oriented colorbar
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    
    plt.gca().set_xlim(190, 200)
    plt.savefig('results/Phase_patterns_zoom_T0.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning/2, durationLearning/2+10)
    plt.savefig('results/Phase_patterns_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning, durationLearning+10)
    plt.savefig('results/Phase_patterns_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationSimulation-10, durationSimulation)
    plt.savefig('results/Phase_patterns_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.savefig('results/Phase_patterns.png', dpi=300, bbox_inches='tight')
    plt.close()
else :  
    plt.gca().set_title('Phase patterns', fontsize=25)
    plt.show() 


#Phases evolution sorted
plt.matshow(np.transpose(phases)[order], cmap=plt.cm.viridis, vmin=-pi, vmax=pi, extent=[0, int(len(phases)*dt), 0, len(phases[0])], origin='lower', aspect='auto')
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.gca().set_ylabel('Neurons', fontsize=20)
plt.gca().xaxis.tick_bottom()
cbar = plt.colorbar(ticks=[-pi, 0, pi])
cbar.ax.set_ylabel('$\\theta$', rotation=180)
cbar.ax.set_yticklabels(['$-\pi$', '0', '$\pi$'])  # vertically oriented colorbar
cbar.ax.tick_params(labelsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    
    plt.gca().set_xlim(190, 200)
    plt.savefig('results/Phase_patterns_sorted_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning/2, durationLearning/2+10)
    plt.savefig('results/Phase_patterns_sorted_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning, durationLearning+10)
    plt.savefig('results/Phase_patterns_sorted_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationSimulation-10, durationSimulation)
    plt.savefig('results/Phase_patterns_sorted_zoom_T3.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.savefig('results/Phase_patterns_sorted.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Phase patterns (sorted)', fontsize=25)
    plt.show() 


#Spikes evolution
colorsCodes = ['red' if inhibitory[i]==1 else 'blue' for i in range(len(inhibitory))]
#plt.eventplot(spikes colors='black', lineoffsets=1, linelengths=0.5, linewidths=1)
plt.eventplot(spikes, colors=colorsCodes, lineoffsets=1, linelengths=0.5, linewidths=1)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.gca().set_ylabel('Neurons', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    
    plt.gca().set_xlim(190, 200)
    plt.savefig('results/spikes_evolution_zoom_T0.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning/2, durationLearning/2+10)
    plt.savefig('results/spikes_evolution_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning, durationLearning+10)
    plt.savefig('results/spikes_evolution_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    
    plt.gca().set_xlim(durationSimulation-10, durationSimulation)
    plt.savefig('results/spikes_evolution_zoom_T3.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.savefig('results/spikes_evolution.png', dpi=300, bbox_inches='tight')
    plt.close()
else : 
    plt.gca().set_title('Spikes of neurons through the time', fontsize=25)
    plt.show()
    

   
#Spikes evolution sorted
plt.eventplot(spikes[order], colors=colorsCodes, lineoffsets=1, linelengths=0.5, linewidths=1)
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
plt.gca().set_ylabel('Neurons', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 

    plt.gca().set_xlim(190, 200)
    plt.savefig('results/spikes_evolution_sorted_zoom_T0.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning/2, durationLearning/2+10)
    plt.savefig('results/spikes_evolution_sorted_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning, durationLearning+10)
    plt.savefig('results/spikes_evolution_sorted_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationSimulation-10, durationSimulation)
    plt.savefig('results/spikes_evolution_sorted_zoom_T3.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.savefig('results/spikes_evolution_sorted.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Spikes of neurons through the time (sorted)', fontsize=25)
    plt.show()



#Firing rates


#Mean firing rate of each neuron
plt.plot(allMeanFrequencies[:, -1], '_')
plt.gca().set_xlabel('Neurons', fontsize=20)
if adimensional :
    plt.gca().set_ylabel('Mean firing rate', fontsize=20)
else:
    plt.gca().set_ylabel('Mean firing rate (Hz)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
if save : 
    plt.savefig('results/mean_firing_rate.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Mean firing rate of each neuron', fontsize=25)
    plt.show()  
   

#Evolution mean firing rates of neurons
fig, ax = plt.subplots()
for i in range(len(allMeanFrequencies)) :
    ax.plot(tpoints, allMeanFrequencies[i], label='Neuron %d' % i) 
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
    plt.gca().set_ylabel('Mean firing rate', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
    plt.gca().set_ylabel('Mean firing rate (Hz)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.gca().legend(fontsize=20)
if save : 
    plt.savefig('results/mean_firing_rates_evolution.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Mean firing rates of neurons', fontsize=25)
    plt.show() 
   

#Evolution mean firing rates of the network
meanFrequenciesNetwork = np.sum(allMeanFrequencies, axis=0)/len(allMeanFrequencies)

allMeanFrequencies1 = allMeanFrequencies[0:40, :]
allMeanFrequencies2 = allMeanFrequencies[40:80, :]
allMeanFrequencies3 = allMeanFrequencies[80:90, :]
allMeanFrequencies4 = allMeanFrequencies[90:100, :]

meanFrequenciesNetwork1 = np.sum(allMeanFrequencies1, axis=0)
meanFrequenciesNetwork1 = meanFrequenciesNetwork1/len(allMeanFrequencies1)

meanFrequenciesNetwork2 = np.sum(allMeanFrequencies2, axis=0)
meanFrequenciesNetwork2 = meanFrequenciesNetwork2/len(allMeanFrequencies2)

meanFrequenciesNetwork3 = np.sum(allMeanFrequencies3, axis=0)
meanFrequenciesNetwork3 = meanFrequenciesNetwork3/len(allMeanFrequencies3)

meanFrequenciesNetwork4 = np.sum(allMeanFrequencies4, axis=0)
meanFrequenciesNetwork4 = meanFrequenciesNetwork4/len(allMeanFrequencies4)

fig, ax = plt.subplots()
#ax.plot(tpoints, meanFrequenciesNetwork, label='Network')
ax.plot(tpoints, meanFrequenciesNetwork1, label='Cluster excitatory 1') 
ax.plot(tpoints, meanFrequenciesNetwork2, label='Cluster excitatory 2')  
ax.plot(tpoints, meanFrequenciesNetwork3, label='Cluster inhibitory 1') 
ax.plot(tpoints, meanFrequenciesNetwork4, label='Cluster inhibitory 2') 
if adimensional :
    plt.gca().set_xlabel('Time', fontsize=20)
    plt.gca().set_ylabel('Mean firing rate', fontsize=20)
else :
    plt.gca().set_xlabel('Time (s)', fontsize=20)
    plt.gca().set_ylabel('Mean firing rate (Hz)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.gca().legend(fontsize=20)
    
if save : 
    plt.gca().set_xlim(190, 200)
    plt.savefig('results/mean_firing_rate_evolution_network_T0.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning/2, durationLearning/2+10)
    plt.savefig('results/mean_firing_rate_evolution_network_zoom_T1.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationLearning, durationLearning+10)
    plt.savefig('results/mean_firing_rate_evolution_network_zoom_T2.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    plt.gca().set_xlim(durationSimulation-10, durationSimulation)
    plt.savefig('results/mean_firing_rate_evolution_network_zoom_T3.png', dpi=300, bbox_inches='tight')
    plt.gca().autoscale()
    
    
    plt.savefig('results/mean_firing_rate_evolution_network.png', dpi=300, bbox_inches='tight')
    plt.close()
else :
    plt.gca().set_title('Mean firing rates of the network', fontsize=25)
    plt.show() 
