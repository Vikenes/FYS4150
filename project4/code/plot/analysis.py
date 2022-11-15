from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
import seaborn as sns 
import pandas as pd 
import os 
import plot 
import analytical 

here = os.path.abspath(".")
data_path = here + "/../../output/data/"




def load(file, folder=data_path, skiprows=1):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, delimiter=",", skiprows=skiprows)
    

def compare_analytical():
    n = load("test_analytical_10000_cycles.txt")[0]
    E, E2, M, M2 = load("test_analytical_10000_cycles.txt")[1:] / n

    eps_avg = E/4 
    m_avg = M/4 

    
    Cv_avg = (E2 - E**2) / 4 
    chi_avg = (M2 - M**2) / 4 

    eps_analytical = analytical.avg_E(1) / 4
    m_analytical   = analytical.avg_M(1) / 4 
    Cv_analytical  = analytical.CV(1) 
    chi_analytical = analytical.chi(1) 

    fig, ax = plt.subplots(2,2)
    ax[0,0].hlines(eps_analytical, n[0], n[-1], color='red')
    ax[0,0].plot(n, eps_avg, '--', color='blue')
    ax[0,1].hlines(m_analytical, n[0], n[-1], color='red')
    ax[0,1].plot(n, m_avg, '--', color='blue')
    ax[1,0].hlines(Cv_analytical, n[0], n[-1], color='red')
    ax[1,0].plot(n, Cv_avg, '--', color='blue')
    ax[1,1].hlines(chi_analytical, n[0], n[-1], color='red')
    ax[1,1].plot(n, chi_avg, '--', color='blue')

    plt.show()


sp = False # only show plots


compare_analytical()