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
    

def compare_analytical(sp):
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


def equilibriation_time(sp):
    T1_unordered = load("equilibriate_L20_T1_unordered.txt") 
    T1_ordered = load("equilibriate_L20_T1_ordered.txt") 

    T2_unordered = load("equilibriate_L20_T2.4_unordered.txt") 
    T2_ordered = load("equilibriate_L20_T2.4_ordered.txt") 


    n_T1_unord = T1_unordered[0]
    n_T1_order = T1_ordered[0]

    n_T2_unord = T2_unordered[0]
    n_T2_order = T2_ordered[0]


    eps_T1_unord, m_T1_unord = T1_unordered[1:] / n_T1_unord
    eps_T1_order, m_T1_order = T1_ordered[1:] / n_T1_order 

    eps_T2_unord, m_T2_unord = T2_unordered[1:] / n_T2_unord
    eps_T2_order, m_T2_order = T2_ordered[1:] / n_T2_order 


    plt.plot(n_T1_order, eps_T1_order, color='red')
    plt.plot(n_T1_unord, eps_T1_unord, '--', color='red')
    plt.plot(n_T2_order, eps_T2_order, color='blue')
    plt.plot(n_T2_unord, eps_T2_unord, '--', color='blue')
    plt.show()


def pdf_histogram(sp):
    eps, eps2 = load("sample_eps_L20_T1_unordered.txt")

    # bins = np.linspace(np.min(eps), np.max(eps), len(eps)) 

    # print(np.min(eps), np.argmax(eps))
    # eps_counts, eps_bins = np.histogram(eps)
    # print(bins)
    # print('--')
    # print(eps_bins)
    # print(eps_counts, eps_bins)
    plt.hist(eps, bins=100)
    # plt.plot(N, eps)
    plt.show()


sp = False # only show plots


compare_analytical(sp)
equilibriation_time(sp)
pdf_histogram(sp)