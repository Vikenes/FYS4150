from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
import seaborn as sns 
import pandas as pd 
import os 
import plot 

here = os.path.abspath(".")
data_path = here + "/../../output/data/"




def load(file, folder=data_path, skiprows=1):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, delimiter=",", skiprows=skiprows)
    



sp = False # only show plots
