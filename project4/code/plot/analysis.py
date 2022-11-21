from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from IPython import embed
import seaborn as sns 
import pandas as pd 
import os 
import plot 
import analytical 

here = os.path.abspath(".")
data_path = here + "/../../output/data/"
latex_path = here + "/../../latex/"




def load(file, folder=data_path, skiprows=1):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, delimiter=",", skiprows=skiprows)
    

def compare_analytical(sp, filename="test_analytical_1000000_cycles.txt"):
    n = load(filename)[0]
    E, E2, M, M2 = load(filename)[1:] / n
    nmax = int(filename.split("_")[2])
    eps_avg = E/4 
    m_avg = M/4 

    
    Cv_avg = (E2 - E**2) / 4 
    chi_avg = (M2 - M**2) / 4 

    eps_analytical = analytical.avg_E(1) / 4
    m_analytical   = analytical.avg_M(1) / 4 
    Cv_analytical  = analytical.CV(1) 
    chi_analytical = analytical.chi(1) 

    # embed()
    indices = np.array([10**i for i in range(int(np.log10(nmax))+1)])-1
    n_array = n[indices]
    eps_avg_array = eps_avg[indices]
    m_avg_array = m_avg[indices]
    Cv_avg_array = Cv_avg[indices]
    chi_avg_array = chi_avg[indices]

    column_names = [r"$\langle \epsilon \rangle$", r"$\langle m \rangle$", r"\langle C_V \rangle$", r"$\langle \chi \rangle$"]
    index_names = [r"$10^{%i}$"%(i) for i in (np.log10(n_array))]
    index_names.append("Analytical")

    data = np.asarray([eps_avg_array, m_avg_array, Cv_avg_array, chi_avg_array])

    data = {
        r"$N$": index_names,
        r"$\langle \epsilon \rangle$": np.append(eps_avg_array, eps_analytical),
        r"$\langle m \rangle$": np.append(m_avg_array, m_analytical),
        r"\langle C_V \rangle$": np.append(Cv_avg_array, Cv_analytical),
        r"$\langle \chi \rangle$": np.append(chi_avg_array, chi_analytical)
    }


    df = pd.DataFrame(data, index=None)
    df.style.format("{:.2f}", subset=column_names).hide(axis="index")#.to_latex(latex_path+"tables/compare_analytical.tex")
    # print(df)
    # st

    N, E_, E2_, M_, M2_ = load("analytical_1000000_cycles_T1.txt")
    # print(N)
    # exit()
    Cv_ = (E2_ - E_**2) / 4 
    chi_ = (M2_ - M_**2) / 4

    fig, ax = plt.subplots(2,2)
    ax[0,0].hlines(eps_analytical, n[0], n[-1], color='red', label='eps')
    ax[0,0].plot(n, eps_avg, '--', color='blue')
    ax[0,0].plot(N, E_/4, 'o', ms=2, color='green')

    ax[0,1].hlines(m_analytical, n[0], n[-1], color='red', label='m')
    ax[0,1].plot(n, m_avg, '--', color='blue')
    ax[0,1].plot(N, M_/4, 'o', ms=2, color='green')
    
    ax[1,0].hlines(Cv_analytical, n[0], n[-1], color='red', label='Cv')
    ax[1,0].plot(n, Cv_avg, '--', color='blue')
    ax[1,0].plot(N, Cv_, 'o', ms=2, color='green')

    ax[1,1].hlines(chi_analytical, n[0], n[-1], color='red', label='chi')
    ax[1,1].plot(n, chi_avg, '--', color='blue')
    ax[1,1].plot(N, chi_, 'o', ms=2, color='green')

    ax[0,0].legend()
    ax[0,1].legend()
    ax[1,0].legend()
    ax[1,1].legend()

    ax[0,0].set_xscale('log')
    ax[0,1].set_xscale('log')
    ax[1,0].set_xscale('log')
    ax[1,1].set_xscale('log')

    plt.show()


def equilibriation_time(sp):
    # OLD 
    n_T1, eps_T1_unord, m_T1_unord = load("equilibriate_L20_T1_unordered.txt")
    n_T1, eps_T1_order, m_T1_order =  load("equilibriate_L20_T1_ordered.txt") 

    n_T2, eps_T2_unord, m_T2_unord = load("equilibriate_L20_T2.4_unordered.txt") 
    n_T2, eps_T2_order, m_T2_order =  load("equilibriate_L20_T2.4_ordered.txt") 

    # NEW 
    N1, E1_unord, M1_unord = load("TESTequilibriate_L20_T1_unordered.txt")
    N1, E1_ord, M1_ord = load("TESTequilibriate_L20_T1_ordered.txt")
    N2, E2_unord, M2_unord = load("TESTequilibriate_L20_T2.4_unordered.txt")
    N2, E2_ord, M2_ord = load("TESTequilibriate_L20_T2.4_ordered.txt")



    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12,7))

    ax[0,0].plot(n_T1, eps_T1_order, color='red')
    ax[0,0].plot(N1, E1_ord, '--', color='blue')

    ax[0,1].plot(n_T1, eps_T1_unord, color='red')
    ax[0,1].plot(N1, E1_unord, '--', color='blue')

    ax[1,0].plot(n_T2, eps_T2_order, color='red')
    ax[1,0].plot(N2, E2_ord, '--', color='blue')

    ax[1,1].plot(n_T2, eps_T2_unord, color='red')
    ax[1,1].plot(N2, E2_unord, '--', color='blue')

    plt.show()


def pdf_histogram(sp):
    epsT1, eps2T1 = load("TESTsample_eps_L20_T1_unordered.txt")
    epsT2, eps2T2 = load("TESTsample_eps_L20_T2.4_unordered.txt")

    epsT1_, eps2T1_ = load("sample_eps_L20_T1_unordered.txt")   ### WRONG!!! 
    epsT2_, eps2T2_ = load("sample_eps_L20_T2.4_unordered.txt") ### WRONG!!! 

    # eps, eps2 = load("sample_eps_L20_T1_unorderedtest.txt")
    # print(len(epsT1))
    # exit()

    # bins = np.linspace(np.min(eps), np.max(eps), len(eps)) 

    # print(np.min(eps), np.argmax(eps))
    # eps_counts, eps_bins = np.histogram(eps)
    # print(bins)
    # print('--')
    # print(eps_bins)
    # print(eps_counts, eps_bins)
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,7))
    ax[0].hist(epsT2, bins=20)
    ax[1].hist(epsT2_, bins=20)
    # plt.plot(N, eps)
    plt.show()


sp = False # only show plots


compare_analytical(sp)
# equilibriation_time(sp)
# pdf_histogram(sp)
# N = 100001
# dn = N // 5

# nn = 10**np.arange(1, 6)
# print(nn)
# exit()
# for i in range(1,N):
    # if i%1000 == 0:
        # print(i)


