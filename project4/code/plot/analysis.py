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

    e1un, m1un, cy1, T1 = load("equil_L20_N100000_T1_unordered.csv", skiprows=0)
    e1o, m1o, cy2, T2 = load("equil_L20_N100000_T1_ordered.csv", skiprows=0)
    e2un, m2un, cy3, T3 = load("equil_L20_N100000_T2.4_unordered.csv", skiprows=0)
    e2o, m2o, cy4, T4 = load("equil_L20_N100000_T2.4_ordered.csv", skiprows=0)

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12,7))

    ax[0,0].plot(n_T1, eps_T1_unord, color='red')
    ax[0,0].plot(cy1, e1un, ':', alpha=0.5, color='blue')
    
    ax[0,1].plot(n_T1, eps_T1_order, color='red')
    ax[0,1].plot(cy2, e1o, ':', alpha=0.5, color='black')

    ax[1,0].plot(n_T2, eps_T2_unord, color='red')
    ax[1,0].plot(cy3, e2un, ':', alpha=0.5, color='blue')
    # ax[1,0].plot(cy3, e2un/cy3, ':', alpha=0.5, color='black')

    ax[1,1].plot(n_T2, eps_T2_order, color='red')
    ax[1,1].plot(cy4, e2o, ':', alpha=0.5, color='blue')
    # ax[1,1].plot(cy4, e2o/cy3, ':', alpha=0.5, color='black')

    

    ax[0,0].set_xscale('log')
    ax[0,1].set_xscale('log')
    ax[1,0].set_xscale('log')
    ax[1,1].set_xscale('log')

    ax[0,0].set_ylim(-2.1,0)
    # ax[0,1].set_ylim(-2.1,0)
    ax[1,0].set_ylim(-2.1,0)
    ax[1,1].set_ylim(-2.1,0)


    

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


def PT(sp):
    cv  = lambda var, L, T: var / L**2 / T**2 
    chi = lambda var, L, T: var / L**2 / T 
    folder = data_path + "parallel/"

    # file1 = load("L20_nT10_NMC50000_Neq15000_serial.csv", folder, skiprows=0)
    # L1 = file1[0]
    # T1, e1, e12, m1, m12, ve1, vm1 = load("L20_nT10_NMC50000_Neq15000_serial.csv", folder, skiprows=0)
    # Cv1 = ve1 / 400 / T1**2
    # chi1 = vm1 / 400 / T1 

    # T2, e2, e22, m2, m22, ve2, vm2 = load("L20_nT10_NMC50000_Neq15000_para.csv", folder, skiprows=0)
    # Cv2 = ve2 / 400 / T2**2
    # chi2= vm2 / 400 / T2 

    T3, e3, e32, m3, m32, ve3, vm3 = load("L20_nT100_NMC100000_Neq15000_para.csv", folder, skiprows=0)
    Cv3 = ve3 / 400 / T3**2 
    chi3 = vm3 / 400 / T3 

    T4 = load("L40_nT100_NMC100000_Neq15000_para.csv", folder, skiprows=0)[0]
    ve4, vm4 = load("L40_nT100_NMC100000_Neq15000_para.csv", folder, skiprows=0)[-2:]
    Cv4 = cv(ve4,40,T4)
    chi4 = chi(vm4,40,T4)

    T5 = load("L60_nT100_NMC100000_Neq15000_para.csv", folder, skiprows=0)[0]
    ve5, vm5 = load("L60_nT100_NMC100000_Neq15000_para.csv", folder, skiprows=0)[-2:]
    Cv5 = cv(ve5,60,T5)
    chi5 = chi(vm5,60,T5)

    # print(e42/400)
    # print(T4.shape)
    # print(ve4.shape, vm4.shape)
    # exit()
    

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,7))
    # ax[0].plot(T1, Cv1, 'o-', label='parall L20')
    # ax[0].plot(T2, Cv2, 'o-', label="serial L20")
    ax[0].plot(T3, Cv3, 'o-', ms=3, lw=0.5, label="L20")
    ax[0].plot(T4, Cv4, 'o-', ms=3, lw=0.5, label="L40")
    ax[0].plot(T5, Cv5, 'o-', ms=3, lw=0.5, label="L60")

    ax[1].plot(T3, chi3, 'o-', ms=3, lw=0.5, label="L20")
    ax[1].plot(T4, chi4, 'o-', ms=3, lw=0.5, label="L40")
    ax[1].plot(T5, chi5, 'o-', ms=3, lw=0.5, label="L60")


    # ax[1].plot(T2, m2**2, 'o-', label='22')
    # ax[1].plot(T2, m22, 'o-', label='22')
    # ax[1].plot(T3, m3**2, 'o-', label='3')
    # ax[1].plot(T3, m32, 'o-', label='32')


    ax[0].legend()
    ax[1].legend()
    plt.show()



sp = False # only show plots


# compare_analytical(sp)
equilibriation_time(sp)
# pdf_histogram(sp)
# PT(sp)

