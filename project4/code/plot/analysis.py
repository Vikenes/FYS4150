from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from IPython import embed
import seaborn as sns 
import pandas as pd 
import os 
import plot 


here = os.path.abspath(".")
data_path = here + "/../../output/data/"
latex_path = here + "/../../latex/"
fig_path = here + "/../../output/plots/temp/"


def analytical(T):
    beta = 1/T 
    Z = 12 + 4*np.cosh(8*beta)
    avg_E = -32*np.sinh(8*beta) / Z
    avg_E2 = 256/Z * np.cosh(8*beta) 
    avg_absM = 8/Z * (2 + np.exp(8*beta))
    avg_M2 = 32/Z * (1 + np.exp(8*beta))
    return avg_E, avg_E2, avg_absM, avg_M2


def load(file, folder=data_path, skiprows=0):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, delimiter=",", skiprows=skiprows)
    

def compare_analytical(sp, filename="anal_Nsamples4_T1.csv"):
    T, E, E2, absM, M2, N = load(filename)

    eps_avg = E/4 
    absm_avg = absM/4
    Cv = (E2 - E**2) / (4 * T**2)
    chi = (M2 - absM**2) / (4 * T) 

    E_avg_anal, E2_avg_anal, absM_avg_anal, M2_avg_anal = analytical(T[0])
    eps_avg_anal = E_avg_anal / 4 
    absm_avg_anal = absM_avg_anal / 4
    Cv_anal = (E2_avg_anal - E_avg_anal**2) / (4 * T[0]**2)
    chi_anal = (M2_avg_anal - absM_avg_anal**2) / (4 * T[0])

    
    

    # eps_analytical = analytical.avg_E(1) / 4
    # m_analytical   = analytical.avg_M(1) / 4 
    # Cv_analytical  = analytical.CV(1) 
    # chi_analytical = analytical.chi(1) 

    # embed()
    # indices = np.array([10**i for i in range(int(np.log10(nmax))+1)])-1
    # n_array = n[indices]
    # eps_avg_array = eps_avg[indices]
    # m_avg_array = m_avg[indices]
    # Cv_avg_array = Cv_avg[indices]
    # chi_avg_array = chi_avg[indices]

    column_names = [r"$\langle \epsilon \rangle$", r"$\langle \abs{m} \rangle$", r"$\langle C_V \rangle$", r"$\langle \chi \rangle$"]
    index_names = [r"$10^{%i}$"%(i) for i in (np.log10(N))]
    index_names.append("Analytical")

    data = np.asarray([eps_avg, absm_avg, Cv, chi])

    data = {
        r"$N$": index_names,
        r"$\langle \epsilon \rangle$": np.append(eps_avg, eps_avg_anal),
        r"$\langle \abs{m} \rangle$": np.append(absm_avg, absm_avg_anal),
        r"$\langle C_V \rangle$": np.append(Cv, Cv_anal),
        r"$\langle \chi \rangle$": np.append(chi, chi_anal)
    }


    df = pd.DataFrame(data, index=None)
    if sp:
        df.style.format("{:.3f}", subset=column_names).hide(axis="index").to_latex(latex_path+"tables/compare_analytical.tex", hrules=True)



def equilibriation_time(sp):

    e1or, m1or, n1, T1   = load("equil_L20_N100000_T1_ordered.csv")     # T=1, aligned initial spins 
    e1un, m1un, n2, T2 = load("equil_L20_N100000_T1_unordered.csv")     # T=1, unordered (random inital spins)
    e2or, m2or, n3, T4   = load("equil_L20_N100000_T2.4_ordered.csv")   # T=2.4, aligned initial spins 
    e2un, m2un, n4, T4 = load("equil_L20_N100000_T2.4_unordered.csv")   # T=2.4, unordered (random inital spins)

    
    ### ENERGY 
    fig1, ax1 = plt.subplots(nrows=1, ncols=2, figsize=(12,7))
    fig1.suptitle(r'$\langle \epsilon \rangle$')

    ax1[0].plot(n1, e1or, '--', color='red', label='Ordered')
    ax1[0].plot(n2, e1un, ':', color='blue', label='Unordered')
    ax1[1].plot(n3, e2or, '--', color='red', label='Ordered')
    ax1[1].plot(n4, e2un, ':', color='blue', label='Unordered')

    
    ax1[0].set_xscale('log')
    ax1[1].set_xscale('log')
    ax1[0].set_ylim(-2.1, 0)
    ax1[1].set_ylim(-2.1, 0)
    ax1[0].legend()
    ax1[1].legend()
    ax1[0].set_title('T=1')
    ax1[1].set_title('T=2.4')

    if sp:
        fig1.savefig(fig_path + "equilibriation_time_energy.png")
        plt.close()



    ### MAGNETIZATION
    fig2, ax2 = plt.subplots(nrows=1, ncols=2, figsize=(12,7))
    fig2.suptitle(r'$\langle \vert m \vert \rangle$')
    
    ax2[0].plot(n1, m1or, '--', color='red', label='Ordered')
    ax2[0].plot(n2, m1un, ':', color='blue', label='Unordered')
    ax2[1].plot(n3, m2or, '--', color='red', label='Ordered')
    ax2[1].plot(n4, m2un, ':', color='blue', label='Unordered')

    ax2[0].set_xscale('log')
    ax2[1].set_xscale('log')
    ax2[0].set_ylim(0, 1.1)
    ax2[1].set_ylim(0, 1.1)
    ax2[0].legend()
    ax2[1].legend()
    ax2[0].set_title('T=1')
    ax2[1].set_title('T=2.4')

    if sp:
        fig2.savefig(fig_path + "equilibriation_time_magnetization.png")
        plt.close()

    if not sp:
        plt.show()


def pdf_histogram(sp):
   
    f1 = load("pdf_T1_unordered.csv", skiprows=0)   # (eps, abs(m), N, T)
    f2 = load("pdf_T2.4_unordered.csv", skiprows=0) # (eps, abs(m), N, T)

    e1 = f1[0]
    e2 = f2[0]
    
    # print((np.max(e2) - np.min(e2))/len(e2))
    bins1 = np.arange(np.min(e1), np.max(e1)+0.01, 0.01)
    bins2 = np.arange(np.min(e2), np.max(e2)+0.01, 0.01)
    N1 = len(e1)
    N2 = len(e2)
    weights1 = np.ones(N1) / N1 
    weights2 = np.ones(N2) / N2

    mean1 = np.mean(e1)
    mean2 = np.mean(e2)

    var1 = np.var(e1)
    var2 = np.var(e2)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,7))
    
    ax[0].hist(e1, bins=bins1, weights=weights1, facecolor='red', edgecolor='black', lw=2)
    ax[1].hist(e2, bins=bins2, weights=weights2, facecolor='red', edgecolor='black', lw=0.5)

    ax[0].set_title(r'T=1, $\langle E \rangle =${:.2f}, $\sigma^2=${:.2e}'.format(mean1, var1))
    ax[1].set_title(r'T=2.4, $\langle E \rangle =${:.2f}, $\sigma^2=${:.2e}'.format(mean2, var2))
    
    ax[0].set_xlabel('E duh')
    ax[1].set_xlabel('E duh')
    ax[0].set_ylabel('P(E) duh')    
    ax[1].set_ylabel('P(E) duh')        

    if sp:
        fig.savefig(fig_path + "histogram.png")
        plt.close()

    else:
        plt.show()


def PT(sp):

    folder = data_path + "parallel/"

    ### (T, E, E2, abs(M), M2, var(E), var(M))
    f1 = load("L40_nT50_NMC100000_Neq15000_para.csv", folder)
    f2 = load("L60_nT50_NMC100000_Neq15000_para.csv", folder)
    f3 = load("L80_nT50_NMC100000_Neq15000_para.csv", folder)
    f4 = load("L100_nT50_NMC100000_Neq15000_para.csv", folder)
    FFF5 = load("L60_nT100_NMC100000_Neq15000_para.csv", folder) # previous sim with nT=100 for comparison

    files = [f1,f2,f3,f4, FFF5]
    Ls = [40, 60, 80, 100, 60]
    colors = ['blue', 'green', 'red', 'black', 'orange']
    fig, ax = plt.subplots(1,2,figsize=(12,7))
    for i, file in enumerate(files):
        T = file[0]
        varE, varM = file[-2:]
        N = (Ls[i])**2 
        Cv = varE / (N * T**2)
        chi = varM / (N * T)
        if i!=4:
            ax[0].plot(T, Cv, 'o--', ms=3, lw=0.5, color=colors[i], label=f'L={Ls[i]}')
            ax[1].plot(T, chi,'o--', ms=3, lw=0.5, color=colors[i], label=f'L={Ls[i]}')
        else:
            ax[0].plot(T, Cv,  lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')
            ax[1].plot(T, chi, lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')


    ax[0].set_xlabel('Temp')
    ax[1].set_xlabel('Temp')

    ax[0].set_ylabel('Cv')
    ax[1].set_ylabel(r'$\chi$')

    ax[0].legend()
    ax[1].legend()
    if sp:
        fig.savefig(fig_path + "phase_transition.png")
        plt.close()
        exit()
    # plt.show()

    
    #### OLD RESULTS WITH nT=100 
    cv  = lambda var, L, T: var / L**2 / T**2 
    chi = lambda var, L, T: var / L**2 / T 


    T2 = load("L20_nT50_NMC100000_Neq15000_para.csv", folder, skiprows=0)[0]
    ve2, vm2 = load("L20_nT50_NMC100000_Neq15000_para.csv", folder, skiprows=0)[-2:]
    Cv2 = cv(ve2,20,T2)
    chi2= cv(vm2,20,T2) 

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

    

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,7))
    ax[0].plot(T2, Cv2, 'o-', ms=3, lw=1, label="L20 n50")
    ax[0].plot(T3, Cv3, 'o-', ms=3, lw=0.5, label="L20")
    ax[0].plot(T4, Cv4, 'o-', ms=3, lw=0.5, label="L40")
    ax[0].plot(T5, Cv5, 'o-', ms=3, lw=0.5, label="L60")

    ax[1].plot(T2, chi2, 'o-', ms=3, lw=1, label="L20 n50")
    ax[1].plot(T3, chi3, 'o-', ms=3, lw=0.5, label="L20")
    ax[1].plot(T4, chi4, 'o-', ms=3, lw=0.5, label="L40")
    ax[1].plot(T5, chi5, 'o-', ms=3, lw=0.5, label="L60")

    ax[0].legend()
    ax[1].legend()

    plt.show()



sp = False # only show plots


compare_analytical(True)
# equilibriation_time(True)
# pdf_histogram(True)
# PT(True)

