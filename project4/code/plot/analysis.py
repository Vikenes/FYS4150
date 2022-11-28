import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.stats import linregress
import os 

plt.style.use("seaborn")

#   Paths
here = os.path.abspath(".")
data_path = here + "/../../output/data/"
latex_path = here + "/../../latex/"
temp_path = here + "/../../output/plots/temp/"
plot_path = here +"/../../output/plots/pdfs/"

#   rc and plot params

TICKLABELSIZE = 20
LABELSIZE = 25
LEGENDSIZE = 20
TITLESIZE = 30

#   Set rc-params

plt.rc("legend", fontsize=LEGENDSIZE, fancybox=True, loc="best", frameon=True, edgecolor="black")
plt.rc("font", size=25)
plt.rc("axes", titlesize=TITLESIZE, labelsize=LABELSIZE)
plt.rc("xtick", labelsize=TICKLABELSIZE)
plt.rc("ytick", labelsize=TICKLABELSIZE)
# plt.rc("tickparams")
plt.rc("lines", linewidth=2.5)
plt.rc('text', usetex=True)
# plt.rc("label", fontsize=20)

plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = 'Times New Roman'

temp = True     # makes temporary .png files instead of .pdf

#   Global setting commands
SAVE = True
PUSH = True
SHOW = False

def save_push(fig, pdf_name, save=SAVE, push=PUSH, show=SHOW, tight=False, png_duplicate=temp):
    """
    This function handles whether you want to show,
    save and/or push the file to git.
    Args:
        fig (matplotlib.figure): Figure you want to handle
        pdfname (string): Name of output pdf file
        args (argparse)
    """
    if tight:
        fig.tight_layout()
    pdfname = pdf_name.replace('.pdf', '').strip() + ".pdf"
    file = plot_path + pdfname
    if save:
        print(f'Saving plot: {file}')
        fig.savefig(file)
        if png_duplicate:
            fig.savefig(temp_path + pdfname.replace('.pdf', '.png'))
    if push:
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")
    if show:
        plt.show()
    if not show:
        plt.close()
    else:
        plt.close()



def set_ax_info(ax, xlabel, ylabel=False, zlabel='none', style='plain', title=None, legend=True):
    """Write title and labels on an axis with the correct fontsizes.
    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel)
    if ylabel != False:
        ax.set_ylabel(ylabel)
    if zlabel != 'none':
        ax.set_zlabel(zlabel)
    ax.set_title(title)
    # ax.tick_params(axis='both', which='major', labelsize=15)
    # ax.yaxis.get_offset_text().set_fontsize(15)
    try:
        ax.ticklabel_format(style=style)
    except AttributeError:
        pass
    if legend:
        ax.legend()
    

# -----------------------------------------------------------------------------
#   General plot code above
#   XXX
#   Specific plot code below
# -----------------------------------------------------------------------------


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
    

def compare_analytical(filename="anal_Nsamples4_T1.csv"):
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
    if SAVE:
        df.style.format("{:.4f}", subset=column_names).hide(axis="index").to_latex(latex_path+"tables/compare_analytical.tex", hrules=True)
        if PUSH:
            os.system(f"git add {latex_path+'tables/compare_analytical.tex'}")
            os.system("git commit -m 'upload plot'")
            os.system("git push")
    if SHOW:
        df.style.format("{:.4f}", subset=column_names).hide(axis="index")
        print(df.to_string())


def equilibriation_time():

    e1or, m1or, n1, T1   = load("equil_L20_N100000_T1_ordered.csv")     # T=1, aligned initial spins 
    e1un, m1un, n2, T2 = load("equil_L20_N100000_T1_unordered.csv")     # T=1, unordered (random initial spins)
    e2or, m2or, n3, T4   = load("equil_L20_N100000_T2.4_ordered.csv")   # T=2.4, aligned initial spins 
    e2un, m2un, n4, T4 = load("equil_L20_N100000_T2.4_unordered.csv")   # T=2.4, unordered (random initial spins)
    
    ### ENERGY 
    fig1, ax1 = plt.subplots(nrows=1, ncols=2, figsize=(12,7), sharey=True)
    fig1.suptitle('Energy')
    ax1[0].plot(n1, e1or, '--', color='red', label='Ordered')
    ax1[0].plot(n2, e1un, ':', color='blue', label='Unordered')
    ax1[1].plot(n3, e2or, '--', color='red', label='Ordered')
    ax1[1].plot(n4, e2un, ':', color='blue', label='Unordered')
    ax1[0].set_xscale('log')
    ax1[1].set_xscale('log')
    ax1[0].set_ylim(-2.1, 0)
    ax1[1].set_ylim(-2.1, 0)
    set_ax_info(ax1[0], xlabel=r'$N$', ylabel=r'$\langle \epsilon \rangle$', title=r'$T=1$')
    set_ax_info(ax1[1], xlabel=r'$N$', title=r'$T=2.4$')
    pdfname1 = "equilibriation_time_energy"

    ### MAGNETIZATION
    fig2, ax2 = plt.subplots(nrows=1, ncols=2, figsize=(12,7), sharey=True)
    fig2.suptitle(r'Magnetization')
    
    ax2[0].plot(n1, m1or, '--', color='red', label='Ordered')
    ax2[0].plot(n2, m1un, ':', color='blue', label='Unordered')
    ax2[1].plot(n3, m2or, '--', color='red', label='Ordered')
    ax2[1].plot(n4, m2un, ':', color='blue', label='Unordered')
    ax2[0].set_xscale('log')
    ax2[1].set_xscale('log')
    ax2[0].set_ylim(0, 1.1)
    ax2[1].set_ylim(0, 1.1)
    set_ax_info(ax2[0], xlabel=r'$N$', ylabel=r'$\langle \vert m \vert \rangle$', title=r'$T=1$')
    set_ax_info(ax2[1], xlabel=r'$N$', title=r'$T=2.4$')
    pdfname2 = "equilibriation_time_magnetization"

    save_push(fig1, pdfname1)
    save_push(fig2, pdfname2)


def pdf_histogram():
   
    f1 = load("pdf_T1_unordered.csv", skiprows=0)   # (eps, abs(m), N, T)
    f2 = load("pdf_T2.4_unordered.csv", skiprows=0) # (eps, abs(m), N, T)

    e1 = f1[0]
    e2 = f2[0]
    
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

    ax[0].set_title(r'$T=1$')
    ax[1].set_title(r'$T=2.4$')
    
    textstr1 = '\n'.join((r'$\langle \epsilon \rangle = $ {:.2f}'.format(mean1), r'$\sigma^2 = $ {:.2e}'.format(var1)))
    textstr2 = '\n'.join((r'$\langle \epsilon \rangle = $ {:.2f}'.format(mean2), r'$\sigma^2 = $ {:.2e}'.format(var2)))

    ax[0].set_xlabel(r'$\epsilon$', fontsize=25)
    ax[1].set_xlabel(r'$\epsilon$', fontsize=25)
    ax[0].set_ylabel(r'$p_\epsilon(\epsilon; T)$', fontsize=25)    
    # ax[1].set_ylabel('P(E) duh')        

    ax[0].text(0.58, 0.95, textstr1, transform=ax[0].transAxes, fontsize=20, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))
    ax[1].text(0.05, 0.95, textstr2, transform=ax[1].transAxes, fontsize=20, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))

    pdfname = "histogram"
    save_push(fig, pdfname)


def PT_NT50():
    folder = data_path + "parallel/"

    ### (T, E, E2, abs(M), M2, var(E), var(M))
    f1 = load("L40_nT50_NMC100000_Neq15000_para.csv", folder)
    f2 = load("L60_nT50_NMC100000_Neq15000_para.csv", folder)
    f3 = load("L80_nT50_NMC100000_Neq15000_para.csv", folder)
    f4 = load("L100_nT50_NMC100000_Neq15000_para.csv", folder)

    files = [f1,f2,f3,f4]
    Ls = [40, 60, 80, 100]
    colors = ['blue', 'green', 'red', 'black', 'orange']

    fig, ax = plt.subplots(1,2,figsize=(12,7))
    for i, file in enumerate(files):
        T = file[0]
        varE, varM = file[-2:]
        N = (Ls[i])**2 
        Cv = varE / (N * T**2)
        chi = varM / (N * T)

        if i!=4:
            ax[0].plot(T, Cv, 'o--', ms=4, lw=1, color=colors[i], label=f'L={Ls[i]}')
            # ax[0].plot(t, cv_fit, '--',color=colors[i])
            # ax[0].plot(*[t_crit, cv_crit], 'o', color='orange')
            ax[1].plot(T, chi,'o--', ms=4, lw=1, color=colors[i], label=f'L={Ls[i]}')
        else:
            ax[0].plot(T, Cv,  lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')
            ax[1].plot(T, chi, lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')

    ax[0].set_title("Heat capacity")
    ax[1].set_title("Susceptibility")

    ax[0].set_xlabel(r'$T$')
    ax[1].set_xlabel(r'$T$')

    ax[0].set_ylabel(r'$C_V$')
    ax[1].set_ylabel(r'$\chi$')

    ax[0].legend()
    ax[1].legend()


    pdfname = "phase_transitionNT50"
    save_push(fig, pdfname)

def PT_NT100():
    folder = data_path + "parallel/"

    ### (T, E, E2, abs(M), M2, var(E), var(M))
    f1 = load("L40_nT100_NMC100000_Neq15000_para.csv", folder)
    f2 = load("L60_nT100_NMC100000_Neq15000_para.csv", folder)
    f3 = load("L80_nT100_NMC100000_Neq15000_para.csv", folder)
    f4 = load("L100_nT100_NMC100000_Neq15000_para.csv", folder)

    files = [f1,f2,f3,f4]
    Ls = [40, 60, 80, 100]
    colors = ['blue', 'green', 'red', 'black', 'orange']

    fig, ax = plt.subplots(1,2,figsize=(12,7))
    for i, file in enumerate(files):
        T = file[0]
        varE, varM = file[-2:]
        N = (Ls[i])**2 
        Cv = varE / (N * T**2)
        chi = varM / (N * T)

        if i!=4:
            ax[0].plot(T, Cv, 'o', ms=4, lw=1, color=colors[i], label=f'L={Ls[i]}')
            # ax[0].plot(t, cv_fit, '--',color=colors[i])
            # ax[0].plot(*[t_crit, cv_crit], 'o', color='orange')
            ax[1].plot(T, chi,'o', ms=4, lw=1, color=colors[i], label=f'L={Ls[i]}')
        else:
            ax[0].plot(T, Cv,  lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')
            ax[1].plot(T, chi, lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')

    ax[0].set_title("Heat capacity")
    ax[1].set_title("Susceptibility")

    ax[0].set_xlabel(r'$T$')
    ax[1].set_xlabel(r'$T$')

    ax[0].set_ylabel(r'$C_V$')
    ax[1].set_ylabel(r'$\chi$')

    ax[0].legend()
    ax[1].legend()


    pdfname = "phase_transitionNT100"
    save_push(fig, pdfname)


def PT_NT101():
    folder = data_path + "parallel/"

    ### (T, E, E2, abs(M), M2, var(E), var(M))
    f1 = load("L40_nT100_NMC1000000_Neq15000_para.csv", folder)
    f2 = load("L60_nT100_NMC1000000_Neq15000_para.csv", folder)
    f3 = load("L80_nT100_NMC1000000_Neq15000_para.csv", folder)
    f4 = load("L100_nT100_NMC1000000_Neq15000_para.csv", folder)

    files = [f1,f2,f3,f4]
    Ls = [40, 60, 80, 100]
    colors = ['blue', 'green', 'red', 'black', 'orange']

    crit_T = []
    crit_CV = []

    fig1, ax1 = plt.subplots(figsize=(12,7))
    for i, file in enumerate(files):
        T = file[0]
        varE, varM = file[-2:]
        N = (Ls[i])**2 
        Cv = varE / (N * T**2)
        chi = varM / (N * T)

        cv_spline = CubicSpline(T, Cv)
        cv_uspline = UnivariateSpline(T, Cv)
        t = np.linspace(T[0], T[-1], 1000)
        cv_fit = cv_uspline(t)
        crit_idx = np.argmax(cv_fit)
        cv_crit = cv_fit[crit_idx]
        t_crit = t[crit_idx]
        crit_T.append(t_crit)
        crit_CV.append(cv_crit)

        if i!=4:
            if i==0:
                ax1.scatter(*[t_crit, cv_crit], marker='X', s=100, c='orange', edgecolors="black", zorder=3, label=r"$T_C$")
            else:
                ax1.scatter(*[t_crit, cv_crit], marker='X', s=100, c='orange', edgecolors="black", zorder=3)
            ax1.plot(T, Cv, 'o', ms=4, lw=1, color=colors[i], zorder=1)
            ax1.plot(t, cv_fit, '--',color=colors[i], label=f'L={Ls[i]}', zorder=2)

            # ax.plot(T, chi,'.', ms=4, lw=1, color=colors[i], label=f'L={Ls[i]}')
        else:
            ax1.plot(T, Cv,  lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')
            ax1.plot(T, chi, lw=3, alpha=0.5, color=colors[i], label=f'L={Ls[i]}, comparing')

    ax1.set_title("Heat capacity")

    ax1.set_xlabel(r'$T$')

    ax1.set_ylabel(r'$C_V$')

    ax1.legend()

    pdfname1 = "phase_transitionNT101"
    save_push(fig1, pdfname1)

    #   Critical temperature analysis

    L_inv = 1/np.array(Ls)
    lreg = linregress(L_inv, crit_T)
    L_inv_array = np.linspace(-0.01, 0.03, 100)
    T_c_infty = lreg.intercept 
    a = lreg.slope 
    print(lreg)
    print(f'a: {a:.4f}')
    print(f"Tc_infty={T_c_infty:.5f}")

    linear_func = lambda xx: T_c_infty + a * xx 

    crit_temp_linreg = linear_func(L_inv)
    crit_temp_act = np.array(crit_T)
    crit_temp_linreg = np.append(crit_temp_linreg, T_c_infty)
    crit_temp_act = np.append(crit_temp_act, np.nan)

    data = {
        r'$L$': np.append(np.array(Ls), r'$\infty$'),
        r'$T_c$ from curve fit': crit_temp_act,
        r'$T_c$ from regression': crit_temp_linreg
    }

    df = pd.DataFrame(data, index=None)
    if SAVE:
        df.style.format("{:.4f}", subset=[r'$T_c$ from curve fit',r'$T_c$ from regression']).hide(axis='index').to_latex(latex_path+"tables/critical_temperatures.tex", hrules=True)
        if PUSH:
            os.system(f"git add {latex_path+'tables/critical_temperatures.tex'}")
            os.system("git commit -m 'upload plot'")
            os.system("git push")
    if SHOW:
        df.style.format("{:.4f}", subset=[r'$T_c$ from curve fit',r'$T_c$ from regression']).hide(axis='index')
        print(df.to_string())

    fig2, ax2 = plt.subplots(figsize=(12,7))
    ax2.plot(L_inv_array, linear_func(L_inv_array), color="blue", label="Linear fit")
    ax2.scatter(1/np.array(Ls), crit_T, marker='X', s=100, c='orange', edgecolors="black", zorder=3, label=r"$T_C$")
    ax2.scatter(*[0, linear_func(0)], marker='X', s=100, c='red', edgecolors="black", zorder=3, label=r"$T_C(L=\infty)$")
    ax2.hlines(T_c_infty, xmin= -0.02, xmax=0, color="red", ls='--', lw=0.7)
    ax2.vlines(0, ymin=np.min(linear_func(L_inv_array))-1, ymax=linear_func(0), colors='red', lw=0.7, ls='--')
    ax2.set_xlim(L_inv_array.min(), L_inv_array.max())
    ax2.set_ylim(linear_func(L_inv_array).min(), linear_func(L_inv_array).max())

    ax2.set_title("Critical temperatures")
    ax2.set_xlabel(r'$L^{-1}$')
    ax2.set_ylabel(r'$T$')
    ax2.legend()

    pdfname2 = 'critical_temperatures'

    save_push(fig2, pdfname2)


if __name__=="__main__":
    compare_analytical()
    equilibriation_time()
    pdf_histogram()
    PT_NT50()
    PT_NT101()
"""
Timing parameters:
T0=2
T1=2.5
N_temperature_steps=20
L=10,20,30,...,80
NMC=100 000
Neq=15 000
Parallel: OMP_NUM_THREADS=5
"""
