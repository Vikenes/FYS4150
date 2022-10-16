import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import seaborn as sns
import pandas as pd
from datetime import datetime

"""
Script used for plotting. 

To show plots without saving: Execute the function "show_plots()" in the end 
"""


# The style we want
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('font', family='DejaVu Sans')

# other rc parameters
plt.rc('figure', figsize=(12,7))

# folder paths
here = os.path.abspath(".")
data_path = here + "/../../output/data/"
plot_path = here + "/../../output/plots/"
latex_path = here + "/../../latex/"


def save_push(fig, pdf_name, save=True, push=False, show=False, tight=True):
    """
    This function handles wether you want to show,
    save and/or push the file to git.

    Args:
        fig (matplotlib.figure): Figure you want to handle
        pdfname (string): Name of output pdf file
        args (argparse)
    """
    if tight:
        fig.tight_layout()
    file = plot_path + pdf_name.replace('.pdf', '').strip() + ".pdf"
    if save:
        print(f'Saving plot: {file}')
        fig.savefig(file)
    if push:
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")
    if show:
        plt.show()
    if not show:
        plt.close()
    else:
        plt.clf()

def set_ax_info(ax, xlabel, ylabel, style='plain', title=None, legend=True):
    """Write title and labels on an axis with the correct fontsizes.

    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.yaxis.get_offset_text().set_fontsize(15)
    try:
        ax.ticklabel_format(style=style)
    except AttributeError:
        pass
    if legend:
        ax.legend(fontsize=15)


def test_single_particle():
    # FIXME

    fig0, ax0 = plt.subplots(layout='constrained')
    fig1, ax1 = plt.subplots(layout='constrained')
    fig2, ax2 = plt.subplots(layout='constrained')

    colors = ['dodgerblue', 'olive', 'darkorange', 'navy']
    markersizes = [16, 13, 9, 4]
    for k in range(1, 5):
        FE = np.loadtxt(data_path + f"tests/FE/single_n{k}.txt", unpack=True, delimiter=",", skiprows=1)
        RK4 = np.loadtxt(data_path + f"tests/RK4/single_n{k}.txt", unpack=True, delimiter=",", skiprows=1)

        t = FE[0] # times
        n = len(t)-1
        zE = FE[3]
        zR = RK4[3]

        omega_z = np.sqrt(2*1/40.078 * 9.65)
        z_anal = 20*np.cos(omega_z * t)


        abserrR = np.abs(z_anal - zR)
        
        errE = np.abs(z_anal - zE)/np.abs(z_anal)
        errR = abserrR/np.abs(z_anal)



        ax0.plot(t, zE, c=colors[k-1], lw=2, ls='--')
        ax0.plot(t, zR, c=colors[k-1], lw=2, ls=':')
        if k == 4:
            ax0.plot(t, z_anal, lw=10, alpha=0.3, color='green', label='analytical')
            
            zero_points = []
            K = int(n/11)
            for j in range(11):
                kj = int(j*K)
                t_vals = t[kj:kj+K+1]
                z_vals = z_anal[kj:kj+K+1]
                idx = np.argmin(np.abs(z_vals))
                zero_points.append(t_vals[idx])

        ax1.plot(t, errE, 'o-', c=colors[k-1], lw=1.5, ms=markersizes[k-1], label=r'$n=%i$'%(n), alpha=.7)
        ax2.plot(t, errR, 'o-', c=colors[k-1], lw=1.5, ms=markersizes[k-1], label=r'$n=%i$'%(n), alpha=.7)
    
    
    ax0.plot(-1, 0, lw=2, ls='--', c='grey', alpha=0.6, label='Forward-Euler')
    ax0.plot(-1, 0, lw=2, ls=':', c='grey', alpha=0.6, label='Runge-Kutta 4')
    ax0.set_xlim(t[0], t[-1])

    for ax in [ax1, ax2]:
        for t0 in zero_points:
            ax.axvline(t0, ls=':', lw=.7, color='r', alpha=.6)
        ax.set_yscale('log')


    set_ax_info(ax0, r"$t$ [$\mu$s]", r"$z$ [$\mu$m]")
    set_ax_info(ax1, r"$t$ [$\mu$s]", r"$\mathbf{r}$", title="Forward-Euler")
    set_ax_info(ax2, r"$t$ [$\mu$s]", r"$\mathbf{r}$", title="Runge-Kutta 4")

    plt.show()



def test_double_particle():

    # FIXME
    FE = np.loadtxt(data_path + "tests/FE/double_with.txt", unpack=True, delimiter=",", skiprows=1)
    RK4 = np.loadtxt(data_path + "tests/RK4/double_with.txt", unpack=True, delimiter=",", skiprows=1)

    Nt = int(len(FE[0])/2)
    t = FE[0,:Nt]
    rE = np.zeros((Nt,3,2))
    rE[:,:,0] = FE[1:4, :Nt].T
    rE[:,:,1] = FE[1:4, Nt:].T
    rR = np.zeros((Nt,3,2))
    rR[:,:,0] = RK4[1:4, :Nt].T
    rR[:,:,1] = RK4[1:4, Nt:].T

    vE = np.zeros((Nt,3,2))
    vE[:,:,0] = FE[4:7, :Nt].T
    vE[:,:,1] = FE[4:7, Nt:].T
    vR = np.zeros((Nt,3,2))
    vR[:,:,0] = RK4[4:7, :Nt].T
    vR[:,:,1] = RK4[4:7, Nt:].T

    cmap = ['copper', 'bone']
    c = ['navy', 'darkorange']
    cp = ['yellow', 'red'] 
    fig, axes = plt.subplots(ncols=2, layout='constrained')
    ax1, ax2 = axes.flat
    p = 0
    coord = 0
    ax1.plot(rR[:,coord,p], rR[:,coord,p], lw=.7, c=c[p], ls='-',  alpha=.5)
    ax1.scatter(rR[:,coord,p], rR[:,coord,p], s=3, marker='o', c=t, cmap=cmap[p], alpha=.7)


    set_ax_info(ax1, r'$x$ [$\mu$m]', r'$y$ [$\mu$m]', title="Without interactions", legend=False)
    set_ax_info(ax2, r'$y$ [$\mu$m]', r'$y$ [$\mu$m]', title="With interactions",    legend=False)


    plt.show()





    fig, axes = plt.subplots(ncols=2, layout='constrained')
    ax1, ax2 = axes.flat

    

    for j, coord in enumerate([0, 2]):
        ax = axes.flat[j]
        for p in [0,1]:
            #ax.plot(rE[:,coord,p], vE[:,coord,p], lw=0.9, c=c[p], ls='--', alpha=.5)
            #ax.scatter(rE[:,coord,p], vE[:,coord,p], s=3, marker='¨', c=t, cmap=cmap[p], alpha=.7)
            ax.plot(rR[:,coord,p], vR[:,coord,p], lw=.7, c=c[p], ls='-',  alpha=.5)
            ax.scatter(rR[:,coord,p], vR[:,coord,p], s=3, marker='o', c=t, cmap=cmap[p], alpha=.7)

        for p in [0,1]:
            ax.plot(rR[0,coord,p],  vR[0,coord,p],  marker="P", ms=12, c=cp[p], alpha=1)
            ax.plot(rR[-1,coord,p], vR[-1,coord,p], marker="*", ms=12, c=cp[p], alpha=1)
            # show start and stop pos!

    #units ={"position":" [μm]", "velocity":" [μm/s]"}
    set_ax_info(ax1, r'$x$ [$\mu$m]', r'$v_x$ [$\mu$m/s]', legend=False)
    set_ax_info(ax2, r'$z$ [$\mu$m]', r'$v_z$ [$\mu$m/s]', legend=False)


    plt.show()


if __name__=="__main__":
    test_single_particle()

    test_double_particle()
