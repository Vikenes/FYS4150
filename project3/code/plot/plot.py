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

def set_ax_info(ax, xlabel, ylabel, zlabel='none', style='plain', title=None, legend=True):
    """Write title and labels on an axis with the correct fontsizes.

    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    if zlabel != 'none':
        ax.set_zlabel(ylabel, fontsize=20)
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
        FE = np.loadtxt(data_path + f"FE/single_n{k}.txt", unpack=True, delimiter=",", skiprows=1)
        RK4 = np.loadtxt(data_path + f"RK4/single_n{k}.txt", unpack=True, delimiter=",", skiprows=1)

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



def test_double_particle(scheme="RK4"):

    tRUa = np.loadtxt(data_path + f"{scheme}/double_without.txt", unpack=True, delimiter=",", skiprows=1)
    tRUb = np.loadtxt(data_path + f"{scheme}/double_with.txt", unpack=True, delimiter=",", skiprows=1)

    Nt = int(len(tRUa[0])/2)
    t = tRUa[0,:Nt]
    # without interactions:
    Ra = np.zeros((Nt,3,2))
    Ra[:,:,0] = tRUa[1:4, :Nt].T
    Ra[:,:,1] = tRUa[1:4, Nt:].T
    Ua = np.zeros((Nt,3,2))
    Ua[:,:,0] = tRUa[4:7, :Nt].T
    Ua[:,:,1] = tRUa[4:7, Nt:].T
    # with interactions
    Rb = np.zeros((Nt,3,2))
    Rb[:,:,0] = tRUb[1:4, :Nt].T
    Rb[:,:,1] = tRUb[1:4, Nt:].T
    Ub = np.zeros((Nt,3,2))
    Ub[:,:,0] = tRUb[4:7, :Nt].T
    Ub[:,:,1] = tRUb[4:7, Nt:].T

    # defining some colours for consistency (fix):
    cmap = ['copper', 'bone'] 
    colourline = ['navy', 'darkorange']
    colourpoint = ['yellow', 'red']

    RR = [Ra, Rb]
    UU = [Ua, Ub] 
    titles = ['Without interactions', 'With interactions']
    
    
    '''Plot in xy plane'''

    def trajectory_plane(axes, coords):
        ii, jj = coords
        if ii == 0:
            x = 'x'
        elif ii == 1:
            x = 'y'
        else:
            x = 'z'
        
        if jj == 0:
            y = 'x'
        elif jj == 1:
            y = 'y'
        else:
            y = 'z'


        for k, ax in enumerate(axes):
            R = RR[k]; U = UU[k]
            for p in [0,1]: 
                ax.scatter(R[:,ii,p], R[:,jj,p], s=3, marker='o', c=t, cmap=cmap[p], alpha=.7)

            for p in [0,1]:
                ax.plot(R[0,ii,p],  R[0,jj,p],  marker="P", ms=12, c=colourpoint[p], alpha=1)
                ax.plot(R[-1,ii,p], R[-1,jj,p], marker="*", ms=12, c=colourpoint[p], alpha=1)
            
            ax.set_aspect('equal')
            set_ax_info(ax, r'$%s$ [$\mu$m]'%x, r'$%s$ [$\mu$m]'%y, title=titles[k], legend=False)

        return axes


    fig, axes = plt.subplots(ncols=2, layout='constrained', sharey=True)
    ax1, ax2 = axes.flat
    axes = trajectory_plane([ax1, ax2], (0,1))

    plt.show() # savepush


    '''Make phase plots'''

    
    def phase_plot(axes, coord):
        if coord == 0:
            x = 'x'
        elif coord == 1:
            x = 'y'
        else:
            x = 'z'
        
        for k, ax in enumerate(axes):
            R = RR[k]; U = UU[k]
            for p in [0,1]:
                #ax.plot(R[:,coord,p], U[:,coord,p], lw=.7, c=colourline[p], ls='-',  alpha=.5)
                ax.scatter(R[:,coord,p], U[:,coord,p], s=3, marker='o', c=t, cmap=cmap[p], alpha=.7)

            for p in [0,1]:
                ax.plot(R[0,coord,p],  U[0,coord,p],  marker="P", ms=12, c=colourpoint[p], alpha=1)
                ax.plot(R[-1,coord,p], U[-1,coord,p], marker="*", ms=12, c=colourpoint[p], alpha=1)
            
            set_ax_info(ax, r'$%s$ [$\mu$m]'%x, r'$v_{%s}$ [$\mu$m/s]'%x, title=titles[k], legend=False)
        return axes


    fig, axes = plt.subplots(ncols=2, layout='constrained', sharey=True)
    axes = phase_plot(axes.flat, 0)
    plt.show()

    fig, axes = plt.subplots(ncols=2, layout='constrained', sharey=True)
    axes = phase_plot(axes.flat, 2)
    plt.show()


    '''3D trajectories'''

    fig, axes = plt.subplots(ncols=2, layout='constrained', subplot_kw={'projection':'3d'})
    ax1, ax2 = axes.flat
    axes = [ax1, ax2]

    for k, ax in enumerate(axes):
        R = RR[k]
        for p in [0, 1]:
            ax.scatter(R[:,0,p], R[:,1,p], R[:,2,p], s=3, marker='o', c=t, cmap=cmap[p], alpha=.7)

        for p in [0,1]:
            ax.plot(R[ 0,0,p], R[ 0,1,p], R[ 0,2,p], marker="P", ms=12, c=colourpoint[p], alpha=1)
            ax.plot(R[-1,0,p], R[-1,1,p], R[-1,2,p], marker="*", ms=12, c=colourpoint[p], alpha=1)
        
        #ax.set_aspect('equal')
        # maybe set xlims etc.
        set_ax_info(ax, xlabel=r'$x$ [$\mu$m]', ylabel=r'$y$ [$\mu$m]', zlabel=r'$z$ [$\mu$m]', title=titles[k], legend=False)
    
    plt.show()



def first_with_timedep():
    # just for testing

    tRU = np.loadtxt(data_path + f"RK4/first.txt", unpack=True, delimiter=",", skiprows=1)

    t_ = tRU[0]
    for i in range(len(t_)):
        if np.abs(t_[i+1]-0)<1e-8:
            idx_stop = i
            break
    
    t = t_[:idx_stop+1]
    Nt = len(t)
    Np = int(len(t_)/Nt)

    R = np.zeros((Nt,3,Np))
    U = np.zeros((Nt,3,Np))
    for p in range(Np):
        j = int(Nt*p)
        j_next = int(Nt*(p+1))
        if p == Np-1:
            R[:,:,p] = tRU[1:4, j:].T
            U[:,:,p] = tRU[4:7, j:].T
        else:
            R[:,:,p] = tRU[1:4, j:j_next].T
            U[:,:,p] = tRU[4:7, j:j_next].T


    cmaps = ['Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    cps = ['purple', 'b', 'g', 'darkorange', 'r', 'orange', 'orangered', 'yellow', 'm', 'pink', 'firebrick']
    fig, ax = plt.subplots(layout='constrained', subplot_kw={'projection':'3d'})
    j = 0
    for p in range(0,Np,Np//3):
        ax.scatter(R[:,0,p], R[:,1,p], R[:,2,p], s=3, marker='o', c=t, cmap=cmaps[j], alpha=.7)
        j+=1
    j=0
    for p in range(0,Np,Np//3):
        ax.plot(R[ 0,0,p], R[ 0,1,p], R[ 0,2,p], marker="P", ms=12, c=cps[j], alpha=1)
        ax.plot(R[-1,0,p], R[-1,1,p], R[-1,2,p], marker="*", ms=12, c=cps[j], alpha=1)
        j+=1
        
    set_ax_info(ax, xlabel=r'$x$ [$\mu$m]', ylabel=r'$y$ [$\mu$m]', zlabel=r'$z$ [$\mu$m]', title="Without interactions", legend=False)
    
    plt.show()


    

def trapped_particles(amplitudes={1:0.1, 2:0.4, 3:0.7}, scheme="RK4"):
    # f1, f2, f3 = 0.1, 0.4, 0.7
    # simulation duration T = 500 micro seconds
    # n_steps = 8000


    fig, ax = plt.subplots(layout="constrained")
    for f_no in amplitudes:
        f = amplitudes[f_no]
        omega, Ntrapped = np.loadtxt(data_path + f"{scheme}/trapped_f{f_no}.txt", unpack=True, delimiter=",")
        ax.plot(omega, Ntrapped, 'o-', lw=2, alpha=.7, label=r"$f_{%i} = %.1f$"%(f_no, f))

    set_ax_info(ax, xlabel=r'$\omega_V$ [MHz]', ylabel=r"$N_\mathrm{trapped}$")

    plt.show()


def trapped_particles_zoom(amplitudes={1:0.1}, scheme="RK4"):
    # f1, f2, f3 = 0.1, 0.4, 0.7
    # simulation duration T = 500 micro seconds
    # n_steps = 8000

    fig, ax = plt.subplots(layout="constrained")
    for f_no in amplitudes:
        f = amplitudes[f_no]
        omega, Ntrapped = np.loadtxt(data_path + f"{scheme}/trapped_f{f_no}_with.txt", unpack=True, delimiter=",")
        ax.plot(omega, Ntrapped, 'o-', lw=2, alpha=.7, label=r"$f_{%i} = %.1f$"%(f_no, f))

    set_ax_info(ax, xlabel=r'$\omega_V$ [MHz]', ylabel=r"$N_\mathrm{trapped}$")

    plt.show()


def check_upper_and_lower_bound():
    # Removed later??
    RK4_ = np.loadtxt(data_path + f"tests/RK4/single_n4.txt", unpack=True, delimiter=",", skiprows=1)

    # print(RK4_.shape)
    t, x, y = RK4_[0:3]
    v0 = RK4_[5][0]
    x0 = x[0]
    m = 40.078 
    B0 = 96.5 
    V_d2 = 9.65 
    omegaz2 = 2*V_d2 /m 
    omega_0 = B0/m 

    omega_minus = omega_0/2 - np.sqrt(omega_0**2 - 2*omegaz2)/2 
    omega_plus  = omega_0/2 + np.sqrt(omega_0**2 - 2*omegaz2)/2 

    A_plus  = (v0 + omega_minus*x0) / (omega_minus - omega_plus)
    A_minus = - (v0 + omega_plus * x0) / (omega_minus - omega_plus)

    R_minus = np.abs(A_plus - A_minus)
    R_plus = A_plus + A_minus

    # R_minus = (2*v0 + x[0]*omega_0) / np.sqrt(omega_0**2 - 2*omegaz2)
    # R_plus  = x[0] 

    theta = np.linspace(0,2*np.pi,int(1e3))
    f_squared = A_plus**2 + A_minus**2 + 2*np.cos((omega_plus-omega_minus)*theta)
    R_t = np.sqrt(f_squared)

    plt.plot(x,y)
    plt.plot(R_plus*np.cos(theta), R_plus*np.sin(theta), label=r'$R_{+}$')
    plt.plot(R_minus*np.cos(theta), R_minus*np.sin(theta), label=r'$R_{-}$')
    plt.axis('equal')
    plt.legend(fontsize=20)
    plt.show()


if __name__=="__main__":
    test_single_particle()

    # test_double_particle()

    #first_with_timedep()
    # trapped_particles()
    # trapped_particles_zoom()
