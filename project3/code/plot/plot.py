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

make_png = True # temporary for saving pngs


def save_push(fig, pdf_name, save=True, push=False, show=False, tight=True, png_duplicate=make_png):
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
    pdfname = pdf_name.replace('.pdf', '').strip() + ".pdf"
    file = plot_path + pdfname
    if save:
        print(f'Saving plot: {file}')
        fig.savefig(file)
        if png_duplicate:
            fig.savefig(plot_path + "png/" + pdfname.replace('.pdf', '.png'))
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

def set_ax_info(ax, xlabel, ylabel, zlabel='none', style='plain', title=None, legend=True, loc='best'):
    """Write title and labels on an axis with the correct fontsizes.

    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel, fontsize=20)
    if ylabel != False:
        ax.set_ylabel(ylabel, fontsize=20)
    if zlabel != 'none':
        ax.set_zlabel(zlabel, fontsize=20)
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.yaxis.get_offset_text().set_fontsize(15)
    try:
        ax.ticklabel_format(style=style)
    except AttributeError:
        pass
    if legend:
        ax.legend(fontsize=15, loc=loc)


def z_analytical(z_anal, z_FE, z_RK, t, savepush=False):

    fig, ax = plt.subplots()
    ax.plot(t, z_anal, lw=10, alpha=0.5, color='green', label='Analytical')
    ax.plot(t, z_FE, color='red', lw=3, ls='solid', label='Forward Euler')
    ax.plot(t, z_RK, color='blue', lw=2, ls='--', label='Runge Kutta 4')
    ax.set_xlim(t[0], t[-1]+12)
    set_ax_info(ax, r"$t\,[\mathrm{\mu s}]$", r"$z\,[\mathrm{\mu m}]$", loc='upper right')
    if savepush:
        save_push(fig, pdf_name="compare_z_analytical", push=True)
    else:
        plt.show()


def two_particles_plane(p_noint, p_int, xlabel, ylabel, fname, savepush=False):
    """
    To be made prettier...
    """
    p1_noint, p2_noint = p_noint 
    p1_int, p2_int = p_int 

    fig, ax = plt.subplots(ncols=2, sharey='all', layout='constrained')
    
    ax[0].scatter(*p1_noint, s=3, marker='o', c="navy", alpha=.7)
    ax[0].scatter(*p2_noint, s=3, marker='o', c="orangered", alpha=.7)
    ax[1].scatter(*p1_int, s=3, marker='o', c="navy", alpha=.7)
    ax[1].scatter(*p2_int, s=3, marker='o', c="orangered", alpha=.7)
    
    ax[0].plot(*p1_noint.T[0], marker="P", ms=12, c="b")
    ax[0].plot(*p2_noint.T[0], marker="P", ms=12, c="y")
    ax[0].plot(*p1_noint.T[-1], marker="*", ms=12, c="b")
    ax[0].plot(*p2_noint.T[-1], marker="*", ms=12, c="y")

    ax[1].plot(*p1_int.T[0], marker="P", ms=12, c="b")
    ax[1].plot(*p2_int.T[0], marker="P", ms=12, c="y")
    ax[1].plot(*p1_int.T[-1], marker="*", ms=12, c="b")
    ax[1].plot(*p2_int.T[-1], marker="*", ms=12, c="y")
    

    ax[0].set_aspect('equal')
    ax[1].set_aspect('equal')
    xlim = ax[0].get_xlim()
    ylim = ax[0].get_ylim()
    ax[0].set_xlim(xlim)
    ax[0].set_ylim(ylim)
    ax[1].set_xlim(xlim)
    ax[1].set_ylim(ylim)

    set_ax_info(ax[0], xlabel, ylabel, title='No interaction', legend=False)
    set_ax_info(ax[1], xlabel, ylabel=False, title='With interaction', legend=False)


    if savepush:
        save_push(fig, fname, push=True)
    else:
        plt.show()


def two_particles_3d(p_noint, p_int, fname, savepush=False):
    p1_noint, p2_noint = p_noint 
    p1_int, p2_int = p_int 

    fig, (ax1, ax2) = plt.subplots(ncols=2,subplot_kw={'projection':'3d'})

    axes = [ax1, ax2]
    
    # plot trajectories:
    ax1.scatter(*p1_noint, s=3, marker='o', c="navy", alpha=.7)
    ax1.scatter(*p2_noint, s=3, marker='o', c="orangered", alpha=.7)
    ax2.scatter(*p1_int, s=3, marker='o', c="navy", alpha=.7)
    ax2.scatter(*p2_int, s=3, marker='o', c="orangered", alpha=.7)

    ax1.plot(*p1_noint.T[0], marker="P", ms=12, c="b")
    ax1.plot(*p2_noint.T[0], marker="P", ms=12, c="y")
    ax1.plot(*p1_noint.T[-1], marker="*", ms=12, c="b")
    ax1.plot(*p2_noint.T[-1], marker="*", ms=12, c="y")

    ax2.plot(*p1_int.T[0], marker="P", ms=12, c="b")
    ax2.plot(*p2_int.T[0], marker="P", ms=12, c="y")
    ax2.plot(*p1_int.T[-1], marker="*", ms=12, c="b")
    ax2.plot(*p2_int.T[-1], marker="*", ms=12, c="y")
    

    xlabel = r"$x\,[\mathrm{\mu m}]$"
    ylabel = r"$y\,[\mathrm{\mu m}]$"
    zlabel = r"$z\,[\mathrm{\mu m}]$"
    
    set_ax_info(ax2, xlabel, ylabel, zlabel, title='With interaction', legend=False)
    set_ax_info(ax1, xlabel, ylabel, zlabel, title='No interaction', legend=False)

    fig.subplots_adjust(left=0.01, bottom=0.02, right=0.96, top=0.98, wspace=0.08, hspace=0.01)
    if savepush:
        save_push(fig, fname, push=True, tight=False)
    else:
        plt.show()


def error_plot(errors, times, fname, title, savepush=False):
    
    colours = ['dodgerblue', 'olive', 'darkorange', 'navy']
    markersizes = [5.5, 5, 4.5, 4]
    
    fig, ax = plt.subplots()
    for i in range(len(errors)):
        n =  len(errors[i]) - 1
        ax.plot(times[i][1:], errors[i][1:], 'o-', ms=markersizes[i], c=colours[i], label=r"$n_{:.0f} = {:.0f}$".format(i+1, n))

    xlabel=r'$t$ [$\mu$s]'
    ylabel=r'Relative error size'# for $n_k$ steps'
    set_ax_info(ax, xlabel, ylabel, title=title)
    ax.legend(fontsize=20)
    ax.set_yscale("log")
    
    if savepush:
        save_push(fig, fname, push=True)
    else:
        plt.show()


def plot_trapped_coarse(N_trapped, omega_V, f_values, fname, title, savepush=False):
    fig, ax = plt.subplots()
    markersizes = [6, 5, 4] 
    colours = ["#4C72B0", "#55A868", "#C44E52"]
    for i in range(len(N_trapped)):
        ax.plot(omega_V, N_trapped[i], 'o-', ms=markersizes[i], c=colours[i], alpha=.7, label=rf"$f=\,${f_values[i]:.1f}")

    xlabel = r"$\omega_V\,[\mathrm{MHz}]$"
    ylabel = r"$N_{\mathrm{trapped}} / N_\mathrm{p}$"
    set_ax_info(ax, xlabel, ylabel, title=title)
    ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1]+ 0.3)
    ax.legend(fontsize=20, loc='lower right')

    
    if savepush:
        save_push(fig, fname, push=True)
    else:
        plt.show()


def plot_trapped_fine(N_trapped, omega_V, legs, fname, title, savepush=False):
    fig, ax = plt.subplots()
    colours = ["#4C72B0", "#55A868", "#C44E52"]
    # assuming f = f1
    ax.plot(omega_V, N_trapped[0],'o-', ms=6, c=colours[0], label=legs[0], alpha=.5)
    ax.plot(omega_V, N_trapped[1],'o--', lw=2, ms=9, c=colours[0], label=legs[1])
    xlabel = r"$\omega_V\,[\mathrm{MHz}]$"
    ylabel = r"$N_{\mathrm{trapped}} / N_\mathrm{p}$"
    set_ax_info(ax, xlabel, ylabel, title=title)
    ax.legend(fontsize=20, loc='lower right')

    if savepush:
        save_push(fig, fname, push=True)
    else:
        plt.show()










if __name__=="__main__":
    print("Use the analysis script. We only hide the ugly stuff here")
    # test_single_particle()

    # test_double_particle()

    #first_with_timedep()
    # trapped_particles()
    # trapped_particles_zoom()
