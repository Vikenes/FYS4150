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

def pt5(infiles, pdf_name="none"):
    N_t, M_t = np.loadtxt(data_path + infiles[0], unpack=True, delimiter=",")
    N_d, M_d = np.loadtxt(data_path + infiles[1], unpack=True, delimiter=",")
    fig, ax = plt.subplots()
    ax.plot(N_t, M_t, 'o', markersize=5, color="b", label=r"tridiagonal matrix $A$")
    ax.plot(N_d[::5], M_d[::5], '^', markersize=9, color="r", label=r"dense matrix $A$")
    title = r"Jacobi rotation method comparison"
    xlabel = r"$N$"
    ylabel = r"$M$"
    set_ax_info(ax, xlabel, ylabel, title=title, legend=True)
    # Option to save, push and show resulting plot
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


def pt6(infiles, pdf_name="none"):
    # infiles: (analytic), (Jacobi)
    xhat, v1, v2, v3 = np.loadtxt(data_path + infiles[0], unpack=True, delimiter=",")
    v_a = [v1, v2, v3]
    xhat, v1, v2, v3 = np.loadtxt(data_path + infiles[1], unpack=True, delimiter=",")
    v_J = [v1, v2, v3]
    
    c = ["b", "r", "g"]
    fig, ax = plt.subplots()

    n = len(xhat)-1
    ms = 25 * n**(-1/3)
    lw = 1.8
    for i, vi in enumerate(v_J):
        ax.plot(xhat, vi, 'o-', lw=lw, markersize=ms, color=c[i], label=r"$\mathbf{v}^{(%i)}$"%(i+1))

    for i, vi in enumerate(v_a):
        ax.plot(xhat, vi, 'o--', lw=lw/2, markersize=ms/2, color='white')

    
    ax.plot(xhat[0]-xhat[-1], 0, 'o-', lw=0.5, markersize=ms/2, color='white', label='analytic')
    ax.set_xlim(xhat[0], xhat[-1])  
    title = r"First three eigenvectors ($n=%i$)"%(n)
    xlabel = r"$\hat{x}$"
    ylabel = r"$v(\hat{x})$"
    set_ax_info(ax, xlabel, ylabel, title=title, legend=True)
    # Option to save, push and show resulting plot
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


if __name__=="__main__":

    infiles = ["transformations_per_tridiag_N_matrix.txt", "transformations_per_dense_N_matrix.txt"]
    pt5(infiles, "jacobi_comparison")

    infiles = ["analytical_solution_10steps.txt", "Jacobi_solution_10steps.txt"]
    pt6(infiles, "solution_10steps")

    infiles = ["analytical_solution_100steps.txt", "Jacobi_solution_100steps.txt"]
    pt6(infiles, "solution_100steps")