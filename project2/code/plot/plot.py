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


def save_push(fig, pdf_name, save=True, push=False, show=False):
    """
    This function handles wether you want to show,
    save and/or push the file to git.

    Args:
        fig (matplotlib.figure): Figure you want to handle
        pdfname (string): Name of output pdf file
        args (argparse)
    """
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

def pt5a(infiles, pdf_name="none"):
    N_t, it_t = np.loadtxt(data_path + infiles[0], unpack=True, delimiter=",")
    N_d, it_d = np.loadtxt(data_path + infiles[1], unpack=True, delimiter=",")
    fig, ax = plt.subplots()
    ax.plot(N_t, it_t, lw=1, ls="-", color="blue", label=r"Tridiagonal matrix $A$")
    ax.plot(N_d, it_d, lw=1, ls="-", color="red", label=r"Dense matrix $A$")
    title = r"Jacobi rotation method comparison"
    xlabel = r"Size of matrix $A$: $N$"
    ylabel = r"Nr. iterations needed for convergence"
    set_ax_info(ax, xlabel, ylabel, title=title, legend=True)
    # Option to save, push and show resulting plot
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


if __name__=="__main__":
    infiles = ["iterations_per_N_tridiag_matrix.txt", "iterations_per_dense_N_matrix.txt"]
    pt5a(infiles, pdf_name="jacobi_comparison")
