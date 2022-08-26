import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import seaborn as sns
import pandas as pd
from datetime import datetime


# The style we want
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('font', family='DejaVu Sans')

# folder paths 
here = os.path.abspath(".")
data_path = here + "/../../output/data/"
plot_path = here + "/../../output/plots/"


def save_push(fig, pdf_name, save=True, push=False, show=False):
    """
    This function handles wether you want to show,
    save and/or push the file to git.

    Args:
        fig (matplotlib.figure): Figure you want to handle
        pdfname (string): Name of output pdf file 
        args (argparse)
    """
    file = plot_path + pdf_name + ".pdf"
    if save:
        print(f'Saving plot: {file}')
        fig.savefig(file)
    if push:
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")
    if show:
        plt.show()
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
    ax.ticklabel_format(style=style)
    if legend:
        ax.legend(fontsize=15)


def line_plot_txt(data_file, pdf_name):
    """Plot a function f(x) from a txt file
    Args:
        datafile: txt file formatted as x, f(x)
        pdf_name: name of the resulting plot
                    saved as: .../plots/pdf_name.pdf
    """

    # Plot 
    x, v = np.loadtxt(data_path + data_file, unpack=True)
    fig, ax = plt.subplots()
    ax.plot(x,v)

    title = r"Analytical solution of $u(x)$"
    xlabel = r"$x$"
    ylabel = r"$u(x)$"
    set_ax_info(ax, xlabel, ylabel, title=title, legend=False)
    
    # Option to save, push and show resulting plot 
    save_push(fig, pdf_name=pdf_name, push=True, save=True)




# line_plot_txt(data_file="x_u.txt", pdf_name='ux')