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


def compare_plots_txt(data_files, pdf_name='none'):
    """Plot a function f(x) from a txt file
    Args:
        datafiles: txt files formatted as x, f(x), first one the analytical solution
        pdf_name: name of the resulting plot
                    saved as: .../plots/pdf_name.pdf
    """

    # Plot 
    fig, axes = plt.subplots(1,2, sharey=True, figsize=(12,7))
    ax1, ax2 = axes
    x, u = np.loadtxt(data_path + data_files[0], unpack=True)    
    ax1.plot(x, u, lw=10, c='powderblue', label=r"$u(x)$")
    ax2.plot(x, u, lw=10, c='powderblue', label=r"$u(x)$")
    
    i = 0 
    colors = ['red', 'magenta', 'darkorange']
    for ax, file in zip([ax1, ax1, ax2], data_files[1:]):
        x, v = np.loadtxt(data_path + file, unpack=True)
        n = x.size-1
        ax.plot(x, v, 'o', markersize=10/n**(1/8), c=colors[i], label=r"$n_{\mathrm{steps}}=%i$"%n)
        i +=1

    xlabel = r"$x$"
    ylabel = r"$u(x)$"
    title = "Numerical solution"
    fig.suptitle(title, fontsize=20)
    set_ax_info(ax1, xlabel, ylabel, legend=True)
    set_ax_info(ax2, xlabel, '', legend=True)
    fig.tight_layout()
    #fig.legend()
    
    # Option to save, push and show resulting plot 
    if pdf_name=="none":
        save_push(fig, pdf_name=pdf_name, show=True, push=False, save=False)
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


def plot_error(data_files, pdf_name='none', relative_error=False):
    fig, ax = plt.subplots(figsize=(12,7))
    i = 0 
    colors = ['red', 'magenta', 'darkorange']
    for file in data_files:
        x, error = np.loadtxt(data_path + file, unpack=True)
        n = x.size+1
        ax.plot(x, error, 'o', markersize=10/n**(1/8), c=colors[i], label=r"$n_{\mathrm{steps}}=%i$"%n)
        i+=1
    xlabel = r"x"
    if relative_error:
        ylabel = r"$\log_{10}(\epsilon)$"
        title = r'Logarithm of relative error $\epsilon$ (except end points)'
    else:
        ylabel = r"$\log_{10}(\Delta)$"
        title = r'Logarithm of absolute error $\Delta$ (except end points)'
    fig.suptitle(title, fontsize=20)
    set_ax_info(ax, xlabel, ylabel, legend=True)
    fig.tight_layout()
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)
    
    
def plot_max_error(data_file, pdf_name='none'):
    fig, ax = plt.subplots(figsize=(12,7))
    steps, error = np.loadtxt(data_path + data_file, unpack=True)
    ax.loglog(steps, error, '--', color='pink', linewidth=2)
    ax.loglog(steps, error, '.', markersize=15, color="red", label=r"$\mathrm{max}(\epsilon)$")
    ax.set_xlabel(r'$n_{\mathrm{steps}}$', fontsize=20)
    ax.set_ylabel(r'Relative error $\epsilon$', fontsize=20)
    title = r'Maximum relative error $\epsilon$'
    fig.suptitle(title, fontsize=20)
    ax.legend(fontsize=20)
    # set_ax_info(ax, xlabel, ylabel) #, title=title) # Need different style!
    fig.tight_layout()
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)

def timed_algorithms(pdf_name='none'):
    filename = "thomas_timed.txt"
    n_steps, g_mean, s_mean, gstd, sstd = np.loadtxt(data_path + filename, unpack=True, delimiter=",")
 

    fig, ax = plt.subplots(figsize=(12,7))
    
    # ax.plot(n_steps, g_mean, '.', markersize=15, color='red', label=r"General Thomas")
    # ax.plot(n_steps, s_mean, '.', markersize=15, color='blue', label=r"Special Thomas")
    ax.fill_between(n_steps, g_mean-gstd, g_mean+gstd, color="red", alpha=.2)
    ax.fill_between(n_steps, s_mean-sstd, s_mean+sstd, color="blue", alpha=.2)

    ax.errorbar(n_steps, g_mean, gstd, marker="o", color="red", label="General Thomas", linestyle="")
    ax.errorbar(n_steps, s_mean, sstd, marker="o", color="blue", label="Special Thomas", linestyle="")

    ax.set_xlabel(r'$n_{\mathrm{steps}}$', fontsize=20)
    ax.set_ylabel(r'Time used per run [s]', fontsize=20)
    ax.set_yscale("log")
    ax.set_xscale("log")
    title = r'Time spend per run (averaged over 500 runs)'
    fig.suptitle(title, fontsize=20)
    ax.legend(fontsize=20)
    fig.tight_layout()
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)

    


# Problem 2

line_plot_txt(data_file="x_u.txt", pdf_name='ux')


# Problem 7

files = ["x_u.txt"]
for n in [10, 100, 1000]:
    files.append(f"num_sol_{n}steps.txt")   

compare_plots_txt(files, "comparison_p7")

# Problem 8

relative_error_files = []
absolute_error_files = []
for n in [10, 100, 1000]:
    relative_error_files.append(f'relative_error{n}steps.txt')
    absolute_error_files.append(f'absolute_error{n}steps.txt')

plot_error(absolute_error_files, pdf_name="absolute_error")
plot_error(relative_error_files, pdf_name="relative_error", relative_error=True)
plot_max_error('max_relative_error.txt', pdf_name="max_relative_error")


# Problem 9 - testing

files = ["x_u.txt"]
for n in [10, 100, 1000]:
    files.append(f"special_num_sol_{n}steps.txt")

compare_plots_txt(files)

# Problem 10
timed_algorithms(pdf_name="thomas_timed")

