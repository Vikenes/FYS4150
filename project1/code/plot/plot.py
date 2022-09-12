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

# other rc parameters
plt.rc('figure', figsize=(12,7))

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
    try:
        ax.ticklabel_format(style=style)
    except AttributeError:
        pass
    if legend:
        ax.legend(fontsize=15)


def line_plot_txt(data_file, pdf_name="none"):
    """Plot a function f(x) from a txt file
    Args:
        datafile: txt file formatted as x, f(x)
        pdf_name: name of the resulting plot
                    saved as: .../plots/pdf_name.pdf
    """

    # Plot
    x, v = np.loadtxt(data_path + data_file, unpack=True, delimiter=",")
    fig, ax = plt.subplots()
    ax.plot(x,v, lw=3, c='powderblue')

    title = r"Analytical solution of $u(x)$"
    xlabel = r"$x$"
    ylabel = r"$u(x)$"
    set_ax_info(ax, xlabel, ylabel, title=title, legend=False)

    # Option to save, push and show resulting plot
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


def compare_plots_txt(data_files, pdf_name='none'):
    """Plot a function f(x) from a txt file
    Args:
        datafiles: txt files formatted as x, f(x), first one the analytical solution
        pdf_name: name of the resulting plot
                    saved as: .../plots/pdf_name.pdf
    """

    # Plot
    fig, axes = plt.subplots(1,2, sharey=True)
    ax1, ax2 = axes
    x, u = np.loadtxt(data_path + data_files[0], unpack=True, delimiter=",")
    ax1.plot(x, u, lw=10, c='powderblue', label=r"$u(x)$")
    ax2.plot(x, u, lw=10, c='powderblue', label=r"$u(x)$")

    i = 0
    colors = ['red', 'magenta', 'darkorange']
    for ax, file in zip([ax1, ax1, ax2], data_files[1:]):
        x, v = np.loadtxt(data_path + file, unpack=True, delimiter=",")
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
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


def plot_error(data_files, pdf_name='none', relative_error=False):
    fig, ax = plt.subplots(figsize=(12,7))
    i = 0
    colors = ['red', 'magenta', 'darkorange']
    for file in data_files:
        x, logeps = np.loadtxt(data_path + file, unpack=True, delimiter=",")
        n = x.size+1
        ax.plot(x, logeps, 'o', markersize=10/n**(1/8), c=colors[i], label=r"$n_{\mathrm{steps}}=%i$"%n)
        i+=1
    xlabel = r"$x$"
    if relative_error:
        ylabel = r"$\log_{10}(\epsilon)$"
        title = r'Logarithm of relative error $\epsilon$ (excluding end points)'
        ax.set_ylim(-5.5, 0.5)
    else:
        ylabel = r"$\log_{10}(\Delta)$"
        title = r'Logarithm of absolute error $\Delta$ (excluding end points)'
    set_ax_info(ax, xlabel, ylabel, title=title, legend=True)
    fig.tight_layout()
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)


def plot_max_error(data_file, pdf_name='none'):
    fig, ax = plt.subplots(figsize=(12,7))
    logsteps, logerror = np.loadtxt(data_path + data_file, unpack=True, delimiter=",")
    ax.plot(logsteps, logerror, '--', color='pink', linewidth=2)
    ax.plot(logsteps, logerror, '.', markersize=20, color="firebrick", label=r"$\mathrm{max}(\epsilon)$")
    xlabel = r'$\log_{10}(n_{\mathrm{steps}})$'
    ylabel = r'$\log_{10}(\epsilon)$'
    title = r'Logarithm of maximum relative error $\epsilon$'
    set_ax_info(ax, xlabel, ylabel, title=title)
    ax.set_yscale("log")
    ax.set_xscale("log")
    fig.tight_layout()
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)

def timed_algorithms(data_file, pdf_name='none', number_of_runs=500):
    n_steps, g_mean, s_mean, g_std, s_std = np.loadtxt(data_path + data_file, unpack=True, delimiter=",")

    fig, ax = plt.subplots()

    ax.fill_between(n_steps, g_mean-g_std, g_mean+g_std, color="navy", alpha=.2)
    ax.fill_between(n_steps, s_mean-s_std, s_mean+s_std, color="olive", alpha=.2)

    ax.errorbar(n_steps, g_mean, g_std, marker="o", color="navy", label="General Thomas", linestyle="")
    ax.errorbar(n_steps, s_mean, s_std, marker="o", color="olive", label="Special Thomas", linestyle="")

    xlabel = r'$n_{\mathrm{steps}}$'
    ylabel = r'$\delta t$ [s]'

    #ax.set_ylabel(r'Time used per run [s]', fontsize=20)
    title = r'Time spent per run $\delta t$ (averaged over %i runs)' %number_of_runs
    set_ax_info(ax, xlabel, ylabel, title=title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.tight_layout()
    if pdf_name=="none":
        plt.show()
    else:
        save_push(fig, pdf_name=pdf_name, push=True, save=True)




# Problem 2

line_plot_txt(data_file="analytical_x_u.txt", pdf_name='ux')


# Problem 7

files = ["analytical_x_u.txt"]
for n in [10, 100, 1000]:
    files.append(f"generalThomas_{n}steps.txt")

compare_plots_txt(files, "comparison_p7")

# Problem 8

rel_error_files = []
abs_error_files = []
for n in [10, 100, 1000]:
    rel_error_files.append(f'rel_error_{n}steps.txt')
    abs_error_files.append(f'abs_error_{n}steps.txt')

plot_error(abs_error_files, pdf_name="absolute_error")
plot_error(rel_error_files, pdf_name="relative_error", relative_error=True)

plot_max_error('max_rel_error.txt', pdf_name="max_relative_error")


# Problem 10
timed_algorithms("thomas_timed.txt", pdf_name="algorithms_timed")


#   SHOW PLOTS ONLY

# # Problem 2

# line_plot_txt(data_file="analytical_x_u.txt")


# # Problem 7

# files = ["analytical_x_u.txt"]
# for n in [10, 100, 1000]:
#     files.append(f"generalThomas_{n}steps.txt")

# compare_plots_txt(files)

# # Problem 8

# rel_error_files = []
# abs_error_files = []
# for n in [10, 100, 1000]:
#     rel_error_files.append(f'rel_error_{n}steps.txt')
#     abs_error_files.append(f'abs_error_{n}steps.txt')

# plot_error(abs_error_files)
# plot_error(rel_error_files, relative_error=True)

# plot_max_error('max_rel_error.txt')


# # Problem 10
# timed_algorithms("thomas_timed.txt")
