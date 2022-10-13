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
    Euler = np.loadtxt(data_path + "Euler_N1.txt", unpack=True, delimiter=",", skiprows=1)
    RK4   = np.loadtxt(data_path + "RK4_N1.txt", unpack=True, delimiter=",", skiprows=1)

    

    t = Euler[0]
    zE = Euler[3]
    zR = RK4[3]

    omega_z = np.sqrt(2*1/40 * 9.65)
    z_anal = 20*np.cos(omega_z * t)


    plt.plot(t, z_anal, lw=10, alpha=0.3, color='green', label='analytical')
    plt.plot(t, zE, lw=2, color='red', label='Euler')
    plt.plot(t, zR, '--', color='blue', label='RK4')
    plt.legend()
    plt.show()

    plt.title('abs difference')
    plt.plot(t, np.abs(z_anal - zE), '--', color='red', label='Euler')
    plt.plot(t, np.abs(z_anal - zR), ':', color='blue', label='RK4')
    plt.legend()
    plt.show()

def test_double_particle():
    Euler = np.loadtxt(data_path + "Euler_N2.txt", unpack=True, delimiter=",", skiprows=1)
    RK4   = np.loadtxt(data_path + "RK4_N2.txt", unpack=True, delimiter=",", skiprows=1)


    Nt = int(len(Euler[0])/2)
    t = Euler[0,:Nt]
    rE = np.zeros((Nt,3,2))
    rE[:,:,0] = Euler[1:4, :Nt].T
    rE[:,:,1] = Euler[1:4, Nt:].T
    rR = np.zeros((Nt,3,2))
    rR[:,:,0] = RK4[1:4, :Nt].T
    rR[:,:,1] = RK4[1:4, Nt:].T

    vE = np.zeros((Nt,3,2))
    vE[:,:,0] = Euler[4:7, :Nt].T
    vE[:,:,1] = Euler[4:7, Nt:].T
    vR = np.zeros((Nt,3,2))
    vR[:,:,0] = RK4[4:7, :Nt].T
    vR[:,:,1] = RK4[4:7, Nt:].T


    fig, axes = plt.subplots(ncols=2)
    ax1, ax2 = axes.flat

    ax1.plot(rE[:,0], vE[:,0], color='red', label='Euler', alpha=.6)
    ax1.plot(rR[:,0], vR[:,0], color='blue', label='RK4')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$v_x$')

    ax2.plot(rE[:,2], vE[:,2], color='red', label='Euler', alpha=.6)
    ax2.plot(rR[:,2], vR[:,2], color='blue', label='RK4')
    ax2.set_xlabel(r'$z$')
    ax2.set_ylabel(r'$v_z$')


    plt.show()


if __name__=="__main__":

    test_double_particle()
