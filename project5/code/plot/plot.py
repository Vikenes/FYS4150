import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import os

plt.style.use("seaborn")

#   Paths
here = os.path.abspath(".")
binfile_path = here + "/../../output/binfiles/"
latex_path = here + "/../../latex/"
temp_path = here + "/../../output/plots/temp/"
plot_path = here +"/../../output/plots/pdfs/"
video_path = here +"/../../output/videos/"

#   rc and plot params
TICKLABELSIZE = 25
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


#   Global setting commands
temp = True     # makes temporary .png files instead of .pdf
SAVE = False
PUSH = False
SHOW = True

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
            png_fname = temp_path + pdfname.replace('.pdf', '.png')
            fig.savefig(png_fname)
            os.system(f"git add {png_fname}")
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
    
