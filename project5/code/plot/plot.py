import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import os
from matplotlib.animation import FuncAnimation
from matplotlib import cm as colourmaps


plt.style.use("seaborn")

### Paths
here = os.path.abspath(".")
binfile_path = here + "/../../output/binfiles/"
latex_path = here + "/../../latex/"
temp_path = here + "/../../output/plots/temp/"
plot_path = here +"/../../output/plots/pdf/"
video_path = here +"/../../output/videos/"

### rc and plot params
# for figures that are small ?
TICKLABELSIZE = 22
LABELSIZE = 25
LEGENDSIZE = 20
TITLESIZE = 30
# for figures that are large ?
SMALLER_TICKLABELSIZE = 16
SMALLER_LABELSIZE = 18
SMALLER_LEGENDSIZE = 16
SMALLER_TITLESIZE = 22

### Set rc-params
plt.rc("legend", fontsize=LEGENDSIZE, fancybox=True, loc="best", frameon=True, edgecolor="black")
plt.rc("font", size=LABELSIZE)
plt.rc("axes", titlesize=TITLESIZE, labelsize=LABELSIZE)
plt.rc("xtick", labelsize=TICKLABELSIZE)
plt.rc("ytick", labelsize=TICKLABELSIZE)
# plt.rc("tickparams")
plt.rc("lines", linewidth=2.5)
plt.rc('text', usetex=True)
# plt.rc("label", fontsize=20)

plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = 'Times New Roman'


### Global setting commands
TEMP = True    # makes temporary .png files instead of .pdf
SAVE = True    # saves .pdf files
PUSH = False    # git stuff
SHOW = False    # show plots


### Functions for saving etc.

def save_push(fig, pdf_name, save=SAVE, push=PUSH, show=SHOW, tight=False, png_duplicate=TEMP):
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


### "Local" functions


def make_colourmap(ax, transposed_data, timestamp, cmap, norm, levels=False, spatial_extent=(0,1,0,1), stamp_colour="white"):
    if levels:
        img = ax.contourf(transposed_data, cmap=plt.get_cmap(cmap), norm=norm, extent=spatial_extent, levels=22)
    else:
        img = ax.imshow(transposed_data, extent=spatial_extent, cmap=plt.get_cmap(cmap), norm=norm)
    
    # img.set_norm(norm)
    time_txt = ax.text(0.95, 0.95, r"$t = %.3f$"%timestamp, color=stamp_colour, ha="right", va="top", fontsize=LEGENDSIZE)
    ax.grid(False)
    ax.set_aspect("equal")
    return img, ax

def default_mapfigure(timepoints, data, c_label=lambda t:"", cmap="gnuplot", num_maps=3, levels=False, vmin=None, vmax=None):
    fig, axes = plt.subplots(nrows=1, ncols=num_maps, sharey=True, figsize=(15, 10))

    axes.flat[0].set_ylabel(r"$y$")

    for j, ax in enumerate(axes.flat):
        v_min = vmin or np.min(data[j])
        v_max = vmax or np.max(data[j])
        norm = colourmaps.colors.Normalize(vmin=v_min, vmax=v_max)
        img, ax = make_colourmap(ax, data[j].T, timepoints[j], cmap, norm, levels)
        ax.set_xlabel(r"$x$")
        ax.tick_params("both", labelsize=SMALLER_TICKLABELSIZE)
        cbar = fig.colorbar(img, ax=ax, location="top", shrink=0.98)#, format="%3.2f")#, ticks=(v_min, (v_max-v_min)/2, v_max))
        # cbar.ax.tick_params(labelsize=SMALLER_TICKLABELSIZE, top=False, bottom=True, length=3, labelbottom=True, labeltop=False)
        cbar.ax.tick_params(labelsize=SMALLER_TICKLABELSIZE, top=True, bottom=False, length=3, labelbottom=False, labeltop=True, rotation=18)
        # cbar.set_label(c_label(timepoints[j]), fontsize=SMALLER_TITLESIZE)
        ax.set_title(c_label(timepoints[j]), fontsize=SMALLER_TITLESIZE)
    
    fig.subplots_adjust(hspace=0.12, wspace=0.12, left=0.06, right=0.98, bottom=0.06, top=0.96)

    return fig, axes.flat

def draw_walls(ax, yc_list, colour="palegreen"):
    xc = 0.5
    h = 0.05
    w = 0.02
    # quick fix
    # ymin = np.min(yc_list)
    # if len(yc_list) > 2:
    #     ymin = np.min(yc_list)+h/2
    #     ymax = np.max(yc_list)+h/2

    #     ax.add_patch(plt.Rectangle((xc-w/2, yc-h/2), w, h, fc=colour, ec=colour, lw=0.2, alpha=0.6, clip_on=False))

    #     for yc in yc_list[1:-1]:
    #         ax.add_patch(plt.Rectangle((xc-w/2, yc-h/2), w, h, fc=colour, ec=colour, lw=0.2, alpha=0.6, clip_on=False))
    # else:
    #     pass
    pass


### Functions to be called in 'analysis.py'

def animate_probability_density(t, P, title=None, wall_y=[], mp4name="animation", spatial_extent=(0,1,0,1), save=SAVE, show=SHOW):

    # Fill p_data
    p_data = []
    for image in P:
        p_data.append(image.T)

    # Create figure
    fig = plt.figure(figsize=(12, 9))
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = colourmaps.colors.Normalize(vmin=0.0, vmax=np.max(p_data[0]))

    # Plot the first frame
    img = ax.imshow(p_data[0], extent=spatial_extent, cmap=plt.get_cmap("gnuplot"), norm=norm)

    # Axis labels
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    if title is not None:
        ax.set_title(title)
    ax.grid(False)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(r"$p(\mathbf{x}; \, t)$")
    # cbar.ax.tick_params(labelsize=TICKLABELSIZE)

    # Add a text element showing the time
    time_txt = ax.text(0.95, 0.95, r"$t = %.3e$"%t[0], color="white", ha="right", va="top", fontsize=LEGENDSIZE)
    

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalise the colour scale to the current frame?
        norm = colourmaps.colors.Normalize(vmin=0.0, vmax=np.max(p_data[i]))
        img.set_norm(norm)

        # Update p data
        img.set_data(p_data[i])

        # Update the time label
        time_txt.set_text(r"$t = %.3e$" %t[i])

        return img
    draw_walls(ax, wall_y)

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(p_data), 2), repeat=True, blit=0)
    
    if show:
        # Run the animation!
        plt.show()
    if save:
        # Save the animation
        anim.save(video_path + mp4name.split(".")[0] + '.mp4', writer="ffmpeg", bitrate=-1, fps=30)

def snapshot_probability_density(t, P, Pmax=None, title=None, wall_y=None, pdfname="snapshot_P", spatial_extent=(0,1,0,1), vline=0.8, num_rows=1, save=SAVE, png_duplicate=TEMP, show=SHOW):

    label = lambda t: r"$p(\mathbf{x}; \, t\!=\!%.3f)$" %t
    fig, axes = default_mapfigure(t, P, label, cmap="gnuplot", num_maps=len(t))#, vmin=0)
    if wall_y is not None:
        for ax in axes:
            draw_walls(ax, wall_y)

    ### plot the screening line:
    ax = axes[-1]
    ax.axvline(vline, ls=':', lw=1.2, c="orangered", alpha=0.8)

    if save:
        fig.savefig(plot_path + pdfname.split(".")[0] + ".pdf")
    if png_duplicate:
        fig.savefig(temp_path + pdfname.split(".")[0] + ".png")
    if show:
        plt.show()
    
def snapshot_real_wavefunction(t, ReU, Umax=None, title=None, wall_y=None, pdfname="snapshot_ReU", spatial_extent=(0,1,0,1), num_rows=1, save=SAVE, png_duplicate=TEMP, show=SHOW):

    label = lambda t: r"$\mathrm{Re}\{u(t\!=\!%.3f, \mathbf{x})\}$" %t
    fig, axes = default_mapfigure(t, ReU, label, cmap="ocean", num_maps=len(t), levels=True)

    if wall_y is not None:
        for ax in axes:
            draw_walls(ax,wall_y)

    if save:
        fig.savefig(plot_path + pdfname.split(".")[0] + ".pdf")
    if png_duplicate:
        fig.savefig(temp_path + pdfname.split(".")[0] + ".png")
    if show:
        plt.show()

def snapshot_imaginary_wavefunction(t, ImU, Ulim=(None, None), title=None, wall_y=None, pdfname="snapshot_ImU", spatial_extent=(0,1,0,1), num_rows=1, save=SAVE, png_duplicate=TEMP, show=SHOW):
    
    label = lambda t: r"$\mathrm{Im}\{u(t\!=\!%.3f, \mathbf{x})\}$" %t
    fig, axes = default_mapfigure(t, ImU, label, cmap="ocean", num_maps=len(t), levels=True)
   
    if wall_y is not None:
        for ax in axes:
            draw_walls(ax, wall_y)


    if save:
        fig.savefig(plot_path + pdfname.split(".")[0] + ".pdf")
    if png_duplicate:
        fig.savefig(temp_path + pdfname.split(".")[0] + ".png")
    if show:
        plt.show()

def total_probability_deviation(t, P_tot_list, title=None, labels=None, pdfname="Ptot_deviation", save=SAVE, png_duplicate=TEMP, show=SHOW):

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    # ax.axhline(0, ls='--', lw=0.9, alpha=0.8, c="orangered")
    colours = ["dodgerblue", "violet", "seagreen"]
    for j, Ptot in enumerate(P_tot_list):
        ax.plot(t, np.abs(1-np.asarray(Ptot)), lw=2.5, c=colours[j], label=labels[j])
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|1-p^\mathrm{tot}(t)|$")
    # ax.set_ylabel(r"$\sum_{i, j} p(\mathbf{x}_{i,j}; \, t)-1$")
    if labels[0] is not None:
        ax.legend()
    
    fig.tight_layout()

    if save:
        fig.savefig(plot_path + pdfname.split(".")[0] + ".pdf")
    if png_duplicate:
        fig.savefig(temp_path + pdfname.split(".")[0] + ".png")
    if show:
        plt.show()

def probability_density_along_screen(y, p, x=0.8, t=0.002, title=None, label=None, xlabel=r"$y$", pdfname="P_along_sceen", save=SAVE, png_duplicate=TEMP, show=SHOW):
    fig, ax = plt.subplots()

    ax.plot(y, np.asarray(p), lw=2.5, c="orangered", label=label)
    ax.set_xlabel(xlabel) 
    ax.set_ylabel(r"$p_{x=%.1f}(y; \, t\!=\!%.3f)$"%(x, t))
    if label is not None:
        ax.legend()
    
    fig.tight_layout()
    if save:
        fig.savefig(plot_path + pdfname.split(".")[0] + ".pdf")
    if png_duplicate:
        fig.savefig(temp_path + pdfname.split(".")[0] + ".png")
    if show:
        plt.show()
    




def show_all():
    plt.show()