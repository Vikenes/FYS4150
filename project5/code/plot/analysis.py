#Please place all analysis here

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import plot as PLOT
import pyarma as pa
import pandas as pd
from IPython import embed
from matplotlib.animation import FuncAnimation

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

# embed()

# U_noslit = pa.cx_cube()
# U_noslit.load("../../output/binfiles/no_slit_arma_cube.bin")
# # # U_noslit.load("../no_slit_arma_cube")
# Uns = np.asarray(U_noslit)

# U_double_slit = pa.cx_cube()
# U_double_slit.load("../../output/binfiles/double_slit_arma_cube.bin")
# # # # U_double_slit.load("../double_slit_arma_cube")
# Uds = np.asarray(U_double_slit)


# # # pns = Uns.conj()*Uns
# pns = np.abs(Uns)**2
# pds = np.abs(Uds)**2
# # pds = Uds.conj()*Uds

# pns = pns[0,:,:]
# print(pns.sum())

# embed()



class Analysis:
    def __init__(self, filename="no_slit_arma_cube", title="No slit"):
        U_arma = pa.cx_cube()
        U_arma.load(binfile_path+filename+".bin")
        self.U = np.asarray(U_arma)
        self.P = np.abs(self.U)**2
        self.title = title


    def animate(self, T=0.008, dt=2.5e-5):
        # Set up a 2D xy grid
        h = 0.005
        x_points = np.arange(0, 1+h, h)
        y_points = np.arange(0, 1+h, h)
        x, y = np.meshgrid(x_points, y_points, sparse=True)

        # Array of time points
        t_points = np.arange(0, T+dt, dt)

        # A function for a Gaussian that is travelling 
        # in the x direction and broadening as time passes
        def z(x,y,t):
            v = 0.5
            x_c = 0.2
            sigma_x = 0.025 + 0.15 * t
            return 1. / (2 * np.pi * np.sqrt(sigma_x)) * np.exp(-0.5 * (x - x_c - v * t)**2 / sigma_x**2)

        # Fill z_data_list with f(x,y,t)
        z_data_list = []

        for image in self.P:
                z_data_list.append(image.T)


        # Some settings
        fontsize = 12
        t_min = t_points[0]
        x_min, x_max = x_points[0], x_points[-1]
        y_min, y_max = y_points[0], y_points[-1]

        # Create figure
        fig = plt.figure(figsize=(12, 7))
        ax = plt.gca()

        # Create a colour scale normalization according to the max z value in the first frame
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

        # Plot the first frame
        img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

        # Axis labels
        plt.xlabel("x", fontsize=LABELSIZE)
        plt.ylabel("y", fontsize=LABELSIZE)
        plt.xticks(fontsize=TICKLABELSIZE)
        plt.yticks(fontsize=TICKLABELSIZE)
        plt.title(self.title, fontsize=TITLESIZE)
        plt.grid(False)

        # Add a colourbar
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label("z(x,y,t)", fontsize=LABELSIZE)
        cbar.ax.tick_params(labelsize=TICKLABELSIZE)

        # Add a text element showing the time
        time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                            horizontalalignment="right", verticalalignment="top", fontsize=LEGENDSIZE)

        # Function that takes care of updating the z data and other things for each frame
        def animation(i):
            # Normalize the colour scale to the current frame?
            norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
            img.set_norm(norm)

            # Update z data
            img.set_data(z_data_list[i])

            # Update the time label
            current_time = t_min + i * dt
            time_txt.set_text("t = {:.3e}".format(current_time))

            return img

        # Use matplotlib.animation.FuncAnimation to put it all together
        anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=True, blit=0)

        # Run the animation!
        plt.show()

        # # Save the animation
        # anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)

NOSLIT = Analysis()
DSLIT = Analysis("double_slit_arma_cube", "Double Slit")

NOSLIT.animate()
DSLIT.animate()
