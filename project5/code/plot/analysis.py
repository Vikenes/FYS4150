#Please place all analysis here

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# import plot as PLOT
import pyarma as pa
import pandas as pd
from IPython import embed
from matplotlib.animation import FuncAnimation


# embed()

U_noslit = pa.cx_cube()
U_noslit.load("../../output/binfiles/no_slit_arma_cube.bin")
# # U_noslit.load("../no_slit_arma_cube")
Uns = np.asarray(U_noslit)

U_double_slit = pa.cx_cube()
U_double_slit.load("../../output/binfiles/double_slit_arma_cube.bin")
# # # U_double_slit.load("../double_slit_arma_cube")
Uds = np.asarray(U_double_slit)


# # pns = Uns.conj()*Uns
pns = np.abs(Uns)**2
pds = np.abs(Uds)**2
# # pds = Uds.conj()*Uds

# pns = pns[0,:,:]
# print(pns.sum())

# embed()




#
# Let's generate a dummy time series for a function z(x,y,t)
#
def animate(ucube=None, T=1, dt=0.005):
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

    if ucube is not None:
        for image in ucube:
            z_data_list.append(image.T)
    else:
        for t in t_points:
            z_data = z(x, y, t)
            z_data_list.append(z_data)


    #
    # Now the list z_data_list contains a series of "frames" of z(x,y,t), 
    # where each frame can be plotted as a 2D image using imshow. Let's
    # animate it!
    #

    # Some settings
    fontsize = 12
    t_min = t_points[0]
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

    # Plot the first frame
    img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("z(x,y,t)", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

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

# animate(pns, T=0.008, dt = 2.5e-5)
animate(pds, T=0.008, dt = 2.5e-5)
