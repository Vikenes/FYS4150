from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
import seaborn as sns 
import pandas as pd 
import os 
import plot 

here = os.path.abspath(".")
data_path = here + "/../../output/data/"
FE_path = data_path + "FE/"
RK_path = data_path + "RK4/"


m = 40.078
B0 = 96.5
V_d2 = 9.65 

omega_z = np.sqrt(2 * V_d2 / m)
omega_0 = B0/m 

omega_plus  = (omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2)) / 2
omega_minus = (omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2)) / 2


def load(file, folder=RK_path, skiprows=1):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, delimiter=",", skiprows=skiprows)
    


def analytical_z(t, z0):
    return z0 * np.cos(omega_z * t)


def analytical_xy(t, x0, v0):
    denom = omega_minus - omega_plus
    A_plus = (v0 + omega_minus*x0) / denom 
    A_minus = - (v0 + omega_plus*x0) / denom 

    f_plus = A_plus * np.exp(-1j*omega_plus*t)
    f_minus = A_minus * np.exp(-1j*omega_minus*t)
    f = f_plus + f_minus 

    x = np.real(f)
    y = np.imag(f)

    return x,y
    




def compare_z_analytical(push=False):
    """
    Sanity check.
    """
    FE = load("single_n4.txt", FE_path)
    RK = load("single_n4.txt")

    t = FE[0]
    z_FE = FE[3]
    z_RK = RK[3]
    z0 = z_FE[0]

    z_anal = analytical_z(t, z0)#z0 * np.cos(omega_z * t)

    plot.z_analytical(z_anal, z_FE, z_RK, t, push)
   
    return None 



def compute_errors():
    RK = load("single_n3.txt")
    t,x,y,z,vx,vy,vz = RK

    v0 = vy[0]
    x0 = x[0]

    x_anal, y_anal = analytical_xy(t, x0, v0)
    plt.plot(x_anal, y_anal, color='red')
    plt.plot(x,y, '--')
    plt.axis('equal')
    plt.show()

    return None 



def xy_plane_movements(push=False):
    # RK4 only 
    p1_no_int = load("double_without_p1.txt")
    p2_no_int = load("double_without_p2.txt")
    p1_int = load("double_with_p1.txt")
    p2_int = load("double_with_p2.txt")

    p1xy = p1_no_int[1:3]
    p2xy = p2_no_int[1:3]

    p1xy_int = p1_int[1:3]
    p2xy_int = p2_int[1:3]

    fname = "xy_two_particles"
    xlabel = r"$x\,[\mathrm{\mu m}]$"
    ylabel = r"$y\,[\mathrm{\mu m}]$"
    plot.two_particles_plane((p1xy, p2xy), (p1xy_int, p2xy_int), xlabel, ylabel, fname=fname, push=push)

def x_phase_plot(push=False):
    p1_no_int = load("double_without_p1.txt")
    p2_no_int = load("double_without_p2.txt")
    p1_int = load("double_with_p1.txt")
    p2_int = load("double_with_p2.txt")



    p1x = p1_no_int[1::3]
    p2x = p2_no_int[1::3]
    p1x_int = p1_int[1::3]
    p2x_int = p2_int[1::3]

    fname = "x_phase_plot"
    xlabel = r"$x\,[\mathrm{\mu m}]$"
    ylabel = r"$v_x\,[\mathrm{\mu m / \mu s}]$"
    plot.two_particles_plane((p1x, p2x), (p1x_int, p2x_int), xlabel, ylabel, fname=fname, push=push)


def z_phase_plot(push=False):
    p1_no_int = load("double_without_p1.txt")
    p2_no_int = load("double_without_p2.txt")
    p1_int = load("double_with_p1.txt")
    p2_int = load("double_with_p2.txt")



    p1x = p1_no_int[1::3]
    p2x = p2_no_int[1::3]
    p1x_int = p1_int[1::3]
    p2x_int = p2_int[1::3]

    fname = "x_phase_plot"
    xlabel = r"$x\,[\mathrm{\mu m}]$"
    ylabel = r"$v_x\,[\mathrm{\mu m / \mu s}]$"
    plot.two_particles_plane((p1x, p2x), (p1x_int, p2x_int), xlabel, ylabel, fname=fname, push=push)

# compare_z_analytical()
# xy_plane_movements()
# compute_errors()
x_phase_plot(True)