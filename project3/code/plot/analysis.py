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



# Analytical solution:
Volt = 9.64852558
m = 40.078
B0 = 10 * Volt
d = 500
V_d2 = 25e-3 * 1e7 * Volt / d**2

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
    




def compare_z_analytical(savepush=False):
    """
    Sanity check.
    """
    FE = load("single_n1.txt", FE_path)
    RK = load("single_n1.txt")

    t = FE[0]
    z_FE = FE[3]
    z_RK = RK[3]
    z0 = z_FE[0]

    z_anal = analytical_z(t, z0)

    plot.z_analytical(z_anal, z_FE, z_RK, t, savepush)
   
    return None 



def compute_errors(method, make_plot=True, savepush=False):
    if method=="FE":
        folder = FE_path 
        title = "Forward Euler"
    elif method=="RK":
        folder = RK_path 
        title = "Runge-Kutta 4"

    abs_errors = []
    rel_errors = []
    times = []
    h_k = []

    for k in range(1,5):
        RK = load(f"single_n{k}.txt", folder)
        t = RK[0]
        r_vec = RK[1:4]
        v_vec = RK[4:]

        x0 = r_vec.T[0,0]
        z0 = r_vec.T[0,-1]
        v0 = v_vec.T[0,1]

        r_anal = np.array([*analytical_xy(t, x0, v0), analytical_z(t, z0)])

        abs_err = np.linalg.norm((r_vec - r_anal), axis=0)
        rel_err = abs_err / np.linalg.norm(r_anal, axis=0)
        
        abs_errors.append(abs_err)
        rel_errors.append(rel_err)
        times.append(t)
        h_k.append(t[-1] / (len(t) - 1))
   
    if make_plot:
        fname = "rel_error_" + method
        plot.error_plot(rel_errors, times, fname, title, savepush)

    else:
        r_err = 0 
        for i in range(1,4):
            r_err += np.log10( np.max(abs_errors[i]) / np.max(abs_errors[i-1]) ) / np.log10(h_k[i] / h_k[i-1])

        print(f"Error convergence rate for {method}: {r_err/3:.5e}")

    return None 



def xy_plane_movements(savepush=False):
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
    plot.two_particles_plane((p1xy, p2xy), (p1xy_int, p2xy_int), xlabel, ylabel, fname=fname, savepush=savepush)

def x_phase_plot(savepush=False):
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
    plot.two_particles_plane((p1x, p2x), (p1x_int, p2x_int), xlabel, ylabel, fname=fname, savepush=savepush)


def z_phase_plot(savepush=False):
    p1_no_int = load("double_without_p1.txt")
    p2_no_int = load("double_without_p2.txt")
    p1_int = load("double_with_p1.txt")
    p2_int = load("double_with_p2.txt")



    p1z = p1_no_int[3::3]
    p2z = p2_no_int[3::3]
    p1z_int = p1_int[3::3]
    p2z_int = p2_int[3::3]

    fname = "z_phase_plot"
    xlabel = r"$z\,[\mathrm{\mu m}]$"
    ylabel = r"$v_z\,[\mathrm{\mu m / \mu s}]$"
    plot.two_particles_plane((p1z, p2z), (p1z_int, p2z_int), xlabel, ylabel, fname=fname, savepush=savepush)


def movement_3d(savepush=False):
    p1_no_int = load("double_without_p1.txt")
    p2_no_int = load("double_without_p2.txt")
    p1_int = load("double_with_p1.txt")
    p2_int = load("double_with_p2.txt")
    
    
    p1 = p1_no_int[1:4]
    p2 = p2_no_int[1:4]

    p1_int = p1_int[1:4]
    p2_int = p2_int[1:4]

    fname = "trajectory_3d"
    plot.two_particles_3d((p1, p2), (p1_int, p2_int), fname=fname, savepush=savepush)



def trapped_without_interaction(savepush=False):
    N0 = 100 
    f_values = [0.1, 0.4, 0.7]
    N_trapped_list = []

    for f in range(1,4):
        omega_V, trapped = load(f"trapped_f{f}_without.txt", skiprows=0)
        N_trapped_list.append(trapped / N0)

    plot.plot_trapped_coarse(N_trapped_list, omega_V, f_values, "trapped_particles_without_interaction", title=None, savepush=savepush)


def trapped_fine_search(savepush=False):
    N0 = 100 
    legs = ["No interactions", "Interactions"]

    f = [0.1, 0.1] 
    n_trapped = []
    omega_V, trapped_int = load(f"trapped_f1_with_fine.txt", skiprows=0)
    omega_V, trapped_noint = load(f"trapped_f1_without_fine.txt", skiprows=0)
    n_trapped.append(trapped_noint / N0)
    n_trapped.append(trapped_int / N0)
    
    plot.plot_trapped_fine(n_trapped, omega_V, legs, "trapped_particles_fine", title=None, savepush=savepush)



# compare_z_analytical(savepush=True)
# compute_errors("RK", savepush=True)
# compute_errors("FE", savepush=True)

# #compute_errors("RK", make_plot=False)
# #compute_errors("FE", make_plot=False)

# # xy_plane_movements(savepush=True)
# # x_phase_plot(savepush=True)
# # z_phase_plot(savepush=True)
# #movement_3d(savepush=True)

# trapped_without_interaction(savepush=True)
# trapped_fine_search(savepush=True)


"""
Uncomment before delivery

sp = False # only show plots

# PROBLEM 8
#   Single-particle case:
compare_z_analytical()
compute_errors("RK", savepush=sp)
compute_errors("FE", savepush=sp)
compute_errors("RK", make_plot=False)
compute_errors("FE", make_plot=False)
#   Double-particle case:
xy_plane_movements(savepush=sp)
x_phase_plot(savepush=sp)
z_phase_plot(savepush=sp)
movement_3d(savepush=sp)

# PROBLEM 9
#   Broad-band scan:
trapped_without_interaction(savepush=sp)
#   Narrow-band scan:
trapped_fine_search(savepush=sp)

"""

print(f"omega_z = {np.sqrt(2/m*V_d2):.2f} MHz")