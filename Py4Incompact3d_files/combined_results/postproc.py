import csv
import math

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.tools.misc import avg_over_axis
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker


mpl.use('TkAgg')

plt.ioff()

# Enable LaTeX for all text
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

####################################################################################
### INPUT FILES

# Reference Values
retau = 180
reb = (np.exp(np.log(180/0.09) / 0.88) / 2) 
rho = 1  # density
nu = 1/reb  # kinematic viscosity
u_tau_ref = 6.37309e-02 #obtained from Lee and Moser
cf_ref = 2*(u_tau_ref)**2
cf_ref_array = np.repeat(cf_ref, 4)

# GCI Values
normres = np.array([0.1982951, 2.5099989, 32])
cf = np.array([0.0082445, 0.0084411, 0.0089732])
normres_ft = 0.0015599452
cf_ft = 0.008124668

#up vs yp filename
upyp_c = np.loadtxt('upyp43125.dat', skiprows=1)
upyp_m = np.loadtxt('upyp52500.dat', skiprows=1)
upyp_f = np.loadtxt('upyp90000.dat', skiprows=1)
upyp_fs = np.loadtxt('upyp240000.dat', skiprows=1)

upyp_ref = np.loadtxt('upyp_reference.dat')


######################################################################################

def main():
    cf_grep = grep(normres,cf)

    print("========================================================================")


    # Insert new values at the beginning
    normres_new = np.insert(normres, 0, normres_ft)
    cf_new = np.insert(cf, 0, cf_ft)

    plot_GREP(normres_new, cf_new, cf_grep, cf_ref_array)

    

    #Calculate u_tau
    u_tau = np.sqrt(cf_new / 2)

    plot_upyp(upyp_c, upyp_m, upyp_f, upyp_fs, upyp_ref, nu, u_tau)

    


    
    
 




    
    







def grep(x, f):
    #Adapted from NASA FORTRAN CODE 
    #https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
    print("GCS AND GREP")
    nd = len(x)

    # Output the number of data sets read
    print("\nNumber of data sets read =", nd)

    print("\n     Grid Size     Quantity\n")
    for n in range(nd):
        print(f"{x[n]:14.6f} {f[n]:14.6f}")

    # Compute the grid refinement ratio, r, between each pair
    r = np.zeros(nd - 1)
    for n in range(nd - 1):
        r[n] = x[n + 1] / x[n]

    # Estimate the order of convergence using the first three data pairs
    p = math.log((f[2] - f[1]) / (f[1] - f[0])) / math.log(r[0])

    print("\nOrder of convergence using first three finest grid")
    print("Order of Convergence, p =", p)

    # Perform Richardson extrapolation to estimate a zero grid value of f
    fexact = f[0] + (f[0] - f[1]) / (r[0] ** p - 1.0)

    if nd >= 3: 
        fsafe = 1.25
    else: 
        fsafe = 3.0

    # Perform Richardson extrapolation to estimate a zero grid value of f
    fexact = f[0] + ( f[0] - f[1] ) / ( r[0]**p - 1.0 )

    print("\nRichardson Extrapolation: Use above order of convergence")
    print("Estimate to zero grid value, f_exact =", fexact)

    print("\nGrid Convergence Index on fine grids. Uses p from above.")
    print("Factor of Safety =", fsafe)
    print("\n  Grid     Refinement")
    print("  Step      Ratio, r       GCI(%)")

    gcif = np.zeros(nd)
    for n in range(nd-1):
        gcif[n] = fsafe * abs(f[n+1] - f[n]) / f[n] / (r[n]**p - 1.0)
        print(f"  {n + 1:2}  {n + 2:3}  {r[n]:14.6f} {gcif[n] * 100.0:14.6f}")

    ratios = np.zeros(nd - 2)
    for n in range(nd - 2):
        if gcif[n + 1] != 0:  # Avoid division by zero
            ratio = r[n] ** p * gcif[n] / gcif[n + 1]
            ratios[n] = ratio
        else:
            ratios[n] = np.nan  # Handle cases where gcif[n + 1] is zero
        print("\nAsymptotic Grid Convergence Ratio", ratio)
    return fexact

def plot_GREP(normres_new, cf_new, cf_grep, cf_ref_array):
    plt.figure()
    plt.rcParams["font.family"] = "serif"
    plt.plot(normres_new, cf_new, marker="o", color = "black", label="GCS Result")
    plt.plot(0, cf_grep, marker='D', linestyle='None', markerfacecolor='none', markeredgecolor='red', markersize=5, label = "GREP")
    plt.plot(normres_new, cf_ref_array, linestyle = "dotted", color = "blue", label="Reference")
    plt.xlabel('Normalised Resolution')
    plt.ylabel('Coefficient of Friction (Cf)')
    plt.legend(frameon=False)
    plt.savefig("GCS"+".png", bbox_inches='tight')
    plt.clf()

def plot_upyp(upyp_c, upyp_m, upyp_f, upyp_fs, upyp_ref, nu, u_tau):
    yp_c = upyp_c[:,0]
    up_c = upyp_c[:,1]/1.45*1.43

    yp_m = upyp_m[:,0]
    up_m = upyp_m[:,1]*1.43

    yp_f = upyp_f[:,0]
    up_f = upyp_f[:,1]*1.43

    yp_fs = upyp_fs[:,0]
    up_fs = upyp_fs[:,1]

    yp_ref = upyp_ref[1:,1]
    up_ref = upyp_ref[1:,2]

    # Plot the data

    # First subplot with original data
    
    plt.subplot(2, 1, 1)
    plt.plot(yp_c, up_c, label='Coarse Grid', color="black")
    plt.plot(yp_m, up_m, label='Medium Grid', color="blue")
    plt.plot(yp_f, up_f, label='Fine Grid', color="green")
    plt.plot(yp_fs, up_fs, label='Finest Grid', color="orange")
    plt.plot(yp_ref, up_ref, label="Reference", color="red", linestyle='dotted')

    # Add labels and title for the first plot
    plt.xlim(1.054557883895367e-02, 1.805584715562184e+02)
    plt.xlabel(r'y$^+$')
    plt.xscale('log')
    plt.ylabel(r'$\langle u \rangle^+$')
    plt.legend(frameon=False)
    plt.text(0.5, 0.9, r'(a)', transform=plt.gca().transAxes, fontsize=10, ha='center')


    # Transform the data
    u_c = up_c * u_tau[3]
    y_c = yp_c * nu / u_tau[3]

    u_m = up_m * u_tau[2]
    y_m = yp_m * nu / u_tau[2]

    u_f = up_f * u_tau[1]
    y_f = yp_f * nu / u_tau[1]

    u_fs = up_fs * u_tau[0]
    y_fs = yp_fs * nu / u_tau[0]

    u_ref = up_ref * u_tau_ref
    y_ref = yp_ref * nu / u_tau_ref

    # Second subplot with transformed data
    plt.subplot(2, 1, 2)
    plt.plot(y_c, u_c, label='Coarse Grid', color="black")
    plt.plot(y_m, u_m, label='Medium Grid', color="blue")
    plt.plot(y_f, u_f, label='Fine Grid', color="green")
    plt.plot(y_fs, u_fs, label='Finest Grid', color="orange")
    plt.plot(y_ref, u_ref, label="Reference", color="red", linestyle='dotted')

    # Add labels and title for the second plot
    plt.xlim(0, 1)
    plt.xlabel(r'y')
    plt.ylabel(r'$\langle u \rangle$')
    plt.legend(frameon=False)
    plt.text(0.5, 0.9, r'(b)', transform=plt.gca().transAxes, fontsize=10, ha='center')

    # Show the plot
    plt.tight_layout()
    plt.savefig('upyp.png', bbox_inches='tight')
    plt.clf()
    ######################################################################################
    

if __name__ == "__main__":
    main()