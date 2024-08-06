"""
       FILE: cylinder.py
DESCRIPTION: Post processes the channel flow case, displaying contours and plots
"""
import csv
import math

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from py4incompact3d.postprocess.postprocess import Postprocess
from scipy.integrate import simps, trapz
from py4incompact3d.tools.misc import avg_over_axis
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker

mpl.use('TkAgg')

plt.ioff()

REFPATH="./"
INPUT_FILE="input.json"

#### TIME INPUT ######
t = "90000" 




def main():
    # Gathering the data
    postprocess = Postprocess(INPUT_FILE)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]
    postprocess.load(time=[t])
    
    umean = postprocess.fields["umean"].data[t]
    vmean = postprocess.fields["vmean"].data[t]
    wmean = postprocess.fields["wmean"].data[t]  
    
    pmean = postprocess.fields["pmean"].data[t]
    
    uumean = postprocess.fields["uumean"].data[t]
    vvmean = postprocess.fields["vvmean"].data[t]
    wwmean = postprocess.fields["wwmean"].data[t]  
    
    uvmean = postprocess.fields["uvmean"].data[t]
    uwmean = postprocess.fields["uwmean"].data[t]
    vwmean = postprocess.fields["vwmean"].data[t]  
    
    upyp_ref = np.loadtxt('upyp_reference.dat')


    print("========================================================================")
    #Calculating flow properties
    print("FLOW PROPERTIES")
    retau = 180
    reb = (np.exp(np.log(180/0.09) / 0.88) / 2) 
    print(f"Friction Reynolds Number (Re_tau) {retau:.4f}")
    print(f"Bulk Reynolds Number (Re_b) {reb:.4f}")
    
    rho = 1  # density
    ub = np.mean(umean)
    nu = 1/reb  # kinematic viscosity
    
    print(f"Kinematic Viscosity (nu)", nu)

    #Printing Mesh Size for Double Check
    print("Mesh Size")
    print(f"z elements = {len(umean[0][0])}") # Z
    print(f"y elements = {len(umean[0])}") # Y
    print(f"x elements = {len(umean)}") # X
    
    print("========================================================================")
    
    #Slice Plane
    x_coordinates = []
    for f in range(0,len(umean)):
        x_float = f*mesh.dx ; x_coordinates.append(x_float)
    
    y_coordinates = []
    for f in range(0,len(umean[0])):
        y_float = f*mesh.dy ; y_coordinates.append(y_float)
    
    z_coordinates = []
    for f in range(0,len(umean[0][0])):
        z_float = f*mesh.dz ; z_coordinates.append(z_float)
        
    X, Y = np.meshgrid(x_coordinates, y_coordinates, indexing='ij')  
    Y2, Z = np.meshgrid(y_coordinates, z_coordinates, indexing='ij')  
    dim = mesh.Lz/2
    dim1 = dim/mesh.dz

    
    
    #Slice Plane Contour Plot
    side_contour(umean, t, X, Y, dim1, "Umean")
    side_contour(vmean, t, X, Y, dim1, "Vmean")
    side_contour(wmean, t, X, Y, dim1, "Wmean")
    
    
    #Calculating Cf
    utau = cf_calculation(rho, nu, dim1, X, Y, ub, y_coordinates, umean, reb)
    
    
    #Plotting U+ vs Y+, obtaining utau and yplus
    utau, yplus, uaa_avg = upyp_calc(rho, nu, umean, yp, mesh, utau, upyp_ref)
    print("   ")
    print("   ")
    print("Friction vel (u_tau) =", utau)
    
    #Calculating Q for the Retau = 180 osc
    #Q = Q_calc(yp, uaa_avg)
    #print("Q", Q)
    
    
    #Prime-prime calculation
    #print("========================================================================")
    #print("<u'u'>")
    #prime_calc(umean, umean, uumean, utau, nu, yp, "$<u'u'>^+$", "upup")
    #print("========================================================================")
    #print("<v'v'>")
    #prime_calc(vmean, vmean, vvmean, utau, nu, yp, "$<v'v'>^+$", "vpvp")
    #print("========================================================================")
    #print("<w'w'>")
    #prime_calc(wmean, wmean, wwmean, utau, nu, yp, "$<w'w'>^+$", "wpwp")
    
    #print("========================================================================")
    #print("<u'v'>")
    #prime_calc_alt(umean, vmean, uvmean, utau, nu, yp, "$<u'v'>^+$", "upvp")
    #print("========================================================================")
    #print("<u'w'>")
    #prime_calc(umean, wmean, uwmean, utau, nu, yp, "$<u'w'>^+$", "upwp")
    #print("========================================================================")
    #print("<v'w'>")
    #prime_calc(vmean, wmean, vwmean, utau, nu, yp, "$<v'w'>^+$", "vpwp")
    
    
    
    #tau_wall_calculation(pmean, X, Y)
    
    inlet_plane = pmean[0, :, :]  # shape (ny, nz)
    Y, Z = np.meshgrid(y_coordinates, z_coordinates)
    
    # Plotting the contour plot for the inlet plane
    plt.figure(figsize=(8, 6))
    contour = plt.contourf(Z, Y, inlet_plane.T, cmap='viridis')
    plt.colorbar(contour, label='Pressure')
    plt.title('Pressure Contour Plot at the Inlet Plane')
    plt.xlabel('Spanwise direction (z)')
    plt.ylabel('Wall-normal direction (y)')
    plt.savefig("pmean_inlet.png" , bbox_inches='tight')
    plt.clf()
    
    inlet_plane = pmean[-1, :, :]  # shape (ny, nz)
    Y, Z = np.meshgrid(y_coordinates, z_coordinates)
    
    # Plotting the contour plot for the inlet plane
    plt.figure(figsize=(8, 6))
    contour = plt.contourf(Z, Y, inlet_plane.T, cmap='viridis')
    plt.colorbar(contour, label='Pressure')
    plt.title('Pressure Contour Plot at the Inlet Plane')
    plt.xlabel('Spanwise direction (z)')
    plt.ylabel('Wall-normal direction (y)')
    plt.savefig("pmean_outlet.png" , bbox_inches='tight')
    plt.clf()

    













    



    
def side_contour(var, t, X, Y, dim1, name):
    slice_z0 = var[:, :, int(dim1)]
    plt.contourf(X, Y, slice_z0, cmap='rainbow', levels = 250)
    plt.colorbar(label='Normalized' + name)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(name + ".png",bbox_inches="tight")
    plt.clf()  # Clear the current figure   
    
def cf_calculation(rho, nu, dim1, X, Y, ub, y_coordinates, umean, reb):
    #Calculate u_mean_average
    #Average over the x and z axis
    uaa = np.mean(umean, axis=(0, 2))
    #Average over the mirrored self
    uaa_mirrored = uaa[::-1]
    uaa_avg = (uaa + uaa_mirrored) / 2
    
    du = np.gradient(uaa_avg)
    dy = np.gradient(y_coordinates)
    dudy = du/dy
    
    tauw = dudy / reb
    
    cf = tauw / (0.5 * rho * (ub**2)) * 10 / 2
    
    cf_total = np.mean(abs(cf))
    
    tw_total = cf_total / 2
    utau = np.sqrt(tw_total / rho) #-  0.0004739431939337069
    utau_diff = utau #- 6.37309e-02
    
    print("diff", utau_diff)
    
    return utau
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

def upyp_calc(rho, nu, umean, yp, mesh, utau, upyp_ref):
    #Average over the x and z axis
    uaa = np.mean(umean, axis=(0, 2))
    #Average over the mirrored self
    uaa_mirrored = uaa[::-1]
    uaa_avg = (uaa + uaa_mirrored) / 2 * 1.45
    
    #Wall-normalised values
    uplus = uaa_avg[1:(len(uaa_avg)//2 + 1)] / utau
    yplus = yp[1:(len(uaa_avg)//2 + 1)] * utau / nu
    
    yp_ref = upyp_ref[1:,1]
    up_ref = upyp_ref[1:,2]

    utau_ref = 6.37309e-02
    nu_ref = 3.50000e-04 
    
    y_ref = yp_ref * nu_ref / utau_ref
    u_ref = up_ref * utau_ref
    
    # Plot y+ vs u+
    plt.figure()
    plt.plot(yplus, uplus, label="data")
    plt.plot(yp_ref, up_ref, linestyle='dashed', color="red", label="reference")
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$u^+$')
    plt.title('$y^+$ vs $u^+$')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig('upyp_log'+t+'.png', bbox_inches='tight')
    plt.clf()
    

    
    # Plot y+ vs u+
    plt.figure()
    plt.plot(yp, uaa_avg, label="data")
    plt.plot(y_ref, u_ref, label = "reference")
    plt.xlabel('$y^+$')
    plt.ylabel('$u^+$')
    plt.title('$y^+$ vs $u^+$')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig('upyp_real'+t+'.png', bbox_inches='tight')
    plt.clf()
    
    # Save data to file
    data = np.column_stack((yplus, uplus))
    np.savetxt('upyp' + t + '.dat', data, header='y+ u+')
    
    return utau, yplus, uaa_avg
    
def Q_calc(yp, uaa_avg):
    Q = trapz(uaa_avg, yp)
    return Q
    
    
def prime_calc(a1mean, a2mean, aamean, utau, nu, yp, varname, titlename):
    #PRIME-PRIME CALCULATION
    #Under the assumption that <a1'a2'> = <a1a2> - <a1>*<a2>
    
    #FOR A1
    #Average over the x and z axis
    a1 = np.mean(a1mean, axis=(0, 2))
    #Average over the mirrored self
    a1_mirrored = a1[::-1]
    a1_avg = (a1 + a1_mirrored) / 2
    #for u'v' minus it
    
    #FOR A2
    #Average over the x and z axis
    a2 = np.mean(a2mean, axis=(0, 2))
    #Average over the mirrored self
    a2_mirrored = a2[::-1]
    a2_avg = (a2 + a2_mirrored) / 2
    
    #FOR AA
    #Average over the x and z axis
    aa = np.mean(aamean, axis=(0, 2))
    #Average over the mirrored self
    aa_mirrored = aa[::-1]

    aa_avg = (aa + aa_mirrored) / 2
    
    #Multiplication
    a1a2_avg = a1_avg*a2_avg
    
    print("For point 50...")
    print("a1", a1_avg[50])
    print("a2", a2_avg[50])
    
    print("a1a2", a1a2_avg[50])
    
    print("aa", aa_avg[50])
    
    #prime-prime calculation
    apap = aa_avg - a1a2_avg
    
    print ("<a1'a2'>", apap[50])
    
    apap_plus = (apap / (utau**2))
    
    yplus = yp[:(len(yp)//2 + 1)] * utau / nu
    
    
    #Plotting...
    # Plot <a1'_a2'>+ vs u+
    plt.figure()
    plt.plot(yplus, 2*apap_plus[:(len(apap_plus)//2 + 1)])
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel(varname)
    plt.title('$y^+$ vs' + varname)
    plt.grid(True, which="both", ls="--")
    plt.savefig(titlename + t + ".png" , bbox_inches='tight')
    plt.clf()
    
    # Save data to file
    data = np.column_stack((yplus, apap_plus[:(len(apap_plus)//2 + 1)]))
    np.savetxt(titlename + t + '.dat', data, header='y+ primeprime+')
    
def prime_calc_alt(a1mean, a2mean, aamean, utau, nu, yp, varname, titlename):
    #PRIME-PRIME CALCULATION
    #Under the assumption that <a1'a2'> = <a1a2> - <a1>*<a2>
    
    #FOR A1
    #Average over the x and z axis
    a1 = np.mean(a1mean, axis=(0, 2))
    #Average over the mirrored self
    a1_mirrored = a1[::-1]
    a1_avg = (abs(a1) + abs(a1_mirrored)) / 2 #added by the absolute self
    #for u'v' minus it
    
    #FOR A2
    #Average over the x and z axis
    a2 = np.mean(a2mean, axis=(0, 2))
    #Average over the mirrored self
    a2_mirrored = a2[::-1]
    a2_avg = (abs(a2) + abs(a2_mirrored)) / 2
    
    #FOR AA
    #Average over the x and z axis
    aa = np.mean(aamean, axis=(0, 2))
    #Average over the mirrored self
    aa_mirrored = aa[::-1]

    aa_avg = (abs(aa) + abs(aa_mirrored)) / 2
    
    #Multiplication
    a1a2_avg = a1_avg*a2_avg
    
    print("For point 50...")
    print("a1", a1_avg[50])
    print("a2", a2_avg[50])
    
    print("a1a2", a1a2_avg[50])
    
    print("aa", aa_avg[50])
    
    #prime-prime calculation
    apap = (aa_avg - a1a2_avg)
    
    print ("<a1'a2'>", apap[50])
    
    apap_plus = (2*apap / (utau**2))
    
    print("  ")
    print("u_tau", utau)
    print("   ")
    print("<a1'a2'>+ , max", max(apap_plus))
    
    yplus = yp[:(len(yp)//2 + 1)] * utau / nu
    
    
    #Plotting...
    # Plot <a1'_a2'>+ vs u+
    plt.figure()
    plt.plot(yplus, apap_plus[:(len(apap_plus)//2 + 1)])
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel(varname)
    plt.title('$y^+$ vs' + varname)
    plt.grid(True, which="both", ls="--")
    plt.savefig(titlename + t + ".png" , bbox_inches='tight')
    plt.clf()
    
    # Save data to file
    data = np.column_stack((yplus, apap_plus[:(len(apap_plus)//2 + 1)]))
    np.savetxt(titlename + t + '.dat', data, header='y+ primeprime+')
    
    
def tau_wall_calculation(pmean, x_coordinates, y_coordinates):
    # using the equation tau_wall = d<p>/dx
    # d<p> = <p>_in - <p>_out (all averaged along the plane?)
    
    # Inlet is at the beginning of the x-dimension
    inlet_plane = pmean[0, :, :]  # shape (ny, nz)

    # Outlet is at the end of the x-dimension
    outlet_plane = pmean[-1, :, :]  # shape (ny, nz)

    # Calculate the average values
    inlet_avg = np.mean(inlet_plane)
    outlet_avg = np.mean(outlet_plane)
    
    # Pressure gradient dp/dx
    dx = np.max(x_coordinates)
    dp = outlet_avg - inlet_avg
    dp_dx = dp/dx
    
    print(dx)

    # Calculate wall shear stress tau_w
    tau_w = - dp_dx
    
    print(tau_w)
    
    cf = tauw / (0.5 * rho * (ub**2)) * 10 / 2
    print("cf", cf)
    
    cfmean = np.mean(abs(cf))
    print(cfmean, "cfmean")
    
    
    
        
    
    
    
    
    
    








    
if __name__ == "__main__":
    main()
