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
from py4incompact3d.tools.misc import avg_over_axis
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker

mpl.use('TkAgg')

plt.ioff()

REFPATH="./"
INPUT_FILE="input.json"

t = "43125" ###TIME MODIFY###




def main():
    
    postprocess = Postprocess(INPUT_FILE)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]
    print(yp)
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


    print("========================================================================")

    retau = 180
    reb = (np.exp(np.log(180/0.09) / 0.88) / 2) 
    print(f"The value of reb is approximately {reb:.4f}")

    rho = 1  # density (kg/m^3)
    ub = np.mean(umean)
    nu = 1/reb  # kinematic viscosity (m^2/s)  

    print("Mesh Size")
    print(f"z elements = {len(umean[0][0])}") # Z
    print(f"y elements = {len(umean[0])}") # Y
    print(f"x elements = {len(umean)}") # X
    
    print("========================================================================")
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
    

        

    du = np.gradient(umean, axis=1)
    dy = np.gradient(y_coordinates)

    dudy = du/dy[:,np.newaxis]
    slice_dudy = dudy[:,:,int(dim1)]
    
    side_contour(umean, t, X, Y, dim1, "Umean")
    side_contour(vmean, t, X, Y, dim1, "Vmean")
    side_contour(wmean, t, X, Y, dim1, "Wmean")
    side_contour(dudy, t, X, Y, dim1, "dudy")
    
    cf_total = cf_calculation(rho, nu, dudy, dim1, X, Y, ub)
    
    utau, yplus = upyp_calc(rho, nu, umean, yp, mesh, cf_total)
    
    print("========================================================================")
    
    print("u_tau", utau)
    
    print("========================================================================")
    print("<u'u'>")
    prime_calc(umean, umean, uumean, utau, nu, yp, "$<u'u'>^+$", "upup")
    print("========================================================================")
    print("<v'v'>")
    prime_calc(vmean, vmean, vvmean, utau, nu, yp, "$<v'v'>^+$", "vpvp")
    print("========================================================================")
    print("<w'w'>")
    prime_calc(wmean, wmean, wwmean, utau, nu, yp, "$<w'w'>^+$", "wpwp")
    
    
    
    
    print("========================================================================")
    print("<u'v'>")
    prime_calc(umean, vmean, uvmean, utau, nu, yp, "$<u'v'>^+$", "upvp")
    print("========================================================================")
    print("<u'w'>")
    prime_calc(umean, wmean, uwmean, utau, nu, yp, "$<u'w'>^+$", "upwp")
    print("========================================================================")
    print("<v'w'>")
    prime_calc(vmean, wmean, vwmean, utau, nu, yp, "$<v'w'>^+$", "vpwp")
    


    













    



    
def side_contour(var, t, X, Y, dim1, name):
    slice_z0 = var[:, :, int(dim1)]
    plt.contourf(X, Y, slice_z0, cmap='rainbow', levels = 250)
    plt.colorbar(label='Normalized' + name)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(name + ".png",bbox_inches="tight")
    plt.clf()  # Clear the current figure   
    
def cf_calculation(rho, nu, dudy, dim1, X, Y, ub):
    tau_x = rho * nu * dudy

    tau_wall = np.mean(tau_x)

    slice_taux = tau_x[:,:,int(dim1)]

    plt.contourf(X, Y, slice_taux, cmap='viridis', levels = 250)
    cbar = plt.colorbar(label='Wall Friction Stress')
    plt.gca().set_aspect('equal', adjustable='box')

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((0, 0))
    cbar.ax.yaxis.set_major_formatter(formatter)

    plt.savefig("tau_wall.png",bbox_inches="tight")
    plt.clf()
    
    tw_bot = tau_x[:,0,:]
    tw_top = tau_x[:, len(tau_x[0])-1, :]
    tw_total = abs(tw_top) + abs(tw_bot)

    mean_tw_top = np.mean(tw_top)
    mean_tw_bot = np.mean(tw_bot)
    mean_tw_total = abs(mean_tw_top) + abs(mean_tw_bot)

    cf_top = mean_tw_top / (0.5 * rho * (ub**2))
    cf_bot = mean_tw_bot / (0.5 * rho * (ub**2))
    cf_total = abs(cf_top) + abs(cf_bot)

    print("Wall Shear Stress")
    print("Top Wall:", mean_tw_top)
    print("Bottom Wall:", mean_tw_bot)
    print("Total:", mean_tw_total)
    print(" ")
    print("Friction Coefficient")
    print("Top Wall:", cf_top)
    print("Bottom Wall:", cf_bot)
    print("Total:", cf_total)
    
    return cf_total

def upyp_calc(rho, nu, umean, yp, mesh, cf_total):
    #Average over the x and z axis
    uaa = np.mean(umean, axis=(0, 2))
    #Average over the mirrored self
    uaa_mirrored = uaa[::-1]
    uaa_avg = (uaa + uaa_mirrored) / 2
    
    #Finding u_tau
    tw_total = cf_total / 2
    utau = np.sqrt(tw_total / rho)
    
    #Wall-normalised values
    uplus = uaa_avg[:(len(uaa_avg)//2 + 1)] / utau
    yplus = yp[:(len(uaa_avg)//2 + 1)] * utau / nu

    
    # Plot y+ vs u+
    plt.figure()
    plt.plot(yplus, uplus)
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$u^+$')
    plt.title('$y^+$ vs $u^+$')
    plt.grid(True, which="both", ls="--")
    plt.savefig('upyp.png', bbox_inches='tight')
    plt.clf()
    
    # Save data to file
    data = np.column_stack((yplus, uplus))
    np.savetxt('upyp' + t + '.dat', data, header='y+ u+')
    
    return utau, yplus
    
def prime_calc(a1mean, a2mean, aamean, utau, nu, yp, varname, titlename):
    #PRIME-PRIME CALCULATION
    #Under the assumption that <a1'a2'> = <a1a2> - <a1>*<a2>
    
    #FOR A1
    #Average over the x and z axis
    a1 = np.mean(a1mean, axis=(0, 2))
    #Average over the mirrored self
    a1_mirrored = a1[::-1]
    a1_avg = (a1 + a1_mirrored) / 2
    
    
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
    
    apap_plus = apap / (utau**2)
    
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
    
    
        
    
    
    
    
    
    








    
if __name__ == "__main__":
    main()
