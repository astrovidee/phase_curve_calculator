"""
Created on Thur Mar 21 18:25:32 2024
@author: Vidya Venkatesan (vidyav1@uci.edu)

Project description: This code is a python version of thermal_phase_curve.pro which calculates phase curve fluxes
based on Koll and Abbot et al., 2014 (https://iopscience.iop.org/article/10.1088/0004-637X/802/1/21)

Required Inputs:
    OLR file for thermal light curve
    RSW file for reflected light curve
Output:
Hopefully a thermal curve 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

############################
##### Load Input Files #####
############################

def load_inputs():
    """
    Load input data files for latitude, OLR, and time.
    
    Returns:
        lat (numpy array): array of latitude values [degrees]
        olr_data (numpy array): OLR data (either 1D or 2D depending on case)
        time (numpy array): time array (used for e>0 case)
    """
    lat = np.loadtxt("phi.txt")
    #olr_data = np.loadtxt("OLR_e0_f50_cs_c2k.txt")
    olr_data=np.loadtxt("OLR_e5_f50_cs_c2k.txt")
    try:
        time = np.loadtxt("n_days.txt")
    except:
        time = np.array([0])  # for e=0 case if time file missing
    return lat, olr_data, time

########################################
##### e=0 Phase Curve Calculation ######
########################################

def compute_phase_curve_e0(lat, olr_data, time, REARTH):
    """
    Compute thermal phase curve for circular orbit (e=0).
    
    Assumes no longitudinal variation in OLR; flux depends only on latitude.
    
    Inputs:
        lat (numpy array): latitude array
        olr_data (numpy array): 1D array of OLR(lat)
        time (numpy array): time array (used for output length)
        REARTH (float): planet radius [meters]

    Returns:
        phase_curve (numpy array): array of disk-integrated flux values (flat)
        tflux (numpy array): total flux in W
    """

    ##### Activity ######
    """
    For each latitude, integrate outgoing longwave flux over the disk using:
        integral = sum(OLR(lat) * cos^2(lat) * cos(view angle))
    Since e=0, no viewing geometry rotation is applied.
    """
    rad = np.pi / 180.
    dtheta = 4 * rad
    dphi = 5 * rad

    # If user accidentally provides 2D OLR for e=0, take first time slice
    if olr_data.ndim == 2:
        olr_data = olr_data[0, :]

    F = 0.0
    D = 0.0
    for x in range(len(lat)):
        LWVAR_temp = olr_data[x] * np.cos(rad * lat[x])**2 * np.cos(rad) * dphi * dtheta
        DENOM_temp = np.cos(rad * lat[x])**2 * np.cos(rad) * dphi * dtheta
        F += LWVAR_temp
        D += DENOM_temp

    avg_flux = F / D
    phase_curve = np.full(len(time)+1, avg_flux)  # flat phase curve
    phase_curve[-1] = phase_curve[0]
    tflux = np.pi * REARTH**2 * phase_curve
    return phase_curve, tflux

########################################
##### e>0 Phase Curve Calculation ######
########################################

def compute_phase_curve_eccentric(lat, olr_data, time, REARTH):
    """
    Compute thermal phase curve for eccentric orbit (e>0).
    
    Rotates observer longitude at each time step and integrates visible hemisphere.

    Inputs:
        lat (numpy array): latitude array
        olr_data (numpy array): 2D array of OLR(time, lat)
        time (numpy array): time array (used for phase angle)
        REARTH (float): planet radius [meters]

    Returns:
        phase_curve (numpy array): array of flux values as function of phase
        tflux (numpy array): total flux in W
    """

    ##### Activity ######
    """
    For each orbital phase (ref_lon), rotate the sub-observer longitude and integrate
    the visible hemisphere using viewing geometry. This includes:
        - visible longitude sector (within +/- 90 deg of sub-observer longitude)
        - limb darkening effects with cos(lat)^2 * cos(tlon)
    """
    rad = np.pi / 180.
    dtheta = 4 * rad
    dphi = 5 * rad

    lon = time
    phase_curve = np.zeros(len(lon)+1)
    F = np.zeros(len(lon)+1)
    D = np.zeros(len(lon)+1)

    for i in range(len(lon)):
        ref_lon = lon[i]
        xview = []
        if ref_lon < 90:
            xview = np.where((lon <= ref_lon + 90.) | (lon >= 360. + (ref_lon - 90.)))[0]
        elif 90 <= ref_lon < 270:
            xview = np.where((lon <= ref_lon + 90.) & (lon >= ref_lon - 90.))[0]
        elif ref_lon >= 270:
            xview = np.where((lon <= ref_lon - 270.) | (lon >= ref_lon - 90.))[0]

        for y in range(len(lat)):
            for x in range(len(xview)):
                tlon = ref_lon - lon[xview[x]]
                LWVAR_temp = olr_data[xview[x], y] * np.cos(rad * lat[y])**2 * np.cos(rad * tlon) * dphi * dtheta
                DENOM_temp = np.cos(rad * lat[y])**2 * np.cos(rad * tlon) * dphi * dtheta
                F[i] += LWVAR_temp
                D[i] += DENOM_temp

        phase_curve[i] = F[i] / D[i]

    phase_curve[-1] = phase_curve[0]
    tflux = np.pi * REARTH**2 * phase_curve
    return phase_curve, tflux

########################################
##### Main driver function #############
########################################

def main():
    """
    Main code driver:
    - Load inputs
    - Ask user for e=0 or e>0 case
    - Compute phase curve accordingly
    - Plot and save results
    """
    REARTH = 6.37122e6 * 0.93  # planet radius in meters

    # Load data
    lat, olr_data, time = load_inputs()

    # Ask user to specify eccentricity case
    ecc_case = input("Is this an e=0 case? (y/n): ").strip().lower()

    if ecc_case == 'y':
        print("Running e=0 phase curve calculation...")
        phase_curve, tflux = compute_phase_curve_e0(lat, olr_data, time, REARTH)
    else:
        print("Running e>0 phase curve calculation...")
        phase_curve, tflux = compute_phase_curve_eccentric(lat, olr_data, time, REARTH)

    # Save output to text file
    np.savetxt("phase_curve_output.txt", phase_curve, fmt='%s')

    # Plot
    plt.figure()
    mpl.rc('font', family='Serif')
    plt.plot(phase_curve)
    plt.xlabel('Orbital phase')
    plt.ylabel('Flux [W/m²]')
    plt.title('Thermal Phase Curve')
    plt.show()

if __name__ == "__main__":
    main()
