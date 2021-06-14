#!/usr/bin/python
'''@package docstring
Call LBLRTM and DISORT
'''

# -*- coding: utf8 -*-

<<<<<<< HEAD
import subprocess
import sys
import os
import numpy             as np
import netCDF4           as nc
import datetime          as dt
sys.path.append("./run_LBLRTM")

import run_LBLRTM


def run_lbldis(wnum_low, wnum_high, sza, resolution=1., scatter=True, mie_ice=1, h2o_self=1, h2o_foreign=1, co2_cont=0, o3_cont=0, o2_cont=0, n2_cont=0, rayleigh=0):
    
    atmospheric_param = np.loadtxt("./input/atmospheric_param.csv", delimiter=",")
    cloud_grid = np.loadtxt("./input/cloud_grid.csv", delimiter=",")
    z_prof = np.loadtxt("./input/z.csv", delimiter=",")
    t_prof = np.loadtxt("./input/T.csv", delimiter=",")
    q_prof = np.loadtxt("./input/h2o.csv", delimiter=",")
    p_prof = np.loadtxt("./input/P.csv", delimiter=",")
    co2_prof = np.loadtxt("./input/co2.csv", delimiter=",")
    o3_prof = np.loadtxt("./input/o3.csv", delimiter=",")
    co_prof = np.loadtxt("./input/co.csv", delimiter=",")
    ch4_prof = np.loadtxt("./input/ch4.csv", delimiter=",")
    n2o_prof = np.loadtxt("./input/n2o.csv", delimiter=",")
    o2_prof = np.loadtxt("./input/o2.csv", delimiter=",")
    with open("./model_paths", "r") as f:
        PATH_TO_LBLRTM = f.readline().rstrip()
        PATH_TO_LBLDIS = f.readline().rstrip()
    '''
    Set up LBLDIS and run LBLRTM/DISORT
    '''
    now = dt.datetime.now()
    outfolder = "{}_{}_{}_{}_{}_{}".format(now.year, now.month, now.day, now.hour, now.minute, now.second)
    os.mkdir(outfolder)
    print(PATH_TO_LBLRTM)
    # Run LBLRTM
    lbldir  = run_LBLRTM.run_LBLRTM(z = z_prof, \
                                    p = p_prof, \
                                    t = t_prof, \
                                    q = q_prof, \
                                    atm = 4, \
                                    hmd_unit='A', \
                                    wnum1 = wnum_low-50.0, \
                                    wnum2 = wnum_high+50.0, \
                                    lbltp5 = 'tp5', \
                                    lbl_home = PATH_TO_LBLRTM, \
                                    path = ".", \
                                    XSELF=np.int(h2o_self), \
                                    XFRGN=np.int(h2o_foreign), \
                                    XCO2C=np.int(co2_cont), \
                                    XO3CN=np.int(o3_cont), \
                                    XO2CN=np.int(o2_cont), \
                                    XN2CN=np.int(n2_cont), \
                                    XRAYL=np.int(rayleigh), \
                                    co2=co2_prof, \
                                    o3=o3_prof, \
                                    co=co_prof, \
                                    ch4=ch4_prof, \
                                    n2o=n2o_prof, \
                                    o2=o2_prof)
            
    n_layer = len(cloud_grid)
        
=======
import os
import subprocess
import sys
import numpy             as np
import netCDF4           as nc
import scipy.signal      as sig
import pandas            as pd

def conv(wavenumber, radiance, atmospheric_param, resolution):
    '''
    Perfom convolution with sinc
    
    Parameter
    ---------
    wavenumber : np.array
        Wavenumber
        
    radiance : np.array
        Radiance
        
    atmospheric_param : list
        Microphysical Cloud Parameters
        
    resolution : float
        Spectral resolution
    
    Returns
    -------
    np.array
        Convoluted radiance
    '''
    nu_out_center = wavenumber - wavenumber[0]
    nu_out = wavenumber
    opd = 0.9/resolution
    if len(atmospheric_param) <= 1:
        radiance = sig.convolve(radiance.flatten(), opd*np.sinc(opd*nu_out_center), \
                                mode="full")[np.arange(0, len(nu_out), 1)] / \
                                np.sum(opd*np.sinc(opd*nu_out_center))
    else:
        for i, dummy in enumerate(atmospheric_param, start=0):#range(len(inp.MCP)):
            radiance[:, i] = sig.convolve(radiance[:, i], opd*np.sinc(opd*nu_out_center), \
                                           mode="full")[np.arange(0, len(nu_out), 1)] / \
                                           np.sum(opd*np.sinc(opd*nu_out_center))
    radiance = np.array(radiance)
   
    return radiance

####################################################################################

def write_lbldis_input(atmospheric_param, t_surf, ssp, path_wdir, path_windows, sza, cloud_grid, scatter, kurucz, sfc_em, log_re, lbldir):
    '''
    Write input file for LBLDIS
    
    Parameter
    ---------
    atmospheric_param : list
        Microphysical cloud parameters
        
    t_surf : float
        Surface temperature. If negative, surface temperature equals temperature of 
        lowermost atmospheric level
        
    ssp : list
        Names of single-scattering databases
        
    path_wdir : str
        Path to output of TCWret
        
    path_windows : str
        Filename of microwindows
        
    sza : float
        Solar Zenith Angle
        
    cloud_grid : list
        Layers of cloud
        
    scatter : bool
        Use scatter in LBLDIS
        
    kurucz : str
        Name of Kurucz database
        
    sfc_em : list
        Surface emissivity
        
    log_re : bool
        Use logarithmic r_eff
        
    lbldir : str
        Path to lblrtm output
        
    '''
    
>>>>>>> 12fbf7266637445bcaba9e1cb3432d7f69d66ce2
    if scatter:
        sign = 1
    else:
        sign = -1
<<<<<<< HEAD

    with open("lbldis.parm", "w") as file_:
=======
    with open("{}/lbldis.parm".format(path_wdir), "w") as file_:
>>>>>>> 12fbf7266637445bcaba9e1cb3432d7f69d66ce2
        file_.write("LBLDIS parameter file\n")
        file_.write("16		Number of streams\n")
        file_.write("{:04.1f} 30. 1.0	Solar ".format(sza))
        file_.write("zenith angle (deg), relative azimuth (deg), solar distance (a.u.)\n")
        file_.write(" 180           Zenith angle (degrees): 0 -> ")
        file_.write("upwelling, 180 -> downwelling\n")
<<<<<<< HEAD
        file_.write("{} {} {}\n".format(wnum_low, wnum_high, resolution))
        file_.write("{}               Cloud parameter option flag: ".format(np.int(sign)))
        file_.write("0: reff and numdens, >=1:  reff and tau\n".format())
        file_.write("{}".format(2 * n_layer))
        file_.write("               Number of cloud layers\n")
        for loop_liq_layer in range(n_layer):
            alt = cloud_grid[loop_liq_layer]
            ind = np.where((z_prof > alt-1e-3) & (z_prof < alt +1e-3))[0]

            temp_of_layer = t_prof[ind]
            if temp_of_layer < 240 + (253 - 240)/2.0:
                mie_liq = 12
            elif temp_of_layer <= 253 + (263 - 253)/2.0:
                mie_liq = 13
            elif temp_of_layer <= 263 + (273 - 263)/2.0:
                mie_liq = 14
            else:
                mie_liq = 0
            file_.write("{} {:5.3f} {:10.8f} -1".format(mie_liq, alt, atmospheric_param[2]))
            tau_liq = atmospheric_param[0]
            file_.write(" {:10.8f}".format(tau_liq/float(n_layer)))
            file_.write("\n")

        for loop_ice_layer in range(n_layer):
            alt = cloud_grid[loop_ice_layer]*1e-3
            file_.write("{} {:5.3f} {:10.8f} -1".format(mie_ice, alt, atmospheric_param[3]))
            tau_ice = atmospheric_param[1]
            file_.write(" {:10.8f}".format(tau_ice/float(n_layer)))
            file_.write("\n")
        file_.write("{}\n".format(lbldir))
        file_.write("solar/solar.kurucz.rad.1cm-1binned.full_disk.asc\n")
        file_.write("15       Number of scattering property databases\n")
        file_.write("ssp/ssp_db.mie_wat.gamma_sigma_0p100\n")
        file_.write("ssp/ssp_db.mie_ice.gamma_sigma_0p100\n")#1
        file_.write("ssp/ssp_db.Aggregate.gamma.0p100\n")#2
        file_.write("ssp/ssp_db.BulletRosette.gamma.0p100\n")#3
        file_.write("ssp/ssp_db.Droxtal.gamma.0p100\n")#4
        file_.write("ssp/ssp_db.HollowCol.gamma.0p100\n")#5
        file_.write("ssp/ssp_db.Plate.gamma.0p100\n")#6
        file_.write("ssp/ssp_db.SolidCol.gamma.0p100\n")#7
        file_.write("ssp/ssp_db.Spheroid.gamma.0p100\n")#8
        file_.write("ssp/ssp_db.mie_gypsum.lognormal_sigma_0p699\n")#9
        file_.write("ssp/ssp_db.mie_kaolinite.lognormal_sigma_0p699\n")#12
        file_.write("ssp/ssp_db.mie_quartz.lognormal_sigma_0p699\n")#11
        file_.write("ssp/ssp_db.mie_wat_zasetsky240.gamma_sigma_0p100\n")#12
        file_.write("ssp/ssp_db.mie_wat_zasetsky253.gamma_sigma_0p100\n")#13
        file_.write("ssp/ssp_db.mie_wat_zasetsky263.gamma_sigma_0p100\n")#14
        file_.write("-1.	Surface temperature (specifying a negative")
        file_.write("value takes the value from profile)\n")
        file_.write("4	Number of surface spectral emissivity lines (wnum, emis)\n")
        file_.write("100 1\n")
        file_.write("700 1\n")
        file_.write("800 1\n")
        file_.write("3000 1\n")
        
    lbldisout_file = 'lbldisout'
    lbldislog = 'lbldislog.txt'
    with open("run_disort.sh", "w") as file_:
        file_.write("#!/bin/bash\n")
        exec_lbldis = '({}/lbldis lbldis.parm 0 {}) >& {}\n'
        file_.write(exec_lbldis.format(PATH_TO_LBLDIS, \
                                       lbldisout_file, lbldislog))
    subprocess.call(["bash", "run_disort.sh".format(outfolder)])
    
    with nc.Dataset("lbldisout.cdf") as disort_out:
        wavenumber = np.array(disort_out.variables['wnum'][:])
        radiance = np.array(disort_out.variables['radiance'][:])
            
    with nc.Dataset("spectral_radiance.nc", "w") as f:
        f.createDimension("const", 1)
        f.createDimension("mcp", 4)
        f.createDimension("level", len(z_prof))
        f.createDimension("wavenumber", len(wavenumber))
        f.createDimension('cgrid', len(cloud_grid))
        
        mcp = f.createVariable("cloud_parameter", "f8", ("mcp", ))
        mcp.units = "[1, 1, um, um]"
        mcp[:] = atmospheric_param
        
        scat = f.createVariable("scatter", "i8", ("const", ))
        scat.units = "boolean"
        scat[:] = scatter
        
        z_prof_out = f.createVariable("z_profile", "f8", ("level", ))
        z_prof_out.units = "km"
        z_prof_out[:] = z_prof     

        t_prof_out = f.createVariable("T_profile", "f8", ("level", ))
        t_prof_out.units = "K"
        t_prof_out[:] = t_prof
        
        q_prof_out = f.createVariable("h2o_profile", "f8", ("level", ))
        q_prof_out.units = "ppmv"
        q_prof_out[:] = q_prof
        
        co2_prof_out = f.createVariable("co2_profile", "f8", ("level", ))
        co2_prof_out.units = "ppmv"
        co2_prof_out[:] = co2_prof
        
        o3_prof_out = f.createVariable("o3_profile", "f8", ("level", ))
        o3_prof_out.units = "ppmv"
        o3_prof_out[:] = o3_prof
        
        co_prof_out = f.createVariable("co_profile", "f8", ("level", ))
        co_prof_out.units = "ppmv"
        co_prof_out[:] = co_prof
        
        ch4_prof_out = f.createVariable("ch4_profile", "f8", ("level", ))
        ch4_prof_out.units = "ppmv"
        ch4_prof_out[:] = ch4_prof
        
        n2o_prof_out = f.createVariable("n2o_profile", "f8", ("level", ))
        n2o_prof_out.units = "ppmv"
        n2o_prof_out[:] = n2o_prof
        
        o2_prof_out = f.createVariable("o2_profile", "f8", ("level", ))
        o2_prof_out.units = "ppmv"
        o2_prof_out[:] = o2_prof
        
        p_prof_out = f.createVariable("p_profile", "f8", ("level", ))
        p_prof_out.units = "hPa"
        p_prof_out[:] = p_prof
        
        wn_out = f.createVariable("wavenumber", "f8", ("wavenumber", ))
        wn_out.units = "cm-1"
        wn_out[:] = wavenumber
        
        rad_out = f.createVariable("spectral_radiance", "f8", ("wavenumber", ))
        rad_out.units = "mW / (sr cm-1 m2)"
        rad_out[:] = radiance
        
        cgrid = f.createVariable("cloud_profile", "f8", ("cgrid", ))
        cgrid.units = "km"
        cgrid[:] = cloud_grid
        
        sza_out = f.createVariable("solar_zenith_angle", "f8", ("const", ))
        sza_out.units = "deg"
        sza_out[:] = sza
        
        mie_ice_out = f.createVariable('mie_ice', 'f8', ('const', ))
        mie_ice_out.units = "1"
        mie_ice_out[:] = mie_ice
        
        XSELF_out = f.createVariable("XSELF", "i8", ("const", ))
        XSELF_out.units = "bool"
        XSELF_out.description = "H2O self broadened continuum absorption multiplicative factor"
        XSELF_out[:] = np.int(h2o_self)
        
        XFRGN_out = f.createVariable("XFRGN", "i8", ("const", ))
        XFRGN_out.units = "bool"
        XFRGN_out.description = "H2O foreign broadened continuum absorption multiplicative fact"
        XFRGN_out[:] = np.int(h2o_foreign)
        
        XCO2C_out = f.createVariable("XCO2C", "i8", ("const", ))
        XCO2C_out.units = "bool"
        XCO2C_out.description = "CO2 continuum absorption multiplicative factor"
        XCO2C_out[:] = np.int(co2_cont)
        
        XO3CN_out = f.createVariable("XO3CN", "i8", ("const", ))
        XO3CN_out.units = "bool"
        XO3CN_out.description = "O3 continuum absorption multiplicative factor"
        XO3CN_out[:] = np.int(o3_cont)
        
        XO2CN_out = f.createVariable("XO2CN", "i8", ("const", ))
        XO2CN_out.units = "bool"
        XO2CN_out.description = "O2 continuum absorption multiplicative factor"
        XO2CN_out[:] = np.int(o2_cont)
        
        XN2CN_out = f.createVariable("XN2CN", "i8", ("const", ))
        XN2CN_out.units = "bool"
        XN2CN_out.description = "N2 continuum absorption multiplicative factor"
        XN2CN_out[:] = np.int(n2_cont)
        
        XRAYL_out = f.createVariable("XRAYL", "i8", ("const", ))
        XRAYL_out.units = "bool"
        XRAYL_out.description = "Rayleigh extinction multiplicative factor"
        XRAYL_out[:] = np.int(rayleigh)
    subprocess.call(["mv", "lbldis.parm", "{}".format(outfolder)])
    subprocess.call(["mv", "lbldislog.txt", "{}".format(outfolder)])
    subprocess.call(["mv", "run_disort.sh", "{}".format(outfolder)])
    subprocess.call(["mv", "tp5", "{}".format(outfolder)])
    subprocess.call(["mv", "lbldisout.cdf", "{}".format(outfolder)])
    subprocess.call(["mv", "spectral_radiance.nc", "{}".format(outfolder)])


    return wavenumber, radiance

if __name__ == '__main__':
    wnum_low = float(sys.argv[1])
    wnum_high = float(sys.argv[2])
    sza = float(sys.argv[3])
    XSELF = int(sys.argv[4])
    XFRGN = int(sys.argv[5])
    XCO2C = int(sys.argv[6])
    XO3CN = int(sys.argv[7])
    XO2CN = int(sys.argv[8])
    XN2CN = int(sys.argv[9])
    XRAYL = int(sys.argv[10])
    MIE_ICE = int(sys.argv[11])
    SCAT = bool(int(sys.argv[12]))
    wavenumber, radiance = run_lbldis(wnum_low, wnum_high, sza, mie_ice=MIE_ICE, resolution=.1, scatter=SCAT,  h2o_self=XSELF, h2o_foreign=XFRGN, co2_cont=XCO2C, o3_cont=XO3CN, o2_cont=XO2CN, n2_cont=XN2CN, rayleigh=XRAYL)
=======
        file_.write("-1 0 0 {}\n".format(path_windows))
        file_.write("{}               ".format(np.int(len(atmospheric_param)*sign)))
        file_.write("Cloud parameter option flag: ")
        file_.write("0: reff and numdens, >=1:  reff and tau\n")
        file_.write("{}".format(len(ssp) * len(cloud_grid)))
        file_.write("               Number of cloud layers\n")
        for loop_liq_layer, dummy in enumerate(cloud_grid, start=0):
            #ii = 0
            alt = cloud_grid[loop_liq_layer]*1e-3
            
            for i, dummy in enumerate(ssp, start=0):
                tau = atmospheric_param[:, i]
                if log_re:
                    file_.write("{} {:5.3f} {:10.8f} -1".format(i, alt, np.exp(atmospheric_param[0, len(atmospheric_param[0])//2+i])))  
                else:
                    file_.write("{} {:5.3f} {:10.8f} -1".format(i, alt, atmospheric_param[0, len(atmospheric_param[0])//2+i]))
                for tau_lay in tau:
                    file_.write(" {:10.8f}".format(tau_lay/len(cloud_grid)))
                file_.write("\n")
        file_.write("{}\n".format(lbldir))
        file_.write("{}\n".format(kurucz))
        num_db = len(ssp)
        file_.write("{}       Number of scattering property databases\n".format(num_db))
        for database in ssp:
            file_.write(database + "\n")
        file_.write("{}	Surface temperature (specifying a negative".format(t_surf))
        file_.write("value takes the value from profile)\n")
        file_.write("{}	Number of surface spectral emissivity lines (wnum, emis)\n".format(len(sfc_em)))
        for eps in sfc_em:
            file_.write("{} {}\n".format(eps[0], eps[1]))

def run_lbldis(atmospheric_param, lblrtm, ssp, wn, atm_grid, path_to_run_lblrtm, path_to_lblrtm, path_to_lbldis, path_wdir, path_windows, sza, cloud_grid, scatter, kurucz, sfc_em, log_re, lbldir, resolution, t_surf=-1, h2o_self=1, h2o_foreign=1, co2_cont=1, o3_cont=1, o2_cont=1, n2_cont=1, rayleigh=1):
    '''
    Set up LBLDIS and run LBLRTM/DISORT
    
    Parameter
    ---------
    atmospheric_param : list
        Microphysical cloud parameters
        
    lblrtm : bool
        If true, run LBLRTM
        
    ssp : list
        Names of single-scattering databases
        
    wn : list
        Spectral limits of calculation
        
    atm_grid : dict
        Atmospheric profile of pressure, altitude, temperature, humidity and trace gases
        
    path_to_lblrtm : str
        Path to binary of lblrtm
        
    path_to_run_lblrtm : 
        Path to source of run_LBLRTM
        
    path_to_lbldis : str
        Path to binary of lbldis
        
    path_wdir : str
        Path to output of TCWret
        
    path_windows : str
        Filename of microwindows
        
    sza : float
        Solar Zenith Angle
        
    cloud_grid : list
        Layers of cloud
        
    scatter : bool
        Use scatter in LBLDIS
        
    kurucz : str
        Name of Kurucz solar database
        
    sfc_em : list
        Surface emissivity
        
    t_surf : float
        Surface temperature. If negative, surface temperature equals temperature of 
        lowermost atmospheric level
        
    log_re : bool
        Use logarithmic r_eff
        
    lbldir : str
        Path to lblrtm output
        
    resolution : float
        Spectral resolution
        
    Returns
    -------
    np.array
        Radiance
    '''
        
    # Run LBLRTM
    if lblrtm:
        sys.path.append(path_to_run_lblrtm)
        import run_LBLRTM
        lbldir = run_LBLRTM.run_LBLRTM(z = np.array(atm_grid['altitude(km)']), \
                                       p = np.array(atm_grid['pressure(hPa)']), \
                                       t = np.array(atm_grid['temperature(K)']), \
                                       q = np.array(atm_grid['humidity(%)']), \
                                       hmd_unit='H', \
                                       wnum1 = wn[0]-50.0, \
                                       wnum2 = wn[1]+50.0, \
                                       lbltp5 = '{}/tp5'.format(path_wdir), \
                                       lbl_home = path_to_lblrtm, \
                                       path = path_wdir, \
                                       XSELF=np.int(h2o_self), \
                                       XFRGN=np.int(h2o_foreign), \
                                       XCO2C=np.int(co2_cont), \
                                       XO3CN=np.int(o3_cont), \
                                       XO2CN=np.int(o2_cont), \
                                       XN2CN=np.int(n2_cont), \
                                       XRAYL=np.int(rayleigh), \
                                       co2 = np.array(atm_grid['co2(ppmv)']), \
                                       o3 = np.array(atm_grid['o3(ppmv)']), \
                                       co = np.array(atm_grid['co(ppmv)']), \
                                       ch4 = np.array(atm_grid['ch4(ppmv)']), \
                                       n2o = np.array(atm_grid['n2o(ppmv)']), \
                                       o2 = np.array(atm_grid['o2(ppmv)']))

    # Write LBLDIS input file
    write_lbldis_input(atmospheric_param, t_surf, ssp, path_wdir, path_windows, sza, cloud_grid, scatter, kurucz, sfc_em, log_re, lbldir)

    # Run LBLDIS
    lbldisout_file = '{}/lbldisout'.format(path_wdir)
    lbldislog = '{}/lbldislog.txt'.format(path_wdir)
    with open("{}/run_disort.sh".format(path_wdir), "w") as file_:
        file_.write("#!/bin/bash\n")
        exec_lbldis = '({}/lbldis {}/lbldis.parm 0 {}) >& {}\n'
        file_.write(exec_lbldis.format(path_to_lbldis, path_wdir, \
                                       lbldisout_file, lbldislog))
    subprocess.call(["bash", "{}/run_disort.sh".format(path_wdir)])
    
    # Read LBLDIS results and perform convolution
    with nc.Dataset("{}/lbldisout.cdf".format(path_wdir)) as disort_out:
        wavenumber = np.array(disort_out.variables['wnum'][:])
        radiance = conv(wavenumber, np.array(disort_out.variables['radiance'][:]), atmospheric_param, resolution)

    return radiance, lbldir

if __name__ == '__main__':
    with open("input.dat", "r") as f:
        num_ssp = int(f.readline())
        atmospheric_param = []
        ssp = []
        for i in range(num_ssp):
            atmospheric_param.append(float(f.readline()))
            atmospheric_param.append(float(f.readline()))
        for i in range(num_ssp):
            ssp.append(f.readline().rstrip())
            
        wn = [float(f.readline()), float(f.readline())]
        path_to_run_lblrtm = f.readline().rstrip()
        path_to_lblrtm = f.readline().rstrip()
        path_to_lbldis = f.readline().rstrip()
        path_wdir = f.readline().rstrip()
        path_windows = f.readline().rstrip()
        sza = float(f.readline())
        num_cld = int(f.readline())
        cloud_grid = []
        for i in range(num_cld):
            cloud_grid.append(float(f.readline()))
        scatter = bool(int(f.readline()))
        kurucz = f.readline().rstrip()
        num_emis = int(f.readline())
        sfc_em = []
        for i in range(num_emis):
            sfc_em.append([float(f.readline()), float(f.readline())])
        log_re = bool(int(f.readline()))
        lbldir = "."
        resolution = float(f.readline())
        h2o_self = int(f.readline())
        h2o_foreign = int(f.readline())
        co2_cont = int(f.readline())
        o3_cont = int(f.readline())
        o2_cont = int(f.readline())
        n2_cont = int(f.readline())
        rayleigh = int(f.readline())

    atm_grid = pd.read_csv("atm_grid.csv")

    if not os.path.exists(path_wdir): os.mkdir(path_wdir)
    if not os.path.exists(lbldir): os.mkdir(lbldir)

    run_lbldis(np.array([atmospheric_param]), True, ssp, wn, atm_grid, path_to_run_lblrtm, path_to_lblrtm, path_to_lbldis, path_wdir, path_windows, sza, cloud_grid, scatter, kurucz, sfc_em, log_re, lbldir, resolution, t_surf=-1, h2o_self=h2o_self, h2o_foreign=h2o_foreign, co2_cont=co2_cont, o3_cont=o3_cont, o2_cont=o2_cont, n2_cont=n2_cont, rayleigh=rayleigh)
>>>>>>> 12fbf7266637445bcaba9e1cb3432d7f69d66ce2
