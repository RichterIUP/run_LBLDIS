#!/usr/bin/python
'''@package docstring
Call LBLRTM and DISORT
'''

# -*- coding: utf8 -*-

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
        
    if scatter:
        sign = 1
    else:
        sign = -1

    with open("lbldis.parm", "w") as file_:
        file_.write("LBLDIS parameter file\n")
        file_.write("16		Number of streams\n")
        file_.write("{:04.1f} 30. 1.0	Solar ".format(sza))
        file_.write("zenith angle (deg), relative azimuth (deg), solar distance (a.u.)\n")
        file_.write(" 180           Zenith angle (degrees): 0 -> ")
        file_.write("upwelling, 180 -> downwelling\n")
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
