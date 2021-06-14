# run\_LBLDIS.py

run\_LBLDIS.py creates input files for LBLRTM and LBLDIS and executes both models. It requires run\_LBLRTM.py, LBLRTM (tested for version 12.8) and LBLDIS (tested for version 3.0)

## Software requirements

run\_LBLDIS.py is written in Python and tested with Python 3.8.5. Necessary libraries: 

- Numpy (creating arrays, calculating sinc)
- Scipy (perform convolution)
- netCDF4 (read output of LBLDIS)
- Pandas (read atmospheric profile)

Python and several scientific libraries are assembled in [Anaconda](https://www.anaconda.com/).

LBLDIS and LBLRTM must be already compiled on your system to use run\_LBLDIS.py

## Usage of run\_LBLDIS.py

```sh
python3 run_lbldis.py
```

run\_LBLDIS.py searches for input.dat and atm\_grid.csv in your run directory. 

## input.dat

input.dat is an ASCII file containing the parameters and paths which are necessary to use run\_LBLDIS.py. It is explained using the example input.dat.example

- Line  1: Must be integer and gives the number of different single scattering databases. For each single scattering database, input.dat must contain two parameters: Optical depth and effective radius
- Line  2: Optical depth corresponding to first single scattering database
- Line  3: Optical depth corresponding to second single scattering database
- Line  4: Effective radius corresponding to first single scattering database
- Line  5: Effective radius corresponding to second single scattering database
- Line  6: First single scattering database. Must be ASCII and readable by LBLDIS
- Line  7: Second single scattering database. Must be ASCII and readable by LBLDIS
- Line  8: Lower limit of spectral interval. Must be at least 50cm-1 below lowest wavenumber in microwindow file 
- Line  9: Upper limit of spectral interval. Must be at least 50cm-1 above highest wavenumber in microwindow file
- Line 10: Path to run\_LBLRTM.py
- Line 11: Path to LBLRTM (root directory, binary must be in PATH\_TO\_LBLRTM/bin/)
- Line 12: Path to LBLDIS (binary)
- Line 13: Output directory
- Line 14: Name of file containing microwindows. Must be readable by LBLDIS
- Line 15: Solar Zenith Angle in degrees. Negative Value deactives solar input
- Line 16: Number of cloud levels
- Line 17 - 21: Cloud levels in meter
- Line 22: Activate (1) or deactivate (0) scattering in LBLDIS
- Line 23: Path to Kurucz solar input. Must be readable by LBLDIS
- Line 24: Number of wavenumbers where surface emissivity is defined. LBLDIS interpolates emissivity to radiation grid
- Line 25: First wavenumber of surface emissivity
- Line 26: First surface emissivity
- Line 27: Second wavenumber of surface emissivity
- Line 28: Second surface emissivity
- Line 29: Third wavenumber of surface emissivity
- Line 30: Third surface emissivity
- Line 31: Effective radii are given in absolute values (0) or as logarithms (1)
- Line 32: Spectral resolution, needed for convolution
- Line 33: Activate (1) or deactivate (0) H2O self continuum
- Line 34: Activate (1) or deactivate (0) H2O foreign continuum
- Line 35: Activate (1) or deactivate (0) CO2 continuum
- Line 36: Activate (1) or deactivate (0) O3 continuum
- Line 37: Activate (1) or deactivate (0) O2 continuum
- Line 38: Activate (1) or deactivate (0) N2 continuum
- Line 39: Activate (1) or deactivate (0) Rayleigh continuum

Line numbers change if differnt numbers of single scattering databases and surface emissivity lines are specified.

## atm\_grid.csv

This file is read using Pandas, therefore the order of columns does not matter. The file must contain following columns (case sensitive):

- altitude(km)
- ch4(ppmv)
- co(ppmv)
- co2(ppmv)
- humidity(%)
- n2o(ppmv)
- o2(ppmv)
- o3(ppmv)
- pressure(hPa)
- temperature(K)
>>>>>>> 12fbf7266637445bcaba9e1cb3432d7f69d66ce2
