# run_LBLDIS.py

run_LBLDIS.py provides an interface to run the radiative transfer models [LBLRTM](http://rtweb.aer.com/lblrtm.html) and [DISORT](http://www.rtatmocn.com/disort/), coupled with [LBLDIS](https://web.archive.org/web/20170508194542/http://www.nssl.noaa.gov/users/dturner/public_html/lbldis/index.html). An older version of DISORT is included in LBLDIS. LBLRTM is used to calculate optical depths of H2O, CH4, CO, CO2, N2O and O3. DISORT calculates the radiative transfer through a plane parallel atmosphere with emitting and scattering objects. All results are in perfect resolution without convolution to a specific instrumental resolution. run_lbldis.py is written in [Python 3](https://www.python.org) and needs the third-party modules NumPy and netCDF4. A python distribution containing a large number of scientific packages is [Anaconda](https://www.anaconda.com/products/individual).

## 1. Download source codes:

- [LBLRTM](http://rtweb.aer.com/lblrtm.html)
- [DISORT (LBLDIS)](https://web.archive.org/web/20170508194542/http://www.nssl.noaa.gov/users/dturner/public_html/lbldis/index.html)

Install python libraries using pip

```sh
> pip install numpy
> pip install netCDF4
```

### Additional files

- [Additional Single scattering databases](https://web.archive.org/web/20170516023452/http://www.nssl.noaa.gov/users/dturner/public_html/lbldis/ADDITIONAL_INFO.html) 

## 2. LBLRTM:

- Extract the source code: 
```sh
> tar -xzvf aer_lblrtm_v12.8.tar.gz
```

LBLRTM requires a spectroscopic database. There are two ways to deal with this:

2.1 Set up the database on your own. To do so, one needs to extract LNFL and the HITRAN database. Next, lnfl needs to be compiled and then the database can be created
```sh
> tar -xzvf lnfl_v3.1.tar.gz
> tar -xzvf aer_v_3.6.tar.gz 
> make -f make_lnfl linuxGNUsgl
> cd ../../aer_v_3.6/line_file
> ../../lnfl/lnfl_v3.1_linux_gnu_sgl
```

2.2 Use the file TAPE3.10-3500cm-1.first_7_molecules

- Create lblrtm/hitran and place the following files:
```sh
> ln -s path/to/TAPE3/file/TAPE3.10-3500.cm-1.first_7_molecules ./tape3.data
> ln -s ../run_examples/xs_files/xs ./xs
> ln -s ../run_examples/xs_files/xs ./x
> ln -s ../run_examples/xs_files/FSCDXS ./FSCDXS
```
- Compile LBLRTM
```sh
> cd lblrtm/build
> make -f make_lblrtm linuxGNUsgl
```
Maybe there could be some errors. Then you should fix the bugs. First bug seems to be in lblatm.f90, line 7967. Insert a whitespace between STOP and the string.
	
- Create the folder lblrtm/bin and link the binary
```sh
> mkdir lblrtm/bin
> ln -s ../lblrtm_v12.8_linux_gnu_sgl ./lblrtm
```

## 3. LBLDIS:

- Create folder lbldis and copy the archive into it. Extract the archive: 
```sh
> mkdir lbldis
> cd lbldis
> tar -xvf lbldis.Release_3_0.tar 
> chmod +w *
```

- By default, LBLDIS wants to be compiled using IFORT. If you want to use GFORTRAN, then change following lines in Makefile:
```sh
- Line 20: Change ifort to gfortran
- Line 26: Remove -nofor_main
- Line 120 and 129: Change -r8 to -fdefault-real-8 -fdefault-double-8
```

- Compile the source
```sh
> make
```

## 3. Input files

- run_lbldis.py reads input files. Atmospheric profiles are in the directory input

### 3.1 model_paths

model_paths contains the paths to LBLRTM and LBLDIS. First line is the path to the root directory of LBLRTM (not to the binary!) and the second line is the path to the binary file of LBLDIS

### 3.2 ./input/atmospheric_param.csv

Contains the cloud optical depth to liquid water, ice water and the effective droplet radii for liquid water and ice water

```sh
1.0,1.0,5.0,30.0
```
### 3.3 ./input/cloud_grid.csv

Contains the height layer of the cloud in km. If this file is empty, clear sky radiances will be calculated

```sh
1.0,2.0,3.0
```

### 3.4 ./input/z.csv

Height layers of the atmosphere in km. Height layers should be defined in a way than the temperature difference between two layers is less then 10 Kelvin!

### 3.5 ./input/P.csv

Pressure layers of the atmosphere in hPa, corresponding to the height layers.

### 3.6 ./input/T.csv

Temperature layers of the atmosphere in K, corresponding to the height layers. Height layers should be defined in a way than the temperature difference between two layers is less then 10 Kelvin!

### 3.7 ./input/h2o.csv

Water Vapour in PPMV.

### 3.8 ./input/ch4.csv

CH4 in PPMV.

### 3.9 ./input/co.csv

CO in PPMV.

### 3.10 ./input/co2.csv

CO2 in PPMV.

### 3.11 ./input/o3.csv

O3 in PPMV

### 3.12 ./input/o2.csv

O2 in PPMV

### 3.13 ./input/n2o.csv

N2O in PPMV

## 4. Run run_lbldis.py

run_lbldis.py needs several command line parameters in following order:
- Lower wavenumber (mininum: 200 cm-1)
- Higher wavenumber (maximum: 3000 cm-1)
- Solar Zenith Angle. If no solar input is wanted, type -1.0
- H2O self broadened continuum absorption multiplicative factor (0 or 1)
- H2O foreign broadened continuum absorption multiplicative factor (0 or 1)
- CO2 continuum absorption multiplicative factor (0 or 1)
- O3 continuum absorption multiplicative factor (0 or 1)
- O2 continuum absorption multiplicative factor (0 or 1)
- N2 continuum absorption multiplicative factor (0 or 1)
- Rayleigh extinction multiplicative factor (0 or 1)
- Ice Particle shape (see table below)
- Scattering in DISORT? (0 or 1)

Spectral radiance interval must be smaller then 2000 cm-1! To run run_lbldis.py, type 

```sh
python3 run_lbldis.py 500.0 1500.0 -1.0 0 0 0 0 0 0 0 1 0
```

This performs a run between 500.0cm-1 and 1500.0cm-1 without solar input and without any continuum for spherical is droplets, without scattering

## 5. Ice particle shapes

| Number | Shape |
|--------|------:|
| 1 | Sphere |
| 2 | Aggregate |
| 3 | Bullet Rosette |
| 4 | Droxtal |
| 5 | Hollow Column |
| 6 | Plate | 
| 7 | Solid Column |
| 8 | Spheroid |