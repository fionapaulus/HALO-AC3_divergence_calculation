# Calculation of divergence and subsidence from dropsondes released during HALO-(AC)$^3$ with Regression Method

Used to process data for Paulus et al. _Airborne measurements of mesoscale divergence at high latitudes during HALO-(AC)³_ with the Regression method.

Author: Fiona Paulus, Institute for Geophysics and Meteorology, University of Cologne, January 2023 (last accessed 18.12.2023) (fpaulus@uni-koeln.de)

Based on George et al. [**JOANNE: Joint dropsonde Observations of the Atmosphere in tropical North atlaNtic meso-scale Environments**](https://essd.copernicus.org/articles/13/5253/2021/essd-13-5253-2021.pdf)


### General
The data processing is performed in 4 different steps (Levels). This processing code starts with Level 1 data, levels 2-4 can be individually performed by turning the `run_Level_x` variable to `True` or `False` in the config file.

The processing can be run by calling 
```
pyhton3 run_processing.py
```

### Level 1: ASPEN Processing
`run_processing.py`starts with Level 1 data, which has been quality controlled using [ASPEN](https://ncar.github.io/aspendocs/pdf/aspendocs.pdf).
The data needs to be stored in directories with a name of the shape `CAMPAIGN_PLATFORM_Dropsondes_DATESTR_FLIGHTSTR` with a folder called `Level_1` with all quality controlled files with a file ending `*QC.nc`. For all further processing steps additional folders `Level_X` will be created in the directory storing the output files from each processing step with file endings `*_LX.nc`

### Level 2: Unit Conversion and Attribute Addition
In Level 2 processing, Non-NAN values of the vertical coordinate are extracted. The vertical coordinate is chosen in `DEFAULT SETTINGS` Section of `config.cfg` and set to `gpsalt` (height above meter) by default. <br />
All variables specified in `DEFAULT SETTINGS` `varname_L2` are included in the Level 2 datasets which are prepared for further processing. If included in `varname_L2`, the units of the following variables will be transformed

| Variable              | Unit (Level 1) | Unit (Level 2) |
|-----------------------|----------------|----------------|
| Relative Humidity     | %              | Fraction       |
| Pressure              | hPa            | Pa             |
| Temperature           | Celsius        | Kelvin         |
| Dew Point Temperature | Celsius        | Kelvin         |

For each Dropsonde a Sonde_ID is added, which correpsonds to the timestamp of the sonde release and global attributes (Platform name, flight time, ...) are added to the dataset. <br />
By default variables `gpsalt [m]`, `lat [°N]`, `lon [°E]`, `time`, `pres [hPa]`, `rh [%]`, `tdry [°C]` and `dp [°C]` are processed.

### Level 3: Interpolation

In Level 3 processing, the files of the quality-controlled individual dropsondes (Level 2) are regridded to match on an evenly spaced grid (height limit and vertical spacing are set in `config.cfg`) and saved as a single output file for each flight. 

Two different interpolation methods are implemented and can be specified in the `GRID SETTINGS` section of `config.cfg` along with the vertical spacing and the maximum height up to which the interpolation is performed. The `linear_interpolate` method uses a direct interpolation algorithm, while the `bin` method first sorts the vertical coordinate data into bins and performs averageing if multiple entries are included in one bin before interpolating. <br />
If the included sondes do not include data between the ground and the defined height limit, they are discarded from the regridding process and this information is stored in a `processing_log.txt` file in the `Level_3` folder.

If not already in the Level 2 ouptut files, specific humidity, saturation specific humidity and water vapor specific humidity are calculated and added to the dataset.

### Level 4: Calculation of Divergence, Vorticity and Advective Tendencies for Dropsonde Patterns



Choose between a carthesian and a spherical coordinate system to calculate divergence and vorticity.

For carthesian coordinates, the distance x and y from 0°N 0°E in meter
```math
x = \text{lon} \cdot 111320\;\text{m} \cdot \cos(\text{lat})
```
```math
y = \text{lat} \cdot 110540\;\text{m}
```
For spherical coordinates, the angles phi and theta are given in radian (computed from lat-lon in degree)
```math 
\theta = \pi \frac{90^\circ - \text{lat}}{180^\circ}
```
```math
\phi = \pi \frac{\text{lon}}{180^\circ}
```

Divergence is defined as
```math
D = \nabla \cdot \vec{v}
```
In carthesian coordinates: D = 	&#8706;<sub>x</sub> u + &#8706;<sub>y</sub> v

In spherical coordinates:  D = 1/(r<sub>Earth</sub> sin(&theta;)) (&#8706;<sub>&phi;</sub> u + &#8706;<sub>&theta;</sub> v)

Vorticity is defined as 
```math
\zeta = \nabla \times \vec{v}
```
In carthesian coordinates:
```math
\zeta = \frac{\mathrm{d}u}{\mathrm{d}y} - \frac{\mathrm{d}v}{\mathrm{d}x}
```
In spherical coordinates: 
```math
D = \frac{1}{R_\text{Earth} \sin(\theta)} \left(\frac{\mathrm{d}v}{\mathrm{d}\phi} - \frac{\mathrm{d}u\sin(\theta)}{\mathrm{d}\theta}\right)
```

## Prerequisites

Phyton Version: Python 3.7.15

Python Packages:
| Package name | Version |
|--------------|---------|
| configparser | 5.0.2 	 |
| numpy        | 1.21.5  |
| pandas       | 1.3.5   |
| xarray       | 0.20.1  |
| tqdm         | 4.64.1  |
| netCDF4      | 1.5.7   |
| metpy        | 1.2.0   |
| geopy        | 2.3.0   |
