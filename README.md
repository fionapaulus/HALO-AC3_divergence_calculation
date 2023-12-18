# Code publication to Paulus et al. - Airborne measurements of mesoscale divergence at high latitudes during HALO-(AC)$^3$

Calculation of divergence and subsidence from dropsondes released during HALO-(AC)$^3$ campaign in circular and rectangular flight patterns with two methods: 

* Line-Integral method
* Regression method

The dropsondes data from the campaign will be published soon and can be accessed already through the [ac3airbrone package](https://github.com/igmk/ac3airborne/). The present code is based on preprocessed (Level 1) data, which has been quality controlled using [ASPEN](https://ncar.github.io/aspendocs/pdf/aspendocs.pdf). The data needs to be stored in directories with a name of the shape `CAMPAIGN_PLATFORM_Dropsondes_DATESTR_FLIGHTSTR` with a folder called `Level_1` with all quality controlled files with a file ending `*QC.nc`. For all further processing steps additional folders `Level_X` will be created in the directory storing the output files from each processing step with file endings `*_LX.nc`.