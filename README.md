# MIRCA-BC-USMX

Tools for generating planted and irrigated area fractions by bias-correcting MIRCA2000 data to match government records in the US and Mexico

This project contains the scripts used to create irrigation parameters by bias-correcting the MIRCA2000 dataset to match state government records in the US and Mexico [(MIRCA-BC-USMX)](https://zenodo.org/record/2630927) (Bohn and Vivoni, 2019a). These parameters are included along with other land surface parameters in input files formatted for use in the VIC model (Liang et al., 1994) release 5.1 (Hamman et al., 2018) extended to handle irrigation (https://github.com/tbohn/VIC/tree/feature/irrig.imperv.deep_esoil). These parameters are intended for use in VIC's "image" driver.

## This project had several goals:
1. Gather data on planted and irrigated areas of major crop types from the US and Mexico from the same set of years so that the evolution of agricultural land use on both sides or the border in the years since the passage of NAFTA could be assessed.
2. Blend the state-wide totals from step 1 with spatially- and seasonally-explicit data to enable spatial analysis of agricultural land use and its evolution.
3. Provide irrigation parameters for the VIC hydrologic model to enable spatial analysis of agricultural water use and its evolution.

The chosen domain was the continental US and Mexico (the USMX domain), at 1/16 degree (6 km) spatial resolution.

## Processing Stages
The processing stages covered by the scripts in this dataset are:
1. Conversion of USDA records to comma-separated value (csv) format files
2. Conversion of SAGARPA records to comma-separated value (csv) format files
3. Conversion of MIRCA2000 data to ascii-format ESRI grid files and resampling to 1/16 degree resolution
4. Reconciliation MIRCA2000 planted and irrigated areas with NLCD_INEGI 2001/2 cropland areas
5. Bias correction of these reconciled 2001/2 areas to 1990, 2000, and 2010 government records

## Inputs:
### Government Records:
 - USDA National Agricultural Statistics Service Agricultural Census (United States, years 1992, 2002, 2012) (USDA, 2014)
 - SADER and SAGARPA Agricultural Statistics (Mexico, years 1992, 2002, 2012) (SADER, 2014; SAGARPA, 2016)

### NLCD_INEGI Harmonized Land Cover Classification aggregated to 1/16 degree:
 - MOD-LSP VIC parameter files (years 1992, 2001, 2011) (Bohn and Vivoni, 2019b), which contain the NLCD_INEGI land cover fractions

### Map of irrigated areas
 - MIRCA2000 global planted and irrigated areas dataset (Portmann et al., 2010)

### Miscellaneous Files:
 - irr_table.USMX.NLCD_INEGI.csv  - table of irrigation parameters for each land cover class in the NLCD_INEGI legend, available as part of this project.
 - county_codes.csv - table of numeric codes for counties in the US, available for download from [Zenodo](https://zenodo.org/record/2630927) (Bohn and Vivoni, 2019a).
 - mun_us.0.01_deg.asc.tgz and mun_mx.0.01_deg.asc.tgz - county maps for the US and Mexico, available for download from [Zenodo](https://zenodo.org/record/2630927) (Bohn and Vivoni, 2019a).

## Outputs:
 - For each year (1992, 2002, 2012), NetCDF-format files containing 12 monthly maps of (a) planted area fractions and (b) irrigated fractions, at 1/16 degree (6 km) resolution over the USMX domain (continental US + Mexico). These area fractions are reconciled with the NLCD_INEGI maps of years 1992, 2001, and 2011 and are named for those years.
 - For each year (1992, 2002, 2012), input parameter files for the VIC land surface model. Each parameter file contains typical input parameters for VIC plus 12 monthly maps of (a) planted fractions and (b) irrigation fractions. The parameters are relative to the NLCD_INEGI maps of 1992, 2001, and 2011 and are named for those years.
Output files are available for download from [Zenodo](https://zenodo.org/record/2630927) (Bohn and Vivoni, 2019a). 

These parameter files were designed for use with VIC 5 (image driver). VIC 5 image driver requires a "domain" file to accompany the parameter file. This domain file is also necessary for disaggregating the daily gridded meteorological forcings to hourly for input to VIC via the disaggregating tool [MetSim](https://github.com/UW-Hydro/MetSim) (Bennett et al., 2018).  We have provided a domain file compatible with the meteorological forcings of Livneh et al (2015) ("L2015" hereafter) and the MIRCA-BC-USMX parameters, on [Zenodo](https://zenodo.org/record/2564019) (Bohn et al., 2019a,b).


## Directory Structure:

./
 - License.txt - GNU public license
 - README.md - this file

docs/
 - [Procedure.md](./docs/Procedure.md) - list of processing steps, the scripts that run them, and their usage

examples/ - batch files containing example commands and arguments

data/ - miscellaneous input data not obtainable elsewhere

tools/ - scripts for processing the MIRCA2000 data and adding the resulting irrigation parameters to VIC parameter files
 - USDA/ - scripts for processing USDA records
 - SAGARPA/ - scripts for processing SAGARPA records
 - MIRCA2000/ - scripts for processing MIRCA2000 data

## References
 - Bennett, A., J. J. Hamman, B. Nijssen, E. A. Clark, and K. M. Andreadis, 2018: UW-Hydro/MetSim: Version 1.1.0 (version 1.1.0). Zenodo, doi:10.5281/zenodo.1256120. http://doi.org/10.5281/zenodo.1256120 (Accessed June 7, 2018).
 - Bohn, T. J., K. M. Whitney, G. Mascaro, and E. R. Vivoni, 2019a: A deterministic approach for approximating the diurnal cycle of precipitation for use in large-scale hydrological modeling. J. Hydromet., 20(2), 297-317.
 - Bohn, T. J., K. M. Whitney, G. Mascaro, and E. R. Vivoni, 2019b: Parameters for PITRI Precipitation Temporal Disaggregation over continental US, Mexico, and southern Canada, 1981-2013 (Version 1) [Data set]. Zenodo, doi:10.5281/zenodo.1402223. https://zenodo.org/record/1402223.
 - Bohn, T. J, and E. R. Vivoni, 2019a: MIRCA-BC.USMX: Irrigated and planted fractions over the Continental United States and Mexico for years 1992, 2002, and 2012 (Version 1.0) [Data set]. Zenodo, doi:10.5281/zenodo.2630927. https://zenodo.org/record/2630927.
 - Bohn, T. J, and E. R. Vivoni, 2019b: MOD-LSP: MODIS-Based Parameters for Variable Infiltration Capacity (VIC) Model over the Continental US, Mexico, and Southern Canada (Version 1.0) [Data set]. Zenodo, doi:10.5281/zenodo.2612560. https://zenodo.org/record/2612560.
 - Hamman, J. J., B. Nijssen, T. J. Bohn, D. R. Gergel, and Y. Mao, 2018: The Variable Infiltration Capacity Model, Version 5 (VIC-5): Infrastructure improvements for new applications and reproducibility. Geosci. Model Dev., 11, 3481–3496.
 - Liang, X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges, 1994: A simple hydrologically based model of land surface water and energy fluxes for general circulation models. J. Geophys. Res. Atmospheres, 99, 14415–14428.
 - Livneh, B., E. A. Rosenberg, C. Lin, B. Nijssen, V. Mishra, K. M. Andreadis, E. P. Maurer, and D. P. Lettenmaier, 2013: A long-term hydrologically based dataset of land surface fluxes and states for the conterminous United States: Update and extensions. J. Clim., 26, 9384–9392.
 - Portmann, F. T., S. Siebert, and Doll, P., 2010: MIRCA2000--Global monthly irrigated and rainfed crop areas around the year 2000: A new high-resolution data set for agricultural and hydrological modeling. Global Biogeochemical Cycles, 24(1), GB1011, doi:10.1029/2008GB003435.
 - SADER, 2014: Estadistica. Secretaria de Agricultura y Desarrollo Rural (SADER), Mexico, DF. https://sader.gob.mx/siap/estadistica.
 - SAGARPA, 2016: Anuario Estadistico de la Produccion Agricola. Secretaria de Agricultura, Ganaderia, Desarrollo Rural, Pesca y Alimentacion (SAGARPA), Mexico, DF. http://infosiap.siap.gob.mx/aagricola_siap_gb/icultivo/index.jsp.
 - USDA, 2014: 2012 Census of Agriculture: United States Summary and State Data. United States Department of Agriculture (USDA), Washington, DC, https://www.agcensus.usda.gov (Accessed November 1, 2016).
