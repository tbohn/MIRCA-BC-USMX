# MIRCA-BC-USMX
Tools for generating planted and irrigated area fractions by bias-correcting MIRCA2000 data to match government records in the US and Mexico

This project contains the scripts used to create irrigation parameters by bias-correcting the MIRCA2000 dataset to match state (county?) government records in the US and Mexico [(MIRCA-BC-USMX)](https://zenodo.org/record/xxxxxxx) (ref). These parameters are included along with other land surface parameters in input files formatted for use in the VIC model (Liang et al., 1994) release 5.1 (Hamman et al., 2018) extended to handle irrigation (ref). These parameters are intended for use in VIC's "image" driver.

## This project had several goals:
1. xxx
2. xxx

The chosen domain was the continental US and Mexico (the USMX domain), at 1/16 degree (6 km) spatial resolution.

## Processing Stages
The processing stages covered by the scripts in this dataset are:
1. Conversion of MIRCA2000 data to ascii-format ESRI grid files and resampling to 1/16 degree resolution
2. Aggregation of NLCD_INEGI data from 30 m pixels to class area fractions at 1/16 degree resolution
3. Reconciliation MIRCA2000 planted and irrigated areas with NLCD_INEGI 2001/2 cropland areas
4. Conversion of USDA records to comma-separated value (csv) format files
5. Conversion of SAGARPA records to comma-separated value (csv) format files
6. Bias correction of these reconciled 2001/2 areas to 1990, 2000, and 2010 government records

## Inputs:
### MIRCA2000 global planted and irrigated areas dataset
 - xxx (Portmann et al., 2010)

### NLCD_INEGI harmonized land cover classification:
 - xxx (years ...) (ref)

### Government Records:
 - xxx (United States, years ...) (ref)
 - xxx (Mexico, years ...) (ref)

### VIC input parameters:
 - xxx (domain files) (ref)
 - xxx (parameter files) (ref)

### Miscellaneous Files:
 - xxx

## Outputs:
 - For each year (xxx, xxx, xxx), 24 ascii-format ESRI grid files at 1/16 degree resolution covering the USMX domain, consisting of 12 monthly maps of (a) planted area fractions and (b) irrigated fractions.
 - For each year (xxx, xxx, xxx), input parameter files for the VIC land surface model, available for download from [Zenodo](https://zenodo.org/record/xxxxxxx) (ref). Each parameter file contains typical input parameters for VIC plus 12 monthly maps of (a) planted fractions and (b) irrigation fractions.

These parameter files were designed for use with VIC 5 (image driver). VIC 5 image driver requires a "domain" file to accompany the parameter file. This domain file is also necessary for disaggregating the daily gridded meteorological forcings to hourly for input to VIC via the disaggregating tool [MetSim](https://github.com/UW-Hydro/MetSim) (Bennett et al., 2018).  We have provided a domain file compatible with the meteorological forcings of Livneh et al (2015) ("L2015" hereafter) and the MIRCA-BC-USMX parameters, on [Zenodo](https://zenodo.org/record/2564019).


## Directory Structure:

./
 - License.txt - GNU public license
 - README.md - this file

docs/
 - [Procedure.md](./docs/Procedure.md) - list of processing steps, the scripts that run them, and their usage

examples/ - batch files containing example commands and arguments

data/ - input files other than the land cover classifications or MODIS observations, e.g., masks, tables of class-specific properties
 - USMX/ - data for USMX domain

tools/ - scripts for processing the MIRCA2000 data and adding the resulting irrigation parameters to VIC parameter files

## References
 - Bennett, A., J. J. Hamman, B. Nijssen, E. A. Clark, and K. M. Andreadis, 2018: UW-Hydro/MetSim: Version 1.1.0 (version 1.1.0). Zenodo, doi:10.5281/zenodo.1256120. http://doi.org/10.5281/zenodo.1256120 (Accessed June 7, 2018).
 - Bohn, T. J., and E. R. Vivoni, 2016: Process-based characterization of evapotranspiration sources over the North American monsoon region. Water Resour. Res., 52, 358–384.
 - Bohn, T. J, and E. R. Vivoni, 2019a: MOD-LSP: MODIS-Based Parameters for Variable Infiltration Capacity (VIC) Model over the Continental US, Mexico, and Southern Canada (Version 1.0) [Data set]. Zenodo, doi:10.528/zenodo.2612560. https://zenodo.org/record/2612560.
 - Bohn, T. J, and E. R. Vivoni, 2019b: NLCD_INEGI: Harmonized US-Mexico Land Cover Change Dataset, 1992/2001/2011 (Version 1.1) [Data set]. Zenodo, doi:10.5281/zenodo.2591501. https://zenodo.org/record/2591501.
 - Hamman, J. J., B. Nijssen, T. J. Bohn, D. R. Gergel, and Y. Mao, 2018: The Variable Infiltration Capacity Model, Version 5 (VIC-5): Infrastructure improvements for new applications and reproducibility. Geosci. Model Dev., 11, 3481–3496.
 - Homer, C. G., and Coauthors, 2015: Completion of the 2011 National Land Cover Database for the conterminous United States - Representing a decade of land cover change information. Photogramm. Eng. Remote Sens., 81, 345–354.
 - INEGI, 2014: Conjunto de datos vectoriales de Uso del Suelo y Vegetación, Escala 1:250 000, Serie V (Capa Unión). http://www.inegi.org.mx/geo/contenidos/recnat/usosuelo/ (Accessed October 1, 2015).
 - Liang, X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges, 1994: A simple hydrologically based model of land surface water and energy fluxes for general circulation models. J. Geophys. Res. Atmospheres, 99, 14415–14428.
 - Livneh, B., E. A. Rosenberg, C. Lin, B. Nijssen, V. Mishra, K. M. Andreadis, E. P. Maurer, and D. P. Lettenmaier, 2013: A long-term hydrologically based dataset of land surface fluxes and states for the conterminous United States: Update and extensions. J. Clim., 26, 9384–9392.
