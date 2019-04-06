# Processing Steps for the MIRCA-BC-USMX Project

To use these scripts, make sure that the paths to the "tools/USDA", "tools/SAGARPA", and "tools/MIRCA" directories are in your `$PATH` environment variable. These scripts require Python 3.6.4 and Perl 5 to be installed on your system.

## Overview of Steps:

1. Install the `AscGridTools` project
2. Download the appropriate input datasets
3. Process USDA files
4. Process SAGARPA files
5. Process MIRCA files

## Step 1. Install the `AscGridTools` project

The scripts in this project depend on tools from the [`AscGridTools` Git repo](https://github.com/tbohn/AscGridTools). Obtain `AscGridTools` by either forking or cloning the repo, and follow the instructions to compile the tools. Add the path to your copy `AscGridTools/tools` to your $PATH environment variable.

## Step 2. Download the appropriate input datasets

You will need the following datasets:
 - USDA Agricultural Census records at the county level for the years 1992, 2002, and 2012. You can obtain them [here](https://www.agcensus.usda.gov). These will be individual pdf files, one per state per year. You should create a top-level directory somewhere (e.g., `my_path/USDA`) that we will refer to as `$USDA_ROOT`. Under `$USDA_ROOT`, create a subdirectory called `orig` and store the pdf files there. The pdf files should be named with the convention `table.$CODE.$YEAR.pdf`, where `$CODE` is the two-digit numerical code for the state (e.g., Alaska's code is 01 and Wyoming's code is 50) and `$YEAR` is the survey year.
 - Maps of counties/municipios  in US and Mexico (mun_us.0.01_deg.asc.tgz and mun_mx.0.01_deg.asc.tgz) along with table linking county codes to county names used in US Census and USDA files (county_codes.csv). You can obtain them from the same [Zenodo location where the outputs of this project are stored](https://zenodo.org/record/2630927). Download these and place them in the `$PROJECT/data` directory.
 - SAGARPA Agricultural records at the county (municipio) level for the years 1980-2012 (from which we will take just 1992, 2002, and 2012). You can obtain them [here](https://sader.gob.mx/siap/estadistica). The files you need are named `180822-agt-cierre-$YEAR1-$YEAR2.csv`, where `$YEAR1` and `$YEAR2` are beginning and ending years for the file. You should create a top-level directory somewhere (e.g., `my_path/SAGARPA`) that we will refer to as `$SAGARPA_ROOT`. Under `$SAGARPA_ROOT`, create a subdirectory called `orig` and store the csv files there. You need to concatenate these into a single file named Agt_cierre_1980-2012.csv, also stored in the `orig` directory.
 - MIRCA2000 dataset, specifically the ascii ESRI grid file versions of the dataset. You can obtain them [here](https://www.uni-frankfurt.de/45218031/data_download). You should create a top-level directory somewhere (e.g., `my_path/MIRCA2000`) that we will refer to as `$MIRCA_ROOT`. Under `$MIRCA_ROOT`, create a subdirectory called `asc.orig` and store the ascii grid files there. The files should be named with the convention `crop_$CODE_$TYPE_$MONTH.asc`, where `$CODE` is the two-digit crop code (from 1 to 26), `$TYPE` is one of `irrigated` or `rainfed`, and `$MONTH` is the 3-digit month number (from 001 to 012).
 - MOD-LSP VIC parameter files based on the NLCD_INEGI land cover classification, at 1/16 degree resolution. You can obtain them [here](https://zenodo.org/record/2612560). Specifically, you'll need the files `params.USMX.NLCD_INEGI.s1992.2000_2016.nc`, `params.USMX.NLCD_INEGI.2001.2000_2016.nc`, and `params.USMX.NLCD_INEGI.2011.2000_2016.nc`. You should create a top-level directory somewhere (e.g., `my_path/MOD-LSP/params`) that we will refer to as `$VIC_PARAM_ROOT`. Store these parameter files in that directory.

## Step 3. Process the USDA files

This step has several sub-steps:
1. Reformat the USDA data from pdf to csv files.

   `wrap_parse_usda_pdfs.pl $USDA_ROOT/orig table $YEAR $PROJECT/data/county_codes.csv $USDA_ROOT/reformat`

   where

   `$USDA_ROOT` = top-level directory where USDA data are stored
   `$YEAR` = the desired census year (1992, 2001, or 2011)
   `$PROJECT` = path to your copy of this Git repo

   Note that this step uses a table of US Census and USDA county codes (county_codes.csv) that was compiled by a combination of perl commands and manual correction of the entries.

2. Combine the data from different states into a single file per year; also sum the irrigated and planted areas to the state level.

   `collate_usda_county_data.pl $USDA_ROOT/reformat $YEAR $USDA_ROOT/summ`

   where

   `$USDA_ROOT` = top-level directory where USDA data are stored
   `$YEAR` = the desired census year (1992, 2001, or 2011)

### Example

An example of how to run these scripts is in `$PROJECT/examples/batch.process_usda.csh`.

To run this batch file, replace all instances of `$PROJECT` in the file with the path to your copy of this GitHub project (and replace any other placeholders beginning with `$` with the appropriate values). Also make sure it is executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file.

## Step 4. Process the SAGARPA files

This step has several sub-steps:
1. Clean up the csv files, which contain some instances of fields containing commas surrounded by quotation marks. Because commas are also the field delimiter, these commas inside quotes are confusing for scripts to process.

   `remove_nondelim_commas_in_csv.pl $SAGARPA_ROOT/orig/Agt_cierre_1980_2012.csv > $SAGARPA_ROOT/no_extra_commas/Agt_cierre_1980_2012.csv`

   where

   `$SAGARPA_ROOT` = top-level directory where SAGARPA data are stored

2. Compute total areas of major irrigated and rainfed crop groups.

   `process_sagarpa_files.pl $SAGARPA_ROOT/no_extra_commas/Agt_cierre_1980_2012.csv > $SAGARPA_ROOT/summ_by_mun/summ_by_mun.csv`

   where

   `$SAGARPA_ROOT` = top-level directory where SAGARPA data are stored

### Example

An example of how to run these scripts is in `$PROJECT/examples/batch.process_sagarpa.csh`.

To run this batch file, replace all instances of `$PROJECT` in the file with the path to your copy of this GitHub project (and replace any other placeholders beginning with `$` with the appropriate values). Also make sure it is executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file.

## Step 5. Process the MIRCA2000 files

This step has several sub-steps:

1. Under `$PROJECT/data/`, extract the contents of `mun_us.0.01_deg.asc.tgz` and `mun_mx.0.01_deg.asc.tgz`. To do this, cd to `$PROJECT/data` and run the commands:

   `tar -xvzf mun_us.0.01_deg.asc.tgz`  
   `tar -xvzf mun_mx.0.01_deg.asc.tgz`  

2. Compute total areas for irrigated, rainfed, and all crops

   `add_crop_areas.pl $MIRCA_ROOT/asc.orig $MIRCA_ROOT/asc.crop_groups` 

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

3. Convert from areas to area fractions

   `convert_areas_to_fractions_of_cell_area.pl $MIRCA_ROOT/asc.crop_groups $MIRCA_ROOT/asc.frac`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

4. Clip to domain and resample to 1/16 degree resolution

   `wrap_gridclip.asc.pl $MIRCA_ROOT/asc.frac crop float 14 50 -125 -65 16 4 $MIRCA_ROOT/USMX/asc.0.08333333_deg.frac_of_cell`
   `wrap_grid_subsample.pl $MIRCA_ROOT/USMX/asc.0.08333333_deg.frac_of_cell crop float 0.0208333333333333 nn 0 float 16 4 $MIRCA_ROOT/USMX/asc.0.02083333_deg.frac_of_cell`
   `wrap_grid_agg.pl $MIRCA_ROOT/USMX/asc.0.02083333_deg.frac_of_cell crop float 3 5 4 $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell.tmp`
   `wrap_gridclip.asc.pl $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell.tmp crop float 14 50 -125 -65 5 4 $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell`
   `rm -rf $MIRCA_ROOT/USMX/asc.0.02083333_deg.frac_of_cell $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell.tmp`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

5. Compute planted fraction (fplant) and irrigated fraction (firr)

   `wrap_compute_firr_of_fplant_from_mirca.pl $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell $MIRCA_ROOT/USMX/asc.0.0625_deg.fplant_and_firr` 

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

6. Convert to netcdf

   `mkdir -p $MIRCA_ROOT/USMX/nc`
   `convert_monthly_single_var_asc_to_nc.py -i $MIRCA_ROOT/USMX/asc.0.0625_deg.fplant_and_firr -p firr_of_fplant_0 -v firr_of_fplant -o $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc`
   `convert_monthly_single_var_asc_to_nc.py -i $MIRCA_ROOT/USMX/asc.0.0625_deg.fplant_and_firr -p fplant_of_cell_0 -v fplant_of_cell -o $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

7. Null out low values

   `null_out_mirca_where_fplant_low.py -f $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -t 0.01 -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc -v firr_of_fplant -o tmp.nc`
   `mv tmp.nc $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc`
   `null_out_mirca_where_fplant_low.py -f $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -t 0.01 -i $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -v fplant_of_cell -o tmp.nc`
   `mv tmp.nc $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc` 

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

8. Compute the planted fraction relative to NLCD_INEGI 2001 cropland areas

   `compute_fplant_of_crop_mirca_lc.py -f $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -p $VIC_PARAM_ROOT/params.USMX.2001.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc`
   `compute_annual_mirca_fplant.py -f $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc -v fplant_of_cell -o $MIRCA_ROOT/USMX/nc/fplant_annual.nc`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored

9. Fill gaps in MIRCA planted areas (where NLCD cropland area of given year > 0)

   `gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc -v fplant_of_crop,fplant_of_cell -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/$VIC_PARAM_FILE -c 12,13 -o $MIRCA_ROOT/USMX/nc/fplant.$YEAR.nc`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored
   `$VIC_PARAM_ROOT` = the directory where VIC parameter files are stored
   `$VIC_PARAM_FILE` = the VIC parameter file corresponding to the desired year
   `$YEAR` = the desired year

10. Fill gaps in MIRCA irrigated areas (where MIRCA gap-filled planted areas > 0)

   `gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc -v firr_of_fplant -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/$VIC_PARAM_FILE -c 12,13 -o $MIRCA_ROOT/USMX/nc/firr_of_fplant.$YEAR.nc`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored
   `$VIC_PARAM_ROOT` = the directory where VIC parameter files are stored
   `$VIC_PARAM_FILE` = the VIC parameter file corresponding to the desired year
   `$YEAR` = the desired year

11. Convert irrigated fractions of planted areas to irrigated fractions of (a) NLCD cropland area and (b) cell area

   `compute_firr_of_crop_and_cell.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.$YEAR.nc -f $MIRCA_ROOT/USMX/nc/fplant.$YEAR.nc -v fplant_of_crop,fplant_of_cell -o $MIRCA_ROOT/USMX/nc/firr.$YEAR.nc -n firr_of_crop,firr_of_cell`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored
   `$VIC_PARAM_ROOT` = the directory where VIC parameter files are stored
   `$VIC_PARAM_FILE` = the VIC parameter file corresponding to the desired year
   `$YEAR` = the desired year

12. Bias correct the planted and irrigated fractions so that their total areas sum to the state-by-state values in government records

   `bias_correct_mirca_fplant_firr.py -v $VIC_PARAM_ROOT/$VIC_PARAM_FILE -c 12,13 -p $MIRCA_ROOT/USMX/nc/fplant.$LCYEAR.nc -i $MIRCA_ROOT/USMX/nc/firr.$LCYEAR.nc -u $USDA_ROOT/summ/summ.$YEAR.csv -s $SAGARPA_ROOT/summ_by_mun.csv -m $PROJECT/data/mun_us.0.01_deg.asc -x $PROJECT/data/mun_mx.0.01_deg.asc -y $YEAR -o $MIRCA_ROOT/USMX/nc/fplant_firr_bc.$LCYEAR.nc`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored
   `$USDA_ROOT` = the top-level directory under which USDA data are stored
   `$SAGARPA_ROOT` = the top-level directory under which SAGARPA data are stored
   `$VIC_PARAM_ROOT` = the directory where VIC parameter files are stored
   `$VIC_PARAM_FILE` = the VIC parameter file corresponding to the desired year
   `$YEAR` = the desired year
   `$LCYEAR` = the desired land cover year (for $YEAR=(1992,2002,2011), $LCYEAR=(s1992,2001,2011))

13. Special Case: Bias correct the planted and irrigated fractions state-by-state in the US and county-by-county in Mexico

   `bias_correct_mirca_fplant_firr.py -v $VIC_PARAM_ROOT/$VIC_PARAM_FILE -c 12,13 -p $MIRCA_ROOT/USMX/nc/fplant.$LCYEAR.nc -i $MIRCA_ROOT/USMX/nc/firr.$LCYEAR.nc -u $USDA_ROOT/summ/summ.$YEAR.csv -s $SAGARPA_ROOT/summ_by_mun.csv -m $PROJECT/data/mun_us.0.01_deg.asc -x $PROJECT/data/mun_mx.0.01_deg.asc -y $YEAR -d -o $MIRCA_ROOT/USMX/nc/fplant_firr_bc.$LCYEAR.mun_mx.nc`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored
   `$VIC_PARAM_ROOT` = the directory where VIC parameter files are stored
   `$VIC_PARAM_FILE` = the VIC parameter file corresponding to the desired year
   `$YEAR` = the desired year
   `$LCYEAR` = the desired land cover year (for $YEAR=(1992,2002,2011), $LCYEAR=(s1992,2001,2011))

14. Insert irrigation parameters into VIC parameter file

   `insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.$LCYEAR.2000_2016.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.$YEAR.nc -o $VIC_PARAM_ROOT/$VIC_PARAM_FILE_WITH_IRRIG`

   where

   `$MIRCA_ROOT` = the top-level directory under which MIRCA2000 data are stored
   `$VIC_PARAM_ROOT` = the directory where VIC parameter files are stored
   `$VIC_PARAM_FILE` = the VIC parameter file corresponding to the desired year
   `$YEAR` = the desired year
   `$LCYEAR` = the desired land cover year (for $YEAR=(1992,2002,2011), $LCYEAR=(s1992,2001,2011))
   `$VIC_PARAM_FILE_WITH_IRRIG` = the VIC parameter file corresponding to the desired year, with irrigation parameters added to it

### Example

An example of how to run these scripts is in `$PROJECT/examples/batch.process_mirca.csh`.

To run this batch file, replace all instances of `$PROJECT` in the file with the path to your copy of this GitHub project (and replace any other placeholders beginning with `$` with the appropriate values). Also make sure it is executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file.

