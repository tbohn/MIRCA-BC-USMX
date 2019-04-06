## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $MIRCA_ROOT with the path to the top-level directory where MIRCA2000 data are stored
## - replace $VIC_PARAM_ROOT with the path to the top-level directory where VIC parameters are stored
## - replace $USDA_ROOT with the path to the top-level directory where USDA files are stored
## - replace $SAGARPA_ROOT with the path to the top-level directory where SAGARPA files are stored
##
## Also, you must run this script from $PROJECT/tools/MIRCA or add $PROJECT/tools/MIRCA to your $PATH environment variable.

## Compute total areas for irrigated, rainfed, and all crops
add_crop_areas.pl $MIRCA_ROOT/asc.orig $MIRCA_ROOT/asc.crop_groups >& log.add_crop_areas.pl.txt

## Convert from areas to area fractions
convert_areas_to_fractions_of_cell_area.pl $MIRCA_ROOT/asc.crop_groups $MIRCA_ROOT/asc.frac >& log.convert_areas_to_fractions_of_cell_area.pl.txt

## Clip to domain and resample to 1/16 degree resolution
wrap_gridclip.asc.pl $MIRCA_ROOT/asc.frac crop float 14 50 -125 -65 16 4 $MIRCA_ROOT/USMX/asc.0.08333333_deg.frac_of_cell >& log.wrap_gridclip.asc.pl.txt
wrap_grid_subsample.pl $MIRCA_ROOT/USMX/asc.0.08333333_deg.frac_of_cell crop float 0.0208333333333333 nn 0 float 16 4 $MIRCA_ROOT/USMX/asc.0.02083333_deg.frac_of_cell >& log.wrap_grid_subsample.pl.txt
wrap_grid_agg.pl $MIRCA_ROOT/USMX/asc.0.02083333_deg.frac_of_cell crop float 3 5 4 $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell.tmp
wrap_gridclip.asc.pl $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell.tmp crop float 14 50 -125 -65 5 4 $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell >& log.wrap_gridclip.asc.pl.txt
rm -rf $MIRCA_ROOT/USMX/asc.0.02083333_deg.frac_of_cell $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell.tmp

## Compute fplant and firr
wrap_compute_firr_of_fplant_from_mirca.pl $MIRCA_ROOT/USMX/asc.0.0625_deg.frac_of_cell $MIRCA_ROOT/USMX/asc.0.0625_deg.fplant_and_firr >& log.wrap_compute_firr_of_fplant_from_mirca.pl.txt

## Convert to netcdf
mkdir -p $MIRCA_ROOT/USMX/nc
convert_monthly_single_var_asc_to_nc.py -i $MIRCA_ROOT/USMX/asc.0.0625_deg.fplant_and_firr -p firr_of_fplant_0 -v firr_of_fplant -o $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc >& log.convert_monthly_single_var_asc_to_nc.py.txt
convert_monthly_single_var_asc_to_nc.py -i $MIRCA_ROOT/USMX/asc.0.0625_deg.fplant_and_firr -p fplant_of_cell_0 -v fplant_of_cell -o $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc >& log.convert_monthly_single_var_asc_to_nc.py.txt

## Null out low values
null_out_mirca_where_fplant_low.py -f $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -t 0.01 -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc -v firr_of_fplant -o tmp.nc
mv tmp.nc $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc >& log.null_out_mirca_where_fplant_low.py.txt
null_out_mirca_where_fplant_low.py -f $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -t 0.01 -i $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -v fplant_of_cell -o tmp.nc
mv tmp.nc $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc >& log.null_out_mirca_where_fplant_low.py.txt

## Compute the planted fraction relative to NLCD_INEGI 2001 cropland areas
compute_fplant_of_crop_mirca_lc.py -f $MIRCA_ROOT/USMX/nc/fplant_of_cell.nc -p $VIC_PARAM_ROOT/params.USMX.2001.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc >& log.compute_fplant_of_crop_mirca_lc.py.txt
compute_annual_mirca_fplant.py -f $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc -v fplant_of_cell -o $MIRCA_ROOT/USMX/nc/fplant_annual.nc >& log.compute_annual_mirca_fplant.py.txt

## Fill gaps in MIRCA planted areas (where NLCD cropland area of given year > 0)
gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc -v fplant_of_crop,fplant_of_cell -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/params.USMX.2001.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/fplant.2001.nc >& log.gap_fill_firr_and_fplant.py.txt
gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc -v fplant_of_crop,fplant_of_cell -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/params.USMX.2011.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/fplant.2011.nc >& log.gap_fill_firr_and_fplant.py.txt
gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/fplant_of_crop.nc -v fplant_of_crop,fplant_of_cell -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/params.USMX.s1992.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/fplant.s1992.nc >& log.gap_fill_firr_and_fplant.py.txt

## Fill gaps in MIRCA irrigated areas (where MIRCA gap-filled planted areas > 0)
gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc -v firr_of_fplant -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/params.USMX.2001.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/firr_of_fplant.2001.nc >& log.gap_fill_firr_and_fplant.py.txt
gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc -v firr_of_fplant -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/params.USMX.2011.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/firr_of_fplant.2011.nc >& log.gap_fill_firr_and_fplant.py.txt
gap_fill_firr_and_fplant.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.nc -v firr_of_fplant -f $MIRCA_ROOT/USMX/nc/fplant_annual.nc -p $VIC_PARAM_ROOT/params.USMX.s1992.2000_2016.nc -c 12,13 -o $MIRCA_ROOT/USMX/nc/firr_of_fplant.s1992.nc >& log.gap_fill_firr_and_fplant.py.txt

## Convert irrigated fractions of planted areas to irrigated fractions of (a) NLCD cropland area and (b) cell area
compute_firr_of_crop_and_cell.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.2001.nc -f $MIRCA_ROOT/USMX/nc/fplant.2001.nc -v fplant_of_crop,fplant_of_cell -o $MIRCA_ROOT/USMX/nc/firr.2001.nc -n firr_of_crop,firr_of_cell >& log.compute_firr_of_crop_and_cell.py.txt
compute_firr_of_crop_and_cell.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.2011.nc -f $MIRCA_ROOT/USMX/nc/fplant.2011.nc -v fplant_of_crop,fplant_of_cell -o $MIRCA_ROOT/USMX/nc/firr.2011.nc -n firr_of_crop,firr_of_cell >& log.compute_firr_of_crop_and_cell.py.txt
compute_firr_of_crop_and_cell.py -i $MIRCA_ROOT/USMX/nc/firr_of_fplant.s1992.nc -f $MIRCA_ROOT/USMX/nc/fplant.s1992.nc -v fplant_of_crop,fplant_of_cell -o $MIRCA_ROOT/USMX/nc/firr.s1992.nc -n firr_of_crop,firr_of_cell >& log.compute_firr_of_crop_and_cell.py.txt

## Bias correct the planted and irrigated fractions so that their total areas sum to the state-by-state values in government records
bias_correct_mirca_fplant_firr.py -v $VIC_PARAM_ROOT/params.USMX.s1992.2000_2016.nc -c 12,13 -p $MIRCA_ROOT/USMX/nc/fplant.s1992.nc -i $MIRCA_ROOT/USMX/nc/firr.s1992.nc -u $USDA_ROOT/summ/summ.1992.csv -s $SAGARPA_ROOT/summ_by_mun.csv -m $PROJECT/data/mun_us.0.01_deg.asc -x $PROJECT/data/mun_mx.0.01_deg.asc -y 1992 -o $MIRCA_ROOT/USMX/nc/fplant_firr_bc.s1992.nc >& log.bias_correct_mirca_fplant_firr.py.1992.txt
bias_correct_mirca_fplant_firr.py -v $VIC_PARAM_ROOT/params.USMX.2001.2000_2016.nc -c 12,13 -p $MIRCA_ROOT/USMX/nc/fplant.2001.nc -i $MIRCA_ROOT/USMX/nc/firr.2001.nc -u $USDA_ROOT/summ/summ.2002.csv -s $SAGARPA_ROOT/summ_by_mun.csv -m $PROJECT/data/mun_us.0.01_deg.asc -x $PROJECT/data/mun_mx.0.01_deg.asc -y 2002 -o $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2001.nc >& log.bias_correct_mirca_fplant_firr.py.2001.txt
bias_correct_mirca_fplant_firr.py -v $VIC_PARAM_ROOT/params.USMX.2011.2000_2016.nc -c 12,13 -p $MIRCA_ROOT/USMX/nc/fplant.2011.nc -i $MIRCA_ROOT/USMX/nc/firr.2011.nc -u $USDA_ROOT/summ/summ.2012.csv -s $SAGARPA_ROOT/summ_by_mun.csv -m $PROJECT/data/mun_us.0.01_deg.asc -x $PROJECT/data/mun_mx.0.01_deg.asc -y 2012 -o $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2011.nc >& log.bias_correct_mirca_fplant_firr.py.2011.txt

## Special Case: Bias correct the planted and irrigated fractions state-by-state in the US and county-by-county in Mexico; only available for years >= 2003
bias_correct_mirca_fplant_firr.py -v $VIC_PARAM_ROOT/params.USMX.2011.2000_2016.nc -c 12,13 -p $MIRCA_ROOT/USMX/nc/fplant.2011.nc -i $MIRCA_ROOT/USMX/nc/firr.2011.nc -u $USDA_ROOT/summ/summ.2012.csv -s $SAGARPA_ROOT/summ_by_mun.csv -m $PROJECT/data/mun_us.0.01_deg.asc -x $PROJECT/data/mun_mx.0.01_deg.asc -y 2012 -d -o $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2011.mun_mx.nc >& log.bias_correct_mirca_fplant_firr.py.2011.mun_mx.txt

## Insert irrigation parameters into VIC parameter file
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.s1992.2000_2016.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.1992.nc -o $VIC_PARAM_ROOT/params/params.USMX.s1992.2000_2016.with_irrig.nc >& log.insert_irrig_params.py.s1992.txt
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.s2001.2000_2016.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2001.nc -o $VIC_PARAM_ROOT/params/params.USMX.s2001.2000_2016.with_irrig.nc >& log.insert_irrig_params.py.s2001.txt
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.2001.2001_2001.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2001.nc -o $VIC_PARAM_ROOT/params/params.USMX.2001.2001_2001.with_irrig.nc >& log.insert_irrig_params.py.2001.txt
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.s2011.2000_2016.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2011.nc -o $VIC_PARAM_ROOT/params/params.USMX.s2011.2000_2016.with_irrig.nc >& log.insert_irrig_params.py.s2011.txt
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.2011.2011_2011.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2011.nc -o $VIC_PARAM_ROOT/params/params.USMX.2011.2011_2011.with_irrig.nc >& log.insert_irrig_params.py.2011.txt
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.s2011.2000_2016.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2011.mun_mx.nc -o $VIC_PARAM_ROOT/params/params.USMX.s2011.2000_2016.with_irrig.mun_mx.nc >& log.insert_irrig_params.py.s2011.mun_mx.txt
insert_irrig_params.py -i $VIC_PARAM_ROOT/params/params.USMX.2011.2011_2011.nc -t $PROJECT/data/irr_table.USMX.NLCD_INEGI.csv -f $MIRCA_ROOT/USMX/nc/fplant_firr_bc.2011.mun_mx.nc -o $VIC_PARAM_ROOT/params/params.USMX.2011.2011_2011.with_irrig.mun_mx.nc >& log.insert_irrig_params.py.2011.mun_mx.txt
