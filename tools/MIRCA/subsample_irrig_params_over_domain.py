#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def subsample_variable(data_in, shape, dtype, llcorner_in,
                       cellsize_in, llcorner, cellsize, mask, fill_value):

    # Replace nans with fill_value
    if ((dtype == 'float' or dtype == 'np.single' or dtype == 'double' or
         dtype == 'int' or dtype == 'np.int32' or dtype == 'np.int64') and
        np.any(np.isnan(data_in))):
        tmp = np.where(np.isnan(data_in),fill_value,data_in)
        data_in = tmp

    # Copy values to output data array
    # Using an explicit loop to account for different grid resolutions
    shape_in = data_in.shape
    data = np.full(shape, fill_value, dtype=dtype)
    if cellsize == cellsize_in:
        y0_off = int((llcorner_in[0] - llcorner[0]) / cellsize)
        y1_off = y0_off + shape_in[-2] - shape[-2]
        y0 = max(y0_off, 0)
        y1 = min((shape[-2] + y1_off), shape[-2])
        y0_in = max(-y0_off, 0)
        y1_in = min((shape_in[-2] - y1_off), shape_in[-2])
        x0_off = int((llcorner_in[1] - llcorner[1]) / cellsize)
        x1_off = x0_off + shape_in[-1] - shape[-1]
        x0 = max(x0_off, 0)
        x1 = min((shape[-1] + x1_off), shape[-1])
        x0_in = max(-x0_off, 0)
        x1_in = min((shape_in[-1] - x1_off), shape_in[-1])
        data[...,y0:y1,x0:x1] = data_in[...,y0_in:y1_in,x0_in:x1_in].astype(dtype)
    else:
        for y in range(shape[-2]):
            ctr_lat = llcorner[0] + (y + 0.5) * cellsize
            y_in = int((ctr_lat - llcorner_in[0]) / cellsize_in)
            if y_in < 0 or y_in >= shape_in[-2]:
                continue
            for x in range(shape[-1]):
                ctr_lon = llcorner[1] + (x + 0.5) * cellsize
                x_in = int((ctr_lon - llcorner_in[1]) / cellsize_in)
                if x_in < 0 or x_in >= shape_in[-1]:
                    continue
                if len(shape) == 2:
                  data[...,y,x] = np.asscalar(data_in[...,y_in,x_in].astype(dtype))
                else:
                  data[...,y,x] = data_in[...,y_in,x_in].astype(dtype)

    # Overwrite values from data_in with fill_value where mask is 0
    data[..., mask != 1] = fill_value

    return data


def main():
    infile = ''
    domainfile = ''
    outfile = ''
    nodata_default = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:d:no:",["infile=","domainfile=","nodata_default","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -d <domainfile> [-n] -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -d <domainfile> [-n] -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-n", "--nodata_default"):
            nodata_default = True
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Define fill values
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-1)
    fill_value_str = ''
    fill_value_float = 9.96920996838687e+36

    # Read domain file
    ds = xr.open_dataset(domainfile)

    # Read relevant dims and vars
    lat = ds['lat']
    lon = ds['lon']
    mask = np.where(ds['mask']==1,1,0).astype(np.int32)
    nLat = len(lat)
    nLon = len(lon)
    minlat = np.asscalar(lat[0])
    maxlat = np.asscalar(lat[-1])
    minlon = np.asscalar(lon[0])
    maxlon = np.asscalar(lon[-1])
    cellsize = ( maxlat - minlat ) / (nLat - 1)
    minlat -= (0.5 * cellsize)
    minlon -= (0.5 * cellsize)
    maxlat += (0.5 * cellsize)
    maxlon += (0.5 * cellsize)
    llcorner = [minlat, minlon]

    # Read infile
    ds_in = xr.open_dataset(infile)

    # Read coordinate variables
    lat_in = ds_in['lat']
    lon_in = ds_in['lon']
    nLat_in = len(lat_in)
    nLon_in = len(lon_in)
    minlat_in = np.asscalar(lat_in[0])
    maxlat_in = np.asscalar(lat_in[-1])
    minlon_in = np.asscalar(lon_in[0])
    maxlon_in = np.asscalar(lon_in[-1])
    cellsize_in = ( maxlat_in - minlat_in ) / (nLat_in - 1)
    minlat_in -= (0.5 * cellsize_in)
    minlon_in -= (0.5 * cellsize_in)
    maxlat_in += (0.5 * cellsize_in)
    maxlon_in += (0.5 * cellsize_in)
    llcorner_in = [minlat_in, minlon_in]

    month = ds_in['month']
    nMonth = len(month)

    # Variables
    varnames_1d = []
    varnames_3d_irr = [
		       'fplant',
                       'firr',
                      ]

    varnames = []
    for varname in varnames_3d_irr:
        varnames.append(varname)

    varnames_int = []
    varnames_str = []

    # Write output file
    ds_out = xr.Dataset(
                        {
                         'mask': (['lat','lon'],mask),
                        },
                        coords={
                                'lat': (['lat'],lat),
                                'lon': (['lon'],lon),
                                'month': (['month'],month),
                               }
                       )

    for varname in ds.coords:
        ds_out[varname].attrs = ds[varname].attrs
        ds_out[varname].encoding = ds[varname].encoding
    ds_out['mask'].attrs = ds['mask'].attrs
    if nodata_default:
        ds_out['mask'].encoding['_FillValue'] = fill_value_mask
    else:
        ds_out['mask'].encoding['_FillValue'] = ds['mask'].encoding['_FillValue']

    # Do the subsampling
    out_dict = {}
    for varname in varnames:

        print('varname',varname)
        if varname in varnames_3d_irr:
            print('3d_irr var:',varname)
            shape = [nMonth,nLat,nLon]
            dimstr = ['month','lat','lon']
        else:
            print('no case for',varname)

        if varname in varnames_int:
            dtype = np.int32
            fill_value = fill_value_int
        elif varname in varnames_str:
            dtype = np.str
            fill_value = fill_value_str
        else:
            dtype = np.single
            fill_value = fill_value_float

        if (varname in varnames_1d):
            out_dict[varname] = ds_in[varname]
        else:
            out_dict[varname] = subsample_variable(ds_in[varname],
                                                   shape,
                                                   dtype,
                                                   llcorner_in,
                                                   cellsize_in,
                                                   llcorner,
                                                   cellsize,
                                                   mask, fill_value)

        ds_out[varname] = (dimstr, out_dict[varname])
        ds_out[varname].attrs = ds_in[varname].attrs
        if nodata_default:
            if varname in varnames_int:
                ds_out[varname].encoding['_FillValue'] = fill_value_int
            elif varname not in varnames_str:
                ds_out[varname].encoding['_FillValue'] = fill_value_float
        else:
            ds_out[varname].encoding = ds_in[varname].encoding

    print('writing to',outfile)
    ds_out.to_netcdf(outfile)

    ds_in.close()
    ds_out.close()
    ds.close()

if __name__ == "__main__":
    main()    
