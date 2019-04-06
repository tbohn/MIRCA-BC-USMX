#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
    fplant_file = ''
    varname = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:v:o:",["fplant_file=","varname=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -f <fplant_file> -v <varname> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -f <fplant_file> -v <varname> -o <outfile>')
            sys.exit()
        elif opt in ("-f", "--fplant_file"):
            fplant_file = arg
        elif opt in ("-v", "--varname"):
            varname = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open fplant file
    ds = xr.open_dataset(fplant_file)

    fplant = np.empty(ds[varname].shape)
    fplant[:] = ds[varname][:]
    nTime = fplant.shape[0]
    nLat = fplant.shape[1]
    nLon = fplant.shape[2]
    lat = np.empty([nLat])
    lon = np.empty([nLon])
    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]

    ds.close()

    tmp1 = np.max(fplant[0:3], axis=0)
    tmp2 = np.max(fplant[5:8], axis=0)
    fplant_annual = tmp1 + tmp2

    # Write to output file
    ds = xr.Dataset(
        coords = {
            'lat': (['lat'], lat),
            'lon': (['lon'], lon),
        }
    )
    ds['lat'].attrs['long_name'] = 'Grid cell center latitude'
    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'degrees North'
    ds['lat'].attrs['axis'] = 'Y'
    ds['lon'].attrs['long_name'] = 'Grid cell center longitude'
    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['units'] = 'degrees East'
    ds['lon'].attrs['axis'] = 'X'
    ds['fplant_annual'] = (['lat','lon'],fplant_annual)
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

