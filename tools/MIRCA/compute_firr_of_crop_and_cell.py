#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
    firr_of_fplant_file = ''
    fplant_file = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:f:v:o:n:",["firr_of_fplant_file=","fplant_file=","varnamelist=","outfile=","outnamelist="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <firr_of_fplant_file> -f <fplant_file> -v <varnamelist> -o <outfile> -n <outnamelist>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <firr_of_fplant_file> -f <fplant_file> -v <varnamelist> -o <outfile> -n <outnamelist>')
            sys.exit()
        elif opt in ("-i", "--firr_of_fplant_file"):
            firr_of_fplant_file = arg
        elif opt in ("-f", "--fplant_file"):
            fplant_of_crop_file = arg
        elif opt in ("-v", "-varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg
        elif opt in ("-n", "-outnamelist"):
            outnamelist = arg
            outnames = outnamelist.lstrip().rstrip().split(',')

    # Open fplant_of_crop file
    ds = xr.open_dataset(fplant_of_crop_file)

    data = {}
    first = True
    for varname in varnames:
        data[varname] = np.empty(ds[varname].shape)
        data[varname][:] = ds[varname][:]
        if first:
            nTime = data[varname].shape[0]
            nLat = data[varname].shape[1]
            nLon = data[varname].shape[2]
        first = False
    lat = np.empty([nLat])
    lon = np.empty([nLon])
    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]

    ds.close()

    # Open firr_of_fplant file
    ds = xr.open_dataset(firr_of_fplant_file)

    firr_of_fplant = np.empty(ds['firr_of_fplant'].shape)
    firr_of_fplant[:] = ds['firr_of_fplant'][:]

    ds.close()

    data_out = {}
    for j in range(len(outnames)):
        data_out[outnames[j]] = np.empty([nTime,nLat,nLon])
        for i in range(nTime):
            data_out[outnames[j]][i] = firr_of_fplant[i] * data[varnames[j]][i]

    month = np.arange(nTime).astype(int)

    # Write to output file
    ds = xr.Dataset(
        coords = {
            'month': (['month'], month),
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
    for outname in outnames:
        ds[outname] = (['month','lat','lon'],data_out[outname])
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

