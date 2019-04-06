#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import pandas as pd

def read_ascii_gridfile(filename, dtype):

    # Read header
    f = open(filename, 'r')
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    header4 = f.readline()
    header5 = f.readline()
    header6 = f.readline()
    ncols = int(header1.partition(' ')[2])
    nrows = int(header2.partition(' ')[2])
    xllcorner = float(header3.partition(' ')[2])
    yllcorner = float(header4.partition(' ')[2])
    cellsize = float(header5.partition(' ')[2])
    nodata = float(header6.partition(' ')[2])
    f.close()

    # Allocate data
    data = np.empty((nrows,ncols), dtype=dtype)

    # Read data
    tmp = np.loadtxt(filename,dtype=dtype,delimiter =' ', skiprows=6)
    for i in range(nrows):
        data[nrows-1-i] = tmp[i].copy()

    llcorner = [yllcorner, xllcorner]

    return (data, nrows, ncols, llcorner, cellsize, nodata)


def main():
    indir = ''
    prefix = ''
    varname = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:p:v:o:",["indir=","prefix=","varname=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <indir> -p <prefix> -v <varname> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <indir> -p <prefix> -v <varname> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--indir"):
            indir = arg
        elif opt in ("-p", "--prefix"):
            prefix = arg
        elif opt in ("-v", "--varname"):
            varname = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Read input files
    for m in range(12):
        monthstr = '{:02d}'.format(m+1)
        file = indir + '/' + prefix + monthstr + '.asc'
        (data_tmp, nrows, ncols, llcorner, cellsize,
            nodata) = read_ascii_gridfile(file, np.float)
        if m == 0:
            data = np.empty([12,nrows,ncols])
            lat = np.empty([nrows])
            lon = np.empty([ncols])
        data[m] = data_tmp

    for y in range(nrows):
        lat[y] = llcorner[0] + (y + 0.5) * cellsize

    for x in range(ncols):
        lon[x] = llcorner[1] + (x + 0.5) * cellsize

    timevar = np.arange(12, dtype=np.int);

    # Write output file
    ds = xr.Dataset(
        coords = {
            'month': (['month'], timevar),
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
    for varname_tmp in ['lat','lon']:
        ds[varname_tmp].encoding['zlib'] = True
    ds[varname] = (['month','lat','lon'], data)
    ds[varname].encoding['zlib'] = True

    print('writing to',outfile)
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()
