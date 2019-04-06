#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def rescale_by_max(x):

    x_max = np.max(x)
    if x_max > 1:
        x /= x_max

    return(x)


def main():
    fplant_file = ''
    paramfile = ''
    cropclassidxlist = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:p:c:o:",["fplant_file=","paramfile=","cropclassidxlist=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -f <fplant_file> -p <paramfile> -c <cropclassidxlist> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -f <fplant_file> -p <paramfile> -c <cropclassidxlist> -o <outfile>')
            sys.exit()
        elif opt in ("-f", "--fplant_file"):
            fplant_file = arg
        elif opt in ("-p", "--paramfile"):
            paramfile = arg
        elif opt in ("-c", "--cropclassidxlist"):
            cropclassidxlist = arg
            cropclassidxs = cropclassidxlist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open fplant file
    ds = xr.open_dataset(fplant_file)

    fplant = np.empty(ds['fplant_of_cell'].shape)
    fplant[:] = ds['fplant_of_cell'][:]
    nTime = fplant.shape[0]
    nLat = fplant.shape[1]
    nLon = fplant.shape[2]
    lat = np.empty([nLat])
    lon = np.empty([nLon])
    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]
    tmp = np.nanmax(fplant, axis=0)
    fplant_mask = np.where(tmp > 0, 1, 0)

    ds.close()

    # Open param file
    ds = xr.open_dataset(paramfile)

    Cv = np.empty(ds['Cv'].shape)
    Cv[:] = ds['Cv'][:]

    ds.close()

    fcrop = np.zeros([nLat,nLon])
    for cropclassidx in cropclassidxs:
        fcrop += Cv[int(cropclassidx)]
    fcrop_mask = np.where(fcrop > 0, 1, 0)
    fplant_of_crop_raw = np.empty([nTime,nLat,nLon])
    fplant_of_crop = np.empty([nTime,nLat,nLon])
    fplant_of_cell = np.empty([nTime,nLat,nLon])
    for i in range(nTime):
        fplant_of_crop_raw[i] = fplant[i] / fcrop
        fplant_of_crop_raw[i] = np.where(fplant_mask > 0, fplant_of_crop_raw[i], np.nan)
        fplant_of_crop_raw[i] = np.where(fcrop_mask > 0, fplant_of_crop_raw[i], np.nan)

    fplant_of_crop[:] = fplant_of_crop_raw[:]
    fplant_of_crop = np.apply_along_axis(rescale_by_max, 0, fplant_of_crop)
    for i in range(nTime):
        fplant_of_cell[i] = fplant_of_crop[i] * fcrop

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
    ds['fcrop'] = (['lat','lon'],fcrop)
    ds['fplant_of_crop_raw'] = (['month','lat','lon'],fplant_of_crop_raw)
    ds['fplant_of_crop'] = (['month','lat','lon'],fplant_of_crop)
    ds['fplant_of_cell'] = (['month','lat','lon'],fplant_of_cell)
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

