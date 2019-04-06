#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
    fplant_file = ''
    thresh = 0
    infile = ''
    varname = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:t:i:v:o:",["fplant_file=","thresh=","infile=","varname=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -f <fplant_file> -t <thresh> -i <infile> -v <varname> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -f <fplant_file> -t <thresh> -i <infile> -v <varname> -o <outfile>')
            sys.exit()
        elif opt in ("-f", "--fplant_file"):
            fplant_file = arg
        elif opt in ("-t", "--thresh"):
            thresh = float(arg)
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-v", "--varname"):
            varname = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open fplant file
    ds = xr.open_dataset(fplant_file)

    fplant = np.empty(ds['fplant_of_cell'].shape)
    fplant[:] = ds['fplant_of_cell'][:]
    fplant_max = np.nanmax(fplant, axis=0)
    mask = np.where(fplant_max < thresh, 0, 1).astype(int)
    mask = np.where(np.isnan(fplant_max), 0, mask).astype(int)

    ds.close()

    # Open input file
    ds = xr.open_dataset(infile)

    data = np.empty(ds[varname].shape)
    data[:] = ds[varname][:]
    dim1 = data.shape[0]
    for i in range(dim1):
        data[i] = np.where(mask == 0, np.nan, data[i])
        data[i] = np.where(data[i] < 0, 0, data[i])
    ds[varname][:] = data[:]

    # Write to output file
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

