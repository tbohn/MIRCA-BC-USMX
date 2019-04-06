#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
    infile = ''
    lc_table = ''
    fplant_firr_file = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:t:f:o:",["infile=","lc_table=","fplant_firr_file=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -t <lc_table> -f <fplant_firr_file> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -t <lc_table> -f <fplant_firr_file> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--lc_table"):
            lc_table = arg
        elif opt in ("-f", "--fplant_firr_file"):
            fplant_firr_file = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open fplant_firr file
    ds = xr.open_dataset(fplant_firr_file)
    fplant_in = np.empty(ds['fplant'].shape)
    fplant_in[:] = ds['fplant'][:]
    firr_in = np.empty(ds['firr'].shape)
    firr_in[:] = ds['firr'][:]
    nTime = fplant_in.shape[0]
    ds.close()

    # Read lc table to get irrig information
    sm_code_dict = {
        'SAT': 0,
        'FC': 1,
        'CR': 2,
        'WP': 3,
    }
    lc_data = np.loadtxt(lc_table, dtype=bytes, delimiter=',').astype(str)
    nClass = lc_data.shape[0]
    ithresh = {}
    itarget = {}
    crop_split = {}
    irr_active = {}
    for v in range(nClass):
        classID = int(lc_data[v,0])
        ithresh[classID] = sm_code_dict[lc_data[v,5]]
        itarget[classID] = sm_code_dict[lc_data[v,6]]
        crop_split[classID] = int(lc_data[v,7])
        irr_active[classID] = int(lc_data[v,8])

    # Open input param file
    ds = xr.open_dataset(infile)

    Cv = np.empty(ds['Cv'].shape)
    Cv[:] = ds['Cv'][:]
    nClass = Cv.shape[0]
    nLat = Cv.shape[1]
    nLon = Cv.shape[2]
    classID = np.empty(ds['veg_class'].shape, dtype=np.int32)
    classID[:] = ds['veg_class'][:]
    lat = np.empty([nLat])
    lon = np.empty([nLon])
    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]
    mask = np.empty([nLat,nLon])
    mask[:] = ds['mask'][:]
    mask_cv = np.where(np.isfinite(Cv), Cv, 0)
    mask_cv = np.where(mask_cv > 0, 1, 0)
    fcrop = np.zeros([nLat,nLon])
    for v in range(nClass):
        if crop_split[classID[v]] == 1:
            fcrop += Cv[v]

    # Normalize fplant and firr by fcrop
    for t in range(nTime):
        fplant_in[t] /= fcrop
        fplant_in[t] = np.where(fcrop <= 0, 0, fplant_in[t])
        fplant_in[t] = np.where(fplant_in[t] > 1, 1, fplant_in[t])
        fplant_in[t] = np.where(mask == 0, np.nan, fplant_in[t])
        firr_in[t] /= fcrop
        firr_in[t] = np.where(fcrop <= 0, 0, firr_in[t])
        firr_in[t] = np.where(firr_in[t] > fplant_in[t], fplant_in[t], firr_in[t])
        firr_in[t] = np.where(mask == 0, np.nan, firr_in[t])

    # Assign to output variables
    fplant_out = np.empty([nClass,nTime,nLat,nLon])
    firr_out = np.empty([nClass,nTime,nLat,nLon])
    for v in range(nClass):
        if crop_split[classID[v]] == 1:
            fplant_out[v] = fplant_in[:]
            firr_out[v] = firr_in[:]
        else:
            fplant_out[v] = 0
            firr_out[v] = 0
        fplant_out[v] = np.where(mask_cv[v] == 1, fplant_out[v], np.nan)
        firr_out[v] = np.where(mask_cv[v] == 1, firr_out[v], np.nan)

    # Check for invalid values
    data = np.sum(firr_out, axis=1)
    mask_data = np.where(np.isfinite(data), 1, 0)
    a = np.where(mask_cv != mask_data)[0]
    b = np.where(mask_cv != mask_data)[1]
    c = np.where(mask_cv != mask_data)[2]
    for v,i,j in zip (a,b,c):
        for t in range(12):
            print(v,i,j,'bad','Cv',Cv[v,i,j],'t',t,'firr',firr_out[v,t,i,j])
        firr_out[v,:,i,j] = 0
    data = np.sum(fplant_out, axis=1)
    mask_data = np.where(np.isfinite(data), 1, 0)
    a = np.where(mask_cv != mask_data)[0]
    b = np.where(mask_cv != mask_data)[1]
    c = np.where(mask_cv != mask_data)[2]
    for v,i,j in zip (a,b,c):
        for t in range(12):
            print(v,i,j,'bad','Cv',Cv[v,i,j],'t',t,'fplant',fplant_out[v,t,i,j])
        fplant_out[v,:,i,j] = 1

    # Other irrig params
    ithresh_out = np.empty([nClass,nLat,nLon])
    itarget_out = np.empty([nClass,nLat,nLon])
    crop_split_out = np.empty([nClass,nLat,nLon])
    irr_active_out = np.empty([nClass,nLat,nLon])
    for v in range(nClass):
        ithresh_out[v] =  ithresh[classID[v]]
        itarget_out[v] =  itarget[classID[v]]
        crop_split_out[v] =  crop_split[classID[v]]
        irr_active_out[v] =  irr_active[classID[v]]

    # Add variables to dataset
    ds['ithresh'] = (['veg_class','lat','lon'], ithresh_out)
    ds['ithresh'].attrs['long_name'] = 'Soil moisture threshold (code) that triggers irrigation'
    ds['ithresh'].attrs['units'] = 'integer code'
    ds['ithresh'].encoding['zlib'] = True
    ds['itarget'] = (['veg_class','lat','lon'], itarget_out)
    ds['itarget'].attrs['long_name'] = 'Soil moisture threshold (code) that ends irrigation'
    ds['itarget'].attrs['units'] = 'integer code'
    ds['itarget'].encoding['zlib'] = True
    ds['crop_split'] = (['veg_class','lat','lon'], crop_split_out)
    ds['crop_split'].attrs['long_name'] = 'Flag (0 or 1) indicating whether to split veg tile into fallow and planted portions'
    ds['crop_split'].attrs['units'] = 'integer code'
    ds['crop_split'].encoding['zlib'] = True
    ds['irr_active'] = (['veg_class','lat','lon'], irr_active_out)
    ds['irr_active'].attrs['long_name'] = 'Flag (0 or 1) indicating whether to run irrgation model in this veg tile'
    ds['irr_active'].attrs['units'] = 'integer code'
    ds['irr_active'].encoding['zlib'] = True
    ds['fcrop'] = (['veg_class','month','lat','lon'], fplant_out)
    ds['fcrop'].attrs['long_name'] = 'Fraction of veg tile that is planted with crops'
    ds['fcrop'].attrs['units'] = 'area fraction'
    ds['fcrop'].encoding['zlib'] = True
    ds['firr'] = (['veg_class','month','lat','lon'], firr_out)
    ds['firr'].attrs['long_name'] = 'Fraction of veg tile that is irrigated'
    ds['firr'].attrs['units'] = 'area fraction'
    ds['firr'].encoding['zlib'] = True

    # Write to output file
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

