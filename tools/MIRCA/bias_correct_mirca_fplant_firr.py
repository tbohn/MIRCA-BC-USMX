#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import pandas as pd
import pysal


def integrate_over_region(data, area, map, regionID):
    mask = np.where(map == int(regionID), 1, 0)
    tmp = data * area * mask
    sum = np.nansum(tmp)
    return sum


def rescale_over_region(data, map, regionID, factor):
    mask = np.where(map == int(regionID), 1, 0)
    tmp = data * mask * factor
    data = np.where(map == int(regionID), tmp, data)
    return data


def rescale_diff_over_region(data, diff, map, regionID, factor):
    diff = np.maximum(diff, 0)
    factor = max(factor, 0)
    factor = min(factor, 1)
    mask = np.where(map == int(regionID), 1, 0)
    tmp = diff * mask * (1-factor)
    data = np.where(map == int(regionID), data+tmp, data)
    return data


def interpolate_months(x):
    # Assumes 24 elements, with original data are at positions 12 to 23
    mbc1 = 1
    mbc2 = 7
    for m in range(12):
        if x[12+mbc2]-x[12+mbc1] != 0:
            x[m] = x[mbc1] + (x[mbc2]-x[mbc1]) * (x[12+m]-x[12+mbc1]) / (x[12+mbc2]-x[12+mbc1])
        else:
            x[m] = x[mbc1]
    return x


def aggregate_space(data_fine, res_ratio):
    nY = int( data_fine.shape[0] / res_ratio )
    nX = int( data_fine.shape[1] / res_ratio )
    data = np.empty([nY,nX])
    for i in range(nY):
        i0 = i * res_ratio
        i1 = i0 + res_ratio
        for j in range(nX):
            j0 = j * res_ratio
            j1 = j0 + res_ratio
            data[i,j] = np.nanmean(data_fine[i0:i1,j0:j1])
    return data


def subsample_variable(data_in, shape, dtype, llcorner_in,
                       cellsize_in, llcorner, cellsize, mask, fill_value):

    if not np.isnan(fill_value):
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


# function to make a gaussian kernel
def make_kernel(ny, nx, yctr, xctr, sigma_space):
    kernel = np.zeros([nx,ny])
    sum = 0
    for i in range(ny):
        for j in range(nx):
          if (sigma_space > 0):
              space_term = ((i-yctr)**2+(j-xctr)**2) / (2*sigma_space**2)
          else:
              space_term = 0
          kernel[i,j] = np.exp( -( space_term ) )
          sum += kernel[i,j]
    for i in range(ny):
        for j in range(nx):
            kernel[i,j] /= sum
    return kernel


# gapfilling with a gaussian kernel
def gapfill (mydata, mycount, gapmask, ss):

    ylen = mydata.shape[0]
    xlen = mydata.shape[1]

    tmpdata = np.empty(mydata.shape)
    tmpdata[:] = mydata[:]

    # Loop over isolated gaps within this time slice
    a = np.where(gapmask == True)[0]
    b = np.where(gapmask == True)[1]
    for i,j in zip(a,b):

        # Make kernel
        kernel = make_kernel(ylen, xlen, i, j, ss)

        # Spatial average coefficients
        weights = kernel * mycount
        sumwts = np.nansum(weights)

        # Weighted data
        wtdata = weights * mydata
        sumwtdata = np.nansum(wtdata)
        if (sumwts > 0):
            tmpdata[i,j] = sumwtdata / sumwts

    return tmpdata


def main():
    vicparam_file = ''
    cropclassidxlist = ''
    fplant_file = ''
    firr_file = ''
    usda_file = ''
    sagarpa_file = ''
    us_map_file = ''
    mx_map_file = ''
    year = -1
    do_mun_mx = False
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hv:c:p:i:u:s:m:x:y:do:",["vicparam_file=","cropclassidxlist=","fplant_file=","firr_file=","usda_file=","sagarpa_file=","us_map_file=","mx_map_file=","year=","do_mun_mx","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -v <vicparam_file> -c <cropclassidxlist> -p <fplant_file> -i <firr_file> -u <usda_file> -s <sagarpa_file> -m <us_map_file> -x <mx_map_file> -y <year> [-d] -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -v <vicparam_file> -c <cropclassidxlist> -p <fplant_file> -i <firr_file> -u <usda_file> -s <sagarpa_file> -m <us_map_file> -x <mx_map_file> -y <year> [-d] -o <outfile>')
            sys.exit()
        elif opt in ("-v", "--vicparam_file"):
            vicparam_file = arg
        elif opt in ("-c", "--cropclassidxlist"):
            cropclassidxlist = arg
            cropclassidxs = cropclassidxlist.lstrip().rstrip().split(',')
        elif opt in ("-p", "--fplant_file"):
            fplant_file = arg
        elif opt in ("-i", "--firr_file"):
            firr_file = arg
        elif opt in ("-u", "--usda_file"):
            usda_file = arg
        elif opt in ("-s", "--sagarpa_file"):
            sagarpa_file = arg
        elif opt in ("-m", "--us_map_file"):
            us_map_file = arg
        elif opt in ("-x", "--mx_map_file"):
            mx_map_file = arg
        elif opt in ("-y", "--year"):
            year = int(arg)
        elif opt in ("-d", "--do_mun_mx"):
            do_mun_mx = True
        elif opt in ("-o", "--outfile"):
            outfile = arg

    mbc1 = int(1) # feb
    mbc2 = int(7) # aug

    # Open vicparam_file
    print('reading',vicparam_file)
    ds = xr.open_dataset(vicparam_file)

    mask = np.empty(ds['mask'].shape, dtype=np.int)
    mask[:] = ds['mask'][:]
    Cv = np.empty(ds['Cv'].shape)
    Cv[:] = ds['Cv'][:]
    nLat = mask.shape[0]
    nLon = mask.shape[1]
    lat = np.empty([nLat])
    lon = np.empty([nLon])
    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]
    fcrop = np.zeros([nLat,nLon])
    for cropclassidx in cropclassidxs:
        fcrop += Cv[int(cropclassidx)]
    cellsize = 0.0625
    llcorner = [lat[0]-0.5*cellsize, lon[0]-0.5*cellsize]
    shape = [nLat, nLon]

    ds.close()

    # Subsample
    print('subsampling mask')
    cellsize_fine = cellsize / 25
    nLat_fine = nLat * 25
    nLon_fine = nLon * 25
    lat_fine = [lat[0] + (i+0.5)*cellsize_fine for i in range(nLat_fine)]
    lon_fine = [lon[0] + (i+0.5)*cellsize_fine for i in range(nLon_fine)]
    llcorner_fine = [lat_fine[0]-0.5*cellsize_fine, lon_fine[0]-0.5*cellsize_fine]
    shape_fine = [nLat_fine, nLon_fine]
    mask_fine = np.empty([nLat_fine,nLon_fine], dtype=np.int)
    for i in range(nLat_fine):
        ii = int(i/25)
        for j in range(nLon_fine):
            jj = int(j/25)
            mask_fine[i,j] = mask[ii,jj]

    print('subsampling fcrop')
    fcrop_fine = subsample_variable(fcrop, shape_fine, np.float, llcorner, cellsize, llcorner_fine, cellsize_fine, mask_fine, np.nan)

    # Open fplant_file
    print('reading',fplant_file)
    ds = xr.open_dataset(fplant_file)

    fplant = np.empty(ds['fplant_of_cell'].shape)
    fplant[:] = ds['fplant_of_cell'][:]
    nTime = ds['fplant_of_cell'].shape[0]

    # Subsample
    print('subsampling fplant')
    fplant_fine = np.empty([nTime,nLat_fine,nLon_fine])
    for t in range(nTime):
        fplant_fine[t] = subsample_variable(fplant[t], shape_fine, np.float, llcorner, cellsize, llcorner_fine, cellsize_fine, mask_fine, np.nan)

    ds.close()

    # Open firr_file
    print('reading',firr_file)
    ds = xr.open_dataset(firr_file)

    firr = np.empty(ds['firr_of_cell'].shape)
    firr[:] = ds['firr_of_cell'][:]

    # Subsample
    print('subsampling firr')
    firr_fine = np.empty([nTime,nLat_fine,nLon_fine])
    for t in range(nTime):
        firr_fine[t] = subsample_variable(firr[t], shape_fine, np.float, llcorner, cellsize, llcorner_fine, cellsize_fine, mask_fine, np.nan)

    ds.close()

    # Read usda_file
    print('reading',usda_file)
    usda = np.loadtxt(usda_file, dtype=bytes, delimiter=',').astype(str)
    nLines_usda = usda.shape[0]

    # Read sagarpa_file
    print('reading',sagarpa_file)
    sagarpa = np.loadtxt(sagarpa_file, dtype=bytes, delimiter=',').astype(str)
    nLines_sagarpa = sagarpa.shape[0]

    print('processing govt areas')
    dictStateIDs = {}
    state_mun_list = {}
    Aplant_gov = {}
    Airr_gov = {}
    for country in ['us', 'mx']:
        dictStateIDs[country] = {}
        state_mun_list[country] = {}
        Aplant_gov[country] = {}
        Airr_gov[country] = {}

    for i in range(nLines_usda):
        state_mun_ID = usda[i,0]
        if state_mun_ID[0] == '#':
            continue
        if state_mun_ID[0] == 'e':
            continue
        stateID = state_mun_ID[0:2]
        dictStateIDs['us'][stateID] = 1
        state_mun_list['us'][stateID] = 1

    for i in range(nLines_sagarpa):
        state_mun_ID = sagarpa[i,0]
        if state_mun_ID[0] == '#':
            continue
        if state_mun_ID[0] == 'e':
            continue
        stateID = state_mun_ID[0:2]
        dictStateIDs['mx'][stateID] = 1
        state_mun_list['mx'][stateID] = 1

    for country in ['us', 'mx']:
        for stateID in dictStateIDs[country].keys():
            Aplant_gov[country][stateID] = {}
            Aplant_gov[country][stateID]['000'] = {}
            Aplant_gov[country][stateID]['000']['a'] = 0
            Airr_gov[country][stateID] = {}
            Airr_gov[country][stateID]['000'] = {}
            state_mun_list[country][stateID] = {}
            Airr_gov[country][stateID]['000']['a'] = 0

    for i in range(nLines_usda):
        state_mun_ID = usda[i,0]
        if state_mun_ID[0] == '#':
            continue
        if state_mun_ID[0] == 'e':
            continue
        stateID = state_mun_ID[0:2]
        munID = state_mun_ID[2:5]
        state_mun_list['us'][stateID][munID] = 1
        Aplant_gov['us'][stateID][munID] = {}
        Airr_gov['us'][stateID][munID] = {}
        Aplant_gov['us'][stateID][munID]['a'] = float(usda[i,3])
        Aplant_gov['us'][stateID]['000']['a'] += Aplant_gov['us'][stateID][munID]['a']
        Airr_gov['us'][stateID][munID]['a'] = float(usda[i,4])
        Airr_gov['us'][stateID]['000']['a'] += Airr_gov['us'][stateID][munID]['a']

    for i in range(nLines_sagarpa):
        state_mun_ID = sagarpa[i,0]
        if state_mun_ID[0] == '#':
            continue
        if state_mun_ID[0] == 'e':
            continue
        if int(sagarpa[i,1]) != year:
            continue
        stateID = state_mun_ID[0:2]
        munID = state_mun_ID[2:5]
        state_mun_list['mx'][stateID][munID] = 1
        Aplant_gov['mx'][stateID][munID] = {}
        Airr_gov['mx'][stateID][munID] = {}
        Airr_gov['mx'][stateID][munID][mbc1] = float(sagarpa[i,13])
        Airr_gov['mx'][stateID][munID][mbc2] = float(sagarpa[i,16])
        Airr_gov['mx'][stateID][munID]['a'] = Airr_gov['mx'][stateID][munID][mbc1] + Airr_gov['mx'][stateID][munID][mbc2]
        Aplant_gov['mx'][stateID][munID][mbc1] = float(sagarpa[i,13]) + float(sagarpa[i,32])
        Aplant_gov['mx'][stateID][munID][mbc2] = float(sagarpa[i,16]) + float(sagarpa[i,35])
        Aplant_gov['mx'][stateID][munID]['a'] = Aplant_gov['mx'][stateID][munID][mbc1] + Aplant_gov['mx'][stateID][munID][mbc2]

    # Read us_map file
    print('reading',us_map_file)
    (us_map, nrows_us, ncols_us, llcorner_us, cellsize_us, nodata_us) = \
        read_ascii_gridfile(us_map_file, np.int32)
    shape_us_map = [nrows_us, ncols_us]

    # Subsample
    print('subsampling us map')
    map_fine = {}
    map_fine['us'] = subsample_variable(us_map, shape_fine, np.int, llcorner_us, cellsize_us, llcorner_fine, cellsize_fine, mask_fine, 0)

    # Read mx_map file
    print('reading',mx_map_file)
    (mx_map, nrows_mx, ncols_mx, llcorner_mx, cellsize_mx, nodata_mx) = \
        read_ascii_gridfile(mx_map_file, np.int32)
    shape_mx_map = [nrows_mx, ncols_mx]

    # Subsample
    print('subsampling mx map')
    map_fine['mx'] = subsample_variable(mx_map, shape_fine, np.int, llcorner_mx, cellsize_mx, llcorner_fine, cellsize_fine, mask_fine, 0)

    # Compute pixel areas (km2)
    print('computing pixel areas')
    area = np.empty([nLat_fine,nLon_fine])
    for i in range(nLat_fine):
        dx = pysal.cg.sphere.arcdist( (lon_fine[0] - (cellsize_fine / 2), lat_fine[i]), (lon_fine[0] + (cellsize_fine / 2), lat_fine[i]))
        dy = pysal.cg.sphere.arcdist( (lon_fine[0], lat_fine[i] - (cellsize_fine / 2)), (lon_fine[0], lat_fine[i] + (cellsize_fine / 2)))
        area[i] = dx * dy
    area = area.copy() * mask_fine
    area = np.where(area == 0, np.nan, area)

    # Bias correction: Compute Crop and Planted Areas
    diff_fplant = {}
    Aplant = {}
    Acrop = {}
    fplant_fine_bc = {}
    firr_fine_bc = {}
    for t in [mbc1, mbc2]:
        fplant_fine_bc[t] = np.empty([nLat_fine,nLon_fine])
        fplant_fine_bc[t][:] = fplant_fine[t,:,:]
        firr_fine_bc[t] = np.empty([nLat_fine,nLon_fine])
        firr_fine_bc[t][:] = firr_fine[t,:,:]
    for country in ['us', 'mx']:
        print('computing crop and planted areas for', country)
        diff_fplant[country] = {}
        Aplant[country] = {}
        Acrop[country] = {}
        for state in state_mun_list[country].keys():
            Acrop[country][state] = {}
            Aplant[country][state] = {}
            Acrop[country][state]['000'] = 0
            Aplant[country][state]['000'] = {}
            for t in [mbc1, mbc2, 'a']:
                Aplant[country][state]['000'][t] = 0
            for mun in state_mun_list[country][state].keys():
                if mun == '000':
                    continue
                smID = state+mun
                Acrop[country][state][mun] = integrate_over_region(fcrop_fine, area, map_fine[country], smID)
                Acrop[country][state]['000'] += Acrop[country][state][mun]
                Aplant[country][state][mun] = {}
                for t in [mbc1, mbc2, 'a']:
                    Aplant[country][state][mun][t] = 0
        for t in [mbc1, mbc2]:
            diff_fplant[country][t] = fcrop_fine - fplant_fine[t]
            for state in state_mun_list[country].keys():
                for mun in state_mun_list[country][state].keys():
                    if mun == '000':
                        continue
                    smID = state+mun
                    Aplant[country][state][mun][t] = integrate_over_region(fplant_fine[t], area, map_fine[country], smID)
                    Aplant[country][state]['000'][t] += Aplant[country][state][mun][t]
        for state in state_mun_list[country].keys():
            Aplant[country][state]['000']['a'] = Aplant[country][state]['000'][mbc1] + Aplant[country][state]['000'][mbc2]
            for mun in state_mun_list[country][state].keys():
                smID = state+mun
                Aplant[country][state][mun]['a'] = Aplant[country][state][mun][mbc1] + Aplant[country][state][mun][mbc2]


    # Bias Correction: adjust fplant (and adjust firr to ensure it doesn't exceed fplant)
    print('adjusting fplant')
    for country in ['us', 'mx']:
        print('country', country)
        for state in state_mun_list[country].keys():
            print('state', state)
            if Acrop[country][state]['000'] <= 0:
                continue
            if country == 'us':  # United States
                for mun in state_mun_list[country][state].keys():
                    if Acrop[country][state][mun] <= 0:
                        continue
                    smID = state+mun
#                    print('debug',country,state,mun,'a','Acrop',Acrop[country][state][mun],'Aplant',Aplant[country][state][mun]['a'],'Aplant_gov',Aplant_gov[country][state][mun]['a'])
                    if Aplant[country][state][mun]['a'] > Aplant_gov[country][state][mun]['a']:
                        if Aplant[country][state][mun]['a'] > 0:
                            factor = Aplant_gov[country][state][mun]['a'] / Aplant[country][state][mun]['a']
                            for t in [mbc1, mbc2]:
                                fplant_fine_bc[t] = rescale_over_region(fplant_fine_bc[t], map_fine[country], smID, factor)
                                firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor)
                    else:
                        # The factor of 2 here is necessary because we are comparing the sum of winter and summer planted areas to the sum of winter and summer crop areas
                        if 2*Acrop[country][state][mun]-Aplant[country][state][mun]['a'] > 0:
                            factor = (2*Acrop[country][state][mun]-Aplant_gov[country][state][mun]['a']) / (2*Acrop[country][state][mun]-Aplant[country][state][mun]['a'])
                            if Aplant[country][state][mun]['a'] > 0:
                                factor2 = Aplant_gov[country][state][mun]['a'] / Aplant[country][state][mun]['a']
                            for t in [mbc1, mbc2]:
                                fplant_fine_bc[t] = rescale_diff_over_region(fplant_fine_bc[t], diff_fplant[country][t], map_fine[country], smID, factor)
                                if Aplant[country][state][mun]['a'] > 0:
                                    firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor2)
            else:  # Mexico
                if year >= 2003 and do_mun_mx:
                    for mun in state_mun_list[country][state].keys():
                        if Acrop[country][state][mun] <= 0:
                            continue
                        smID = state+mun
                        for t in [mbc1, mbc2]:
                            if Aplant[country][state][mun][t] > Aplant_gov[country][state][mun][t]:
                                if Aplant[country][state][mun][t] > 0:
                                    factor = Aplant_gov[country][state][mun][t] / Aplant[country][state][mun][t]
                                    fplant_fine_bc[t] = rescale_over_region(fplant_fine_bc[t], map_fine[country], smID, factor)
                                    firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor)
                            else:
                                if Acrop[country][state][mun]-Aplant[country][state][mun][t] > 0:
                                    factor = (Acrop[country][state][mun]-Aplant_gov[country][state][mun][t]) / (Acrop[country][state][mun]-Aplant[country][state][mun][t])
                                    if Aplant[country][state][mun][t] > 0:
                                        factor2 = Aplant_gov[country][state][mun][t] / Aplant[country][state][mun][t]
                                        fplant_fine_bc[t] = rescale_diff_over_region(fplant_fine_bc[t], diff_fplant[country][t], map_fine[country], smID, factor)
                                        if Aplant[country][state][mun][t] > 0:
                                            firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor2)
                else:
                    if Acrop[country][state]['000'] <= 0:
                        continue
                    for t in [mbc1, mbc2]:
                        if Aplant[country][state]['000'][t] > Aplant_gov[country][state]['000'][t]:
                            if Aplant[country][state]['000'][t] > 0:
                                factor = Aplant_gov[country][state]['000'][t] / Aplant[country][state]['000'][t]
                                for mun in state_mun_list[country][state].keys():
                                    smID = state+mun
                                    fplant_fine_bc[t] = rescale_over_region(fplant_fine_bc[t], map_fine[country], smID, factor)
                                    firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor)
                        else:
                            if Acrop[country][state]['000']-Aplant[country][state]['000'][t] > 0:
                                factor = (Acrop[country][state]['000']-Aplant_gov[country][state]['000'][t]) / (Acrop[country][state]['000']-Aplant[country][state]['000'][t])
                                if Aplant[country][state]['000'][t] > 0:
                                    factor2 = Aplant_gov[country][state]['000'][t] / Aplant[country][state]['000'][t]
                                for mun in state_mun_list[country][state].keys():
                                    smID = state+mun
                                    fplant_fine_bc[t] = rescale_diff_over_region(fplant_fine_bc[t], diff_fplant[country][t], map_fine[country], smID, factor)
                                    if Aplant[country][state]['000'][t] > 0:
                                        firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor2)


    # Bias correction: Compute Irrigated Areas
    diff_firr = {}
    Airr = {}
    for country in ['us', 'mx']:
        print('computing irrigated areas for', country)
        Airr[country] = {}
        diff_firr[country] = {}
        for state in state_mun_list[country].keys():
            Airr[country][state] = {}
            Airr[country][state]['000'] = {}
            for t in [mbc1, mbc2, 'a']:
                Airr[country][state]['000'][t] = 0
            for mun in state_mun_list[country][state].keys():
                if mun == '000':
                    continue
                smID = state+mun
                Airr[country][state][mun] = {}
                for t in [mbc1, mbc2, 'a']:
                    Airr[country][state][mun][t] = 0
        for t in [mbc1, mbc2]:
            diff_firr[country][t] = fplant_fine_bc[t] - firr_fine_bc[t]
            for state in state_mun_list[country].keys():
                for mun in state_mun_list[country][state].keys():
                    if mun == '000':
                        continue
                    smID = state+mun
                    Airr[country][state][mun][t] = integrate_over_region(firr_fine_bc[t], area, map_fine[country], smID)
                    Airr[country][state]['000'][t] += Airr[country][state][mun][t]
        for state in state_mun_list[country].keys():
            Airr[country][state]['000']['a'] = Airr[country][state]['000'][mbc1] + Airr[country][state]['000'][mbc2]
            for mun in state_mun_list[country][state].keys():
                smID = state+mun
                Airr[country][state][mun]['a'] = Airr[country][state][mun][mbc1] + Airr[country][state][mun][mbc2]

    # Bias Correction: adjust firr
    print('adjusting firr')
    for country in ['us', 'mx']:
        print('country', country)
        for state in state_mun_list[country].keys():
            print('state', state)
            if Acrop[country][state]['000'] <= 0:
                continue
            if country == 'us':  # United States
                for mun in state_mun_list[country][state].keys():
                    if Acrop[country][state][mun] <= 0:
                        continue
                    smID = state+mun
                    if Airr[country][state][mun]['a'] > Airr_gov[country][state][mun]['a']:
                        if Airr[country][state][mun]['a'] > 0:
                            factor = Airr_gov[country][state][mun]['a'] / Airr[country][state][mun]['a']
                            for t in [mbc1, mbc2]:
                                firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor)
                    else:
                        for t in [mbc1, mbc2]:
                            if Aplant[country][state][mun][t]-Airr[country][state][mun][t] > 0:
                                factor = (Aplant[country][state][mun][t]-Airr_gov[country][state][mun]['a']) / (Aplant[country][state][mun][t]-Airr[country][state][mun][t])
                                firr_fine_bc[t] = rescale_diff_over_region(firr_fine_bc[t], diff_firr[country][t], map_fine[country], smID, factor)
            else:  # Mexico
                if year >= 2003:
                    for mun in state_mun_list[country][state].keys():
                        smID = state+mun
                        for t in [mbc1, mbc2]:
                            if Aplant[country][state][mun][t] <= 0:
                                continue
                            if Airr[country][state][mun][t] > Airr_gov[country][state][mun][t]:
                                if Airr[country][state][mun][t] > 0:
                                    factor = Airr_gov[country][state][mun][t] / Airr[country][state][mun][t]
                                    firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor)
                            else:
                                if Aplant[country][state][mun][t]-Airr[country][state][mun][t] > 0:
                                    factor = (Aplant[country][state][mun][t]-Airr_gov[country][state][mun][t]) / (Aplant[country][state][mun][t]-Airr[country][state][mun][t])
                                    firr_fine_bc[t] = rescale_diff_over_region(firr_fine_bc[t], diff_firr[country][t], map_fine[country], smID, factor)
                else:
                    for t in [mbc1, mbc2]:
                        if Aplant[country][state]['000'][t] <= 0:
                            continue
                        if Airr[country][state]['000'][t] > Airr_gov[country][state]['000'][t]:
                            if Airr[country][state]['000'][t] > 0:
                                factor = Airr_gov[country][state]['000'][t] / Airr[country][state]['000'][t]
                                for mun in state_mun_list[country][state].keys():
                                    smID = state+mun
                                    firr_fine_bc[t] = rescale_over_region(firr_fine_bc[t], map_fine[country], smID, factor)
                        else:
                            if Aplant[country][state]['000'][t]-Airr[country][state]['000'][t] > 0:
                                factor = (Aplant[country][state]['000'][t]-Airr_gov[country][state]['000'][t]) / (Aplant[country][state]['000'][t]-Airr[country][state]['000'][t])
                                for mun in state_mun_list[country][state].keys():
                                    smID = state+mun
                                    firr_fine_bc[t] = rescale_diff_over_region(firr_fine_bc[t], diff_firr[country][t], map_fine[country], smID, factor)


    # Aggregate to 0.0625 deg
    print('aggregating to 0.0625 deg')
    fplant_bc = np.empty([nTime,nLat,nLon])
    firr_bc = np.empty([nTime,nLat,nLon])
    for t in [mbc1, mbc2]:
        fplant_bc[t] = aggregate_space(fplant_fine_bc[t],25)
        firr_bc[t] = aggregate_space(firr_fine_bc[t],25)

    # Interpolation to other months
    print('interpolating to other months')
    tmpdata = np.empty([2*nTime,nLat,nLon])
    tmpdata[mbc1] = fplant_bc[mbc1]
    tmpdata[mbc2] = fplant_bc[mbc2]
    tmpdata[12:24] = fplant
    tmpdata = np.apply_along_axis(interpolate_months, 0, tmpdata)
    fplant_bc = tmpdata[0:12]
    tmpdata = np.empty([2*nTime,nLat,nLon])
    tmpdata[mbc1] = firr_bc[mbc1]
    tmpdata[mbc2] = firr_bc[mbc2]
    tmpdata[12:24] = firr
    tmpdata = np.apply_along_axis(interpolate_months, 0, tmpdata)
    firr_bc = tmpdata[0:12]

    # Make sure data stay within valid bounds
    for t in range(nTime):
        fplant_bc[t] = np.minimum(fplant_bc[t], fcrop)
        fplant_bc[t] = np.maximum(fplant_bc[t], 0)
        firr_bc[t] = np.minimum(firr_bc[t], fplant_bc[t])
        firr_bc[t] = np.maximum(firr_bc[t], 0)

    # Fill gaps with spatial interpolation
    mask_crop = np.where(fcrop > 0, 1, 0)
    sigma = 5
    for t in range(nTime):
        mask_fplant = np.where(np.isfinite(fplant_bc[t]), 1, 0)
        mask_firr = np.where(np.isfinite(firr_bc[t]), 1, 0)
        mask_fplant_firr = mask_fplant * mask_firr
        gapmask = np.where(mask_crop - mask_fplant_firr > 0, 1, 0)
        fplant_bc[t, gapmask == 1] = np.nan
        firr_bc[t, gapmask == 1] = np.nan
        fplant_bc[t] = gapfill(fplant_bc[t], mask_fplant_firr, gapmask, sigma)
        firr_bc[t] = gapfill(firr_bc[t], mask_fplant_firr, gapmask, sigma)

    # Make sure gap-filled data stay within valid bounds
    for t in range(nTime):
        fplant_bc[t] = np.where(np.isnan(fplant_bc[t]), 0, fplant_bc[t])
        fplant_bc[t] = np.minimum(fplant_bc[t], fcrop)
        fplant_bc[t] = np.maximum(fplant_bc[t], 0)
        firr_bc[t] = np.where(np.isnan(firr_bc[t]), 0, firr_bc[t])
        firr_bc[t] = np.minimum(firr_bc[t], fplant_bc[t])
        firr_bc[t] = np.maximum(firr_bc[t], 0)

    # Write output file
    month = np.arange(nTime, dtype=np.int) + 1
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

    for varname in ['lat','lon']:
        ds[varname].encoding['zlib'] = True

    ds['fplant'] = (['month','lat','lon'], fplant_bc)
    ds['firr'] = (['month','lat','lon'], firr_bc)

    print('writing to',outfile)
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()
