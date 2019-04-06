#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt


# function to make a gaussian kernel
def make_kernel(radius, sigma):
    kxlen = 2*radius+1
    kylen = 2*radius+1
    kernel = np.zeros([kxlen,kylen])
    sum = 0
    for i in range(kxlen):
        for j in range(kylen):
            kernel[i,j] = np.exp( -((i-radius)**2+(j-radius)**2) / (2*sigma**2) )
            sum += kernel[i,j]
    for i in range(kxlen):
        for j in range(kylen):
            kernel[i,j] /= sum
    return kernel


# gapfilling with a gaussian kernel
def gapfill (mydata, opmask, data_with_buffer, data_count_with_buffer, radius, sigma):

    # Note: ideally we could simplify this function

    xlen = mydata.shape[1]
    ylen = mydata.shape[0]
    xlen_wbuf = data_with_buffer.shape[1]
    ylen_wbuf = data_with_buffer.shape[0]
    buf = int(0.5 * (ylen_wbuf - ylen))

    # make kernel
    kernel = make_kernel(radius, sigma)
    [kylen,kxlen] = kernel.shape

    # loop for filling isolated gaps with kernel-weighted spatial average
    a = np.where(opmask == True)[0]
    b = np.where(opmask == True)[1]
    for i,j in zip(a,b):

            # Fill tmpdata and tmpcount arrays with data in neighborhood surrounding this gap
            i0 = i-radius
            i1 = i+radius+1
            j0 = j-radius
            j1 = j+radius+1
            i0_wbuf = i0+buf
            if i0_wbuf < 0:
                i0_offset = -i0_wbuf
                i0_wbuf = 0
            else:
                i0_offset = 0
            i1_wbuf = i1+buf
            if i1_wbuf > ylen_wbuf:
                i1_wbuf = ylen_wbuf
            j0_wbuf = j0+buf
            if j0_wbuf < 0:
                j0_offset = -j0_wbuf
                j0_wbuf = 0
            else:
                j0_offset = 0
            j1_wbuf = j1+buf
            if j1_wbuf > xlen_wbuf:
                j1_wbuf = xlen_wbuf
            ilen = i1_wbuf-i0_wbuf
            jlen = j1_wbuf-j0_wbuf

            tmpdata = np.empty(kernel.shape)
            tmpdata[i0_offset:i0_offset+ilen,j0_offset:j0_offset+jlen] = data_with_buffer[i0_wbuf:i1_wbuf,j0_wbuf:j1_wbuf]
            tmpcount = np.empty(kernel.shape)
            tmpcount[i0_offset:i0_offset+ilen,j0_offset:j0_offset+jlen] = data_count_with_buffer[i0_wbuf:i1_wbuf,j0_wbuf:j1_wbuf]

            # Mask out any invalid data points
            tmpmask = 1 - np.isnan(tmpdata).astype(int)
            tmpdata = tmpdata * tmpmask.astype(float)
            tmpcount = tmpcount * tmpmask

            # If we have valid points to interpolate between, do interpolation
            weights = kernel * tmpcount
            wtdata = weights * tmpdata
            sumwts = np.nansum(weights)
            sumwtdata = np.nansum(wtdata)
            if (sumwts > 0):
                mydata[i,j] = sumwtdata / sumwts

    return mydata


def main():
    infile = ''
    varnamelist = ''
    fplant_annual_file = ''
    paramfile = ''
    cropclassidxlist = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:v:f:p:c:o:",["infile=","varnamelist=","fplant_annual_file=","paramfile=","cropclassidxlist=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -v <varnamelist> -f <fplant_annual_file> -p <paramfile> -c <cropclassidxlist> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -v <varnamelist> -f <fplant_annual_file> -p <paramfile> -c <cropclassidxlist> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-f", "--fplant_annual_file"):
            fplant_annual_file = arg
        elif opt in ("-p", "--paramfile"):
            paramfile = arg
        elif opt in ("-c", "--cropclassidxlist"):
            cropclassidxlist = arg
            cropclassidxs = cropclassidxlist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg

    data = {}

    # Open infile file
    ds = xr.open_dataset(infile)

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

    # Open fplant_annual file
    ds = xr.open_dataset(fplant_annual_file)

    fplant_annual = np.empty(ds['fplant_annual'].shape)
    fplant_annual[:] = ds['fplant_annual'][:]
    fplant_mask = np.where(fplant_annual > 0, 1, 0)
    wt = np.where(fplant_annual > 1, 1, fplant_annual)

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
    tmp = fcrop_mask - fplant_mask
    missing_mask = np.where(tmp == 1, True, False)
    radius = 16
    sigma = 3

    # Add buffer around data
    nLat_wbuf = nLat + 2 * radius
    nLon_wbuf = nLon + 2 * radius
    data_wbuf = np.full([nTime,nLat_wbuf,nLon_wbuf],np.nan)
    data_count_wbuf = np.zeros([nTime,nLat_wbuf,nLon_wbuf])
    first = True
    for varname in varnames:
        data_wbuf[:,radius:nLat+radius,radius:nLon+radius] = data[varname][:]
        for t in range(nTime):
            data_count_wbuf[t,radius:nLat+radius,radius:nLon+radius] = wt[:]

        # Fill gaps
        for t in range(nTime):
            data[varname][t] = gapfill(data[varname][t],missing_mask,data_wbuf[t],data_count_wbuf[t],radius,sigma)
        if first:
            data_mask = np.where(data[varname][0] >= 0, 1, 0)
            tmp = fcrop_mask - data_mask
            missing_mask_final = np.where(tmp == 1, True, False)
            a = np.where(missing_mask_final == True)[0]
            b = np.where(missing_mask_final == True)[1]
            for i,j in zip(a,b):
                data[varname][:,i,j] = 0
            first = False

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
    for varname in varnames:
        ds[varname] = (['month','lat','lon'],data[varname])
    ds['missing_initial'] = (['lat','lon'],missing_mask)
    ds['missing_final'] = (['lat','lon'],missing_mask_final)
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

