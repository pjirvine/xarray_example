"""
This file contains useful functions for processing netcdf files.

get_time_slice() - This function returns a time slice from the CMIP6 CEDA archive, first concatenating the files in the selected data folder together.
get_seasonal_mean_std() - This function returns 2 xarray datasets for the time-mean and standard deviation calculated over the selected years for the selected season.
get_fixed() - This function returns the fixed area and land fraction variables

"""

import numpy as np
import xarray as xr
import os

def get_time_slice(dates, model, centre, var, domain, exp, project, run, grid, time_files=0):
    """
    This function returns a time slice from the CMIP6 CEDA archive, first concatenating the files in the selected data folder together.
    USAGE:
    dates [list] = ['2070-01-01','2100-01-01'] - selects the time-slice to calculate over
    model, centre, var, domain, exp, project, run [string] - specifies which files to 
    time_files [integer] = 0, 1, 2... - by default [0] all files will be concatenated before the time-slice is extracted. 
        this may be time-consuming for long experiments. If you know that your time slice spans, e.g. only the last 2 files then enter 2.
    """
    
    """
    Collect the files needed from ceda_dir 
    """
    # define the data directory to search for files using lists
    ceda_dir='/badc/cmip6/data/CMIP6/{project}/{centre}/{model}/{exp}/{run}/{domain}/{var}/{grid}/latest/'.format(project=project, centre=centre, var=var, domain=domain, model=model, exp=exp, run=run, grid=grid)
    
    try:
        ceda_dir_files = os.listdir(ceda_dir)
    except FileNotFoundError as error:
        print(error)
        return
    
    if time_files == 0: # select all files
        file_list = ceda_dir_files
    elif time_files < len(ceda_dir_files): # select number of time files chosen
        file_list = ceda_dir_files[len(ceda_dir_files)-time_files:]
    else:
        print("ERROR: time_files >= length of file list")
        return
    
    """
    Concatenate files, if necessary, then select time
    """
    if len(file_list) == 1: # if there's just 1 file open it.
        ds = xr.open_dataset(ceda_dir + file_list[0])
    elif len(file_list) > 1: # if there's a list, open all of them and concatenate them together.
        ds_list = [ xr.open_dataset(ceda_dir + idx) for idx in file_list ]
        ds = xr.concat(ds_list, 'time')
    else:
        print("ERROR: file_list:", file_list)
        return
    
    # select time period and return
    try:
        ds_slice = ds.sel(time=slice(dates[0],dates[1]))
        return ds_slice # this line returns the output. 
    except ValueError as error:
        print("error in dates:",dates)
        print(error)
        return
# end def.

def get_seasonal_mean_std(season, dates, data_dir, model, centre, var, domain, exp, project, run, grid, time_files=0, over_write=False):
    """
    This function returns 2 xarray datasets for the time-mean and standard deviation calculated over the selected years for the selected season.
    A netcdf copy of this processes file will be saved in the ~/data/ folder. The function will check to see if the file is already present before 
    processing the raw files again, saving time.
    
    USAGE:
    season [string] = 'ANN', 'DJF', 'MAM', 'JJA', or 'SON' - selects the season to calculate
    dates [list] = ['2070-01-01','2100-01-01'] - selects the time-slice to calculate over
    data_dir = '/home/users/<USERNAME>/data/' - enter your username here!
    model, centre, var, domain, exp, project, run [string] - specifies which files to 
    time_files [integer] = 0, 1, 2... - by default [0] all files will be concatenated before the time-slice is extracted. 
        this may be time-consuming for long experiments. If you know that your time slice spans, e.g. only the last 2 files then enter 2.
    """
    
    """
    First define directories and filenames
    """
    # define the data directory to search for files using lists
    ceda_dir='/badc/cmip6/data/CMIP6/{project}/{centre}/{model}/{exp}/{run}/{domain}/{var}/{grid}/latest/'.format(project=project, centre=centre, var=var, domain=domain, model=model, exp=exp, run=run, grid=grid)
    # define the base of the filename for the output file(s) using lists, note this matches the file format in the data directory except for the date_4_file section and the missing .nc. 
    out_base='{var}_{domain}_{model}_{exp}_{run}_{grid}_{time_range}'.format(var=var, domain=domain, model=model, exp=exp, run=run, grid=grid, time_range=dates[0]+'_'+dates[1])
    # define a simplified folder structure for our output data directory.
    data_dir_full=data_dir + '{model}/{exp}/{var}/'.format(model=model, exp=exp, var=var)
    
    # Format the output filenames
    fname_mean = out_base + '_' + season + '_mean.nc'
    fname_std = out_base + '_' + season + '_std.nc'
    
    # specify the full output file paths
    fpath_mean = os.path.join(data_dir_full,fname_mean)
    fpath_std = os.path.join(data_dir_full,fname_std)
    
    """
    Check if processed files already exist, if so return those and end.
    """
    if os.path.isfile(fpath_mean) and os.path.isfile(fpath_std) and not over_write: # and not over_write:
        ds_mean = xr.open_dataset(fpath_mean)
        ds_std = xr.open_dataset(fpath_std)
        print("loading existing files", out_base, season)
        return ds_mean, ds_std
    
    print("processing files", out_base, season)
    
    # make directory to store output netcdfs
    os.makedirs(data_dir_full, exist_ok=True)
    
    """
    Use get_time_slice() to collect the needed data and take the time-slice needed.
    """
    args=[dates,model,centre,var,domain,exp,project,run,grid,time_files]
    ds_tslice = get_time_slice(*args) # *args passes the list of values in as arguments to the get_time_slice function.
    
    """
    Take seasonal mean and standard deviation
    """
    season_list = ['DJF','MAM','JJA','SON']
    if season == 'ANN':
        # calculate annual mean and standard deviation
        ds_yearly = ds_tslice.groupby('time.year').mean(dim='time') # take mean over every year
        ds_seas_mean = ds_yearly.mean(dim='year')
        ds_seas_std = ds_yearly.std(dim='year')

    elif season in season_list:
        ds_seasonal = list(ds_tslice.groupby('time.season')) # split into seasons
        # take mean over each season to produce a timeseries
        ds_seasonal_series = [ idx[1].groupby('time.year').mean('time') for idx in ds_seasonal if idx[0] == season ]
        # calculate mean and standard deviation across this seasonal timeseries.
        ds_seas_mean = ds_seasonal_series[0].mean('year')
        ds_seas_std = ds_seasonal_series[0].std('year')
    else:
        print("only ANN, DJF, MAM, JJA, or SON allowed, you entered:", season)
        return
    
    """
    Finally, save output to netcdf to use again later and also return result
    """
    # output datasets to netcdf
    ds_seas_mean.to_netcdf(fpath_mean)
    ds_seas_std.to_netcdf(fpath_std)
    # return seasonal (or annual) mean and standard deviation
    return ds_seas_mean, ds_seas_std

# end def

def get_ens_seasonal_mean_std(season, dates, data_dir, model, centre, var, domain, exp, project, runs, grid, time_files=0, over_write=False):
    """
    This function returns 2 xarray datasets for the time-mean and standard deviation calculated over the selected years and runs for the selected season.
    it is the same as get_seasonal_mean_std except that a list of runs must be provided rather than a single 
    
    USAGE:
    season [string] = 'ANN', 'DJF', 'MAM', 'JJA', or 'SON' - selects the season to calculate
    dates [list] = ['2070-01-01','2100-01-01'] - selects the time-slice to calculate over
    data_dir = '/home/users/<USERNAME>/data/' - enter your username here!
    model, centre, var, domain, exp, project [string] - specifies which files to work on
    runs [list of strings] - specifies the runs to combine for the ensemble-mean
    time_files [integer] = 0, 1, 2... - by default [0] all files will be concatenated before the time-slice is extracted. 
        this may be time-consuming for long experiments. If you know that your time slice spans, e.g. only the last 2 files then enter 2.
    """
    
    """
    First define directories and filenames
    """
    # define the data directory but leave {run} undefined
    ceda_dir='/badc/cmip6/data/CMIP6/{project}/{centre}/{model}/{exp}/{run}/{domain}/{var}/{grid}/latest/'.format(project=project, centre=centre, var=var, domain=domain, model=model, exp=exp, run='{run}', grid=grid)
    # define the base of the filename for the output ensemble-mean file(s).
    out_base='{var}_{domain}_{model}_{exp}_ens-mean_{grid}_{time_range}'.format(var=var, domain=domain, model=model, exp=exp, grid=grid, time_range=dates[0]+'_'+dates[1])
    # define a simplified folder structure for our output data directory.
    data_dir_full=data_dir + '{model}/{exp}/{var}/'.format(model=model, exp=exp, var=var)
    
    # Format the output filenames
    fname_mean = out_base + '_' + season + '_mean.nc'
    fname_std = out_base + '_' + season + '_std.nc'
    
    # specify the full output file paths
    fpath_mean = os.path.join(data_dir_full,fname_mean)
    fpath_std = os.path.join(data_dir_full,fname_std)
    
    """
    Check if processed files already exist, if so return those and end.
    """
    if os.path.isfile(fpath_mean) and os.path.isfile(fpath_std) and not over_write: # and not over_write:
        ds_mean = xr.open_dataset(fpath_mean)
        ds_std = xr.open_dataset(fpath_std)
        print("loading existing files", out_base, season)
        return ds_mean, ds_std
    
    print("processing files", out_base, season)
    
    # make directory to store output netcdfs
    os.makedirs(data_dir_full, exist_ok=True)
    
    """
    Use get_time_slice() to collect the needed data and take the time-slice needed.
    """
    
    # test if runs are a list.
    if type(runs) is not list:
        print('runs must be a list. you entered:',runs)
        return 
    
    # Run get_time_slice() once for each run.
    args_list=[[dates,model,centre,var,domain,exp,project,IDX,grid,time_files] for IDX in runs] # make X lists of input arguments
    tslices = [get_time_slice(*ARGS) for ARGS in args_list] # make a list of X timeslice outputs
    
    # Combine the X runs into an ensemble dataset along a new "run" dimension
    ds_ens = xr.concat(tslices, "run")
    
    """
    Take seasonal mean and standard deviation
    """
    season_list = ['DJF','MAM','JJA','SON']
    if season == 'ANN':
        # calculate annual mean and standard deviation
        ds_yearly = ds_ens.groupby('time.year').mean(dim='time') # take mean over every year
        ds_seas_mean = ds_yearly.mean(dim=['year','run']) # take mean over years and runs
        ds_seas_std = ds_yearly.std(dim=['year','run'])

    elif season in season_list:
        ds_seasonal = list(ds_ens.groupby('time.season')) # split into seasons
        # take mean over each season to produce a timeseries
        ds_seasonal_series = [ idx[1].groupby('time.year').mean('time') for idx in ds_70_100_seasons if idx[0] == season ]
        # calculate mean and standard deviation across this seasonal timeseries.
        ds_seas_mean = ds_seasonal_series[0].mean(dim=['year','run']) # take mean over years and runs
        ds_seas_std = ds_seasonal_series[0].std(dim=['year','run'])
    else:
        print("only ANN, DJF, MAM, JJA, or SON allowed, you entered:", season)
        return
    
    """
    Finally, save output to netcdf to use again later and also return result
    """
    # output datasets to netcdf
    ds_seas_mean.to_netcdf(fpath_mean)
    ds_seas_std.to_netcdf(fpath_std)
    # return seasonal (or annual) mean and standard deviation
    return ds_seas_mean, ds_seas_std

# end def

def get_fixed(centre, model, run, grid='gn'):
    """
    This function returns the fixed area and land fraction variables:
    areacella - area of gridcells
    sftlf - fraction of each grid cell that is land
    
    You'll need to find which run is the pre-industrial baseline. Search here: https://esgf-index1.ceda.ac.uk/search/cmip6-ceda/
    search for 'areacella' and the model name.
    """
    
    project='CMIP'
    exp='piControl'
    domain='fx'
    
    # define the folder and filenames
    ceda_dir='/badc/cmip6/data/CMIP6/{project}/{centre}/{model}/{exp}/{run}/{domain}/{var}/{grid}/latest/'
    fname='{var}_{domain}_{model}_{exp}_{run}_{grid}.nc'
    
    # Get area
    var = 'areacella'
    dir_path = ceda_dir.format(project=project, centre=centre, var=var, domain=domain, model=model, exp=exp, run=run, grid=grid)
    f_path = fname.format(var=var, domain=domain, model=model, exp=exp, run=run, grid=grid)
    ds_area = xr.open_dataset(dir_path + f_path)
    
    # Get land fraction
    var = 'sftlf'
    dir_path = ceda_dir.format(project=project, centre=centre, var=var, domain=domain, model=model, exp=exp, run=run, grid=grid)
    f_path = fname.format(var=var, domain=domain, model=model, exp=exp, run=run, grid=grid)
    ds_land = xr.open_dataset(dir_path + f_path)
    
    return ds_area, ds_land
# end def