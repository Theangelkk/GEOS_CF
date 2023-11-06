# conda activate GEOS_CF_Google

# Libreries
import os
from datetime import datetime, timedelta
import xarray as xr
from copy import copy
import argparse
import warnings

warnings.filterwarnings("ignore")

def joinpath(rootdir, targetdir):
    return os.path.join(os.sep, rootdir + os.sep, targetdir)

# Path of EEA_data
path_main_dir_GEOS_CF_data = os.environ['GEOS_CF_data']

if path_main_dir_GEOS_CF_data == "":
    print("Error: set the environmental variables of GEOS_CF_data")
    exit(-1)

def valid_datetime(dt):
    for fmt in ('%Y-%m-%dT%H:%M', '%Y-%m-%dT%H:%M:%S',
                '%Y-%m-%d %H:%M', '%Y-%m-%d %H:%M:%S'):
        try:
            return datetime.strptime(dt, fmt)
        except ValueError:
            pass
    raise argparse.ArgumentTypeError("Invalid date: '{0}'.".format(dt))

def valid_date(d):
    t = 'T00:00'
    return valid_datetime(d + t)

parser = argparse.ArgumentParser(description='Merge of NetCDF of GEOS CF data')
parser.add_argument('-s_date', '--start_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
parser.add_argument('-e_date', '--end_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
args = vars(parser.parse_args())

# Temporal interval to consider
start_datetime = args["start_date"]
end_datetime = args["end_date"]

# List of air pollutants to consider
list_air_pollutants = ["NO2", "O3", "CO", "SO2", "PM25_RH35_GCC", "PM25_RH35_GOCART"]

path_dir_data = joinpath(path_main_dir_GEOS_CF_data, "Google_data")
path_dir_formula_data = joinpath(path_dir_data, "Formula_denisity")
path_dir_fixed_data = joinpath(path_dir_data, "Fixed_denisity")

# For each air pollutant
for air_chem in list_air_pollutants:

    current_path_dir_air_formula_data = joinpath(path_dir_formula_data, air_chem)
    current_path_dir_air_fixed_data = joinpath(path_dir_fixed_data, air_chem)

    esito = False

    current_date = start_datetime.date()
    delta = timedelta(hours=24)

    if current_date.month + 1 == 13:
        end_current_date = datetime(current_date.year + 1, 1, 1, 0, 0).date()
    else:
        end_current_date = datetime(current_date.year, current_date.month + 1, 1, 0, 0).date()

    while esito == False:

        current_netcdf = None

        string_current_date = str(current_date.year) + "-" + str(current_date.month).zfill(2)

        PATH_CURRENT_DAY_FORMULA_DIR = joinpath(current_path_dir_air_formula_data, string_current_date) 
        PATH_CURRENT_DAY_FIXED_DIR = joinpath(current_path_dir_air_fixed_data, string_current_date) 

        diff_dates = end_current_date - current_date
        diff_dates_days = int(diff_dates.total_seconds() / (60*60*24))

        current_date_day = current_date
    
        for time in range(diff_dates_days):
        
            path_current_formula_netcdf_file = joinpath(PATH_CURRENT_DAY_FORMULA_DIR, str(current_date_day.day).zfill(2) + ".nc")
            path_current_fixed_netcdf_file = joinpath(PATH_CURRENT_DAY_FIXED_DIR, str(current_date_day.day).zfill(2) + ".nc")

            if time == 0:
                current_formula_netcdf = copy(xr.open_dataset(path_current_formula_netcdf_file))
                current_fixed_netcdf = copy(xr.open_dataset(path_current_fixed_netcdf_file))

                print(current_date_day.isoformat() + " analysed")
            
            else:
                xr_formula_ds = xr.open_dataset(path_current_formula_netcdf_file)
                xr_fixed_ds = xr.open_dataset(path_current_fixed_netcdf_file)

                current_formula_netcdf = xr.merge([current_formula_netcdf, xr_formula_ds])
                current_fixed_netcdf = xr.merge([current_fixed_netcdf, xr_fixed_ds])

            print(current_date_day.isoformat() + " analysed")

            current_date_day += delta

        path_formula_netcdf = joinpath(PATH_CURRENT_DAY_FORMULA_DIR, string_current_date + ".nc")
        path_fixed_netcdf = joinpath(PATH_CURRENT_DAY_FIXED_DIR, string_current_date + ".nc")
    
        if os.path.exists(path_formula_netcdf):
            os.remove(path_formula_netcdf)

        if os.path.exists(path_fixed_netcdf):
            os.remove(path_fixed_netcdf)

        current_formula_netcdf.to_netcdf(path_formula_netcdf)
        current_fixed_netcdf.to_netcdf(path_fixed_netcdf)

        current_date = end_current_date

        if current_date >= end_datetime.date():
            current_date = datetime(end_datetime.year, end_datetime.month, 1, 0, 0).date()
            end_current_date = end_datetime.date()
            esito = True
        else:
            if current_date.month + 1 == 13:
                end_current_date = datetime(current_date.year + 1, 1, 1, 0, 0).date()
            else:
                end_current_date = datetime(current_date.year, current_date.month + 1, 1, 0, 0).date()