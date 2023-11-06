# conda activate GEOS_CF_Google

# Important information of installation "ee":
# /Users/angelocasolaro/miniconda3/envs/GEOS_CF_Google/lib/python3.8/site-packages/ee
# l'import di io --> from io import StringIO

# Libreries
import os
import sys
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import argparse
import warnings

# Library of Google Earth Engine
import ee

warnings.filterwarnings("ignore")

# Path of GEOS_CF_data
path_main_dir_GEOS_CF_data = os.environ['GEOS_CF_data']

if path_main_dir_GEOS_CF_data == "":
    print("Error: set the environmental variables of GEOS_CF_data")
    exit(-1)

if not os.path.exists(path_main_dir_GEOS_CF_data):
  os.mkdir(path_main_dir_GEOS_CF_data)

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

parser = argparse.ArgumentParser(description='Script for downloading GEOS CF data set')
parser.add_argument('-s_date', '--start_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
parser.add_argument('-e_date', '--end_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
args = vars(parser.parse_args())

# Temporal interval to consider
start_datetime = args["start_date"]
end_datetime = args["end_date"]

# List of all air pollutants (Volume Mixing ratio)
# NO2, SO2, CO, O3: mol mol^-1
# PM25_RH35_GCC: ug m^-3
# PM25_RH35_GOCART: kg m^-3
list_air_pollutants = ["NO2", "O3", "CO", "SO2", "PM25_RH35_GCC", "PM25_RH35_GOCART"]

# Italy coordinates
lat_italy_bnds, lon_italy_bnds = [32.0,50.0], [5.0,21.0]

# Dictionary of molar mass (g/mol) for each air pollutants considered
# http://stcorp.github.io/harp/doc/html/algorithms/definitions.html
dict_molar_mass = {}

dict_molar_mass["AIR"] = 28.9644
dict_molar_mass["CO"] = 28.0101
dict_molar_mass["NO2"] = 46.00550
dict_molar_mass["O3"] = 47.99820
dict_molar_mass["SO2"] = 64.0638

# Trigger the authentication flow.
ee.Authenticate()

# Initialize the library.
ee.Initialize()

geos_cf_ee = ee.ImageCollection('NASA/GEOS-CF/v1/rpl/tavg1hr')

def joinpath(rootdir, targetdir):
    return os.path.join(os.sep, rootdir + os.sep, targetdir)

def ee_array_to_df(arr, list_of_bands):
    
    """Transforms client-side ee.Image.getRegion array to pandas.DataFrame."""
    df = pd.DataFrame(arr)

    # Rearrange the header.
    headers = df.iloc[0]
    df = pd.DataFrame(df.values[1:], columns=headers)

    # Remove rows without data inside.
    df = df[['longitude', 'latitude', 'time', *list_of_bands]]

    # Convert the data to numeric values.
    for band in list_of_bands:
        df[band] = pd.to_numeric(df[band], errors='coerce')

    # Convert the time field into a datetime.
    df['datetime'] = pd.to_datetime(df['time'], unit='ms')

    # Keep the columns of interest.
    df = df[['datetime', 'longitude', 'latitude', *list_of_bands]]
    df = df.set_index(['datetime', 'longitude', 'latitude'])

    return df

# Extraction of Italy measures from ImageCollection
def ext_dataset_ee_italy_for_day(geos_cf_ee, measure, start_day):

    global lat_italy_bnds, lon_italy_bnds

    # Dimensions of GEOS-CF cell in meters
    scale_factor = 27750

    geometry = ee.Geometry.Polygon(
            [[  [lon_italy_bnds[0], lat_italy_bnds[1]],
                [lon_italy_bnds[0], lat_italy_bnds[0]],
                [lon_italy_bnds[1], lat_italy_bnds[0]],
                [lon_italy_bnds[1], lat_italy_bnds[1]]]], None, False)

    delta_day_hours = timedelta(hours=24)

    # Day after of current day considered
    end_day = start_day + delta_day_hours
    dataset_day_ee = geos_cf_ee.filterDate(ee.Date(start_day.isoformat()), ee.Date(end_day.isoformat()))

    # Extraction of measure specified
    measure_ee = dataset_day_ee.select(measure)
    
    measure_italy_ee = measure_ee.getRegion(geometry, scale=scale_factor).getInfo()

    return measure_italy_ee

def E_s_T(temp):
  
    eso = 6.1078

    c4 = 0.43884187 * 1e-8

    c = [   
            0.99999683, -0.90826951*1e-2, 0.78736169*1e-4, -0.61117958*1e-6, \
            c4, -0.29883885*1e-10, 0.21874425*1e-12, -0.17892321*1e-14, \
            0.11112018*1e-16, -0.30994571*1e-19
        ]

    p = (c[0]+temp*(c[1]+temp*(c[2]+temp* \
        (c[3]+temp*(c[4]+temp*(c[5]+temp*(c[6]+temp* \
        (c[7]+temp*(c[8]+temp*(c[9]))))))))))

    # Polynomial function
    E_s = eso / pow(p,8)

    return E_s

def compute_current_air_density(temperature_ds, rh_ds):

    temperature_ds_C = temperature_ds - 273.15
    
    # Specific constant of dry air: 287.05 J/(kgK)
    constant_dry_air = 287.05

    # Specific constant of water vapor: 461.495 J/(kgK)
    constant_h2o = 461.495

    # In order to compute the density of water vapor is
    # necessary to compute the partial pressure of water vapor, defined with 
    # the term P_v
    # P_v = Relative_umidity * E_s(T)
    P_v = rh_ds * E_s_T(temperature_ds_C)

    # Compute the density of water vapor (kg/m^3)
    # Formula: Atmosphere_pressure(hPa) / (constant_dry_air(J/(khK)) * Air_temp(K))
    h2o_density = P_v / (constant_h2o * temperature_ds)

    # Before to compute the air density, it is require to compute the pressure
    # of air dry that can be retrive as the difference of total pressure,
    # expresse in Pa, minus the pressure of water vapor.

    # Value of atmosphere pressure at model level 55 --> level 72 (GEOS CF) - about 288 meters
    surface_pressure_Pa = 985.15 * 100

    P_d = surface_pressure_Pa - P_v

    # Calcolo della consistenza dell'aria secca

    # Compute the density of air dry (kg/m^3)
    # Formula: Atmosphere_pressure(hPa) /(constant_dry_air(J/(khK)) * Air_temp(K))
    dry_air_density = P_d / (constant_dry_air * temperature_ds)

    # Formula of air density (kg/m^3)
    air_density = dry_air_density + h2o_density

    return air_density

path_dir_data = joinpath(path_main_dir_GEOS_CF_data, "Google_data")

if not os.path.exists(path_dir_data):
    os.mkdir(path_dir_data)

path_dir_formula_data = joinpath(path_dir_data, "Formula_denisity")

if not os.path.exists(path_dir_formula_data):
      os.mkdir(path_dir_formula_data)

path_dir_fixed_data = joinpath(path_dir_data, "Fixed_denisity")

if not os.path.exists(path_dir_fixed_data):
      os.mkdir(path_dir_fixed_data)

for air_chem in list_air_pollutants:

    path_dir_air_formula_data = joinpath(path_dir_formula_data, air_chem)
    path_dir_air_fixed_data = joinpath(path_dir_fixed_data, air_chem)

    if not os.path.exists(path_dir_air_formula_data):
        os.mkdir(path_dir_air_formula_data)
    
    if not os.path.exists(path_dir_air_fixed_data):
        os.mkdir(path_dir_air_fixed_data)

# Compute the difference between start and end time in seconds
diff_dates = end_datetime - start_datetime

# Compute the number of days to load by Google Engine
diff_dates_days = int(diff_dates.total_seconds() / (60*60*24))

delta = timedelta(hours=24)

# For each air pollutant to observe
for air_chem in list_air_pollutants:

    current_date = start_datetime

    res_fin = False
    idx_time = 0

    current_path_dir_air_formula_data = joinpath(path_dir_formula_data, air_chem)
    current_path_dir_air_fixed_data = joinpath(path_dir_fixed_data, air_chem)

    while res_fin == False:

        if idx_time >= diff_dates_days:
            res_fin = True
        else:
            try:
                currentDay = current_date.day
                currentMonth = current_date.month
                currentYear = current_date.year

                string_current_date = str(currentYear) + "-" + str(currentMonth).zfill(2)
                
                # It is a new day to load by Google Engine
                current_air_pollutant_italy_ee = ext_dataset_ee_italy_for_day(geos_cf_ee, air_chem, current_date)
                
                current_day_path_dir_air_formula_data = joinpath(current_path_dir_air_formula_data, string_current_date)
                current_day_path_dir_air_fixed_data = joinpath(current_path_dir_air_fixed_data, string_current_date)

                if not os.path.exists(current_day_path_dir_air_formula_data):
                    os.mkdir(current_day_path_dir_air_formula_data)
                
                if not os.path.exists(current_day_path_dir_air_fixed_data):
                    os.mkdir(current_day_path_dir_air_fixed_data)
            
                current_air_pollutant_ds = ee_array_to_df(current_air_pollutant_italy_ee, [air_chem]).to_xarray()

                # -------------------- Temperature in Kelvin ---------------
                current_temperature = ee_array_to_df(ext_dataset_ee_italy_for_day(geos_cf_ee, 'T', current_date), ['T']).to_xarray()

                # -------------------- Relative umidity --------------------
                current_relative_umidity = ee_array_to_df(ext_dataset_ee_italy_for_day(geos_cf_ee, 'RH', current_date), ['RH']).to_xarray()

                if air_chem == "NO2" or air_chem == "O3" or air_chem == "SO2" or air_chem == "CO":
                    
                    # --------------- Conversion from VMR (mol mol^-1) to MMR (kg kg^-1) ------------------
                    # Formula: MMR (kg kg^-1) = vmr (mol mol^-1) * (mass_molar_chm / mass_molar_air)
                    current_air_pollutant_mmr = current_air_pollutant_ds[air_chem] * \
                                                (dict_molar_mass[air_chem] / dict_molar_mass["AIR"])

                    # --------------- Conversion from MMR (kg kg^-1) to concentration (kg m^-3) -----------
                    # Formula: MMR * Densità dell'aria umità --> (kg kg^-1) * (kg m^-3) --> (kg m^-3)
                    current_surface_air_density_italy = compute_current_air_density(current_temperature['T'], current_relative_umidity['RH'])

                    current_air_pollutant_formula_kg_m3 = current_air_pollutant_mmr * current_surface_air_density_italy

                    air_density_fixed = 1.191403
                    current_air_pollutant_fixed_kg_m3 = current_air_pollutant_mmr * air_density_fixed
            
                elif air_chem == "PM25_RH35_GOCART":
                    current_air_pollutant_formula_kg_m3 = current_air_pollutant_ds[air_chem]
                    current_air_pollutant_fixed_kg_m3 = current_air_pollutant_ds[air_chem]
            
                elif air_chem == "PM25_RH35_GCC":
                    current_air_pollutant_formula = current_air_pollutant_ds[air_chem]
                    current_air_pollutant_fixed = current_air_pollutant_ds[air_chem]

                if  air_chem == "NO2" or air_chem == "O3" or air_chem == "SO2" \
                    or air_chem == "PM25_RH35_GOCART":

                    # -------------------- Conversion (kg m^-3) --> (ug m^-3) --------------------
                    current_air_pollutant_formula = current_air_pollutant_formula_kg_m3 * 1e+9
                    current_air_pollutant_fixed = current_air_pollutant_fixed_kg_m3 * 1e+9
        
                elif air_chem == "CO":

                    # -------------------- Conversion (kg m^-3) --> (mg m^-3) --------------------
                    current_air_pollutant_formula = current_air_pollutant_formula_kg_m3 * 1e+6
                    current_air_pollutant_fixed = current_air_pollutant_fixed_kg_m3 * 1e+6
            
                # Save NetCDF file
                path_netcdf_formula = joinpath(current_day_path_dir_air_formula_data, str(currentDay).zfill(2) + ".nc")
                path_netcdf_fixed = joinpath(current_day_path_dir_air_fixed_data, str(currentDay).zfill(2) + ".nc")

                current_air_pollutant_formula.to_netcdf(path_netcdf_formula)

                print(air_chem + " FORMULA " + current_date.isoformat() + " saved")

                current_air_pollutant_fixed.to_netcdf(path_netcdf_fixed)

                print(air_chem + " FIXED " + current_date.isoformat() + " saved")

                current_date += delta
                idx_time +=1

            except:
                print("Error download " + current_date.isoformat() + " retry...")