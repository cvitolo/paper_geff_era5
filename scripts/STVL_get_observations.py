# This is also available as notebook here:
# https://github.com/cvitolo/GEFF_notebooks/blob/master/GetObservations.ipynb

import datetime as dt
import numpy
import pandas as pd
import metview
import calendar
from ecmp.MediaPy.StvlRetrieve import StvlRetrieveToGeopoints

# Define parameters of interest and related accumulation periods
# https://confluence.ecmwf.int/display/VER/stvl+service
# Parameters:
# total precipitation in the last 24h (mm), wind speed (m/s), 2m temperature (K), 2m dewpoint temperature (K)
parameters = ['tp', '10ff', '2t', '2d']
# Accumulation period (in seconds)
obs_period = [24*3600, 0, 0, 0]

# The STVL service is limited to 2000 requests, I had to split the total number of requests into monthly chunks
for month in range(1, 13):
    # Find last day of month to set the time window to 1 calendar month
    month_as_string = str(format(month, '02d'))
    last_of_month = str(calendar.monthrange(2017, month)[1])
    period = pd.date_range('2017-' + month_as_string + '-01T00:00:00.000Z',
                           '2017-' + month_as_string + '-' + last_of_month + 'T23:00:00.000Z',
                           freq = 'H')
    print(period)

    # Initialise empty list
    appended_data = []

    # Retrieve observations from STVL
    for param, obst in zip(parameters, obs_period):
        print(param, obst)
        for geo in StvlRetrieveToGeopoints(table = 'observation',
                                           parameter = param,
                                           period = obst,
                                           reference_datetimes = period):
            data = geo.to_dataframe()
            data['param'] = param
            # store DataFrame in list
            appended_data.append(data)

    # Concatenate list to dataframe
    appended_data = pd.concat(appended_data)
    # Write DataFrame to a csv file
    appended_data.to_csv('/hugetmp/SYNOP/Observations_2017' + month_as_string + '.csv')
