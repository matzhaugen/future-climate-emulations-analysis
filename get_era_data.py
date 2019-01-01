#!/usr/bin/env python
import calendar
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
 
ERA_DIR = '/project/moyer/mahaugen/reanalysis/era_interim_land'
TIME_WINDOW = "{}-01-01/to/{}-12-31"
def retrieve_interim():
    """      
       A function to demonstrate how to iterate efficiently over several years and months etc    
       for a particular interim_request.        
    """
    year_range = (1979, 2016)
    north_west_south_east = (50, -125, 30, -70)
    
    target = "{}/interim_daily_latlon{:2d}-{:3d}-{:2d}-{:2d}_years{:4d}-{:4d}.grb".format(ERA_DIR, 
        *(north_west_south_east + year_range))
    interim_request(
        TIME_WINDOW.format(year_range[0], year_range[1]), 
        north_west_south_east,
        target)
    
def interim_request(requestDates, north_west_south_east, target):
    """      
        An ERA interim request for analysis pressure level data.
        Change the keywords below to adapt it to your needs.
        (eg to add or to remove  levels, parameters, times etc)
        Request cost per day is 112 fields, 14.2326 Mbytes
    """
    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": requestDates,
        "expver": "1",
        # North West South East
        "area": "{:3d}/{:3d}/{:2d}/{:2d}".format(*north_west_south_east),
        "grid": "0.75/0.75",
        "levtype": "sfc",
        "param": "167.128",
        "step": "0",
        "stream": "oper",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "type": "an",
        "target": target,
    })
if __name__ == '__main__':
    retrieve_interim()