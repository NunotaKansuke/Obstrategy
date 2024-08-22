#!/usr/bin/env python
# coding: utf-8

#import requests
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
#from astroplan import Observer
import time
from astropy.coordinates import AltAz, SkyCoord
#import re
#from astroplan.plots import plot_sky
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib.dates import HourLocator, MinuteLocator, DateFormatter
import pytz
aflica = pytz.timezone('Africa/Johannesburg')
from astropy.coordinates import get_sun
from GetTobedone import get_observable_grid
import pandas as pd
import sys

#------------------------------------------------------------------------------------#

def time2degrees(times):
    return times[0] * 360/24 + times[1] * 360/24/60 + times[2] * 360/24/60/60

def dms2degrees(dms):
    return dms[0] + (dms[1] / 60) + (dms[2] / 60/60)

def str2float(string):
    num  = string.split(":")
    num = [float(n) for n in num]
    if num[0]<0:
        num[1],num[2] = -num[1], -num[2] 
    return tuple(num)

def convert_ra_hour_to_deg(ra):
    [ra_hour,ra_min,ra_sec] = ra.split(":")
    ra_hour,ra_min,ra_sec = float(ra_hour),float(ra_min),float(ra_sec)
    ra_deg = (ra_hour + ra_min / 60 + ra_sec / 3600) * 15
    return ra_deg

def convert_dec_dms_to_deg(dec):
    [dec_deg, dec_min, dec_sec] = dec.split(":")
    dec_deg, dec_min, dec_sec = float(dec_deg),float(dec_min),float(dec_sec)
    _dec_deg = abs(dec_deg)  # Ensure positive value for calculation
    result = _dec_deg + dec_min / 60 + dec_sec / 3600
    if dec_deg < 0:
        result *= -1
    return result

def convert_altaz_xy(alt,az):
    alt_rad, az_rad = np.deg2rad(alt), np.deg2rad(az)
    x = np.abs(np.cos(alt_rad))*np.cos(az_rad+np.pi/2)
    y = np.abs(np.cos(alt_rad))*np.sin(az_rad+np.pi/2)
    return x, y, alt_rad>0

#------------------------------------------------------------------------------------#

class MakeGridScript():
    def __init__(self,margins=10):
        #----deffine constant-------------------#
        self.ll, self.ul = 25, 87
        self.obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m) 
        self.now = datetime.now(aflica).replace(tzinfo=None)+timedelta(minutes=margins)
        self.utc_now = self.now + timedelta(hours=-2)

    def get_altaz(self,ra,dec,time): #引数はutc
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        altaz_frame = AltAz(obstime=time, location=self.obs_loc)
        altaz_coords = coords.transform_to(altaz_frame)
        alt = altaz_coords.alt.degree
        az = altaz_coords.az.degree
        return alt, az
    
    def make_order_grid(self,hour=4,alt_min=30,alt_max=87):
        start_time = self.utc_now
        end_time = start_time + timedelta(hours=hour)
        current_time = start_time
        num_list = []
        time_list = []
        alt_list = []
        renzoku = 0
        renzoku_2 = 0
        while True:
            if renzoku == 0:
                grids = get_observable_grid(current_time,alt_min=alt_min,alt_max=alt_max)[3:]
            if renzoku>=len(grids) and renzoku_2 == 0:
                renzoku = 0
                renzoku_2 += 1 #renzoku_2がないとここで無限ループに入る時がある。それは、その時間に何もとるものがない時
                continue
            if renzoku>=len(grids) and renzoku_2 == 1: #その時間に何もとるものがない時
                print(f"there is no observable object at {(current_time+timedelta(hours=2)).strftime('%H:%M')}")
                current_time += timedelta(minutes=4.5)
                renzoku = 0
                continue
            grid = grids[renzoku]
            alt, az = self.get_altaz(grid["ra"],grid["dec"],current_time)
            if (alt < alt_min) or (alt > alt_max):
                renzoku = 0
                continue
            renzoku += 1
            if not grid["num"] in num_list: 
                num_list.append(grid["num"])
                time_list.append(current_time+timedelta(hours=2))
                alt_list.append(alt)
            current_time += timedelta(minutes=4.5)
            if current_time > end_time:
                break

        return num_list, time_list, alt_list
    
    def make_script_grid(self,path,hour=4,alt_min=30,alt_max=87):
        order, obs_time, alt_list = self.make_order_grid(hour=hour,alt_min=alt_min,alt_max=alt_max)
        grid_script = pd.read_csv("data/obsable_all_sky_grid-3.csv")
        output = pd.DataFrame(columns=grid_script.columns)
        for i,time,alt in zip(order,obs_time,alt_list):
            grid_script.loc[i - 1, "Comment1"] = time.strftime('%H:%M')
            grid_script.loc[i - 1, "Comment2"] = alt
            row = grid_script.iloc[i-1]
            output = pd.concat([output, row.to_frame().T], ignore_index=True)

        output.to_csv(path,index=False)


#----------------------------------------------------------#
if __name__ == "__main__":
    tmp = MakeGridScript(margins=10)
    today = datetime.now(aflica).replace(tzinfo=None).strftime('%Y_%m_%d')
    tmp.make_script_grid(f"./all_sky_grid/all_sky_grid{today}.csv",alt_min=25,hour=5)