#!/usr/bin/env python
# coding: utf-8

#import requests
import numpy as np
import pandas as pd
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
from astropy.coordinates import get_body
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

class CalcPoint():
    def __init__(self, path):
        #----deffine constant-------------------#
        self.ll, self.ul = 25, 87
        self.sun_night, self.sun_twillight= 0, -15
        self.obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m) 
        self.now = datetime.now(aflica).replace(tzinfo=None)
        self.utc_now = self.now + timedelta(hours=-2)
        self.noon = self.now.replace(hour=12,minute=0,second=0,microsecond=0)

        #----input script-----------------------#
        script = pd.read_csv(path)
        self.script = script
        self.script["RA"] = script["RA"].map(convert_ra_hour_to_deg)
        self.script["DEC"] = script["DEC"].map(convert_dec_dms_to_deg)

    def calc_moon_altaz(self,time):
        time = Time(time, format='datetime', scale='utc')
        moon = get_body("moon",time).transform_to(AltAz(obstime=time,location=self.obs_loc)) 

        return moon.alt, moon.az
    
    def get_altaz(self,ra,dec,time): #引数はutc
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        altaz_frame = AltAz(obstime=time, location=self.obs_loc)
        altaz_coords = coords.transform_to(altaz_frame)
        alt = altaz_coords.alt.degree
        az = altaz_coords.az.degree
        return alt, az


if __name__ == "__main__":

    test = CalcPoint("/Users/nunotahiroshisuke/Desktop/iral/work/Obstrategy/input_exa/input_example.csv")
    print(test.script["RA"])
    print(test.script["DEC"])