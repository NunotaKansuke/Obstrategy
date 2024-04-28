#!/usr/bin/env python
# coding: utf-8

import requests
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
from astroplan import Observer
import time
from astropy.coordinates import AltAz, SkyCoord
import re
from astroplan.plots import plot_sky
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib.dates import HourLocator, MinuteLocator, DateFormatter
import pytz
aflica = pytz.timezone('Africa/Johannesburg')

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

class Obstrategy():
    def __init__(self,path):
        gb_field_path = "./data//PRIME_LB_20230719_nunota.dat"
        self.gb_field = np.genfromtxt(gb_field_path,names=["name","ra","dec","l","b","type"],dtype=[("num", int), ("ra", float), ("dec", float), ("l", float), ("b", float), ("type", int)])
        
        self.list = np.genfromtxt(path,usecols=[0,1,2],names=["name","ra","dec"],encoding="utf-8",dtype=None)
        self.now = datetime.now(aflica).replace(tzinfo=None)
        self.obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m)

    def get_altaz(self,ra,dec,time):
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        altaz_frame = AltAz(obstime=time, location=self.obs_loc)
        altaz_coords = coords.transform_to(altaz_frame)
        alt = altaz_coords.alt.degree
        az = altaz_coords.az.degree
        return alt, az


tmp = Obstrategy("./data/test_list.dat")
print(tmp.get_altaz(18.989, -78.922,tmp.now))






