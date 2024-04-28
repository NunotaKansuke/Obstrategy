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
from astropy.coordinates import get_sun

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

#------------------------------------------------------------------------------------#

class Obstrategy():
    def __init__(self,path):
        self.ll, self.ul = 25, 87
        gb_field_path = "./data//PRIME_LB_20230719_nunota.dat"
        self.gb_field = np.genfromtxt(gb_field_path,names=["name","ra","dec","l","b","type"],dtype=[("num", int), ("ra", float), ("dec", float), ("l", float), ("b", float), ("type", int)])
        self.list = np.genfromtxt(path,usecols=[0,1,2],names=["name","ra","dec"],encoding="utf-8",dtype=None)
        self.now = datetime.now(aflica).replace(tzinfo=None)
        self.utc_now = self.now + timedelta(hours=-2)
        self.noon = self.now.replace(hour=12,minute=0,second=0,microsecond=0)
        self.obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m) 

    def get_altaz(self,ra,dec,time): #引数はutc
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        altaz_frame = AltAz(obstime=time, location=self.obs_loc)
        altaz_coords = coords.transform_to(altaz_frame)
        alt = altaz_coords.alt.degree
        az = altaz_coords.az.degree
        return alt, az
    
    def calc_night_start(self): #utcで返す
        start_time = Time(self.noon, format='datetime', scale='utc')
        ind_time = start_time

        while True:
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if sun_alt < -15:
                disappear_hour = ind_time
                break
            ind_time +=  timedelta(hours=1)

        while ind_time >= disappear_hour+timedelta(hours=-1):
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if not sun_alt < -15:
                appear_10minuite = ind_time
                break
            ind_time -= timedelta(minutes=10)
        
        while ind_time <= appear_10minuite + timedelta(minutes=10):
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if sun_alt < -15:
                return ind_time
            ind_time +=  timedelta(minutes=1)
        

    def which_appear(self,ra,dec,time):
        alt, az = self.get_altaz(ra,dec,time)
        return alt > self.ll and alt < self.ul

    def calc_appear(self,ra,dec,start_time,end_time): #utcで返す
        ind_time = start_time
        appear = False

        if self.which_appear(ra,dec,start_time):
            return start_time
        
        #二分探索的なやつ
        while ind_time < end_time: #あとで、このループの終わり方を考え直す
            if self.which_appear(ra,dec,ind_time):
                appear_hour = ind_time
                appear = True
                break
            ind_time +=  timedelta(hours=1)
        
        if not appear:
            return end_time

        while ind_time >= appear_hour+timedelta(hours=-1):
            if not self.which_appear(ra,dec,ind_time):
                appear_10minuite = ind_time
                break
            ind_time -= timedelta(minutes=10)
        
        while ind_time <= appear_10minuite + timedelta(minutes=10):
            if self.which_appear(ra,dec,ind_time):
                return ind_time
            ind_time +=  timedelta(minutes=1)


tmp = Obstrategy("./data/test_list.dat")
ra_array = np.array([10,20,30])
dec_array = np.array([10,20,30])
print(tmp.calc_appear(-50,-50,tmp.now,tmp.now+timedelta(hours=12)).replace(microsecond=0,second=0))
print(tmp.utc_now)
print(tmp.now)
print(tmp.calc_night_start())