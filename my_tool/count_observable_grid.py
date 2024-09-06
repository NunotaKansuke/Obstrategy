#!/usr/bin/env python
# coding: utf-8

from matplotlib import rcParams
rcParams["font.size"] = 15
rcParams["axes.linewidth"] = 3
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.major.size'] = 10
rcParams['xtick.major.width'] = 2
rcParams['xtick.minor.size'] = 5
rcParams['xtick.minor.width'] = 2
rcParams['ytick.major.size'] = 10
rcParams['ytick.major.width'] = 2
rcParams['ytick.minor.size'] = 5
rcParams['ytick.minor.width'] = 2

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
from datetime import datetime, timedelta, timezone
from matplotlib.dates import HourLocator, MinuteLocator, DateFormatter
import pytz
aflica = pytz.timezone('Africa/Johannesburg')
from astropy.coordinates import get_sun
from source.GetTobedone import get_observable_grid
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m) 
now = datetime(2024, 9, 4, 1, 0)
utc_now = datetime(2024, 9, 4, 23, 0, 0, tzinfo=timezone.utc)

tmp = get_observable_grid(utc_now,alt_min=25,alt_max=87)

current_time = utc_now 

count = []
time_list = []
for i in range(100):
   tmp = get_observable_grid(current_time,alt_min=25,alt_max=87)
   count.append(len(tmp))
   time_list.append(current_time+timedelta(hours=2))
   current_time += timedelta(minutes=3)

plt.plot(time_list,count)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

plt.gcf().autofmt_xdate()

plt.xlabel("Time")
plt.ylabel("Count")
plt.minorticks_on()
plt.show()

print(count)