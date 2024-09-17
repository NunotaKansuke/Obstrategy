import requests
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.coordinates import AltAz, SkyCoord
from datetime import datetime
import re

import pytz
aflica = pytz.timezone('Africa/Johannesburg')

url1 = "http://primeenv.saao.ac.za/tobedoneHL.dat"
response = requests.get(url1) 

data = response.text

tmp =  data.split("\n")
del tmp[0:2]
del tmp[len(tmp)-1]


field_num_list, field_ra_list, field_dec_list = [], [], []
for line in tmp:
    line = line.strip().split()
    line[0] = re.sub(r'[^0-9]', '', line[0])
    field_num_list.append(int(line[0]))
    field_ra_list.append(float(line[1]))
    field_dec_list.append(float(line[2]))

field_num_list, field_ra_list, field_dec_list = np.array(field_num_list), np.array(field_ra_list), np.array(field_dec_list)
field_array = np.empty(field_num_list.shape[0], dtype=[("num", int), ("ra", float), ("dec", float),("alt",float),("az",float)])
field_array["num"], field_array["ra"], field_array["dec"] = field_num_list, field_ra_list, field_dec_list
field_array["ra"] *= 360/24
obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m)

def get_observable_grid_HL(time,alt_min=30,alt_max=87):#時間はutc
    
    obs_time = Time(time, format='datetime', scale='utc')
    
    coords = SkyCoord(field_array["ra"], field_array["dec"], frame='icrs', unit="deg")
    altaz_frame = AltAz(obstime=obs_time, location=obs_loc)
    altaz_coords = coords.transform_to(altaz_frame)
    
    field_array["alt"] = altaz_coords.alt.degree
    field_array["az"] = altaz_coords.az.degree
    
    cond_az1, cond_az2 = field_array["alt"] > alt_min, field_array["alt"] <alt_max
    ind = np.where(cond_az1&cond_az2)
    return field_array[ind]

if __name__ == "__main__":
    import datetime
    tmp = get_observable_grid_HL(datetime.datetime.now())
    print(len(tmp))
    print(tmp[np.argsort(tmp["alt"])][::-1])