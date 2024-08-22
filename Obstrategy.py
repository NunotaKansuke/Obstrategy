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

class Obstrategy():
    def __init__(self,path):
        #----load data--------------------------#
        gb_field_path = "./data/PRIME_LB_20230719_nunota.dat"
        self.gb_field = np.genfromtxt(gb_field_path,names=["name","ra","dec","l","b","type","x","y"]
                                      ,dtype=[("num", int),("ra", float),("dec", float),("l", float),("b", float),("type", int),("x", int),("y", int)])
        self.list = np.genfromtxt(path,usecols=[0,1,2],names=["name","ra","dec"],encoding="utf-8",dtype=None)
        #----deffine constant-------------------#
        self.ll, self.ul = 25, 87
        self.sun_night, self.sun_twillight= 0, -15
        self.obs_loc = EarthLocation(lat=-32.3763*u.deg, lon=20.8107*u.deg, height=1798*u.m) 
        self.now = datetime.now(aflica).replace(tzinfo=None)
        self.utc_now = self.now + timedelta(hours=-2)
        self.noon = self.now.replace(hour=12,minute=0,second=0,microsecond=0)
        self.pre = False #wether prepare for plot
        self.main = np.where(self.gb_field["type"]==1) #gbのメインフィールド
        self.sub = np.where(self.gb_field["type"]==2) #gbのメインの両端のフィールド
        self.moa = np.where(self.gb_field["type"]==3) #moaのフィールド
        self.other = np.where(self.gb_field["type"]==0) #その他のフィールド
        self.m_s = np.union1d(self.main, self.sub) #メインと両端
        self.m_s_m = np.union1d(self.m_s, self.moa) #メインと両端とmoa
        self.gb_mode = self.m_s #何を使いますかってモード。今はメインと両端とmoaにしてる。引数で変えれるようにする。
        #----run function-----------------------#
        self.calc_obs_start()
        self.calc_obs_end()
        #self.set_observable_time()
        self.set_observable_time_gb()
        #---------------------------------------#

    def get_altaz(self,ra,dec,time): #引数はutc
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        altaz_frame = AltAz(obstime=time, location=self.obs_loc)
        altaz_coords = coords.transform_to(altaz_frame)
        alt = altaz_coords.alt.degree
        az = altaz_coords.az.degree
        return alt, az
    
    def calc_obs_start(self): #utcで返す
        ind_time = Time(self.noon, format='datetime', scale='utc')

        while True:
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if sun_alt < self.sun_twillight:
                disappear_hour = ind_time
                break
            ind_time +=  timedelta(hours=1)

        while ind_time >= disappear_hour+timedelta(hours=-1):
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if not sun_alt < self.sun_twillight:
                appear_10minuite = ind_time
                break
            ind_time -= timedelta(minutes=10)
        
        while ind_time <= appear_10minuite + timedelta(minutes=10):
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if sun_alt < self.sun_twillight:
                self.obs_start = ind_time.value
            ind_time +=  timedelta(minutes=1)

    def calc_obs_end(self): #utcで返す
        ind_time = Time(self.obs_start+timedelta(hours=+3), format='datetime', scale='utc')

        while True:
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if sun_alt > self.sun_twillight:
                disappear_hour = ind_time
                break
            ind_time +=  timedelta(hours=1)

        while ind_time >= disappear_hour+timedelta(hours=-1):
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if not sun_alt > self.sun_twillight:
                appear_10minuite = ind_time
                break
            ind_time -= timedelta(minutes=10)
        
        while ind_time <= appear_10minuite + timedelta(minutes=10):
            sun_altaz = get_sun(ind_time).transform_to(AltAz(obstime=ind_time,location=self.obs_loc)) 
            sun_alt = sun_altaz.alt.value
            if sun_alt > self.sun_twillight:
                self.obs_end = ind_time.value
            ind_time +=  timedelta(minutes=1)
        
    def which_appear(self,ra,dec,time):
        alt, az = self.get_altaz(ra,dec,time)
        return alt > self.ll and alt < self.ul

    def calc_appear(self,ra,dec,start_time,end_time): #utcで返す
        if self.which_appear(ra,dec,start_time):
            return start_time
        if not self.which_appear(ra,dec,end_time):
            return end_time
        
        #二分探索
        time1, time2 = start_time, end_time
        while abs(time1-time2)>timedelta(minutes=1):
            time3 = time1 + (time2-time1)/2
            if self.which_appear(ra,dec,time3):
                time2 = time3
            else:
                time1 = time3
        return time3
    
    def calc_disappear(self,ra,dec,start_time,end_time): #wutcで返す
        if not self.which_appear(ra,dec,start_time):
            return start_time
        if self.which_appear(ra,dec,end_time):
            return end_time
        
        #二分探索
        time1, time2 = start_time, end_time
        while abs(time1-time2)>timedelta(minutes=1):
            time3 = time1 + (time2-time1)/2
            if not self.which_appear(ra,dec,time3):
                time2 = time3
            else:
                time1 = time3
        return time3

    def set_observable_time(self):
        start_time_list, end_time_list = [], []

        for i in range(self.list.shape[0]):
            ra, dec = self.list["ra"][i], self.list["dec"][i]
            appear_time = self.calc_appear(ra,dec,self.obs_start,self.obs_end)
            disappear_time = self.calc_disappear(ra,dec,self.obs_start,self.obs_end)
            start_time_list.append(appear_time)
            end_time_list.append(disappear_time)

        new_data = np.zeros(self.list.shape[0], dtype=self.list.dtype.descr + [("start", "datetime64[s]"), ("end", "datetime64[s]")])
        for name in self.list.dtype.names:
            new_data[name] = self.list[name]
        new_data["start"] = start_time_list
        new_data["end"] = end_time_list

        self.list = new_data

    def set_observable_time_gb(self):
        start_time_list, end_time_list = [], []

        for i in range(self.gb_field.shape[0]):
            ra, dec = self.gb_field["ra"][i], self.gb_field["dec"][i]
            appear_time = self.calc_appear(ra,dec,self.obs_start,self.obs_end)
            disappear_time = self.calc_disappear(ra,dec,self.obs_start,self.obs_end)
            start_time_list.append(appear_time)
            end_time_list.append(disappear_time)

        new_data = np.zeros(self.gb_field.shape[0], dtype=self.gb_field.dtype.descr + [("start", "datetime64[s]"), ("end", "datetime64[s]")])
        for name in self.gb_field.dtype.names:
            new_data[name] = self.gb_field[name]
        new_data["start"] = start_time_list
        new_data["end"] = end_time_list

        self.gb_field = new_data
        self.gb_start = np.min(self.gb_field["start"][self.gb_mode]).astype(datetime)
        self.gb_end = np.max(self.gb_field["end"][self.gb_mode]).astype(datetime)

    def prepare_for_plot(self):
        time_list = []
        alt_array, az_array = [], []
        alt_array_gb, az_array_gb = [], []
        current_time = self.obs_start
        while current_time <= self.obs_end:
            time_list.append(current_time)
            alt, az = self.get_altaz(self.list["ra"],self.list["dec"],current_time)
            alt_gb, az_gb = self.get_altaz(self.gb_field["ra"],self.gb_field["dec"],current_time)
            current_time += timedelta(minutes=30)
            alt_array.append(alt)
            az_array.append(az)
            alt_array_gb.append(alt_gb)
            az_array_gb.append(az_gb)

        self.time_list = np.array(time_list) #utc時間
        self.alt_array, self.az_array = np.array(alt_array).T, np.array(az_array).T
        self.alt_array_gb, self.az_array_gb = np.array(alt_array_gb).T, np.array(az_array_gb).T
        self.pre = True

    def plot_alt(self):
        if not self.pre:
            self.prepare_for_plot()

        plt.figure(figsize=(6,5))

        for i in range(self.alt_array.shape[0]):
            plt.plot(self.time_list+timedelta(hours=+2),self.alt_array[i],lw=2,c="blue",alpha=0.5) #南ア時間でプロット
        for i in self.other[0]:
           plt.plot(self.time_list+timedelta(hours=+2),self.alt_array_gb[i],lw=2,c="orange",alpha=0.5) #南ア時間でプロット
        for i in self.m_s_m:
            plt.plot(self.time_list+timedelta(hours=+2),self.alt_array_gb[i],lw=2,c="red",alpha=0.5) #南ア時間でプロット
        
        plt.hlines(y=self.ll,xmin=self.time_list[0]+timedelta(hours=+2),xmax=self.time_list[-1]+timedelta(hours=+2),colors="cyan",linestyles="--",label="low limit= 22")
        plt.gca().xaxis.set_major_locator(HourLocator())
        plt.gca().xaxis.set_minor_locator(MinuteLocator(interval=60)) 
        plt.gca().xaxis.set_major_formatter(DateFormatter('%H')) 
        plt.ylabel("Altitude")
        plt.xlabel(f"Time, starting at {self.time_list[0].strftime('%Y-%m-%d')}")
        plt.ylim(0,90)
        plt.grid(True, linestyle='--', color='gray')
        plt.minorticks_on()
        plt.legend(loc="best")
        plt.show()

    def plot_altaz(self):
        if not self.pre:
            self.prepare_for_plot()
        
        plt.figure(figsize=(5,5))
        x, y, _ = convert_altaz_xy(self.alt_array.T,self.az_array.T)
        plt.plot(x,y,lw=1,c="blue",alpha=0.5)
        x, y, _ = convert_altaz_xy(self.alt_array_gb.T,self.az_array_gb.T)
        plt.plot(x.T[self.other].T,y.T[self.other].T,lw=1,c="orange",alpha=0.5)
        plt.plot(x.T[self.m_s_m].T,y.T[self.m_s_m].T,lw=1,c="red",alpha=0.5)

        #--------------plot circle----------------#
        theta = np.linspace(0, 2*np.pi, 100)
        x = np.cos(theta)
        y = np.sin(theta)
        plt.plot(x, y,c="black")
        x = np.cos(np.deg2rad(self.ul))*np.cos(theta)
        y = np.cos(np.deg2rad(self.ul))*np.sin(theta)
        plt.plot(x, y,c="C1")
        x = np.cos(np.deg2rad(self.ll))*np.cos(theta)
        y = np.cos(np.deg2rad(self.ll))*np.sin(theta)
        plt.plot(x, y, c="C1")
        #------------------------------------------#

        plt.axis('equal')
        plt.xticks([]) 
        plt.yticks([])
        plt.xlim(-1.1,1.1)
        plt.ylim(-1.1,1.1)
        plt.show()

    def set_grid(self):
        grids = get_observable_grid(self.obs_start) #self.obs_startの時点で上に上がっているgridだけを取得
        unique_values, indices = np.unique(grids["dec"], return_inverse=True) #decが同じやつ同士でまとめる

        grouped_indices = [np.where(indices == i)[0] for i in range(len(unique_values))if len(np.where(indices == i)[0]) >= 10] #decが同じグループの中身が10こ以上のものだけ採用
        last_disappear =  np.array([self.calc_disappear(grids["ra"][grouped_indices[i][-1]],grids["dec"][grouped_indices[i][-1]],
                                                        self.obs_start,self.gb_start) for i in range(len(grouped_indices))]) #そのグループで一番最後に観測されるgridが沈む時間を計算
        sort_ind = np.argsort(last_disappear) #その時間でソートする
        obs_orders = [grouped_indices[i] for i in sort_ind]
        obs_orders = np.array(np.concatenate([arr.flatten() for arr in obs_orders]))
        return grids[obs_orders]
    
    def set_gb(self):
        obs_gb = self.gb_field[self.gb_mode]
        num_obs = np.zeros(obs_gb.shape)
        num_obs[0] = 1
        self.gb_dt = timedelta(minutes=4) #gbのひとつの観測にかかる時間(後から調節する)
        obs_gb = obs_gb[np.argsort(obs_gb["start"])]
        order = [obs_gb["name"][0]]
        current = obs_gb[0] 
        ind_time = self.gb_start+self.gb_dt
        self.nanshume = [] #スクリプトの各行が何周目かを示す。あとで、1周目はHbandで2周目はJbandで撮るみたいな時に使う
        nanshume = 1
        while ind_time < self.obs_end: #観測できるかつ距離が短いかつできるだけ左上その中でまだ観測していない
            can_obs = np.where(ind_time>obs_gb["start"]) #まだ、gbが沈むのは考慮に入れていない
            next_cand = np.where(num_obs == np.min(num_obs[can_obs]))[0] #観測できるfieldの中で最小の観測回数と同様の観測回数をもつfield
            dist = np.abs(current["x"]-obs_gb["x"]) + np.abs(current["y"]-obs_gb["y"])
            dist_min = np.min(dist[np.intersect1d(can_obs, next_cand)])
            dist_ind = np.where(dist == dist_min)
            comb_ind = np.intersect1d(dist_ind,np.intersect1d(can_obs, next_cand))
            order.append(obs_gb["name"][comb_ind[0]])
            self.nanshume.append(nanshume)
            num_obs[comb_ind[0]] +=1
            ind_time += self.gb_dt
            if np.all(num_obs > 0): 
                current = obs_gb[0]
                num_obs=np.zeros(obs_gb.shape)
                nanshume += 1
            else:
                 current = obs_gb[comb_ind[0]]
        self.gb_order = order
        self.nanshume = np.array(self.nanshume)

    def make_script_gb(self,path,band=[]):
        self.set_gb()

        if len(band) > np.max(self.nanshume):
            sys.stderr.write(f"len(band) must be smaller than {np.max(self.nanshume)+1}\n")
            sys.exit(1)

        while len(band) < np.max(self.nanshume):
            band.append("H")

        gb_script = pd.read_csv("data/gb_script.csv")
        output = pd.DataFrame(columns=gb_script.columns)

        for i in self.gb_order:
            row = gb_script.iloc[i-1] 
            output = pd.concat([output, row.to_frame().T], ignore_index=True)

        for i in range(np.max(self.nanshume)):
            ind = np.where(self.nanshume==i+1)[0]
            output["Filter2"].iloc[ind] = band[i]
        
        output.to_csv(path,index=False)

    def make_script_grid(self,path):
        order = self.set_grid()
        grid_script = pd.read_csv("data/obsable_all_sky_grid-3.csv")
        output = pd.DataFrame(columns=grid_script.columns)
        for i in order["num"]:
            row = grid_script.iloc[i-1]
            output = pd.concat([output, row.to_frame().T], ignore_index=True)

        output.to_csv(path,index=False)


    def print_time(self):
        print(f"Start time for observation {self.obs_start+timedelta(hours=+2)}")
        print(f"Start time for bulge {self.gb_start+timedelta(hours=+2)}")
        print(f"End time {self.obs_end+timedelta(hours=+2)}")
            

    
#----------------------------------------------------------#
if __name__ == "__main__":
    tmp = Obstrategy("./data/test_list.dat")
    #print(tmp.calc_appear(-50,-50,tmp.now,tmp.now+timedelta(hours=12)).replace(microsecond=0,second=0))
    #print(tmp.calc_disappear(-50,-50,tmp.now+timedelta(hours=5),tmp.now+timedelta(hours=12)).replace(microsecond=0,second=0))
    #print(tmp.utc_now)
    #print(tmp.now)
    #print(tmp.obs_start,tmp.obs_end)
    #print(tmp.list)
    #print(tmp.gb_field[tmp.m_s_m])
    #tmp.plot_altaz()
    #tmp.plot_alt()
    #print(tmp.calc_disappear(16.159 ,-68.668,tmp.obs_start,tmp.obs_end))
    #print(tmp.calc_appear(16.159 ,-68.668,tmp.obs_start,tmp.obs_end))
    #tmp.make_script_grid("tmp.csv")
    #tmp.make_script_gb("tmp_gb.csv",["H","Y"])
    tmp.make_script_grid_from_now("tmp.csv")