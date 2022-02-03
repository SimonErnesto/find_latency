# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import mne
import matplotlib.pyplot as plt 
from scipy.stats import mode

os.chdir("")

evokeds = mne.read_evokeds("example_p3b_ave.fif")[0]


#the code for computing the fractional area is based on this post: https://lindeloev.net/hej-verden/

#find latencies with fractional peak or fractional area methods from mne evokeds and start and end points of time-window to search
def find_latency(evokeds, tmin, tmax, method, percents=None, out=None):
    sfreq = evokeds.info['sfreq'] 
    basetime = int(abs(evokeds.baseline[0]*sfreq))
    times = np.array(evokeds.times)[basetime:]
    elecs = np.array(evokeds.info.ch_names)
    samps = evokeds.get_data()[:,basetime:]*1e6 #scale to microvolts
    nelecs, nsamps = samps.shape
    start = round(tmin*sfreq)
    end = round(tmax*sfreq)
    if percents == None:
        percents = [0.25, 0.75]
    onset = []
    offset = []
    peak = []
    amp = []
    for e in range(nelecs):
        s = samps[e][start:end]
        if method=='frac': #compute fractional peak (50%)
            peak_idx = start + np.absolute(s).argmax() #index of maximum peak (global max of preselcted time-window)
            if peak_idx == start:
                start = start-2
            if peak_idx == end:
                end = end+2
            half_peak = np.absolute(samps[e][start:end]).max()*0.5 #half of the maximum peak
            onset_idx = start + np.absolute(samps[e][start:peak_idx]-half_peak).argmin() #index of latency onset (from left to right, so it need to be added to start)
            offset_idx = peak_idx + np.absolute(samps[e][peak_idx:end]-half_peak).argmin() #index of latency offset (from left to right, so it need to be added to peak latency)
            onset.append(times[onset_idx]) 
            offset.append(times[offset_idx])
            peak.append(times[peak_idx])
            amp.append(samps[e][peak_idx]) #amplitude at peak (global minimum)
        if (method=='area_pos') or (method=='area' and s.mean() > 0): #if the most common amplitude of the ERP time-window is positive, the cumulative sum is corrected by samples - min(smaples)
            peak_idx = start + np.absolute(np.cumsum(s-min(s))-sum(s-min(s))*0.5).argmin() #index of area midpoint obtained via cumulative sum method
            onset_idx = start + (np.absolute(np.cumsum(s-min(s))-sum(s-min(s))*percents[0])).argmin() #index of area lower percent
            offset_idx = start + (np.absolute(np.cumsum(s-min(s))-sum(s-min(s))*percents[1])).argmin() #index of area upper percent
            onset.append(times[onset_idx]) 
            offset.append(times[offset_idx])
            peak.append(times[peak_idx])
            amp.append(samps[e][peak_idx]) #amplitude at peak (fractional area)
        if (method=='area_neg') or (method=='area' and s.mean() < 0): #if the most common amplitude of the ERP time-window is negative, the cumulative sum is corrected by samples - max(smaples)
            peak_idx = start + np.absolute(np.cumsum(s-max(s))-sum(s-max(s))*0.5).argmin() #index of area midpoint obtained via cumulative sum method
            onset_idx = start + (np.absolute(np.cumsum(s-max(s))-sum(s-max(s))*percents[0])).argmin() #index of area lower percent
            offset_idx = start + (np.absolute(np.cumsum(s-max(s))-sum(s-max(s))*percents[1])).argmin() #index of area upper percent
            onset.append(times[onset_idx]) 
            offset.append(times[offset_idx])
            peak.append(times[peak_idx])
            amp.append(samps[e][peak_idx]) #amplitude at peak (fractional area)
        if method=='area+frac' and s.mean() > 0: #if the most common amplitude of the ERP time-window is positive, the cumulative sum is corrected by samples - min(smaples)
            peak_idx = start + np.absolute(np.cumsum(s-min(s))-sum(s-min(s))*0.5).argmin() #peak index from whole epoch
            half_peak = abs(samps[e][peak_idx])*0.5 #half of the maximum peak
            ons = start + np.absolute(samps[e][start:peak_idx]-half_peak).argmin() #index of latency onset (from left to right, so it need to be added to start)
            off = peak_idx + np.absolute(samps[e][peak_idx:end]-half_peak).argmin() #index of latency offset (from left to right, so it need to be added to peak latency)
            onset.append(times[ons]) 
            offset.append(times[off])
            peak.append(times[peak_idx])
            amp.append(samps[e][peak_idx]) #amplitude at peak (global maximum)
        if method=='area+frac' and s.mean() < 0: #if the most common amplitude of the ERP time-window is negative, the cumulative sum is corrected by samples - max(smaples)
            peak_idx = start + np.absolute(np.cumsum(s-max(s))-sum(s-max(s))*0.5).argmin() #peak index from whole epoch
            half_peak = abs(samps[e][peak_idx])*0.5 #half of the maximum peak
            ons = start + np.absolute(samps[e][start:peak_idx]-half_peak).argmin() #index of latency onset (from left to right, so it need to be added to start)
            off = peak_idx + np.absolute(samps[e][peak_idx:end]-half_peak).argmin() #index of latency offset (from left to right, so it need to be added to peak latency)
            onset.append(times[ons]) 
            offset.append(times[off])
            peak.append(times[peak_idx])
            amp.append(samps[e][peak_idx]) #amplitude at peak (global maximum)
    onset = np.array(onset)
    offset = np.array(offset)
    peak = np.array(peak)
    amp = np.array(amp)
    if out == 'dataframe':
       return  pd.DataFrame({'electrode':elecs, 'onset':onset, 'offset':offset, 'peak':peak, 'peak_amp':amp})
    else:
        return elecs,onset,offset,peak,amp 


####Create dataframes with results

latencies_area =  find_latency(evokeds, tmin=0.25, tmax=0.75, method='area', out='dataframe')

latencies_frac = find_latency(evokeds, tmin=0.25, tmax=0.75, method='frac', out='dataframe')

latencies_area_pos = find_latency(evokeds, tmin=0.25, tmax=0.75, method='area_pos', out='dataframe')

latencies_area_frac = find_latency(evokeds, tmin=0.25, tmax=0.75, method='area+frac', out='dataframe')

   
###Vizualise results    
start = 0.25
end = 0.75         
sfreq = evokeds.info['sfreq'] 
basetime = int(abs(evokeds.baseline[0]*sfreq))
times = np.array(evokeds.times)[basetime:]
elecs = np.array(evokeds.info.ch_names)
samps = evokeds.get_data()[:,basetime:]*1e6
nelecs, nsamps = samps.shape
start = round(start*sfreq)
end = round(end*sfreq)
erp = samps[12][start:end] #12 corresponds to Pz

pz_times = latencies_area[latencies_area.electrode=='Pz']
pk = int(round(pz_times['peak']*256).values[0])
ons = int(round(pz_times['onset']*256).values[0])
off = int(round(pz_times['offset']*256).values[0])
plt.axvline(0, color="k", linestyle="-", alpha=0.5)
plt.axhline(0, color="k", linestyle="-", alpha=0.5)
plt.plot(evokeds.times, evokeds.get_data()[12]*1e6, color="m", label="Evoked")
plt.fill_between(times[start:end], erp[:], erp[0], color='tan', alpha=0.3, label="Analysis time-window")
plt.vlines(pz_times['peak'].values[0], ymin=erp[0], ymax=samps[12][pk], color="k", linestyle="--", alpha=0.5, label="Area midpoint (fractional area peak)")
plt.vlines(pz_times['onset'].values[0], ymin=erp[0], ymax=samps[12][ons], color="k", linestyle=":", alpha=0.3, label="25% onset/offset")
plt.vlines(pz_times['offset'].values[0], ymin=erp[0], ymax=samps[12][off], color="k", linestyle=":", alpha=0.3)
plt.title("Fractional Area (Pz)")
plt.legend()
plt.ylabel("Amplitude (\u03BCV)")
plt.xlabel("Time (seconds)")
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.savefig('frac_area.png', dpi=300)
plt.close()

pz_times = latencies_frac[latencies_frac.electrode=='Pz']
pk = int(round(pz_times['peak']*256).values[0])
ons = int(round(pz_times['onset']*256).values[0])
off = int(round(pz_times['offset']*256).values[0])
plt.axvline(0, color="k", linestyle="-", alpha=0.5)
plt.axhline(0, color="k", linestyle="-", alpha=0.5)
plt.plot(evokeds.times, evokeds.get_data()[12]*1e6, color="c", label="Evoked")
plt.fill_between(times[start:end], erp[:], erp[0], color='tan', alpha=0.3, label="Analysis time-window")
plt.vlines(pz_times['peak'].values[0], ymin=erp[0], ymax=samps[12][pk], color="k", linestyle="--", alpha=0.5, label="Maximum amplitude (peak)")
plt.vlines(pz_times['onset'].values[0], ymin=erp[0], ymax=samps[12][ons], color="k", linestyle=":", alpha=0.3, label="Half-peak onset/offset")
plt.vlines(pz_times['offset'].values[0], ymin=erp[0], ymax=samps[12][off], color="k", linestyle=":", alpha=0.3)
plt.title("Fractional Peak (Pz)")
plt.legend()
plt.ylabel("Amplitude (\u03BCV)")
plt.xlabel("Time (seconds)")
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.savefig('frac_peak.png', dpi=300)
plt.close()