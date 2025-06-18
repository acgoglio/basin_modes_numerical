
# -------------------------------------------
# Extract the main 8 modes at each grid point
# -------------------------------------------

import glob
import sys
import xarray as xr
import netCDF4 as nc
import numpy as np
from scipy.signal import periodogram
import heapq
from functools import partial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.signal import detrend
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
mpl.use('Agg')

def amp_main_modes(lat_idx,lon_idx,ssh_ts_all,dt):

   ########
   
   # Flag and threshold [h] for filtering the spectrum the threshold is also used as plot minimum 
   flag_filter='true'
   th_filter=39

   ###################
   
   # Convert SSH time series to NumPy array and remove NaNs
   time_series_point = np.array(ssh_ts_all)
   valid_indices = np.logical_not(np.isnan(time_series_point))
   time_series_clean = time_series_point[valid_indices]
   
   if len(time_series_clean) == 0:
        return np.full(8, np.nan), np.full(8, np.nan)

   #### SPECTRUM ANALYSIS

   # Compute FFT
   spt_len = len(time_series_clean)
   fft = np.fft.fft(time_series_clean)
   freq = np.fft.fftfreq(spt_len, d=dt)

   # Select only positive frequencies (first half of FFT, excluding Nyquist frequency)
   half_spt_len = spt_len // 2
   freq_positive = freq[:half_spt_len]  # Only positive frequencies
   fft_positive = fft[:half_spt_len]  # Only corresponding FFT values

   # Ensure we are not including zero frequency
   mask = freq_positive > 0
   freq_positive = freq_positive[mask]
   fft_positive = fft_positive[mask]

   # Compute Power Spectrum Density (normalized)
   spt = (np.abs(fft_positive) ** 2) / spt_len #**2

   # Compute Periods in hours
   periods = 1 / freq_positive / 3600  # Convert periods to hours

   # Compute Amplitudes (correct scaling)
   amplitudes = (2 / spt_len) * np.abs(fft_positive)

   if flag_filter=='true':
      #print ('Filter = true')
      # Apply high-pass filter: Set frequencies below the threshold to zero
      high_pass_threshold = 1 / (th_filter * 3600)  # Corresponding to 48 hours in Hz
      fft_positive[freq_positive < high_pass_threshold] = 0  # Filter out frequencies below the threshold

      # Recompute the power spectrum and amplitudes after filtering
      spt = (np.abs(fft_positive) ** 2) / spt_len #**2
      amplitudes = (2 / spt_len) * np.abs(fft_positive)

   spt_smooth = spt
   amp_smooth = amplitudes

   # Found peaks in spt and in amplitude:
   amp_peaks, _ = find_peaks(amp_smooth)
   amp_peak_frequencies = freq_positive[amp_peaks]
   amp_peak_amplitudes = amp_smooth[amp_peaks]

   # Order based on amplitudes
   sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
   amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

   # Remove close modes:
   period_sorted = 1/amp_peak_frequencies_sorted/3600  # periodi dei modi
   amplitude_sorted = amp_peak_amplitudes_sorted  # ampiezze dei modi

   final_periods = period_sorted
   final_amplitudes = amplitude_sorted

#   final_periods = []
#   final_amplitudes = []
#
#   for i in range(len(period_sorted)):
#       if period_sorted[i] > 28: #30:
#           tolerance = 10  
#       elif 25 <= period_sorted[i] <= 28: #30:
#           tolerance = 5 #10  
#       elif 12 <= period_sorted[i] < 25:
#           tolerance = 2  
#       elif 6 <= period_sorted[i] < 12:
#           tolerance = 0.5   
#       else:
#           tolerance = 0.1  
# 
#   keep_mode = True  
#   for j in range(len(final_periods)):
#        if abs(period_sorted[i] - final_periods[j]) <= tolerance:
#            if amplitude_sorted[i] <= final_amplitudes[j]:
#                keep_mode = False  
#                break
#            else:
#                final_periods[j] = period_sorted[i]
#                final_amplitudes[j] = amplitude_sorted[i]
#                keep_mode = False  
#                break
#    
#   if keep_mode:
#        final_periods.append(period_sorted[i])
#        final_amplitudes.append(amplitude_sorted[i])

   # Final arrays
   amp_peak_amplitudes_sorted=np.array(final_amplitudes) 
   amp_peak_period_sorted=np.array(final_periods)

   #############
   # Return the values

   amp_peak_periods_main=[]
   amp_peak_amplitudes_main=[]

   for i in range(0,8):
     try:
       amp_peak_periods_main.append(amp_peak_period_sorted[i])
       amp_peak_amplitudes_main.append(amp_peak_amplitudes_sorted[i])
     except:
       amp_peak_periods_main.append(np.nan)
       amp_peak_amplitudes_main.append(np.nan)
   
   return np.array(amp_peak_periods_main),np.array(amp_peak_amplitudes_main)

